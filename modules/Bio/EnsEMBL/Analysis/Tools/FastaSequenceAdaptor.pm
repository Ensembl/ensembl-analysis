=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::FastaSequenceAdaptor

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Tools::FastaSequenceAdaptor;

use strict;
use warnings;

use Bio::DB::HTS::Faidx;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence  qw(reverse_comp);

use parent qw(Bio::EnsEMBL::DBSQL::SequenceAdaptor);


=head2 fetch_by_Slice_start_end_strand

  Arg  [1]   : Bio::EnsEMBL::Slice slice
               The slice from which you want the sequence
  Arg  [2]   : (optional) int startBasePair
               The start base pair relative to the start of the slice. Negative
               values or values greater than the length of the slice are fine.
               default = 1
  Arg  [3]   : (optional) int endBasePair
               The end base pair relative to the start of the slice. Negative
               values or values greater than the length of the slice are fine,
               but the end must be greater than or equal to the start
               count from 1
               default = the length of the slice
  Arg  [4]   : (optional) int strand
               1, -1
               default = 1
  Example    : $dna = $seq_adptr->fetch_by_Slice_start_end_strand($slice, 1,
                                                                  1000, -1);
  Description: Retrieves from db the sequence for this slice
               uses AssemblyMapper to find the assembly
  Returntype : string
  Exceptions : endBasePair should be less or equal to length of slice
  Caller     : Bio::EnsEMBL::Slice::seq(), Slice::subseq()
  Status     : Stable

=cut

sub fetch_by_Slice_start_end_strand {
   my ( $self, $slice, $start, $end, $strand ) = @_;

   if(!ref($slice) || !($slice->isa("Bio::EnsEMBL::Slice") or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
     throw("Slice argument is required.");
   }

   $start = 1 if(!defined($start));


   if ( ( !defined($end) || $start > $end || $start < 0 || $end < 0 || $slice->start> $slice->end ) && $slice->is_circular ) {

       if ( !defined($end) || ($start > $end ) ) {
     return $self->_fetch_by_Slice_start_end_strand_circular( $slice, $start, $end, $strand );
       }

       if ( defined($end) && ($end < 0) ) {
     $end += $slice->seq_region_length;
       }

       if ($start < 0) {
           $start += $slice->seq_region_length;
       }

       if($slice->start> $slice->end) {
           return $self->_fetch_by_Slice_start_end_strand_circular( $slice, $slice->start, $slice->end, $strand );
       }
  }

  if ( ( !defined($end) ) && (not $slice->is_circular) ) {
           $end = $slice->end() - $slice->start() + 1;
  }

  if ( $start > $end ) {
      throw("Start must be less than or equal to end.");
  }

   $strand ||= 1;
   #get a new slice that spans the exact region to retrieve dna from
   my $right_expand  = $end - $slice->length(); #negative is fine
   my $left_expand   = 1 - $start; #negative is fine

   if($right_expand || $left_expand) {
     $slice = $slice->expand($left_expand, $right_expand);
   }

   #retrieve normalized 'non-symlinked' slices
   #this allows us to support haplotypes and PARs
   my $slice_adaptor = $slice->adaptor();
   my @symproj=@{$slice_adaptor->fetch_normalized_slice_projection($slice)};

   if(@symproj == 0) {
     throw('Could not retrieve normalized Slices. Database contains ' .
           'incorrect assembly_exception information.');
   }

   #call this method again with any slices that were 'symlinked' to by this
   #slice
   if(@symproj != 1 || $symproj[0]->[2] != $slice) {
     my $seq;
     foreach my $segment (@symproj) {
       my $symlink_slice = $segment->[2];
       #get sequence from each symlinked area
       $seq .= ${$self->fetch_by_Slice_start_end_strand($symlink_slice,
                                                        1,undef,1)};
     }
     if($strand == -1) {
       reverse_comp(\$seq);
     }
     return \$seq;
   }

   # we need to project this slice onto the sequence coordinate system
   # even if the slice is in the same coord system, we want to trim out
   # flanking gaps (if the slice is past the edges of the seqregion)
   my $seq = ${$self->_fetch_seq($slice->seq_region_name, $slice->start, $slice->length())};

   #if the sequence is too short it is because we came in with a seqlevel
   #slice that was partially off of the seq_region.  Pad the end with Ns
   #to make long enough
   if(length($seq) != $slice->length()) {
     $seq .= 'N' x ($slice->length() - length($seq));
   }

   if(defined($self->{_rna_edits_cache}) and defined($self->{_rna_edits_cache}->{$slice->get_seq_region_id})){
     $self->_rna_edit($slice,\$seq);
   }

   #if they asked for the negative slice strand revcomp the whole thing
   reverse_comp(\$seq) if($strand == -1);

   return \$seq;
}

=head2 _fetch_raw_seq

 Arg [1]    : String $id, accession of the sequence
 Arg [2]    : Int $start, 5' position of the sequence to fetch on the forward strand
 Arg [3]    : Int $length, length of the sequence to fetch
 Description: Retrieve the DNA sequence of the region specified in the arguments
 Returntype : StringRef $seq, the sequence could be a whole chromosome
 Exceptions : Throws if the fetcher object is not Bio::DB::HTS::Faidx or Bio::DB::Fasta

=cut

sub _fetch_raw_seq {
  my ($self, $id, $start, $length) = @_;

  my $fasta_db = $self->{_fasta_db};
  throw("No fasta file found for ".$self->dbc->dbname) unless ($fasta_db);
  my $fa_length = $fasta_db->length($id);
  my $seq;
  my $end = $start+$length-1;
  if ($fa_length and $fa_length > 0) {
    my $location_string = "$id:$start-$end";
    if($fasta_db->isa('Bio::DB::HTS::Faidx')) {
      ($seq, $length) = $fasta_db->get_sequence($location_string);
    }
    elsif($fasta_db->isa('Bio::DB::Fasta')) {
      $seq = $fasta_db->seq($location_string);
    }
    else {
      throw("ERROR: Don't know how to fetch sequence from a ".ref($fasta_db));
    }
    return \$seq;
  }
  else {
    throw("Could not find sequence for $id $start $end from $fasta_db $fa_length");
  }
}


=head2 fasta

 Arg [1]    : String $fasta, path to a fasta file
 Description: Create a fetcher object for a fasta file. It is prefered to use
              an indexed fasta file.
 Returntype : Bio::DB::HTS::Faidx or Bio::DB::Fasta depending on the presence
              of an index file
 Exceptions : Throws if the file extension is not '.fa' or '.fasta'

=cut

sub fasta {
  my ($self, $fasta) = @_;

  if ($fasta) {
    if ($fasta =~ /\.fa$/ or $fasta =~ /\.fasta$/) {
      if (-e "$fasta.fai") {
        $self->{_fasta_db} = Bio::DB::HTS::Faidx->new($fasta);
      }
      else {
        $self->{_fasta_db} = Bio::DB::Fasta->new($fasta);
      }
    }
    else {
      throw("'$fasta' should have a '.fa' or '.fasta' extension");
    }
  }
  return $self->{_fasta_db};
}


=head2 _populate_seq_region_edits

 Arg [1]    : None
 Description: Empty method to avoid any connection to the MySQL server to check
              for the presence of sequence edit. Only LRGs and maybe some other
              obscure species uses sequence edit.
 Returntype : None
 Exceptions : None

=cut

sub _populate_seq_region_edits {
  my ($self) = @_;

  return ;
}

1;

