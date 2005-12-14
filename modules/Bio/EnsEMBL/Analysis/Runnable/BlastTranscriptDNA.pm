# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA->
  new(
      -transcript => $transcript,
      -program => 'wublastn',
      -database => 'embl_vertrna',
      -options => 'hitdist=40 -cpus=1',
      -parser => $bplitewrapper,
      -filter => $featurefilter,
     );
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

This module acts as an intermediate between blast and transcripts. It
is primarily used by the runnabledb BlastGenscanDNA to blast
the protein sequence of a genscan or any other ab initio prediction 
againsts a dna database. It instantiates a blast runnable passing it the
Bio::Seq of the transcript translation as the query sequence. This
module expects all the same arguments as a standard blast with the 
exception or a query sequence as these must be passed into the blast 
runnable it instantiates

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Blast);


=head2 new

  Arg [1]         : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA
  Arg [Transcript]: Bio::EnsEMBL::Transcript
  Function        : create a BlastTranscriptDNA runnable 
  Returntype      : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA
  Exceptions      : none 
  Example         : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($pt) = rearrange(['TRANSCRIPT'], @args);
  $self->transcript($pt);
  return $self;
}



=head2 transcript

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA
  Arg [2]   : Bio::EnsEMBL::Transcript
  Function  : container for transcript
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: throws if not passed a Bio::EnsEMBL::Transcript
  Example   : 

=cut


sub transcript{
  my ($self, $pt) = @_;
  if($pt){
    throw("BlastGenscanDNA:transcript must be a ".
          "Bio::EnsEMBL::Transcript ") 
      unless($pt->isa("Bio::EnsEMBL::Transcript"));
    $self->{'prediction_transcript'} = $pt;
  }
  return $self->{'prediction_transcript'};
}


=head2 output

  Function  : override the output method to allow its
   resetting to empty

=cut

sub output {
  my ($self, $arr_ref, $reset) = @_;

  if(!$self->{'output'}){
    $self->{'output'} = [];
  }
  if($arr_ref){
    throw("Must pass Runnable:output an arrayref not a ".$arr_ref)
      unless(ref($arr_ref) eq 'ARRAY');
    if ($reset) {
      $self->{'output'} = $arr_ref;
    } else {
      push(@{$self->{'output'}}, @$arr_ref);
    }
  }
  return $self->{'output'};
}


=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA
  Arg [2]   : string, working directory
  Function  : instantiations Blast runnable and runs it then converts
  blast hits back into genomic coords
  Returntype: none
  Exceptions: none
  Example   : 

=cut

sub run{
  my ($self, $dir) = @_;

  $self->workdir($dir) if($dir);

  my $pep = $self->transcript->translate;
  $pep->id($self->transcript->dbID);

  if($pep->length <= 3){
    #transcripts this length cause problems for blast
    return;
  }

  my $query = $self->query;
  $self->query($pep);  
  $self->SUPER::run($dir);
  $self->query($query);

  my $out = $self->output;
  $self->output([], 1);
  $self->align_hits_to_query($out);
}



=head2 align_hits_to_query

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA
  Arg [2]   : arrayref of Bio::EnsEMBL::BaseAlignFeatures
  Function  : convert the features from blast from peptide coodinates
  to genomic coordinates
  Returntype: none
  Exceptions: 
  Example   : 

=cut


sub align_hits_to_query {
  my ( $self, $features )  = @_;
  
  # for each feature
  my @features = sort{ $a->start <=> $b->start} @$features;
  for my $feature ( @features ) {
    my %exon_hash = ();
    # for each ungapped piece in it
    my @ungapped = $feature->ungapped_features;
    for my $ugFeature ( @ungapped ) {
      my @split = $self->transcript->pep2genomic($ugFeature->start(),
                                              $ugFeature->end());
      
      my $cdna_total = 1;
      foreach my $gcoord ( @split ) {
        if($gcoord->isa('Bio::EnsEMBL::Mapper::Gap')) {
          $cdna_total += $gcoord->end - $gcoord->start + 1;
          next;
        }
        
        my $cdna_start = $cdna_total;
        my $gstart  = $gcoord->start;
        my $gend    = $gcoord->end;
        my $gstrand = $gcoord->strand;
        my $cdna_end = $gend - $gstart + $cdna_start;
        $cdna_total += $gend - $gstart + 1;
        
        #determine which exon this genomic coordinate overlaps
        my $exon;
        foreach my $e (@{$self->transcript->get_all_Exons}) {
          if($gstart >= $e->start && $gend <= $e->end) {
            $exon = $e;
            last;
          }
        }


        # first, eat away non complete codons from start
        while(( $cdna_start - 1 ) % 3 != 0 ) {
          $cdna_start++;
          if( $gstrand == 1 ) {
            $gstart++;
          } else {
            $gend--;
          }
        }
        
        # and from end
        while( $cdna_end  % 3 != 0 ) {
          $cdna_end--;
          if( $gstrand == 1 ) {
            $gend--;
          } else {
            $gstart++;
          }
        }
        
        if( $cdna_end <= $cdna_start ) {
          next;
        }
        #recalculate the hit coordinates.  They may be split up into
        # seperate features by introns
        my ($hstrand, $hstart, $hend);
        $hstrand = $feature->hstrand();
        if($hstrand == 1) {
          $hstart = $cdna_start - 1 + $ugFeature->hstart();
          $hend = $cdna_end - 1 + $ugFeature->hstart();
        } else {
          $hend = $ugFeature->hend() - $cdna_start + 1;
          $hstart = $ugFeature->hend() - $cdna_end + 1;
        }
        my $fp = $self->feature_factory->create_feature_pair
          ($gstart, $gend, $gstrand, $feature->score, $hstart,
           $hend, $hstrand, $feature->hseqname, $feature->percent_id, 
           $feature->p_value);
        #store generated feature pairs, hashed on exons
        push( @{$exon_hash{$exon}}, $fp );
      }
    }
    # Take the pieces for each exon and make gapped feature
    foreach my $ex ( keys %exon_hash ) {
      my $dna_align_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new
        (-features => $exon_hash{$ex});
      $self->output([$dna_align_feature]);
    }
  }

}
