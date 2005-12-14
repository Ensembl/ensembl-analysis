# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep->
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
is primarily used by the runnabledb BlastGenscanPep to blast
the protein sequence of a genscan or any other ab initio prediction 
againsts a protein database. It instantiates a blast runnable passing it 
the Bio::Seq of the transcript translation as the query sequence. This
module expects all the same arguments as a standard blast with the 
exception or a query sequence as these must be passed into the blast 
runnable it instantiates

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Blast);

=head2 new

  Arg [1]         : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
  Arg [Transcript]: Bio::EnsEMBL::Transcript
  Function        : create a BlastTranscriptPep runnable 
  Returntype      : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
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

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
  Arg [2]   : Bio::EnsEMBL::Transcript
  Function  : container for transcript
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: throws if not passed a Bio::EnsEMBL::Transcript
  Example   : 

=cut

sub transcript{
  my ($self, $pt) = @_;
  if($pt){
    throw("BlastGenscanPep:transcript must be a ".
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

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
  Arg [2]   : string, working directory
  Function  : instantiations Blast runnable and runs it then converts
  blast hits back into genomic coords
  Returntype: none
  Exceptions: none
  Example   : 

=cut

sub run{
  my ($self, $dir) = @_;

  my $pep = $self->transcript->translate;
  $pep->id($self->transcript->dbID);
  $self->query($pep);

  if($pep->length <= 3){
    #transcripts this length cause problems for blast
    return;
  }

  $self->SUPER::run($dir);

  my $out = $self->output;
  $self->output([], 1);
  $self->align_hits_to_query($out);
}



=head2 align_hits_to_query

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
  Arg [2]   : arrayref of Bio::EnsEMBL::BaseAlignFeatures
  Function  : convert the features from blast from peptide coodinates
  to genomic coordinates
  Returntype: none
  Exceptions: 
  Example   : 

=cut

sub align_hits_to_query {
  my ($self, $features)  = @_;

  my $ff = $self->feature_factory();
  my @output;
  for my $feature ( @$features ) {
    my %exon_hash = ();
    # for each ungapped piece in it
    foreach my $ugFeature ( $feature->ungapped_features() ) {
      my $cdna_total = 1;
      #convert peptide coords to genomic coords
      my @split = $self->transcript->pep2genomic($ugFeature->start(),
                                                 $ugFeature->end());
      foreach my $gcoord ( @split ) {
        if($gcoord->isa('Bio::EnsEMBL::Mapper::Gap')) {
          $cdna_total += $gcoord->end - $gcoord->start + 1;
          next;
        }
        
        my $gstart = $gcoord->start;
        my $gend   = $gcoord->end;
        my $gstrand = $gcoord->strand;
        my $cdna_start = $cdna_total;
        my $cdna_end = $cdna_start + $gend - $gstart;
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

        my $hstart = (($cdna_start+2)/3) + $ugFeature->hstart() - 1;
        my $hend = ($cdna_end / 3) + $ugFeature->hstart() - 1;
        my $fp = $ff->create_feature_pair($gstart, $gend, $gstrand, 
                                          $feature->score, $hstart, 
                                          $hend, 1, 
                                          $feature->hseqname, 
                                          $feature->percent_id, 
                                          $feature->p_value,
                                          undef,
                                          $self->query,
                                          $feature->analysis);
        push( @{$exon_hash{$exon}}, $fp );
      }
    }
    
    # Take the pieces for each exon and make gapped feature
    foreach my $ex ( keys %exon_hash ) {
      my $dna_align_feature = Bio::EnsEMBL::DnaPepAlignFeature->new
        (-features => $exon_hash{$ex});
      push(@output, $dna_align_feature);
    }
  }
  $self->output(\@output);
}
