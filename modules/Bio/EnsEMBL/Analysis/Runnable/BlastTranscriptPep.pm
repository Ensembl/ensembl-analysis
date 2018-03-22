=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep - 

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

=cut


package Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable::Blast');


=head2 new

 Arg [Transcript]: Bio::EnsEMBL::Transcript
 Description     : create a BlastTranscriptPep runnable 
 Returntype      : Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
 Exceptions      : None 

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($pt) = rearrange(['TRANSCRIPT'], @args);
  $self->transcript($pt);
  return $self;
}


=head2 transcript

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: container for transcript
 Returntype : Bio::EnsEMBL::Transcript
 Exceptions : throws if not passed a Bio::EnsEMBL::Transcript

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


=head2 run

  Arg [1]    : string, working directory
  Description: instantiations Blast runnable and runs it then converts
               blast hits back into genomic coords
  Returntype : none
  Exceptions : none

=cut

sub run{
  my ($self, $dir) = @_;

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

  my $out = $self->clean_output;
  $self->align_hits_to_query($out);
}


=head2 align_hits_to_query

 Arg [1]    : Arrayref of Bio::EnsEMBL::BaseAlignFeatures
 Description: convert the features from blast from peptide coodinates
              to genomic coordinates
 Returntype : None
 Exceptions : None

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
        (-features => $exon_hash{$ex}, -align_type => 'ensembl');
      push(@output, $dna_align_feature);
    }
  }
  $self->output(\@output);
}

1;
