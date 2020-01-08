# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::BlastPep

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::BlastPep->
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
is primarily used to blast the protein sequence of transcripts
againsts a protein database. It instantiates a blast runnable passing it 
the Bio::Seq of the transcript translation as the query sequence. This
module expects all the same arguments as a standard blast with the 
exception or a query sequence as these must be passed into the blast 
runnable it instantiates

Based on BlastTranscriptPep and BlastGenescanPep.

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut


package Bio::EnsEMBL::Analysis::Runnable::BlastPep;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;
use Bio::EnsEMBL::ProteinFeature;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep);

=head2 new

  Arg [1]         : Bio::EnsEMBL::Analysis::Runnable::BlastPep
  Arg [Transcript]: Bio::EnsEMBL::Transcript
  Function        : create a BlastTranscriptPep runnable
  Returntype      : Bio::EnsEMBL::Analysis::Runnable::BlastPep
  Exceptions      : none
  Example         :

=cut

sub new {
  my ($class,@args) = @_;
  my $self = new Bio::EnsEMBL::Analysis::Runnable::Blast @args;
  bless $self,$class;
  return $self;
}

=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastPep
  Arg [2]   : string, working directory
  Function  : instantiations Blast runnable and runs it then converts
  blast hits back into ProteinFeatures
  Returntype: none
  Exceptions: none
  Example   : 

=cut

sub run {
  my ($self, $dir) = @_;

  $self->Bio::EnsEMBL::Analysis::Runnable::Blast::run($dir);
  my $out = $self->output;
  $self->output([], 1);
  $self->align_hits_to_query($out);
}



=head2 align_hits_to_query

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BlastPep
  Arg [2]   : arrayref of Bio::EnsEMBL::BaseAlignFeatures
  Function  : convert the features to ProteinFeatures
  Returntype: arrayref of Bio::EnsEMBL::ProteinFeatures
  Exceptions: 
  Example   : 

=cut

sub align_hits_to_query {
  my ($self, $features)  = @_;

  my @output;
  for my $feature ( @$features ) {

      my $_feature = Bio::EnsEMBL::ProteinFeature->new(
	      				  -start      => $feature->start,
					  -end        => $feature->end, 
					  -hstart     => $feature->hstart,
					  -hend       => $feature->hend,
                                          -percent_id => $feature->percent_id, 
                                          -score      => $feature->score, 
                                          -p_value    => $feature->p_value,
					  -hseqname   => $feature->hseqname, 
                                          -seqname    => $self->query->id,
                                          -analysis   => $feature->analysis);
      push(@output, $_feature);
    }
  $self->output(\@output);
}
1;
