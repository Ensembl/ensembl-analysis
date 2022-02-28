
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;


=head1 DESCRIPTION

fetches sequence data from database an instantiates and runs the
augustus runnable, this inherits from the Genscan runnableDB an as such doesnt
implement much itself

=head1 CONTACT

Post questions to : anacode-people@sanger.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAugustus;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input{
  my ($self) = @_;

  my $repeat_masking = $self->param('repeat_masking_logic_names');

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  $self->hrdb_set_con($dba,'target_db');

  my $input_id = $self->param('iid');
  my $slice = $self->fetch_sequence($input_id,$dba,$repeat_masking);
  $self->query($slice);


  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                              -program_file => $self->param('augustus_path'),
                                            );

  $self->param('analysis',$analysis);


  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus->new
                 (
                   -query => $self->query,
                   -program => $analysis->program_file,
                   -analysis => $self->param('analysis'),
                   -species => $self->param('species'),
                 );

  $self->runnable($runnable);


  my $seq = $self->query->seq;
  if ($seq =~ /[CATG]{3}/) {
     $self->input_is_void(0);
  } else {
     $self->input_is_void(1);
     $self->warning("Need at least 3 nucleotides - maybe your sequence was fully repeatmasked ...");
  }

  return 1;
}

=head2 runnable_path

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus
  Function  : return the runnable path
  Returntype: string
  Exceptions:
  Example   : my $runnable = $self->runnable_path->new
                               (
                                -query    => $self->query,
                                -program  => $self->analysis->program_file,
                                -analysis => $self->analysis,
                                %parameters,
                               );

=cut


sub runnable_path{
  my ($self);
  return "Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus";
}

sub run {
  my ($self) = @_;

  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    $self->output($runnable->output);
  }

  return 1;
}


sub write_output {
  my ($self) = @_;

  my $pt_adaptor = $self->hrdb_get_con('target_db')->get_PredictionTranscriptAdaptor;
  my $analysis = $self->param('analysis');

  my $ff = $self->feature_factory;
  foreach my $pt(@{ $self->output()}){
    $pt->analysis($self->analysis);
    $pt->slice($self->query) if(!$pt->slice);
    print STDERR "Validate transcript ".$pt->seqname."\n";
    # dismiss transcript with invalid translation
    eval {
      $ff->validate_prediction_transcript($pt, 1);
    };

    if($@) {
      $self->warning($@);
      next;
    }
    $pt_adaptor->store($pt);
  }

}

1;
