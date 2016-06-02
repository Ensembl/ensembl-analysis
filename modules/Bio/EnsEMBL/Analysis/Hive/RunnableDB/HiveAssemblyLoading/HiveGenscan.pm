=head1 LICENSE
# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Genscan - 

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Genscan->
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
genscan runnable


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Runnable::Genscan;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genscan
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   :

=cut



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
                                              -program_file => $self->param('genscan_path'),
                                              -db => 'HumanIso.smat',
                                              -db_file => $self->param('genscan_matrix_path'),
                                            );

  $self->param('analysis',$analysis);

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Genscan->new
                 (
                   -query => $self->query,
                   -program => $analysis->program_file,
                   -analysis => $self->param('analysis'),
                   -matrix => $analysis->db_file,
                   %parameters,
                 );

  $self->runnable($runnable);


  my $seq = $self->query->seq;
  if ($seq =~ /[CATG]{3}/) {
     $self->input_is_void(0);
  } else {
    $self->input_is_void(1);
    $self->warning("Need at least 3 nucleotides - maybe your sequence was fully repeatmasked ...");
    $self->input_job->autoflow(0);
    $self->complete_early('No input seq to process');
  }

  return 1;
}


sub run {
  my ($self) = @_;

  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    $self->output($runnable->output);
  }

  return 1;
}

=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genscan
  Function  : writes the prediction transcripts back to the database
  after validation
  Returntype: none
  Exceptions:
  Example   :

=cut



sub write_output{
  my ($self) = @_;

  my $genscan_adaptor = $self->hrdb_get_con('target_db')->get_PredictionTranscriptAdaptor;
  my $analysis = $self->param('analysis');

  foreach my $feature ( @{ $self->output() } ) {
    $feature->analysis($analysis);

    if ( !defined( $feature->slice() ) ) {
      $feature->slice( $self->query() );
    }

    $self->feature_factory->validate_prediction_transcript($feature, 1);

    eval { $genscan_adaptor->store($feature); };
    if ($@) {
      $self->throw("RunnableDB::write_output() failed:failed to store '".$feature."' into database '".$genscan_adaptor->dbc()->dbname()."': ".$@);
    }
  }

  my $output_hash = {};
  $output_hash->{'iid'} = $self->param('iid');
  $self->dataflow_output_id($output_hash,4);
  $self->dataflow_output_id($output_hash,1);

  return 1;

}


#sub runnable_path{
#  my ($self);
#  return "Bio::EnsEMBL::Analysis::Runnable::Genscan";
#}

1;
