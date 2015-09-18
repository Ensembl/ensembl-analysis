=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker - 

=head1 SYNOPSIS

  my $repeat_masker = Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $repeat_masker->fetch_input;
  $repeat_masker->run;
  $repeat_masker->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable RepeatMasker which wraps the program RepeatMasker

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the repeat_feature and repeat_consensus tables 

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Runnable::RepeatMasker;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');




=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  $self->hrdb_set_con($dba,'target_db');

  my $input_id = $self->param('iid');
  my $slice = $self->fetch_sequence($input_id,$dba);
  $self->query($slice);

  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                              -program_file => $self->param('repeatmasker_path'),
                                              -parameters => $self->param('commandline_params'),
                                            );

  $self->param('analysis',$analysis);

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::RepeatMasker->new
    (
     -query => $self->query,
     -program => $analysis->program_file,
     -analysis => $self->param('analysis'),
     %parameters,
    );
  $self->runnable($runnable);
  say "Runnable ref: ".ref($self->runnable);
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

sub write_output {
  my ($self) = @_;
  my $rf_adaptor = $self->hrdb_get_con('target_db')->get_RepeatFeatureAdaptor;
  my $analysis = $self->param('analysis');

  # Keep track of the analysis_id we have here, because the store()
  # method might change it if the analysis does not already exist in
  # the output database, which will make running the next job difficult
  # or impossible (because the analysis tables weren't in sync).  If we
  # find that the analysis ID changes, throw() so that the genebuilder
  # realizes that she has to sync thhe analysis tables.


  foreach my $feature ( @{ $self->output() } ) {
    $feature->analysis($analysis);

    if ( !defined( $feature->slice() ) ) {
      $feature->slice( $self->query() );
    }

    $self->feature_factory->validate($feature);

    eval { $rf_adaptor->store($feature); };
    if ($@) {
      $self->throw("RunnableDB::write_output() failed:failed to store '".$feature."' into database '".$rf_adaptor->dbc()->dbname()."': ".$@);
    }
  }

  my $output_hash = {};
  $output_hash->{'iid'} = $self->param('iid');
  $self->dataflow_output_id($output_hash,4);
  $self->dataflow_output_id($output_hash,1);

  return 1;

}

1;
