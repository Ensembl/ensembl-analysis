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

Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker -

=head1 SYNOPSIS

  my $repeat_masker = Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker->
  new(
      -input_id => ['contig::AL805961.22.1.166258:1:5000:1','contig::AL805961.22.1.166258:5001:10000:1'],
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
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(parse_timer);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

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

  # This adaptor MUST be set before attaching the dna_db, as the core api currently
  # forces the dna_db to be the db used for storing the features, even if a separate
  # output db is specified
#  my $rfa = $dba->get_RepeatFeatureAdaptor;
#  $self->get_adaptor($rfa);

  if($self->param('use_genome_flatfile')) {
    say "Ingoring dna table and using fasta file for sequence fetching";
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } elsif($self->param('dna_db')) {
    say "Attaching dna db to target";
    my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $dba->dnadb($dna_dba);
  } else {
    say "Assuming the target db has dna";
  }

#  if($self->param_is_defined('dna_db')) {
#    say "Attaching dna_db to output db adaptor";
#    my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
#    $dba->dnadb($dna_dba);
#  } else {
#    say "No dna_db param defined, so assuming target_db has dna";
#  }

  $self->hrdb_set_con($dba,'target_db');

  my $slice_array = $self->param('iid');
  unless(ref($slice_array) eq "ARRAY") {
    $self->throw("Expected an input id to be an array reference. Type found: ".ref($slice_array));
  }

  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                              -program_file => $self->param('repeatmasker_path'),
                                              -parameters => $self->param('commandline_params'),
                                            );

  $self->analysis($analysis);

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  foreach my $slice_name (@{$slice_array}) {
    my $slice = $self->fetch_sequence($slice_name,$dba);
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::RepeatMasker->new
                   (
                     -query => $slice,
                     -program => $analysis->program_file,
                     -analysis => $analysis,
                      %parameters,
                   );
    $self->runnable($runnable);
  }
  return 1;
}


=head2 run

 Arg [1]    : None
 Description: For each runnable, execute the run method of the runnable and store the output
              in output. If the parameter 'disconnect_jobs' is set to 1, it will disconnect
              the worker from the database. Only use this is your job is longer than 15-30 minutes
 Returntype : Arrayref of object to be stored
 Exceptions : None

=cut

sub run {
  my ($self) = @_;
  $self->dbc->disconnect_if_idle() if ($self->param('disconnect_jobs'));

  # If timer_batch is defined then use this to set the timer for the runnables. For the first
  # runnable the timer will be the value of timer_batch. For the next runnable it will be the
  # time original value of timer_batch minus the time it took to run the first runnable, which
  # is stored in remaining time and so on. The batch will fail whenever the total time spent in
  # runnables reaches the value of timer_batch
  if($self->param_is_defined('timer_batch')) {
    my $remaining_time = parse_timer($self->param('timer_batch'));
    say "You have defined a time limit on your batch: ".$self->param('timer_batch')." (".$remaining_time." secs)";
    foreach my $runnable(@{$self->runnable}) {
      $runnable->timer($remaining_time);
      eval {
        $runnable->run;
        $remaining_time = $runnable->remaining_time;
        $self->output($runnable->output);
        $self->slices_processed($runnable->query->name);
        say "\nCompleted a runnable, remaining time: ".$remaining_time." secs";
      };

      if($@) {
        my $except = $@;
        if($except =~ /still running after your timer/) {
          $self->batch_failed(1);
          $self->param('_branch_to_flow_to_on_fail',-2);
          last;
        }
      }
    } # end foreach
  } else {
    foreach my $runnable (@{$self->runnable}){
      $runnable->run;
      $self->output($runnable->output);
    }
  }
  return $self->output;
}



=head2 write_output

 Arg [1]    : None
 Description: Write the objects stored in output in the output database. It fetches the adaptor
              using the get_adaptor method
              It uses a Bio::EnsEMBL::Analysis::Tools::FeatureFactory to validate the objects
 Returntype : Integer 1
 Exceptions : Throws when storing a feature fails

=cut

sub write_output {
  my ($self) = @_;

  my $adaptor  = $self->get_adaptor();
  my $analysis = $self->analysis();

  # if a batch fails store only the slice names that worked here and output them
  my $slices_processed = $self->slices_processed();
  my $output_ids_to_remove = {};
  my $batch_failed = $self->batch_failed;
  my $failure_branch_code = $self->param('_branch_to_flow_to_on_fail');

  # Store the features on the appropriate slice. For any slice that has
  foreach my $feature ( @{ $self->output() } ) {
    $feature->analysis($analysis);

    if ( !defined( $feature->slice() ) ) {
      $feature->slice( $self->query() );
    }

    $self->feature_factory->validate($feature);

    eval { $adaptor->store($feature); };
    if ($@) {
      $self->throw("RunnableDB::write_output() failed: failed to store '".$feature."' into database '".
                   $adaptor->dbc->dbname."': ".$@);
    }
  }

  # If the batch is marked as failed, loop through all the slices that were in the original input id list
  # and compare them to the list of slices that had features stored. The ones that are not in the stored
  # list are put on the output array and then flow occurs on the failure branch (-2)
  if($batch_failed) {
    my $input_ids = $self->param_required("iid");
    my $output_ids = [];
    foreach my $input_id (@{$input_ids}) {
      unless($slices_processed->{$input_id}) {
        say "Slice ".$input_id." was not processed so will pass on for rebatching";
        push(@{$output_ids},$input_id);
      } else {
        say "Completed store for slice ".$input_id." so will not flow into failed id set for rebatching";
      }
    }

    unless(scalar(@{$output_ids})) {
      $self->throw("The batch was classed as failed, but there were no slices to process in the output array. Something is wrong");
    }

    my $output_hash = {};
    $output_hash->{'iid'} = $output_ids;
    $self->dataflow_output_id($output_hash,$failure_branch_code);
  }

  return 1;
} ## end sub write_output



=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker
  Arg [2]   : Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
  Function  : get/set repeatfeature adaptor
  Returntype: Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
  Exceptions: none
  Example   :

=cut


sub get_adaptor {
#  my ($self,$rfa) = @_;
  my ($self) = @_;
  my $rf_adaptor = $self->hrdb_get_con('target_db')->get_RepeatFeatureAdaptor;
  #if($rfa) {
  #  $self->param('_rfa',$rfa);
  #}

  #return($self->param('_rfa'));
  return($rf_adaptor);
}


sub batch_failed {
  my ($self,$batch_failed) = @_;
  if($batch_failed) {
    $self->param('_batch_failed',$batch_failed);
  }
  return ($self->param('_batch_failed'));
}


sub slices_processed {
  my ($self,$slice_name) = @_;
  unless($self->param_is_defined('_slices_processed')) {
    $self->param('_slices_processed',{});
  }
  if($slice_name) {
    $self->param('_slices_processed')->{$slice_name} = 1;
  }
  return ($self->param('_slices_processed'));
}

1;
