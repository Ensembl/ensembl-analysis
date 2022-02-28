=head1 LICENSE

# Copyright [2017-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscan -

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscan->
  new(
      -analysis => $analysis,
      -db => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

 This module inherits from the Blast runnable and instantiates
 BlastTranscriptPep or BlastTranscriptDNA passing in prediction transcript

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan;

use strict;
use warnings;
use feature 'say';

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlast');

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscan
  Function  : fetch sequence and prediction transcripts of database,
  read config files instantiate the filter, parser and finally the blast
  runnables
  Returntype: none
  Exceptions: none
  Example   :

=cut

sub fetch_input{
  my ($self) = @_;

  unless($self->param('blast_db_path')) {
    $self->throw("You did not pass in the blast_db_path parameter. This is required to locate the blast db");
  }

  my $sequence_type = $self->param('sequence_type');
  my $runnable_class = "";
  if($sequence_type eq 'peptide') {
    use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;
    $runnable_class = "Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep";
  } elsif($sequence_type eq 'dna') {
    use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA;
    $runnable_class = "Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA";
  } else {
    $self->throw("No sequence_type param provided. Must be 'dna' or 'peptide'");
  }

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  my $prediction_transcript_dba = $self->hrdb_get_dba($self->param('prediction_transcript_db'));
  my $pta = $prediction_transcript_dba->get_PredictionTranscriptAdaptor();
  $self->get_adaptor($pta);

  if($self->param_is_defined('dna_db')) {
    say "Attaching dna_db to output db adaptor";
    my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $dba->dnadb($dna_dba);
  } else {
    say "No dna_db param defined, so assuming target_db has dna";
  }

  $self->hrdb_set_con($dba,'target_db');

  $self->create_analysis();
  $self->analysis->program($self->param('blast_program')) if ($self->param_is_defined('blast_program'));
  $self->analysis->program_file($self->param('blast_exe_path')) if ($self->param_is_defined('blast_exe_path'));
  $self->analysis->parameters($self->param('commandline_params')) if ($self->param_is_defined('commandline_params'));
  $self->analysis->db_file($self->param('blast_db_path')) if ($self->param_is_defined('blast_db_path'));
  $self->analysis->db($self->param('blast_db_name')) if ($self->param_is_defined('blast_db_name'));

  my $input_ids = $self->param('iid');
  my $input_id_type = $self->param('iid_type');

  my %blast = %{$self->BLAST_PARAMS};
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }

  # submit blast module to use via analysis_parameters column of analysis table
  my $options_string ;
  my %options = %{$self->PARSER_PARAMS};

  if ( $blast{-type}=~m/ncbi/ ) {
    if ( $options{-query_type}=~m/pep/ ) {
      if ( $options{-database_type}=~m/pep/ ) {
           $options_string = '-p blastp' ;
	 } elsif ( $options{-database_type}=~m/dna/ ) {
         $options_string = '-p tblastn' ;
       }
    }

    if ( $options{-query_type}=~m/dna/ ) {
      if ( $options{-database_type}=~m/dna/ ) {
           $options_string = '-p blastn' ;
	 }elsif ( $options{-database_type}=~m/pep/ ) {
           $options_string = '-p blastx' ;
	 }
    }
  }

  my @prediction_transcripts = ();
  if($input_id_type eq 'slice') {
    foreach my $slice_name (@{$input_ids}) {
      my $slice = $self->fetch_sequence($slice_name,$dba);
      @prediction_transcripts = @{$pta->fetch_all_by_Slice($slice)};

      # Remove anything crossing a right hand slice boundry as it will be picked up on another slice
      if($prediction_transcripts[$#prediction_transcripts]->end > $slice->end) {
        pop(@prediction_transcripts);
      }

      foreach my $prediction_transcript (@prediction_transcripts) {
        $self->create_runnable($prediction_transcript,$slice,$parser,$filter,$options_string,\%blast,$runnable_class);
      }
    }
  } elsif($input_id_type eq 'feature_id') {
    foreach my $feature_id (@{$input_ids}) {
      my $prediction_transcript = $pta->fetch_by_dbID($feature_id);
      my $slice = $prediction_transcript->slice();
      $self->create_runnable($prediction_transcript,$slice,$parser,$filter,$options_string,\%blast,$runnable_class);
    }
  } else {
    $self->throw("The input_id type you have specified is not currently supported by this module\ninput_id_type: ".$input_id_type);
  }

}

sub create_runnable {
  my ($self,$prediction_transcript,$slice,$parser,$filter,$options_string,$blast,$runnable_class) = @_;
  my $runnable = $runnable_class->
                   new(
                        -transcript => $prediction_transcript,
                        -query => $slice,
                        -program => $self->analysis->program_file,
                        -parser => $parser,
                        -filter => $filter,
                        -database => $self->analysis->db_file,
                        -analysis => $self->analysis,
                        -options => $options_string,
                         %{$blast},
                       );
  $runnable->timer($self->param('timer'));
  $self->runnable($runnable);
}


=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscan
  Arg [2]   : Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor
  Function  : get/set prediction transcript adaptor
  Returntype: Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor
  Exceptions: none
  Example   :

=cut

sub get_adaptor {
  my ($self,$pta) = @_;
  if($pta) {
    $self->param('_pta',$pta);
  }

  return($self->param('_pta'));
}

1;
