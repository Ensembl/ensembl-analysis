=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep->
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
 BlastTranscriptPep passing in prediction transcript

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanPep;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlast');

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep
  Function  : fetch sequence and prediction transcripts of database,
  read config files instantiate the filter, parser and finally the blast
  runnables
  Returntype: none
  Exceptions: none
  Example   :

=cut

sub fetch_input{
  my ($self) = @_;

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  $self->hrdb_set_con($dba,'target_db');


  unless($self->param('blast_db_path')) {
    $self->throw("You did not pass in the blast_db_path parameter. This is required to locate the blast db");
  }

 # Make an analysis object and set it, this will allow the module to write to the output db
  $self->create_analysis;
  $self->analysis->program($self->param('blast_program')) if ($self->param_is_defined('blast_program'));
  $self->analysis->program_file($self->param('blast_exe_path')) if ($self->param_is_defined('blast_exe_path'));
  $self->analysis->parameters($self->param('commandline_params')) if ($self->param_is_defined('commandline_params'));
  $self->analysis->db_file($self->param('blast_db_path')) if ($self->param_is_defined('blast_db_path'));
  $self->analysis->db($self->param('blast_db_name')) if ($self->param_is_defined('blast_db_name'));

  my $pta = $dba->get_PredictionTranscriptAdaptor;

  my $input_id = $self->param('iid');
  my $input_id_type = $self->param('iid_type');
  my @pts ;
  if($input_id_type eq 'slice') {
    my $logic_names = $self->param('prediction_transcript_logic_names');
    if ( !ref($logic_names) || scalar(@$logic_names) == 0 ) {
      $logic_names = ['genscan'];
    }

    my $slice = $self->fetch_sequence($input_id,$dba);
    $self->query($slice);

    foreach my $logic_name (@$logic_names) {
      my $pt = $pta->fetch_all_by_Slice($self->query, $logic_name);
      push @pts, @$pt ;
    }
  } elsif($input_id_type eq 'feature_id') {
    my $pt_feature = $pta->fetch_by_dbID($input_id);
    my $pt = [$pt_feature];
    push @pts, @$pt ;

    my $feature_slice = $pt_feature->feature_Slice();
    my $cs = $feature_slice->coord_system;
    my $cs_name = $cs->name;
    my $cs_version = $cs->version;
    my $slice_name = $pt_feature->seq_region_name;
    my $slice_start = $pt_feature->seq_region_start;
    my $slice_end = $pt_feature->seq_region_end;

    my $slice_id = $cs_name.":".$cs_version.":".$slice_name.":".$slice_start.":".$slice_end.":1";

    my $slice = $self->fetch_sequence($slice_id,$dba);
    $self->query($slice);

  } else {
    $self->throw("The input_id type you have specified is not currently supported by this module\ninput_id_type: ".$input_id_type);
  }

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


  foreach my $t(@pts){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep->
      new(
          -transcript => $t,
          -query => $self->query,
          -program => $self->analysis->program_file,
          -parser => $parser,
          -filter => $filter,
          -database => $self->analysis->db_file,
          -analysis => $self->analysis,
          -options => $options_string,
          %blast,
         );
    $self->runnable($runnable);
  }
}

1;


