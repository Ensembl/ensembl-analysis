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

#!/usr/bin/env perl

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  if ( $self->param('slice') || $self->param('slice_to_feature_ids') ) {
    my $dba = $self->hrdb_get_dba( $self->param('target_db') );
    $self->hrdb_set_con( $dba, 'target_db' );
  }

  return 1;
}

sub run {
  my $self = shift;

  if ( !( $self->param('slice') ) &&
       !( $self->param('single') )         &&
       !( $self->param('file') )           &&
       !( $self->param('translation_id') ) &&
       !( $self->param('hap_pair') )       &&
       !( $self->param('chunk') )          &&
       !( $self->param('slice_to_feature_ids') ) )
  {
    $self->throw( "Must define input as either contig, slice, file, translation_id " .
                  "single, seq_level, top_level, hap_pair, chunk or slice_to_feature_ids" );
  }

  if ( $self->param('slice') && $self->param('chunk') ) {
    $self->throw("You have selected both the slice and the chunk file, select one or the other");
  }

  if ( $self->param('slice') ) {
    my $input_id_factory = new Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory(
                                  -db                    => $self->hrdb_get_con('target_db'),
                                  -slice                 => $self->param('slice'),
                                  -single                => $self->param('single'),
                                  -file                  => $self->param('file'),
                                  -translation_id        => $self->param('translation_id'),
                                  -seq_level             => $self->param('seq_level'),
                                  -top_level             => $self->param('top_level'),
                                  -include_non_reference => $self->param('include_non_reference'),
                                  -dir                   => $self->param('dir'),
                                  -regex                 => $self->param('regex'),
                                  -single_name           => 'genome',                                # Don't know why this is set this way
                                  -logic_name            => $self->param('logic_name'),
                                  -input_id_type         => $self->param('input_id_type'),
                                  -coord_system          => $self->param('coord_system_name'),
                                  -coord_system_version  => $self->param('coord_system_version'),
                                  -slice_size            => $self->param('slice_size'),
                                  -slice_overlaps        => $self->param('slice_overlap'),
                                  -seq_region_name       => $self->param('seq_region_name'),
                                  -hap_pair              => $self->param('hap_pair'), );

    $input_id_factory->generate_input_ids;

    $self->output_ids( $input_id_factory->input_ids );
  } ## end if ( $self->param('slice'...))
  elsif ( $self->param('chunk') ) {
    if ( $self->param_is_defined('num_chunk') || $self->param_is_defined('seqs_per_chunk') ) {
      $self->make_chunk_files();
    }
    $self->create_chunk_ids();
  }
  elsif ( $self->param('slice_to_feature_ids') ) {
    $self->convert_slice_to_feature_ids();
  }
  return 1;
} ## end sub run

sub make_chunk_files {
  my $self = shift;

  my $input_file;
  my $chunk_dir = $self->param('chunk_output_dir');
  my $chunk_num;

  if ( $self->param_is_defined('input_file_path') ) {
    $input_file = $self->param('input_file_path');
  }
  elsif ( $self->param_is_defined('rechunk_dir_path') && $self->param_is_defined('rechunk') ) {
    if ( $self->param('rechunk') ) {
      $input_file = $self->param('rechunk_dir_path') . "/" . $self->input_id;
    }
  }

  else {
    $input_file = $self->input_id;
  }

  unless ( -e $input_file ) {
    $self->throw( "Your input file '" . $input_file . "' does not exist!!!" );
  }

  unless ( -e $chunk_dir ) {
    `mkdir -p $chunk_dir`;
  }

  unless ( $self->param_is_defined('fastasplit_random_path') ) {
    $self->throw( "You haven't defined a path to fastasplit_random. Please define this using the fastasplit_random_path " .
                  " flag in your pipeline config" );
  }

  my $fastasplit_random_path = $self->param('fastasplit_random_path');
  unless ( -e $fastasplit_random_path ) {
    $self->throw(
      "The path provided to the fastasplit_random exe does not exist. Please check the path in the config:\n" . $fastasplit_random_path );
  }

  if ( $self->param_is_defined('seqs_per_chunk') ) {
    my $num_seqs = `grep -c '>' $input_file`;
    $chunk_num = int( $num_seqs/$self->param('seqs_per_chunk') );
  }

  say "Chunking input file to " . $chunk_num . " output files";
  my $fastasplit_command   = $fastasplit_random_path . " " . $input_file . " " . $chunk_num . " " . $chunk_dir;
  my $fastasplit_exit_code = system($fastasplit_command);
  unless ( $fastasplit_exit_code == 0 ) {
    $self->throw( $fastasplit_random_path . " returned an error code:\n" . $fastasplit_exit_code );
  }

} ## end sub make_chunk_files

sub create_chunk_ids {
  my $self = shift;

  my $input_file;
  my $chunk_dir = $self->param('chunk_output_dir');

  if ( $self->param_is_defined('input_file_path') ) {
    $input_file = $self->param('input_file_path');
  }
  else {
    $input_file = $self->input_id;
  }

  # Get the name without the extension as fastasplit_random cuts off the extension
  $input_file =~ /[^\/]+$/;
  $input_file = $&;
  $input_file =~ s/\.[^\.]+$//;

  my @chunk_array = glob $chunk_dir . "/" . $input_file . "_chunk_*";

  unless ( scalar(@chunk_array) ) {
    $self->throw(
       "Found no files in chunk dir using glob. Chunk dir:\n" . $chunk_dir . "/" . "\nChunk generic name:\n" . $input_file . "_chunk_*" );
  }

  for ( my $i = 0; $i < scalar(@chunk_array); $i++ ) {
    $chunk_array[$i] =~ /[^\/]+$/;
    $chunk_array[$i] = $&;
  }
  $self->output_ids( \@chunk_array );
} ## end sub create_chunk_ids

sub write_output {
  my $self = shift;

  my $output_ids = $self->output_ids();

  unless ( scalar( @{$output_ids} ) ) {
    $self->warning("No input ids generated for this analysis!");
  }

  foreach my $output_id ( @{$output_ids} ) {
    if ( $self->param_is_defined('skip_mito') &&
         ( $self->param('skip_mito') == 1 || $self->param('skip_mito') eq 'yes' ) &&
         $self->param_is_defined('slice') &&
         ( $self->param('slice') == 1 || $self->param('slice') eq 'yes' ) &&
         $output_id =~ /^.+\:.+\:MT\:/ )
    {
      next;
    }

    my $output_hash = {};
    $output_hash->{'iid'} = $output_id;
    $self->dataflow_output_id( $output_hash, 4 );
    $self->dataflow_output_id( $output_hash, 1 );
  }

  return 1;
} ## end sub write_output

sub convert_slice_to_feature_ids {
  my ($self) = @_;

  unless ( $self->param('iid') ) {
    $self->throw("Failed to provide an input id. Expected to find a slice input id using \$self->param('iid')");
  }

  unless ( $self->param('feature_type') ) {
    $self->throw( "You're trying to convert a slice to a set of feature ids but haven't provided a feature type. " .
                  "Expected \$self->param('feature_type')" );
  }

  my $output_id_array = [];

  my $input_id = $self->param('iid');
  my $slice = $self->fetch_sequence( $input_id, $self->hrdb_get_con('target_db') );
  $self->query($slice);

  if ( $self->param('feature_type') eq 'prediction_transcript' ) {
    my $pta         = $self->hrdb_get_con('target_db')->get_PredictionTranscriptAdaptor;
    my $logic_names = $self->param('prediction_transcript_logic_names');
    if ( !ref($logic_names) || scalar(@$logic_names) == 0 ) {
      $logic_names = ['genscan'];
    }
    my @pts;
    foreach my $logic_name (@$logic_names) {
      my $pt = $pta->fetch_all_by_Slice( $self->query, $logic_name );
      foreach my $pt_feature ( @{$pt} ) {
        push( @{$output_id_array}, $pt_feature->dbID() );
      }
    }
  }
  else {
    $self->throw("The feature_type you provided is not currently supported by the code.\nfeature_type: " . $self->param('feature_type') );
  }

  $self->output_ids($output_id_array);
} ## end sub convert_slice_to_feature_ids

sub input_id_factory {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    unless ( $value->isa('Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory') ) {
      $self->throw(
               "To set an input id factory object it must be of type Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory, not a " . $value );
    }
    $self->param( '_input_id_factory', $value );
  }

  return self->param('_input_id_factory');
}

sub output_ids {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->param( '_output_ids', $value );
  }

  return $self->param('_output_ids');
}

1;
