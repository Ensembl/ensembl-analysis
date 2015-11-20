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
use Data::Dumper;

use Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  if($self->param('slice') ||
     $self->param('slice_to_feature_ids') ||
     $self->param('feature_region')) {
     my $dba = $self->hrdb_get_dba($self->param('target_db'));
     $self->hrdb_set_con($dba,'target_db');
  }

  return 1;
}

sub run {
  my $self = shift;

#  if (!($self->param('slice')) && !($self->param('single')) && !($self->param('file')) &&
#      !($self->param('translation_id')) && !($self->param('hap_pair')) && !($self->param('chunk')) &&
#      !($self->param('slice_to_feature_ids')) && !($self->param('split_slice')) && !($self->param('uniprot_accession'))
#     ) {
#    $self->throw("Must define input as either contig, slice, file, translation_id ".
#                 "single, seq_level, top_level, hap_pair, chunk or slice_to_feature_ids");
#  }

  if($self->param('slice') && $self->param('chunk')) {
    $self->throw("You have selected both the slice and the chunk file, select one or the other");
  }

  if($self->param('slice')) {
    $self->create_slice_ids();
  } elsif($self->param('chunk')) {
    $self->create_chunk_ids();
  } elsif($self->param('slice_to_feature_ids')) {
    $self->convert_slice_to_feature_ids();
  } elsif($self->param('split_slice')) {
    $self->split_slice();
  } elsif($self->param('uniprot_accession')) {
    $self->uniprot_accession();
  } elsif($self->param('rechunk_uniprot_accession')) {
    $self->rechunk_uniprot_accession();
  } elsif($self->param('feature_region')) {
    $self->feature_region();
  } elsif($self->param('feature_id')) {
    $self->feature_id();
  } else {
    $self->throw('You have not specified one of the recognised operation types');
  }
  return 1;
}

sub create_slice_ids {
  my ($self) = @_;

  my $input_id_factory = new Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory
     (
       -db => $self->hrdb_get_con('target_db'),
       -slice => $self->param('slice'),
       -single => $self->param('single'),
       -file => $self->param('file'),
       -translation_id => $self->param('translation_id'),
       -seq_level => $self->param('seq_level'),
       -top_level => $self->param('top_level'),
       -include_non_reference => $self->param('include_non_reference'),
       -dir => $self->param('dir'),
       -regex => $self->param('regex'),
       -single_name => 'genome', # Don't know why this is set this way
       -logic_name => $self->param('logic_name'),
       -input_id_type => $self->param('input_id_type'),
       -coord_system => $self->param('coord_system_name'),
       -coord_system_version => $self->param('coord_system_version'),
       -slice_size => $self->param('slice_size'),
       -slice_overlaps => $self->param('slice_overlap'),
       -seq_region_name => $self->param('seq_region_name'),
       -hap_pair => $self->param('hap_pair'),
     );

  $input_id_factory->generate_input_ids;
  $self->output_ids($input_id_factory->input_ids);
}


sub create_chunk_ids {
  my $self = shift;

  if($self->param_is_defined('num_chunk') || $self->param_is_defined('seqs_per_chunk')) {
      $self->make_chunk_files();
  }

  my $input_file;
  my $chunk_dir = $self->param('chunk_output_dir');
  # If this is passed in then the idea is that chunk_output_dir is a parent dir and this dir name
  # is the specific one for the input id
  if($self->param('chunk_dir_name')) {
    $chunk_dir .= "/".$self->param('chunk_dir_name');
  }


  if($self->param_is_defined('input_file_path')) {
      $input_file = $self->param('input_file_path');
    } else {
      $input_file = $self->input_id;
  }

  # Get the name without the extension as fastasplit_random cuts off the extension
  $input_file =~ /[^\/]+$/;
  $input_file = $&;
  $input_file =~ s/\.[^\.]+$//;

  my @chunk_array = glob $chunk_dir."/".$input_file."_chunk_*";

  unless(scalar(@chunk_array)) {
    $self->throw("Found no files in chunk dir using glob. Chunk dir:\n".
                 $chunk_dir."/"."\nChunk generic name:\n".$input_file."_chunk_*");
  }

  for(my $i=0; $i < scalar(@chunk_array); $i++) {
    $chunk_array[$i] =~ /[^\/]+$/;
    $chunk_array[$i] = $&;
    if($self->param('chunk_dir_name')) {
      $chunk_array[$i] = $self->param('chunk_dir_name')."/".$chunk_array[$i];
    }
  }
  $self->output_ids(\@chunk_array);
}


sub make_chunk_files {
  my $self = shift;

  my $input_file;
  my $chunk_dir = $self->param('chunk_output_dir');
  # If this is passed in then the idea is that chunk_output_dir is a parent dir and this dir name
  # is the specific one for the input id
  if($self->param('chunk_dir_name')) {
    $chunk_dir .= "/".$self->param('chunk_dir_name');
  }
  my $chunk_num;

  if($self->param_is_defined('input_file_path')) {
      $input_file = $self->param('input_file_path');
  } elsif($self->param_is_defined('rechunk_dir_path') && $self->param_is_defined('rechunk')) {
    if($self->param('rechunk')) {
      $input_file = $self->param('rechunk_dir_path')."/".$self->input_id;
    }
  }

  else {
      $input_file = $self->input_id;
  }

  unless(-e $input_file) {
      $self->throw("Your input file '".$input_file."' does not exist!!!");
  }

  unless(-e $chunk_dir) {
    `mkdir -p $chunk_dir`;
  }

  unless($self->param_is_defined('fastasplit_random_path')) {
    $self->throw("You haven't defined a path to fastasplit_random. Please define this using the fastasplit_random_path ".
                 " flag in your pipeline config");
  }

  my $fastasplit_random_path = $self->param('fastasplit_random_path');
  unless(-e $fastasplit_random_path) {
    $self->throw("The path provided to the fastasplit_random exe does not exist. Please check the path in the config:\n".
                 $fastasplit_random_path);
  }

  if($self->param_is_defined('seqs_per_chunk')) {
    my $num_seqs = `grep -c '>' $input_file`;
    $chunk_num = int($num_seqs / $self->param('seqs_per_chunk'));
  }

  say "Chunking input file to ".$chunk_num." output files";
  my $fastasplit_command = $fastasplit_random_path." ".$input_file." ".$chunk_num." ".$chunk_dir;
  my $fastasplit_exit_code = system($fastasplit_command);
  unless($fastasplit_exit_code == 0){
    $self->throw($fastasplit_random_path." returned an error code:\n".$fastasplit_exit_code);
  }

}


sub write_output {
  my $self = shift;

  my $output_ids = $self->output_ids();

  unless(scalar(@{$output_ids})) {
    $self->warning("No input ids generated for this analysis!");
  }

  foreach my $output_id (@{$output_ids}) {

    # Skip over mito slices, needs updating/refining
    if($self->param_is_defined('skip_mito') && ($self->param('skip_mito') == 1 || $self->param('skip_mito') eq 'yes') &&
       $self->param_is_defined('slice') && ($self->param('slice') == 1 || $self->param('slice') eq 'yes') &&
       $output_id =~ /^.+\:.+\:MT\:/) {
       next;
    }

    # Yet to implement this, will allow for slices that don't have any of the specified features to be skipped
    if($self->param('check_slice_for_features')) {
      unless($self->check_slice_for_features()) {
        next;
      }
    }

    my $output_hash = {};
    $output_hash->{'iid'} = $output_id;

    $self->dataflow_output_id($output_hash,4);
    $self->dataflow_output_id($output_hash,1);
  }

  return 1;
}

sub convert_slice_to_feature_ids {
  my ($self) = @_;

  unless($self->param('iid')) {
    $self->throw("Failed to provide an input id. Expected to find a slice input id using \$self->param('iid')");
  }

  unless($self->param('feature_type')) {
    $self->throw("You're trying to convert a slice to a set of feature ids but haven't provided a feature type. ".
                 "Expected \$self->param('feature_type')");
  }

  my $output_id_array = [];

  my $input_id = $self->param('iid');
  my $slice = $self->fetch_sequence($input_id,$self->hrdb_get_con('target_db'));
  $self->query($slice);

  if($self->param('feature_type') eq 'prediction_transcript') {
    my $pta = $self->hrdb_get_con('target_db')->get_PredictionTranscriptAdaptor;
    my $logic_names = $self->param('prediction_transcript_logic_names');
    if (!ref($logic_names) || scalar(@$logic_names) == 0 ) {
      $logic_names = ['genscan'];
    }
    my @pts ;
    foreach my $logic_name (@$logic_names) {
      my $pt = $pta->fetch_all_by_Slice($self->query, $logic_name);
      foreach my $pt_feature (@{$pt}) {
        push(@{$output_id_array},$pt_feature->dbID());
      }
    }
  } elsif($self->param('feature_type') eq 'gene') {
    my $ga = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
    my $logic_names = $self->param('gene_logic_names');

    foreach my $logic_name (@$logic_names) {
      my $genes = $ga->fetch_all_by_Slice($self->query, $logic_name);
      foreach my $gene_feature (@{$genes}) {
        push(@{$output_id_array},[$gene_feature->dbID()]);
      }
    }
  } else {
    $self->throw("The feature_type you provided is not currently supported by the code.\nfeature_type: ".$self->param('feature_type'));
  }

  $self->output_ids($output_id_array);
}

sub split_slice {
  my ($self) = @_;

  unless($self->param('iid')) {
    $self->throw("Failed to provide an input id. Expected to find a slice input id using \$self->param('iid')");
  }

  unless($self->param('slice_size')) {
    $self->throw("You selected the split_slice option, but did not specific 'slice_size'. Need a size to split into");
  }

  my $output_id_array = [];

  my $slice_name = $self->param('iid');
  my $target_slice_size = $self->param('slice_size');

  # 'scaffold:PapAnu2.0:JH684492.1:1:489941:1'

  my @slice_array = split(':',$slice_name);
  my $slice_length = $slice_array[4]-$slice_array[3] + 1;
  if($slice_length <= $target_slice_size) {
    push(@{$output_id_array},$slice_name);
  } else {
    my $remainder = $slice_length % $target_slice_size;
    my $loop_count = int($slice_length / $target_slice_size);
    my $slice_start = $slice_array[3];
    my $slice_end = $slice_start + $target_slice_size - 1;
    my $new_slice = $slice_array[0].':'.$slice_array[1].':'.$slice_array[2].':'.$slice_start.':'.$slice_end.':'.$slice_array[5];
    push(@{$output_id_array},$new_slice);
    my $i=0;
    for($i=1; $i<$loop_count; $i++) {
      $slice_start += $target_slice_size;
      $slice_end += $target_slice_size;
      $new_slice = $slice_array[0].':'.$slice_array[1].':'.$slice_array[2].':'.$slice_start.':'.$slice_end.':'.$slice_array[5];
      push(@{$output_id_array},$new_slice);
    }
    if($remainder) {
      $slice_start += $target_slice_size;
      $slice_end += $remainder;
      $new_slice = $slice_array[0].':'.$slice_array[1].':'.$slice_array[2].':'.$slice_start.':'.$slice_end.':'.$slice_array[5];
      push(@{$output_id_array},$new_slice);
    }
  }

  $self->output_ids($output_id_array);
}

sub uniprot_accession {
  my ($self) = @_;

  my $output_id_array = [];

  unless($self->param('uniprot_batch_size')) {
    $self->throw("You've select to batch uniprot ids but haven't passed in a batch size using 'uniprot_batch_size'");
  }

  unless($self->param('uniprot_table_name')) {
    $self->throw("You've select to batch uniprot ids but haven't passed the name of the uniprot table in ".
                 "the pipeline database using 'uniprot_table_name'");
  }

  my $batch_size = $self->param('uniprot_batch_size');
  my $table_name = $self->param('uniprot_table_name');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($table_name);

  $table_adaptor->column_set();
  my $accessions = $table_adaptor->fetch_all(undef,undef,undef,'accession');
  my $accession_array = [];
  foreach my $accession (@{$accessions}) {
    my $size = scalar(@{$accession_array});
    if($size == $batch_size) {
      push(@{$output_id_array},$accession_array);
      $accession_array = [];
    }
    push(@{$accession_array},$accession);
  }

  if(scalar(@{$accession_array})) {
    push(@{$output_id_array},$accession_array);
  }

  $self->output_ids($output_id_array);
}


sub rechunk_uniprot_accession {
  my ($self) = @_;

  my $output_id_array = [];

  unless($self->param('uniprot_batch_size')) {
    $self->throw("You've select to batch uniprot ids but haven't passed in a batch size using 'uniprot_batch_size'");
  }

  unless($self->param('iid')) {
    $self->throw("You've select to rechunk uniprot ids but haven't passed in an input_id using 'iid'");
  }

  my $batch_size = $self->param('uniprot_batch_size');

  my $input_accession_array = $self->param('iid');
  my $output_accession_array = [];
  foreach my $accession (@{$input_accession_array}) {
    my $size = scalar(@{$output_accession_array});
    if($size == $batch_size) {
      push(@{$output_id_array},$output_accession_array);
      $output_accession_array = [];
    }
    push(@{$output_accession_array},$accession);
  }

  if(scalar(@{$output_accession_array})) {
    push(@{$output_id_array},$output_accession_array);
  }

  $self->output_ids($output_id_array);
}



sub feature_region {
  my ($self) = @_;

  unless($self->param('feature_type')) {
    $self->throw("You're trying to convert a slice to a set of feature ids but haven't provided a feature type. ".
                 "Expected \$self->param('feature_type')");
  }

  my $output_id_array = [];
  if($self->param('feature_type') eq 'gene') {
    my $ga = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
    my $logic_names = $self->param('gene_logic_names');
    my $padding = $self->param('region_padding');

    unless($self->param('region_padding')) {
      $self->warning("You didn't pass in any value for padding. Defaulting to 10000");
      $padding = 10000;
    }

    unless($logic_names) {
      $self->throw("You didn't pass in an arrayref of logic names for the genes. Pass this in using the 'gene_logic_names' param");
    }

    foreach my $logic_name (@$logic_names) {
      my $genes = $ga->fetch_all_by_logic_name($logic_name);
      foreach my $gene (@{$genes}) {
        my $transcripts = $gene->get_all_Transcripts();
        foreach my $transcript (@{$transcripts}) {
          my $start = $transcript->seq_region_start;
          my $end = $transcript->seq_region_end;
          my $strand = $transcript->strand;
          my $slice = $transcript->slice();
          my $slice_length = $slice->length();
          if($padding) {
            $start = $start - $padding;
            if($start < 1) {
              $start = 1;
            }
            $end = $end + $padding;
            my $slice = $transcript->slice();
            my $slice_length = $slice->length();
            if($end > $slice_length) {
              $end = $slice_length;
            }
          }

          my @slice_array = split(':',$slice->name());
          $slice_array[3] = $start;
          $slice_array[4] = $end;
          $slice_array[5] = $strand;
          my $new_slice_name = join(':',@slice_array);

          my @transcript_supporting_evidence = @{$transcript->get_all_supporting_features};
          my $slice_hash;
          foreach my $tse (@transcript_supporting_evidence) {
            my $hit_name = $tse->hseqname();
            if($hit_name =~ /\:/) {
              $self->throw("The hit name for the supporting evidence has a colon in it and this will break the output id structure. ".
                           "Transcript dbID: ".$transcript->dbID.", hit_name: ".$hit_name);
            }
            my $hit_slice_name = $new_slice_name.":".$hit_name;
            $slice_hash->{$hit_slice_name} = 1;
          }

          foreach my $output_slice (keys(%$slice_hash)) {
            push(@{$output_id_array},$output_slice);
          }
        }
      }
    }
  } else {
    $self->throw("The feature_type you provided is not currently supported by the code.\nfeature_type: ".$self->param('feature_type'));
  }

  $self->output_ids($output_id_array);
}


sub feature_id {
  my ($self) = @_;

  my $dba = $self->hrdb_get_dba($self->param('target_db'));

  my $output_id_array = [];
  unless($self->param('feature_type')) {
    $self->throw("You're trying to output a set of feature ids but haven't provided a feature type. ".
                 "Expected \$self->param('feature_type')");
  }

  my $type = $self->param('feature_type');
  my $logic_names = $self->param('feature_logic_names');
  my $feature_adaptor;

  if($type eq 'transcript') {
    $feature_adaptor = $dba->get_TranscriptAdaptor;
  } else {
    $self->throw("The feature type you requested is not supported in the code yet. Feature type:\n".$type);
  }

  if($logic_names) {
    foreach my $logic_name (@$logic_names) {
      my $features = $feature_adaptor->fetch_all_by_logic_name($logic_name);
      foreach my $feature (@{$features}) {
        my $db_id = $feature->dbID();
        push(@{$output_id_array},$db_id);
      }
    }
  } else {
    $self->warning("No logic names passed in using 'feature_logic_names' param, so will fetch all features");
    my $features = $feature_adaptor->fetch_all();
    foreach my $feature (@{$features}) {
      my $db_id = $feature->dbID();
      push(@{$output_id_array},$db_id);
    }
  }

  $self->output_ids($output_id_array);
}


sub check_slice_for_features {
  my ($self) = @_;

  my $features_present = 0;
  my $feature_type = $self->param('feature_type');
  return $features_present;

}


sub input_id_factory {
 my ($self,$value) = @_;

  if (defined $value) {
    unless($value->isa('Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory')) {
      $self->throw("To set an input id factory object it must be of type Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory, not a ".$value);
    }
    $self->param('_input_id_factory',$value);
  }

  return self->param('_input_id_factory');
}


sub output_ids {
 my ($self,$value) = @_;

  if (defined $value) {
    $self->param('_output_ids',$value);
  }

  return $self->param('_output_ids');
}

1;
