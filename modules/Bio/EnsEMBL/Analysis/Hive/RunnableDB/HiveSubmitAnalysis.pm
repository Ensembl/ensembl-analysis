# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Hive::Utils qw(destringify);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(hrdb_get_dba);
# We will have to choose if we want to really use a InputIDFactory or if this will be the InputIdFactory
use Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');


=head2 param_defaults

 Description: It allows the definition of default parameters for all inherting module.
              These are the default values:
               seq_level => 0,
               top_level => 1,
               include_non_reference => 0,
               hap_pair => 0,
               mitochondrion => 0,
               slice_size => 0,
               slice_overlaps => 0,
               coord_system_name => 'toplevel'
 Returntype : Hashref, containing all default parameters
 Exceptions : None

=cut

sub param_defaults {
    my $self = shift;

    return {
        %{$self->SUPER::param_defaults},
        seq_level => 0,
        top_level => 1,
        include_non_reference => 0,
        hap_pair => 0,
        mitochondrion => 0,
        slice_size => 0,
        slice_overlaps => 0,
        coord_system_name => 'toplevel'
    }
}


=head2 fetch_input

 Arg [1]    : None
 Description: It will create input ids for the next analysis using 'iid_type' and
              other parameters
 Returntype : Integer, 1
 Exceptions : None

=cut

sub fetch_input {
  my $self = shift;

  if($self->param('iid_type') eq 'chunk') {
    $self->create_chunk_ids();
  } elsif($self->param('iid_type') eq 'split_slice') {
    # Maybe we could use the method in Core: Bio::EnsEMBL::Utils::Slice qw(split_Slices);
    $self->split_slice();
  } elsif($self->param('iid_type') eq 'uniprot_accession') {
    $self->uniprot_accession();
  } elsif($self->param('iid_type') eq 'rechunk_uniprot_accession') {
    $self->rechunk_uniprot_accession();
  } elsif($self->param('cdna_accession')) {
    $self->cdna_accession();
  } else {
      $self->param('target_db', destringify($self->param('target_db'))) if (ref($self->param('target_db')) ne 'HASH');
      my $dba = hrdb_get_dba($self->param('target_db'));
      if($self->param('iid_type') eq 'slice') {
        $self->param('slice', 1);
        $self->create_slice_ids($dba);
      } elsif($self->param('iid_type') eq 'patch_slice') {
        $self->param('slice', 1);
        $self->create_patch_ids($dba);
      } elsif($self->param('iid_type') eq 'slice_to_feature_ids') {
        $self->convert_slice_to_feature_ids($dba);
      } elsif($self->param('iid_type') eq 'feature_region') {
        if ($self->param_is_defined('use_annotation')) {
          $self->throw('Could not find annotation file '.$self->param('annotation_file')) if ($self->param('use_annotation') and !-e $self->param('annotation_file'));
          $self->feature_region_annotation($dba);
        }
        else {
          $self->feature_region($dba);
        }
      } elsif($self->param('iid_type') eq 'feature_id') {
        $self->feature_id($dba);
      } else {
        $self->throw('You have not specified one of the recognised operation types');
      }
  }

  $self->param('column_names', ['iid']);
  return 1;
}


=head2 create_slice_ids

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Create input ids based on slices using a Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : None

=cut

sub create_slice_ids {
  my ($self, $dba) = @_;

  my $input_id_factory = new Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory
     (
       -db => $dba,
       -slice => $self->param('slice'),
       -seq_level => $self->param('seq_level'),
       -top_level => $self->param('top_level'),
       -include_non_reference => $self->param('include_non_reference'),
       -coord_system => $self->param('coord_system_name'),
       -coord_system_version => $self->param_is_defined('coord_system_version') ? $self->param('coord_system_version') : undef,
       -slice_size => $self->param('slice_size'),
       -slice_overlaps => $self->param('slice_overlaps'),
#       -seq_region_name => $self->param('seq_region_name'),
       -hap_pair => $self->param('hap_pair'),
       -mt => $self->param_is_defined('mitochondrion') ? $self->param('mitochondrion') : 0,
       -min_slice_length => $self->param('min_slice_length'),
     );

  $input_id_factory->generate_input_ids;
  my $input_ids = $input_id_factory->input_ids;

  if($self->param('feature_constraint')) {
    my $filtered_ids = $self->filter_slice_on_features($input_ids,$dba);
    $input_ids = $filtered_ids;
  }

  if($self->param('batch_slice_ids')) {
    my $batched_ids = $self->batch_slice_ids($input_ids);
    $input_ids = $batched_ids;
  }

  $self->param('inputlist', $input_ids);

}

=head2 create_patch_ids

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Create input ids based on slices for patch_novel and patch_fix assembly exceptions.
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : None

=cut

sub create_patch_ids {
  my ($self,$dba) = @_;
  my $code_patch_novel = 'patch_novel';
  my $code_patch_fix = 'patch_fix';
  my @pt;
  push(@pt,$code_patch_novel,$code_patch_fix);

  my $input_ids = [];

  my $sa = $dba->get_SliceAdaptor();

  #get non-ref but not duplicates
  my @slices = @{$sa->fetch_all('toplevel',undef,1)};
  print scalar(@slices)."\n";
  foreach my $slice (@slices) {
    foreach my $type (@pt) {
      my @slice_attributes = @{$slice->get_all_Attributes($type)};
      if (scalar(@slice_attributes) > 0) {
        push(@{$input_ids},$slice->name());
      }
    }
  }

  if ($self->param('feature_constraint')) {
    my $filtered_ids = $self->filter_slice_on_features($input_ids,$dba);
    $input_ids = $filtered_ids;
  }

  if ($self->param('batch_slice_ids')) {
    my $batched_ids = $self->batch_slice_ids($input_ids);
    $input_ids = $batched_ids;
  }
  $self->param('inputlist',$input_ids);
}

=head2 create_chunk_ids

 Arg [1]    : None
 Description: Creates input ids using 'input_file_path' as input and writes in 'chunk_output_dir'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : None

=cut

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
      $input_file = $self->param('iid');
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
  $self->param('inputlist', \@chunk_array);
}


=head2 make_chunk_files

 Arg [1]    : None
 Description: Create the files when using create_chunk_ids
 Returntype : None
 Exceptions : Throws if 'fastasplit_random_path' is not defined
              Throws if it cannot find 'fastasplit_random_path'
              Throws if 'fastasplit_random_path' failed

=cut

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
      $input_file = $self->param('rechunk_dir_path')."/".$self->param('iid');
    }
  }

  else {
      $input_file = $self->param('iid');
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


=head2 convert_slice_to_feature_ids

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Create input ids using the dbID of feature on the slice provided in 'iid'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'iid' is not defined
              Throws if 'feature_type' is not defined
              Throws if 'feature_type' is not supported

=cut

sub convert_slice_to_feature_ids {
  my ($self, $dba) = @_;

  unless($self->param('iid')) {
    $self->throw("Failed to provide an input id. Expected to find a slice input id using \$self->param('iid')");
  }

  unless($self->param('feature_type')) {
    $self->throw("You're trying to convert a slice to a set of feature ids but haven't provided a feature type. ".
                 "Expected \$self->param('feature_type')");
  }

  my $output_id_array = [];

  my $input_id = $self->param('iid');
  my $slice = $dba->get_SliceAdaptor->fetch_by_name($input_id);

  if($self->param('feature_type') eq 'prediction_transcript') {
    if ($self->param_is_defined('create_stable_ids')) {
        my $sqlquery = 'UPDATE prediction_transcript SET stable_id = CONCAT("'.$self->param('stable_id_prefix').'X", LPAD(prediction_transcript_id, 11, 0)) WHERE seq_region_id = '.$slice->get_seq_region_id;
        my $sth = $dba->dbc->prepare($sqlquery);
        $sth->execute();
    }
    my $pta = $dba->get_PredictionTranscriptAdaptor;
    my $logic_names = $self->param('logic_name');
    if (!ref($logic_names) || scalar(@$logic_names) == 0 ) {
      $logic_names = ['genscan'];
    }
    my @pts ;
    foreach my $logic_name (@$logic_names) {
      my $pt = $pta->fetch_all_by_Slice($slice, $logic_name);
      foreach my $pt_feature (@{$pt}) {
        push(@{$output_id_array},$pt_feature->dbID());
      }
    }
  } elsif($self->param('feature_type') eq 'gene') {
    if ($self->param_is_defined('create_stable_ids')) {
        my $sqlquery = 'UPDATE gene SET stable_id = CONCAT("'.$self->param('stable_id_prefix').'G", LPAD(gene_id, 11, 0)) WHERE seq_region_id = '.$slice->get_seq_region_id;
        my $sth = $dba->dbc->prepare($sqlquery);
        $sth->execute();
    }
    my $ga = $dba->get_GeneAdaptor;
    my $logic_names = $self->param('logic_name');

    foreach my $logic_name (@$logic_names) {
      my $genes = $ga->fetch_all_by_Slice($slice, $logic_name);
      foreach my $gene_feature (@{$genes}) {
        push(@{$output_id_array}, ($self->param_is_defined('use_stable_ids') ? $gene_feature->stable_id : $gene_feature->dbID()));
      }
    }
  } else {
    $self->throw("The feature_type you provided is not currently supported by the code.\nfeature_type: ".$self->param('feature_type'));
  }

  $self->param('inputlist', $output_id_array);
}


=head2 split_slice

 Arg [1]    : None
 Description: Split the slice given in 'iid' in smaller chunks given by 'slice_size'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'iid' is not defined
              Throws if 'slice_size' is not defined

=cut

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

  $self->param('inputlist', $output_id_array);
}


=head2 uniprot_accession

 Arg [1]    : None
 Description: Create input ids based on the custom table in the hive pipeline database specified
              by 'uniprot_table_name'. You need to specify a batch size with 'uniprot_batch_size'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'uniprot_batch_size' is not defined
              Throws if 'uniprot_table_name' is not defined

=cut

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
      push(@{$output_id_array},[$accession_array]);
      $accession_array = [];
    }
    push(@{$accession_array},$accession);
  }

  if(scalar(@{$accession_array})) {
    push(@{$output_id_array},[$accession_array]);
  }

  $self->param('inputlist', $output_id_array);
}


=head2 cdna_accession

 Arg [1]    : None
 Description: Create input_ids based on the sequence accessions stored in Hive
 Returntype : None
 Exceptions : Throws if 'cdna_batch_size' is defined but not set
              Throws if 'cdna_table_name' is defined but not set

=cut

sub cdna_accession {
  my ($self) = @_;

  my $output_id_array = [];

  unless($self->param('cdna_batch_size')) {
    $self->throw("You've selected to batch cdna ids but haven't passed in a batch size using 'cdna_batch_size'");
  }

  unless($self->param('cdna_table_name')) {
    $self->throw("You've selected to batch cdna ids but haven't passed the name of the cdna table in ".
                 "the pipeline database using 'cdna_table_name'");
  }

  my $batch_size = $self->param('cdna_batch_size');
  my $table_name = $self->param('cdna_table_name');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($table_name);

  $table_adaptor->column_set();
  my $accessions = $table_adaptor->fetch_all(undef,undef,undef,'accession');
  my $accession_array = [];
  foreach my $accession (@{$accessions}) {
    my $size = scalar(@{$accession_array});
    if($size == $batch_size) {
      push(@{$output_id_array},[$accession_array]);
      $accession_array = [];
    }
    push(@{$accession_array},$accession);
  }

  if(scalar(@{$accession_array})) {
    push(@{$output_id_array},[$accession_array]);
  }

  $self->param('inputlist', $output_id_array);
}

=head2 rechunk_uniprot_accession

 Arg [1]    : None
 Description: Create smaller chunks of uniprot accession using 'uniprot_batch_size'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'uniprot_batch_size' is not defined
              Throws if 'iid' is not defined

=cut

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
      push(@{$output_id_array},[$output_accession_array]);
      $output_accession_array = [];
    }
    push(@{$output_accession_array},$accession);
  }

  if(scalar(@{$output_accession_array})) {
    push(@{$output_id_array},[$output_accession_array]);
  }

  $self->param('inputlist', $output_id_array);
}

sub feature_region_annotation {
  my ($self, $dba) = @_;

  my %accessions;
  my @input_ids;
  my $use_annotation = $self->param_is_defined('annotation_file');
  if ($use_annotation) {
    open(F, $self->param('annotation_file')) or $self->throw('Could not open supplied annotation file for reading');
    while (<F>) {
      my @fields = split;
      $accessions{$fields[0]} = 1;
    }
    close(F) || $self->throw('Could not close annotation_file: '.$self->param('annotation_file'));
  }

  my $slices;
  if ($self->param_is_defined('iid')) {
    $slices = [$dba->get_SliceAdaptor->fetch_by_name($self->param('iid'))];
  }
  else {
    $slices = $dba->get_SliceAdaptor->fetch_all('toplevel');
  }
  foreach my $slice (@$slices) {
    foreach my $gene (@{$slice->get_all_Genes}) {
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
        my $evidence = $transcript->get_all_supporting_features->[0];
        last if ($use_annotation and !exists $accessions{$evidence->hseqname});
        my $slice = $transcript->slice;
        my $start = $transcript->seq_region_start-$self->param('region_padding');
        my $end = $transcript->seq_region_end+$self->param('region_padding');
        if ($self->param_is_defined('region_padding')) {
          $start = 1 if ($start < 1);
          $end = $slice->end if ($end > $slice->end);
        }
        push(@input_ids, join(':', $transcript->coord_system_name,
                                $slice->coord_system->version,
                                $slice->seq_region_name,
                                $start,
                                $end,
                                1,
                                ':'.$evidence->hseqname));
      }
    }
  }
  $self->param('inputlist', \@input_ids);
}

=head2 feature_region

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Creates input ids based on the presence of feature 'feature_type' in a region
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'feature_type' is not defined
              Throws if 'logic_name' is not defined
              Throws if the name of the supporting evidence contains ':'

=cut

sub feature_region {
  my ($self, $dba) = @_;

  unless($self->param('feature_type')) {
    $self->throw("You're trying to convert a slice to a set of feature ids but haven't provided a feature type. ".
                 "Expected \$self->param('feature_type')");
  }

  my $output_id_array = [];
  if($self->param('feature_type') eq 'gene') {
    my $ga = $dba->get_GeneAdaptor;
    my $logic_names = $self->param('logic_name');
    my $padding = $self->param('region_padding');

    unless($self->param('region_padding')) {
      $self->warning("You didn't pass in any value for padding. Defaulting to 10000");
      $padding = 10000;
    }

    unless($logic_names) {
      $self->throw("You didn't pass in an arrayref of logic names for the genes. Pass this in using the 'logic_name' param");
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

  $self->param('inputlist', $output_id_array);
}


=head2 feature_id

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Creates input ids based on features 'feature_type'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'feature_type' is not defined
              Throws if 'feature_type' is not supported

=cut

sub feature_id {
  my ($self, $dba) = @_;

  my $output_id_array = [];
  unless($self->param('feature_type')) {
    $self->throw("You're trying to output a set of feature ids but haven't provided a feature type. ".
                 "Expected \$self->param('feature_type')");
  }

  my $type = $self->param('feature_type');

#  my $logic_names = $self->param('logic_name');
  my $logic_names = $self->param('feature_logic_names');
  # feature_restriction is a way to do specific restrictions that aren't easy to model. One example is 'protein_coding'
  # which will check if the feature hash a translation
  my $feature_restriction = $self->param('feature_restriction');
  my $feature_adaptor;

  if($type eq 'transcript') {
    $feature_adaptor = $dba->get_TranscriptAdaptor;
  } elsif($type eq 'gene') {
    $feature_adaptor = $dba->get_GeneAdaptor;
  } else {
    $self->throw("The feature type you requested is not supported in the code yet. Feature type:\n".$type);
  }

  my $features= [];
  if($logic_names) {
    foreach my $logic_name (@$logic_names) {
      push(@{$features},@{$feature_adaptor->fetch_all_by_logic_name($logic_name)});
    }
  } else {
    $self->warning("No logic names passed in using 'feature_logic_names' param, so will fetch all features");
    push(@{$features},@{$feature_adaptor->fetch_all()});
  }

  foreach my $feature (@{$features}) {
    unless($self->feature_restriction($feature,$type,$feature_restriction)) {
      my $db_id = $feature->dbID();
      push(@{$output_id_array},$db_id);
    }
  }

  $self->param('inputlist', $output_id_array);
}


sub feature_restriction {
  my ($self,$feature,$type,$restriction) = @_;
  my $feature_restricted = 0;

  if($restriction) {
    # Transcript restrictions go here
    if($type eq 'gene' || $type eq 'transcript') {
      if($restriction eq 'protein_coding') {
        # Note initially this is based on translation and not biotype, but for projection I've switched to to biotype temporarily
        unless($feature->biotype() eq 'protein_coding') {
          $feature_restricted = 1;
        }
      } elsif($restriction eq 'biotype') {
        # future code with a new param for a biotype array should go here
      }
    } # End if type eq gene or transcript
  } # End if restriction

  return($feature_restricted);

}


sub filter_slice_on_features {
  my ($self,$slice_names,$dba) = @_;

  my $feature_dbs = [];
  my @output_slices;
  my $feature_slices = {};

  unless($self->param('feature_type')) {
    $self->throw("You have selected to use a feature_constraint, but haven't passed in the feature_type parameter");
  }

  unless($self->param('feature_dbs')) {
    $self->warning("You have selected to use a feature_constraint, but haven't passed in a feature dbs array ref using ".
                   "the feature_dbs param. Defaulting to the target db");
   push(@{$feature_dbs},$dba);
  } else {
    foreach my $feature_db (@{$self->param('feature_dbs')}) {
      $dba = hrdb_get_dba($feature_db);
      push(@{$feature_dbs},$dba);
    }
  }

  my $feature_type = $self->param('feature_type');
  if($feature_type eq 'gene') {
    foreach my $slice_name (@$slice_names) {
      foreach my $db (@$feature_dbs) {
        if ($db->get_GeneAdaptor->count_all_by_Slice($db->get_SliceAdaptor->fetch_by_name($slice_name))) {
          push(@output_slices, $slice_name);
          last;
        }
      }
    }
  } else {
    $self->throw("The feature type you have selected to constrain the slices on is not currently implemented in the code. ".
                 "Feature type selected: ".$feature_type);
  }
  return \@output_slices;
}

sub batch_slice_ids {
  my ($self,$slice_names) = @_;
  my $batch_target_size = 1000000;
  if($self->param('batch_target_size')) {
    $batch_target_size = $self->param('batch_target_size');
  }

  # "chromosome:Mmul_8.0.1:1:1:998163:1"
  my %length_hash;
  foreach my $slice_name (@{$slice_names}) {
    $slice_name =~ /(\d+)\:(\d+)\:\d+$/;
    my $start = $1;
    my $end = $2;
    my $length = $end - $start + 1;
    $length_hash{$slice_name} = $length;
  }

  my $all_batches = [];
  my $single_batch_array = [];
  my $total_length = 0;
  # Sort shortest first
  foreach my $slice_name (sort { $length_hash{$a} <=> $length_hash{$b} } keys %length_hash) {
    my $length = $length_hash{$slice_name};
    if($length + $total_length > $batch_target_size) {
      push(@{$all_batches},[$single_batch_array]);
      $single_batch_array = [];
      $total_length = 0;
      push(@{$single_batch_array},$slice_name);
      $total_length += $length;
    } else {
      push(@{$single_batch_array},$slice_name);
      $total_length += $length;
    }
  }
  push(@{$all_batches},[$single_batch_array]);
  return($all_batches);
}


sub input_id_factory {
 my ($self,$value) = @_;

  if (defined $value) {
    unless($value->isa('Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory')) {
      $self->throw("To set an input id factory object it must be of type Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory, not a ".$value);
    }
    $self->param('_input_id_factory',$value);
  }

  return $self->param('_input_id_factory');
}

1;
