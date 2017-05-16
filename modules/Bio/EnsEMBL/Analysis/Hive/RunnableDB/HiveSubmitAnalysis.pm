# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');


=head2 param_defaults

 Description: It allows the definition of default parameters for all inherting module.
              These are the default values:
               include_non_reference => 0,
               include_duplicates => 0,
               include_lrg => 0, #This one should never be used but it's here for compatibility
               mitochondrion => 0,
               slice_size => 0,
               slice_overlaps => 0,
               coord_system_name => 'toplevel',
               coord_system_version => undef,
               stable_id_prefix => 'ENS',
               region_padding => 10000,
               use_annotation => 0,
               batch_target_size => 1000000,
               table_column_name => 'accession',
               column_names => ['iid'],
 Returntype : Hashref, containing all default parameters
 Exceptions : None

=cut

sub param_defaults {
    my $self = shift;

    return {
        %{$self->SUPER::param_defaults},
        include_non_reference => 0,
        include_duplicates => 0,
        include_lrg => 0, #This one should never be used but it's here for compatibility
        mitochondrion => 0,
        use_annotation => 0,
        batch_slice_ids => 0,
        slice_size => 0,
        slice_overlaps => 0,
        coord_system_name => 'toplevel',
        coord_system_version => undef,
        stable_id_prefix => 'ENS',
        region_padding => 10000,
        batch_target_size => 1000000,
        table_column_name => 'accession',
        column_names => ['iid'],
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

  my $iid_type = $self->param_required('iid_type');
  if($iid_type eq 'chunk') {
    $self->create_chunk_ids();
  } elsif($iid_type eq 'sequence_accession') {
    $self->sequence_accession();
  } elsif($iid_type eq 'rechunk') {
    $self->rechunk_input_ids();
  } else {
      $self->param('target_db', destringify($self->param('target_db'))) if (ref($self->param('target_db')) ne 'HASH');
      my $dba = hrdb_get_dba($self->param('target_db'));
      if($iid_type eq 'slice') {
        $self->create_slice_ids($dba);
      } elsif($iid_type eq 'split_slice') {
        # Not sure it's the correct call but Core has a split_Slice method
        # so it's better to use it
        $self->split_slice($dba);
      } elsif($iid_type eq 'patch_slice') {
        $self->param('include_non_reference',1);
        $self->create_slice_ids($dba);
      } elsif($iid_type eq 'slice_to_feature_ids') {
        $self->convert_slice_to_feature_ids($dba);
      } elsif($iid_type eq 'feature_region') {
        $self->feature_region($dba);
      } elsif($iid_type eq 'feature_id') {
        $self->feature_id($dba);
      } else {
        $self->throw('You have not specified one of the recognised operation types');
      }
  }
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

  if ($self->param('slice_size') < 0) {
    $self->throw('Slice size must be >= 0. Currently '.$self->param('slice_size'));
  }

  my $sa = $dba->get_SliceAdaptor;

  my $slices = $sa->fetch_all($self->param('coord_system_name'),
                              $self->param('coord_system_version'),
                              $self->param('include_non_reference'),
                              $self->param('include_duplicates'),
                              $self->param('include_lrg'));

  if (!$self->param('mitochondrion')) {
    my $mt = $sa->fetch_by_region('toplevel', 'MT');
    if ($mt) {
      my @ids = grep {$_->seq_region_name ne $mt->seq_region_name} @$slices;
      $slices = \@ids;
    }
  }

  if ($self->param('iid_type') eq 'patch_slice') {
    my @pt = ('patch_novel', 'patch_fix');
    my @tmp_slices;
    foreach my $slice (@$slices) {
      foreach my $type (@pt) {
        if (scalar(@{$slice->get_all_Attributes($type)}) > 0) {
          push(@tmp_slices, $slice);
          last;
        }
      }
    }
    $slices = \@tmp_slices;
  }

  if($self->param('slice_size') > 0) {
    $slices = split_Slices($slices, $self->param_required('slice_size'), $self->param('slice_overlaps'));
  }

  if($self->param('feature_constraint')) {
    $slices = $self->filter_slice_on_features($slices, $dba);
  }

  if($self->param('batch_slice_ids')) {
    $slices = $self->batch_slice_ids($slices);
    $self->param('inputlist', $slices);
  }
  else {
    my @input_ids = map {$_->name} @$slices;
    $self->param('inputlist', \@input_ids);
  }

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

  my $output_id_array = [];

  my $slice = $dba->get_SliceAdaptor->fetch_by_name($self->param_required('iid'));
  my $feature_type = $self->param_required('feature_type');
  my $template_sql = 'UPDATE %s SET stable_id = CONCAT("'.$self->param('stable_id_prefix').'%s", LPAD(%s_id, 11, 0)) WHERE seq_region_id = %d';

  if($feature_type eq 'prediction_transcript') {
    if ($self->param_is_defined('create_stable_ids')) {
#        my $sqlquery = 'UPDATE prediction_transcript SET stable_id = CONCAT("'.$self->param('stable_id_prefix').'X", LPAD(prediction_transcript_id, 11, 0)) WHERE seq_region_id = '.$slice->get_seq_region_id;
        my $sqlquery = sprintf($template_sql, $feature_type, 'X', $feature_type, $slice->get_seq_region_id);
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
  } elsif($feature_type eq 'gene') {
    if ($self->param_is_defined('create_stable_ids')) {
#        my $sqlquery = 'UPDATE gene SET stable_id = CONCAT("'.$self->param('stable_id_prefix').'G", LPAD(gene_id, 11, 0)) WHERE seq_region_id = '.$slice->get_seq_region_id;
        my $sqlquery = sprintf($template_sql, $feature_type, 'G', $feature_type, $slice->get_seq_region_id);
        my $sth = $dba->dbc->prepare($sqlquery);
        $sth->execute();
    }
    my $ga = $dba->get_GeneAdaptor;
    my $logic_names = $self->param('logic_name');

    my $use_stable_ids = $self->param_is_defined('use_stable_ids');
    foreach my $logic_name (@$logic_names) {
      my $genes = $ga->fetch_all_by_Slice($slice, $logic_name);
      foreach my $gene_feature (@{$genes}) {
        push(@{$output_id_array}, ($use_stable_ids ? $gene_feature->stable_id : $gene_feature->dbID()));
      }
    }
  } else {
    $self->throw("The feature_type you provided is not currently supported by the code.\nfeature_type: $feature_type");
  }

  $self->param('inputlist', $output_id_array);
}


=head2 split_slice

 Arg [1]    : None
 Description: Split the slice given in 'iid' in smaller chunks given by 'slice_size'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'iid' is not defined
              Throws if 'slice_size' is negative

=cut

sub split_slice {
  my ($self, $dba) = @_;

  if ($self->param('slice_size') < 0) {
    $self->throw('Slice size must be >= 0. Currently '.$self->param('slice_size'));
  }

  my $sa = $dba->get_SliceAdaptor;
  my $slices = split_Slices([$sa->fetch_by_name($self->param_required('iid'))], $self->param('slice_size'), $self->param('slice_overlaps'));
  my @input_ids = map {$_->name} @$slices;

  $self->param('inputlist', \@input_ids);
}


=head2 sequence_accession

 Arg [1]    : None
 Description: Create input_ids based on the sequence accessions stored in Hive
 Returntype : None
 Exceptions : Throws if 'batch_size' is defined but not set
              Throws if 'sequence_table_name' is defined but not set

=cut

sub sequence_accession {
  my ($self) = @_;

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));

  $table_adaptor->column_set();
  my $accessions = $table_adaptor->fetch_all(undef,undef,undef, $self->param('table_column_name'));

  $self->param('inputlist', $self->_chunk_input_ids($self->param_required('batch_size'), $accessions));
}


=head2 rechunk_input_ids

 Arg [1]    : None
 Description: Create smaller chunks of input ids using 'batch_size'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'batch_size' is not defined
              Throws if 'iid' is not defined

=cut

sub rechunk_input_ids {
  my ($self) = @_;

  $self->param('inputlist', $self->_chunk_input_ids($self->param_required('batch_size'), $self->param_required('iid')));
}


=head2 feature_region

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Creates slice input ids with an accession to the evidence supporting the models.
              It is used to realign evidences on a smaller region of the genome so we can use
              more expensive search types.
              If you specify an arrayref of logic names, it will fetch models from these logic
              names. Otherwise it will fetch all genes.
              If you set 'use_annotation' and provide an annotation file with 'annotation_file'
              it will ony realign models which are in the annotation file.
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if 'annotation_file' does not exists or fails to open/close when using 'use_annotation'
              Throws if the accession has a ':' in it's name as it will disrupt the input_id
              Throws if the feature type cannot be processed

=cut

sub feature_region {
  my ($self, $dba) = @_;

  my @input_ids;
  my $feature_type = $self->param_required('feature_type');
  if ($feature_type eq 'gene') {
    my %accessions;
    my $use_annotation = $self->param('use_annotation');
    if ($use_annotation) {
      open(F, $self->param_required('annotation_file')) || $self->throw('Could not open annotation file for reading '.$self->param('annotation_file'));
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
      $slices = $dba->get_SliceAdaptor->fetch_all($self->param('coord_system_name'));
    }
    foreach my $slice (@$slices) {
      my $logic_names = [''];
      if ($self->param_is_defined('logic_name')) {
        $logic_names = $self->param('logic_name');
      }
      my $base_name = join(':', $slice->coord_system_name,
                                $slice->coord_system->coord_system_version,
                                $slice->seq_region_name);
      my %slice_hash;
      foreach my $logic_name (@$logic_names) {
        foreach my $gene (@{$slice->get_all_Genes($logic_name)}) {
          foreach my $transcript (@{$gene->get_all_Transcripts}) {
            my $start = $transcript->seq_region_start-$self->param('region_padding');
            my $end = $transcript->seq_region_end+$self->param('region_padding');
            $start = 1 if ($start < 1);
            $end = $slice->end if ($end > $slice->end);
            my $strand = $use_annotation ? 1 : $transcript->strand;
            my $tbase_name = join(':', $base_name,
                                       $start,
                                       $end,
                                       $strand);
            foreach my $tse (@{$transcript->get_all_supporting_features}) {
              my $hit_name = $tse->hseqname();
              last if ($use_annotation and !exists $accessions{$hit_name});
              if($hit_name =~ /\:/) {
                $self->throw("The hit name for the supporting evidence has a colon in it and this will break the output id structure. ".
                             "Transcript dbID: ".$transcript->dbID.", hit_name: ".$hit_name);
              }
              my $hit_slice_name = $tbase_name.':'.$hit_name;
              $slice_hash{$hit_slice_name} = 1;
            }
          }
        }
      }
      push(@input_ids,keys %slice_hash);
    }
  }
  else {
    $self->throw("The feature_type you provided is not currently supported by the code.\nfeature_type: $feature_type");
  }
  $self->param('inputlist', \@input_ids);
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
  my $type = $self->param_required('feature_type');

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


=head2 filter_slice_on_features

 Arg [1]    : Arrayref of Bio::EnsEMBL::Slice
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor (optional)
 Description: Creates input ids for a region if any object of type 'feature_type'
              can be found for this region.
              If an arrayref of database hash cannot be found in 'feature_dbs', you
              will need to have 'target_db' set.
 Returntype : Arrayref of Bio::EnsEMBL::Slice
 Exceptions : Throws if the 'feature_type' cannot be processed

=cut

sub filter_slice_on_features {
  my ($self,$slices,$dba) = @_;

  my @feature_dbs;
  my @output_slices;

  if($self->param_is_defined('feature_dbs')) {
    foreach my $feature_db (@{$self->param('feature_dbs')}) {
      push(@feature_dbs, hrdb_get_dba($feature_db));
    }
  } else {
    $self->warning("You have selected to use a feature_constraint, but haven't passed in a feature dbs array ref using ".
                   "the feature_dbs param. Defaulting to the target db");
    @feature_dbs = ($dba);
  }

  my $feature_type = $self->param_required('feature_type');
  if($feature_type eq 'gene') {
    foreach my $slice (@$slices) {
      foreach my $db (@feature_dbs) {
        if ($db->get_GeneAdaptor->count_all_by_Slice($slice)) {
          push(@output_slices, $slice);
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


=head2 batch_slice_ids

 Arg [1]    : Arrayref Bio::EnsEMBL::Slice
 Description: Creates batch of slices to cover at least 'batch_target_size' bp.
              It's useful when you have many small scaffolds to avoid having too
              many jobs
 Returntype : Arrayref of arrayref
 Exceptions : None

=cut

sub batch_slice_ids {
  my ($self, $slices) = @_;
  my $batch_target_size = $self->param('batch_target_size');

  my $all_batches = [];
  my $single_batch_array = [];
  my $total_length = 0;
  # Sort shortest first
  foreach my $slice (sort { $a->length <=> $b->length } $slices) {
    my $length = $slice->length;
    if($length + $total_length > $batch_target_size) {
      push(@{$all_batches},[$single_batch_array]);
      $single_batch_array = [];
      $total_length = 0;
    }
    push(@{$single_batch_array}, $slice->name);
    $total_length += $length;
  }
  push(@{$all_batches},[$single_batch_array]);
  return($all_batches);
}


=head2 _chunk_input_ids

 Arg [1]    : Int $batch_size, size of the batch to create
 Arg [2]    : Arrayref Object, an array ref of object
 Description: Split Arg[2] into smaller array of length Arg[1]
 Returntype : Arrayref of Arrayref
 Exceptions : None

=cut

sub _chunk_input_ids {
  my ($self, $batch_size, $values) = @_;
  my @output_id_array;

  my $output_accession_array = [];
  foreach my $accession (@$values) {
    if(scalar(@{$output_accession_array}) == $batch_size) {
      push(@output_id_array, [$output_accession_array]);
      $output_accession_array = [];
    }
    push(@{$output_accession_array},$accession);
  }

  if(scalar(@{$output_accession_array})) {
    push(@output_id_array,[$output_accession_array]);
  }

  return \@output_id_array;
}

1;
