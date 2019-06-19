# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Hive::Utils qw(destringify);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(hrdb_get_dba is_slice_name);
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');


=head2 param_defaults

 Description: It allows the definition of default parameters for all inherting module.
              These are the default values:
               include_non_reference => 0,
               feature_id_include_non_reference => 1,
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
        feature_id_include_non_reference => 1,
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
        logic_name => [],
        column_names => ['iid'],
        feature_restriction => undef,
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

  if($self->param('skip_analysis')) {
    $self->input_job->autoflow(0);
    $self->complete_early('Skip analysis flag is enabled, so no input ids will be generated');
  }

  my $iid_type = $self->param_required('iid_type');
  if($iid_type eq 'chunk') {
    $self->create_chunk_ids();
  } elsif($iid_type eq 'sequence_accession') {
    $self->sequence_accession();
  } elsif($iid_type eq 'rechunk') {
    $self->rechunk_input_ids();
  } elsif($iid_type eq 'fastq_range') {
    $self->fastq_range($self->param_required('fastq_file'),$self->param_required('batch_size'));
  } else {
      $self->param('target_db', destringify($self->param('target_db'))) if (ref($self->param('target_db')) ne 'HASH');
      my $dba = hrdb_get_dba($self->param('target_db'));
      if($iid_type eq 'slice') {
        $self->create_slice_ids($dba);
      } elsif($iid_type eq 'split_slice') {
        # Not sure it's the correct call but Core has a split_Slice method
        # so it's better to use it
        $self->split_slice($dba);
      } elsif($iid_type eq 'stranded_slice') {
        $self->create_stranded_slice_ids($dba);
      } elsif($iid_type eq 'rebatch_and_resize_slices') {
        $self->rebatch_and_resize_slices($dba);
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
        $self->throw('You have not specified one of the recognised operation types: '.$iid_type);
      }
  }
}


=head2 create_slice_ids

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Create input ids based on slices
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
my @input = map {$_->name} @$slices;
say "slices from db are ", scalar(@input);
  if (!$self->param('mitochondrion')) {
    my $mt = $sa->fetch_by_region('toplevel', 'MT');
    if ($mt) {
say "mt is ", $mt;
      my @ids = grep {$_->seq_region_name ne $mt->seq_region_name} @$slices;
      $slices = \@ids;
    }
say "inside MT";
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
say "slice > 0";
  }

  if($self->param('min_slice_length')) {
    $slices = $self->filter_slice_on_size($slices);
say "slice length is ", $self->param('min_slice_length');
  }

  if($self->param('feature_constraint')) {
    $slices = $self->filter_slice_on_features($slices, $dba);
say "constraint is ", $self->param('feature_constraint');
  }

  if($self->param('batch_slice_ids')) {
my @inpt = map {$_->name} @$slices;
say "slices from db are ", scalar(@inpt);
    $slices = $self->batch_slice_ids($slices);
#my @inp = map {$_->name} @$slices;
#say "slices from db are ", scalar(@inp);
    say "batch is ", Dumper($slices);
    $self->param('inputlist', $slices);
  } else {
    my @input_ids = map {$_->name} @$slices;
   say "slices from db are ", Dumper(@input_ids);
    $self->param('inputlist', \@input_ids);
  }

}


=head2 create_stranded_slice_ids

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: Create input ids based on slices, creates two sets, one for each strand
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : None

=cut

sub create_stranded_slice_ids {
  my ($self,$dba) = @_;

  $self->create_slice_ids($dba);
  my @slice_names = @{$self->param('inputlist')};
  my @stranded_slice_names = ();
  foreach my $slice_name (@slice_names) {
    push(@stranded_slice_names,$slice_name);
    my $stranded_slice_name = $slice_name;
    $stranded_slice_name =~ s/\:[^:]+$/\:-1/;
    push(@stranded_slice_names,$stranded_slice_name);
  }
  $self->param('inputlist',\@stranded_slice_names);
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

  my $feature_type = $self->param_required('feature_type');
  my $template_sql = 'UPDATE %s SET stable_id = CONCAT("'.$self->param('stable_id_prefix').'%s", LPAD(%s_id, 11, 0)) WHERE seq_region_id = %d';

  my $adaptor;
  my $feature_denominator;
  my $id_method;
  if ($self->param_is_defined('use_stable_ids') and $self->param('use_stable_ids')) {
    $id_method = 'stable_id';
  }
  else {
    $id_method = 'dbID';
  }
  my $batch_size = 0;
  my $logic_names = $self->param('logic_name');
  if($feature_type eq 'prediction_transcript') {
    $adaptor = $dba->get_PredictionTranscriptAdaptor;
    $feature_denominator = 'X';
    $batch_size = $self->param_required('batch_size');
    if (!ref($logic_names) || scalar(@$logic_names) == 0 ) {
      $logic_names = ['genscan'];
    }
  } elsif($feature_type eq 'gene') {
    $adaptor = $dba->get_GeneAdaptor;
    $feature_denominator = 'G';
    if (!ref($logic_names) || scalar(@$logic_names) == 0 ) {
      $logic_names = [undef];
    }
  } else {
    $self->throw("The feature_type you provided is not currently supported by the code.\nfeature_type: $feature_type");
  }
  my $slice_names = $self->param_required('iid');
  if (!ref($slice_names)) {
    $slice_names = [$slice_names];
  }
  my $sa = $dba->get_SliceAdaptor();
  foreach my $slice_name (@{$slice_names}) {
    my $slice = $sa->fetch_by_name($slice_name);
    if ($self->param_is_defined('create_stable_ids') and $self->param('create_stable_ids')) {
      my $sqlquery = sprintf($template_sql, $feature_type, $feature_denominator, $feature_type, $slice->get_seq_region_id);
      my $sth = $dba->dbc->prepare($sqlquery);
      $sth->execute();
    }

    foreach my $logic_name (@$logic_names) {
      foreach my $feature (@{$adaptor->fetch_all_by_Slice($slice, $logic_name)}) {
        push(@{$output_id_array}, $feature->$id_method);
      }
    }
  }
  if ($batch_size) {
    $self->param('inputlist', $self->_chunk_input_ids($batch_size, $output_id_array));
  }
  else {
    $self->param('inputlist', $output_id_array);
  }
}


=head2 batch_feature_ids

 Arg [1]    : None
 Description: Batch a set of input ids
 Returntype : None
 Exceptions : Throws if 'iid' is not defined
              Throws if 'slice_size' is negative

=cut


=head2 split_slice

 Arg [1]    : None
 Description: Split the slice given in 'iid' in smaller chunks given by 'slice_size'
              If 'iid' is an arrayref, it considers that it is small enough and does
              not split the slices.
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

  my $iid = $self->param_required('iid');
  if (ref($iid) eq 'ARRAY' and @$iid > 1) {
    $self->param('inputlist', [[$iid]]);
  }
  else {
    my $sa = $dba->get_SliceAdaptor;
    if (ref($iid) eq 'ARRAY') {
      $iid = $iid->[0];
    }
    my $slices = split_Slices([$sa->fetch_by_name($iid)], $self->param('slice_size'), $self->param('slice_overlaps'));
    my @input_ids = map {$_->name} @$slices;

    $self->param('inputlist', \@input_ids);
  }
}


=head2 rebatch_and_resize_slices

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Description: This can be used to rebatch a list of slices towards a new target
              batch size. So you could rebatch an arrayref of slice names from a
              1mb target size to a 100kb target size. This can be useful if an
              assembly has a huge number of tiny toplevel sequences, sometimes
              the overhead associated with having many runnables is large. You
              can also ask to split the individual slices themselves, so if for example
              you had 1mb slices you could ask for 100kb slice sizes instead. This
              dual functionality should cover decreasing slices in all scenarios.
              You can do things like set target batch size to 100kb and the slice
              size to 100kb so that if you had a 1mb batch to begin with, consisting
              of a 400kb slice and a 600kb slice you would end up with 10 output
              jobs, each with a 100kb slice. You can do split the input batch into
              individual jobs by setting the target batch size to 1, but then having
              slize size set to whatever you want in terms of splitting stuff. This
              would work well for analyses that get stuck on certain regions
 Returntype : None
 Exceptions : Throws if 'iid' is not defined
              Throws if 'slice_size' defined and it is negative or zero

=cut

sub rebatch_and_resize_slices {
  my ($self, $dba) = @_;

  my $sa = $dba->get_SliceAdaptor;

  my $slices = [];
  foreach my $slice_name (@{$self->param_required('iid')}) {
    my $slice = $sa->fetch_by_name($slice_name);
    push(@{$slices},$slice);
  }

  if ($self->param_is_defined('slice_size')) {
    unless($self->param('slice_size') > 0) {
      $self->throw("Slice size must be > 0. Currently ".$self->param('slice_size'));
    }
    $slices = split_Slices($slices, $self->param('slice_size'), $self->param('slice_overlaps'));
  }

  $self->param('inputlist',$self->batch_slice_ids($slices));
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
  my $constraint = undef;
  if ($self->param_is_defined('constraint')) {
    $constraint = $self->param('constraint');
  }
  my $accessions = $table_adaptor->fetch_all($constraint, undef, undef, $self->param('table_column_name'));

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
  my $slices;
  if ($self->param_is_defined('iid')) {
    if (is_slice_name($self->param('iid'))) {
      $slices = [$dba->get_SliceAdaptor->fetch_by_name($self->param('iid'))];
    }
  }
  else {
    $slices = $dba->get_SliceAdaptor->fetch_all($self->param('coord_system_name'));
  }
  my $logic_names = $self->param('logic_name');
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

    foreach my $slice (@$slices) {
      my $base_name = join(':', $slice->coord_system_name,
                                $slice->coord_system->version,
                                $slice->seq_region_name);
      my %slice_hash;
      foreach my $logic_name (@$logic_names) {
        foreach my $gene (@{$slice->get_all_Genes($logic_name)}) {
          foreach my $transcript (@{$gene->get_all_Transcripts}) {
            my ($start, $end) = _padded_slice_coord($transcript, $slice, $self->param('region_padding'));
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
  elsif ($feature_type eq 'protein_align_feature') {
    if ($slices) {
      foreach my $slice (@$slices) {
        my $base_name = join(':', $slice->coord_system_name,
                                  $slice->coord_system->version,
                                  $slice->seq_region_name);
        my %slice_hash;
        foreach my $logic_name (@$logic_names) {
          foreach my $paf (@{$slice->get_all_ProteinAlignFeatures($logic_name)}) {
            my ($start, $end) = _padded_slice_coord($paf, $slice, $self->param('region_padding'));
            my $tbase_name = join(':', $base_name,
                                       $start,
                                       $end,
                                       1,
                                       $logic_name);
            my $hit_name = $paf->hseqname();
            if($hit_name =~ /\:/) {
              $self->throw("The supporting evidence has a colon in its name and this will break the output id structure. ".
                           "dbID: ".$paf->dbID.", name: $hit_name");
            }
            my $hit_slice_name = $tbase_name.':'.$hit_name;
            $slice_hash{$hit_slice_name} = 1;
          }
        }
        push(@input_ids,keys %slice_hash);
      }
    }
    else {
# this code is only for BlastMiniGenewise when recovering from genblast
      my $paf_adaptor = $dba->get_ProteinAlignFeatureAdaptor;
      my %slice_hash;
      foreach my $iid (@{$self->param('iid')}) {
        if($iid =~ /\:/) {
          $self->throw('The supporting evidence has a colon in its name and this will break the output id structure. '.
                       "name: $iid");
        }
        foreach my $logic_name (@$logic_names) {
          foreach my $paf (@{$paf_adaptor->fetch_all_by_hit_name($iid, $logic_name)}) {
            my $slice = $paf->slice;
            my $hit_slice_name = join(':', $slice->coord_system_name,
                                      $slice->coord_system->version,
                                      $slice->seq_region_name,
                                      $slice->start,
                                      $slice->end,
                                      1,
                                      $logic_name,
                                      $iid);
            $slice_hash{$hit_slice_name} = 1;
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

  my $logic_names = [undef];
  if ($self->param_is_defined('feature_logic_names')) {
    $logic_names = $self->param('feature_logic_names');
  }
  else {
    $self->warning("No logic names passed in using 'feature_logic_names' param, so will fetch all features");
  }
  # feature_restriction is a way to do specific restrictions that aren't easy to model. One example is 'protein_coding'
  # which will check if the feature hash a translation
  my $feature_restriction = $self->param_is_defined('feature_restriction') ? $self->param('feature_restriction') : undef;
  my $feature_adaptor;

  if($type eq 'transcript') {
    $feature_adaptor = $dba->get_TranscriptAdaptor;
  } elsif($type eq 'gene') {
    $feature_adaptor = $dba->get_GeneAdaptor;
  } else {
    $self->throw("The feature type you requested is not supported in the code yet. Feature type:\n".$type);
  }
  my $slices;
  if ($self->param_is_defined('iid') and is_slice_name($self->param('iid'))) {
    $slices = [$dba->get_SliceAdaptor->fetch_by_name($self->param('iid'))];
  }
  else {
    $slices = $dba->get_SliceAdaptor->fetch_all($self->param('coord_system_name'));
  }
  foreach my $slice (@$slices) {
    foreach my $logic_name (@$logic_names) {
        foreach my $feature (@{$feature_adaptor->fetch_all_by_Slice($slice, $logic_name)}) {
          if($self->param_is_defined('exclude_biotype')) {
            foreach my $biotype (@{$self->param('exclude_biotype')}){
               if ($feature->biotype eq $biotype) {
                   $self->warning("You've defined a biotype that is not allowed to be copied. Something is wrong");
               }
               else{
                  push(@$output_id_array, $feature->dbID) unless ($self->feature_restriction($feature, $type, $feature_restriction));
               }
           }
         }
         else{
                  push(@$output_id_array, $feature->dbID) unless ($self->feature_restriction($feature, $type, $feature_restriction));
               }
      }
    }
  }

  if($self->param_is_defined('batch_size')) {
    unless($self->param('batch_size') > 0) {
      $self->throw("You've defined a batch size param of less than 1. Something is wrong");
    }
    $output_id_array = $self->_chunk_input_ids($self->param('batch_size'),$output_id_array);
  }

  $self->param('inputlist', $output_id_array);
}


sub feature_restriction {
  my ($self,$feature,$type,$restriction) = @_;
  my $feature_restricted = 0;

  if($restriction) {
    # Transcript restrictions go here
    if($type eq 'gene' || $type eq 'transcript') {
      if($restriction eq 'has_translation') {
        # Note initially this is based on translation and not biotype, but for projection I've switched to to biotype temporarily
        unless($feature->translation) {
          $feature_restricted = 1;
          return($feature_restricted);
        }
      } elsif($restriction eq 'biotype') {
        unless($self->param_required('biotypes')->{$feature->biotype()}) {
          $feature_restricted = 1;
        }
      } elsif($restriction eq 'projection') {
        return($self->assess_projection_transcript($feature));
      } else {
        $self->throw("You've selected a features restriction type that is not recognised: ".$restriction);
      }
    } # End if type eq gene or transcript
  } # End if restriction

  return($feature_restricted);

}



sub assess_projection_transcript {
  my ($self,$current_transcript) = @_;

  unless(ref($current_transcript) eq "Bio::EnsEMBL::Transcript") {
    $self->throw("The assess_projection_transcript subroutine expects a transcript object. Found object type: ".ref($current_transcript));
  }

  my $biotypes = $self->param_required('biotypes');
  my $feature_restricted = 0;
  unless($biotypes->{$current_transcript->biotype()}) {
    $feature_restricted = 1;
    return($feature_restricted);
  }

  my $attribs = $current_transcript->get_all_Attributes();
  my $readthrough = 0;
  my $cds_incomplete = 0;
  foreach my $attrib (@{$attribs}) {
    my $code = $attrib->code();
    if($code eq 'cds_start_NF' || $code eq 'cds_end_NF') {
      $cds_incomplete = 1;
    }

    # Remove readthroughs
    if($code eq 'readthrough_tra') {
      $feature_restricted = 1;
      return($feature_restricted);
    }
  }

  if($cds_incomplete) {
    my $parent_gene = $current_transcript->get_Gene();
    my $all_transcripts = $parent_gene->get_all_Transcripts();
    # If it's the only transcript, then keep it, so return 0
    if(scalar(@{$all_transcripts}) == 1) {
      return($feature_restricted);
    }

    my $current_cds_length = $current_transcript->translation->length();

    # Ignore tiny transcripts
    my $min_cds_length = 30;
    if($current_cds_length < $min_cds_length) {
      $feature_restricted = 1;
      return($feature_restricted);
    }

    # If any other coding transcript has a larger CDS then restrict this transcript
    foreach my $transcript (@{$all_transcripts}) {
      if($transcript->translation) {
        my $cds_length = $transcript->translation->length();
        if($cds_length > $current_cds_length) {
          $feature_restricted = 1;
          return($feature_restricted);
        }
      }
    }
  } # end if($cds_incomplete)

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


=head2 filter_slice_on_size

 Arg [1]    : Arrayref of Bio::EnsEMBL::Slice
 Description: Filters the slices based on the the min_slice_length param. Useful if an assembly contains huge amounts of tiny
              toplevel slices. While these small slices can be batched together, the overhead can sometimes be extremely large
              for little to no benefit in terms of annotation
 Returntype : Arrayref of Bio::EnsEMBL::Slice
 Exceptions : throw if there are no slices after filtering

=cut

sub filter_slice_on_size {
  my ($self,$slices) = @_;

  my @output_slices = ();
  my $min_size = $self->param('min_slice_length');
  say "Slice count pre length filter: ".scalar(@{$slices});
  foreach my $slice (@$slices) {
    if($slice->length >= $min_size) {
      push(@output_slices, $slice);
    }
  }

  if(scalar(@output_slices) == 0) {
    $self->throw('You set the min_slice_length param to '.$min_size." but after filtering there were no input ids. Something is probably wrong");
  }

  say "Slice count post length filter: ".scalar(@output_slices);
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
  # Sort shortest firs
  foreach my $slice (sort { $a->length <=> $b->length } @$slices) {
    my $length = $slice->length;
    if($length + $total_length > $batch_target_size) {
      #this check handles rare cases where the shortest slice is > than the batch_target_size. if so, skip pushing to $all_batches array until $single_batch_array is populated with slice name
      #this is to avoid creating an extra job without input id
      if (scalar(@$single_batch_array)){
         push(@{$all_batches},[$single_batch_array]);
         $single_batch_array = [];
         $total_length = 0;
       }
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

=head2 _padded_slice_coord

 Arg [1]    : Bio::EnsEMBL::Feature, the object to get the start and end of the slice from
 Arg [2]    : Bio::EnsEMBL::Slice, the toplevel slice where Arg[1] is
 Arg [3]    : Int, the size of the padding
 Description: It calculate the start and end of the slice to be fetched for next analysis
              It makes sure that the slice is not out of bound
 Returntype : Array of Int, start and end of the slice
 Exceptions : None

=cut

sub _padded_slice_coord {
  my ($object, $slice, $padding) = @_;

  my $start = $object->seq_region_start-$padding;
  my $end = $object->seq_region_end+$padding;
  $start = 1 if ($start < 1);
  $end = $slice->end if ($end > $slice->end);
  return $start, $end;
}


=head2 fastq_range

 Arg [1]    : Path to the fastq file
 Arg [2]    : Batch size to split the sequences into
 Description: Take in a fastq file and output ranges based on a batch size
 Returntype : Array of Int, start and end of the index range
 Exceptions : Fasta file doesn't exist
              Fasta file has no headers

=cut

sub fastq_range {
  my ($self,$fastq_file,$batch_size) = @_;

  my $batch_array = [];
  unless(-e $fastq_file) {
    $self->throw("You have selected to generate a fastq range, but the fasta file doesn't exist. Path specified:\n".$fastq_file);
  }

  my $seq_count = `wc -l $fastq_file`;
  $seq_count =~ /^(\d+)/;
  $seq_count = $1;
  $seq_count = $seq_count / 4;

  unless($seq_count > 0) {
    $self->throw("You have selected to generate a fastq range, but the fastq file doesn't have any headers. Path specified:\n".$fastq_file);
  }

  my $start = 0;
  my $end = $start + $batch_size - 1;

  if($end > $seq_count - 1) {
    $end = $seq_count - 1;
  }

  push(@$batch_array,[$start,$end]);
  while($end + $batch_size + 1 < $seq_count) {
    $start = $end + 1;
    $end = $start + $batch_size - 1;
    push(@$batch_array,[$start,$end]);
  }

  if($end < $seq_count - 1) {
    $start = $end + 1;
    $end = $seq_count - 1;
    push(@$batch_array,[$start,$end]);
  }

  $self->param('inputlist', $self->_chunk_input_ids(1, $batch_array));
}

1;
