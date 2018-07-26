=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSequences

=head1 SYNOPSIS


=head1 DESCRIPTION

Base module to load data into Hive custom tables. Modules should override
create_row_data which will be called before loading data into the table.
the accession field will be return in an arrayref on branch 2

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::CalculateLowCoverageSlices;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;
use File::Basename;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    window_jump          => 10000,
    window_size          => 250000,
    min_slice_split_size => 5000000,
    min_depth            => 50,
    min_cov              => 500,
    base_length          => 1000000,
    max_slice_size       => 5000000,
    slice_overlap        => 2500000,
    sample_count         => 1,
    data_flow_branch     => 2,
  }
}



sub fetch_input {
  my ($self) = @_;

  my $bam_file = $self->param_required('alignment_bam_file');
  unless(-e $bam_file) {
    $self->throw("Could not find the BAM file on the following path:\n".$bam_file);
  }

  my $dna_dba = $self->hrdb_get_dba($self->param_required('dna_db'));
  $self->hrdb_set_con($dna_dba,'dna_db');
  my $slice_adaptor = $dna_dba->get_SliceAdaptor();
  my $slices = [];
  my $slice_names = $self->param_required('iid');
  foreach my $slice_name (@$slice_names) {
    say "Initial slice name: ".$slice_name;
    my $slice = $slice_adaptor->fetch_by_name($slice_name);
    push(@{$slices},$slice);
  }

  unless(scalar(@$slices)) {
    $self->warning("Found no slices from the input array");
    $self->input_job->autoflow(0);
    $self->complete_early('There are no slices to process');
  }

  $self->input_slices($slices);
}


sub run {
  my ($self) = @_;

  my $max_slice_size = $self->param('max_slice_size');
  my $slice_overlap = $self->param('slice_overlap');
  my $slice_adaptor = $self->hrdb_get_con('dna_db')->get_SliceAdaptor();

  my $output_slices = [];
  my $slices = $self->input_slices();
  foreach my $slice (@$slices) {
    say "Processing slice: ".$slice->name;
    my $slice_names = $self->process_slice($slice);
    foreach my $slice_name (@$slice_names) {
      my $new_slice = $slice_adaptor->fetch_by_name($slice_name);
      if($new_slice->length < $max_slice_size) {
        push(@$output_slices,$new_slice);
      } else {
        my $split_slices = split_Slices([$new_slice],$max_slice_size,$slice_overlap);
        unless(scalar(@$split_slices)) {
          $self->throw("Tried to split a slice but the split array was empty. Original slice name: ".$new_slice->name);
        }
        foreach my $split_slice (@$split_slices) {
          push(@$output_slices,$split_slice);
        }
      }
    }

    foreach my $output_slice (@$output_slices) {
      say "Created: ".$output_slice->name;
    }
    push(@$output_slices,$slice);
  } # end foreach my $slice (@$slices) {


    foreach my $output_slice (@$output_slices) {
      say "Pre output: ".$output_slice->name;
    }

  $self->output($output_slices);
}


sub write_output {
  my ($self) = @_;

  foreach my $output_slice (@{$self->output}) {
    my $output_hash->{'iid'} = $output_slice->name;
    $output_hash->{'alignment_bam_file'} = $self->param('alignment_bam_file');
    $self->dataflow_output_id($output_hash, $self->param('_branch_to_flow_to'));
  }

}


sub process_slice {
  my ($self,$slice) = @_;

  my $bam_file = $self->param('alignment_bam_file');
  my $window_jump = $self->param('window_jump');
  my $window_size = $self->param('window_size');
  my $min_slice_split_size = $self->param('min_slice_split_size');
  my $min_depth = $self->param('min_depth');
  my $min_cov = $self->param('min_cov');
  my $base_length = $self->param('base_length');
  my $sample_count= $self->param('sample_count');

  my $slice_names = [];
  my $assembly_name = $slice->coord_system->version;
  my $slice_type = $slice->coord_system->name;
  my $slice_seq_region_name = $slice->seq_region_name;
  my $slice_end = $slice->length;
  my $current_slice_start = 1;

  say "Getting depth: ";
  my @depth_array;
  open(SAMTOOLS_DEPTH,"samtools depth -Q 30 -r ".$slice_seq_region_name.":1-".$slice_end." ".$bam_file." |");
  while(my $depth_line = <SAMTOOLS_DEPTH>) {
    chomp($depth_line);
    my ($region,$base,$cov) = split('\t',$depth_line);
    $depth_array[$base] = $cov;
  }
  close SAMTOOLS_DEPTH;
  say "Got depth";

  for(my $i=$base_length; $i<($slice_end - $window_size); $i += $window_jump) {
    my $window_start = $i;
    my $window_end = $i + $window_size;
    if($window_end > $slice_end) {
      $window_end = $slice_end;
    }

    if($self->low_coverage_window($slice_seq_region_name,$window_start,$window_end,\@depth_array)) {
      my $midpoint = int(($window_start + $window_end) / 2);
      my $new_slice_name = $slice_type.":".$assembly_name.":".$slice_seq_region_name.":".$current_slice_start.":".$midpoint.":1";
      push(@$slice_names,$new_slice_name);
      $current_slice_start = $midpoint + 1;
      if($current_slice_start >= $slice_end) {
        $current_slice_start = 0;
        last;
      } elsif($current_slice_start + $base_length >= $slice_end) {
        $new_slice_name = $slice_type.":".$assembly_name.":".$slice_seq_region_name.":".$current_slice_start.":".$slice_end.":1";
        push(@$slice_names,$new_slice_name);
        $current_slice_start = 0;
        last;
      } else {
        $i = $midpoint + $base_length - $window_jump;
      }
    }
  } # End for(my $i=1;

  if($current_slice_start) {
    my $new_slice_name = $slice_type.":".$assembly_name.":".$slice_seq_region_name.":".$current_slice_start.":".$slice_end.":1";
    push(@$slice_names,$new_slice_name);
  }

  undef(@depth_array);
  return($slice_names);
}

sub low_coverage_window {
  my ($self,$slice_seq_region_name,$window_start,$window_end,$depth_array) = @_;

  my $min_depth = $self->param('min_depth');
  my $min_cov = $self->param('min_cov');
  my $base_length = $self->param('base_length');

  my $passing_bases = 0;
  for (my $i=$window_start-1; $i<$window_end-1; $i++) {
    if(${$depth_array}[$i] && ${$depth_array}[$i] > $min_depth) {
      $passing_bases++;
      if($passing_bases >= $min_cov) {
        say "Region passed threshold, moving to next window";
        return(0);
      }
    }
  }
  say "Region is low coverage/depth: ".$window_start.":".$window_end;
  return(1);
}


sub input_slices {
  my ($self,$slices) = @_;

  if($slices) {
    $self->param('_input_slices',$slices);
  }

  return($self->param('_input_slices'));
}


1;
