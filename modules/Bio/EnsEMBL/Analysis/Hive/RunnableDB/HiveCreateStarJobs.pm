=head1 LICENSE
 

Copyright [2019-2020] EMBL-European Bioinformatics Institute

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

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateStarJobs;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;
use Data::Dumper;
use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my $self = shift;

  return {
    %{$self->SUPER::param_defaults},
    'compression_ratio' => 3,
    'target_batch_size' => 10000000000,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Creates input id based on a custom table 'csvfile_table' in the hive database
              It will generate the parameters for STAR based on the data for each file
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my $self = shift;

  my $table_adaptor = $self->db->get_NakedTableAdaptor;
  $table_adaptor->table_name($self->param('csvfile_table'));

  my $results = $table_adaptor->fetch_all();
  my @output_ids;
  my $column_hash = $self->param('column_names');
  my $samples_hash = {};

  foreach my $result (@$results) {
    my $sample_id = $result->{$self->param('sample_id_column')};
    my $sample_name = $result->{$self->param('sample_column')};
    if ($result->{$self->param('filename_column')} !~ /_split/){
      if(exists $samples_hash->{$sample_id}) {
        push(@{$samples_hash->{$sample_id}->{'files'}},$result->{$self->param('filename_column')});
      } else {
        $samples_hash->{$sample_id}->{'files'} = [$result->{$self->param('filename_column')}];
        $samples_hash->{$sample_id}->{$self->param('sample_column')} = $sample_name;
        $samples_hash->{$sample_id}->{$self->param('sample_id_column')} = $sample_id;
      }
    }
  }

#  my $output_ids = $self->batch_samples($samples_hash);
  $self->param('inputlist',$self->batch_samples($samples_hash));
}


sub batch_samples {
  my ($self,$samples_hash) = @_;

  say Dumper($samples_hash);
  my $fastq_dir = $self->param_required('input_dir');

  my $file_sizes = {};
  foreach my $sample_id (keys(%$samples_hash)) {
    my $sample = $samples_hash->{$sample_id};
    my $file_name = ${$sample->{'files'}}[0];
    my $file_path = catfile($fastq_dir,$file_name);
    my $file_size = `stat -c %s $file_path`;
    chomp($file_size);

    unless($file_name =~ /\.gz$/) {
      $file_size = $file_size / $self->param_required('compression_ratio');
    }
    $file_sizes->{$sample_id} = $file_size;
  }

  my $output_ids = $self->build_sorted_batches($samples_hash,$file_sizes);
  return($output_ids);
}


sub build_sorted_batches {
  my ($self,$samples_hash,$file_sizes) = @_;

  my $all_batches = [];
  my $single_batch_array = [];
  my $target_batch_size = $self->param('target_batch_size');
  my $total_size = 0;
  foreach my $sample_id (sort { $file_sizes->{$a} <=> $file_sizes->{$b} } keys %$file_sizes) {
    my $sample = $samples_hash->{$sample_id};
    my $size_multiplier = 1;
    if(scalar(@{$sample->{'files'}}) > 1) {
      $size_multiplier = 2;
    }

    my $file_size = $file_sizes->{$sample_id} * $size_multiplier;
    if($file_size + $total_size > $target_batch_size) {
      if(scalar(@$single_batch_array)){
        push(@{$all_batches},[$single_batch_array]);
        $single_batch_array = [];
        $total_size = 0;
      }
    }
    push(@{$single_batch_array}, $sample);
    $total_size += $file_size;
  }
  push(@{$all_batches},[$single_batch_array]);
  return($all_batches);
}

1;
