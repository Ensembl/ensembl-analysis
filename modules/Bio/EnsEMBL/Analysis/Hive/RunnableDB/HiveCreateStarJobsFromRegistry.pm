=head1 LICENSE
 
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateStarJobsFromRegistry;

use strict;
use warnings;
use File::Spec::Functions;
use feature 'say';
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
                'compression_ratio' => 3,
                'target_batch_size' => 10000000000,
 Returntype : Hashref
 Exceptions : None

=cut

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
  my $registry_adaptor = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
        -user   => $self->param('user'),
        -dbname => $self->param('registry_db'),
        -host   => $self->param('registry_host'),
        -port   => $self->param('registry_port'),
        -pass   => '',,
        -driver => 'mysql',,
    );
  my @results = $registry_adaptor->fetch_short_reads_by_species_id($self->param('taxon_id'));
  my @output_ids;
  my $column_hash = $self->param('column_names');
  my $samples_hash = {};
  foreach my $result (@results) {
    my @arr = split(/\t/,$result);
    $arr[1] =~ s/^\s+|\s+$//g;
    my $sample_id = $arr[1];
    my $sample_name = $arr[2];
    $arr[3] =~ s/^\s+|\s+$//g;
    #Separate paired fastq files to display in right format for Star alignment
    my @reads = split(/;/,$arr[3]);
    if (scalar(@reads) > 1){#if paired ended reads found
      foreach my $input_read(@reads){
        #Add read to input parameter
        if(exists $samples_hash->{$sample_id}) {
          push(@{$samples_hash->{$sample_id}->{'files'}},$input_read);
        } 
        else {
          $samples_hash->{$sample_id}->{'files'} = [$input_read];
          $samples_hash->{$sample_id}->{$self->param('sample_column')} = $sample_name;
          $samples_hash->{$sample_id}->{$self->param('sample_id_column')} = $sample_id;
        }
      }
    }
  }
  $self->param('inputlist',$self->batch_samples($samples_hash));
}


=head2 batch_samples

 Arg [1]    : Hashref, containing the information from the 'csvfile_table'
 Description: Batch files depending on their size
 Returntype : Arrayref of array containing sample names to be batched together
 Exceptions : None

=cut

sub batch_samples {
  my ($self,$samples_hash) = @_;

  my $fastq_dir = $self->param_required('input_dir');
  $fastq_dir =~ s/^\s+|\s+$//g;
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


=head2 build_sorted_batches

 Arg [1]    : Hashref, containing the information from the 'csvfile_table'
 Arg [2]    : Hashref, containing the the size of each file
 Description: Build the size sorted batches
 Returntype : Arrayref of array containing sample names to be batched together
 Exceptions : None

=cut

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
        push(@{$all_batches},[$single_batch_array,1]);
        $single_batch_array = [];
        $total_size = 0;
      }
    }
    push(@{$single_batch_array}, $sample);
    $total_size += $file_size;
  }
  push(@{$all_batches},[$single_batch_array,1]);
  return($all_batches);
}

1;
