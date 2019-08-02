=head1 LICENSE
 
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute
 
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateStarTissueJobs;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;
use Data::Dumper;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my $self = shift;

  return {
    %{$self->SUPER::param_defaults},
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Creates input id based on a custom table 'csvfile_table' in the hive database
              It will generate the parameters for BWA based on the data for each file
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
  my $output_dir = $self->param('output_dir');
  foreach my $result (@$results) {
    my $sample_name = $result->{$self->param('sample_column')};
    my $sample_id = $result->{$self->param('sample_id_column')};
    if(exists $samples_hash->{$sample_name}->{$sample_id}) {
      next;
    }

   $samples_hash->{$sample_name}->{$sample_id} = 1;

    say "FERGAL SAMPLE RES: ".$sample_name;
#    $samples_hash->{$sample_name} = 1;
    if(exists $samples_hash->{$sample_name}) {
      push(@{$samples_hash->{$sample_name}->{'files'}},catfile($output_dir,$sample_id.'.bam'));
    } else {
      $samples_hash->{$sample_name}->{'files'} = [catfile($output_dir,$sample_id.'.bam')];
    }
  }

  $self->param('_samples_hash',$samples_hash);
#  $self->param('inputlist',$self->batch_samples($samples_hash));
}


sub run {
  my ($self) = @_;
}


sub write_output {
  my ($self) = @_;

  my $samples_hash = $self->param('_samples_hash');
  foreach my $sample_name (keys(%$samples_hash)) {
#    my $files = $samples_hash->{$sample_name}->{'files'};
#    foreach my $file (@$files) {
#    $self->dataflow_output_id({'samplename' => {'files' => $files,'name' => $sample_name}},1);
     say "FERGAL DEBUG: ".$sample_name;
     say Dumper($samples_hash);
     $self->dataflow_output_id({'tissue_data' => {$sample_name => $samples_hash->{$sample_name}->{'files'} }},2);
#  }
  }
}

1;
