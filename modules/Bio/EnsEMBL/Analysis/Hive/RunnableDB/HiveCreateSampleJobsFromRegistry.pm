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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateSampleJobsFromRegistry;

use strict;
use warnings;
use File::Spec::Functions;
use feature 'say';
use Bio::EnsEMBL::Analysis::Hive::DBSQL::TranscriptomicRegistryAdaptor;
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
  my $transcriptomic_adaptor = new Bio::EnsEMBL::Analysis::Hive::DBSQL::TranscriptomicRegistryAdaptor(
        -user   => $self->param('user'),
        -dbname => $self->param('registry_db'),
        -host   => $self->param('registry_host'),
        -port   => $self->param('registry_port'),
        -pass   => '',
        -driver => 'mysql',
    );
  my @results = $transcriptomic_adaptor->fetch_sample_by_taxon_id($self->param('taxon_id'));
  my @output_ids;
  my $column_hash = $self->param('column_names');
  my @samples;
  foreach my $result (@results) {
    push (@samples, [$result,'ENA']);
  }
  $self->param('inputlist',\@samples);
}



1;
