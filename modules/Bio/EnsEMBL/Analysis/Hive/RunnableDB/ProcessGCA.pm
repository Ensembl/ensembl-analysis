=head1 LICENSE

 Copyright [2021] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessGCA

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessGCA->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module dumps cds coord data from Ensembl dbs to file

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessGCA;

use warnings;
use strict;
use feature 'say';
use JSON;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
  };
}

=head2 run
p
 Arg [1]    : None
 Description: Run a python script to collect info from registry.


=cut

sub run {
  my ($self) = @_;

  my $python_script =  catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'info_from_registry', 'start_pipeline_from_registry.py');

  my $cmd = "python $python_script";
  say "Running: $cmd";


  my $json_output = `$cmd`;
  if ($? != 0) {
    $self->throw("Python script failed with exit code $? - output:\n$json_output");
  }

  # Convert JSON string to Perl hash
  my $result_hash = decode_json($json_output);

  # Store results for write_output
  $self->param('output_params', $result_hash);

  say "Assembly Accession: " . $result_hash->{'assembly_accession'};
  say "Clade: " . $result_hash->{'clade'};

}


=head2 write_output
p
 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $output_params = $self->param('output_params');
  $self->dataflow_output_id( $output_params, 1 );
}


1;
