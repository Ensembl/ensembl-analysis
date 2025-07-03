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
use IPC::Open3;
use File::Spec::Functions qw(catfile);
use Symbol qw(gensym);
use Data::Dumper;


use Bio::EnsEMBL::Utils::Exception qw(throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
      %{$self->SUPER::param_defaults},  # inherit defaults from HiveBaseRunnableDB or parent
  };
}

=head2 run
p
 Arg [1]    : None
 Description: Run a python script to collect info from registry.


=cut


sub run {
  my ($self) = @_;

  my $python_script =  catfile(
    $self->param('enscode_root_dir'),
    'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'info_from_registry', 'start_pipeline_from_registry.py'
  );

  # Get all params as hashref - this method is not allowed :(((((((
  my $all_param = $self->param;
  say Dumper($all_param);  # returns a hashref of all parameters

  # Encode all params to JSON
  my $input_json = encode_json($all_param);
  say "Input JSON to python:\n$input_json";

  my $cmd = "python $python_script";

  say "Running: $cmd with input JSON: $input_json";

  my $err = gensym;  # for capturing stderr

  my $pid = open3(my $in, my $out, $err, $cmd);

  # Write JSON input to python stdin and close
  print $in $input_json;
  close $in;

  # Read stdout and stderr from Python
  my $json_output = do { local $/; <$out> };
  my $stderr_output = do { local $/; <$err> };

  waitpid($pid, 0);
  my $exit_code = $? >> 8;

  if ($exit_code != 0) {
    $self->throw("Python script failed with exit code $exit_code\nSTDERR:\n$stderr_output\nSTDOUT:\n$json_output");
  }

  my $result_hash;
  eval {
    $result_hash = decode_json($json_output);
  };
  if ($@) {
    $self->throw("Failed to decode JSON output from Python script: $@\nOutput was:\n$json_output\nSTDERR:\n$stderr_output");
  }

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
