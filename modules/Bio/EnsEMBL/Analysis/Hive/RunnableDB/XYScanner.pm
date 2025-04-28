=head1 LICENSE

 Copyright [2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::XYScanner

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::XYScanner->new();

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::XYScanner;

use warnings;
use strict;
use feature 'say';

use File::Spec::Functions;
use Bio::EnsEMBL::Hive::Utils qw(destringify);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
  }
}


sub run {
  my ($self) = @_;
  my $xy_scanner_path = $self->param_required('xy_scanner_path');
  my $x_marker_fasta_path = $self->param_required('x_marker_fasta_path');
  my $y_marker_fasta_path = $self->param_required('y_marker_fasta_path');
  my $target_genome_index = $self->param_required('target_genome_index');
  my $output_dir = $self->param_required('output_dir');
  my $cmd = "python ".$xy_scanner_path." --genome_index ".$target_genome_index.
                                          " --x_marker_file ".$x_marker_fasta_path.
                                          " --y_marker_file ".$y_marker_fasta_path.
                                          " --output_dir ".$output_dir;

  say "Running command:";
  say $cmd;

  my $result = system($cmd);
  if($result) {
    $self->throw("Non-zero exit code encountered from the command");
  }
  my $output_file = "$output_dir/xy_scanner.out";
  my $output_params = destringify($self->input_job->input_id);
  if (-e $output_file) {  # Check if the file exists
    open(IN, $output_file) or $self->throw("Failed to open $output_file: $!");
    #open(IN,$output_dir."/xy_scanner.out");
    my $xy_result = <IN>;
    chomp $xy_result;
    close IN;

    #   my $output_params = destringify($self->input_job->input_id);
    $output_params->{'xy_scanner'} = $xy_result;
     } else {
    $self->warning("Output file $output_file does not exist. Setting xy_scanner to empty string.");
    $output_params->{'xy_scanner'} = undef;
  }
  $self->input_job->input_id($output_params);
}


sub write_output {
  my ($self) = @_;
}

1;
