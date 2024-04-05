=pod

=head1 NAME

    Bio::EnsEMBL::Analysis::Hive::RunnableDB::FileFactory;

=head1 SYNOPSIS

    standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::FileFactory \
                    --inputfile 'core_genes_for_deletion_filename' \
                    --flow_into "{ 2 => { 'mysql://ensadmin:${ENSADMIN_PSW}@127.0.0.1:2914/lg4_compara_families_70/meta' => {'meta_key'=>'module_name','meta_value'=>'#_0#'} } }""

=head1 DESCRIPTION

    This is a generic RunnableDB module for creating batches of similar jobs using dataflow mechanism
    (a fan of jobs is created in one branch and the funnel in another).
    Make sure you wire this building block properly from outside using the 'file' parameter.

    You can supply as parameter one source of ids from which the batches will be generated:

        param('inputfile');  The list is contained in a file whose name is supplied as parameter: 'inputfile' => 'myfile.txt'

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2024] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    Please subscribe to the Hive mailing list:  http://listserver.ebi.ac.uk/mailman/listinfo/ehive-users  to discuss Hive-related questions or to be notified of our updates

=cut


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::FileFactory;

use strict;
use warnings;

use File::Spec::Functions qw(tmpdir);
use base ('Bio::EnsEMBL::Hive::Process');


=head2 param_defaults

    Description : Implements param_defaults() interface method of Bio::EnsEMBL::Hive::Process that defines module defaults for parameters.

=cut

sub param_defaults {

    return {
        'max_chunk_size'  => 500,
        'output_prefix'     => 'my_chunk_',
        'output_suffix'     => '.txt',

        'inputfile'         => undef,
        'output_dir'        => tmpdir(),

        'fan_branch_code'   => 2,
    };
}


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).
                    Because we want to stream the data more efficiently, all functionality is in write_output();

=cut

sub run {
}


=head2 write_output

    Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
                    The main bulk of this Runnable's functionality is here.
                    Iterates through all output ids in output_ids, splits them into separate files ("chunks") using a cut-off length and dataflows one job per chunk.

=cut

sub write_output {
  my $self = shift @_;

  my $fan_branch_code         = $self->param('fan_branch_code');

  my $max_chunk_size    = $self->param('max_chunk_size');
  my $output_prefix       = $self->param('output_prefix');
  my $output_suffix       = $self->param('output_suffix');

  my $chunk_number = 1;   # counts the chunks
  my $chunk_size   = 0;   # number of sequences in the current chunk
  my $chunk_filename   = $self->param('output_dir')."/".$output_prefix.$chunk_number.$output_suffix;

  open(my $chunk_file,'>',$chunk_filename) or die "Could not open file '$chunk_filename' for writing $!";

  open INPUTFILE, $self->param('inputfile') or die $!;
  while (my $output_id = <INPUTFILE>) {
    chomp($output_id);
    print $chunk_file $output_id."\n";
    $chunk_size++;

    if ($chunk_size >= $max_chunk_size) {
      # dataflow the current chunk
      $self->dataflow_output_id( {
                'file' => $chunk_filename,
                'chunk_number' => $chunk_number,
                'chunk_size' => $chunk_size
              },
              $fan_branch_code);

      close $chunk_file;

      # start writing to the next one
      $chunk_size = 0;
      $chunk_number++;
      $chunk_filename   = $self->param('output_dir')."/".$output_prefix.$chunk_number.$output_suffix;
      open($chunk_file,'>',$chunk_filename) or die "Could not open file '$chunk_filename' for writing $!";
    }
  }

  if ($chunk_size) { # flush the last chunk
    $self->dataflow_output_id( {
                'file' => $chunk_filename,
                'chunk_number' => $chunk_number,
                'chunk_size' => $chunk_size
              },
              $fan_branch_code);
      close $chunk_file;
    } else {
      unlink $chunk_filename unless (stat($chunk_filename))[7];
    }
}


1;
