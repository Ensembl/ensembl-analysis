=head1 LICENSE

 Copyright [2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::Star

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::Minimap2;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for minimap2
 Returntype : None
 Exceptions : Throws if 'genome_file' does not exist
              Throws if 'input_file' does not exist

=cut

sub fetch_input {
  my ($self) = @_;

  my $output_dba = $self->hrdb_get_dba($self->param('target_db'));

  if($self->param('dna_db')) {
    my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'target_db');

  my $genome_file = $self->param('genome_file');
  unless(-e $genome_file) {
    $self->throw("Could not find the genome file. Path used:\n".$genome_file);
  }

  my $input_file = $self->param('input_file');
  unless(-e $input_file) {
    $self->throw("Could not find the input file. Path used:\n".$input_file);
  }

  my $program = $self->param('minimap2_path');
  my $paftools =  $self->param('paftools_path');

  unless($program) {
    $program = "minimap2";
  }

  unless($paftools) {
    $paftools = "paftools.js";
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis       => $self->create_analysis,
       -program        => $program,
       -paftools_path  => $paftools,
#       -options        => $self->param('minimap2_options'),
       -genome_file    => $genome_file,
       -input_file     => $input_file,
       -database_adaptor => $output_dba,
    );

  $self->runnable($runnable);
}

=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $target_dba = $self->hrdb_get_con("target_db");
  my $gene_adaptor = $target_dba->get_GeneAdaptor();

  my $genes = $self->output();
  foreach my $gene (@$genes) {
    $gene_adaptor->store($gene);
  }

}

1;
