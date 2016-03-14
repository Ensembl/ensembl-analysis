# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BWA




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::BWA->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses BWA to align fastq to a genomic sequence

=head1 CONTACT

Post general queries to B<http://lists.ensembl.org/mailman/listinfo/dev>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA;

use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::Runnable::BWA;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my ($self) = @_;
  my $filename =  $self->param('wide_input_dir') ."/" .$self->param('filename');
  $self->throw("Fastq file  $filename not found\n") unless ( -e $filename );
  my $program = $self->param('wide_short_read_aligner');
  $self->throw("BWA program not defined in analysis \n") unless (defined $program);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BWA->new
    (
     -analysis => $self->create_analysis,
     -program  => $program,
     -options  => $self->param('short_read_aligner_options'),
     -outdir   => $self->param('wide_output_dir'),
     -genome   => $self->param('wide_genome_file'),
     -fastq    => $filename,
    );
  $self->runnable($runnable);
}

# override write output as we have nothing for the db
sub write_output {
  my ($self) = @_;

  $self->dataflow_output_id({fastq => {filename => $self->param('filename'), is_mate_1 => $self->param('is_mate_1')}}, 1);
}

1;
