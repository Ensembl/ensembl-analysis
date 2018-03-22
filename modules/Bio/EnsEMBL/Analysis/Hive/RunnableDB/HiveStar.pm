=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
 Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::Star->new( );

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar;

use warnings;
use strict;

use Bio::EnsEMBL::Analysis::Runnable::Star;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for STAR
 Returntype : None
 Exceptions : Throws if 'filename' does not exist
              Throws if 'fastqpair' does not exist
              Throws if 'wide_short_read_aligner' is not defined

=cut

sub fetch_input {
  my ($self) = @_;
  my $filename = $self->param('wide_input_dir').'/'.$self->param('filename');
  $self->throw("Fastq file $filename not found\n") unless ( -e $filename );
  my $fastqpair = $self->param('wide_input_dir').'/'.$self->param('fastqpair');
  $self->throw("Fastq file $fastqpair not found\n") unless ( -e $fastqpair );
  my $program = $self->param('wide_short_read_aligner');
  $self->throw("Star program not defined in analysis\n") unless (defined $program);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Star->new
    (
     -analysis       => $self->create_analysis,
     -program        => $program,
     -options        => $self->param('short_read_aligner_options'),
     -outdir         => $self->param('output_dir'),
     -workdir        => $self->param('temp_dir'),
     -genome         => $self->param('wide_genome_file'),
     -fastq          => $filename,
     -fastqpair      => $fastqpair,
     -sam_attributes => $self->param('sam_attributes'),
     -decompress     => $self->param('decompress'),
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

# I need to overwrite branch 1 as I'm using the accu table and the next analysis is on branch 1
  $self->dataflow_output_id({filename => $self->output->[0]}, 1);
}

1;
