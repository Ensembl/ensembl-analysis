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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM;

use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::Runnable::BWA2BAM;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my ($self) = @_;
  my $program = $self->param('short_read_aligner');
  $self->throw("BWA program not defined in analysis\n") unless (defined $program);
  my $filename = $self->param('input_dir').'/'.$self->input_id;
  $self->throw("Fastq file $filename not found\n") unless (-e $filename);
  my $fastqpair;

  my $method;
  if ( $self->param('paired') ) {
    $method = ' sampe '.$self->param('sampe_options');
    $fastqpair = $self->param('input_dir') .'/'.$$self->param('fastq_pair');
    $self->throw('You are missing the paired reads from '.$fastqpair."\n") unless (-e $fastqpair);
  } else {
    $method = ' samse '.$self->param('samse_options');
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BWA2BAM->new
    (
     -analysis   => $self->create_analysis,
     -program    => $program,
     -fastq      => $filename,
     -fastqpair  => $fastqpair,
     -options    => $method,
     -outdir     => $self->param('output_dir'),
     -genome     => $self->param('genomefile'),
	 -samtools => $self->param('samtools_path'),
     -header => $self->param('header_file'),
     -min_mapped => $self->param('min_mapped'),
     -min_paired => $self->param('min_paired'),
    );
  $self->runnable($runnable);
}
