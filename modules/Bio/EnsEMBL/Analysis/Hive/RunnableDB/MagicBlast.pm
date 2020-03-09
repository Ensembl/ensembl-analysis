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

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::MagicBlast

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::MagicBlast->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses MagicBlast to align fastq to a genomic sequence

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::MagicBlast;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::MagicBlast;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _branch_to_flow_to => 2,
    threads => 1,
    create_sorted_bam => 1,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for MagicBlast
 Returntype : None
 Exceptions : Throws if 'filename' does not exist
              Throws if 'fastqpair' does not exist
              Throws if 'wide_short_read_aligner' is not defined

=cut

sub fetch_input {
  my ($self) = @_;

  my $input_ids = $self->param('SM');
  say "Found ".scalar(@$input_ids)." input ids";
  foreach my $input_id (@$input_ids) {
    my $sample_id = $input_id->{'ID'};
    say "Processing sample: ".$sample_id;
    my $files = $input_id->{'files'};
    my $file1 = ${$files}[0];
    my $file2 = ${$files}[1];

    say "Found file: ".$file1;
    my $filepath1 = $self->param('input_dir').'/'.$file1;
    $self->throw("Fastq file ".$filepath1." not found\n") unless ( -e $filepath1 );

    my $filepath2 = "";
    if($file2) {
      say "Found paired file: ".$file2;
      $filepath2 = $self->param('input_dir').'/'.$file2;
      $self->throw("Fastq file ".$filepath2." not found\n") unless ( -e $filepath2 );
    }

    my $program = $self->param('short_read_aligner');
    $self->throw("MagicBlast program not defined in analysis\n") unless (defined $program);

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::MagicBlast->new
    (
     -analysis       => $self->create_analysis,
     -program        => $program,
     -options        => $self->param('short_read_aligner_options'),
     -outdir         => $self->param('output_dir'),
     -genome         => $self->param('genome_file'),
     -sample_name    => $sample_id,
     -fastq          => $filepath1,
     -fastqpair      => $filepath2,
     -threads        => $self->param('num_threads'),
     -create_sorted_bam => $self->param('create_sorted_bam'),
    );
    $self->runnable($runnable);
  }

}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $output_files = $self->output;
  foreach my $output_file (@$output_files) {
    say "Output file: ".$output_file;
    $self->dataflow_output_id([{'iid' => $output_file}], $self->param('_branch_to_flow_to'));
  }
}

1;
