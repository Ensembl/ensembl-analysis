=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar->new( );

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


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
 Returntype : Hashref
                threads => 1,
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    threads => 1,
  }
}


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

  my $input_ids = $self->param('SM');
  $self->say_with_header('Found '.scalar(@$input_ids).' input ids');
  foreach my $input_id (@$input_ids) {
    my $sample_id = $input_id->{'ID'};
    $self->say_with_header("Processing sample: $sample_id");
    my $files = $input_id->{'files'};
    my $file1 = ${$files}[0];
    my $file2 = ${$files}[1];

    $self->say_with_header("Found file: $file1");
    my $filepath1 = $self->param('input_dir').'/'.$file1;
    $self->throw("Fastq file ".$filepath1." not found\n") unless ( -e $filepath1 );

    my $filepath2 = "";
    if($file2) {
      $self->say_with_header("Found paired file: $file2");
      $filepath2 = $self->param('input_dir').'/'.$file2;
      $self->throw("Fastq file ".$filepath2." not found\n") unless ( -e $filepath2 );
    }

    my $program = $self->param('short_read_aligner');
    $self->throw("Star program not defined in analysis\n") unless (defined $program);

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Star->new
    (
     -analysis       => $self->create_analysis,
     -program        => $program,
     -options        => $self->param('short_read_aligner_options'),
     -outdir         => $self->param('output_dir'),
     -genome_dir     => $self->param('genome_dir'),
     -genome     => $self->param('genome_dir')."/Genome",
     -sample_name    => $sample_id,
     -fastq          => $filepath1,
     -fastqpair      => $filepath2,
     -threads        => $self->param('num_threads'),
    );
    if ($self->param_is_defined('rg_lines')) {
      $runnable->rg_lines($self->param('rg_lines'));
    }
    else {
      $runnable->rg_lines("ID:$sample_id\tSM:".$input_id->{SM}."\tPL:ILLUMINA");
    }
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
    $self->say_with_header("Output file: ".$output_file);
    $self->dataflow_output_id([{'iid' => $output_file}], $self->param('_branch_to_flow_to'));
  }
}

1;
