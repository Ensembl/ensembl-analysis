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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::GenerateGTFLoadingJobs

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::GenerateGTFLoadingJobs->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module takes a gtf file (or a set or dir of gtf files) and generates input ids based
on the number of genes in the each file. It generates a set of ranges based on a batch size
and then outputs the ranges (paired with the appropriate file). This is just a way to speed
up loading into a db since loading large gtfs (which often involves calculating things such
as the translation) can be really slow

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::GenerateGTFLoadingJobs;

use warnings;
use strict;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
                transcripts_per_batch => 5000,
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    transcripts_per_batch => 5000,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch input expects either a dir containing gtfs or an arrayref of
              gtf files
 Returntype : None
 Exceptions : Throws if it doesn't find input files

=cut

sub fetch_input {
  my ($self) = @_;

  my $gtf_files = [];
  if($self->param('gtf_dir')) {
    my $gtf_dir = $self->param('gtf_dir');
    unless(-e $gtf_dir) {
      $self->throw("A path to a gtf dir was porvided, but the dir does not exist. Path provided:\n".$gtf_dir);
    }
    $gtf_files = [glob($gtf_dir."/*.gtf")];
  } else {
    $gtf_files = $self->param('iid');
  }

  unless(scalar(@$gtf_files)) {
    $self->throw("Found no gtf files. This module expects either a dir containing gtf files or an iid with an arrayref of gtf files");
  }

  $self->param('gtf_files',$gtf_files);
}


=head2 run

 Arg [1]    : None
 Description: Run will go through the gtf files, count the genes and then make batches
 Returntype : None
 Exceptions : Throws if no output has been created has no files would have been find

=cut

sub run {
  my ($self) = @_;

  my $gtf_files = $self->param('gtf_files');
  my @output_ids;
  foreach my $gtf_file (@$gtf_files) {
    my $logic_name = (split '/', $gtf_file)[-1];
    $logic_name =~ s/.gtf$//;
    foreach my $chunk (@{$self->process_gtf_file($gtf_file)}){
      push(@output_ids,{gtf_array => $chunk, logic_name => $logic_name."_rnaseq_gene"});
    }
  }

  unless(scalar(@output_ids)) {
    $self->throw("No output ids created. Something went wrong");
  }
  $self->dataflow_output_id(\@output_ids, 2);
}


=head2 write_output

 Arg [1]    : None
 Description: Writes the output ids on branch '_branch_to_flow_to' usually 2
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $output_ids = $self->output();
  foreach my $output_id (@$output_ids) {
    $self->dataflow_output_id([{'iid' => $output_id}], $self->param('_branch_to_flow_to'));
  }
}


=head2 process_gtf_file

 Arg [1]    : String, path to the GTF file
 Description: Count the number of transcripts in the file and create a data structure to
              store the file path and the range to process
 Returntype : Arrayref of array, the array has three elements, filename, start, end
 Exceptions : Throws if the number of transcrits is 0

=cut

sub process_gtf_file {
  my ($self,$gtf_file) = @_;

  unless(-e $gtf_file) {
    $self->throw("The GTF file specified does not exist. Path checked:\n".$gtf_file);
  }

  my $batch_size = $self->param_required('transcripts_per_batch');
  my $batch_array = [];

  my $count_transcripts_command = "awk '{print \$3}' ".$gtf_file." | grep 'transcript' | wc -l";
  my $transcript_count = `$count_transcripts_command`;
  chomp($transcript_count);

  unless($transcript_count > 0) {
    $self->throw("Transcript count in file is < 1. File checked:\n".$gtf_file."\nCommandline used:\n".$count_transcripts_command);
  }

  my $start = 0;
  my $end = $start + $batch_size - 1;

  if($end > $transcript_count - 1) {
    $end = $transcript_count - 1;
  }

  push(@$batch_array,[$gtf_file,$start,$end]);
  while($end + $batch_size + 1 < $transcript_count) {
    $start = $end + 1;
    $end = $start + $batch_size - 1;
    push(@$batch_array,[$gtf_file,$start,$end]);
  }

  if($end < $transcript_count - 1) {
    $start = $end + 1;
    $end = $transcript_count - 1;
    push(@$batch_array,[$gtf_file,$start,$end]);
  }
  return($batch_array);
}

1;
