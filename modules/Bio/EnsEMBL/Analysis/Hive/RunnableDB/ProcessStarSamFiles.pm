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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessStarSamFiles;

use warnings;
use strict;
use feature 'say';
use File::Basename;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    '_branch_to_flow_to' => 1,
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

  my $sam_file = $self->param('iid');
  unless(-e $sam_file) {
    $self->throw("The sam file does not appear to exist. Path used:\n".$sam_file);
  }

  $self->param('_sam_file',$sam_file);
}


sub run {
  my ($self) = @_;

  my $sam_file = $self->param('_sam_file');
  $self->process_sam($sam_file);
}


sub write_output {
  my ($self) = @_;

  my $bam_file = $self->param('_bam_file');
  $self->dataflow_output_id({filename => $bam_file}, 1);
}


sub process_sam {
  my ($self,$sam_file) = @_;

  my $unsorted_bam_file = $sam_file;
  if($unsorted_bam_file =~ /_Aligned\.out\.sam/) {
    $unsorted_bam_file =~ s/_Aligned\.out\.sam/\.unsorted\.bam/;
  } elsif($unsorted_bam_file =~ /\.sam/) {
    $unsorted_bam_file =~ s/\.sam/\.unsorted\.bam/;
  } else {
    $self->throw("Unexpected sam filename pattern. Filename: ".$unsorted_bam_file);
  }

  my $bam_command = 'samtools view -S -b '.$sam_file.' > '.$unsorted_bam_file;
  my $result = system($bam_command);
  if($result) {
    $self->throw("Error running bam conversion. Commandline used:\n".$bam_command);
  }

  my $sorted_bam_file = $unsorted_bam_file;
  $sorted_bam_file =~ s/\.unsorted//;
  my $sorted_bam_command = 'samtools sort '.$unsorted_bam_file.' -o '.$sorted_bam_file;
  $result = system($sorted_bam_command);
  if($result) {
    $self->throw("Error running sorting bam. Commandline used:\n".$sorted_bam_command);
  }
  system("rm ".$unsorted_bam_file);

  my $index_bam_command = 'samtools index '.$sorted_bam_file;
  $result = system($index_bam_command);
  if($result) {
    $self->throw("Error running indexing bam. Commandline used:\n".$index_bam_command);
  }

  $self->param('_bam_file',$sorted_bam_file);
}


1;
