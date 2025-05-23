=head1 LICENSE

 Copyright [2020-2024] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::Scallop

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::Scallop->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses Scallop to construct a gtf file from a sorted bam file

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::Scallop;

use warnings;
use strict;
use feature 'say';
use File::Basename;

use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Analysis::Runnable::Scallop;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    input_file_extensions => ['_Aligned.sortedByCoord.out.bam','.mb.sorted.bam','.bam'],
    num_threads => 1,
    delete_input_file => 0,
	skip_analysis => 0,
  }
}



=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for Stringtie
 Returntype : None
 Exceptions : Throws if no input ids

=cut

sub fetch_input {
  my ($self) = @_;

  if ($self->param('skip_analysis')) {
    $self->complete_early('I was asked to skip this analysis');
  }

  my $input_ids = [$self->param('iid')];
  unless(scalar(@$input_ids)) {
    $self->throw("Found no input ids, something is wrong");
  }

  my $program = $self->param('scallop_path');

  unless($program) {
    $program = "scallop";
  }

  my $output_dir = $self->param_required('output_dir');
  unless(-e $output_dir) {
    my $result = system('mkdir -p '.$output_dir);
    if($result) {
      $self->throw("Could not create an output dir. Path used:\n".$output_dir);
    }
  }

  foreach my $input_file (@$input_ids) {
    my $analysis = $self->create_analysis;
    my $sample_name = $self->get_sample_name($input_file,$self->param_required('csv_summary_file'));
    unless($sample_name) {
      $sample_name = $self->get_sample_name($input_file,$self->param_required('csv_summary_file_genus'));
      unless($sample_name) {
        $self->throw("Checked both the normal and genus csv files and failed to match the input file to a line, so could not retrieve the sample name");
      }
    }

    $analysis->logic_name($sample_name."_scallop");
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Scallop->new(
         -analysis          => $analysis,
         -program           => $program,
         -input_file        => $input_file,
         -output_dir        => $output_dir,
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


#  my $target_dba = $self->hrdb_get_con("target_db");
#  my $gene_adaptor = $target_dba->get_GeneAdaptor();

#  my $genes = $self->output();
#  say "Total genes to output: ".scalar(@$genes);
#  foreach my $gene (@$genes) {
#    $gene_adaptor->store($gene);
#  }

}


sub get_sample_name {
  my ($self,$input_file_path,$csv_file) = @_;

  my $input_file = basename($input_file_path);
  chomp($input_file);

  my $file_extensions_to_remove = $self->param('input_file_extensions');

  foreach my $extension (@$file_extensions_to_remove) {
    $input_file =~ s/$extension//;
  }

  say "Searching for string ".$input_file." in the following csv file:\n".$csv_file;

  my $sample_name;

  unless(open(IN,$csv_file)) {
    $self->throw("Could not find the summary file to get the sample name. Path used:\n".$csv_file);
  }


  while(<IN>) {
    my $line = $_;
    unless($line =~ /$input_file/) {
      next;
    }

    # Here we should have found the line, so the sample name should be in the first column
    unless($line =~ /^([^\t]+)\t/) {
      $self->throw("Tried to parse the sample name from corresponding line in the csv file but failed (parsing regex issue?). Line:\n".$line);
    }

    $sample_name = $1;
    last;
  }
  close IN;

  return($sample_name);
}

1;
