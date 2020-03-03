=head1 LICENSE

 Copyright [2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses Stringtie2 to construct a gtf file from a sorted bam file

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2;

use warnings;
use strict;
use feature 'say';
use File::Basename;

use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Analysis::Runnable::Stringtie2;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    input_file_extension => '_Aligned.sortedByCoord.out.bam',
    num_threads => 1,
    delete_input_file => 0,
  }
}



=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for minimap2
 Returntype : None
 Exceptions : Throws if 'genome_file' does not exist
              Throws if 'input_file' does not exist

=cut

sub fetch_input {
  my ($self) = @_;

#  my $target_dba = $self->hrdb_get_dba($self->param('target_db'));

#  my $dna_dba;
#  if($self->param('use_genome_flatfile')) {
#    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
#      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
#    }
#    setup_fasta(
#                 -FASTA => $self->param_required('genome_file'),
#               );
#  } else {
#    $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
#    $target_dba->dnadb($dna_dba);
#  }

#  $self->hrdb_set_con($target_dba,'target_db');

  my $genome_file = $self->param_required('genome_file');
  unless(-e $genome_file) {
    $self->throw("Could not find the genome file. Path used:\n".$genome_file);
  }

  my $input_ids = [$self->param('iid')];
  unless(scalar(@$input_ids)) {
    $self->throw("Found no input ids, something is wrong");
  }

  my $program = $self->param('stringtie2_path');

  unless($program) {
    $program = "stringtie2";
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

    $analysis->logic_name($sample_name."_stringtie");
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Stringtie2->new(
         -analysis          => $analysis,
         -program           => $program,
         -input_file        => $input_file,
         -num_threads       => $self->param('num_threads'),
         -output_dir        => $output_dir,
#         -database_adaptor  => $target_dba,
         -delete_input_file => $self->param('delete_input_file'), # NB!! only set this when creating ranged files, not when using the original input file
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

  my $file_extension_to_remove = $self->param('input_file_extension');
  if($file_extension_to_remove) {
    $input_file =~ s/$file_extension_to_remove//;
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
