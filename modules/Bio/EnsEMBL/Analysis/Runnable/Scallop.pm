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

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Scallop

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::Scallop->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses Scallop to take a sorted bam file and construct
a gtf file of transcripts based on the intron data

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Scallop;

use warnings;
use strict;
use feature 'say';

use File::Spec::Functions;
use File::Basename;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

use parent ('Bio::EnsEMBL::Analysis::Runnable');


=head2 new

 Arg [DECOMPRESS]           : String as a command like 'gzip -c -'
 Arg [EXPECTED_ATTRIBUTES]  : String specify the attribute expected for the output, see STAR manual
 Description                : Creates a  object to align reads to a genome using STAR
 Returntype                 :
 Exceptions                 : Throws if WORKDIR does not exist
                              Throws if the genome has not been indexed

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($input_file, $output_dir, $delete_input_file) = rearrange([qw (INPUT_FILE OUTPUT_DIR DELETE_INPUT_FILE)],@args);
  $self->input_file($input_file);
  $self->output_dir($output_dir);
  $self->delete_input_file($delete_input_file);
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run scallop to construct transripts from short reads. The resulting output file will be stored in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $input_file    = $self->input_file;
  my $output_file = basename($input_file);
  $output_file =~ s/\.bam/\.gtf/;
  my $output_file_path = catfile($self->output_dir(),$output_file);

  if($self->delete_input_file) {
    $self->files_to_delete($input_file);
  }

  my $options = $self->options;

  # run stringtie
  my $scallop_command = $self->program." -i ".$input_file." -o ".$output_file_path;

  $self->warning("Command:\n".$scallop_command."\n");
  execute_with_wait($scallop_command);

  # check that the output GTF file does not have zero size
  if (-s $output_file_path) {
    $self->output([$output_file]);
  } else {
    $self->throw("Error running scallop\nThe output GTF file has zero size: $output_file_path\n");
  }
}


sub genome_index {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_genome_index} = $val;
  }

  return $self->{_genome_index};
}


sub input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_input_file} = $val;
  }

  return $self->{_input_file};
}


sub output_dir {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_output_dir} = $val;
  }
  return $self->{_output_dir};
}


sub database_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_database_adaptor} = $val;
  }

  return $self->{_database_adaptor};
}


sub delete_input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_delete_input_file} = $val;
  }

  return $self->{_delete_input_file};
}

1;
