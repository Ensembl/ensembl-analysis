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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2Merge

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2Merge->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses Stringtie2 to merge sample gtf files

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2Merge;

use warnings;
use strict;
use feature 'say';
use File::Spec::Functions;
use File::Basename;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    input_file_extension => '_Aligned.sortedByCoord.out.gtf',
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

  my $merged_sample_hash = {};
  $merged_sample_hash->{'merged'} = [];

  my @all_gtf_files = ();
  my $gtf_dirs = $self->param_required('input_gtf_dirs');
  foreach my $gtf_dir (@$gtf_dirs) {
    my @gtf_files = glob($gtf_dir."/*.gtf");
    unless(scalar(@gtf_files)) {
      $self->throw("Did not find any sample gtf files in the stringtie gtf dir. Path used:\n".$gtf_dir);
    }
    push(@all_gtf_files,@gtf_files);
  }

  unless(scalar(@all_gtf_files)) {
    $self->throw("No gtf dirs were provided");
  }

  foreach my $input_file (@all_gtf_files) {
    my $sample_name = $self->get_sample_name($input_file,$self->param_required('csv_summary_file'));
    unless($sample_name) {
      $sample_name = $self->get_sample_name($input_file,$self->param_required('csv_summary_file_genus'));
      unless($sample_name) {
        $self->throw("Checked both the normal and genus csv files and failed to match the input file to a line, so could not retrieve the sample name");
      }
    }

    unless($merged_sample_hash->{$sample_name}) {
      $merged_sample_hash->{$sample_name} = [];
    }

    push(@{$merged_sample_hash->{$sample_name}},$input_file);
    push(@{$merged_sample_hash->{'merged'}},$input_file);
  }

  $self->param('merged_sample_hash',$merged_sample_hash);
}


=head2 run

 Arg [1]    : None
 Description: Run will merge all the samples based on tissue name via stringtie merge
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $stringtie_path = $self->param_required('stringtie2_path');
  my $merged_sample_hash = $self->param('merged_sample_hash');
  my $stringtie_merge_dir = $self->param_required('stringtie_merge_dir');
  my $result = system('mkdir -p '.$stringtie_merge_dir);
  if($result) {
    $self->throw("Could not create stringtie merge output dir. Path used:\n".$stringtie_merge_dir);
  }

  my @tissues = keys(%{$merged_sample_hash});
  unless(scalar(@tissues)) {
    $self->throw("There were no keys in the merged sample hash. This is wrong");
  }

  foreach my $tissue (@tissues) {
    my @tissue_files = @{$merged_sample_hash->{$tissue}};
    unless(scalar(@tissue_files)) {
      $self->throw("There were no files for the tissue. This should not happen. Tissue affected: ".$tissue);
    }

    my $merged_file_path = catfile($stringtie_merge_dir,$tissue.".gtf");

    if(scalar(@tissue_files) == 1) {
      say_with_header("For ".$tissue." only found one file so will just copy it to merged dir");
      my $copy_command = "cp ".$tissue_files[0]." ".$merged_file_path;
      $result = system($copy_command);
      if($result) {
        $self->throw("Issue copying the gtf file to the merged dir. Commandline used:\n".$copy_command);
      }
      next;
    }

    my $tissue_file_list_path = $merged_file_path.".list";
    unless(open(OUT,">".$tissue_file_list_path)) {
      $self->throw("Failed to open output file to list gtf files for merge in. Path used:\n".$tissue_file_list_path);
    }

    foreach my $file_path (@tissue_files) {
      say OUT $file_path;
    }
    close(OUT) || $self->throw('Could not close '.$tissue_file_list_path);

    my $stringtie_merge_command = $stringtie_path." --merge -p ".$self->param('num_threads')." -o ".$merged_file_path." ".$tissue_file_list_path;
    $result = system($stringtie_merge_command);
    if($result) {
      $self->throw("Stringtie merge returned a non-zero exit code. Commandline used:\n".$stringtie_merge_command);
    }
  }
}


sub write_output {
  my ($self) = @_;
}


sub get_sample_name {
  my ($self,$input_file_path,$csv_file) = @_;

  my $input_file = basename($input_file_path);
  chomp($input_file);

  my $file_extension_to_remove = $self->param('input_file_extension');
  if($file_extension_to_remove) {
    $input_file =~ s/$file_extension_to_remove//;
  }

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
  close(IN) || $self->throw('Could not close '.$csv_file);

  return($sample_name);
}

1;
