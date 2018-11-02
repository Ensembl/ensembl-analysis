=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSequences

=head1 SYNOPSIS


=head1 DESCRIPTION

Base module to load data into Hive custom tables. Modules should override
create_row_data which will be called before loading data into the table.
the accession field will be return in an arrayref on branch 2

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::SplitFastQFiles;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;
use File::Basename;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

 Arg [1]    : None
 Description: Check that sequence_file and table_name exists. Load the module
              needed to parse the file
 Returntype :
 Exceptions :

=cut

sub fetch_input {
  my $self = shift;

  $self->param_required('max_reads_per_split');
  $self->param_required('max_total_reads');
  $self->param_required('rnaseq_summary_file');
  $self->param_required('fastq_dir');
  $self->param_required('iid');
}


sub run {
  my ($self) = @_;

  my $fastq_file_name = $self->param('iid');
  my $max_lines_per_split = $self->param('max_reads_per_split') * 4;

  my $fastq_file = catfile($self->param('fastq_dir'),$fastq_file_name);
  unless(-e $fastq_file) {
    $self->throw("Could not find fastq file on path. Path:\n".$fastq_file);
  }

  say "Counting the reads in ".$fastq_file_name;
  my $count_cmd = 'zcat '.$fastq_file.' | echo $((`wc -l`/4))';
  my $count_reads = `$count_cmd`;
  chomp($count_reads);
  unless($count_reads) {
    $self->throw("Failed to count reads in the following file:\n".$fastq_file);
  }
  say "Found ".$count_reads." reads";

  my $csv_file = $self->param('rnaseq_summary_file');
  unless(-e $csv_file) {
    $self->throw("Could not find the csv file. Path used:\n".$csv_file);
  }

  my $grep_cmd = "grep '".$fastq_file_name."' ".$csv_file;
  my $csv_entry = `$grep_cmd`;
  chomp($csv_entry);
  unless($csv_entry) {
    $self->throw("Could not find a corresponding line for ".$fastq_file_name." in the csv file:\n".$csv_file);
  }


  if($count_reads <= $self->param('max_reads_per_split')) {
    $self->for_csv_output($csv_entry);
    return;
  }

  my $decompressed_name = $fastq_file_name;
  unless($decompressed_name =~ s/\.gz$//) {
    $self->throw("Did not remove a .gz extension when trying to generate the decompressed version of the name. Name used:\n".$decompressed_name);
  }

  say "Preparing to split the fastq file";
  my $pairing_regex = '(\S+)(\_\d\.\S+)';
  unless($fastq_file_name =~ /$pairing_regex/) {
    $self->throw('Failed name did not match the regex. Regex used: '.$pairing_regex.' filename: '.$fastq_file_name);
  }

  my $sample_name = $1;
  my $remainder = $2;
  my $split_file_paths = $self->split_fastq($fastq_file_name,$self->param('fastq_dir'),$max_lines_per_split,$sample_name,$remainder);
  say "Fastq file has been split";

  my $decompressed_glob_path = catfile($self->param('fastq_dir'),$sample_name."_*_split");

  my @split_files = glob($decompressed_glob_path);
  unless(scalar(@split_files)) {
    $self->throw("Failed to glob split files in input dir. Glob path used:\n".$decompressed_glob_path);
  }

  foreach my $split_file (@$split_file_paths) {
    say "Processing split file: ".$split_file;
    my $split_file_name = basename($split_file);
    my $split_pairing_regex = '(\S+\_\d+)\_\d\.\S+';
    unless($split_file_name =~ /$split_pairing_regex/) {
      $self->throw('Failed name did not match the split regex. Regex used: '.$split_pairing_regex.' filename: '.$split_file_name);
    }
    my $split_sample_name = $1;

    my $new_csv_entry = $csv_entry;
    say "FERGAL SAMPLE: ".$sample_name;
    say "FERGAL SPLIT SAMPLE: ".$split_sample_name;
    unless($new_csv_entry =~ s/\t$sample_name\t/\t$split_sample_name\t/) {
      $self->throw("Failed to update the original sample name in the csv (".$fastq_file_name.") to the split sample name (".$split_file_name.")");
    }
    say "FERGAL NEW: ".$new_csv_entry;
    unless($new_csv_entry =~ s/$fastq_file_name/\t$split_file_name\t/) {
      $self->throw("Failed to update the original file name in the csv (".$fastq_file_name.") to the split file name (".$split_file_name.")");
    }

    $new_csv_entry =~ s/\t+/\t/g;
    $self->for_csv_output($new_csv_entry);
  }
}


sub write_output {
  my ($self) = @_;

  my $output_file = catfile($self->param('fastq_dir'),$self->param('iid')."_new.csv");
  unless(open(OUT,">".$output_file)) {
    $self->throw("Failed to open the new csv file for writing. Path used:\n".$output_file);
  }

  foreach my $line (@{$self->for_csv_output()}) {
    say OUT $line;
  }
  close OUT;
}


sub split_fastq {
  my ($self,$file_name,$dir,$max_lines_per_split,$sample_name,$remainder) = @_;

  my $split_fastq_files = [];
  my $count = 0;
  my $increment = 0;
  my $file_path = catfile($dir,$file_name);
  say "FERGAL 1";
  say "FERGAL 1 PATH: ".$file_path;

  if($file_name =~ /\.gz/) {
    unless(open(IN_FASTQ, "gunzip -c ".$file_path." |")) {
      $self->throw("Couldn't open the fastq file to read. Path used: ".$file_path);
    }
  } else {
    unless(open(IN_FASTQ, $file_path)) {
      $self->throw("Couldn't open the fastq file to read. Path used: ".$file_path);
    }
  }

  my $new_output_file = $sample_name.'_'.$increment.$remainder.'_split';
  my $max_total_reads = $self->param('max_total_reads');
  open(OUT_FASTQ, '>'.catfile($dir, $new_output_file)) || $self->throw('Could not open output file for writing. Path used: '.catfile($dir, $new_output_file));
  push(@{$split_fastq_files}, $new_output_file);
  while(my $line = <IN_FASTQ>) {
    print OUT_FASTQ $line;
    ++$count;
    if($count == $max_total_reads) {
      # If it has gone over the read limit at this point then we just cut here. Not a great system in general, sub sampling would be better
      # but should work fine in most cases
      $self->warning("The read count from the fastq has passed the max allowed number of reads. Will not read further.\nMax total reads: $max_total_reads");
      last;
    }
    elsif($count && $count % $max_lines_per_split == 0) {
      close(OUT_FASTQ) || $self->throw("Could not close the fastq file ".catfile($dir, $sample_name.'_'.$increment.$remainder.'_split'));
      ++$increment;
      my $new_output_file = $sample_name.'_'.$increment.$remainder.'_split';
      open(OUT_FASTQ, '>'.catfile($dir, $new_output_file)) || $self->throw('Could not open output file for writing. Path used: '.catfile($dir, $new_output_file));
      push(@{$split_fastq_files},$new_output_file);
    }
  }
  close IN_FASTQ;
  close(OUT_FASTQ) || $self->throw("Could not close the fastq file");

  if(-z catfile($dir, $split_fastq_files->[-1])) {
    my $empty_file = pop(@$split_fastq_files);
    unlink catfile($dir, $empty_file) || $self->throw('Could not delete empty file '.$empty_file);
  }
  return($split_fastq_files);
}


sub for_csv_output {
  my ($self,$line) = @_;
  unless($self->param_is_defined('_for_csv_output')) {
    $self->param('_for_csv_output',[]);
  }

  if($line) {
    push(@{$self->param('_for_csv_output')},$line);
  }

  return($self->param('_for_csv_output'));
}

1;
