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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessStarIntrons;

use warnings;
use strict;
use feature 'say';
use File::Basename;

use Bio::EnsEMBL::Analysis::Runnable::Star;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _branch_to_flow_to => 2,
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
  $self->process_introns($sam_file);
}


sub write_output {
  my ($self) = @_;

  my $bam_file = $self->param('_bam_file');
  $self->dataflow_output_id({intron_filename => $bam_file}, 1);
}


sub process_introns {
  my ($self,$sam_file) = @_;

  my $intron_sam_file = $sam_file;
  if($intron_sam_file =~ /_Aligned\.out\.sam/) {
    $intron_sam_file =~ s/_Aligned\.out\.sam/\.introns\.sam/;
  } elsif($intron_sam_file =~ /\.sam/) {
    $intron_sam_file =~ s/\.sam/\.introns\.sam/;
  } else {
    $self->throw("Unexpected sam filename pattern. Filename: ".$intron_sam_file);
  }

  my $unsorted_intron_bam_file = $intron_sam_file;
  $unsorted_intron_bam_file =~ s/\.sam/\.unsorted\.bam/;

  open(OUT,">".$intron_sam_file);

#  open(IN,"samtools view -S $sam_file|");
  unless(open(IN,$sam_file)) {
    $self->throw("Could not open the input sam file for reading. Filepath:\n".$sam_file);
  }

  while(<IN>) {
    my $line = $_;
    if($line =~ /^\@/) {
      print OUT $line;
    } else {
      my @ele = split("\t",$line);
      my $cigar = $ele[5];
      my $read_length = length($ele[9]);

#      my $ensembl_cigar = $self->convert_cigar($cigar);
#      unless($cigar eq $ensembl_cigar) {
#        unless($line =~ s/$cigar/$ensembl_cigar/) {
#          say "Issue with updating the cigar string in the sam line. Line processed:\n".$line;
#        }
#        say "CIGAR string updated: ".$cigar."->".$ensembl_cigar;
#      }

      if($self->process_read($cigar)) {
        print OUT $line;
      } else {
        say "Rejecting cigar: ".$cigar;
      }
    }
  }
  close IN;
  close OUT;

  my $unsorted_bam_command = 'samtools view -S -b '.$intron_sam_file.' > '.$unsorted_intron_bam_file;
  my $result = system($unsorted_bam_command);
  if($result) {
    $self->throw("Error running bam conversion. Commandline used:\n".$unsorted_bam_command);
  }

  system("rm ".$intron_sam_file);

  my $sorted_intron_bam_file = $unsorted_intron_bam_file;
  $sorted_intron_bam_file =~ s/\.unsorted//;
  my $sorted_bam_command = 'samtools sort '.$unsorted_intron_bam_file.' -o '.$sorted_intron_bam_file;
  $result = system($sorted_bam_command);
  if($result) {
    $self->throw("Error running sorting bam. Commandline used:\n".$sorted_bam_command);
  }
  system("rm ".$unsorted_intron_bam_file);

  my $index_bam_command = 'samtools index '.$sorted_intron_bam_file;
  $result = system($index_bam_command);
  if($result) {
    $self->throw("Error running indexing bam. Commandline used:\n".$index_bam_command);
  }

  $self->param('_bam_file',$sorted_intron_bam_file);
}


sub process_read {
  my ($self,$cigar) = @_;

  my $keep = 1;
  unless($cigar =~ /N/) {
    $keep = 0;
    return($keep);
  }

  my $min_alignment_anchor = 15;
  my $min_intron_size = 75;

  while($cigar =~ /(\d+)M(\d+)N(\d+)M/) {
    my $match_left = $1;
    my $intron = $2;
    my $match_right = $3;
    if($match_left < $min_alignment_anchor || $match_right < $min_alignment_anchor || $intron < $min_intron_size) {
      $keep = 0;
    }
    $cigar =~ s/$&//;
  }
  return($keep);
}


#sub convert_cigar {
#  my ($self,$cigar) = @_;

#  $cigar =~ s/^\d+S//;
#  $cigar =~ s/\d+S$//;

#  if($cigar =~ /^(\d+)S(\d+)M/) {
#    my $clip_count = $1;
#    my $match_count = $2;
#    my $new_match = ($clip_count + $match_count)."M";
#    $cigar =~ s/$&/$new_match/;
#  }

#  if($cigar =~ /(\d+)M(\d+)S$/) {
#    my $clip_count = $1;
#    my $match_count = $2;
#    my $new_match = ($clip_count + $match_count)."M";
#    $cigar =~ s/$&/$new_match/;
#  }

#  $cigar =~ s/N/I/g;
#  return($cigar);
#}


1;
