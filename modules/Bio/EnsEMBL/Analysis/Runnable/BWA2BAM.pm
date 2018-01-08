=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::BWA2BAM

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::BWA2BAM->new();

 $runnable->run;
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module uses BWA to align fastq to a genomic sequence


=cut

package Bio::EnsEMBL::Analysis::Runnable::BWA2BAM;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable::BWA;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
$| = 1;
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BWA);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ($header, $fastqpair, $samtools, $min_mapped, $min_paired, $bam_prefix) = rearrange([qw(HEADER FASTQPAIR SAMTOOLS MIN_MAPPED MIN_PAIRED BAM_PREFIX)],@args);
  $self->fastqpair($fastqpair);
  $self->throw("You must define a path to samtools cannot find $samtools\n")  
    unless $samtools && -e $samtools;
  $self->samtools($samtools);
  $self->header($header);
  $self->min_mapped($min_mapped);
  $self->min_paired($min_paired);
  $self->bam_prefix($bam_prefix);

  return $self;
}

=head2 run

  Args       : none
  Description: Create the BAM files using BWA to first join the pairs
               If we have paired reads. Then sort using samtools.
               It deletes the temporary files.
               It stores the name of the resulting BAM file in $self->output
  Returntype : none

=cut 

sub run {
  my ($self) = @_;
  # get a list of files to use

  my $fastq = $self->fastq;
  my $fastqpair = $self->fastqpair;
  my $outdir = $self->outdir;
  my $program = $self->program;
  my $method = $self->options;
  my $samtools = $self->samtools;
  my $header = $self->header;
  my $readgroup;
  my $filename;
  my $outfile;
  my $pairfilename;
  my $total_reads = 0;
  # count how many reads we have in the fasta file to start with
  my $command;
  my $fh;
  if (-B $fastq) {
    $command = "gunzip -c $fastq | wc -l";
  }
  else {
    $command = "wc -l $fastq";
  }
  print STDERR "$command\n";
  open  ( $fh,"$command 2>&1 |" ) ||
    $self->throw("Error counting reads");
  while (<$fh>){
    chomp;
    if ( $_ =~ /^\s*(\d+)/ ) {
        if ( $fastqpair ) {
        $total_reads = $1 / 2;
  } else {
    $total_reads = $1 / 4;
  }
  print STDERR "Starting with $total_reads reads\n";
    }
  }
  close($fh) || $self->throw("Failed counting reads");

  unless ($total_reads) {
    $self->throw("unable to count the reads in the fastq file $fastq\n");
  }
  my @tmp = split(/\//,$fastq);
  $filename = pop @tmp;
  $outfile = $self->bam_prefix || $filename if ($self->bam_prefix);
  print "Filename $outfile\n";
  if ( $fastqpair ) {
    my @tmp = split(/\//,$fastqpair);
    $pairfilename = pop @tmp;
    print "Filename $pairfilename\n";
  }
  # add the header line if there is one
  if ( $header ) {
    open (HEAD,"$header") or $self->throw("No read group file found $header\n");
    $readgroup = "-r ";
    while (<HEAD>){
      chomp;
      $readgroup .= "\"$_\"";
    }
    print "using readgroup line $readgroup\n";
    close(HEAD) || $self->throw("Could not close $header\n");
  }

  # run bwa
  my $sai_fastq_files = "$outdir/$filename.sai $fastq";
  $sai_fastq_files = "$outdir/$filename.sai $outdir/$pairfilename.sai $fastq $fastqpair" if ($fastqpair);

  $command = join(' ', $program, $method, $readgroup, $self->genome, $sai_fastq_files, '|', $samtools, 'view - -b -S -o', "$outdir/$outfile.bam");
  
  print STDERR "Command: $command\n";
  open  ( $fh,"$command 2>&1 |" ) ||
    $self->throw("Error processing alignment $@\n");
    while (<$fh>){
  chomp;
  print STDERR "VIEW: $_\n";
  if ( $_ =~ /truncated file/ ) {
    $self->throw("Error converting file\n");
  }
    }
  close($fh) || $self->throw("Failed processing alignment");

  my $sorted_bam = $outdir.'/'.$outfile.'_sorted';
  # sort the bam
  $command = "$samtools sort $outdir/$outfile.bam $sorted_bam";
  print STDERR "Sort: $command\n";
    open  ( $fh,"$command 2>&1 |" ) ||
      $self->throw("Error sorting bam $@\n");
    while (<$fh>){
      chomp;
      print STDERR "SORT: $_\n";
    }
  close($fh) || $self->throw("Failed sorting bam");
  # index the bam
  $command = "$samtools index $sorted_bam.bam";
  print STDERR "Index: $command\n";
    open  ( $fh,"$command 2>&1 |" ) ||
      $self->throw("Error indexing bam $@\n");
    while (<$fh>){
      chomp;
      print STDERR "INDEX: $_\n";
    }
  close($fh) || $self->throw("Failed indexing bam");

 # check the reads with flagstat
  $command = "$samtools flagstat $sorted_bam.bam";
  print STDERR "Got $total_reads to check\nCheck: $command\n";
    open  ( $fh,"$command 2>&1 |" ) ||
      $self->throw("Error checking alignment $@\n");
    while (<$fh>){
      print STDERR "$_";
      if ( $_ =~ /EOF marker is absent/ or /invalid BAM binary header/ ) {
	$self->throw("EOF marker absent, something went wrong\n");
      }
      if ( $_ =~ /(\d+) \+ \d+ in total \(/ ) {
	$self->throw("Got $1 reads in flagstat  rather than $total_reads in fastq - something went wrong\n")
	  unless ( $total_reads - $1 ) <=1 ;
      }
      elsif (/^\s*\d+.*mapped\s+\(([0-9.]+)/) {
          warning("$sorted_bam.bam is below the threshold of ".$self->min_mapped.": $1") if ($self->min_mapped > $1);
      }
      elsif (/\s*\d+.*properly paired \(([0-9.]+)/) {
          warning("$sorted_bam.bam is below the threshold of ".$self->min_paired.": $1") if ($self->min_paired > $1);
      }
    }
  close($fh) || $self->throw("Failed checking alignment");

  #if the BAM file has been sorted and indexed OK delete the original BAM file that was generated
  $command = "rm $outdir/$outfile".".bam";
  print STDERR "Delete: $command\n";
  open  ( $fh,"$command 2>&1 |" ) || $self->throw("Error deleting unsorted bam $@\n");
  while (<$fh>)
  {
      chomp;
      print STDERR "DELETE: $_\n";
  }
  close($fh) || $self->throw("Failed deleting bam");
  $self->output([$sorted_bam.'.bam']);
}






#Containers
#=================================================================

sub samtools {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_samtools'} = $value;
  }
  
  if (exists($self->{'_samtools'})) {
    return $self->{'_samtools'};
  } else {
    return undef;
  }
}

sub fastqpair {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_fastqpair'} = $value;
  }
  
  if (exists($self->{'_fastqpair'})) {
    return $self->{'_fastqpair'};
  } else {
    return undef;
  }
}


sub header {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_header'} = $value;
  }
  
  if (exists($self->{'_header'})) {
    return $self->{'_header'};
  } else {
    return undef;
  }
}

sub min_mapped {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_min_mapped'} = $value;
  }

  if (exists($self->{'_min_mapped'})) {
    return $self->{'_min_mapped'};
  } else {
    return undef;
  }
}

sub min_paired {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_min_paired'} = $value;
  }

  if (exists($self->{'_min_paired'})) {
    return $self->{'_min_paired'};
  } else {
    return undef;
  }
}

sub bam_prefix {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_bam_prefix'} = $value;
  }

  if (exists($self->{'_bam_prefix'})) {
    return $self->{'_bam_prefix'};
  } else {
    return undef;
  }
}

1;
