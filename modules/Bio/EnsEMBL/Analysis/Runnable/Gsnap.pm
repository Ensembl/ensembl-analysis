=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Gsnap

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::Gsnap->new();

 $runnable->run;
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module uses Gsnap to align fastq to a genomic sequence

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Gsnap;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
$| = 1;
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ($options, $fastq,$paired,$outdir,$indir,$genome,$genome_name,$samtools,$header) = 
    rearrange([qw (OPTIONS FASTQ PAIRED OUTDIR INDIR GENOME GENOMENAME SAMTOOLS HEADER)],@args);
  $self->throw("You must defne a fastq file(s)\n")  unless $fastq ;
  $self->fastq($fastq);
  $self->paired($paired);
  $self->options($options);
  $self->throw("You must defne an output dir\n")  unless $outdir;
  $self->outdir($outdir);
  $self->throw("You must defne an input dir\n")  unless $indir;
  $self->indir($indir);
  $self->throw("You must defne a genome dir\n")  unless $genome;
  $self->genome($genome);
  $self->throw("You must defne a genome name \n")  unless $genome_name;
  $self->genome_name($genome_name);
  $self->throw("You must defne a path to samtools cannot find $samtools\n")  
    unless $samtools && -e $samtools;
  $self->samtools($samtools);
  $self->header($header);
  return $self;
}

=head2 run

  Args       : none
  Description: Merges Sam files defined in the config using Samtools
  Returntype : none

=cut 

sub run {
  my ($self) = @_;
  # get a list of files to use
  my @files;

  my $fastq = $self->fastq;
  my $options = $self->options;
  # default behaviour
  $options = " --nthreads 1 --novelsplicing=1 --format=sam " unless $options;
  my $outdir = $self->outdir;
  my $program = $self->program;
  my $samtools = $self->samtools;
  my $header = $self->header;
  my $readgroup;
  my $filename;
  my $pair_file;
  my $fh;
  my $string;
  my $success = 0;  
  #  count how many reads we have in the fasta file to start with
  my $total_reads = 0;
  my $command = "wc -l $fastq";
  print STDERR "$command\n";
  eval  {
    open  ( my $fh,"$command 2>&1 |" ) || 
      $self->throw("Error counting reads");
    while (<$fh>){
      chomp;
      if ( $_ =~ /(\d+) $fastq/ ) {
      	if ( $pair_file ) {
   	  $total_reads = $1 / 2;
	} else {
	  $total_reads = $1 / 4;
	}
	print STDERR "Starting with $total_reads reads\n";
      }
    }
  }; if($@){
    $self->throw("Error processing alignment \n$@\n");
  }
  
  my @tmp = split(/\//,$fastq);
  $filename = pop @tmp;
  print "Filename $filename\n";
  if ( $self->paired ) {
    # if it is a paired lane there will be 2 filenames separated by a colon
    # need to replace the : with a space
    if ( $filename =~ /(\S+):(\S+)/ ) {
      $filename = $1;
      $pair_file = $2;
      $self->throw("Fastq file  $pair_file not found\n")
	unless ( -e $self->indir . "/$pair_file" );
    } else {
      $self->throw("Have paired flag set but cannot find paired file names in $filename\n");
    }
  } 
  $self->throw("Fastq file  $filename not found\n")
    unless ( -e $self->indir . "/$filename" );
  
  # run gsnap
  my @indexcommand = split(/\//,$self->program);
  pop @indexcommand;
  my $indexcommand = join("/",@indexcommand);
  print "Remember to index your genome file something like:\n$indexcommand/gmap_build -d genome -D " . 
    $self->genome ." -d " . 
      $self->genome_name ."\n";
  
  
  # add the header line if there is one
  if ( $header ) {
    open (HEAD,"$header") or $self->throw("No read group file found $header\n");
    my %kvp;
    my @tags;
    while (<HEAD>){
      chomp;
      foreach my $tag ( split (/\t/,$_) ){
	if ( $tag =~ /(\S+):(.*)$/ ) {
	  $kvp{$1} = $2;
	}
      }
    }
    # now put together the read group options
    $readgroup .= ' --read-group-id "'. $kvp{'ID'} .'" ' 
      if $kvp{'ID'};
    $readgroup .= ' --read-group-name "'. $kvp{'SM'} .'" ' 
      if $kvp{'SM'};
    $readgroup .= ' --read-group-library "'. $kvp{'LB'} .'" ' 
      if $kvp{'LB'};
    $readgroup .= ' --read-group-platform "'. $kvp{'PL'} .'" ' 
      if $kvp{'PL'};
  }
  
  
  $command = "$program  $options $readgroup -D " . $self->genome .
    " -d " . $self->genome_name . " " . $self->indir . "/$filename " ;
  if ( $pair_file ) {
    $command .=  $self->indir . "/$pair_file " ;
  }
  $command .= " \| $samtools  view - -b -S -o $outdir/$filename.bam ";
  
  print STDERR "Command: $command\n";
 # exit;
  eval {
     open  ( my $fh,"$command 2>&1 |" ) || 
      $self->throw("Error running gsnap $@\n");  
    while (<$fh>){
      chomp;
      print STDERR "GSNAP: $_\n";
     }
   }; if($@){
     $self->throw("Error aligning $filename \n$@\n");
   } 
  # sort the bam
  $command = "$samtools  sort $outdir/$filename.bam $outdir/$filename"."_sorted";
  print STDERR "Sort: $command\n";
  eval {
    open  ( my $fh,"$command 2>&1 |" ) || 
      $self->throw("Error sorting bam $@\n");  
    while (<$fh>){
      chomp;
      print STDERR "SORT: $_\n";
    }
  }; if($@){
    $self->throw("Error sorting bam \n$@\n");
  }  
  # index the bam
  $command = "$samtools  index $outdir/$filename"."_sorted.bam";
  print STDERR "Index: $command\n";
  eval {
    open  ( my $fh,"$command 2>&1 |" ) || 
      $self->throw("Error indexing bam $@\n");
    while (<$fh>){
      chomp;
      print STDERR "INDEX: $_\n";
    }
  }; if($@){
    $self->throw("Error indexing bam \n$@\n");
  } 
  
  # check the reads with flagstat
  $command = "$samtools  flagstat $outdir/$filename"."_sorted.bam";
  print STDERR "Got $total_reads to check/nCheck: $command\n";
  eval {
    open  ( my $fh,"$command 2>&1 |" ) || 
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
    }
  }; if($@){
    $self->throw("Error checking alignment \n$@\n");
  }
  
  
}




#Containers
#=================================================================

sub fastq {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_fastq'} = $value;
  }
  
  if (exists($self->{'_fastq'})) {
    return $self->{'_fastq'};
  } else {
    return undef;
  }
}

sub paired {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_paired'} = $value;
  }
  
  if (exists($self->{'_paired'})) {
    return $self->{'_paired'};
  } else {
    return undef;
  }
}

sub options {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_options'} = $value;
  }
  
  if (exists($self->{'_options'})) {
    return $self->{'_options'};
  } else {
    return undef;
  }
}

sub outdir {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_outdir'} = $value;
  }
  
  if (exists($self->{'_outdir'})) {
    return $self->{'_outdir'};
  } else {
    return undef;
  }
}

sub indir {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_indir'} = $value;
  }
  
  if (exists($self->{'_indir'})) {
    return $self->{'_indir'};
  } else {
    return undef;
  }
}


sub genome {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_genome'} = $value;
  }
  
  if (exists($self->{'_genome'})) {
    return $self->{'_genome'};
  } else {
    return undef;
  }
}

sub genome_name { 
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_genome_name'} = $value;
  }
  
  if (exists($self->{'_genome_name'})) {
    return $self->{'_genome_name'};
  } else {
    return undef;
  }
}

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
