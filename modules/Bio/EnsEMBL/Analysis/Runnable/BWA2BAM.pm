=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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
  my ($header,$method,$samtools) = rearrange([qw(HEADER METHOD SAMTOOLS)],@args);
  $self->throw("You must defne an alignment processing method not $method\n")  unless $method ;
  $self->method($method);
  $self->throw("You must define a path to samtools cannot find $samtools\n")  
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
  my $fastqpair = $self->fastqpair;
  my $options = $self->options;
  my $outdir = $self->outdir;
  my $program = $self->program;
  my $method = $self->method;
  my $samtools = $self->samtools;
  my $header = $self->header;
  my $readgroup;
  my $filename;
  my $outfile;
  my $pairfilename;
  my $total_reads = 0;
  # count how many reads we have in the fasta file to start with
  my $command = "wc -l $fastq";
  print STDERR "$command\n";
  eval  {
    open  ( my $fh,"$command 2>&1 |" ) || 
      $self->throw("Error counting reads");
    while (<$fh>){
      chomp;
      if ( $_ =~ /(\d+) $fastq/ ) {
      	if ( $fastqpair ) {
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
  
  unless ($total_reads) {
    $self->throw("unable to count the reads in the fastq file $fastq\n");
  }
  my @tmp = split(/\//,$fastq);
  $filename = pop @tmp;
  print "Filename $filename\n";
  $outfile = $filename;
  #($outfile) = $filename =~ /^(.+)\.[^.]+$/;
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
  }

  # run bwa
  unless ( -e (  $self->genome.".ann" ) ) {
    $self->throw("Genome file must be indexed \ntry " . $self->program . " index " . $self->genome ."\n"); 
  }
  
  $command = "$program $method $readgroup " . $self->genome .
    " $outdir/$filename.sai $fastq  \| $samtools  view - -b -S -o $outdir/$outfile.bam ";
  if ( $fastqpair ) {
    $command = "$program $method $readgroup " . $self->genome .
      " $outdir/$filename.sai $outdir/$pairfilename.sai $fastq $fastqpair \| $samtools  view - -b -S -o $outdir/$outfile.bam ";
  }
  
  print STDERR "Command: $command\n";
  eval {
    open  ( my $fh,"$command 2>&1 |" ) || 
      $self->throw("Error processing alignment $@\n");
      while (<$fh>){
	chomp;
	print STDERR "VIEW: $_\n";
	if ( $_ =~ /truncated file/ ) {
	  $self->throw("Error converting file\n");
	}
      }
   }; if($@){
    $self->throw("Error processing alignment \n$@\n");
  }

  # sort the bam
  $command = "$samtools  sort $outdir/$outfile.bam $outdir/$outfile"."_sorted";
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
  $command = "$samtools  index $outdir/$outfile"."_sorted.bam";
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
  $command = "$samtools  flagstat $outdir/$outfile"."_sorted.bam";
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


sub method {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_method'} = $value;
  }
  
  if (exists($self->{'_method'})) {
    return $self->{'_method'};
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



