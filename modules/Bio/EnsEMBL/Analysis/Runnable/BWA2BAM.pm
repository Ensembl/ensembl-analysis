=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::Samtools;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

use parent ('Bio::EnsEMBL::Analysis::Runnable::BWA');


=head2 new

 Arg [HEADER]    : String
 Arg [FASTQPAIR] : String
 Arg [SAMTOOLS]  : String
 Arg [MIN_MAPPED]: Integer
 Arg [MIN_PAIRED]: Integer
 Arg [BAM_PREFIX]: String
 Description     : Creates a new Bio::EnsEMBL::Analysis::Runnable::BWA2BAM object
 Returntype      : Bio::EnsEMBL::Analysis::Runnable::BWA2BAM
 Exceptions      : Throws if it cannot create a Bio::EnsEMBL::Analysis::Runnable::Samtools object

=cut
  
sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ($header, $fastqpair, $samtools, $min_mapped, $min_paired, $bam_prefix) = rearrange([qw(HEADER FASTQPAIR SAMTOOLS MIN_MAPPED MIN_PAIRED BAM_PREFIX)],@args);
  $self->fastqpair($fastqpair);
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
               It also store the % of mapped read andof paired reads if it is below the threshold
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
  my $header = $self->header;
  my $samtools = $self->samtools;
  my $readgroup;
  my $filename;
  my $outfile;
  my $pairfilename;
  my $total_reads = 0;
  # count how many reads we have in the fasta file to start with
  print "fastq pair is $fastqpair and fastq is $fastq\n";
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
      $readgroup =~ s/\s+(\w+:)/\\t$1/g;
    }
    print "using readgroup line $readgroup\n";
    close(HEAD) || $self->throw("Could not close $header\n");
  }

  # run bwa
  my $sai_fastq_files = "$outdir/$filename.sai $fastq";
  $self->files_to_delete("$outdir/$filename.sai");
  if ($fastqpair) {
    $sai_fastq_files = "$outdir/$filename.sai $outdir/$pairfilename.sai $fastq $fastqpair";
#    $self->files_to_delete("$outdir/$pairfilename.sai");
  }

  $command = join(' ', $program, $method, $readgroup, $self->genome, $sai_fastq_files, '|', $samtools->make_commandline('view', '-b -S', $outdir.'/'.$outfile.'_unsorted.bam', '-'));
 print "command is $command\n"; 
  execute_with_wait($command, 'Failed processing alignment: '.$command."\n$?\n");

  my $sorted_bam = $outdir.'/'.$outfile;
  # sort the bam
  $samtools->sort($sorted_bam, $sorted_bam.'_unsorted.bam');
  # index the bam
  $samtools->index($sorted_bam.'.bam');

 # check the reads with flagstat
  my $flagstat_results = $samtools->flagstat($sorted_bam.'.bam', 1);
  $self->throw('Got '.$flagstat_results->[0]." reads in flagstat rather than $total_reads in fastq - something went wrong\n")
    unless ($flagstat_results->[0] == $total_reads);
  my @output = ($sorted_bam.'.bam');
  if (($fastqpair and $self->min_paired > $flagstat_results->[2] and $self->min_mapped > $flagstat_results->[1])
     or (!$fastqpair and $self->min_mapped > $flagstat_results->[1])) {
    warning("Mappings for $sorted_bam.bam is bad:\n".
            'Min mapped '.$self->min_mapped.': '.$flagstat_results->[1]."\n".
            'Min paired '.$self->min_paired.': '.$flagstat_results->[2]."\n"
            );
    push(@output, $flagstat_results->[1], $flagstat_results->[2]);
  }
  $self->output(\@output);

  #if the BAM file has been sorted and indexed OK delete the original BAM file that was generated
  $self->files_to_delete("$outdir/$outfile"."_unsorted.bam");
  $self->delete_files;
}






#Containers
#=================================================================

=head2 samtools

 Arg [1]    : String (optional)
 Description: Getter/setter, creates an Bio::EnsEMBL::Analysis::Runnable::Samtools object
              when Arg[1] is defined and set the usage of threads.
 Returntype : Bio::EnsEMBL::Analysis::Runnable::Samtools
 Exceptions : Trows if it cannot create the Bio::EnsEMBL::Analysis::Runnable::Samtools object

=cut

sub samtools {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_samtools} = Bio::EnsEMBL::Analysis::Runnable::Samtools->new($value);
    $self->{_samtools}->use_threading($1) if ($self->analysis->parameters =~ /-use_threads\s+=>\s+(\d+)/);
  }
  
  if (exists($self->{'_samtools'})) {
    return $self->{'_samtools'};
  } else {
    return;
  }
}


=head2 fastqpair

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub fastqpair {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_fastqpair'} = $value;
  }
  
  if (exists($self->{'_fastqpair'})) {
    return $self->{'_fastqpair'};
  } else {
    return;
  }
}


=head2 header

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub header {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_header'} = $value;
  }
  
  if (exists($self->{'_header'})) {
    return $self->{'_header'};
  } else {
    return;
  }
}


=head2 min_mapped

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub min_mapped {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_min_mapped'} = $value;
  }

  if (exists($self->{'_min_mapped'})) {
    return $self->{'_min_mapped'};
  } else {
    return;
  }
}


=head2 min_paired

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub min_paired {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_min_paired'} = $value;
  }

  if (exists($self->{'_min_paired'})) {
    return $self->{'_min_paired'};
  } else {
    return;
  }
}


=head2 bam_prefix

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub bam_prefix {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_bam_prefix'} = $value;
  }

  if (exists($self->{'_bam_prefix'})) {
    return $self->{'_bam_prefix'};
  } else {
    return;
  }
}

1;
