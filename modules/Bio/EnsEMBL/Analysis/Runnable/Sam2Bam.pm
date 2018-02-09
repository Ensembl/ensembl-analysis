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

Bio::EnsEMBL::Analysis::Runnable::Sam2Bam

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::Sam2Bam->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses samtools to convert a directory containing SAM
files into a single sorted indexed merged BAM file


=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Sam2Bam;

use warnings ;
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Analysis::Runnable::Samtools;

use parent ('Bio::EnsEMBL::Analysis::Runnable');


=head2 new

 Arg [HEADER]  : String
 Arg [SAMFILES]: String
 Arg [BAMFILE] : String
 Arg [GENOME]  : String
 Description   : Creates a new Bio::EnsEMBL::Analysis::Runnable::Sam2Bam object to convert SAM files
                 to a BAM file
 Returntype    : Bio::EnsEMBL::Analysis::Runnable::Sam2Bam
 Exceptions    : Throws if you have no files to work on
                 Throws if BAMFILE or GENOME is not set
              
=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ($header, $samfiles, $bamfile, $genome) = rearrange([qw(HEADER SAMFILES BAMFILE GENOME)],@args);
  $self->throw("You have no files to work on\n") unless (scalar(@$samfiles));
  $self->samfiles($samfiles);
  $self->headerfile($header);
  $self->throw("You must define an output file\n")  unless $bamfile;
  $self->bamfile($bamfile);
  $self->throw("You must define a genome file\n")  unless $genome;
  $self->genome($genome);
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
  my $bamfile = $self->bamfile;
  my $program = $self->program;
  my $samtools = Bio::EnsEMBL::Analysis::Runnable::Samtools->new(
                        -program => $program,
                        );
  my $count = 0;
  my @fails;
  # next make all the sam files into one big sam flie
  open (BAM ,">$bamfile.sam" ) or $self->throw("Cannot open sam file for merging $bamfile.sam");
  foreach my $file ( @{$self->samfiles} ) {
    my $line;
    open ( SAM ,"$file") or $self->throw("Cannot open file '$file'\n");
    my $line_count = 0;
    while (<SAM>) {
      # 1st file copy the header all the others just copy the data
      if ($_ =~ /^(\@\w+)/) {
        $line = $_;
      }
      else {
        chomp;
        print BAM "$_\n";
        $line_count++;
      }
    }
    close(SAM) || $self->throw("Failed closing '$file'\n");
    $count++;
    push @fails,$file  unless ( $line eq '@EOF' or $line_count == 0 );
    #last if $count >= 100;
  }
  close(BAM) || $self->throw("Failed closing sam file for merging $bamfile.sam");
  print "Merged $count files\n";
  if ( scalar(@fails) > 0 ) {
    $self->throw(join("\n", 'The following sam files failed to complete, you need to run them again', @fails));
  }
  my $fh;
  # now call sam tools to do the conversion.
  # is the genome file indexed?
  # might want to check if it's already indexed first
  my $command = "$program faidx " . $self->genome ;
  unless ( -e (  $self->genome.".fai" ) ) {
    print "Indexing genome file\n";
    print STDERR "$command \n";
    system("$command 2> /tmp/sam2bam_index.err");
    open  ( $fh,"/tmp/sam2bam_index.err" ) or die ("Cannot find STDERR from fasta indexing\n");
    # write output
    while(<$fh>){
      print STDERR "INDEX $_";
      $self->throw('Samtools failed to index '.$self->genome) if ($_ =~ /fail/ or $_ =~ /abort/ or $_ =~ /truncated/ )
     }
     close($fh) || $self->throw("Cannot close STDERR from fasta indexing");
     $self->files_to_delete("/tmp/sam2bam/index.err");
  }

  $command = "$program view -b -h -S -T " . $self->genome ." $bamfile.sam >  $bamfile" . "_unsorted.bam ";
  print STDERR "$command \n";
  system("$command 2> /tmp/sam2bam_view.err");
  open  ( $fh,"/tmp/sam2bam_view.err" ) or die ("Cannot find STDERR from samtools view\n");
  # write output
  while(<$fh>){
    print STDERR "IMPORT $_";
    $self->throw('Samtools failed to created an unsorted bam file '.$bamfile.'_unsorted.bam') if ($_ =~ /fail/ or $_ =~ /abort/ or $_ =~ /truncated/ )
  }
  close($fh) || $self->throw("Cannot close STDERR from samtools view");
  $self->files_to_delete($bamfile.'.sam');
  $self->files_to_delete("/tmp/sam2bam_view.err");

  # add readgroup info if there is any
  if ( $self->headerfile ) {
    # dump out the sequence header
     $command = "$program view  -H   $bamfile" . "_unsorted.bam > $bamfile.header";
     print STDERR "$command \n";
     system("$command");
     $command = " cat " . $self->headerfile ." >>  $bamfile.header  ";
     print STDERR "$command \n";
     system("$command");
     $command = "$program reheader   $bamfile.header $bamfile" . "_unsorted.bam > $bamfile.fixed_header.bam";
     print STDERR "$command \n";
     system("$command");
     $self->files_to_delete("/tmp/sam2bam_reheader.err");
     $command = "mv  $bamfile.fixed_header.bam $bamfile"."_unsorted.bam ";
     print STDERR "$command \n";
     system("$command");
  }

  $samtools->sort($bamfile, $bamfile.'_unsorted.bam');
  $samtools->index($bamfile.'.bam');

  $self->files_to_delete($bamfile.'_unsorted.bam');
  $self->delete_files();
}




#Containers
#=================================================================


=head2 headerfile

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : Throws if file doesn't exist

=cut

sub headerfile {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->throw("Cannot find read group header file $value\n")
      unless (-e $value);
    $self->{'_headerfile'} = $value;
  }

  if (exists($self->{'_headerfile'})) {
    return $self->{'_headerfile'};
  } else {
    return undef;
  }
}



=head2 samfiles

 Arg [1]    : (optional) Arrayref of String
 Description: Getter/setter
 Returntype : Arrayref of String
 Exceptions : None

=cut

sub samfiles {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_samfiles'} = $value;
  }

  if (exists($self->{'_samfiles'})) {
    return $self->{'_samfiles'};
  } else {
    return undef;
  }
}


=head2 bamfile

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub bamfile {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_bamfile'} = $value;
  }

  if (exists($self->{'_bamfile'})) {
    return $self->{'_bamfile'};
  } else {
    return undef;
  }
}


=head2 genome

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : Throws if the file doesn't exists

=cut

sub genome {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->throw($value.' does not exist!') unless (-e $value);
    $self->{'_genome'} = $value;
  }

  if (exists($self->{'_genome'})) {
    return $self->{'_genome'};
  } else {
    return undef;
  }
}

1;
