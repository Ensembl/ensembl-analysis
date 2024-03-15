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

  Bio::EnsEMBL::Analysis::Runnable::RNAfold

=head1 SYNOPSIS
 
 my $runnable = Bio::EnsEMBL::Analysis::Runnable::RNAfold->new
    (
     -analysis  => $analysis,
     -sequence  => $seq,
     -structure => $str,
    );
  $self->runnable($runnable);

=head1 DESCRIPTION

Runnble to wrap RNAfold, part of the Vienna RNA package.
Takes a Bio::Seq object and an optional structure string as parameters.
If a structure is provided it uses the structure to constrain the folding
prediction (RNAFold -C).
The resulting structure string is run-length encoded.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::RNAFold;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

my $verbose = "";

=head2 new

  Title      : new
  Usage      :   my $runnable = Bio::EnsEMBL::Analysis::Runnable::RNAfold->new
             :    (
             :     -sequence  => $seq,
             :     -structure => $structure,
             :     -analysis  => $self->analysis,
             :    );
  Function   : Instantiates new RNAfold runnable
  Returns    : Bio::EnsEMBL::Analysis::Runnable::RNAFold
  Exceptions : Thows if no sequence is provided
  Args       : sequence (Bio::PrimarySeq)
             : analysis (Bio::EnsEMBL::Analysis)
             : structure (String)

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($sequence,$structure) = rearrange(['SEQUENCE','STRUCTURE'], @args);
  $self->throw("RNAFold: Cannot work without a sequence object\n") unless $sequence;
  $self->sequence($sequence);
  $self->structure_constraint($structure) if $structure;
  return $self;
}

=head2 new

  Title      : run
  Usage      : my $runnable->run;
  Function   : Run method.
  Returns    : None
  Exceptions : None
  Args       : None

=cut

sub run{
  my ($self) = @_;
  print STDERR "Run analysis\n" if $verbose;
  my $filename = $self->write_seq($self->sequence);
  $self->RNAfold($self->sequence,$filename);
  my $structure = $self->parse;
  $self->structure($structure);
  $self->encoded_str($self->encode_str($structure));
#  print STDERR "delete temp files\n"  if $verbose;
  $self->delete_files;
}

=head2 RNAfold

  Title      : RNAfold
  Usage      : $self->RNAfold($seq,$filename);
  Function   : Wrapper for RNAfold
  Returns    : none
  Exceptions : Throws if RNAfold fails to run
  Args       : sequence (Bio::PrimarySeq)
             : filename (String)

=cut

sub RNAfold{
  my ($self,$seq,$filename)=@_;
  my $command  = "RNAfold ";
  my $options = "";
  if ($self->structure_constraint){
    $options = " -C ";
  }
  my $results_file = $self->create_filename("RNAfold","txt");
  $self->files_to_delete($results_file);
  # delete the postcript file that RNAfold generates
  $self->files_to_delete(catfile($self->workdir, $seq->display_id."_ss.ps"));
  $self->resultsfile($results_file);
  $command .= "$options < $filename  2>&1 > ".$results_file;
  print STDERR "Running RNAfold ".$command."\n";
  chdir($self->workdir); # RNAfold create files in the current directory
  open(my $fh, "$command |") || 
    $self->throw("Error opening RNAfold cmd <$command>");
  # this loop reads the STDERR from the blast command
  while(<$fh>){
    if(/FATAL:(.+)/){
      my $match = $1;
      $self->throw("miRNA: RNAfold failed to run: $match$@\n");
    }
  }
  close $fh;
  return 1;
}

=head2 parse

  Title      : parse
  Usage      : my $structure = $self->parse;
  Function   : Parses the results of RNAfold
  Returns    : structure (String)
  Exceptions : Throws if results file cannot be opened or closed
  Args       : None

=cut

sub parse{
  my ($self)=@_;  
  my $results = $self->resultsfile;
  my $structure;
  my $score;
  open(RNAFOLD, $results) or $self->throw("FAILED to open ".$results.
					  " miRNA:parse_results\n$@\n");
 LINE: while(<RNAFOLD>){
    chomp;
    if ($_ =~ /([().]+)\s\(\s*(-*\d+.\d+)\)$/){
      $structure = $1;
      $score = $2;
    }
  }
  close(RNAFOLD) or $self->throw("FAILED to close ".$results.
				 " miRNA:parse_results\n$@\n");
  $self->score($score);
  if ($structure){
    return $structure;
  } else {
    return 0;
  }
}

=head2 write_seq

  Title      : write_seq
  Usage      : my $filename = $self->write_seq($daf);
  Function   : Writes the dna sequence file 
  Returns    : filename (String)
  Exceptions : Throws if it cannot write to the file
  Args       : Bio::PrimarySeq

=cut

sub write_seq{
  my ($self,$seq)=@_;
  my $filename = $self->create_filename("miRNA","seq");
  # have to write file so the sequence is all on a single line 
  # cos thats the way RNAfold likes it
  $self->files_to_delete($filename);
  eval{
    open (FILE,">$filename");
    print FILE ">".$seq->display_id."\n";
    print FILE $seq->seq."\n";
    print FILE $self->structure_constraint if ($self->structure_constraint);
    close FILE;
  };
  if ($@){
    $self->throw
      ("RNAFold: Error writing to seq file $@\n");
  };
  $self->files_to_delete($filename);
  return $filename;
}


=head2 encode_str

  Title      : encode_str
  Usage      : my $encoded_str = $runnable->encode_string($string)
  Function   : Does string length encoding to reduce size of structure string
             : splits strings if they are longer then 200 characters so they 
             : will fit in the transcript attribute table, gives a range value
             : at the start of the string indicating the start and stop positions of the 
             : structure on the transcript
  Returns    : String
  Exceptions : None
  Args       : String

=cut

sub encode_str{
  my ($self,$string)= @_;
  my @codes;
  my $start = 1;
  my $count=0;
  my $code;
  my @elements = split //,$string;
  my $last_chr = "";
  my @array =[];
  foreach my $chr (@elements){
    $count++;
    if ($chr eq $last_chr){
	push @array,$chr;
      }
    else {
      if ($code && length($code) > 200 && scalar(@array) == 1){
	push @codes,"$start:$count\t$code";
	$code = undef;
	$start = $count+1;
      }
      # Character has changed print STDERR it and the associated array length
      if (scalar(@array) > 1){
	$code .= scalar(@array);
	@array = [];
      }
      $code .= "$chr";
      $last_chr = $chr;
    }
  }
# last element
  if (scalar(@array) > 1){
    $code .= scalar(@array);
  }
  push @codes,"$start:$count\t$code";
  return \@codes;
}



##################################################################################
# Containers


=head2 sequence

  Title      : sequence
  Usage      : my $seq = $runnable->sequence
  Function   : Get/ set for the sequence object
  Returns    : Bio::PrimarySeq
  Exceptions : Throws if not passed a Bio::PrimarySeq object
  Args       : BioPrimarySeq

=cut

sub  sequence{
  my ($self, $sequence) = @_;
  if ($sequence){
    $self->throw("RNAFold: Sequence needs to be a Bio::PrimarySeq object not a $sequence\n")
		 unless ($sequence->isa("Bio::PrimarySeq"));
    $self->{'_sequence'} = $sequence;
  }
  return $self->{'_sequence'};
}

=head2 structure_constraint

  Title      : structure_constraint
  Usage      : my $str = $runnable->structure_constraint
  Function   : Get/ set for the structure string used by RNAfold to constrain
             : the folding prediction 
  Returns    : String
  Exceptions : None
  Args       : String

=cut

sub  structure_constraint{
  my ($self, $structure_constraint) = @_;
  if ($structure_constraint){
    $self->{'_structure_constraint'} = $structure_constraint;
  }
  return $self->{'_structure_constraint'};
}

=head2 structure

  Title      : structure
  Usage      : my $str = $runnable->structure
  Function   : Get/ set for the structure string predicted by RNAfold
  Returns    : String
  Exceptions : None
  Args       : String

=cut

sub  structure{
  my ($self, $structure) = @_;
  if ($structure){
    $self->{'_structure'} = $structure;
  }
  return $self->{'_structure'};
}

=head2 encoded_str

  Title      : encoded_str
  Usage      : my $str = $runnable->encoded_str
  Function   : Get/ set for the run-length encoded structure string predicted by RNAfold
  Returns    : String
  Exceptions : None
  Args       : String

=cut

sub  encoded_str{
  my ($self, $encoded_str) = @_;
  if ($encoded_str){
    $self->{'_encoded_str'} = $encoded_str;
  }
  return $self->{'_encoded_str'};
}

=head2 score

  Title      : score
  Usage      : my $score = $runnable->score
  Function   : Get/ set for the score predicted by RNAfold
  Returns    : String
  Exceptions : None
  Args       : String

=cut

sub  score{
  my ($self, $score) = @_;
  if ($score){
    $self->{'_score'} = $score;
  }
  return $self->{'_score'};
}
1;
