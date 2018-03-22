=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Dust - 

=head1 SYNOPSIS

  my $dust = Bio::EnsEMBL::Analysis::Runnable::Dust->
  new(
      -query => $slice,
      -program => 'tcdust',
     );
  $dust->run;
  my @repeats = @{$dust->output};

=head1 DESCRIPTION

Dust is a wrapper for the tcdust program which runs the dust algorithm
to identify and mask simple repeats

=cut


package Bio::EnsEMBL::Analysis::Runnable::Dust;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent qw(Bio::EnsEMBL::Analysis::Runnable);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Dust
  Arg [2]   : int, masking threshold
  Arg [3]   : int, word size
  Arg [4]   : int, window size
  Function  : create a new  Bio::EnsEMBL::Analysis::Runnable::Dust
  Returntype: Bio::EnsEMBL::Analysis::Runnable::Dust
  Exceptions: 
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($level, $word_size, $window_size) = rearrange(['MASKING_THRESHOLD',
                                                     'WORD_SIZE', 
                                                     'WINDOW_SIZE',
                                                    ], @args);
  ##################
  #SETTING DEFAULTS#
  ##################
  $self->program('tcdust') if(!$self->program);
  ##################

  $self->level($level) if($level);
  $self->word_size($word_size) if($word_size);
  $self->window_size($window_size) if($window_size);
  return $self;
}

#container methods

=head2 level

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Dust
  Arg [2]   : int for specified value (for more info see tcdust)
  Function  : container for specified variable. This pod refers to the
  three methods below level, window size and word size. These are simple 
  containers which dont do more than hold and return an given value
  nothing is defined
  Returntype: int
  Exceptions: 
  Example   : 

=cut


sub level{
  my $self = shift;
  $self->{'level'} = shift if(@_);
  return $self->{'level'};
}

=head2 window_size

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Dust
  Arg [2]   : int for specified value (for more info see tcdust)
  Function  : container for specified variable. This pod refers to the
  three methods below level, window size and word size. These are simple 
  containers which dont do more than hold and return an given value
  nothing is defined
  Returntype: int
  Exceptions: 
  Example   : 

=cut
sub window_size{
  my $self = shift;
  $self->{'window_size'} = shift if(@_);
  return $self->{'window_size'};
}

=head2 word_size

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Dust
  Arg [2]   : int for specified value (for more info see tcdust)
  Function  : container for specified variable. This pod refers to the
  three methods below level, window size and word size. These are simple 
  containers which dont do more than hold and return an given value
  nothing is defined
  Returntype: int
  Exceptions: 
  Example   : 

=cut
sub word_size{
  my $self = shift;
  $self->{'word_size'} = shift if(@_);
  return $self->{'word_size'};
}



#utility methods




=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Dust
  Arg [2]   : string, program name
  Function  : constructs a commandline and runs the program passed
  in, the generic method in Runnable isnt used as Dust doesnt
  fit this module
  Returntype: none
  Exceptions: throws if run failed because system doesnt
  return 0
  Example   : 

=cut

sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable Dust::run_analysis ") 
    unless($program && -x $program);
  my $command = $self->program;
  $command .= " -l ".$self->level if($self->level);
  $command .= " -k ".$self->word_size if($self->word_size);
  $command .= " -w ".$self->window_size if($self->window_size);
  $command .= " -x ".$self->queryfile." > ".$self->resultsfile;
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
}


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Dust
  Arg [2]   : string, filename
  Function  : open and parse the results file into repeat
  features
  Returntype: none 
  Exceptions: throws on failure to open or close output file
  Example   : 

=cut

sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  my $ff = $self->feature_factory;
  my @output;
  open(DUST, $results) or throw("FAILED to open ".$results);
 LINE:while(<DUST>){
    chomp;
    next LINE if(/^>/);
    if (/(\d+)\D+(\d+)/) {
      my ($start, $end) = ($1, $2);
      $start++;
      $end++;
      my $rc = $ff->create_repeat_consensus('dust', 'dust', 'simple', 'N');
      my $rf = $ff->create_repeat_feature($start, $end, 0, 0, $start,
                                          $end, $rc, $self->query->name,
                                          $self->query);
      push(@output, $rf);
    }
  }
  $self->output(\@output,1);
  close(DUST) or throw("FAILED to close ".$results);
}

1;
