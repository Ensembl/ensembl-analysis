# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Dust
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Dust

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

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Runnable::Dust;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);



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

=head2 containers

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
  $self->{'min_length'} = shift if(@_);
  return $self->{'min_length'};
}

sub window_size{
  my $self = shift;
  $self->{'min_gc_content'} = shift if(@_);
  return $self->{'min_gc_content'};
}

sub word_size{
  my $self = shift;
  $self->{'min_oe'} = shift if(@_);
  return $self->{'min_oe'};
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
    if (/(\d+)\.\.(\d+)/) {
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
  $self->output(\@output);
  close(DUST) or throw("FAILED to close ".$results);
}
