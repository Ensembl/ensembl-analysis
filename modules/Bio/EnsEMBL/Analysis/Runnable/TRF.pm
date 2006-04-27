# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::TRF
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::TRF

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::TRF->
  new(
      -query => $slice,
      -program => 'trf',
     );
  $runnable->run;
  my @repeats = @{$runnable->output};

=head1 DESCRIPTION

TRF expects to run the program trf
and produce RepeatFeatures which can be stored in the repeat_feature
and repeat_consensus tables in the core database


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Runnable::TRF;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::TRF
  Arg [2]   : int, match
  Arg [3]   : int, mismatch
  Arg [4]   : int, delta
  Arg [5]   : int, match probablilty
  Arg [6]   : int, indel probablilty
  Arg [7]   : int, minimun score
  Arg [8]   : int, maximun period size
  Function  : contruct a new Bio::EnsEMBL::Analysis::Runnable::TRF
  runnable
  Returntype: Bio::EnsEMBL::Analysis::Runnable::TRF
  Exceptions: none
  Example   : 

=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($match, $mismatch, $delta, $pm, $pi, $minscore,
      $maxperiod) = rearrange(['MATCH', 'MISMATCH', 'DELTA', 'PM', 'PI',
                               'MINSCORE', 'MAX_PERIOD'], @args);

  ######################
  #SETTING THE DEFAULTS#
  ######################
  $match = 2 unless(defined $match);
  $mismatch = 5 unless(defined $mismatch);
  $delta = 7 unless(defined $delta);
  $pm = 80 unless(defined $pm);
  $pi = 10 unless(defined $pi);
  $minscore = 40 unless(defined $minscore);
  $maxperiod = 500 unless(defined $maxperiod);
  ######################
  $self->program('trf') if(!$self->program);
  $self->match     ($match)     if defined $match;
  $self->mismatch  ($mismatch)  if defined $mismatch;
  $self->delta     ($delta)     if defined $delta;
  $self->pm        ($pm)        if defined $pm;
  $self->pi        ($pi)        if defined $pi;
  $self->minscore  ($minscore)  if defined $minscore;
  $self->maxperiod ($maxperiod) if defined $maxperiod;
  return $self;
}


# Please use: trf File Match Mismatch Delta PM PI Minscore MaxPeriod 
# [options]
# Where: (all scores are positive)
#   File = sequences input file
#   Match  = matching weight
#   Mismatch  = mismatching penalty
#   Delta = indel penalty
#   PM = match probability (whole number)
#   PI = indel probability (whole number)
#   Minscore = minimum alignment score to report
#   MaxPeriod = maximum period size to report
#   [options] = one or more of the following :
#                -m    masked sequence file
#                -f    flanking sequence
#                -h    avoid html output 



=head2 containers

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::CPG
  Arg [2]   : int,
  Function  : container for specified variable. This pod refers to the
  7 methods below match, mismatch, delta, pm, pi, minscore and maxperiod. 
  These are simple containers which dont do more than hold and return an 
  given value nothing is defined
  Returntype: int
  Exceptions: 
  Example   : 

=cut

sub match {
    my ($self, $arg) = @_;

    $self->{'_match'} = $arg if defined $arg;
    
    return $self->{'_match'};
}

sub mismatch {
    my ($self, $arg) = @_;

    $self->{'_mismatch'} = $arg if defined $arg;
    return $self->{'_mismatch'};
}

sub delta {
    my ($self, $arg) = @_;

    $self->{'_delta'} = $arg if defined $arg;

    return $self->{'_delta'};
}

sub pm {
    my ($self, $arg) = @_;

    $self->{'_pm'} = $arg if defined $arg;

    return $self->{'_pm'};
}

sub pi {
    my ($self, $arg) = @_;

    $self->{'_pi'} = $arg if defined $arg;

    return $self->{'_pi'};
}

sub minscore {
    my ($self, $arg) = @_;

    $self->{'_minscore'} = $arg if defined $arg;

    return $self->{'_minscore'};
}

sub maxperiod {
    my ($self, $arg) = @_;

    $self->{'_maxperiod'} = $arg if defined $arg;

    return $self->{'_maxperiod'};
}



=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::TRF
  Arg [2]   : string, program name
  Function  : create and open a commandline for the program trf
  Returntype: none
  Exceptions: throws if the program in not executable or if the results
  file doesnt exist
  Example   : 

=cut



sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  my @arg_list = ($self->match, $self->mismatch, $self->delta, $self->pm,
                  $self->pi, $self->minscore, $self->maxperiod);
	$self->options(join(' ', @arg_list));
  throw($program." is not executable TRF::run_analysis ") 
    unless($program && -x $program);
  # Use this command if you are using TRF version 4.0.0
  my $command = $program." ".$self->queryfile." ".$self->options." -d -h";

  # Use this command if your TRF program is previous to the 4.0.0 version
  # my $command = $program." ".$self->queryfile." ".$self->options." -d";

  print "Running analysis ".$command."\n";
  my $exit = system($command);
  foreach my $file(glob $self->queryfile."*"){
    $self->files_to_delete($file);
  }
  $self->resultsfile($self->queryfile.".".join('.', @arg_list).".dat");
  if(!-e $self->resultsfile){
    throw("FAILED to run TRF on ".$self->queryfile." ".
          $self->resultsfile." has not been produced ".
          "TRF:run_analysis");
  }
}


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::TRF
  Arg [2]   : string, results filename
  Function  : parse the specifed file and produce RepeatFeatures
  Returntype: nine
  Exceptions: throws if fails to open or close the results file
  Example   : 

=cut


sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  my $ff = $self->feature_factory;
  open(OUT, "<".$results) or throw("FAILED to open ".$results.
                                   "TRF:parse_results");
  my  @output;
 LINE:while(<OUT>){
    my $seqname;
    if (/Sequence: (\w+)/) {
      $seqname = $1;
    }
    next LINE unless (/^\d/); # ignore introductory lines
    my ($start, $end, $period_size,$copy_number, $consensus_size, 
        $percent_matches, $percent_indels, $score, $A, $C, $G, $T,
        $entropy, $mer) = split;
    if (($score < 50 && $percent_matches >= 80 && $copy_number > 2 && $period_size < 10) || ($copy_number >= 2 && $percent_matches >= 70 && $score >= 50)){
      my $rc = $ff->create_repeat_consensus('trf', 'trf', '', $mer);
      my $rf = $ff->create_repeat_feature($start, $end, 0, $score, 1,
                                         ($end - $start +1), $rc);
      push(@output, $rf);
    }
  }
  print "No tandem repeats found" if(@output == 0);
  $self->output(\@output);
  close(OUT) or throw("FAILED to close ".$results.
                      "TRF:parse_results");
}

1;
