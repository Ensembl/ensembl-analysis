# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::CPG
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::CPG

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::CPG->new
  (
   -query => $slice,
   -program => 'cpg',
  );
  $runnable->run;
  my @simple_features = @{$runnable->output};

=head1 DESCRIPTION

CPG expects to run the program CPG and produces SimpleFeature which can be
stored in the simple_feature table in the core database


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::CPG;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::CPG
  Arg [2]   : int, minimun length of hit
  Arg [3]   : int, minimun gc content of hit
  Arg [4]   : int, minimum oe of hit
  Function  : create a Bio::EnsEMBL::Analysis::Runnable::CPG
  Returntype: Bio::EnsEMBL::Analysis::Runnable::CPG
  Exceptions: none
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($min_length, $min_gc, $min_oe) = rearrange(['MIN_LENGTH',
                                                  'MIN_GC_CONTENT',
                                                  'MIN_OE',
                                                 ], @args);
  ##################
  #SETTING DEFAULTS#
  ##################
  $self->min_length(1000);
  $self->min_gc_content(50);
  $self->min_oe(0.6);
  #################

  $self->min_length($min_length);
  $self->min_gc_content($min_gc);
  $self->min_oe($min_oe);
  return $self;
}



=head2 containers

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::CPG
  Arg [2]   : int, minimun value
  Function  : ontainer for specified variable. This pod refers to the
  three methods below min_length, min_gc_content, min_oe. These are simple 
  containers which dont do more than hold and return an given value
  nothing is defined
  Returntype: int
  Exceptions: 
  Example   : 

=cut


sub min_length{
  my $self = shift;
  $self->{'min_length'} = shift if(@_);
  return $self->{'min_length'};
}

sub min_gc_content{
  my $self = shift;
  $self->{'min_gc_content'} = shift if(@_);
  return $self->{'min_gc_content'};
}

sub min_oe{
  my $self = shift;
  $self->{'min_oe'} = shift if(@_);
  return $self->{'min_oe'};
}



=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::CPG
  Arg [2]   : string, filename
  Function  : to open and parse the results file
  Returntype: none
  Exceptions: throws on failure to open or close the results file
  or if the results file doesnt exist
  Example   : 

=cut



sub parse_results{
  my ($self, $results) = @_;

  if(!$results){
    $results = $self->resultsfile;
  }
  if(!-e $results){
    throw("Can't parse an no existance results file ".$results.
          " CPG:parse_results");
  }
  my @output;
  open(CPG, $results) or throw("FAILED to open ".$results.
                               " CPG:parse_results");
  LINE:while(<CPG>){
    if (/\d+/){ #ignore introductory lines
      chomp;
      my @elements = split;
      my ($name, $start, $end, $score, $gc_content, $oe) 
        = @elements[0, 1, 2, 3, 6, 7];
      if($oe eq "-"){ 
        $oe = 0; 
      }
      my $length = $end - $start + 1;
      next LINE unless($length >= $self->min_length && 
                       $gc_content >= $self->min_gc_content &&
                       $oe >= $self->min_oe);
      my $sf = Bio::EnsEMBL::SimpleFeature->new();
      $sf->start($start);
      $sf->end($end);
      $sf->strand(0);
      $sf->score($score);
      $sf->display_label("oe = ".$oe);
      $sf->seqname($name);
      push(@output, $sf);
    }
  }
  $self->output(\@output);
  close(CPG) or throw("FAILED to close ".$results.
                      " CPG:parse_results");
}


1;
