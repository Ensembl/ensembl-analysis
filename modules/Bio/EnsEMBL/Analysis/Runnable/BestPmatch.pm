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

Bio::EnsEMBL::Analysis::Runnable::BestPmatch - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::BestPmatch;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($pafs, $min_coverage) = rearrange(['PROTEIN_HITS', 'MIN_COVERAGE'], @args);
  
  ###SETTING DEFAULTS###
  $self->min_coverage(25);
  #######################

  $self->protein_hits($pafs);
  $self->min_coverage($min_coverage);
  return $self;
}


sub protein_hits{
  my ($self, $ref) = @_;
  if($ref){
    if(ref($ref) ne "ARRAY"){
      throw("Runnable::BestPmatch Must pass protein hits an array ref not ".$ref);
    }
    $self->{protein_hits} = $ref;
  }
  return $self->{protein_hits};
}


sub min_coverage{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'min_coverage'} = $arg;
  }
  return $self->{'min_coverage'};
}

sub protein_hash{ 
  my ($self, $hash) = @_;
  if(!$self->{'_proteins'}){
    $self->{'_proteins'} = {};
  }
  if($hash){
    throw("Must pass protein_hash a hashref not ".$hash) 
      if(ref($hash) ne 'HASH');
    $self->{'_proteins'} = $hash;
  }
  return $self->{'_proteins'};
}


sub run{
  my ($self) = @_;
  my %prots;
  foreach my $hit(@{$self->protein_hits}){
    push (@{$prots{$hit->hseqname}}, $hit);
  }
  $self->protein_hash(\%prots);
  my $hits = $self->prune_hits;
  my %unique;
  foreach my $hit(@$hits){
    my $string = id($hit)."-".coord_string($hit);
    if(!$unique{$string}){
      $hit->analysis($self->analysis);
      $unique{$string} = $hit;
    }
  }
  my @output = values(%unique);
  $self->output(\@output);
}


sub prune_hits{
  my ($self, $hits) = @_;

  $hits = $self->protein_hits if(!$hits);
  my @chosen;
  my %prots = %{$self->protein_hash};
 PROTEIN:foreach my $p(keys(%{$self->protein_hash})){
    my $allhits = $prots{$p};
    my @sorted = sort {$b->score <=> $a->score} @$allhits;

    my $first = shift(@sorted);
    my $score_boundary = $first->score() - 2; 
    my $lower_threshold = $self->min_coverage;
    next PROTEIN if $first->score < $lower_threshold;
    push(@chosen, $first);
  PRUNE:foreach my $hit(@sorted){
      last PRUNE if($hit->score < $score_boundary);
      last PRUNE if($hit->score < $lower_threshold);
      push(@chosen, $hit);
    }
  }
  return \@chosen;
}

1;
