=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

#
#
#
=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::ProteinAnnotation::PIRSF

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::Runnable::ProteinAnnotation::PIRSF->
    new ( 
    -query      => $query,
    -analysis   => $analysis,
    -database   => $db,
    -datfile    => $datfile,
    );
    

=head1 DESCRIPTION

=head1 CONTACT

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);


sub multiprotein{
  my ($self) = @_;
  return 1;
}

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);
  
  my ($threshold_file) = rearrange(['THRESHOLDS'], @args);

  if (not defined $threshold_file) {
    throw("You must supply a threshold file for PIRSF");
  } elsif (not -e $threshold_file) {
    throw("Threshold file '$threshold_file' could not be found");
  }

  $self->thresholds($threshold_file);

  throw("You must supply this runnable with a sequence object")
      if not ref($self->query) or not $self->query->isa("Bio::PrimarySeqI");

  return $self;
}


###############################
sub run_analysis {
  my ($self) = @_;

  my $dat = $self->process_thresholds_file($self->thresholds);

  my ($main_hmm, $sf_hmm) = split(/;/, $self->database);

  my $main_run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam->
      new(-query => $self->query,
          -analysis => $self->analysis,
          -database => $main_hmm,
          -options  => $self->options);
  $main_run->run;

  my $best_hit = $self->filter_main($main_run->output, $dat);

  if (defined $best_hit and exists($dat->{$best_hit->hseqname}->{children})) {
    my $sf_run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam->
        new(-query => $self->query,
            -analysis => $self->analysis,
            -database => $sf_hmm,
            -options  => $self->options);
    $sf_run->run;

    my @relevant_hits;
    foreach my $hit (@{$sf_run->output}) {
      if ($dat->{$best_hit->hseqname}->{children}->{$hit->hseqname}) {
        push @relevant_hits, $hit;
      }
    }

    if (@relevant_hits) {
      my $best_sf_hit = $self->filter_sub_family(\@relevant_hits, $dat);
      if (defined $best_sf_hit) {
        $best_hit = $best_sf_hit;
      }
    }
  }

  my $out = [];
  if (defined $best_hit) {
    push @$out, $best_hit;
  }

  $self->output($out);
}


sub parse_results {
  my ($self) = @_;
  
  # nothing to do
  return;
}



###############################
sub filter_main {
  my ($self, $hits, $dat) = @_;

  my @pass_hits;

  foreach my $hit (@$hits) {
    my $this_dat = $dat->{$hit->hseqname};
    
    if (not defined $this_dat->{mean_score}) {
      # this must be a sub-family hit; filter it out, hoping 
      # that the superfamily will be hit, and hence this sub-family 
      # hit in the refined search
      next;
    }

    my $len_prop = $hit->length / $self->query->length;
    my $len_diff = abs($self->query->length - $this_dat->{mean_length});
    my $score_diff = $hit->score - $this_dat->{mean_score};

    if ($len_prop >= 0.8 and
        ($score_diff >= -1*$this_dat->{std_score} or 
         ($hit->score >= $this_dat->{min_score} and 
          $score_diff >= -2.5 * $this_dat->{std_score})) and
        ($len_diff < 3.5 * $this_dat->{std_length})) {
      push @pass_hits, $hit;
    }
  }

  if (@pass_hits) {
    @pass_hits = sort { $a->p_value <=> $b->p_value} @pass_hits;
    return $pass_hits[0];
  } else {
    return undef;
  }
}


sub filter_sub_family {
  my ($self, $hits, $dat);
  
  my @pass_hits;
    foreach my $hit (@$hits) {
    my $this_dat = $dat->{$hit->hseqname};

    if ($hit->score >= $this_dat->{min_score}) {
      push @pass_hits, $hit;
    }
  }

  if (@pass_hits) {
    @pass_hits = sort { $a->pvalue <=> $b->pvalue} @pass_hits;
    return $pass_hits[0];
  } else {
    return undef;
  }

}


##############################
sub process_thresholds_file {
  my ($self, $file) = @_;

  my (%data);

  open(DAT, $file); 
  while(<DAT>) {
    /^\>(\S+)(\s*.+)?/ and do {
      my $acc = $1; my $rest = $2;
      
      my $line = <DAT>; $line = <DAT>;
      
      my @fields = split /\s+/, $line;

      $data{$acc} = {
        acc => $acc,
        mean_length => $fields[0],
        std_length  => $fields[1],
        min_score   => $fields[2],
        mean_score  => @fields > 3 ? $fields[3] : undef,
        std_score   => @fields > 3 ? $fields[4] : undef,
      };

      if ($rest and $rest =~ /child:\s+(.+)/) {
        foreach my $cacc (split /\s+/, $1) {
          $data{$acc}->{children}->{$cacc} = 1;
        }
      }
    }
  }

  return \%data;
}


sub thresholds {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_thresholds} = $val;
  }

  return $self->{_thresholds};
}
