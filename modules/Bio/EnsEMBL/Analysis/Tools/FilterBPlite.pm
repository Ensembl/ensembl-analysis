package Bio::EnsEMBL::Analysis::Tools::FilterBPlite;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use Bio::EnsEMBL::Analysis::Tools::FeatureFilter;
use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Tools::BPliteWrapper);


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($threshold_type, $threshold, $coverage, $filter) = rearrange
    (['THRESHOLD_TYPE', 'THRESHOLD', 'COVERAGE', 'FILTER'], @args);

  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->coverage(10);
  $self->filter(1);
  ######################

  $self->threshold_type($threshold_type);
  $self->threshold($threshold);
  $self->coverage($coverage) if(defined($coverage));
  $self->filter($filter) if(defined($filter));
  return $self;
}


sub threshold_type{
  my $self = shift;
  $self->{'threshold_type'} = shift if(@_);
  return $self->{'threshold_type'};
}


sub threshold{
  my $self = shift;
  $self->{'threshold'} = shift if(@_);
  return $self->{'threshold'};
}

sub coverage{
  my $self = shift;
  $self->{'coverage'} = shift if(@_);
  return $self->{'coverage'};
}

sub filter{
  my $self = shift;
  $self->{'filter'} = shift if(@_);
  return $self->{'filter'};
}

sub get_hsps{
  my ($self, $parser) = @_;
  my $regex = $self->regex;
  my @output;
  my $ids;
  if($self->filter){
    $ids = $self->filter_hits($parser);
  }
  my $second = $self->get_parser($self->filename);
 NAME: while(my $sbjct = $second->nextSbjct){
    if($self->filter && !($ids->{$sbjct->name})){
      next NAME;
    }
    my ($name) = $sbjct->name =~ /$regex/;
    throw("Error parsing name from ".$sbjct->name." check your ".
          "blast setup and blast headers") unless($name);
  HSP: while (my $hsp = $sbjct->nextHSP) {
      if($self->is_hsp_valid($hsp)){
        push(@output, $self->split_hsp($hsp, $name));
      }
    }
  }
  $self->output(\@output);
}

sub filter_hits{
  my ($self, $parser) = @_;
  my %ids;
  my @features;
 SUB:while(my $sbjct = $parser->nextSbjct){
    my $name = $sbjct->name;
  HSP:while (my $hsp = $sbjct->nextHSP) {
      if($self->is_hsp_valid($hsp)){
        my $qstart = $hsp->query->start();
        my $hstart = $hsp->subject->start();
        
        my $qend   = $hsp->query->end();
        my $hend   = $hsp->subject->end();      
        
        my $qstrand = $hsp->query->strand();
        my $hstrand = $hsp->subject->strand();
        
        my $score  = $hsp->score;
        my $p_value = $hsp->P;
        my $percent = $hsp->percent;

        my $fp = $self->feature_factory->create_feature_pair
          ($qstart, $qend, $qstrand, $score, $hstart,
           $hend, $hstrand, $name, $percent, $p_value);

        push(@features,$fp);
      }
    }
  }
  my $search = Bio::EnsEMBL::Analysis::Tools::FeatureFilter->new
    (
     -coverage => $self->coverage
    );
  my @newfeatures = @{$search->filter_results(\@features)};
  foreach my $f (@newfeatures) {
    my $id = $f->hseqname;
    $ids{$id} = 1;
  }
  
  return \%ids;
}

sub is_hsp_valid{
  my ($self, $hsp) = @_;
  if($self->threshold_type){
    if ($self->threshold_type eq "PID") {
      return 0 if ($hsp->percent < $self->threshold);
    } elsif ($self->threshold_type eq "SCORE") {
      return 0 if ($hsp->score < $self->threshold);
    } elsif ($self->threshold_type eq "PVALUE") {
      return 0 if($hsp->P > $self->threshold);
    } 
  }
  return $hsp;
}
