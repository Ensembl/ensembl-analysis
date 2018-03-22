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
  my $lower_threshold = $self->min_coverage;
 PROTEIN:foreach my $p(keys(%{$self->protein_hash})){
    my $allhits = $prots{$p};
    my @sorted = sort {$b->score <=> $a->score} @$allhits;

    my $first = shift(@sorted);
    my $score_boundary = $first->score() - 2; 
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
