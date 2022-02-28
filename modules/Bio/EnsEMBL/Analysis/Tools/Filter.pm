=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Filter

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Tools::Filter;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


=head2 new

 Arg [COVERAGE]                 : Int $cov, minimum coverage
 Arg [PERCENT_ID]               : Int $pid, minimum percent identity
 Arg [BEST_IN_GENOME]           : Boolean, only getting the best model
 Arg [REJECT_PROCESSED_PSEUDOS] : Boolean, reject models if a spliced alignment already exist
 Arg [VERBOSITY]                : Int $verbose, the higher the value the more talkative is the filter
 Description                    : Constructor
 Returntype                     : Bio::EnsEMBL::Analysis::Tools::Filter object
 Exceptions                     : None

=cut

sub new {
  my ($class,@args) = @_;

  my $self = bless {}, $class;

  my ($min_coverage,
      $min_percent,
      $best_in_genome,
      $rpp,
      $verbosity) =
        rearrange([
                   'COVERAGE',
                   'PERCENT_ID',
                   'BEST_IN_GENOME',
                   'REJECT_PROCESSED_PSEUDOS',
                   'VERBOSITY',
                  ], @args);

  $self->min_coverage($min_coverage) if (defined $min_coverage);
  $self->min_percent($min_percent) if (defined $min_percent);
  $self->best_in_genome($best_in_genome) if (defined $best_in_genome);
  $self->reject_processed_pseudos($rpp) if (defined $rpp);
  $self->verbosity($verbosity) if (defined $verbosity);
  return $self;
}


=head2 filter_results

 Arg [1]    : Arrayref of Bio::EnsEMBL::Feature
 Description: Only this method should be used. All inheriting class must implement it
              You should not pass arguments in the methods, they should all have been set
              at the creation of the object through a hashref
 Returntype : Arrayref of Bio::EnsEMBL::Transcript
 Exceptions : Throws if Arg[1] is not an arrayref

=cut

sub filter_results {
  my ($self, $objects) = @_;

  throw('You should give an arrayref of objects') unless (ref($objects) eq 'ARRAY');
  throw('You should implement the filter_results method');
}


=head2 min_coverage

 Arg [1]    : Int $cov
 Description: Getter/setter for the minimum hit coverage
 Returntype : Int or undef if not set
 Exceptions : None

=cut

sub min_coverage {
  my $self = shift;
  $self->{'_min_coverage'} = shift if(@_);

  return exists($self->{'_min_coverage'}) ? $self->{'_min_coverage'} : undef;
}


=head2 min_percent

 Arg [1]    : Int pid
 Description: Getter/setter for the minimum percentage identity allowed
 Returntype : Int or undef if not set
 Exceptions : None

=cut

sub min_percent {
  my $self = shift;
  $self->{'_min_percent'} = shift if(@_);

  return exists($self->{'_min_percent'}) ? $self->{'_min_percent'} : undef;
}


=head2 best_in_genome

 Arg [1]    : Boolean $best
 Description: Getter/setter for allowing only the best hits or any hits
 Returntype : Boolean
 Exceptions : None

=cut

sub best_in_genome {
  my $self = shift;
  $self->{'_best_in_genome'} = shift if(@_);

  return exists($self->{'_best_in_genome'}) ? $self->{'_best_in_genome'} : 0;
}


=head2 reject_processed_pseudos

 Arg [1]    : Boolean
 Description: Getter/setter for rejecting the single exon hits which have a spliced
              alignment in the genome
 Returntype : Boolean
 Exceptions : None

=cut

sub reject_processed_pseudos {
  my $self = shift;
  $self->{'_reject_processed_pseudos'} = shift if(@_);

  return exists($self->{'_reject_processed_pseudos'}) ? $self->{'_reject_processed_pseudos'} : 0;
}


=head2 verbosity

 Arg [1]    : Int
 Description: Getter/setter for the verbosity of the filter. The higher Arg[1] the more verbose it is.
 Returntype :
 Exceptions :

=cut

sub verbosity {
  my $self = shift;
  $self->{'_verbosity'} = shift if(@_);

  return exists($self->{'_verbosity'}) ? $self->{'_verbosity'} : 0;
}

1;
