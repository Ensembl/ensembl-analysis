=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Bam2Introns

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::RunnableDB::Bam2Introns->new(
    );

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

Takes an input id of a rough transcript and fetches reads associated with that model
from a BAM file, then runs a spliced exonerate alignment against genomic DNA or the
rough transcript. Writes output as SAM files.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::Bam2Introns;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::SeqFeature::Lite;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::DB::Sam;

$| = 1;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ($percent_id, $coverage, $missmatch) = rearrange( [qw(PERCENT_ID COVERAGE MISSMATCH)],@args );
  $self->PERCENT_ID($percent_id);
  $self->COVERAGE($coverage);
  $self->MISSMATCH($missmatch);
  return $self;
}

sub run  {
  my ( $self) = @_;
  # set up the output files
#  my $query_seq = $self->create_filename("B2I_reads","fa");
#  my $genomic_seq = $self->create_filename("B2I_transcript","fa");
#  $self->files_to_delete($query_seq);
#  $self->files_to_delete($genomic_seq);
#  $self->query_file($query_seq);
#  $self->target_file($genomic_seq);
#  $self->write_seq_file($self->query, $genomic_seq);
#  $self->write_seq_file($self->query_seqs, $query_seq);
  # now to run the runnable
  $self->SUPER::run();  # attach the read seq to the output features
  $self->process_features;
}

sub process_features {
  my ($self) = @_;
  my $features = $self->clean_output;
  my @new_features;
  foreach my $feat ( @$features ) {
    # filter on coverage and percent id
    next unless $feat->percent_id >= $self->PERCENT_ID;
    next unless $feat->hcoverage >= $self->COVERAGE;
    next unless $feat->{"_intron"};
    # check missmatches
    if ( $self->MISSMATCH ) {
      my $aligned_length = abs($feat->hend - $feat->hstart) +1;
      my $matches = $aligned_length *  $feat->percent_id / 100;
      my $missmatches = ( $aligned_length - $matches) / $aligned_length * 100;
      next if ($missmatches > $self->MISSMATCH);
    }
    push @new_features, $feat;
  }

  $self->output(\@new_features);
}



###########################################################
# containers

sub MISSMATCH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MISSMATCH'} = $value;
  }

  if (exists($self->{'_CONFIG_MISSMATCH'})) {
    return $self->{'_CONFIG_MISSMATCH'};
  } else {
    return undef;
  }
}

sub PERCENT_ID {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PERCENT_ID'} = $value;
  }

  if (exists($self->{'_CONFIG_PERCENT_ID'})) {
    return $self->{'_CONFIG_PERCENT_ID'};
  } else {
    return undef;
  }
}

sub COVERAGE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_COVERAGE'} = $value;
  }

  if (exists($self->{'_CONFIG_COVERAGE'})) {
    return $self->{'_CONFIG_COVERAGE'};
  } else {
    return undef;
  }
}

1;
