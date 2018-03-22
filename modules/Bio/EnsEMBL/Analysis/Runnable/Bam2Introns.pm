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
use strict;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature');

=head2 new

 Arg [PERCENT_ID]: Integer
 Arg [COVERAGE]  : Integer
 Arg [MISSMATCH] : Integer
 Description     : Create a new Bio::EnsEMBL::Analysis::Runnable::Bam2Introns object
 Returntype      : Bio::EnsEMBL::Analysis::Runnable::Bam2Introns
 Exceptions      : None

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ($percent_id, $coverage, $missmatch) = rearrange( [qw(PERCENT_ID COVERAGE MISSMATCH)],@args );
  $self->PERCENT_ID($percent_id);
  $self->COVERAGE($coverage);
  $self->MISSMATCH($missmatch);
  return $self;
}


=head2 run

 Arg [1]    : None
 Description: Align the reads to a region of the genome to find candidate introns
 Returntype : None
 Exceptions : None

=cut

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


=head2 process_features

 Arg [1]    : None
 Description: Filter the result from Exonerate based on coverage and percent of identity
 Returntype : None
 Exceptions : None

=cut

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

=head2 MISSMATCH

 Arg [1]    : (optional) Integer
 Description: Getter/setter
 Returntype : Integer
 Exceptions : None

=cut

sub MISSMATCH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MISSMATCH'} = $value;
  }

  if (exists($self->{'_CONFIG_MISSMATCH'})) {
    return $self->{'_CONFIG_MISSMATCH'};
  } else {
    return;
  }
}


=head2 PERCENT_ID

 Arg [1]    : (optional) Integer
 Description: Getter/setter
 Returntype : Integer
 Exceptions : None

=cut

sub PERCENT_ID {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PERCENT_ID'} = $value;
  }

  if (exists($self->{'_CONFIG_PERCENT_ID'})) {
    return $self->{'_CONFIG_PERCENT_ID'};
  } else {
    return;
  }
}


=head2 COVERAGE

 Arg [1]    : (optional) Integer
 Description: Getter/setter
 Returntype : Integer
 Exceptions : None

=cut

sub COVERAGE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_COVERAGE'} = $value;
  }

  if (exists($self->{'_CONFIG_COVERAGE'})) {
    return $self->{'_CONFIG_COVERAGE'};
  } else {
    return;
  }
}

1;
