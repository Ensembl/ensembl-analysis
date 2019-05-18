# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
#
# Object for storing sequence dna_align_feature_history details
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory.pm - Stores details of an dna_align_feature_history run

=head1 SYNOPSIS
    my $obj    = new Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory(
        -id                     => $id,
        -seq_region_id          => $seq_region_id,
        -analysis_id            => $analysis_id,
        -align_feature_id_start => $align_feature_id_start,
        -align_feature_id_end   => $align_feature_id_end,
        -db_version             => $db_version,
        -date                   => $date,
        );

  The dna_align_feature_history table in Otter databases is not part of the 
  core Ensembl schema. This module provides an API for accessing the data
  from this table by creating a dna_align_feature_history object that can
  be attached to a dna_align_feature.

=head1 DESCRIPTION


=head1 CONTACT

Post questions to the EnsEMBL dev mailing list: <http://lists.ensembl.org/mailman/listinfo/dev>

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [..]   :  Takes a set of named arguments
  Example    : $dna_align_feature_history = new Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory(
        -id                     => $id,
        -seq_region_id          => $seq_region_id,
        -analysis_id            => $analysis_id,
        -align_feature_id_start => $align_feature_id_start,
        -align_feature_id_end   => $align_feature_id_end,
        -db_version             => $db_version,
        -date                   => $date,
        );

  Description: Creates a new DnaAlignFeatureHistory object
  Returntype : Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $seq_region_id, $analysis, $align_feature_id_start, 
      $align_feature_id_end, $db_version, $date, $adaptor) = 

	  rearrange([qw(ID
                        SEQ_REGION_ID
                        ANALYSIS
                        ALIGN_FEATURE_ID_START
                        ALIGN_FEATURE_ID_END
                        DB_VERSION
                        DATE
                        ADAPTOR
                     )],@args);


  $self->dbID           ($id);
  $self->adaptor        ($adaptor);
  $self->seq_region_id  ($seq_region_id);

  if($analysis) {
    if(!ref($analysis) || !$analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw('-ANALYSIS argument must be a Bio::EnsEMBL::Analysis not '.  $analysis);
    }
  } else {
    throw("Please attach an analysis");
  }
  $self->analysis ( $analysis);
  $self->align_feature_id_start ($align_feature_id_start);
  $self->align_feature_id_end ($align_feature_id_end);
  $self->db_version     ($db_version);
  $self->date           ($date);
  return $self; # success - we hope!
}


=head2 seq_region_id

  Arg [1]    : string $seq_region_id
  Example    : none
  Description: get/set for attribute seq_region_id
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_id {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{_seq_region_id} = $arg;
    }

    return $self->{_seq_region_id};
}

=head2 analysis

  Arg [1]    : (optional) Bio::EnsEMBL::Analysis $analysis
  Example    : $feature->analysis(new Bio::EnsEMBL::Analysis(...))
  Description: Getter/Setter for the analysis that is associated with
               this feature.  The analysis describes how this feature
               was derived.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : thrown if an invalid argument is passed
  Caller     : general
  Status     : Stable

=cut

sub analysis {
  my $self = shift;

  if(@_) {
    my $an = shift;
    if(defined($an) && (!ref($an) || !$an->isa('Bio::EnsEMBL::Analysis'))) {
      throw('analysis argument must be a Bio::EnsEMBL::Analysis');
    }
    $self->{'analysis'} = $an;
  }

  return $self->{'analysis'};
}


=head2 align_feature_id_start

  Arg [1]    : string $align_feature_id_start
  Example    : none
  Description: get/set for attribute align_feature_id_start
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub align_feature_id_start {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{_align_feature_id_start} = $arg;
    }

    return $self->{_align_feature_id_start};
}

=head2 align_feature_id_end

  Arg [1]    : string $align_feature_id_end
  Example    : none
  Description: get/set for attribute align_feature_id_end
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub align_feature_id_end {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{_align_feature_id_end} = $arg;
    }

    return $self->{_align_feature_id_end};
}

=head2 db_version

  Arg [1]    : string $db_version
  Example    : none
  Description: get/set for attribute db_version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub db_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_db_version} = $arg;
    }

    return $self->{_db_version};
}

=head2 date

  Arg [1]    : string $date
  Example    : none
  Description: get/set for attribute date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub date {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{_date} = $arg;
    }

    return $self->{_date};
}

1;
