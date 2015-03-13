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

Bio::EnsEMBL::Analysis::RunnableDB::RNASeq_UTRCorrection

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::RunnableDB::RNASeq_UTRCorrection->new(
    );

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

Take RNASeq gene models as input to correct the UTR. It tries to find reads
with polyA tail to define the 3' UTR. In any other cases it search for a drop
in the read coverage. When there is no UTR it tries to create one.
Because we overpredict the length of the UTR, if a model already have UTR,
it can only be shorten but never extended. We need better data like 3' pulldown
data to be able to clearly define the UTR.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::RNASeq_UTRCorrection;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::RNASeq_UTRCorrection;
use Bio::EnsEMBL::Analysis::Runnable::RNASeq_UTRCorrection;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->db->disconnect_when_inactive(1);
  $self->read_and_check_config($RNASEQ_UTRCORRECTION_CONFIG_BY_LOGIC);
  return $self;
}

# fetch input
# get the transcript in question
# get all the reads that overlap
# do we need to hold them in memory or can we stream them into the file?

sub fetch_input {
  my ($self) = @_;
  my $gene_db = $self->get_dbadaptor($self->GENES_DB);
  my $slice_adaptor =  $gene_db->get_SliceAdaptor;
  my $slice = $slice_adaptor->fetch_by_name($self->input_id);
  my $genes = $slice->get_all_Genes($self->LOGIC_NAME, undef, 1, $self->SOURCE, $self->BIOTYPE);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::RNASeq_UTRCorrection->new
    (
     -analysis     => $self->analysis,
     -genes      => $genes,
     -bam_file => $self->BAM_FILE,
     -slice => $slice,
    );
  $self->runnable($runnable);
}

=head2 get_adaptor

  Function  : Return the proper DBAdaptor to write the features using write_output from RunnableDB
  Returntype: Bio::EnsEMBL::DBSQL::GeneAdaptor
  Exceptions: None

=cut

sub get_adaptor {
    my $self = shift;

    if (!exists $self->{outputdb_adaptor}) {
       $self->{outputdb_adaptor} = $self->get_dbadaptor($self->OUTPUT_DB)->get_GeneAdaptor;
    }
    return $self->{outputdb_adaptor};
}


###########################################################
# containers

sub GENES_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_GENES_DB'} = $value;
  }

  if (exists($self->{'_CONFIG_GENES_DB'})) {
    return $self->{'_CONFIG_GENES_DB'};
  } else {
    return undef;
  }
}

sub BAM_FILE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BAM_FILE'} = $value;
  }

  if (exists($self->{'_CONFIG_BAM_FILE'})) {
    return $self->{'_CONFIG_BAM_FILE'};
  } else {
    return undef;
  }
}

sub LOGIC_NAME {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_LOGIC_NAME'} = $value;
  }

  if (exists($self->{'_CONFIG_LOGIC_NAME'})) {
    return $self->{'_CONFIG_LOGIC_NAME'};
  } else {
    return undef;
  }
}

sub SOURCE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_SOURCE'} = $value;
  }

  if (exists($self->{'_CONFIG_SOURCE'})) {
    return $self->{'_CONFIG_SOURCE'};
  } else {
    return undef;
  }
}

sub BIOTYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BIOTYPE'} = $value;
  }

  if (exists($self->{'_CONFIG_BIOTYPE'})) {
    return $self->{'_CONFIG_BIOTYPE'};
  } else {
    return undef;
  }
}

sub OUTPUT_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUTPUT_DB'} = $value;
  }

  if (exists($self->{'_CONFIG_OUTPUT_DB'})) {
    return $self->{'_CONFIG_OUTPUT_DB'};
  } else {
    return undef;
  }
}

1;
