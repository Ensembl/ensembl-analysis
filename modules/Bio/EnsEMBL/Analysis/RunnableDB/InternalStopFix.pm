=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild -

=head1 SYNOPSIS


=head1 DESCRIPTION

This class removes the internal stop codon from models. If the model has one stop codon,
the model's biotype is changed to the value specified in STOP_CODON_BIOTYPE and a new model
is created with a 3bp intron with the biotype EDITED_BIOTYPE. If the model have more than
1 stop codon, it only changes the biotype to the value specified in STOP_CODON_BIOTYPE


=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::InternalStopFix;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::InternalStopFix;
use Bio::EnsEMBL::Analysis::Config::InternalStopFix;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($INTERNALSTOP_FIX_CONFIG_BY_LOGIC);

  return $self;
}

sub fetch_input {
    my $self = shift;

    my $db = $self->get_dbadaptor($self->GENES_DB);
    my $slice = $db->get_SliceAdaptor->fetch_by_name($self->input_id);
    my $genes = $slice->get_all_Genes($self->LOGIC_NAME, undef, 1, $self->SOURCE, $self->BIOTYPE);
    $self->input_is_void(1) if (scalar(@$genes) == 0);
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::InternalStopFix->new(
        -analysis => $self->analysis,
        -stop_codon_biotype => $self->STOP_CODON_BIOTYPE,
        -edited_biotype => $self->EDITED_BIOTYPE,
        -genes => $genes,
    );
    $self->runnable($runnable);
}

sub write_output {
    my $self = shift;

    my $genes = $self->output;
    my $gene_adaptor = $self->get_dbadaptor($self->OUTPUT_DB)->get_GeneAdaptor();
    my $update_gene_adaptor = $self->get_dbadaptor($self->GENES_DB)->get_GeneAdaptor();
    foreach my $gene (@$genes) {
        if ($gene->dbID) {
            $update_gene_adaptor->update($gene);
        }
        else {
            $gene_adaptor->store($gene);
        }
    }
}

sub GENES_DB {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_genes_db} = $val;
  }
  return $self->{_genes_db};
}

sub OUTPUT_DB {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_output_db} = $val;
  }
  return $self->{_output_db};
}

sub EDITED_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_edited_biotype} = $val;
  }
  return $self->{_edited_biotype};
}

sub STOP_CODON_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_stop_codon_biotype} = $val;
  }
  return $self->{_stop_codon_biotype};
}

sub LOGIC_NAME {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_logic_name} = $val;
  }
  return $self->{_logic_name};
}

sub SOURCE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_source} = $val;
  }
  return $self->{_source};
}

sub BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_biotype} = $val;
  }
  return $self->{_biotype};
}

1;
