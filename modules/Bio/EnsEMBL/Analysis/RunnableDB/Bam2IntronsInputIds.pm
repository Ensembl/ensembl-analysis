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

Bio::EnsEMBL::Analysis::RunnableDB::Bam2IntronsInputIds -

=head1 SYNOPSIS

  [bam2intronsinputids]
  module=Bam2IntronsInputIds
  parameters=logic_name=>bam2introns, genes_db=>ROUGHDB, input_id_type=>RNASEQID

=head1 DESCRIPTION

  Merge the bam files to create the "super" bam file needed by
  Bam2Genes, the rough model step.
  It uses samtools to check but it uses picard for merging,
  no config file as everything is fetched from the analysis table
  The module is looking for files named: *_sorted.bam

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Bam2IntronsInputIds;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new{
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my $hash = $self->parameters_hash;
  throw('You need to specify a genes database') unless (exists $hash->{genes_db});
  $self->genes_db($hash->{genes_db});
  throw('You need to specify a logic_name') unless (exists $hash->{logic_name});
  $self->logic_name($hash->{logic_name});
  throw('You need to specify a input_id_type') unless (exists $hash->{input_id_type});
  $self->input_id_type($hash->{input_id_type});
  $self->pipeline_db($hash->{pipeline_db}) if (exists $hash->{pipeline_db});
  $self->genes_logic_name($hash->{genes_logic_name}) if (exists $hash->{genes_logic_name});
  return $self;
}

sub fetch_input {
    my ($self) = @_;

    my $db = $self->get_dbadaptor($self->genes_db);
    my $seq_region_id = $self->db->get_SliceAdaptor->fetch_by_name($self->input_id)->get_seq_region_id;
    my $sql_query = 'SELECT stable_id FROM gene WHERE seq_region_id = '.$seq_region_id;
    if ($self->genes_logic_name) {
        my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($self->genes_logic_name);
        $sql_query .= ' AND analysis_id = '.$analysis->dbID;
    }
    my $sth = $db->dbc->prepare($sql_query);
    $sth->execute;
    my $stable_id;
    my @stable_ids;
    $sth->bind_columns(\$stable_id);
    while ($sth->fetch) {
        push (@stable_ids, $stable_id);
    }
    $self->output(\@stable_ids);
}

sub run {
    my ($self) = @_;

    return 1;
}

sub write_output {
    my ($self) = @_;

    my $outdb = $self->db;
    $outdb = $self->get_dbadaptor($self->pipeline_db) if ($self->pipeline_db);
    my $analysis = $outdb->get_AnalysisAdaptor->fetch_by_logic_name($self->logic_name);
    my $sic = $outdb->get_StateInfoContainer;
    foreach my $stable_id (@{$self->output}) {
        $sic->store_input_id_analysis($stable_id, $analysis, $self->input_id_type);
    }
    return 1;
}



sub genes_db {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_genes_db'} = $value;
    }

    if (exists($self->{'_genes_db'})) {
        return $self->{'_genes_db'};
    } else {
        return undef;
    }
}


sub pipeline_db {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_pipeline_db'} = $value;
    }

    if (exists($self->{'_pipeline_db'})) {
        return $self->{'_pipeline_db'};
    } else {
        return undef;
    }
}

sub genes_logic_name {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_genes_logic_name'} = $value;
    }

    if (exists($self->{'_genes_logic_name'})) {
        return $self->{'_genes_logic_name'};
    } else {
        return undef;
    }
}


sub logic_name {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_logic_name'} = $value;
    }

    if (exists($self->{'_logic_name'})) {
        return $self->{'_logic_name'};
    } else {
        return undef;
    }
}


sub input_id_type {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_input_id_type'} = $value;
    }

    if (exists($self->{'_input_id_type'})) {
        return $self->{'_input_id_type'};
    } else {
        return undef;
    }
}


1;
