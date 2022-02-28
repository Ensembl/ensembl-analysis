=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::AddStableIds -

=head1 SYNOPSIS

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
It is meant to provide the interface for mapping ESTs to the genome
sequence and writing the results as genes. By the way Exonerate is run
we do not cluster transcripts into genes and only write one transcript per gene.
we then create a dbadaptor for the target database.


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::AddStableIds;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::AddStableIds qw (ADD_STABLEIDS_BY_LOGIC);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->read_and_check_config($ADD_STABLEIDS_BY_LOGIC);

    return $self;
}


sub fetch_input {
    my( $self) = @_;
    if (defined $self->PREFIX) {
        $self->prefix($self->PREFIX);
    }
    else {
        my $meta_adptor = $self->db->get_MetaContainer;
        $self->prefix($meta_adptor->single_value_by_key('species.stable_id_prefix'));
    }
    throw("Could not get a prefix for the stable_id") unless (defined $self->prefix);

    return 1;
}

sub run{
    my ($self) = @_;

    logger_info("Everything is done in the write_output method");
    return 1;
}


sub write_output{
    my ($self,@output) = @_;

    my $outdb = $self->get_dbadaptor($self->GENES_DB);
    my $seq_region_id = $outdb->get_SliceAdaptor->fetch_by_name($self->input_id)->get_seq_region_id;
    my $analysis_id;
    $analysis_id = $outdb->get_AnalysisAdaptor->logic_name($self->LOGIC_NAME) if (defined $self->LOGIC_NAME);

    foreach my $table ('gene', 'transcript', 'exon') {
        my $sql_query = 'UPDATE '.$table.' SET stable_id = CONCAT("'.$self->prefix.'", "'.$self->one_letter($table).'", LPAD('.$table.'_id, 11, 0)) WHERE seq_region_id = '.$seq_region_id;
        $sql_query .= ' AND analysis_id = '.$analysis_id if (defined $analysis_id);
        my $sth = $outdb->dbc->prepare($sql_query);
        $sth->execute();
    }
    my $sql_query = 'UPDATE translation tln, transcript t SET tln.stable_id = CONCAT("'.$self->prefix.'", "'.$self->one_letter('translation').'", LPAD(tln.translation_id, 11, 0)) WHERE t.transcript_id = tln.transcript_id AND t.seq_region_id = '.$seq_region_id;
    $sql_query .= ' AND t.analysis_id = '.$analysis_id if (defined $analysis_id);
    my $sth = $outdb->dbc->prepare($sql_query);
    $sth->execute();
}



sub prefix {
    my ($self, $value) = @_;

    if (defined $value) {
        $self->{prefix} = $value;
    }
    return $self->{prefix};
}


sub one_letter {
    my ($self, $table) = @_;

    if ($table eq 'gene') {
        return 'G';
    }
    elsif ($table eq 'transcript') {
        return 'T';
    }
    elsif ($table eq 'translation') {
        return 'P';
    }
    elsif ($table eq 'exon') {
        return 'E';
    }
    else {
        throw('You should specify a proper MySQL table');
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


sub PREFIX {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_PREFIX'} = $value;
    }

    if (exists($self->{'_CONFIG_PREFIX'})) {
        return $self->{'_CONFIG_PREFIX'};
    } else {
        return undef;
    }
}


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


1;
