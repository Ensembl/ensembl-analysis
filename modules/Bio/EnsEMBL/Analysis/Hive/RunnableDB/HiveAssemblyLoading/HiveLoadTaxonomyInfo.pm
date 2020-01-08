=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo;

use strict;
use warnings;

use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Returns the defaults parameters
               ensembl_release => $ENV{ENSEMBL_RELEASE},
 Returntype : None
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    ensembl_release => $ENV{ENSEMBL_RELEASE},
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch the taxonomy id from the aprameter 'taxon_id' or from the 'target_db'
              meta table. It will check that species.scientific_name and species.common_name
              are set (it does not check their correctness).
              The analysis can be skipped with the 'skip_analysis' flag.
 Returntype : None
 Exceptions : Throws if 'taxonomy_db' is not set
              Throws if 'target_db' is not set

=cut

sub fetch_input {
  my $self = shift;

  if ($self->param('skip_analysis')) {
    $self->complete_early('Skipping the analysis');
  }
  else {
    my $target_db = $self->get_database_by_name('target_db');
    if (!($self->param_is_defined('taxon_id') and defined $self->param('taxon_id'))) {
      my $taxon_id = $target_db->get_MetaContainer->get_taxonomy_id;
      if ($taxon_id) {
        $self->param('taxon_id', $taxon_id);
      }
      else {
        $self->throw('No taxonomy id found in the Core database or in the parameter "taxon_id"');
      }
    }
    my $value = $target_db->get_MetaContainer->get_scientific_name;
    if ($value) {
      $self->param('scientific_name', $value);
    }
    $value = $target_db->get_MetaContainer->get_common_name;
    if ($value) {
      $self->param('common_name', $value);
    }
    $self->say_with_header('Looking for taxonomy data for '.$self->param('taxon_id'));
    my $taxonomy_db = $self->get_database_by_name('taxonomy_db', undef, 'Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor');
    $self->say_with_header('Using '.$taxonomy_db->dbc->dbname.'@'.$taxonomy_db->dbc->host);
    $self->hrdb_set_con($target_db, 'target_db');
    $self->hrdb_set_con($taxonomy_db, 'taxonomy_db');
  }
}

=head2 run

 Arg [1]    : None
 Description: Fetch the classification information from the 'taxonomy_db' database.
              If 'species.scientific_name' or 'species.common_name' are not set, we
              will attempt to fetch them.
 Returntype : None
 Exceptions : None

=cut

sub run {
  my $self = shift;

  my $node_adaptor = $self->hrdb_get_con('taxonomy_db')->get_TaxonomyNodeAdaptor;
  my $node = $node_adaptor->fetch_by_taxon_id($self->param('taxon_id'));
  if (!$self->param_is_defined('scientific_name')) {
    my $scientific_name = $node->name('scientific name');
    if ($scientific_name) {
      $self->say_with_header($scientific_name);
      $self->output([['species.scientific_name', $scientific_name]]);
    }
  }
  if (!$self->param_is_defined('common_name')) {
    my $common_name = $node->name('genbank common name');
    if ($common_name) {
      $self->say_with_header($common_name);
      $self->output([['species.common_name', $common_name]]);
    }
  }
  foreach my $ancestor ( @{ $node_adaptor->fetch_ancestors($node)}){
    next if ($ancestor->rank eq 'genus');
    $self->say_with_header($ancestor->name);
    $self->output([['species.classification', $ancestor->name]]);
    last if ($ancestor->rank eq 'superkingdom');
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Write the classification in the meta table
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

  my $meta_adaptor = $self->hrdb_get_con('target_db')->get_MetaContainer;
  foreach my $data (@{$self->output}) {
    $meta_adaptor->store_key_value(@$data);
  }
}

1;
