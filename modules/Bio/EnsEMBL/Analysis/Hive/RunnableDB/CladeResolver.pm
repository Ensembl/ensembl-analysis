# Copyright [2018-2020] EMBL-European Bioinformatics Institute
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

package  Bio::EnsEMBL::Analysis::Hive::RunnableDB::CladeResolver;

use strict;
use warnings;
use JSON;
use Carp qw(croak);
use Exporter 'import';
our @EXPORT_OK = qw(get_clade_from_taxon_id);

my %taxonomy_rank_order = (
    'species' => 1, 'genus' => 2, 'family' => 3,
    'order' => 4, 'class' => 5, 'phylum' => 6, 'kingdom' => 7
);

sub get_clade_from_taxon_id {
    my ($registry_dba, $taxon_id, $clade_json_path) = @_;
    croak "Missing registry_dba" unless $registry_dba;
    croak "Missing taxon_id" unless defined $taxon_id;
    croak "Missing JSON file path" unless $clade_json_path;

    # Load JSON
    open(my $fh, '<', $clade_json_path) or croak "Cannot open $clade_json_path: $!";
    local $/;
    my $json_text = <$fh>;
    close($fh);
    my $clade_config = decode_json($json_text);

    # Query taxonomy hierarchy
    my $sql = qq{
        SELECT taxon_class_id, taxon_class
        FROM taxonomy
        WHERE lowest_taxon_id = ?
    };
    my $sth = $registry_dba->dbc->prepare($sql);
    $sth->execute($taxon_id);

    my @lineage;
    while (my ($class_id, $rank) = $sth->fetchrow_array) {
        next unless exists $taxonomy_rank_order{$rank};
        push @lineage, { taxon_id => $class_id, rank => $rank };
    }

    @lineage = sort {
        $taxonomy_rank_order{$a->{rank}} <=> $taxonomy_rank_order{$b->{rank}}
    } @lineage;

    # Match to JSON
    foreach my $entry (@lineage) {
        my $tid = $entry->{taxon_id};
        foreach my $clade_name (keys %$clade_config) {
            return $clade_name if $clade_config->{$clade_name}->{taxon_id} == $tid;
        }
    }

    croak "No clade match found in JSON for taxon_id $taxon_id";
}

1;