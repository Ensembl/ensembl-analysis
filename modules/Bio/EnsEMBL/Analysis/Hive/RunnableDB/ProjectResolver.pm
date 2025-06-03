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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProjectResolver;

use strict;
use warnings;
use JSON;
use Exporter 'import';
use Carp;
use File::Slurp;

our @EXPORT_OK = qw(get_bioproject_group);

# Path to the JSON file
my $json_path = catfile( $self->o('enscode_root_dir'),'ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Config/bioproject_mapping.json');

# Internal storage for the parsed data
my $bioproject_data;

# Load JSON only once
BEGIN {
    eval {
        my $json_text = read_file($json_path);
        $bioproject_data = decode_json($json_text);
    };
    if ($@) {
        die "Failed to load BioProject JSON from $json_path: $@";
    }
}

# Lookup BioProject ID and return group name
sub get_bioproject_group {
    my ($bioproject_id) = @_;

    croak "BioProject ID is required" unless defined $bioproject_id;
    return $bioproject_data->{$bioproject_id};
}

1;