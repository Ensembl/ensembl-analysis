# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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


package Bio::EnsEMBL::Analysis::Tools::SoftwareConfigLoad;

use strict;
use warnings;
use JSON;
use File::Basename;
use File::Spec;
use Exporter 'import';
our @EXPORT_OK = qw(get_software_path);

my $config_file = File::Spec->catfile(dirname(__FILE__), 'SoftwareConfig.json');  # Find suitable location

sub get_software_path {
  my ($software_type, $tool) = @_;

  open my $fh, '<', $config_file or die "Could not open config file: $!";
  my $json_text = do { local $/; <$fh> };
  close $fh;

  my $config = decode_json($json_text);
  

  # Validate inputs
  unless ($software_type && exists $config->{software_paths}{$software_type}) {
    die "Software type '$software_type' not found in config. Available types: "
      . join(", ", keys %{ $config->{software_paths} }) . "\n";
  }

  unless (exists $config->{software_paths}{$software_type}{$tool}) {
    die "Tool '$tool' not found for software type '$software_type'. Available tools: "
      . join(", ", keys %{ $config->{software_paths}{$software_type} }) . "\n";
  }

  return $config->{software_paths}{$software_type}{$tool};

}
