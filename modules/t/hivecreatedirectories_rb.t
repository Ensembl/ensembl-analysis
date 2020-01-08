#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
use strict;
use warnings;

use File::Path qw(remove_tree);
use File::Spec::Functions qw(catdir);

use Test::More;

use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);

use_ok('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDirectories');

my $directory = catdir($ENV{PWD}, 'test_directory');
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDirectories', # module
	{ # input param hash
    paths => [$directory],
	},
);

my @stat = stat($directory);
cmp_ok(sprintf("%04o", $stat[2] & 07777), 'eq', 2775, 'Checking permissions for default production directory');
remove_tree($directory);

my $directory_2 = catdir($directory, 'directory_test');
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDirectories', # module
	{ # input param hash
    paths => [{
      path => $directory_2,
      mode => 0755,
      }],
	},
);

@stat = stat($directory);
cmp_ok(sprintf("%04o", $stat[2] & 07777), 'eq', '0755', 'Checking permissions for directory');
@stat = stat($directory_2);
cmp_ok(sprintf("%04o", $stat[2] & 07777), 'eq', '0755', 'Checking permissions for subdirectory');
remove_tree($directory);

done_testing();
