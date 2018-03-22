#!/usr/bin/env perl
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

use strict;
use warnings;

use Test::More;
use File::Spec::Functions qw(devnull);

use Bio::EnsEMBL::Test::TestUtils;

use_ok('Bio::EnsEMBL::Analysis::Tools::Utilities');

my $db = Bio::EnsEMBL::Analysis::Tools::Utilities::get_database_from_registry('sus_scrofa', 'Core');
cmp_ok(ref($db), 'eq', 'Bio::EnsEMBL::DBSQL::DBAdaptor', 'Checking object type for "get_database_from_registry"');
cmp_ok($db->dbc->host, 'eq', 'ensembldb.ensembl.org', 'Checking host name for "get_database_from_registry"');
cmp_ok($db->dbc->port, '==', 3306, 'Checking port for "get_database_from_registry"');
cmp_ok($db->dbc->user, 'eq', 'anonymous', 'Checking user name for "get_database_from_registry"');
cmp_ok($db->dbc->dbname, '=~', 'sus_scrofa_core_', 'Checking dbname for "get_database_from_registry"');

cmp_ok(Bio::EnsEMBL::Analysis::Tools::Utilities::is_slice_name('chromosome:GRCh38:1:1:199099:1'), '==', 1, 'Cheking "is_slice_name" is correct');
cmp_ok(Bio::EnsEMBL::Analysis::Tools::Utilities::is_slice_name('chromosome:GRCh38:1:1:199099:1:A7S8A7.1:1'), '==', 0, 'Cheking "is_slice_name" is wrong');
cmp_ok(Bio::EnsEMBL::Analysis::Tools::Utilities::is_slice_name('chromosome:GRCh38::1:199099:1'), '==', 0, 'Cheking "is_slice_name" is wrong');


cmp_ok(Bio::EnsEMBL::Analysis::Tools::Utilities::parse_timer('2m'), '==', 120, 'Cheking "parse_timer" with 2m');
cmp_ok(Bio::EnsEMBL::Analysis::Tools::Utilities::parse_timer('1h2m'), '==', 3720, 'Cheking "parse_timer" with 1h2m');
cmp_ok(Bio::EnsEMBL::Analysis::Tools::Utilities::parse_timer('1h2m3'), '==', 3723, 'Cheking "parse_timer" with 1h2m3');
cmp_ok(Bio::EnsEMBL::Analysis::Tools::Utilities::parse_timer('54'), '==', 54, 'Cheking "parse_timer" with 54 seconds');

cmp_ok(Bio::EnsEMBL::Analysis::Tools::Utilities::execute_with_wait('sleep 2'), '==', 1, 'Cheking "execute_with_wait" with sleep');
my $time = time;
eval {
  Bio::EnsEMBL::Analysis::Tools::Utilities::execute_with_wait('sleep t > '.devnull.' 2>&1', undef, 2);
};
$time = time()-$time;
cmp_ok($time, '>=', 2, 'Checking "execute_with_wait" fails and wait using the time in parameters');
cmp_ok($time, '<', 10, 'Checking "execute_with_wait" fails and wait using the time in parameters');
$time = time;
eval {
  Bio::EnsEMBL::Analysis::Tools::Utilities::execute_with_wait('sleep t > '.devnull.' 2>&1', 'Doh', 2);
};
$time = time()-$time;
cmp_ok($time, '>=', 2, 'Checking "execute_with_wait" fails and wait');
cmp_ok($time, '<', 10, 'Checking "execute_with_wait" fails and wait using the time and message in parameters');
cmp_ok($@, '=~', 'Doh', 'Checking "execute_with_wait" fails and wait using the time and message in parameters');

done_testing();
