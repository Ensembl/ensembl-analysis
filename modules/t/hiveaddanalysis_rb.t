#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_database_from_registry);

use_ok('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses');

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor('core');
my %target_db = (
  -dbname => $db->dbc->dbname,
  -host   => $db->dbc->host,
  -port   => $db->dbc->port,
  -user   => $db->dbc->user,
  -pass   => $db->dbc->pass,
  -driver => $db->dbc->driver,
);

my $logic_name = 'test_analysis';
my $db_logic_name = 'repeatmask_repbase_sus_scrofa';
my $db_dbid = 10000;
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses', # module
	{ # input param hash
    source_type => 'list',
    analyses => [{'logic_name' => $logic_name, '-module' => 'AddAnalysis'}],
    target_db => \%target_db,
	},
);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
cmp_ok($analysis->logic_name, 'eq', $logic_name, 'Checking analysis present with source_type "list"');

my $core_db = get_database_from_registry('sus_scrofa', 'Core');
my %source_db = (
  -dbname => $core_db->dbc->dbname,
  -host   => $core_db->dbc->host,
  -port   => $core_db->dbc->port,
  -user   => $core_db->dbc->user,
  -pass   => $core_db->dbc->pass,
  -driver => $core_db->dbc->driver,
);

standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses', # module
	{ # input param hash
    source_type => 'db',
    analyses => [$db_logic_name],
    source_db => \%source_db,
    target_db => \%target_db,
    analysis_start => $db_dbid,
	},
);

$analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($db_logic_name);
cmp_ok($analysis->logic_name, 'eq', $db_logic_name, 'Checking analysis present with source_type "db"');
cmp_ok($analysis->dbID, 'eq', $db_dbid, 'Checking dbID when using "analysis_start"');

done_testing();
