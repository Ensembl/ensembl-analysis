#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);

use_ok('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatCoverage');

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

standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatCoverage', # module
	{ # input param hash
    source_db => \%target_db,
    repeat_logic_names => ['RepeatMask'],
    coord_system_version => 'NCBI33',
	},
	[ # list of events to test for (just 1 event in this case)
		[ # start event
			'WARNING', # event to test for (could be WARNING)
			$db->dbc->dbname . "\nAnalyses: RepeatMask\nTotal bases = 62842997\nTotal masked = 504576\t( 0.80% masked)\n", # expected data flowed out
		], # end event
		[ # start event
			'DATAFLOW', # event to test for (could be WARNING)
			{repeat_mask_coverage => 0.802915239704433574}, # expected data flowed out
			2 # dataflow branch
		], # end event
	]
);

done_testing();
