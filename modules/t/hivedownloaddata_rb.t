#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

use File::Spec::Functions qw(catfile);
use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);
use Test::Most;

my $expected_dataflow = [{filename => catfile($ENV{PWD}, 'README')}];

sub cleaning {
  unlink catfile($ENV{PWD}, 'README');
}

use_ok('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData');

standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData', # module
	{ # input param hash
		'target'     => 'ftp://ftp.ensembl.org/pub/release-90/README',
		'download_method'     => 'ftp',
    'md5sum' => '2a3ea6c67b0bf0a1cc72c5b15c73b931',
    'output_dir' => $ENV{PWD},
	},
	[ # list of events to test for (just 1 event in this case)
		[ # start event
			'DATAFLOW', # event to test for (could be WARNING)
			$expected_dataflow, # expected data flowed out
			2 # dataflow branch
		], # end event
	]
);

standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData', # module
	{ # input param hash
		'target'     => 'ftp://ftp.ensembl.org/pub/release-90/README',
		'download_method'     => 'ftp',
    'output_dir' => $ENV{PWD},
	},
	[ # list of events to test for (just 1 event in this case)
		[ # start event
			'DATAFLOW', # event to test for (could be WARNING)
			$expected_dataflow, # expected data flowed out
			2 # dataflow branch
		], # end event
	]
);

done_testing();

cleaning();
