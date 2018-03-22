#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use File::Spec::Functions qw(catfile);
use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);
use Test::Most;
use Net::FTP;

my $version = 0;
my $client = Net::FTP->new('ftp.ebi.ac.uk');
$client->login;
$client->cwd('pub/databases/embl/release/doc') || die('could not go into some dir');
foreach my $file ($client->ls) {
  if ($file =~ /Release_(\d+)/) {
    $version = $1;
    last;
  }
}
$client->quit;
die('Cannot get ENA version, cannot test the module ') unless ($version);


my $base_url = 'ftp://ftp.ebi.ac.uk/pub/databases/embl';
my @expected_dataflow;
foreach my $index (1..9) {
  push(@expected_dataflow, {url => $base_url.'/release/std/rel_std_hum_0'.$index."_r$version.dat.gz"});
}
push(@expected_dataflow, {url => $base_url."/release/std/rel_htc_hum_01_r$version.dat.gz"});

use_ok('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFTPURLs');

standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFTPURLs',
	{
		'base_url'     => $base_url,
		'file_list'     => ['release/std/rel_std_hum_0*', 'release/std/rel_htc_hum_0*'],
	},
	[ # list of events to test for (just 1 event in this case)
		[ # start event
			'DATAFLOW', # event to test for (could be WARNING)
			\@expected_dataflow, # expected data flowed out
			2 # dataflow branch
		], # end event
	]
);

done_testing();
