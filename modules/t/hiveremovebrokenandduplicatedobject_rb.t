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

use File::Spec::Functions qw(catfile);

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

my $standalone_path = 'standaloneJob.pl';
if (-e catfile('ensembl-hive', 'scripts', $standalone_path)) {
  $standalone_path = catfile('ensembl-hive', 'scripts', $standalone_path);
}

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor('core');

# First checking that it doesn't delete any genes if there is nothing to remove
my $num_genes = count_rows($db, 'gene');
my $num_transcripts = count_rows($db, 'transcript');
system($standalone_path.' Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects -target_db \'{"-host" => "'.$db->dbc->host.'", "-user" => "'.$db->dbc->user.'", "-port" => "'.$db->dbc->port.'", "-pass" => "'.$db->dbc->password.'", "-dbname" => "'.$db->dbc->dbname.'"}\'');
ok($? == 0, 'Checking the job ran successfully');
cmp_ok(count_rows($db, 'gene'), '==', $num_genes, 'Checking no change done in gene');
cmp_ok(count_rows($db, 'transcript'), '==', $num_transcripts, 'Checking no change done in transcript');
$multi->save('core', 'gene', 'transcript', 'translation', 'exon', 'exon_transcript', 'dna_align_feature', 'protein_align_feature', 'transcript_supporting_feature', 'supporting_feature');

# I will delete a transcript to pretend a broken gene
my $ta = $db->get_TranscriptAdaptor;
my $t = $ta->fetch_by_stable_id('ENST00000262651');
$ta->remove($t);
system($standalone_path.' Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects -target_db \'{"-host" => "'.$db->dbc->host.'", "-user" => "'.$db->dbc->user.'", "-port" => "'.$db->dbc->port.'", "-pass" => "'.$db->dbc->password.'", "-dbname" => "'.$db->dbc->dbname.'"}\'');
ok($? == 0, 'Checking the job ran successfully');
cmp_ok(count_rows($db, 'gene'), '==', $num_genes-1, 'Checking deletion in gene');
cmp_ok(count_rows($db, 'transcript'), '==', $num_transcripts-1, 'Checking deletion in transcript');

# Now we reload and delete the transcript to check parameter 'iid'
$multi->restore;
$t = $ta->fetch_by_stable_id('ENST00000262651');
$ta->remove($t);

# First testing on a region without problem
$num_genes = count_rows($db, 'gene', 'WHERE seq_region_id = 965899');
$num_transcripts = count_rows($db, 'transcript', 'WHERE seq_region_id = 965899');
system($standalone_path.' Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects -target_db \'{"-host" => "'.$db->dbc->host.'", "-user" => "'.$db->dbc->user.'", "-port" => "'.$db->dbc->port.'", "-pass" => "'.$db->dbc->password.'", "-dbname" => "'.$db->dbc->dbname.'"}\' -iid chromosome:NCBI33:MT_NC_001807:1:16571:1');
ok($? == 0, 'Checking the job ran successfully');
cmp_ok(count_rows($db, 'gene', 'WHERE seq_region_id = 965899'), '==', $num_genes, 'Checking no change in gene on MT_NC_001807');
cmp_ok(count_rows($db, 'transcript', 'WHERE seq_region_id = 965899'), '==', $num_transcripts, 'Checking no change in transcript on MT_NC_001807');

# Now testing on a region with a problem
$num_genes = count_rows($db, 'gene', 'WHERE seq_region_id = 469283');
$num_transcripts = count_rows($db, 'transcript', 'WHERE seq_region_id = 469283');
system($standalone_path.' Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects -target_db \'{"-host" => "'.$db->dbc->host.'", "-user" => "'.$db->dbc->user.'", "-port" => "'.$db->dbc->port.'", "-pass" => "'.$db->dbc->password.'", "-dbname" => "'.$db->dbc->dbname.'"}\' -iid chromosome:NCBI33:20:1:62842997:1');
ok($? == 0, 'Checking the job ran successfully');
cmp_ok(count_rows($db, 'gene', 'WHERE seq_region_id = 469283'), '==', $num_genes-1, 'Checking no change in gene on 20');
cmp_ok(count_rows($db, 'transcript', 'WHERE seq_region_id = 469283'), '==', $num_transcripts, 'Checking no change in transcript on 20');

done_testing();
