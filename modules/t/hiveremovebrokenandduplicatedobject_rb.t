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

use Test::More;
#use Test::Warnings;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

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

# First checking that it doesn't delete any genes if there is nothing to remove
my $num_genes = count_rows($db, 'gene');
my $num_transcripts = count_rows($db, 'transcript');
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects',
	{
    target_db => \%target_db,
	},
);
cmp_ok(count_rows($db, 'gene'), '==', $num_genes, 'Checking no change done in gene');
cmp_ok(count_rows($db, 'transcript'), '==', $num_transcripts, 'Checking no change done in transcript');
$multi->save('core', 'gene', 'transcript', 'translation', 'exon', 'exon_transcript', 'dna_align_feature', 'protein_align_feature', 'transcript_supporting_feature', 'supporting_feature');

# I will delete a transcript to pretend a broken gene
my $ta = $db->get_TranscriptAdaptor;
my $t = $ta->fetch_by_stable_id('ENST00000262651');
$ta->remove($t);
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects',
	{
    target_db => \%target_db,
	},
);
cmp_ok(count_rows($db, 'gene'), '==', $num_genes-1, 'Checking deletion in gene');
cmp_ok(count_rows($db, 'transcript'), '==', $num_transcripts-1, 'Checking deletion in transcript');

# Now we reload and delete the transcript to check parameter 'iid'
$multi->restore;
$t = $ta->fetch_by_stable_id('ENST00000262651');
$ta->remove($t);

my $num_genes_mt = count_rows($db, 'gene', 'WHERE seq_region_id = 965899');
my $num_transcripts_mt = count_rows($db, 'transcript', 'WHERE seq_region_id = 965899');
my $num_genes_20 = count_rows($db, 'gene', 'WHERE seq_region_id = 469283');
my $num_transcripts_20 = count_rows($db, 'transcript', 'WHERE seq_region_id = 469283');

# First testing on a region without problem
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects',
	{
    target_db => \%target_db,
    iid => 'chromosome:NCBI33:MT_NC_001807:1:16571:1',
	},
);
cmp_ok(count_rows($db, 'gene', 'WHERE seq_region_id = 965899'), '==', $num_genes_mt, 'Checking no change in gene on MT_NC_001807');
cmp_ok(count_rows($db, 'transcript', 'WHERE seq_region_id = 965899'), '==', $num_transcripts_mt, 'Checking no change in transcript on MT_NC_001807');
cmp_ok(count_rows($db, 'gene', 'WHERE seq_region_id = 469283'), '==', $num_genes_20, 'Checking no change in gene on 20');
cmp_ok(count_rows($db, 'transcript', 'WHERE seq_region_id = 469283'), '==', $num_transcripts_20, 'Checking no change in transcript on 20');

# Now testing on a region with a problem
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects',
	{
    target_db => \%target_db,
    iid => 'chromosome:NCBI33:20:1:62842997:1',
	},
);
cmp_ok(count_rows($db, 'gene', 'WHERE seq_region_id = 965899'), '==', $num_genes_mt, 'Checking no change in gene on MT_NC_001807');
cmp_ok(count_rows($db, 'transcript', 'WHERE seq_region_id = 965899'), '==', $num_transcripts_mt, 'Checking no change in transcript on MT_NC_001807');
cmp_ok(count_rows($db, 'gene', 'WHERE seq_region_id = 469283'), '==', $num_genes_20-1, 'Checking deletion in gene on 20');
cmp_ok(count_rows($db, 'transcript', 'WHERE seq_region_id = 469283'), '==', $num_transcripts_20, 'Checking no change in transcript on 20');

# Now we reload and delete the transcript to check parameter 'iid'
$multi->restore;
$num_genes = count_rows($db, 'gene', 'WHERE seq_region_id = 469283');
$num_transcripts = count_rows($db, 'transcript', 'WHERE seq_region_id = 469283');
my $num_exons = count_rows($db, 'exon', 'WHERE seq_region_id = 469283');
my $g = $db->get_GeneAdaptor->fetch_by_stable_id('ENSG00000174873');
empty_Gene($g);
$db->get_GeneAdaptor->store($g);
cmp_ok(count_rows($db, 'gene', 'WHERE seq_region_id = 469283'), '==', $num_genes+1, 'Checking added gene on 20');
cmp_ok(count_rows($db, 'transcript', 'WHERE seq_region_id = 469283'), '==', $num_transcripts+1, 'Checking added transcript on 20');
cmp_ok(count_rows($db, 'exon', 'WHERE seq_region_id = 469283'), '==', $num_exons+6, 'Checking added exons on 20');
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects',
	{
    target_db => \%target_db,
	},
);
cmp_ok(count_rows($db, 'gene', 'WHERE seq_region_id = 469283'), '==', $num_genes, 'Checking removed gene on 20');
cmp_ok(count_rows($db, 'transcript', 'WHERE seq_region_id = 469283'), '==', $num_transcripts, 'Checking removed transcript on 20');
cmp_ok(count_rows($db, 'exon', 'WHERE seq_region_id = 469283'), '==', $num_exons, 'Checking added exons on 20');

done_testing();
