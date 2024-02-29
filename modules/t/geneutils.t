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


use_ok('Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils');

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);

my $stable_id = 'ENSPAGT00005008336';
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('pararge_aegeria');
my $db = $multi->get_DBAdaptor('core');
my $original_transcript = $db->get_TranscriptAdaptor->fetch_by_stable_id($stable_id);
my $original_translation = $original_transcript->translation;
my $original_protein = $original_translation->seq;
my $cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $original_transcript->get_all_Exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', $original_protein, 'Same exons, same sequence');
cmp_ok($cds_transcript->translation->start, '==', $original_translation->start, 'Same exons, same start');
cmp_ok($cds_transcript->translation->end, '==', $original_translation->end, 'Same exons, same end');
my @new_exons;
foreach my $exon (@{$original_transcript->get_all_Exons}) {
  my %tmp_exon = %$exon;
  my $new_exon = Bio::EnsEMBL::Exon->new_fast(\%tmp_exon);
  $new_exon->phase(-1);
  $new_exon->end_phase(-1);
  push(@new_exons, $new_exon);
}
$new_exons[0]->start($new_exons[0]->start+1);
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => \@new_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', $original_protein, 'Final exon shorten, same translation');
cmp_ok($cds_transcript->translation->start, '==', $original_translation->start, 'Final exon shorten, same start');
cmp_ok($cds_transcript->translation->end, '==', $original_translation->end, 'Final exon shorten, same end');
my $cds_exons = $original_transcript->get_all_translateable_Exons;
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', $original_protein, 'CDS exons, same translation');
cmp_ok($cds_transcript->translation->start, '==', 1, 'CDS exons, translation start is 1');
cmp_ok($cds_transcript->translation->end, '==', $original_translation->end, 'CDS exons, same translation end');

$cds_exons = $original_transcript->get_all_translateable_Exons;
# Shifting coordinates is ok here because the start and end exons are not the same objects as the originals
$cds_exons->[-1]->end($cds_exons->[-1]->end-1);
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', $original_protein, 'Stop removed, same sequence');
cmp_ok($cds_transcript->translation->start, '==', 1, 'Stop removed, translation start is 1');
cmp_ok($cds_transcript->translation->end, '==', 75, 'Stop removed, end translation shorter');

$cds_exons = $original_transcript->get_all_translateable_Exons;
# Shifting coordinates is ok here because the start and end exons are not the same objects as the originals
$cds_exons->[-1]->end($cds_exons->[-1]->end-4);
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', substr($original_protein, 0, length($original_protein)-1), 'Last amino acid removed frame 1, truncated sequence');
cmp_ok($cds_transcript->translation->start, '==', 1, 'Last amino acid removed frame 1, translation start is 1');
cmp_ok($cds_transcript->translation->end, '==', 72, 'Last amino acid removed frame 1, translation end upstream');

$cds_exons = $original_transcript->get_all_translateable_Exons;
# Shifting coordinates is ok here because the start and end exons are not the same objects as the originals
$cds_exons->[-1]->end($cds_exons->[-1]->end-5);
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', substr($original_protein, 0, length($original_protein)-1), 'Last amino acid removed frame 2, truncated sequence');
cmp_ok($cds_transcript->translation->start, '==', 1, 'Last amino acid removed frame 2, translation start is 1');
cmp_ok($cds_transcript->translation->end, '==', 72, 'Last amino acid removed frame 2, translation end upstream');

$cds_exons = $original_transcript->get_all_translateable_Exons;
# Shifting coordinates is ok here because the start and end exons are not the same objects as the originals
$cds_exons->[0]->start($cds_exons->[0]->start-1);
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', $original_protein, 'First amino acid removed frame 1, truncated sequence');
cmp_ok($cds_transcript->translation->start, '==', 2, 'First amino acid removed frame 1, translation start downstream');
cmp_ok($cds_transcript->translation->end, '==', $original_translation->end, 'First amino acid removed frame 1, same translation end');

$cds_exons = $original_transcript->get_all_translateable_Exons;
# Shifting coordinates is ok here because the start and end exons are not the same objects as the originals
$cds_exons->[0]->start($cds_exons->[0]->start-2);
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', $original_protein, 'First amino acid removed frame 1, truncated sequence');
cmp_ok($cds_transcript->translation->start, '==', 3, 'First amino acid removed frame 1, translation start downstream');
cmp_ok($cds_transcript->translation->end, '==', $original_translation->end, 'First amino acid removed frame 1, same translation end');

$cds_exons = $original_transcript->get_all_translateable_Exons;
# Shifting coordinates is ok here because the start and end exons are not the same objects as the originals
$cds_exons->[0]->start($cds_exons->[0]->start+4);
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', substr($original_protein, 35), 'First amino acid removed frame 1, truncated sequence');
cmp_ok($cds_transcript->translation->start, '==', 66, 'First amino acid removed frame 1, translation start downstream');
cmp_ok($cds_transcript->translation->end, '==', $original_translation->end, 'First amino acid removed frame 1, same translation end');

$cds_exons = $original_transcript->get_all_translateable_Exons;
# Shifting coordinates is ok here because the start and end exons are not the same objects as the originals
$cds_exons->[0]->start($cds_exons->[0]->start+5);
$cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
compute_translation($cds_transcript);
cmp_ok($cds_transcript->translation->seq, 'eq', substr($original_protein, 35), 'First amino acid removed frame 2, truncated sequence');
cmp_ok($cds_transcript->translation->start, '==', 66, 'First amino acid removed frame 2, translation start downstream');
cmp_ok($cds_transcript->translation->end, '==', $original_translation->end, 'First amino acid removed frame 2, same translation end');

done_testing();
