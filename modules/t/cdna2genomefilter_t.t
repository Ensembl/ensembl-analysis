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

use Bio::EnsEMBL::Analysis::Tools::Filter::cDNA2GenomeTranscriptFilter;

my %params = (
  -coverage => 50,
  -percent_id => 87,
  -reject_processed_pseudos => 1,
  -best_in_genome => 1,
  -verbosity => 2,
);
my $classic_cigar = '45M432378I786M1D3M987I78M';
my $codon_cigar = '45M432378I78M543C987I435C67M678263I78M';
my $codon_cigar_check = '45M432378I621M987I502M678263I78M';

my $filter = new_ok('Bio::EnsEMBL::Analysis::Tools::Filter::cDNA2GenomeTranscriptFilter');
is($filter->min_coverage, undef, 'Checking default min_coverage');
is($filter->min_percent, undef, 'Checking default min_percent');
cmp_ok($filter->reject_processed_pseudos, '==', 0, 'Checking default reject_processed_pseudos');
cmp_ok($filter->best_in_genome, '==', 0, 'Checking default best_in_genome');
cmp_ok($filter->verbosity, '==', 0, 'Checking default verbosity');

$filter->min_coverage(90);
$filter->min_percent(97);
$filter->reject_processed_pseudos(1);
$filter->best_in_genome(1);
$filter->verbosity(1);
cmp_ok($filter->min_coverage, '==', 90, 'Checking min_coverage');
cmp_ok($filter->min_percent, '==', 97, 'Checking min_percent');
cmp_ok($filter->reject_processed_pseudos, '==', 1, 'Checking reject_processed_pseudos');
cmp_ok($filter->best_in_genome, '==', 1, 'Checking best_in_genome');
cmp_ok($filter->verbosity, '==', 1, 'Checking verbosity');

$filter = Bio::EnsEMBL::Analysis::Tools::Filter::cDNA2GenomeTranscriptFilter->new(%params);
cmp_ok($filter->min_coverage, '==', 50, 'Checking min_coverage');
cmp_ok($filter->min_percent, '==', 87, 'Checking min_percent');
cmp_ok($filter->reject_processed_pseudos, '==', 1, 'Checking reject_processed_pseudos');
cmp_ok($filter->best_in_genome, '==', 1, 'Checking best_in_genome');
cmp_ok($filter->verbosity, '==', 2, 'Checking verbosity');

eval{
  $filter->filter_results;
};
like($@, qr/You should give an arrayref of Bio::EnsEMBL::Transcript/, 'Checking fails on empty');

cmp_ok($filter->_update_cigar_string($classic_cigar), 'eq', $classic_cigar, 'Check no substitution on normal cigar');
cmp_ok($filter->_update_cigar_string($codon_cigar), 'eq', $codon_cigar_check, 'Check substitution on cdna2genome cigar');
done_testing();
