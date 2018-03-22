#!/usr/env perl
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

use Test::More;

use Bio::EnsEMBL::Analysis::Tools::Filter;

my %params = (
  -coverage => 50,
  -percent_id => 87,
  -reject_processed_pseudos => 1,
  -best_in_genome => 1,
  -verbosity => 2,
);

my $filter = new_ok('Bio::EnsEMBL::Analysis::Tools::Filter');
ok(!defined($filter->min_coverage), 'Checking default min_coverage');
ok(!defined($filter->min_percent), 'Checking default min_percent');
ok($filter->reject_processed_pseudos == 0, 'Checking default reject_processed_pseudos');
ok($filter->best_in_genome == 0, 'Checking default best_in_genome');
ok($filter->verbosity == 0, 'Checking default verbosity');

$filter->min_coverage(90);
$filter->min_percent(97);
$filter->reject_processed_pseudos(1);
$filter->best_in_genome(1);
$filter->verbosity(1);
ok($filter->min_coverage == 90, 'Checking min_coverage');
ok($filter->min_percent == 97, 'Checking min_percent');
ok($filter->reject_processed_pseudos == 1, 'Checking reject_processed_pseudos');
ok($filter->best_in_genome == 1, 'Checking best_in_genome');
ok($filter->verbosity == 1, 'Checking verbosity');

$filter = Bio::EnsEMBL::Analysis::Tools::Filter->new(%params);
ok($filter->min_coverage == 50, 'Checking min_coverage');
ok($filter->min_percent == 87, 'Checking min_percent');
ok($filter->reject_processed_pseudos == 1, 'Checking reject_processed_pseudos');
ok($filter->best_in_genome == 1, 'Checking best_in_genome');
ok($filter->verbosity == 2, 'Checking verbosity');

eval{
  $filter->filter_results;
};
ok($@ && $@ =~ /You should give an arrayref of objects/, 'Checking fails on empty');
eval{
  $filter->filter_results([]);
};
ok($@ && $@ =~ /You should implement the filter_results method/, 'Checking fails on not implemented');
done_testing();
