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

use Test::More;

#use Bio::EnsEMBL::Test::TestUtils;
#use Bio::EnsEMBL::Test::MultiTestDB;
#
#use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);
#
#use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_database_from_registry);

use_ok('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition');

TODO: {
  local $TODO = 'Proper tests needed';
  note($TODO);
}

####
# The code and data below is the test case which was in the module. A proper test needs to be done
####
#
#my %source_db = (
#);
#my %dna_db = (
#);
#my %target_db = (
#);
#
#standaloneJob(
#	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition', # module
#	{ # input param hash
#    dna_db => \%dna_db,
#    donor_dbs => $self->o('utr_donor_dbs'),
#    acceptor_dbs => $self->o('utr_acceptor_dbs'),
#    utr_biotype_priorities => $self->o('utr_biotype_priorities'),
#    target_db => \%target_db,
#    iid => $input_id,
#    iid_type => 'slice',
#	},
#);
#
#
#
#sub donor_test_cases {
#  my ($self) = @_;
#
#  my $dba = $self->hrdb_get_con('target_db');
#  my $slice_adaptor = $dba->get_SliceAdaptor();
#  my $slice = $slice_adaptor->fetch_by_name('chromosome:Mmul_8.0.1:9:1:129882849:1');
#
#  my $exon_1 = new Bio::EnsEMBL::Exon(
#                                       -START     => 100,
#                                       -END       => 199,
#                                       -STRAND    => 1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_2 = new Bio::EnsEMBL::Exon(
#                                       -START     => 300,
#                                       -END       => 399,
#                                       -STRAND    => 1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_3 = new Bio::EnsEMBL::Exon(
#                                       -START     => 500,
#                                       -END       => 599,
#                                       -STRAND    => 1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_4 = new Bio::EnsEMBL::Exon(
#                                       -START     => 700,
#                                       -END       => 799,
#                                       -STRAND    => 1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_5 = new Bio::EnsEMBL::Exon(
#                                       -START     => 900,
#                                       -END       => 999,
#                                       -STRAND    => 1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my  @exons_set_1 = ($exon_1,$exon_2,$exon_3,$exon_4,$exon_5);
#
#  my $transcript_1 = new Bio::EnsEMBL::Transcript(
#                                                   -EXONS => \@exons_set_1,
#                                                   -STRAND    => 1,
#                                                   -SLICE     => $slice,
#                                                   -ANALYSIS  => $self->analysis);
#
#
#  my $exon_6 = new Bio::EnsEMBL::Exon(
#                                       -START     => 900,
#                                       -END       => 999,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_7 = new Bio::EnsEMBL::Exon(
#                                       -START     => 1100,
#                                       -END       => 1199,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_8 = new Bio::EnsEMBL::Exon(
#                                       -START     => 1300,
#                                       -END       => 1399,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_9 = new Bio::EnsEMBL::Exon(
#                                       -START     => 1500,
#                                       -END       => 1599,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_10 = new Bio::EnsEMBL::Exon(
#                                       -START     => 1700,
#                                       -END       => 1799,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_11 = new Bio::EnsEMBL::Exon(
#                                       -START     => 1900,
#                                       -END       => 1999,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_12 = new Bio::EnsEMBL::Exon(
#                                       -START     => 2100,
#                                       -END       => 2199,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my $exon_13 = new Bio::EnsEMBL::Exon(
#                                       -START     => 2300,
#                                       -END       => 2399,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => -1);
#
#  my  @exons_set_2 = ($exon_12,$exon_11,$exon_10,$exon_9,$exon_8,$exon_7,$exon_6);
#
#  my $transcript_2 = new Bio::EnsEMBL::Transcript(
#                                                   -EXONS => \@exons_set_2,
#                                                   -STRAND    => -1,
#                                                   -SLICE     => $slice,
#                                                   -ANALYSIS  => $self->analysis);
#
#  my  @exons_set_3 = ($exon_13,$exon_12,$exon_11,$exon_10,$exon_9,$exon_8,$exon_7,$exon_6);
#
#  my $transcript_3 = new Bio::EnsEMBL::Transcript(
#                                                   -EXONS => \@exons_set_3,
#                                                   -STRAND    => -1,
#                                                   -SLICE     => $slice,
#                                                   -ANALYSIS  => $self->analysis);
#
#  $transcript_1->biotype('cdna');
#  $transcript_2->biotype('cdna');
#  $transcript_3->biotype('cdna_predicted');
#
#  say "Created the following test donor transcripts: ";
#  say "DONOR T1: (".$transcript_1->start.":".$transcript_1->end.":".$transcript_1->strand.")";
#  my $exons = $transcript_1->get_all_Exons();
#  foreach my $exon (@{$exons}) {
#    print "(".$exon->start."..".$exon->end.")";
#  }
#  print "\n";
#
#  say "DONOR T2: (".$transcript_2->start.":".$transcript_2->end.":".$transcript_2->strand.")";
#  $exons = $transcript_2->get_all_Exons();
#  foreach my $exon (@{$exons}) {
#    print "(".$exon->start."..".$exon->end.")";
#  }
#  print "\n";
#
#  say "DONOR T3: (".$transcript_3->start.":".$transcript_3->end.":".$transcript_3->strand.")";
#  $exons = $transcript_3->get_all_Exons();
#  foreach my $exon (@{$exons}) {
#    print "(".$exon->start."..".$exon->end.")";
#  }
#  print "\n";
#
#  return([$transcript_1,$transcript_2,$transcript_3]);
#}
#
#sub acceptor_test_cases {
#  my ($self) = @_;
#
#  my $dba = $self->hrdb_get_con('target_db');
#  my $slice_adaptor = $dba->get_SliceAdaptor();
#  my $slice = $slice_adaptor->fetch_by_name('chromosome:Mmul_8.0.1:9:1:129882849:1');
#
#  my $exon_1 = new Bio::EnsEMBL::Exon(
#                                       -START     => 550,
#                                       -END       => 599,
#                                       -STRAND    => 1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => 0,
#                                       -END_PHASE => 0);
#
#  my $exon_2 = new Bio::EnsEMBL::Exon(
#                                      -START     => 700,
#                                      -END       => 759,
#                                      -STRAND    => 1,
#                                      -SLICE     => $slice,
#                                      -ANALYSIS  => $self->analysis,
#                                      -PHASE     => 0,
#                                      -END_PHASE => 0);
#
#  my $translation_1 = new Bio::EnsEMBL::Translation(
#                                                     -START_EXON => $exon_1,
#                                                     -END_EXON => $exon_2,
#                                                     -SEQ_START => 1,
#                                                     -SEQ_END => 49,
#                                                   );
#  my  @exons_set_1 = ($exon_1,$exon_2);
#
#  my $transcript_1 = new Bio::EnsEMBL::Transcript( -DBID  => 1,
#                                                   -EXONS => \@exons_set_1,
#                                                   -STRAND    => 1,
#                                                   -SLICE     => $slice,
#                                                   -ANALYSIS  => $self->analysis);
#
#  $transcript_1->translation($translation_1);
#
#  my $exon_3 = new Bio::EnsEMBL::Exon(
#                                       -START     => 1350,
#                                       -END       => 1399,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => 0,
#                                       -END_PHASE => 0);
#
#  my $exon_4 = new Bio::EnsEMBL::Exon(
#                                      -START     => 1500,
#                                      -END       => 1599,
#                                      -STRAND    => -1,
#                                      -SLICE     => $slice,
#                                      -ANALYSIS  => $self->analysis,
#                                      -PHASE     => 0,
#                                      -END_PHASE => 0);
#
#  my $exon_5 = new Bio::EnsEMBL::Exon(
#                                      -START     => 1700,
#                                      -END       => 1759,
#                                      -STRAND    => -1,
#                                      -SLICE     => $slice,
#                                      -ANALYSIS  => $self->analysis,
#                                      -PHASE     => 0,
#                                      -END_PHASE => 0);
#
#  my $translation_2 = new Bio::EnsEMBL::Translation(
#                                                     -START_EXON => $exon_5,
#                                                     -END_EXON => $exon_3,
#                                                     -SEQ_START => 1,
#                                                     -SEQ_END => 49,
#                                                   );
#  my  @exons_set_2 = ($exon_5,$exon_4,$exon_3);
#
#  my $transcript_2 = new Bio::EnsEMBL::Transcript( -DBID  => 2,
#                                                   -EXONS => \@exons_set_2,
#                                                   -STRAND    => -1,
#                                                   -SLICE     => $slice,
#                                                   -ANALYSIS  => $self->analysis);
#
#  $transcript_2->translation($translation_2);
#
#
#  my $exon_6 = new Bio::EnsEMBL::Exon(
#                                       -START     => 1910,
#                                       -END       => 1989,
#                                       -STRAND    => -1,
#                                       -SLICE     => $slice,
#                                       -ANALYSIS  => $self->analysis,
#                                       -PHASE     => 0,
#                                       -END_PHASE => 0);
#
#  my $translation_3 = new Bio::EnsEMBL::Translation(
#                                                     -START_EXON => $exon_6,
#                                                     -END_EXON => $exon_6,
#                                                     -SEQ_START => 1,
#                                                     -SEQ_END => 49,
#                                                   );
#
#  my  @exons_set_3 = ($exon_6);
#
#  my $transcript_3 = new Bio::EnsEMBL::Transcript( -DBID  => 3,
#                                                   -EXONS => \@exons_set_3,
#                                                   -STRAND    => -1,
#                                                   -SLICE     => $slice,
#                                                   -ANALYSIS  => $self->analysis);
#
#  $transcript_3->translation($translation_3);
#
#  say "ACCEPTOR T1: (".$transcript_1->start.":".$transcript_1->end.":".$transcript_1->strand.")";
#  my $exons = $transcript_1->get_all_Exons();
#  foreach my $exon (@{$exons}) {
#    print "(".$exon->start."..".$exon->end.")";
#  }
#  print "\n";
#  say "ACCEPTOR TN1: ".$transcript_1->translation()->seq();
#
#
#  say "ACCEPTOR T2: (".$transcript_2->start.":".$transcript_2->end.":".$transcript_2->strand.")";
#  $exons = $transcript_2->get_all_Exons();
#  foreach my $exon (@{$exons}) {
#    print "(".$exon->start."..".$exon->end.")";
#  }
#  print "\n";
#  say "ACCEPTOR TN2: ".$transcript_2->translation()->seq();
#
#
#  say "ACCEPTOR T3: (".$transcript_3->start.":".$transcript_3->end.":".$transcript_3->strand.")";
#  $exons = $transcript_3->get_all_Exons();
#  foreach my $exon (@{$exons}) {
#    print "(".$exon->start."..".$exon->end.")";
#  }
#  print "\n";
#  say "ACCEPTOR TN3: ".$transcript_3->translation()->seq();
#
#
#  return([$transcript_1,$transcript_2,$transcript_3]);
#}
done_testing();
