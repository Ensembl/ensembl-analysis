use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../";
use Test::More;

use Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic;
use Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic;

my $layer_config = Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic->new
  ->_master_config('zebrafish_basic');

is(scalar(@$layer_config), 10, 'zebrafish has ten layers');
is_deeply(
  [map { $_->{ID} } @$layer_config],
  [map { "LAYER$_" } 1 .. 10],
  'layer IDs are unique and ascending',
);

for my $index (0 .. $#$layer_config) {
  my $expected = [map { "LAYER$_" } 1 .. $index];
  is_deeply(
    $layer_config->[$index]{FILTER_AGAINST} || [],
    $expected,
    "$layer_config->[$index]{ID} filters against all earlier layers",
  );
}

my %expected_projection_layer = (
  projection_1          => 'LAYER2',
  projection_2          => 'LAYER3',
  projection_3          => 'LAYER3',
  projection_4          => 'LAYER4',
  projection_1_noncanon => 'LAYER7',
  projection_2_noncanon => 'LAYER7',
  projection_3_noncanon => 'LAYER7',
  projection_4_noncanon => 'LAYER7',
  projection_1_pseudo   => 'LAYER7',
  projection_2_pseudo   => 'LAYER7',
  projection_3_pseudo   => 'LAYER7',
  projection_4_pseudo   => 'LAYER7',
);

my %projection_occurrences;
for my $layer (@$layer_config) {
  for my $biotype (@{ $layer->{BIOTYPES} }) {
    $projection_occurrences{$biotype}++ if exists $expected_projection_layer{$biotype};
    is($layer->{ID}, $expected_projection_layer{$biotype}, "$biotype is in the expected layer")
      if exists $expected_projection_layer{$biotype};
  }
}
is($projection_occurrences{$_} || 0, 1, "$_ occurs exactly once")
  for sort keys %expected_projection_layer;

my %expected_evidence_layer = (
  cdna2genome => 'LAYER2', edited => 'LAYER2', gw_gtag => 'LAYER2', gw_exo => 'LAYER2',
  rnaseq_tissue_4 => 'LAYER5', rnaseq_tissue_5 => 'LAYER5',
  cdna_4 => 'LAYER5', cdna_5 => 'LAYER5',
  rnaseq_tissue_6 => 'LAYER6', cdna_6 => 'LAYER6',
  cdna => 'LAYER10', rnaseq_tissue => 'LAYER10',
);
$expected_evidence_layer{"rnaseq_tissue_$_"} = 'LAYER3' for 1 .. 3;
$expected_evidence_layer{"cdna_$_"} = 'LAYER3' for 1 .. 3;
for my $layer (@$layer_config) {
  my %biotypes = map { $_ => 1 } @{ $layer->{BIOTYPES} };
  is($layer->{ID}, $expected_evidence_layer{$_}, "$_ is in the expected layer")
    for grep { $biotypes{$_} } keys %expected_evidence_layer;
}

my $genebuilder_config = Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic->new
  ->_master_config('zebrafish_basic');
my %genebuilder_biotypes = map { $_ => 1 } @$genebuilder_config;
ok($genebuilder_biotypes{$_}, "genebuilder retains $_")
  for keys %expected_projection_layer;

my $transcript_selection = "$Bin/../Bio/EnsEMBL/Analysis/Hive/Config/TranscriptSelection.pm";
open my $fh, '<', $transcript_selection or die "Cannot read $transcript_selection: $!";
local $/;
my $transcript_selection_source = <$fh>;
close $fh;
my ($utr_donor_block) = $transcript_selection_source =~ /(utr_donor_dbs\s*=>\s*\[.*?\])/s;
like(
  $utr_donor_block || '',
  qr/utr_donor_dbs\s*=>\s*\[\s*\$self->o\('cdna_db'\).*?\$self->o\('rnaseq_for_layer_db'\).*?\$self->o\('long_read_final_db'\).*?\]/s,
  'UTR donors are cDNA, RNA layer, and long-read databases',
);
unlike($utr_donor_block || '', qr/selected_projection_db/, 'selected projections are not UTR donors');

done_testing();
