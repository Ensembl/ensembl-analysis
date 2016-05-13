
use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan tests => 33; }

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::ReadBaseDataFromGff;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

ok(1);

my $gff = ReadBaseDataFromGff->new();
$gff->file("gff_example.gff");
my @example_genes = @{ $gff->get_Genes() };
my $gene1         = $example_genes[0];
my $gene2         = $example_genes[1];
my $gene3         = $example_genes[2];
my $gene4         = $example_genes[3];
my $gene5         = $example_genes[4];

$gene1->biotype("set1");
$gene2->biotype("set2");
$gene3->biotype("set2");
$gene4->biotype("set2");
$gene5->biotype("set2");

my @gene1 = ($gene1);
my @gene2 = ($gene2);
my @gene3 = ($gene3);
my @gene4 = ($gene4);
my @gene5 = ($gene5);

my $types_hash = make_types_hash( \@gene1, \@gene2, "set1", "set2" );
my @genes1 = ( $gene1, $gene2 );
print "Testing cluster genes \n";
my ( $clustered1, $unclustered1 ) = cluster_Genes( \@genes1, $types_hash );
if ( scalar( @{$clustered1} ) > 0 ) {
  print "Looking in the clustered\n";
  my @new_genes1 = @{ @{$clustered1}[0]->get_Genes() };
  all_in_one( \@new_genes1, \@genes1 );
}
else {
  throw("Didn't cluster\n");
}
print "\n";

print "Testing cluster genes without strand\n";
my ( $clustered2, $unclustered2 ) = cluster_Genes_without_strand( \@genes1, $types_hash );
if ( scalar( @{$clustered2} ) > 0 ) {
  print "Looking in the clustered\n";
  my @new_genes2 = @{ @{$clustered2}[0]->get_Genes() };
  all_in_one( \@new_genes2, \@genes1 );
}
else {
  throw("Didn't cluster!");
}
print "\n";

print "Testing cluster genes on coding\n";
my ( $clustered3, $unclustered3 ) = cluster_Genes_by_coding_exon_overlap( \@genes1, $types_hash );
if ( scalar( @{$unclustered3} ) > 0 ) {
  print "Looking in the unclustered\n";
  my @new_genes3;
  foreach my $noncluster ( @{$unclustered3} ) {
    push @new_genes3, @{ $noncluster->get_Genes() };
  }
  all_in_one( \@new_genes3, \@genes1 );
}
else {
  throw("Shouldn't cluster !");
}
print "\n";

print "Testing cluster genes on coding after adding translation\n";
$gff->add_Translation($gene1);
$gff->add_Translation($gene2);

my ( $clustered4, $unclustered4 ) = cluster_Genes_by_coding_exon_overlap( \@genes1, $types_hash );
if ( scalar( @{$clustered4} ) > 0 ) {
  print "Looking in the clustered\n";
  my @new_genes4 = @{ @{$clustered4}[0]->get_Genes() };
  all_in_one( \@new_genes4, \@genes1 );
}
else {
  throw("Didn't cluster !");
}
print "\n";

print "Testing cluster genes with the make_types_hash_with_genes method directly\n";
my ( $types_hash2, $genes5 ) = @{ make_types_hash_with_genes( \@gene1, \@gene2, "set1", "set2" ) };
my ( $clustered5, $unclustered5 ) = cluster_Genes( $genes5, $types_hash2 );
if ( scalar( @{$clustered5} ) > 0 ) {
  print "Looking in the clustered\n";
  my @new_genes5 = @{ @{$clustered5}[0]->get_Genes() };
  all_in_one( \@new_genes5, \@genes1 );
}
else {
  throw("Didn't cluster !");
}
print "\n";

print "Testing simple_cluster_Genes\n";
my ( $clustered6, $unclustered6 ) = @{ simple_cluster_Genes( \@gene1, "set1", \@gene2, "set2" ) };
if ( scalar( @{$clustered6} ) > 0 ) {
  print "Looking in the clustered\n";
  my @new_genes6 = @{ @{$clustered6}[0]->get_Genes() };
  all_in_one( \@new_genes6, \@genes1 );
}
else {
  throw("Didn't cluster !");
}
print "\n";

print "Testing get_twoway_clusters\n";
my $cluster7   = get_twoway_clusters($clustered6);
my @new_genes7 = @{ @{$cluster7}[0]->get_Genes() };
all_in_one( \@new_genes7, \@genes1 );
print "\n";

print "Testing get_twoway_clustering_genes_of_set\n";
print "Twoway clusters with set 1 from twoway clusters\n";
my @new_genes8 = @{ get_twoway_clustering_genes_of_set( $clustered6, "set1" ) };
ok( scalar(@new_genes8), 1 );
ok( $new_genes8[0],      $gene1 );
print "Twoway clusters with set 2 from all clusters\n";
my @new_genes9 = @{ get_twoway_clustering_genes_of_set( $clustered6, "set2" ) };
ok( scalar(@new_genes9), 1 );
ok( $new_genes9[0],      $gene2 );
print "\n";

print "Testing cluster genes on new set\n";
my @genes10 = ( $gene1, $gene3 );
my ( $clustered10, $unclustered10 ) = cluster_Genes( \@genes10, $types_hash );
if ( scalar( @{$unclustered10} ) > 0 ) {
  print "Looking in the unclustered\n";
  my @new_genes10;
  foreach my $noncluster ( @{$unclustered10} ) {
    push @new_genes10, @{ $noncluster->get_Genes() };
  }
  all_in_one( \@new_genes10, \@genes10 );
}
else {
  throw("Shouldn't cluster !");
}
print "\n";

print "Testing cluster genes without strand\n";
my ( $clustered11, $unclustered11 ) = cluster_Genes_without_strand( \@genes10, $types_hash );
if ( scalar( @{$clustered11} ) > 0 ) {
  print "Looking in the clustered\n";
  my @new_genes11 = @{ @{$clustered11}[0]->get_Genes() };
  all_in_one( \@new_genes11, \@genes10 );
}
else {
  throw("Didn't cluster !");
}
print "\n";

print "Testing cluster genes on new set\n";
my @genes12 = ( $gene1, $gene4 );
my ( $clustered12, $unclustered12 ) = cluster_Genes_without_strand( \@genes12, $types_hash );
if ( scalar( @{$unclustered12} ) > 0 ) {
  print "Looking in the unclustered\n";
  my @new_genes12;
  foreach my $noncluster ( @{$unclustered12} ) {
    push @new_genes12, @{ $noncluster->get_Genes() };
  }
  all_in_one( \@new_genes12, \@genes12 );
}
else {
  throw("Shouldn't cluster !");
}
print "\n";

print "Testing get_single_cluster\n";
my @genes13 = ( $gene2, $gene3 );
my ( $clustered13, $unclustered13 ) = cluster_Genes_without_strand( \@genes13, $types_hash );
my $cluster13   = get_single_clusters($clustered13);
my @new_genes13 = @{ @{$cluster13}[0]->get_Genes() };
all_in_one( \@new_genes13, \@genes13 );
print "\n";

print "Testing get_oneway_clustering_genes_of_set\n";
my @genes_set14 = get_oneway_clustering_genes_of_set( $clustered13, "set2" );
all_in_one( \@genes13, @genes_set14 );
print "\n";

print "Testing cluster genes on coding with non coding exon overlap\n";
$gff->add_Translation($gene5);
my @genes15 = ( $gene1, $gene5 );

my ( $clustered15, $unclustered15 ) = cluster_Genes_by_coding_exon_overlap( \@genes15, $types_hash );
if ( scalar( @{$unclustered15} ) > 0 ) {
  print "Looking in the unclustered\n";
  my @new_genes15;
  foreach my $noncluster ( @{$unclustered15} ) {
    push @new_genes15, @{ $noncluster->get_Genes() };
  }
  all_in_one( \@new_genes15, \@genes15 );
}
else {
  throw("Shouldn't cluster !");
}
print "\n";

print "Checking if it clusters on genomic overlap\n";
my ( $clustered16, $unclustered16 ) = cluster_Genes( \@genes15, $types_hash );
if ( scalar( @{$clustered16} ) > 0 ) {
  print "Looking in the clustered\n";
  my @new_genes16 = @{ @{$clustered16}[0]->get_Genes() };
  all_in_one( \@new_genes16, \@genes15 );
}
else {
  throw("Didn't cluster !");
}
print "\n";

sub all_in_one {
  my ( $cluster, $genes ) = @_;
  my @new_genes = sort { $a <=> $b } @$cluster;
  my $count = 0;
  foreach my $gene ( sort { $a <=> $b } @$genes ) {
    ok( $new_genes[$count], $gene );
    $count++;
  }

}

