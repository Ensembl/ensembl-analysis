#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

=head1

Simply cluster the merged gene set with the ensembl only gene set
and find ensembl-only clusters

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Getopt::Long;

# ensembl genes
my $ensemblhost;
my $ensembluser;
my $ensembldbname;
my $ensemblport   = 3306;

# merged genes
my $mergedhost;
my $mergeduser;
my $mergeddbname;
my $mergedport   = 3306;

# output
my $outfile = 'stdout';
my $verbose;

# genome bits
my @default_ensembl_genetypes = ('IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'Mt_rRNA', 'Mt_tRNA', 'processed_pseudogene', 'protein_coding', 'pseudogene');
my @genetypes = ();
my @default_merged_genetypes = ('3prime_overlapping_ncrna', 'antisense', 'IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene','IG_V_gene', 'IG_V_pseudogene', 'lincRNA', 'miRNA', 'misc_RNA', 'Mt_rRNA', 'Mt_tRNA', 'non_coding', 'polymorphic_pseudogene', 'processed_pseudogene', 'processed_transcript', 'protein_coding', 'pseudogene', 'rRNA', 'sense_intronic', 'sense_overlapping', 'snoRNA', 'snRNA', 'TR_V_gene', 'TR_V_pseudogene');
my @merged_genetypes = ();


# set names
my $ensembl_setname = 'ensembl';
my $merged_setname = 'merged';

# assembly
my $coord_system_version;
my $autotype = 0;

#  ~~~~~~
#  No changes required below this point
#  ~~~~~~

# where do we get the ensembl structures from:

&GetOptions(
  'ensemblhost:s'          => \$ensemblhost,
  'ensembluser:s'          => \$ensembluser,
  'ensembldbname:s'        => \$ensembldbname,
  'ensemblport:n'          => \$ensemblport,
  'mergedhost:s'           => \$mergedhost,
  'mergeduser:s'           => \$mergeduser,
  'mergeddbname:s'         => \$mergeddbname,
  'ensemblhost:s'             => \$ensemblhost,
  'ensembluser:s'             => \$ensembluser,
  'ensembldbname:s'           => \$ensembldbname,
  'ensemblport:n'             => \$ensemblport,
  'coord_system_version:s' => \$coord_system_version,
  'mergedport:n'           => \$mergedport,
  'merged_genetypes:s'     => \@merged_genetypes,
  'genetypes:s'            => \@genetypes,
  'outfile:s'              => \$outfile,
  'verbose'                => \$verbose,
  'mergedsource=s'         => \$merged_setname,
  'ensemblsource=s'           => \$ensembl_setname,
  'autotype!'           => \$autotype,
);

print STDERR "Ensembl database: name $ensembldbname host $ensemblhost port $ensemblport, \nMerged database: name $mergeddbname host $mergedhost port $mergedport\n";

if (scalar(@merged_genetypes)) {
  @merged_genetypes = split(/,/,join(',',@merged_genetypes));
} else {
  @merged_genetypes = @default_merged_genetypes;
}
if (scalar(@genetypes)) {
  @genetypes = split(/,/,join(',',@genetypes));
} else {
  @genetypes = @default_ensembl_genetypes;
}

# connect to dbs and get adaptors
my $ensembldb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $ensemblhost,
  -user   => $ensembluser,
  -port   => $ensemblport,
  -dbname => $ensembldbname
);
my $ensemblsa = $ensembldb->get_SliceAdaptor();

my $mergeddb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $mergedhost,
  -user   => $mergeduser,
  -port   => $mergedport,
  -dbname => $mergeddbname
);
my $mergedsa   = $mergeddb->get_SliceAdaptor();

# open outfile
my $fh;
if ($outfile && $outfile ne "stdout") {
  open FH,">$outfile" or die "couldn't open file ".$outfile." $!";
  $fh = \*FH;
} else {
  $fh = \*STDOUT;
}


if ($autotype) {
    my $sqlquery = 'SELECT DISTINCT(biotype) FROM gene';
    my $sth = $ensembldb->dbc->prepare($sqlquery);
    $sth->execute();
    my $rows = $sth->fetchall_arrayref;
    @genetypes = ();
    foreach my $row (@$rows) {
        push(@genetypes, $row->[0]);
    }
    $sth = $mergeddb->dbc->prepare($sqlquery);
    $sth->execute();
    @merged_genetypes = ();
    $rows = $sth->fetchall_arrayref;
    foreach my $row (@$rows) {
        push(@merged_genetypes, $row->[0]);
    }
}

# # #
# OK, now we begin to do stuff
# # #

my $toplevels = $ensemblsa->fetch_all('toplevel',$coord_system_version,0,1);

foreach my $ensemblslice (sort @$toplevels) {
  print "ENSEMBL ".$ensemblslice->name."\n";

  my $count;
  my $mergedslice = $mergedsa->fetch_by_name($ensemblslice->name);
  print "MERGED ".$mergedslice->name."\n\n";

  # # #
  # now we only fetch the merged gene biotypes that we are interested in
  # # #
  my ($mergedgenes, $all_merged_biotypes) = do_gene_fetch($mergedslice, $merged_setname, \@merged_genetypes, undef);
  print STDERR "\nGot " . scalar(@$mergedgenes) . " $merged_setname genes\n" if $verbose;

  # # #
  # now we only fetch the ensembl gene biotypes that we are interested in
  # # #
  my ($ensemblgenes,$all_ensembl_biotypes) = do_gene_fetch($ensemblslice, $ensembl_setname, \@genetypes, undef);
  print STDERR "\nGot " . scalar(@$ensemblgenes) . " $ensembl_setname genes\n" if $verbose;

  # # #
  # Now cluster the genes and look for clusters with > 1 gene
  # # #
  # make typeshash
  my @genes = (@$mergedgenes,@$ensemblgenes);
  my %types_hash;
  $types_hash{$merged_setname} = [@$all_merged_biotypes];
  $types_hash{$ensembl_setname} = [@$all_ensembl_biotypes];


  # cluster genes
  print STDERR "\nClustering genes on all exons...\n" if $verbose;
  # cluster on all exons, not just coding exons
  my ($clusters, $unclustered) = cluster_Genes(\@genes, \%types_hash, undef);

  # # #
  # loop thru clusters
  # # #
  print STDERR "Have clusters ".scalar(@$clusters)." and ".scalar(@$unclustered)." unclustered\n";

  UNCLUSTERED: foreach my $uncl (@$unclustered) {
    my $genomic_location = get_genomic_location($ensemblslice->name, $uncl);

    my @ensemblgenes = @{$uncl->get_Genes_by_Set($ensembl_setname)};

    foreach my $g (@ensemblgenes) {
      # we only want to report the unclustered ensembl genes
      print "Missing gene ".$g->stable_id." biotype ".$g->biotype." slice $genomic_location\n";
    }
  } # clusyer
} # chr



sub do_gene_fetch {
  my ($slice, $setname, $biotypes) = @_;
  my @genes;
  my %biotypes;

  foreach my $genetype (@$biotypes) {
    my @tmp = @{$slice->get_all_Genes(undef, undef, 1, undef, $genetype)};
    print STDERR "Got " . scalar(@tmp) . " $setname $genetype genes\n" if $verbose;
    foreach my $g (@tmp) {
      $g->biotype($setname."_".$g->biotype);
      push @genes,$g;
      $biotypes{$g->biotype} = 1;
    }
  }
  my @bio = keys %biotypes;
  print "Got the following $setname gene biotypes : @bio" if $verbose;
  return \@genes, \@bio;
}

sub get_genomic_location {
  my ($name, $cluster) = @_;

  my @genomic_fields = split(":",$name);
  my $genomic_location = $genomic_fields[0].":".$genomic_fields[1].":".
                         $genomic_fields[2].":".$cluster->start.":".
                         $cluster->end.":".$genomic_fields[5];

  return $genomic_location;
}
