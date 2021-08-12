#!/usr/bin/env perl

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

## Script to scan a core DB and replace the Ensembl stable IDs with custom
## braker style stable ids where the gene in question has been derived from
## running braker 2 pipeline (determined by source field).
## So, for example, gene ENSVATG00005000003 would become BRAKER_VATG_5000003
## Requires a core DB and connection details, plus admin user to update the
## stable ids.
## Example:
## perl stable_id_convert.pl -u ensro -h $GBS6 -dbname jma_vanessa_atalanta_gca905147705v1_core_test -p $GBP6

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $help = 0;
my $port = 3306;
my ( $host, $dbname, $user, $pass );

&GetOptions(
  'dbhost|host|h:s'   => \$host,
  'dbuser|user|u:s'   => \$user,
  'dbpass|pass|p:s'   => \$pass,
  'dbname|db|D:s'     => \$dbname,
  'dbport|port|P:n'   => \$port,
);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -dbname => $dbname,
  -pass   => $pass,
  ) or throw( "can not connect to $dbname@$host\n");

my $slice_adaptor = $db->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all( 'toplevel', undef, 1 );
my $ga = $db->get_GeneAdaptor();
my $braker = 0;

while ( my $slice = shift @{$slices} ) {
  my $genes = $slice->get_all_Genes();

  while ( my $gene = shift @{$genes} ) {

    ## Check if gene model is derived from braker, if so update stable id
    ## for gene and all associated transcripts/exons/tranlsations.
    if ( $gene->source =~ m/^braker_/ ) {
      $braker = 1;
#      print "\n\n", $gene->dbID, "\t", $gene->stable_id, "\t", $gene->source;
      my $stable_g = update_stable($gene->stable_id);
#      print "\t$stable_g";
      update_core($gene->stable_id, $stable_g, "gene");

      my @transcripts = @{ $gene->get_all_Transcripts };

      while ( my $transcript = shift @transcripts ) {
#        print "\n", $transcript->dbID, "\t", $transcript->stable_id, "\t";
        my $stable_t = update_stable( $transcript->stable_id );
#        print "\t$stable_t";
        update_core($transcript->stable_id, $stable_t, "transcript");

        my @exons = @{ $transcript->get_all_Exons() };

        while ( my $exon = shift @exons ) {
#          print "\n", $exon->dbID, "\t", $exon->stable_id, "\t";
          my $stable_e = update_stable( $exon->stable_id );
#          print "\t$stable_e";
          update_core($exon->stable_id, $stable_e, "exon");
        }

        if ( defined( $transcript->translation() ) ) {
#          print "\n", $transcript->translation()->stable_id();
          my $stable_p = update_stable( $transcript->translation()->stable_id );
#          print "\t$stable_p";
          update_core($transcript->translation()->stable_id, $stable_p, "translation");
        }
      }
    }
  }
}

if ( $braker == 0 ) {
  print "No braker derived genes were found in $dbname, no changes were made.\n";
}

## accepts ensembl style stable ID and returns custom braker stable id.
sub update_stable {

  my $stable = shift;
  if ( $stable =~ m/^ENS[A-Z]+[0-9]+/ ) {
    $stable =~ s/ENS/BRAKER_/;
    $stable =~ s/0+(?=[0-9])/_/;
  }
  ## Note: this will probably trigger on exons that are part of multiple transcripts
  ## and can safely be ignored in that case
  elsif ( $stable =~ m/^BRAKER_/ ) {
    print "Looks like $stable is already a braker id.\n";
  }
  else {
    throw( "No Ensembl like stable IDs found, perhaps run_stable_ids did not complete.\n" );
  }
  return $stable;
}

sub update_core {
  my ($old_stable, $new_stable, $table) = (@_);
  my $dbc = $db->dbc;

  my $sql = "UPDATE " . $table
            . " SET stable_id = '" . $new_stable . "'"
            . " WHERE stable_id = '" . $old_stable . "'";
#  print $sql, "\n";

  my $update_sth = $dbc->prepare($sql);
  $update_sth->execute;
}
