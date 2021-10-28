#!/usr/bin/env perl
# Copyright [2021] EMBL-European Bioinformatics Institute
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

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $dbname;
my $user;
my $host;
my $port;
my $pass;

my $options = GetOptions ("user|dbuser|u=s"   => \$user,
                          "host|dbhost|h=s"   => \$host,
                          "port|dbport|P=i"   => \$port,
                          "dbname|db|D=s"     => \$dbname,
                          "dbpass|pass|p=s"   => \$pass);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

my $gene_adaptor = $db->get_GeneAdaptor();
my $genes = $gene_adaptor->fetch_all();

say "Looping through genes...";
foreach my $gene (@$genes) {
  my $biotype = $gene->biotype;
  unless($gene->biotype eq 'lncRNA') {
    next;
  }

  my $transcripts = $gene->get_all_Transcripts;
  my $remove = remove_small_or_overlapping($gene);
  if($remove) {
    $gene_adaptor->remove($gene);
  }

}

sub remove_small_or_overlapping {
  my ($gene) = @_;

  my $overlap_okay = 0;
  my $overlapping_genes = $gene->get_overlapping_Genes();
  foreach my $overlapping_gene (@$overlapping_genes) {
    if($overlapping_gene->biotype eq 'protein_coding') {
      say "Removing gene due to overlap";
      return(1);
    }
  }

  # check min length
  my $length_okay = 0;
  my $transcripts = $gene->get_all_Transcripts;
  foreach my $transcript (@$transcripts) {
    if($transcript->length >= 200) {
      $length_okay = 1;
    }
  }

  unless($length_okay) {
    say "Removing gene due to length";
    return(1);
  }

  return(0);
}


#sub remove_gene {
#  my ($transcript) = @_;

#  my $min_seq_size = 1000;
#  my $exons = $transcript->get_all_Exons;

#  if($transcript->length < 200) {
#    say "Found small gene: ".$transcript->length;
#    return(1);
#  }

#  if(scalar(@$exons) >= 3 ) {
#    return(0);
#  } elsif($transcript->length >= $min_seq_size) {
#    my $overlapping_genes = $transcript->get_overlapping_Genes(1);
#    if($overlapping_genes) {
#      foreach my $gene (@$overlapping_genes) {
#        if($gene->biotype eq 'protein_coding') {
#          say "Found overlapping protein_coding gene on same strand";
#          return(2);
#	}
#      }
#    }
#    return(0);
#  } else {
#    return(1);
#  }
#}
