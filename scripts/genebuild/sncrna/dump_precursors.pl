#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my ($dbname, $dbhost, $dbport, $dbuser, $dna_dbname, $dna_dbhost, $dna_dbport, $dna_dbuser, $dna_dbpass, $working_dir);


GetOptions( 'dbhost|host|h:s'        => \$dbhost,
            'dbport|port|p:n'        => \$dbport,
            'dbname|d:s'        => \$dbname,
            'dbuser|user|u:s'        => \$dbuser,
            'dnadbname|D:s' => \$dna_dbname,
            'dnadbport|P:n' => \$dna_dbport,
            'dnadbhost|H:s' => \$dna_dbhost,
            'dnadbuser|U:s' => \$dna_dbuser,
            'password|e:s' => \$dna_dbpass,
            'output_dir|o:s' => \$working_dir);


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-DBNAME => $dbname,
  	-HOST => $dbhost,
  	-PORT => $dbport,
  	-USER => $dbuser,
	-DRIVER => 'mysql',
);

my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -DBNAME => $dna_dbname,
  -HOST => $dna_dbhost,
  -PORT => $dna_dbport,
  -USER => $dna_dbuser,
  -DRIVER => 'mysql',
);

if ($dnadb){
  $db->dnadb($dnadb);
}

# dump putative stem-loops
my $gene_adaptor = $db->get_GeneAdaptor();
my @genes = @{ $gene_adaptor->fetch_all_by_biotype('miRNA')};

my $fn = $working_dir . "/identified_mirnas.fa";

open(FH, '>', $fn) or die "Could not write to $fn";

foreach my $gene (@genes){
    my $strand = $gene->strand() > 0 ? "+" : "-";

    print FH ">" .
          $gene->dbID(), ":",
          $gene->seq_region_name(), ":",
          $gene->seq_region_start(), "-",
          $gene->seq_region_end(), ":",
          $strand, "\n",
          $gene->seq(), "\n";

}

close(FH);
