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

=head1 NAME
  get_transcripts.pl
=head1 DESCRIPTION
Retrieves the DNA for all transcripts from a given database and dumps out to fasta format file.
=head1 OPTIONS
=head1
=head2 Database connection options
  -dbhost   Host name for database.
  -dbport   What port to connect to.
  -dbname   What name to connect to.
  -dbuser   What username to connect as.
  -dbpass   What password to use.

  -dna_dbhost   Host name for database.
  -dna_dbport   What port to connect to.
  -dna_dbname   What name to connect to.
  -dna_dbuser   What username to connect as.
  -dna_dbpass   What password to use.

=head2 Other options
  --coords coordinate system name to use, default toplevel
=head1 EXAMPLES
  ./get_transcripts.pl --dbhost genebuild5 \
    --dbuser ensro \
    --dbname <database_name> \
    --dna_dbhost genebuild6 \
    --dna_dbname <dna_database_name> \
    > output.fasta
=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;

## Some defaults

my $coord_system = 'toplevel';

my $dbname = '';
my $user   = 'ensro';
my $host = $ENV{GBS5};
my $port = $ENV{GBP5};
my $pass = '';

my $dna_dbname = '';
my $dna_user   = 'ensro';
my $dna_host   = $ENV{GBS6};
my $dna_port   = $ENV{GBP6};
my $dna_pass  = '';

my $options = GetOptions ("user|dbuser|u=s"         => \$user,
                          "host|dbhost|h=s"         => \$host,
                          "port|dbport|P=i"         => \$port,
                          "dbname|db|D=s"           => \$dbname,
                          "dbpass|pass|p=s"         => \$pass,

                          "dna_user|dna_dbuser|u=s" => \$dna_user,
                          "dna_host|dna_dbhost|h=s" => \$dna_host,
                          "dna_port|dna_dbport|P=i" => \$dna_port,
                          "dna_dbname|dna_db|D=s"   => \$dna_dbname,
                          "dna_dbpass|dna_pass|p=s" => \$dna_pass,

                          "coords:s"                => \$coord_system,);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dna_port,
  -user    => $dna_user,
  -host    => $dna_host,
  -dbname  => $dna_dbname,
  -pass    => $dna_pass);

$db->dnadb($dna_db);

my $slice_adaptor = $db->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all('toplevel', undef, 1 );
my $gene_adaptor = $db->get_GeneAdaptor();

while ( my $slice = shift @{$slices} ) {
  my $genes = $slice->get_all_Genes();
  while ( my $gene = shift @{$genes} ) {
    my @transcripts = @{$gene->get_all_Transcripts};
    if (!@transcripts) {
      throw("No transcripts found in $dna_dbname\n");
    }
    while ( my $transcript = shift @transcripts ) {
      my $name =  ">". $transcript->seq_region_name. ":" .$transcript->seq_region_start. ":" .$transcript->seq_region_end. ":" .$transcript->biotype. ":" .$transcript->dbID;
      print $name, "\n", $transcript->seq->seq, "\n";
    }
  }
}
