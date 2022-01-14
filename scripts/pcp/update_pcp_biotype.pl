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

# Updates biotype in pcp db if BOTH protein coding prediction tools agree that
# the gene is protein coding

=head1 NAME
  update_pcp_biotype.pl
=head1 DESCRIPTION
Reads the output of both RNASamba and CPC2 and where they both agree that a transcript is protein coding the biotype
is updated to pcp_protein_coding in the selected database.
Throws if results files have non matching gene models/ number of gene models
=head1 OPTIONS
=head1
=head2 Database connection options
  -dbhost   Host name for database.
  -dbport   What port to connect to.
  -dbname   What name to connect to.
  -dbuser   What username to connect as.
  -dbpass   What password to use.
=head2 Other options
  --coords coordinate system name to use, default toplevel
  --rnaSamba RNASamba output file to be used
  --cpc2 CPC2 output file to be used
=head1 EXAMPLES
  ./update_pcp_biotype.pl --dbhost genebuild5 \
    --dbname <database_name> \
    --dbhost GBS5
    --dport GBP5
    --rnaSamba <RNASamba output file> \
    --cpc2 <CPC2 output file> \
=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

## Some defaults
my $coord_system = 'toplevel';
my $dbname = '';
my $user   = 'ensadmin';
my $host = '';
my $port = '';
my $pass = '';
my ($cpc2_file, $rnasamba_file ) = ( '', '' );

my $options = GetOptions ("user|dbuser|u=s"     => \$user,
                          "host|dbhost|h=s"     => \$host,
                          "port|dbport|P=i"     => \$port,
                          "dbname|db|D=s"       => \$dbname,
                          "dbpass|pass|p=s"     => \$pass,
                          "cpc2=s"              => \$cpc2_file,
                          "rnaSamba|rnas=s"     => \$rnasamba_file,
                          "coords:s"            => \$coord_system,);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

my %cpc_results = parse_results($cpc2_file, 8);
my %rnasamba_results = parse_results($rnasamba_file, 2);
my %selected_genes;

if (scalar(keys %cpc_results) != scalar(keys %rnasamba_results)) {
  throw("Results files have different number of gene models\n");
}

foreach my $key (keys %cpc_results) {
  throw("Results files contain different gene models\n") unless (exists $rnasamba_results{$key});
  ## Currently using AND to account for relative sensitivity  / specifity differences
  ## between the algorithms.
  if ($rnasamba_results{$key} eq 'coding' and $cpc_results{$key} eq 'coding') {
      $selected_genes{$key} = 1;
  }
}

my $slice_adaptor = $db->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all('toplevel', undef, 1 );
my $gene_adaptor = $db->get_GeneAdaptor();

while (my $slice = shift @{$slices}) {
  my $genes = $slice->get_all_Genes();
  while (my $gene = shift @{$genes}) {
    my $check = $gene->dbID;
    if ($selected_genes{$check}) {
      my $current = $gene->biotype;
      $gene->biotype('pcp_protein_coding');
      $gene_adaptor->update($gene);
    }
  }
}

=head2
Parses results from the classification program output, akes the filename and
column where protein coding classifier stored. Throws if files not found
=cut

sub parse_results {
  my ($in_file, $column_num) = @_;
  my %results;
  if (!$in_file) {
    throw("Input file is missing!\n");
  }
  open(my $FILE,'<',$in_file) or die("Could not open file ".$in_file);
  my $header = <$FILE>;

  while(<$FILE>) {
    chomp( my @row = split'\t', $_ );
    my @temp = split':', $row[0];
    $results{$temp[-1]} = $row[$column_num];
  }
  
  close($FILE) or die("Could not close file ".$in_file);
  return %results;
}
