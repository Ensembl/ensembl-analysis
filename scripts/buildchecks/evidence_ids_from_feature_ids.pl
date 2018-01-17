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

=head1 NAME

evidence_ids_from_feature_ids.pl

=head1 SYNOPSIS

perl evidence_ids_from_feature_ids.pl  -host host -user ensro -port 3306 
-dbname gwdb -feature_id 4 

This would return the protein based evidence ids for gene with dbID 4

=head1 DESCRIPTION

This script will fetch either the ids of protein or dna evidence supporting
a the gene  with the given dbID

=head1 OPTIONS

  -dbhost database host
  -dbuser database user
  -dbpass database password
  -dbport database port
  -dbname database name
  -feature_id dbID or stable id of feature to get ids for
  -evidence_type protein_align_feature or dna_align_feature, protein by
  default
  -id_type gene or transcript gene by default
  -gene_type type from gene table
  -primary_evidence, only consider primary evidence
  -id_list_file file containing a list of dbIDs
  -info print out descriptions of the ids
  -verbose print out information about the run
  -pfetch indicate you have pfetch to fetch descriptions. This is on
  by default but can be switched off with -nopfetch
  -help print perl docs

=head1 EXAMPLES

perl evidence_ids_from_feature_ids.pl -host host -user ensro -port 3306 -dbname gw_db -feature_id 1

This will find the protein evidence associated with the gene with the dbID
1

perl evidence_ids_from_feature_ids.pl -host host -user ensro -port 3306 -dbname gw_db -feature_id 4 -info -feature_type transcript -evidence_type
dna_align_feature

This will find the dna evidence associated with the transcript with the 
dbID 1 and print out the descriptions from pfetch

=head1 NOTES

By default this script gets the evidence associated with all exons which
are part of the gene or transcript id passed in. If the -primary_evidence
tag is used only the transcript_supporting_evidence table is queryed so
only primary evidence is returned

Note this uses a module which is found in this directory so you must
either run this script in your directory or put this directory in your
PERL5LIB

Also by default these scripts consider protein evidence and gene ids

=head1 CONTACT

the ensembl-dev mailing list <http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use warnings ;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use SupportingEvidenceInfo;


my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my $protein_id;
my $table_name = 'protein_align_feature';
my $gene_type;
my $protein_file;
my $info;
my $id_type = 'gene';
my $verbose;
my $pfetch = 1;
my $help;
my $only_primary;
GetOptions( 
            'dbhost|host|h=s'      => \$dbhost,
            'dbname|db|D=s'      => \$dbname,
            'dbuser|user|u=s'      => \$dbuser,
            'dbpass|pass|p=s'      => \$dbpass,
            'dbport|port|P=s'      => \$dbport,
            'feature_id=s' => \$protein_id,
            'evidence_type=s' => \$table_name,
            'id_type=s' => \$id_type,
            'gene_type=s' => \$gene_type,
            'primary_evidence!' => \$only_primary, 
            'id_list_file=s' => \$protein_file,
            'info!' => \$info,
            'verbose!' => \$verbose,
            'pfetch!' => \$pfetch,
            'help!' => \$help,
           ) or perldocs("Failed to get opts");

perldocs() if($help);

if(!$dbhost || !$dbname || !$dbuser){
  throw("Must pass in host, user and dbname with -host $dbhost -user $dbuser".
        " -dbname $dbname ");

}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $dbhost,
   -user   => $dbuser,
   -dbname => $dbname,
   -pass   => $dbpass,
   -port   => $dbport,
  );



my $evidence_info = SupportingEvidenceInfo->new();
$evidence_info->db($db);
$evidence_info->verbose($verbose);
$evidence_info->evidence_type($table_name);
$evidence_info->id_type($id_type);
$evidence_info->info($info);
$evidence_info->have_pfetch($pfetch);
$evidence_info->primary_evidence($only_primary);
my $ids = [];
if($protein_file){
  $ids = $evidence_info->read_id_file($protein_file);
}else{
  push(@$ids, $protein_id);
}

foreach my $id(@$ids){
  if(!($id =~ /^\d+/)){
    my $temp = $id;
    $id = $evidence_info->dbID_from_stable_id($temp);
  }
  $evidence_info->evidence_ids_from_feature_id($id);
}


sub perldocs{
  my ($msg) = @_;
  print $msg."\n" if($msg);
  exec('perldoc', $0);
  exit(0);
}
