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

evidence_ids_from_evidence_ids.pl

=head1 SYNOPSIS

perl feature_ids_from_evidence.pl  -host host -user ensro -port 3306 
-dbname gwdb -evidence_id Q60D5

This would return the gene ids supported by the protein id Q60D5

=head1 DESCRIPTION

This script will fetch the gene or transcripts supported by the supplied
evidence identified

=head1 OPTIONS

  -dbhost database host
  -dbuser database user
  -dbpass database password
  -dbport database port
  -dbname database name
  -evidence_id dbID of feature to get ids for
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

perl feature_ids_from_evidence.pl -host yourhost -user ensro -port 3306 -dbname gw_db -evidence_id Q9CW79-1

This will find the gene ids associatied with Q9CW79-1

perl feature_ids_from_evidence.pl -host yourhost -user ensro -port 3306 -dbname gw_db -evidence_id NM_082790.1 -id_type transcript -evidence_type
dna_align_feature

This will find the transcript ids associated with dna evidence NM_082790.1

=head1 NOTES

By default this script joins to the transcript table from the align_feature
tables via the exon and supporting_feature tables but if the primary
evidence tag is used it will only look to the transcript_supporting_feature
table so only find features which have the evidence supplied as primary
evidence

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
my $help;
my $only_primary;
GetOptions( 
            'dbhost|host|h=s'      => \$dbhost,
            'dbname|db|D=s'      => \$dbname,
            'dbuser|user|u=s'      => \$dbuser,
            'dbpass|pass|p=s'      => \$dbpass,
            'dbport|port|P=s'      => \$dbport,
            'evidence_id=s' => \$protein_id,
            'evidence_type=s' => \$table_name,
            'id_type=s' => \$id_type,
            'primary_evidence!' => \$only_primary, 
            'gene_type=s' => \$gene_type,
            'id_list_file=s' => \$protein_file,
            'info!' => \$info,
            'verbose!' => \$verbose,
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
$evidence_info->primary_evidence($only_primary);
my $ids = [];
if($protein_file){
  $ids = $evidence_info->read_id_file($protein_file);
}else{
  push(@$ids, $protein_id);
}

foreach my $id(@$ids){
  $evidence_info->feature_ids_from_evidence_id($id, $gene_type);
}


sub perldocs{
  my ($msg) = @_;
  print $msg."\n" if($msg);
  exec('perldoc', $0);
  exit(0);
}
