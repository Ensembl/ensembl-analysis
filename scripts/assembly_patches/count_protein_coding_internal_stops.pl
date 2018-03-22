#!/usr/bin/env perl

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use strict ;
use warnings ;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils 
  qw(contains_internal_stops);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils 
  qw(Gene_info);
	
my $host   = '';
my $port   = 3306;
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my $biotype = '';
my $dnadbname = 'sf7_patch_core_62';
my $dnadbhost = 'genebuild5';
my $gene_outfile = '';

my $internal_stop = 0;
my $processed = 0;

&GetOptions(
            'dbhost:s'   => \$host,
            'dbport:n'   => \$port,
            'dbname:s'   => \$dbname,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'dnadbname:s'=> \$dnadbname,
            'dnadbhost:s'=> \$dnadbhost,
	    'biotype:s'   => \$biotype,
	    'gene_outfile:s' => \$gene_outfile,
            );
if ($biotype eq '' || $gene_outfile eq '') {
  die "Outfile and biotype must be specified\n";
}
	   
open(GENE, ">", $gene_outfile) || die("can't open gene output file: $!");

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor
(
 -host   => $host,
 -user   => $dbuser,
 -port   => $port,
 -dbname => $dbname,
 -pass => $dbpass,
);   

# The DNA is in the ref db
my $dnadb =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new( -dbname => $dnadbname,
                                       -host   => $dnadbhost,
                                       -user   => $dbuser,
                                       -pass   => $dbpass,
                                       -port   => $port );

#print "Connected to ".$db->dbc->dbname." on ".$db->dbc->host."\n";
$db->dnadb($dnadb);

my @genes = @{ $db->get_GeneAdaptor->fetch_all_by_biotype($biotype) };

print "Number of genes retrieved: ".scalar(@genes)."\n";

foreach my $gene (@genes) {
	my $contains_stop_translations = 0;
	TRAN: foreach my $tr ( @{ $gene->get_all_Transcripts } ) {
    #print $tr->dbID."\n";
    next TRAN if $tr->biotype ne 'protein_coding';
    #print $tr->dbID." has biotype = protein_coding\n";
    $contains_stop_translations = contains_internal_stops($tr);
		if($contains_stop_translations){
			$internal_stop++;
      print GENE $contains_stop_translations." ".Gene_info($gene)."\n";
			#$contains_stop_translations = 1;
		}
		#print "Tested for internal stops\n";
		$processed++;
	}
}

print "Total transcripts processed: ".$processed."\n";
print "Transcripts with internal stops: ".$internal_stop."\n";

