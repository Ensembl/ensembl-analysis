#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;

$| = 1;

my $host         = '';
my $user         = '';
my $pass         = '';
my $port         = '';
my $dbname       = '';
my $gene_host    = '';
my $gene_user    = '';
my $gene_pass    = '';
my $gene_port    = '';
my $gene_dbname  = '';
my $out_host    = '';
my $out_user    = '';
my $out_pass    = '';
my $out_port    = '';
my $out_dbname  = '';
my $attrib_id = 347;#attrib_type id for attrib proj_alt_seq
my @patch_types = ('PATCH_FIX','PATCH_NOVEL');
my $out_dir = "";

&GetOptions( 'mapdbhost:s'       => \$host,
             'mapdbuser:s'       => \$user,
             'mapdbname:s'       => \$dbname,
             'mapdbpass:s'       => \$pass,
             'mapdbport:n'       => \$port,
             'genedbhost:s'      => \$gene_host,
             'genedbuser:s'      => \$gene_user,
             'genedbname:s'      => \$gene_dbname,
             'genedbpass:s'      => \$gene_pass,
             'genedbport:n'      => \$gene_port,
             'outdbhost:s'       => \$out_host,
             'outdbuser:s'       => \$out_user,
             'outdbname:s'       => \$out_dbname,
             'outdbpass:s'       => \$out_pass,
             'outdbport:n'       => \$out_port,
             'outdir:s'          => \$out_dir,
             'attrib_id:n'       => \$attrib_id);

my $mapdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                             -user   => $user,
                                             -pass   => $pass,
                                             -port   => $port,
                                             -dbname => $dbname );

my $genedb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -dnadb => $mapdb,
                                             -host   => $gene_host,
                                             -user   => $gene_user,
                                             -pass   => $gene_pass,
                                             -port   => $gene_port,
                                             -dbname => $gene_dbname);

my $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -dnadb => $mapdb,
                                             -host   => $out_host,
                                             -user   => $out_user,
                                             -pass   => $out_pass,
                                             -port   => $out_port,
                                             -dbname => $out_dbname);


my $ga = $genedb->get_GeneAdaptor();
my $outga =$outdb->get_GeneAdaptor();
my $geneta = $genedb->get_TranscriptAdaptor();

my $alt_pep_file = $out_dir."/altered_translations.sql";
my $alt_tscript_file = $out_dir."/altered_transcripts.sql";
open (PEP, ">", $alt_pep_file);
open (TSCRIPT, ">", $alt_tscript_file);

my @proj_genes = @{$outga->fetch_all};
print scalar(@proj_genes)." genes in projection db\n";

#for each projected gene
foreach my $pg (@proj_genes){
  my $gene_file = $out_dir."/".$pg->stable_id.".fa"; 
  open (GENE, ">", $gene_file);
  my $gene_transcript_changed = 0;
  my $gene_translation_changed = 0;
  my @transcripts = @{$pg->get_all_Transcripts};
  #for each transcript
  foreach my $t (@transcripts){
    my $stable_id = $t->stable_id;
    my $t_seq = $t->seq->seq;
    my $ref_seq = $geneta->fetch_by_stable_id($stable_id)->seq->seq;
    print $stable_id." patch\n";
    print $t_seq."\n";
    print $stable_id." ref\n";
    print $ref_seq."\n";
    if($t_seq eq $ref_seq){
      print $stable_id." ".$t->biotype." seqs are the same\n";
    }
    else{
      $gene_transcript_changed++;
      my $db_id = $t->dbID;
      print TSCRIPT "insert into transcript_attrib(transcript_id, attrib_type_id, value) values(".$db_id.", ".$attrib_id.",'');\n";
      print $stable_id." ".$t->biotype." seqs are different\n";
    }
    #check translation
    if($t->translation){
      my $t_pep = $t->translation->seq;
      my $ref_pep = $geneta->fetch_by_stable_id($stable_id)->translation->seq;
      print $t->translation->stable_id." patch\n";
      print $t_pep."\n";
      print $t->translation->stable_id." ref\n";
      print $ref_pep."\n";
      if($t_pep eq $ref_pep){
        print $t->translation->stable_id." ".$t->biotype." pep seq is the same\n";
      }
      else{
        $gene_translation_changed++;
        print $t->translation->stable_id." ".$t->biotype." pep seq is different\n";
        my $file = $out_dir."/".$t->translation->stable_id.".fa";
        my $db_id = $t->translation->dbID;
        print PEP "insert into translation_attrib(translation_id, attrib_type_id, value) values(".$db_id.", ".$attrib_id.",'');\n";
        open (OUT, ">", $file);
        print OUT ">ref_".$t->translation->stable_id."\n";
        print OUT $ref_pep."\n";
        print OUT ">patch_".$t->translation->stable_id."\n";
        print OUT $t_pep."\n";
        close(OUT);
        system("muscle -clw -in $file -out $file.aln");
        print GENE ">ref_".$t->translation->stable_id."\n";
        print GENE $ref_pep."\n";
        print GENE ">patch_".$t->translation->stable_id."\n";
        print GENE $t_pep."\n";
      }
    } else {
      print $stable_id." ".$t->biotype." has no pep\n";
    }

  }
close(GENE);
if($gene_transcript_changed){
  print $pg->stable_id." has $gene_transcript_changed altered transcript sequences\n"; 
}
if($gene_translation_changed){
  print $pg->stable_id." has $gene_translation_changed altered translations\n";
  system("muscle -clw -in $gene_file -out $gene_file.aln");  
}

}

print "Completed\n";
