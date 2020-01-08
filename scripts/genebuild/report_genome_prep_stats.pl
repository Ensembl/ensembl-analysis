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

use warnings;
use strict;
use feature 'say';

use Getopt::Long qw(:config no_ignore_case);

my $host;
my $port;
my $user;
my $dbname;
my $type;

GetOptions( 'host|dbhost|h=s'       => \$host,
            'port|dbport|P=s'       => \$port,
            'user|dbuser|u=s'       => \$user,
            'dbname|db|D=s'         => \$dbname,
            'report_type=s'         => \$type);

unless($type) {
  $type = 'all';
}

my $sql_start = "mysql -h".$host." -D".$dbname." -u".$user." -P".$port." -NB -e ";

if ($type eq "assembly_loading" || $type eq "all") {
  my $assembly_name = "'select version from coord_system where name=\"chromosome\"'";
  my $seq_region_count = "'select count(*) from seq_region'";
  my $chromosome_count = "'select count(*) from seq_region where coord_system_id=(select coord_system_id from coord_system where name=\"chromosome\")'";
  my $scaffold_count = "'select count(*) from seq_region where coord_system_id=(select coord_system_id from coord_system where name=\"scaffold\")'";
  my $contig_count = "'select count(*) from seq_region where coord_system_id=(select coord_system_id from coord_system where name=\"contig\")'";
  my $dna_count = "'select count(*) from dna'";
  my $toplevel_count = "'select count(*) from seq_region_attrib where attrib_type_id=6'";
  my $assembly_count = "'select count(*) from assembly'";
  my $repeatmasker_count = "'select count(*) from repeat_feature join analysis using (analysis_id) where logic_name like \"repeatmask%\"'";
  my $dust_count = "'select count(*) from repeat_feature where analysis_id=(select analysis_id from analysis where logic_name=\"dust\")'";
  my $trf_count = "'select count(*) from repeat_feature where analysis_id=(select analysis_id from analysis where logic_name=\"trf\")'";
  my $eponine_count = "'select count(*) from simple_feature where analysis_id=(select analysis_id from analysis where logic_name=\"eponine\")'";
  my $firstef_count = "'select count(*) from simple_feature where analysis_id=(select analysis_id from analysis where logic_name=\"firstef\")'";
  my $cpg_count = "'select count(*) from simple_feature where analysis_id=(select analysis_id from analysis where logic_name=\"cpg\")'";
  my $trnascan_count = "'select count(*) from simple_feature where analysis_id=(select analysis_id from analysis where logic_name=\"trnascan\")'";
  my $genscan_count = "'select count(*) from prediction_transcript'";
  my $vertrna_count = "'select count(*) from dna_align_feature where analysis_id=(select analysis_id from analysis where logic_name=\"vertrna\")'";
  my $unigene_count = "'select count(*) from dna_align_feature where analysis_id=(select analysis_id from analysis where logic_name=\"unigene\")'";
  my $uniprot_count = "'select count(*) from protein_align_feature where analysis_id=(select analysis_id from analysis where logic_name=\"uniprot\")'";


  say "==============================";
  say "Assembly Info: ";
  say "==============================";
  print "assembly name:      ".`$sql_start$assembly_name`;
  print "seq region count:   ".`$sql_start$seq_region_count`;
  print "chromosome count:   ".`$sql_start$chromosome_count`;
  print "scaffold count:     ".`$sql_start$scaffold_count`;
  print "contig count:       ".`$sql_start$contig_count`;
  print "dna count:          ".`$sql_start$dna_count`;
  print "toplevel count:     ".`$sql_start$toplevel_count`;
  print "assembly count:     ".`$sql_start$assembly_count`;

  print "\n";
  say "==============================";
  say "Feature Info: ";
  say "==============================";
  print "repeatmasker count: ".`$sql_start$repeatmasker_count`;
  print "dust count:         ".`$sql_start$dust_count`;
  print "trf count:          ".`$sql_start$trf_count`;
  print "eponine count:      ".`$sql_start$eponine_count`;
  print "firstef count:      ".`$sql_start$firstef_count`;
  print "cpg count:          ".`$sql_start$cpg_count`;
  print "trnascan count:     ".`$sql_start$trnascan_count`;
  print "genscan count:      ".`$sql_start$genscan_count`;
  print "vertrna count:      ".`$sql_start$vertrna_count`;
  print "unigene count:      ".`$sql_start$unigene_count`;
  print "unprot count:       ".`$sql_start$uniprot_count`;
}

if ($type eq "annotation" || $type eq "all") {
  my $gene_count = "'select count(*) from gene'";
  my $transcript_count = "'select count(*) from transcript'";
  my $exon_count = "'select count(*) from exon'";
  my $translation_count = "'select count(*) from translation'";
  my $exon_supporting_feature_count = "'select count(*) from supporting_feature'";
  my $gene_biotype_count = "'select count(*),biotype from gene group by biotype'";
  my $transcript_biotype_count = "'select count(*),biotype from transcript group by biotype'";

  print "\n";
  say "==============================";
  say "Annotation Info: ";
  say "==============================";
  print "gene count:                    ".`$sql_start$gene_count`;
  print "transcript count:              ".`$sql_start$transcript_count`;
  print "exon count:                    ".`$sql_start$exon_count`;
  print "translation count:             ".`$sql_start$translation_count`;
  print "exon supporting feature count: ".`$sql_start$exon_supporting_feature_count`;
  print "gene biotype count:\n".`$sql_start$gene_biotype_count`;
  print "transcript biotype count:\n".`$sql_start$transcript_biotype_count`;
}

exit;

