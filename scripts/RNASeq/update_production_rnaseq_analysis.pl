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
=head1 NAME

  update_production_rnaseq_analysis.pl

=head1 SYNOPSIS

   perl update_production_rnaseq_analysis.pl \
    --dbhost hostname      \
    --dbname database_name \
    --version assembly version \
    --species genus_species    \
    --id your_production_database_id \
    --institute "RNASeq Institute" \ 
    --infile tab_delim_anlyses_file \
    --sqlfile output.sql 


=head1 DESCRIPTION

 Script to provide sql files to insert rows into the following 
 ensembl-production tables:

   1. output.sql -> analysis_description
   2. output.sql.web_data -> web_data
   3. output.sql.analysis_web_data -> analysis_web_data 

  It's important that the database is updated with the cotrolled vocab tables,
  analysis_description and web_data, before using the sql in output.sql.analysis_web_data
  since this last file performs sub selects on the 2 CV tables. 

=head1 EXAMPLE

  perl $ENSEMBL_ANAYSIS/scripts/RNASeq/update_production_rnaseq_analysis.pl --dbhost genebuild2 --dbname db8_anolis_carolinensis_rnaseq_71_2 \
    --version AnoCar2.0 --species anolis_carolinensis --id 102717 --institute "Broad Institute" --infile ~/anole_analyses.tab \
    --sqlfile ~/anole_add_rnaseq_to_production.sql 

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.
    
=cut


use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


# Connection to the target DB
my $host   = '';
my $port   = '3306';
my $user   = 'ensro';
my $dbname = '';

my $ephost   = 'mysql-ens-sta-1';
my $epport   = '4519';
my $epuser   = 'ensro';
my $epdbname = 'ensembl_production';

my $institute;
my $version;
my $write;
my $file;
my $sqlfile;
my $species_path;
my $creator;
my $help;

&GetOptions (
        'dbhost=s'    => \$host,
        'dbport=s'    => \$port,
        'dbuser=s'    => \$user,
        'dbname=s'    => \$dbname,
        'ephost=s'    => \$ephost,
        'epport=s'    => \$epport,
        'epuser=s'    => \$epuser,
        'epdbname=s'    => \$epdbname,      
        'infile=s'    => \$file,
        'species=s'   => \$species_path,
        'institute=s' => \$institute,
        'version=s'   => \$version,
        'id=i'        => \$creator, # ask web if you have never added anything to ensembl_production before
        'sqlfile=s'   => \$sqlfile,
        'help!'       => \$help,
        );

&Usage() if ( !($creator and $species_path and $file) or $help);

my $description_sqlquery = ''; # this will be seperate insert statements
# the other 2 tables are each populated with a single multi-row insert statement:
my $analysis_webdata_sqlquery = 'INSERT INTO analysis_web_data (analysis_description_id, web_data_id, species_id, db_type, displayable , created_by, created_at) VALUES';
my $webdata_sqlquery = 'INSERT INTO web_data (data, created_by, created_at) VALUES';

my $db_version;
my $web_data_id;

my %paired = (
  1 => 'Paired-end',
  0 => 'Single-end',
);

my %comments = (
  22 => 'Use for RNAseq introns',
  62 => 'Use for gene models built by our RNAseq pipeline',
  92 => 'Use for RNAseq BAM files',
);

my $species = ucfirst($species_path);
$species =~ s/_/ /g;
my $epdbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
  -host => $ephost,
  -port => $epport,
  -user => $epuser,
  -dbname => $epdbname,
  -driver => 'mysql'
);

my $epsth = $epdbc->prepare('SELECT species_id, web_name FROM species WHERE is_current = 1 AND production_name = ?');
$epsth->bind_param(1, $species_path);
$epsth->execute();
my ($species_id, $common_name) = $epsth->fetchrow_array();
$common_name = lcfirst($common_name);
my $sthcomment = $epdbc->prepare('SELECT comment FROM web_data WHERE web_data_id = ?');
foreach my $wdi (92, 62, 22) {
  $sthcomment->bind_param(1, $wdi);
  $sthcomment->execute();
  my ($comment) = $sthcomment->fetchrow_array();
  if ($comments{$wdi} ne $comment) {
    print STDERR 'Check the web_data for web_data_id ', $wdi, "\n";
  }
}

my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
  -host => $host,
  -port => $port,
  -user => $user,
  -dbname => $dbname,
  -driver => 'mysql'
);

my $sth = $dbc->prepare('SELECT df.name FROM data_file df LEFT JOIN analysis a ON df.analysis_id = a.analysis_id WHERE df.file_type = "BAMCOV" AND a.logic_name = ?');

my %seen_B;
my %seen_I;
my %seen_G;

open(RF, $file) || die("Could not open $file");

while (my $line = <RF>) {
  print "******************\n".$line;
  my ($logic_name, $is_paired, $length, $group, $row, $group_order) = $line =~ /(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/;
  my $group_ = $group;
  $group =~ s/_/ /g;
  my $row_ = $row;
  $row =~ s/_/ /g;
  my $display_label = ucfirst($logic_name);
  $display_label =~ s/_/ /g;
  my $description;
  

  # check for existing web_data or analysis_description.
  my $rows = check_for_existing_rows(sprintf 'SELECT web_data_id FROM web_data WHERE data LIKE \'%%"BAM files"%%"group" => "%s"%%"row" => "%s"%%\'',$group, $row);
  if ($rows != 0) {
    warn("Existing web_data for $logic_name - Not adding web_data.");
    $seen_B{$group_.$row_} = 1;
    $seen_I{$group_.$row_} = 1;
    $seen_G{$group_.$row_} = 1;
  }   

  $rows = check_for_existing_rows(sprintf 'SELECT analysis_description_id FROM analysis_description WHERE logic_name = "%s"', $logic_name);
  if ($rows != 0) {
    throw("Existing analysis description for $logic_name - check what already exists in the production database.");
  }

  if ($display_label =~ s/\s+bwa\s+/ /) {
    $display_label .= ' RNAseq alignments';
    $sth->bind_param(1, $logic_name);
    $sth->execute();
    my ($bamfile) = $sth->fetchrow_array();
    $description = sprintf('%s %i bp RNAseq reads from a %s. Data were obtained from %s and aligned to the %s assembly using BWA. Download this BAM file from the <a href="ftp://ftp.ensembl.org/pub/data_files/%s/%s/rnaseq/">Ensembl FTP</a> (file:<a href="ftp://ftp.ensembl.org/pub/data_files/%s/%s/rnaseq/%s.bam">%s.bam</a>).', $paired{$is_paired}, $length, $common_name, $institute, $version, $species_path, $version, $species_path, $version, $bamfile, $bamfile);
    my ($tissue) = $logic_name =~ /[a-z]+_bwa_(\S+)/;
    $tissue =~ s/_/ /g;

    $analysis_webdata_sqlquery .=  sprintf('((SELECT analysis_description_id FROM analysis_description WHERE logic_name = "%s"), (SELECT web_data_id FROM web_data WHERE data LIKE \'%%"BAM files"%%"group" => "%s"%%"row" => "%s"%%\'), %i, "rnaseq", 1, %i, NOW()),', $logic_name, $group, $row, $species_id, $creator)." \\\n";
    $webdata_sqlquery .= sprintf('(\'{"colour_key" => "bam", "label_key" => "RNASeq [biotype]", "matrix" => {"column" => "BAM files", "group_order" => "%i", "group" => "%s", "menu" => "rnaseq", "row" => "%s"}, "type" => "rnaseq", "zmenu" => "RNASeq_bam"}\', %d, NOW()),', $group_order, $group, $row, $creator)." \\\n" unless $seen_B{$group_.$row_};
    $seen_B{$group_.$row_} = 1;
  }                                                                                                                              
  elsif ($logic_name =~ /[a-z]+_(\w+)_rnaseq_introns/) {                                                                         
    $description = sprintf('Spliced RNASeq read support for %s for %s.', $common_name, $1);                                    
    $display_label =~ s/\s+rnaseq introns/ intron-spanning reads/;                                                                                
    $analysis_webdata_sqlquery .= sprintf('((SELECT analysis_description_id FROM analysis_description WHERE logic_name = "%s"), NULL, %i, "rnaseq", 1, %i, NOW()),', $logic_name, $species_id, $creator)." \\\n";
  }                                                                                                                              
  elsif ($logic_name =~ /[a-z]+_(\w+)_introns/) {                                                                                
    my $tissue = $1;
    $description = sprintf('Spliced RNASeq read support for %s in %s', $common_name, $tissue);                                      
    $display_label =~ s/\s+introns/ intron-spanning reads/;                                                                                
    $analysis_webdata_sqlquery .= sprintf('((SELECT analysis_description_id FROM analysis_description WHERE logic_name = "%s"), (SELECT web_data_id FROM web_data WHERE data LIKE \'%%"Intron-spanning reads"%%"group" => "%s"%%"row" => "%s"%%\'), %i, "rnaseq", 1, %i, NOW()),', $logic_name, $group, $row, $species_id, $creator)." \\\n";
    $tissue =~ s/_/ /g;
    $webdata_sqlquery .= sprintf('(\'{"additional_renderers" => ["histogram", "Variable height"],"colour_key" => "intron_support", "matrix" => {"column" => "Intron-spanning reads", "group" => "%s", "menu" => "rnaseq", "row" => "%s"}, "type" => "rnaseq", "zmenu" => "Supporting_alignment"}\', %d, NOW()),', $group, $row, $creator)." \\\n" unless $seen_I{$group_.$row_};
    $seen_I{$group_.$row_} = 1;
  }                                                                                                                              
  elsif ($logic_name =~ /[a-z]+_(\w+)_rnaseq/) {                                                                                 
    my $tissue = $1;
    $description = sprintf('Annotation generated using only RNASeq data from %s %s tissue.', $common_name, $tissue);                
    $display_label =~ s/\s+rnaseq/ RNASeq/;                                                                                               
    $analysis_webdata_sqlquery .= sprintf('((SELECT analysis_description_id FROM analysis_description WHERE logic_name = "%s"), (SELECT web_data_id FROM web_data WHERE data LIKE \'%%"Gene models"%%"group" => "%s"%%"row" => "%s"%%\'), %i, "rnaseq", 1, %i, NOW()),', $logic_name, $group, $row, $species_id, $creator)." \\\n";
    $tissue =~ s/_/ /g;
    $webdata_sqlquery .= sprintf('(\'{"colour_key" => "human_rnaseq", "label_key" => "RNASeq [biotype]", "matrix" => {"column" => "Gene models", "group" => "%s", "menu" => "rnaseq", "row" => "%s"}, "type" => "rnaseq", "zmenu" => "RNASeq"}\', %d, NOW()),', $group, $row, $creator)." \\\n" unless $seen_G{$group_.$row_};
    $seen_G{$group_.$row_} = 1;
  }
  $description_sqlquery .= sprintf('INSERT INTO analysis_description (logic_name, description, display_label, db_version, is_current, created_by, created_at) VALUES (\'%s\', \'%s\', \'%s\', 0, 1, \'%s\', NOW() );',
                                                                     $logic_name, $description, $display_label, $creator)."\n";
}

open(WF, ">$sqlfile") || die("Could not open $sqlfile");
print WF $description_sqlquery;
close(WF) || die("Could not close $sqlfile");

$analysis_webdata_sqlquery =~ s/, \\$/;/;
open(WF, ">$sqlfile.analysis_web_data") || die("Could not open $sqlfile.analysis_web_data");
print WF $analysis_webdata_sqlquery;
close(WF) || die("Could not close $sqlfile.analysis_web_data");

$webdata_sqlquery =~ s/, \\$/;/;
open(WF, ">$sqlfile.web_data") || die("Could not open $sqlfile.web_data");
print WF $webdata_sqlquery;
close(WF) || die("Could not close $sqlfile.web_data");


# when adding web_data or analyis_description check to see if any rows are already existing 
sub check_for_existing_rows {
  my $sth = $epdbc->prepare(shift);
  my $rows = $sth->execute();
  return $rows;
}



sub Usage {
    print <<EOF
 $0 -dbhost <dbhost> -dbuser <user> -dbname <bdname> [-dbport <port>] -id <your_id> -species <production db name> -infile <file> -institute <institute name> -version <assembly version> -sqlfile <file name>

     -dbhost         Host name of your rnaseq database
     -dbport         Host port, default is 3306
     -dbuser         User name, default is ensro
     -dbname         Name of your rnaseq database
     -infile         File with 3 columns: logic_name  paired-end(0/1)  length of the reads
     -species        Species name formatted like the production_name, i.e anolis_carolinensis
     -institute      Name of the institute/consortium who provided the data set
     -version        Assembly version, i.e AnoCar2.0
     -id             Your id in the production database. You can find it by looking at data you've added in the production database
     -sqlfile        The output file containing a MySQL query to insert your analysis description.
                      It creates 3 files: <file name> for the analysis_description table
                                          <file name>.web_data for the analysis_web_data table
                                          <file name>.analysis_web_data for linking analysis_description and web_data
     -help
EOF
;
  exit(0) if ($help);
  exit(1);
}
