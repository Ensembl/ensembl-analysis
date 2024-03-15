#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

=pod

=head1 NAME

load_assembly.pl

=head1 DESCRIPTION

This script downloads all the primary assembly data for the given species, creates the reference database for them and runs the load_seq_region script for each sequence level.

To be implemented:
Checks after loading seq regions.

Read ensembl-personal/genebuilders/sap_assembly_loading file for further information.

=head1 OPTIONS

-ftp_species_path	Species path to assembly main directory at NCBI FTP site (it usually starts with "/genbank/genomes/").

-local_path		Local path where the 'ftp_species_path' contents will be downloaded and the dumped out tables after creating the reference db will be saved. The script will create it if necessary.

-dbhost    		host name where the reference database will be created

-dbport    		what port to connect

-dbname    		name for the new reference db that will be created 

-dbuser    		what username to connect as

-dbpass    		what password to use

-prod_dbhost  production db host name

-prod_dbport  production db port

-prod_dbname  production db name 

-prod_dbuser  production db username

-prod_dbpass  production db password

-enscode	    path to the directories containing ensembl, ensembl-pipeline, ensembl-analysis. If not specified ENSCODE environment variable will be used

-ass_name	    Unique name that identifies the assembly. Check the genebuilders' sap_assembly_loading document if you do not know how to choose it.

-ena			    Use this option if you want the get_ENA_sequences.pl script to fetch the contigs from the ENA ftp site. If this is not used, the contigs will be fetched from the NCBI ftp site.

-ena_tax		  the taxonomy code (mam, rod, mus, hum or vrt)

-ena_id			  the four letter code of the WGS project (http://www.ncbi.nlm.nih.gov/genome/assembly/, search your assembly name)

-genebuild_id overwrite the lookup of genebuilder_id and use the one passed as an argument here

-skip			    Number of steps to be skipped. 0 by default. For example: if skip=2, the script will begin by running the step 3 (download contigs from NCBI).

-only			    Run only the specified step. For example, it can be set like "-only 4" to run only step number 4. The 'skip' parameter will be ignored if set.

-until        Run up to and including this step number

=head2 Output options:

  -verbose    Use this option to get more print statements to follow
                        the script. Set to 0 (not verbose) by default to get
                        only the final summary.

	-help		    Show usage.

=head1 EXAMPLE USAGE

=head1

bsub -M 1000000 -R 'select[mem>1000] rusage[mem=1000]' -o dog_load_assembly.out -e dog_load_assembly.err "perl load_assembly.pl -ftp_species_path /genbank/genomes/Eukaryotes/vertebrates_mammals/Canis_lupus/CanFam3.1 -local_path /lustre/scratch101/sanger/cgg/CanFam3.1_test -dbhost genebuild1 -dbname cgg_dog_ref_test -dbuser ensadmin -dbpass ***  -enscode /nfs/users/nfs_c/cgg/src -ass_name CanFam3.1 -ena_tax mam -ena_id aaex -verbose"

=cut


use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use File::Copy;
use File::Basename;
use List::Util qw(sum);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::ConfigWriter;

my $ftp_species_path;
my $local_path;
my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;
my $prod_dbhost;
my $prod_dbuser;
my $prod_dbpass;
my $prod_dbport;
my $prod_dbname;
my $dir_enscode ;
my $dir_ensgbscripts ;
my $ass_name;
my $ena_tax;
my $ena_id;
my $genebuild_id;
my $skip = 0; # skip 0 steps by default
my $only = 0;
my $help = 0;
my $verbose = 0;
my $ftp_server = "https://ftp.ncbi.nlm.nih.gov";
my $primary_assembly_dir_name = "Primary_Assembly"; # primary assembly dir name both the remote and the local one
my $contigs_dir_name = "contigs"; # contigs local dir name (it will be created in the local_path)
my $contigs_accessions_filename = "contig_accs.txt";
my $contigs_increment = "10000"; # number of contigs that will be downloaded at a time by using the interval query
my $configdir;
my $until = 0;
my $ena = 0;

my $all_contigs_filename = "contig.fa"; # filename for the fasta file containing all the downloaded contigs once concatenated (IT MUST BE THE SAME AS THE ONE DEFINED IN ensembl-personal/genebuilders/scripts/get_ENA_sequences.pl)
my $DOWNLOAD_RETRIES = 4; # number of retries for a wget command

GetOptions(
  'ftp_species_path:s' => \$ftp_species_path, # species path to assembly main directory at NCBI FTP site (it usually starts with "/genbank/genomes/")
  'local_path:s'    => \$local_path, # local path where the 'ftp_species_path' contents will be downloaded together with the dumpped out tables after creating the reference db. The script will create it if necessary.
  'dbhost|host|h=s' => \$dbhost,
  'dbname|db|D=s'   => \$dbname,
  'dbuser|user|u=s' => \$dbuser,
  'dbpass|pass|p=s' => \$dbpass,
  'dbport|port|P=s' => \$dbport,
  'prod_dbhost=s'   => \$prod_dbhost,
  'prod_dbname=s'   => \$prod_dbname,
  'prod_dbuser=s'   => \$prod_dbuser,
  'prod_dbpass=s'   => \$prod_dbpass,
  'prod_dbport=s'   => \$prod_dbport,
  'enscode=s'       => \$dir_enscode,  
  'ass_name=s'      => \$ass_name, # unique name that identifies the assembly
  'ena_tax=s'       => \$ena_tax,  # the taxonomy code (mam, rod, mus, hum or vrt), always lower case
  'ena_id=s'        => \$ena_id,   # the four letter code of the WGS project, always lower case
  'genebuild_id=s'  => \$genebuild_id, # overwrite any attem,pt to lookup your genebuild id
  'skip=s'          => \$skip,     # number of steps to be skipped
  'only=s'          => \$only,     # step number to be executed
  'until=s'         => \$until,
  'configdir=s'     => \$configdir,
  'ena!'            => \$ena,      # option to fetch contigs from ENA ftp site (if not set=default=NCBI ftp)
  'help'            => \$help,
  'verbose'         => \$verbose
);

print $0, "\n";

if (!$ftp_species_path or !$local_path or
    !$dbhost or !$dbname or !$dbuser or !$dbpass or !$dbport or
    !$prod_dbhost or !$prod_dbname or !$prod_dbuser or !$prod_dbpass or !$prod_dbport or
    !$ass_name or !$ena_tax or !$ena_id or $help) {
    &usage;
    exit(1);
}

if( !$dir_enscode )
{
    $dir_enscode = $ENV{'ENSCODE'} ;
    if( !$dir_enscode )
    {
        print "Please specify ENSCODE directory\n" ;
        &usage;
        exit(1);
    }
    $dir_enscode =~ s'/$'' ;
}

my $dir_assembly_scripts = $dir_enscode."/ensembl-analysis/scripts/assembly_loading/";

#-------
# BEGIN
#-------

#add / at the end of the paths if it cannot be found to avoid possible errors
if (!($local_path =~ /\/$/)) {
  $local_path .= "/";
}

if (!($ftp_species_path =~ /\/$/)) {
  $ftp_species_path .= "/";
}

# if parameter 'only' is set, make sure that we 'skip' anything else
if ($only > 0) {
  $skip = 99999; # hopefully, we wont have 99999 steps in a genebuild
}

if (($skip >= 0) and ($until == 0)) {
  $until = 99999;
}

# always lower case for these values
$ena_tax = lc($ena_tax);
$ena_id = lc($ena_id);

# print the number of skipped steps
print("The first step will be skipped.\n") if ($verbose and $skip == 1 and $only == 0);
print("The $skip first steps will be skipped.\n") if ($verbose and $skip > 1 and $only == 0);
print "Everything after step $until will be omitted\n" if ($until >= 1);

# print the number of the only step executed
print("Only the step number $only will be executed.\n") if ($verbose and $only > 0);

# download the entire primary assembly directory
print("\n===== STEP 1: DOWNLOAD PRIMARY ASSEMBLY DIRECTORY FROM FTP =====\n") if ($verbose);
if (($skip < 1 and $until >= 1) or $only == 1) {
    download_ftp_dir($ftp_server,
                     $ftp_species_path.$primary_assembly_dir_name,
                     $local_path.$primary_assembly_dir_name,
                     $verbose);
}

# unzip every compressed file
print("\n===== STEP 2: UNZIP PRIMARY ASSEMBLY ZIPPED FILES =====\n") if ($verbose);
unzip($local_path.$primary_assembly_dir_name,$verbose) if (($skip < 2 and $until >= 2) or $only == 2);

# download contigs
if (($skip < 3 and $until >= 3) or $only == 3) {
	print("\n===== STEP 3: DOWNLOAD CONTIGS FROM FTP =====\n") if ($verbose);
	download_ftp_contigs($ena,
				$ena_tax,
				$ena_id,
				$dir_assembly_scripts,
				$local_path.$primary_assembly_dir_name,
				$local_path.$contigs_dir_name,
				$all_contigs_filename,
				$verbose);

	# if the contig.fa file contains sequences from the ENA ftp site,
	# it is necessary to do the following headers fix
	# it is always necessary to fix it, even if the ncbi flag is used because we always try to
	# fetch some missing sequences at the last stage of get_ENA_sequences.pl script
	print("\n===== FIX ENA CONTIG HEADERS =====\n") if ($verbose);
	&fix_ENA_contig_headers($local_path.$contigs_dir_name."/".$all_contigs_filename,
				  $verbose);
}

# create the reference db where the downloaded assembly will be loaded
print("\n===== STEP 4: CREATE REFERENCE DB =====\n") if ($verbose);
create_ref_db($dbhost,$dbport,$dbuser,$dbpass,$dbname,
              $prod_dbhost,$prod_dbport,$prod_dbuser,$prod_dbpass,$prod_dbname,
              $dir_enscode,$local_path,$verbose) if (($skip < 4 and $until >= 4) or $only == 4);

print("\n===== STEP 5: LOAD SEQ_REGIONS =====\n") if ($verbose);
my @agp_filepaths = ();
@agp_filepaths =
load_seq_regions($dbhost,$dbport,$dbuser,$dbpass,$dbname,
                 $dir_enscode,
                 $local_path.$primary_assembly_dir_name,
                 $local_path.$contigs_dir_name."/".$all_contigs_filename,
                 $ass_name,
                 $verbose)
                 if (($skip < 5 and $until >= 5) or $only == 5);

print("\n===== STEP 6: LOAD ASSEMBLY INFORMATION =====\n") if ($verbose);
# load assembly information
load_assembly_info($dbhost,$dbport,$dbuser,$dbpass,$dbname,
                   $dir_enscode,
                   $verbose,
                   $local_path.$primary_assembly_dir_name,
                   @agp_filepaths)
                   if (($skip < 6 and $until >= 6) or $only == 6);

# run top level script
print("\n===== STEP 7: RUN TOP LEVEL SCRIPT =====\n") if ($verbose);
&top_level($dbhost,$dbport,$dbuser,$dbpass,$dbname,
          $dir_enscode,
          $verbose)
          if (($skip < 7 and $until >= 7) or $only == 7);

# Dump the top level sequence and compare with the downloaded version
print("\n===== STEP 8: CHECK: COMPARE THE TOP LEVEL SEQUENCE WITH THE DOWNLOADED VERSION =====\n") if ($verbose);
&check_top_level($dbhost,$dbport,$dbuser,$dbpass,$dbname,
                $dir_enscode,
                $local_path.$primary_assembly_dir_name,
                $verbose)
                if (($skip < 8 and $until >= 8) or $only == 8);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host => $dbhost,
    -port => $dbport,
    -user => $dbuser,
    -pass => $dbpass,
    -dbname => $dbname);
if (($skip < 9 and $until >= 9) or $only == 9) {
    print("\n===== STEP 9: UPDATING META TABLE, SYNONYMS AND KARYOTYPE RANK =====\n") if ($verbose);
    add_meta_table_values($genebuild_id);
    add_seq_region_synonym();
}

if (($skip < 10 and $until >= 10) or $only == 10) {
    if ($configdir) {
        print("\n===== STEP 10: COPYING CONFIG FILES =====\n") if ($verbose);
        copying_core_config_files();
    }
}

#load_taxonomy();
#load_meta_data();
#nt_transform_check(); # if 3 coord systems;

#-------
# END
#-------

=head2 download_ftp_dir

  Arg [1]   : string, ftp server
  Arg [2]   : string, ftp path to assembly's main directory at NCBI FTP site
  Arg [3]   : string, local path where [2] contents will be downloaded 
  Arg [4]   : verbose flag
  Function  : Connects to the ftp in Arg [2], downloads all the contents in Arg [2] to [3] by using wget.
  Exceptions: some thrown if given the wrong types
  Example   : download_ftp_dir("ftp.ncbi.nlm.nih.gov","/genbank/genomes/Eukaryotes/vertebrates_mammals/Canis_lupus/CanFam3.1","/nfs/ensembl/cgg/dog/CanFam3.1/download_test/Primary_Assembly",1);

=cut

sub download_ftp_dir {
  my ($ftp_server,$ftp_path,$local_dir,$verbose) = @_;

  my $wget_verbose = "-nv";

  my @dirs = split(/\//,$ftp_path); # get array of names of subdirectories
  my @cleaned_dirs = map {$_ ? $_ : ()} @dirs; # delete empty elements in the array
  my $numDirs = @cleaned_dirs; # count the number of dirs we need to skip to be able to put the files in the local path specified as parameter

  if (system("wget --no-proxy $wget_verbose -r -nH --cut-dirs=$numDirs --reject *.rm.out.gz -P $local_dir $ftp_server$ftp_path")) {
      throw("could not download the AGP, FASTA and info files\n");
  }

  # check if no file was downloaded
  my $num_files = int(`find $local_dir -type f | wc -l`);
  if ($num_files == 0) {
    throw("No file was downloaded from $ftp_server$ftp_path to $local_dir. Please, check that both paths are valid.");
  } else {
    print("$num_files files were downloaded.\n") if ($verbose);
  }
  
  # get the name of the directory where the assembly report file is
  my @ftp_path_array = split('/',$ftp_path);
  my $ass_report_dir = pop(@ftp_path_array); # Primary assembly dir
  $ass_report_dir = pop(@ftp_path_array);
  $ass_report_dir = pop(@ftp_path_array);
  
  if (system("wget $wget_verbose -nH -P $local_dir $ftp_server$ftp_path/../../$ass_report_dir"."_assembly_report.txt -O $local_dir/assembly_report.txt")) {
    throw("Could not download *_assembly_report.txt file from $ftp_server$ftp_path/../../ to $local_dir. Please, check that both paths are valid.");
  }
  else {
    print("Assembly report file was downloaded.\n") if ($verbose);
  }
}

sub unzip {
  my ($path,$verbose) = @_;
  my $gunzip_verbose = "";
  $gunzip_verbose = "-v" if ($verbose);
  print ("Unzipping the compressed files...\n") if ($verbose);
  throw("gunzip operation failed. Please check your error log file.") if (system("gunzip $gunzip_verbose -r $path") == 1);
  print ("Unzipping finished!\n") if ($verbose);
}

sub download_ftp_contigs {
  my ($ena,$ena_tax,$ena_id,$dir_assembly_scripts,$ass_path,$output_path,$all_contigs_filename,$verbose) = @_;

  my $ena_option = '-ncbi';
  $ena_option = '-ena' if ($ena);

  print("The contigs will be downloaded from the $ena_option- ftp site.\n") if ($verbose);

  # check if the output dir for contigs exists; otherwise, create it
  if (-e "$output_path") {
    print("Output path $output_path found.\n") if ($verbose);
  } else {
    `mkdir -p $output_path`;
    if (-e "$output_path") {
      print("Output path $output_path not found.\n$output_path created successfully.\n") if ($verbose);
    } else {
      throw("Cannot create output path for contigs $output_path\n");
    }
  }

  my $get_sequence_script = '/get_ENA_sequences.pl';
  if( -e $dir_assembly_scripts.$get_sequence_script ) {
    $get_sequence_script = $dir_assembly_scripts.$get_sequence_script;
  } else {
    throw("Couldn't find ".$get_sequence_script." in ".$dir_assembly_scripts."!");
  }
  print("Running script $get_sequence_script....\n");

  `perl $get_sequence_script -tax $ena_tax -id $ena_id -agp $ass_path $ena_option -outdir $output_path -no_desc -verbose -n_correction > $output_path/get_ENA_sequences.out`;
  throw($get_sequence_script.' died') if ($?);

  print("The output file of the get_ENA_sequences.pl script was written to $output_path/get_ENA_sequences.out\n");

  # check that some contigs were downloaded
  my $num_contigs = int(`grep -c ">" $output_path/$all_contigs_filename`);
  if ($num_contigs <= 0) {
    throw("No contig has been downloaded. Please, check your parameters (especially the ENA-related ones) and the ENA FTP site.");
  } else {
    print("The contigs FASTA file $all_contigs_filename was saved in $output_path. $num_contigs contigs were downloaded.\n") if ($verbose);
  }
}

sub fix_ENA_contig_headers() {
  my ($all_contigs_filepath,$verbose) = @_;

  if (-e $all_contigs_filepath) {
    print("File $all_contigs_filepath found. Fixing headers... ") if ($verbose);

    my $contigs_before_mv = int(`grep -c ">" $all_contigs_filepath`);

    # 1st condition (if): fixes headers fetched from ENA SVA like '>ENA|AC011295|AC011295.3' by replacing them with 'AC011295.3'
    # 2nd condition (elsif): fixes old-style headers, which don't appear these days (they got fixed in get_ENA_sequences.pl) but it doesn't harm
    `cat $all_contigs_filepath | perl -e 'while (<>) {if ( /\\>\\S+\\|\\S+\\|(\\S+).*/ ) {print ">\$1\n";} elsif ( /\\>(\\S+)/ ) {print ">\$1\n";} else {print;}}' > $all_contigs_filepath.fixed`;

    print("Finished!\n") if ($verbose);

    my $contigs_after_mv = int(`grep -c ">" $all_contigs_filepath.fixed`);

    if ($contigs_after_mv != $contigs_before_mv) {
      throw("The number of contigs before fixing the headers is different from the number of contigs after the fix. Please check your files $all_contigs_filepath and $all_contigs_filepath.fixed");
    }
    else {
      # Replace the all_dog_contigs.fa file with the *_header_fixed.fa file so we don't keep both versions (takes up a lot of space!):
      print("The number of contigs before fixing headers is $contigs_before_mv and after fixing it is $contigs_after_mv\n");
      `mv $all_contigs_filepath.fixed $all_contigs_filepath`;

      print("File $all_contigs_filepath.fixed moved to $all_contigs_filepath\n") if ($verbose);
    }
  } else {
    throw("File $all_contigs_filepath cannot be found.");
  }
}

=head2 create_ref_db

  Arg [1]   : string, db host
  Arg [2]   : string, db port
  Arg [3]   : string, db user
  Arg [4]   : string, db pass
  Arg [5]   : string, db name for the new ref db
  Arg [6]   : string, path to the dir where your local ensembl code root dir is
  Arg [7]   : string, path where the tables will be dumped out
  Function  : create reference db according to the parameters specified above
  Exceptions: some thrown if given the wrong types and cannot execute the system and mysql commands
  Example   : create_ref_db("genebuild1","3306","user","pass","cgg_dog_ref","/nfs/users/nfs_c/cgg/enscode",$local_path);

=cut
sub create_ref_db {
  my ($host,$port,$user,$pass,$dbname,
      $prod_host,$prod_port,$prod_user,$prod_pass,$prod_dbname,
      $enscode,$dump_path,$verbose) = @_;

  my $prod_db_check = `mysql -h$prod_host -P$prod_port -u$prod_user -p$prod_pass -NB -e"show databases like '$prod_dbname'"`;
  if ($prod_db_check =~ /$dbname/) {
    throw("Cannot find production database. The database $prod_dbname does not exist on host $prod_host port $prod_port");
  }

  # CREATE DB--------------------------------------------------------------
  # check if db exists
  my $dbname_check = `mysql -h$host -P$port -u$user -p$pass -NB -e"show databases like '$dbname'"`;
  if ($dbname_check =~ /$dbname/) {
    throw("Cannot create reference database. The database $dbname already exists on host $host port $port");
  }

  # create db
  `mysql -h$host -P$port -u$user -p$pass -e'create database $dbname'`;

  # check if the db was created
  $dbname_check = `mysql -h$host -P$port -u$user -p$pass -NB -e"show databases like '$dbname'"`;
  if ($dbname_check =~ /$dbname/) {
    print("The reference database $dbname was successfully created on host $host, port $port\n") if ($verbose);
  } else {
    throw("The reference database $dbname could not be created on host $host, port $port");
  }

  # CREATE MAIN TABLES---------------------------------------------------
  # check if file table.sql can be found
  if (!-e $enscode.'/ensembl/sql/table.sql') {
    throw("Cannot create reference database main tables. File $enscode/ensembl/sql/table.sql cannot be found.");
  }

  # create main tables
  `mysql -h$host -P$port -u$user -p$pass -D$dbname < $enscode/ensembl/sql/table.sql`;

  # check if the main tables were created successfully
  my $num_tables_created = int(`mysqlshow --count -h$host -P$port -u$user -p$pass $dbname | grep rows | awk '{print \$1}'`);
  # get the number of tables from the last line of the query
  # it should be a line like "78 rows in set.", so doing "int" will return only the number
#  my $num_tables_created = int($brand_new_tables[@brand_new_tables-2]);
  my $num_tables_to_create = int(`grep -ci "CREATE TABLE" $enscode/ensembl/sql/table.sql`);

  if ($num_tables_created != $num_tables_to_create) {
    throw("The reference database main tables could not be created. The number of created tables differs from the number of tables to be created. Number of created tables=$num_tables_created; Number of tables to be created=$num_tables_to_create");
  } else {
    print("$num_tables_created reference database main tables were created successfully.\n") if ($verbose);
  }

  # CREATE PIPELINE-RELATED TABLES-----------------------------------------
  # check if file table.sql can be found
  if (!-e $enscode.'/ensembl-pipeline/sql/table.sql') {
    throw("Cannot create reference database pipeline-related tables. File $enscode/ensembl-pipeline/sql/table.sql cannot be found.");
  }

  # create pipeline-related tables
  `mysql -h$host -P$port -u$user -p$pass -D$dbname < $enscode/ensembl-pipeline/sql/table.sql`;

  # check if the pipeline-related tables were created successfully
  my $num_pip_tables_created = `mysqlshow --count -h$host -P$port -u$user -p$pass $dbname | grep rows | awk '{print \$1}'`-$num_tables_created;

  # get the number of tables from the last line of the query
  # it should be a line like "78 rows in set.", so doing "int" will return only the number
  my $num_pip_tables_to_create = int(`grep -ci "CREATE TABLE" $enscode/ensembl-pipeline/sql/table.sql`);

  if ($num_pip_tables_created != $num_pip_tables_to_create) {
    throw("The reference database pipeline-related tables could not be created. The number of created tables differs from the number of tables to be created. Number of created tables=$num_pip_tables_created; Number of tables to be created=$num_pip_tables_to_create");
  } else {
    print("$num_pip_tables_created reference database pipeline-related tables were created successfully.\n") if ($verbose);
  }

  # POPULATE PRODUCTION DB TABLES
  # create dump_path
  `mkdir -p $dump_path`;
  if (!-e $dump_path) {
    throw("Cannot create dump_path $dump_path before populating production db tables");
  }

  # populate production db tables
  my $populate_production_tables_script = '/ensembl-production/scripts/production_database/populate_production_db_tables.pl';
  if (-e $enscode.$populate_production_tables_script) {
    `perl $enscode$populate_production_tables_script -h $host -P $port -u $user -p $pass -d $dbname -mh $prod_host -md $prod_dbname -mu $prod_user -mp $prod_pass -mP $prod_port -dp $dump_path`;
  } else {
    throw("Could not find the production script: ".$populate_production_tables_script." in\n ".$enscode."\n");
  }

  my $num_populated_tables = int(`mysqlshow --count -h$host -P$port -u$user -p$pass $dbname | awk '{if (int(\$6) > 0) print \$6}' | wc -l`);

  my $NUM_PRODUCTION_DB_TABLES = 5;
  if ($num_populated_tables == $NUM_PRODUCTION_DB_TABLES) {
    print("The production db tables were populated successfully.\n") if ($verbose);
  } else {
    throw("The production db tables could not be populated. Number of populated tables=$num_populated_tables. Expected number of populated tables=$NUM_PRODUCTION_DB_TABLES.");
  }
}

=head2 load_seq_regions

  Arg [1]   : string, db host
  Arg [2]   : string, db port
  Arg [3]   : string, db user
  Arg [4]   : string, db pass
  Arg [5]   : string, db name for the new ref db
  Arg [6]   : string, path to the dir where your local ensembl code root dir is
  Arg [7]   : string, path to the "Primary_Assembly" local path
  Arg [8]   : string, path to the all_contigs FASTA file (including full filename)
  Arg [9]   : string, name for the assembly which will be used as coord_system_version when loading the regions.
  Function  : load assembly seq regions into the given database
  Returns   : array of strings, filepaths to all of the AGP files needed for the load_agp.pl script later
  Exceptions: some thrown if given the wrong types
  Example   : load_seq_regions("genebuild1","3306","user","pass","cgg_dog_ref","/nfs/users/nfs_c/cgg/enscode","/nfs/ensembl/cgg/dog/CanFam3.1/download/assembly/Primary_Assembly","/nfs/ensembl/cgg/dog/CanFam3.1/download/assembly/contigs/contig.fa","CanFam3.1");

=cut
sub load_seq_regions {
  my ($host,$port,$user,$pass,$dbname,$enscode,$ass_path,$contigs_path,$ass_name,$verbose) = @_;

  # remove possibly (previous executions) created agp files to avoid inappropriate concatenation
  `rm $ass_path/*_all.agp`;
  `rm $ass_path/*_all.comp.agp`;

  print("\nLooking for and concatenating non-comp AGP files...\n");
  # for each subdir in Primary_Assembly (there should be some of assembled_chromosomes,placed_scaffolds,
  # unplaced_scaffolds,unlocalized_scaffolds) look for any .agp file excluding those whose name contains
  # substring "comp" concatenate them into assembled_chromosomes_all.agp, scaffolds_all.agp. These will
  # be saved into the $ass_path given

  # assuming bourne shell (sh)
  # first chromosomes
  system(
<<COMMAND
    for ass_subdir in `ls -d $ass_path/*chromosomes/`; do
      echo "ass path: "$ass_path
      echo "ass sub dir: "\$ass_subdir
      for agpfile in `ls \$ass_subdir"AGP"/*.agp | grep -v comp`; do
        echo "agpfile: "\$agpfile
        cat \$agpfile >> $ass_path/`basename \$ass_subdir`_all.agp
      done;
      grep '>' \$ass_subdir"FASTA"/*.fna | wc -l > $ass_path/`basename \$ass_subdir`_all.agp.fastawc
    done;
COMMAND
  );
  # and now scaffolds
  system(
<<COMMAND
    for ass_subdir in `ls -d $ass_path/*_scaffolds/`; do
      echo "ass path: "$ass_path
      echo "ass sub dir: "\$ass_subdir
      for agpfile in `ls \$ass_subdir"AGP"/*.agp | grep -v comp`; do
        echo "agpfile: "\$agpfile
        cat \$agpfile >> $ass_path/scaffolds_all.agp
      done;
      grep '>' \$ass_subdir"FASTA"/*.fna | wc -l > $ass_path/scaffolds_all.agp.fastawc
    done;
COMMAND
  );


  # load *_all.agp files into the database
  # Note that 'ls' is executed with -r (reverse) option so that the 'assembled_chromosomes' dir is the last one to be processed in the foreach loop. In this way, we will get the seq_regions loaded in the desired order as stated in the 'sop_assembly_loading' document after the genebuilder convention meeting. Note that the loading order will change if the dir names do.

  my $agp_filenames = `ls -r $ass_path/*_all.agp`; # note that the comp file is not required here so it is being excluded
  my @agp_filenames_arr = split("\n",$agp_filenames);

  if (@agp_filenames_arr <= 0) {
    throw("The *_all.agp files could not be created in $ass_path");
  } elsif ($verbose) {
    print("\nThe following *_all.agp files were created in $ass_path\n\n");
    print("$agp_filenames\n");
  }

  # decide which will be the rank 1
  # the chromosome if available; otherwise, the scaffold
  my $rank = 2;
  if (basename($agp_filenames_arr[@agp_filenames_arr-1]) =~ /chr/) {
    # the last string in the array will be the filename for the chromosomes if available.
    # chromosomes agp file found, so rank 1 for chromosome level
    $rank = 3;
    print("Chromosome agp file found. Chromosome level will be of rank 1, scaffold level will be of rank 2 and contig level will be of rank 3.\n") if ($verbose);
  } else {
    # chromosome agp file not found, so rank 1 for scaffold level
    print("Chromosome agp file not found. Scaffold level will be of rank 1 and contig level will be of rank 2.\n") if ($verbose);
  }

  # use the load seq_regions script to load the contigs
  print("\nLoading the contigs...\n") if ($verbose);
  print("Processing file $contigs_path\n") if ($verbose);

  my $contigs_path_basename = basename($contigs_path);
  `perl $enscode/ensembl-pipeline/scripts/load_seq_region.pl -dbhost $host -dbuser $user -dbpass $pass -dbport $port -dbname $dbname -coord_system_name contig -coord_system_version $ass_name -rank $rank -default_version -fasta_file $contigs_path -sequence_level -noverbose > $ass_path/load_seq_region_$contigs_path_basename.out`;

  print("The output file $ass_path/load_seq_region_$contigs_path_basename.out was written.\n");
  # maybe check for exceptions in the output file here

  # check contigs
  my $num_contigs = int(`mysql -h$host -P$port -u$user -p$pass -D$dbname -NB -e'select count(*) from seq_region'`);

  my $num_dna = int(`mysql -h$host -P$port -u$user -p$pass -D$dbname -NB -e'select count(*) from dna'`);

  if ($num_contigs != $num_dna) {
    throw("The number of 'seq_region' table loaded contigs ($num_contigs) differs from the number of 'dna' table rows ($num_dna)");
  } else {
    print("The number of 'seq_region' table loaded contigs ($num_contigs) is the same as the number of 'dna' table rows ($num_dna). Great!\n") if ($verbose);
  }

  $rank--;

  # use the load seq_regions script to load the scaffolds
  # for each existing scaffold agp file (unplaced, placed, unlocalized)
  print("\nLoading the scaffolds...\n") if ($verbose);
  for (my $i = 0; $i < @agp_filenames_arr-1; $i++) {
    my $curr_agp_path = $agp_filenames_arr[$i];
    my $curr_agp_basename = basename($curr_agp_path);

    print("Processing file $curr_agp_path\n") if ($verbose);

    `perl $enscode/ensembl-pipeline/scripts/load_seq_region.pl -dbhost $host -dbuser $user -dbpass $pass -dbport $port -dbname $dbname -coord_system_name scaffold -coord_system_version $ass_name -rank $rank -default_version -agp_file $curr_agp_path > $ass_path/load_seq_region_$curr_agp_basename.out`;

    print("The output file $ass_path/load_seq_region_$curr_agp_basename.out was written.\n");
    # maybe check for exceptions in the output file here

    # check scaffolds
    # maybe check fasta length here
  }

  $rank--;

  my @comp_agp_filenames_arr;
  if ($rank > 0) {
    # if we have chromosomes


    print("\nLooking for and concatenating comp AGP files...\n");
    # for .comp do the same as above (there should only be comp files for chromosome level)
    system(
<<COMMAND
      for ass_subdir in `ls -d $ass_path/*/`; do
        for agpfile in `ls \$ass_subdir"AGP"/*.agp | grep comp`; do
          cat \$agpfile >> $ass_path/`basename \$ass_subdir`_all.comp.agp
        done;
        grep '>' \$ass_subdir"FASTA"/*.fna | wc -l > $ass_path/`basename \$ass_subdir`_all.comp.agp.fastawc
      done;
COMMAND
    );

    # check if the comp agp file was created
    my $comp_agp_filenames = `ls -r $ass_path/*_all.comp.agp`;
    @comp_agp_filenames_arr = split("\n",$comp_agp_filenames);

    if (@comp_agp_filenames_arr <= 0) {
      throw("The *_all.comp_agp files could not be created in $ass_path");
    } elsif ($verbose) {
      print("\nThe following *_all.comp_agp files were created in $ass_path\n");
      print("$comp_agp_filenames\n");
    }

    my $chr_path = $agp_filenames_arr[@agp_filenames_arr-1];
    my $chr_path_basename = basename($chr_path);

    # use the load seq_regions script to load the chromosomes
    print("Loading the chromosomes...\n") if ($verbose);
    print("Processing file $chr_path\n") if ($verbose);

    `perl $enscode/ensembl-pipeline/scripts/load_seq_region.pl -dbhost $host -dbuser $user -dbpass $pass -dbport $port -dbname $dbname -coord_system_name chromosome -coord_system_version $ass_name -rank $rank -default_version -agp_file $chr_path > $ass_path/load_seq_region_$chr_path_basename.out`;

    print("The output file $ass_path/load_seq_region_$chr_path_basename.out was written.\n");
    # maybe check for exceptions in the output file here

    # maybe check fasta lengths here
    # maybe check warnings about non-ATGCN (RYKMSWBDHV) bases
  }

  # ought to update contig version to NULL
  `mysql -h$host -P$port -u$user -p$pass -D$dbname -e'update coord_system set version=NULL where name = "contig"'`;
  
  # check that it was updated
  my $contig_version = `mysql -h$host -P$port -u$user -p$pass -D$dbname -e'select version from coord_system where name = "contig"'`;

  if ($contig_version =~ /NULL/) {
    print("Column 'version' in 'coord_system' table set to NULL for contig level as required.\n") if ($verbose); 
  } else {
    throw("Column 'version' in 'coord_system' table could not be set to NULL for contig level.\n");
  }

  return @agp_filenames_arr,@comp_agp_filenames_arr;
}

sub load_assembly_info {
  my ($host,$port,$user,$pass,$dbname,$enscode,$verbose,$ass_path,@agp_filepaths) = @_;

  my $chr_level_exists = 0; # flag to check easily the meta table new rows at the end
  my $assembled;
  my $component;
  my @num_items_file;

  print("\nLoading assembly info...\n");
                                               
  # if we use the flags 'skip', 'until' and/or 'only', we need to check this just in case we had stopped the previous step and then resumed from this step
  # after editing the AGP files (this usually happens when dealing with the PAR in a new human assembly)
  if (scalar(@agp_filepaths) <= 0) {
    my $agp_filenames = `ls -r $ass_path/*_all*agp`; # note that the comp file is included here
    @agp_filepaths = split("\n",$agp_filenames);

    if (@agp_filepaths <= 0) {
      throw("The *_all*agp files could not be found in $ass_path");
    } elsif ($verbose) {
      print("\nThe following *_all*agp files were found in $ass_path\n\n");
      print("$agp_filenames\n");
    }
  }

  my $last_item_db = 0;
  # use the load_agp script to load the assembly info
  for (my $i = 0; $i < @agp_filepaths; $i++) {
    my $curr_agp_path = $agp_filepaths[$i];
    my $curr_agp_basename = basename($curr_agp_path);
    print("- Current AGP basename: $curr_agp_basename\n") if ($verbose);
    my $curr_agp_path_outfile = $curr_agp_path.".out";

    print("- Current AGP file: $curr_agp_path\n") if ($verbose);
    
    if ($curr_agp_basename =~ /scaf/) { # scaffold to contig
      $assembled = "scaffold";
      $component = "contig";
      print("Scaffold to contig detected.\n") if ($verbose);
    } elsif ($curr_agp_basename =~ /comp/) { # chromosome to contig
      $assembled = "chromosome";
      $component = "contig";
      $chr_level_exists = 1; # flag to check easily the meta table new rows at the end
      print("Chromosome to contig detected.\n") if ($verbose);
    } else { # chromosome to scaffold
      $assembled = "chromosome";
      $component = "scaffold";
      $chr_level_exists = 1; # flag to check easily the meta table new rows at the end
      print("Chromosome to scaffold detected.\n") if ($verbose);
    }
    `perl $enscode/ensembl-pipeline/scripts/load_agp.pl -dbhost $host -dbuser $user -dbpass $pass -dbport $port -dbname $dbname -assembled_name $assembled -component_name $component -agp_file $curr_agp_path >& $curr_agp_path_outfile`;

    print("$curr_agp_path_outfile output file written.\n") if ($verbose);

    # maybe look for exceptions and/or warnings in the output file here

    # check number of loaded items
    $num_items_file[$i] = int(`grep -v "#" $curr_agp_path | grep -v "\\WN\\W" | grep -c -v "\\WU\\W"`); #skipping gaps as not stored in db like in load_agp.pl
    my $num_items_db = int(`mysql -h$host -P$port -u$user -p$pass -D$dbname -NB -e'select count(*) from assembly'`);

#    my $index = $i-1 < 0 ? 0 : $i-1;
#    print STDERR $index, " hhh", $num_items_file[$i],"\n";
#    my $sum_items_file = sum($num_items_file[0..$index]);
#    my $num_items_loaded = $num_items_db-$sum_items_file;
print STDERR 'LOADED: ', $num_items_file[$i], ' ALREADY: ', $num_items_db, ' LAST: ', $last_item_db, "\n";
    my $num_items_loaded = $num_items_db-$last_item_db;
    $last_item_db = $num_items_db;
    if ($num_items_loaded != $num_items_file[$i]) {
      throw("num:$num_items_db sum:$last_item_db , The number of items loaded into the assembly table ($num_items_loaded) differs from the number of items in file ($num_items_file[$i]) for file $curr_agp_path");
    } else {
      print("The number of items loaded into the assembly table ($num_items_loaded) is the same as the number of items in file ($num_items_file[$i]) for file $curr_agp_path. Great!\n") if ($verbose);
    }

    # check warning messages in output file
    my $numWarnings = int(`grep "WARN" $curr_agp_path_outfile | wc -l`);
    if ($numWarnings > 0) {
      warning("There are some warning messages in $curr_agp_path_outfile that should be checked (grep 'WARN').\n");
    }

    my $numAlready = int(`grep "You are already using" $curr_agp_path_outfile | wc -l`);
    if ($numAlready > 0) {
      warning("There are some warning messages in $curr_agp_path_outfile that should be checked (grep 'You are already using').\n");
    }
  }

  # check meta table new rows for assembly.mapping
  my $numMetaRows = int(`mysql -h$host -P$port -u$user -p$pass -D$dbname -NB -e'select count(*) from meta where meta_key = "assembly.mapping"'`);
  if ($chr_level_exists) {
    if ($numMetaRows == 3) {
      print("3 new rows found in 'meta' table for 'assembly.mapping' meta_key. Great!\n") if ($verbose);
    } else {
      warning("$numMetaRows new rows found in 'meta' table for 'assembly.mapping' meta_key. 3 expected.");
    }
  } else { # if there is no chromosome level
    if ($numMetaRows == 1) {
      print("1 new row found in 'meta' table for 'assembly.mapping' meta_key. Great!\n") if ($verbose);
    } else {
      throw("$numMetaRows new rows found in 'meta' table for 'assembly.mapping' meta_key. 1 expected.");
    }
  }

  # print number of rows in assembly table
  my $numAssemblyRows = int(`mysql -h$host -P$port -u$user -p$pass -D$dbname -NB -e'select count(*) from assembly'`);
  print("\n$numAssemblyRows row(s) have been inserted into 'assembly' table\n");
}

sub top_level() {
  my ($host,$port,$user,$pass,$dbname,$enscode,$verbose) = @_;

  `perl $enscode/ensembl-pipeline/scripts/set_toplevel.pl -dbhost $host -dbuser $user -dbpass $pass -dbport $port -dbname $dbname`;

  my $numTopLevel = int(`mysql -h$host -P$port -u$user -p$pass -D$dbname -NB -e'select count(*) from seq_region_attrib where attrib_type_id = 6'`); # 6 corresponds to Top Level
  if ($numTopLevel > 0) {
    print("$numTopLevel Top Level Non-Redundant Sequence Regions found after running ensembl-pipeline/scripts/set_toplevel.pl\n");
  } else {
    throw("$numTopLevel Top Level Non-Redundant Sequence Regions found after running ensembl-pipeline/scripts/set_toplevel.pl");
  }
}

# dump the toplevel sequence and compare with the downloaded version (chromosomes, unplaced and unlocalized scaffolds)
sub check_top_level() {
  my ($host,$port,$user,$pass,$dbname,$enscode,$ass_path,$verbose) = @_;

  # remove to avoid problems
  # there should not be any additional .fa file
  # apart from the files in FASTA directories
  # otherwise, we would have problems with the 'find' command below
  `rm $ass_path/downloaded_toplevel.fa` if (-e "$ass_path/downloaded_toplevel.fa");
  `rm $ass_path/toplevel.fa` if (-e "$ass_path/toplevel.fa");
  `rm $ass_path/combined_toplevel.fa` if (-e "$ass_path/combined_toplevel.fa");

  # combine chromosomes and unplaced and unlocalized scaffolds which were downloaded
  print("\nCombining downloaded chromosomes and unplaced and unlocalized scaffolds (if available) into $ass_path/downloaded_toplevel.fa...\n") if ($verbose);

  # only the "non-placed" fasta files that are located in FASTA subdirs
  `find $ass_path -name '*.fna' -type f ! -name '*.placed.*' | grep FASTA | xargs cat >> $ass_path/downloaded_toplevel.fa`;

  my $numTopLevelDownloaded = int(`grep -c ">" $ass_path/downloaded_toplevel.fa`);

  print("\nDumping top level sequence...\n");
  # dump toplevel
  `perl $enscode/ensembl-analysis/scripts/sequence_dump.pl -dbhost $host -dbuser $user -dbpass $pass -dbport $port -dbname $dbname -toplevel -onefile -output_dir $ass_path/ > $ass_path/sequence_dump.out`;

  print("The output file $ass_path/sequence_dump.out was written.\n");
  
#  # check output file
#  my $sequence_dump_success = int(`wc -l $ass_path/sequence_dump.out`); # it should be empty if successfully completed
#  if ($sequence_dump_success == 0) {
#    print("The top level sequence was dumped successfully. File $ass_path/toplevel.fa created.\n") if ($verbose);
#  } else {
#    throw("It seems that the top level sequence was not dumped successfully. Please check the output file $ass_path/sequence_dump.out");
#  }

  # get number of top level sequences
  my $numTopLevel = int(`grep -c ">" $ass_path/toplevel.fa`);
  if ($numTopLevel == 0) {
    throw("$ass_path/toplevel.fa contains no sequence after dumping.");
  }

  if ($numTopLevel != $numTopLevelDownloaded) {
    throw("The number of dumped top level sequences ($numTopLevel) does not match the number of downloaded top level sequences ($numTopLevelDownloaded).");
  } else {
    print("The number of dumped top level sequences ($numTopLevel) matches the number of downloaded top level sequences ($numTopLevelDownloaded). Great!\n") if ($verbose);
  }

  # compare using fastanrdb
  # combine files
  print("Combining dumped top level file and downloaded top level file into $ass_path/combined_toplevel.fa...\n") if ($verbose);
  `cat $ass_path/toplevel.fa $ass_path/downloaded_toplevel.fa > $ass_path/combined_toplevel.fa`;

  my $numTopLevelCombined = int(`grep -c ">" $ass_path/combined_toplevel.fa`);
  my $sumTopLevel = $numTopLevel+$numTopLevelDownloaded;
  if ($numTopLevelCombined != $sumTopLevel) {
    throw("The number of combined top level sequences ($numTopLevelCombined) does not match the sum of dumped top level and downloaded top level sequences ($sumTopLevel).");
  } else {
    print("The number of combined top level sequences ($numTopLevelCombined) matches the sum of dumped and downloaded top level sequences ($sumTopLevel). Great!\n") if ($verbose);
  }

  # running fastanrdb to check that both sets of sequences are equivalent
  print("\nRunning fastanrdb...\n") if ($verbose);

  throw("fastanrdb died unexpectedly") if (system("fastanrdb $ass_path/combined_toplevel.fa > $ass_path/fastanrdb.out"));
  print("The output file $ass_path/fastanrdb.out was written.\n");

  my $numTopLevelFastanrdb = int(`grep -c ">" $ass_path/fastanrdb.out`);
  if ($numTopLevelFastanrdb != $numTopLevel) {
    warning("The number of sequences in file $ass_path/fastanrdb.out is $numTopLevelFastanrdb but it should be $numTopLevel. It means that we have duplicates.\n");
    print("The number of sequences in file $ass_path/fastanrdb.out is $numTopLevelFastanrdb but it should be $numTopLevel. It means that we have duplicates.\n");

    # maybe add more checks to inform about the sequences that does not match
    # cat $SCRA/fastanrdb.log | perl -ne 'next unless ($a,$b) = /^>(\w+\.\w+)\s*(\w+\.\w+)/; if ($a eq $b) {print "$a\n"}' | wc -l
    # to look at cases where NCBI and Ref DB agree
    # bsub -o $SCRA/fastanrdb_minus_r.log fastanrdb temp.fa -r
    # to compare sequences in reverse complement too

  } else {
    print("The number of sequences in file $ass_path/fastanrdb.out is $numTopLevelFastanrdb, which is good since it is the expected number, but there might be odd cases where duplicates can be present. More checks will be run.\n");
  }

  # check that the number of identifiers in headers matches the expected number of top level sequences*2
  my $numIds = int(`grep '>' $ass_path/fastanrdb.out | wc -w`);
  if ($numIds != $numTopLevel*2) {
    throw("The number of identifiers in headers in $ass_path/fastanrdb.out ($numIds) does not match the number_of_top_level_sequences*2 ($numTopLevel*2)");
  } else {
    print("The number of identifiers in headers in $ass_path/fastanrdb.out ($numIds) matches the number_of_top_level_sequences*2 ($numTopLevel*2). Great!\n") if ($verbose);
  }

  # check that there is no header with one identifier only
  my $numNonMatchingSeqs = int(`grep '>' $ass_path/fastanrdb.out | grep -c -v '\\w*[ ]\\w*'`);
  if ($numNonMatchingSeqs > 0) {
    throw("There are $numNonMatchingSeqs sequences with one identifier only in $ass_path/fastanrdb.out\nfastanrdb check FAILED!");
  } else {
    print("There is no one-identifier-only header in $ass_path/fastanrdb.out\nThe dumped top level sequence and the downloaded top level sequence are equivalent.\nfastanrdb check passed!\n");
  }
}
sub add_meta_table_values {
  my $genebuilder_id = shift;
  my $meta_adaptor = $db->get_MetaContainerAdaptor;
  my $genebuild_id_file;
  if (-e $dir_enscode.'/ensembl-production_private/release_coordination/docs/genebuilder_id.txt'){
    print $dir_enscode, "\n";
    $genebuild_id_file = $dir_enscode.'/ensembl-production_private/release_coordination/docs/genebuilder_id.txt';
  }
  if (defined $genebuilder_id) {
    $meta_adaptor->store_key_value('genebuild.id', $genebuilder_id);
  } elsif ($genebuild_id_file) {
    # this isn't going to work if your EBI and Sanger usernames differ - in that case you should have used the script arg -genebuild_id
    $genebuilder_id = `gawk '/$ENV{USER}/ {print \$1}' $genebuild_id_file`;
    $genebuilder_id =~ s/\s*(\d+)\s*/$1/;
    if ($genebuilder_id eq "") {
      warning("Could not find your genebuilder ID in genebuilder_id.txt, please add your id in the meta table!");
    } else {
      $meta_adaptor->store_key_value('genebuild.id', $genebuilder_id);
    }
  } else {
    warning "Could not find /ensembl-production_private/release_coordination/docs/genebuilder_id.txt, please add your id in the meta table!";
  } 
  $meta_adaptor->store_key_value('marker.priority', 1);
  $meta_adaptor->store_key_value('assembly.coverage_depth', 'high');
  if (-e $local_path.$primary_assembly_dir_name.'/assembly_report.txt') {
    open(RF, $local_path.$primary_assembly_dir_name.'/assembly_report.txt') || die("Could not open $local_path$primary_assembly_dir_name/assembly_report.txt");
    my $description_defined = 0;
    my $assembly_name;
    while (my $line = <RF>) {
      next if ($line !~ /^#/);
      if ($line =~ /^#\s*Date:\s*(\d+)-(\d+)-\d+/) {
        $meta_adaptor->store_key_value('assembly.date', $1.'-'.$2);
      }
      elsif ($line =~ /^#\s*Description:\s*(\S+)/) {
        $description_defined = 1;
        $meta_adaptor->store_key_value('assembly.name', $1);
      }
      elsif ($line =~ /^#\s*Assembly Name:\s*(\S+)/) {
        $assembly_name = $1;
        $meta_adaptor->store_key_value('assembly.default', $assembly_name);
      }
      elsif ($line =~ /^#\s*GenBank Assembly ID:\s*(\S+)/) {
        $meta_adaptor->store_key_value('assembly.accession', $1);
        $meta_adaptor->store_key_value('assembly.web_accession_source', 'NCBI');
        $meta_adaptor->store_key_value('assembly.web_accession_type', 'GenBank Assembly ID');
      }
    }
    if (!$description_defined) {
      $meta_adaptor->store_key_value('assembly.name', $assembly_name);
    }
  }
  else {
    print STDERR 'Could not find ', $local_path.$primary_assembly_dir_name.'/ASSEMBLY_INFO', "\n ==> No meta information will be inserted.\n";
  }
  print STDERR "Please check the meta table: SELECT * FROM meta ORDER BY meta_key; to be sure that it is correct!\n";
}

sub add_seq_region_synonym {
    #I'm not using the API here as I will modify the seq_region table if there is chromosomes
    if (-e $local_path.$primary_assembly_dir_name.'/assembled_chromosomes/chr2acc') {
        open(RF, $local_path.$primary_assembly_dir_name.'/assembled_chromosomes/chr2acc') || die("Could not open ".$local_path.$primary_assembly_dir_name.'/assembled_chromosomes/chr2acc');
        my $sth_select = $db->dbc->prepare('SELECT sr.seq_region_id FROM seq_region sr, coord_system cs WHERE cs.coord_system_id = sr.coord_system_id AND sr.name = ? AND cs.rank = 1');
        my $sth_insdc = $db->dbc->prepare('SELECT external_db_id FROM external_db WHERE db_name = "INSDC"');
        $sth_insdc->execute();
        my ($insdc_db_id) = $sth_insdc->fetchrow_array;
        # This is crap to hard code the external_db_id but it will more than probably never change as it is INSDC...
        my $sth_insert = $db->dbc->prepare('INSERT INTO seq_region_synonym (seq_region_id, synonym, external_db_id) VALUES(?, ?, ?)');
        my $sth_update = $db->dbc->prepare('UPDATE seq_region set name = ? WHERE seq_region_id = ?');
        while (my $line = <RF>) {
            next if ($line =~ /^#/);
            my ($synonym, $seq_region_name) = $line =~ /(\w+)\s+(\S+)/;
            $sth_select->bind_param(1, $seq_region_name);
            $sth_select->execute();
            my ($seq_region_id) = $sth_select->fetchrow_array();
            $sth_insert->bind_param(1, $seq_region_id);
            $sth_insert->bind_param(2, $seq_region_name);
            $sth_insert->bind_param(3, $insdc_db_id);
            $sth_insert->execute();
            $sth_update->bind_param(1, $synonym);
            $sth_update->bind_param(2, $seq_region_id);
            $sth_update->execute();
        }
        close(RF);
    }
    else {
        print STDERR 'Could not find ', $local_path.$primary_assembly_dir_name.'/assembled_chromosomes/chr2acc', "\n ==> No synonym will be inserted in seq_region_synonym for the chromosomes.\n";
    }
    foreach my $file ('component_localID2acc', 'scaffold_localID2acc') {
        next unless (-e $local_path.$primary_assembly_dir_name.'/'.$file);
        print STDERR "You will need to update the external_db_id for the synonyms of scaffold or contigs!\n";
        open(RF, $local_path.$primary_assembly_dir_name.'/'.$file) || die("Could not open ".$local_path.$primary_assembly_dir_name.'/'.$file);
        my $sth_select = $db->dbc->prepare('SELECT seq_region_id FROM seq_region WHERE name = ?');
        my $sth_insert = $db->dbc->prepare('INSERT INTO seq_region_synonym (seq_region_id, synonym) VALUES(?, ?)');
        while (my $line = <RF>) {
            next if ($line =~ /^#/);
            my ($synonym, $seq_region_name) = $line =~ /(\S+)\s+(\S+)/;
            if ($synonym eq 'na') {
                $synonym = $seq_region_name;
	    }
            $sth_select->bind_param(1, $seq_region_name);
            $sth_select->execute();
            my ($seq_region_id) = $sth_select->fetchrow_array();
            $sth_insert->bind_param(1, $seq_region_id);
            $sth_insert->bind_param(2, $synonym);
            $sth_insert->execute();
        }
    }
}

sub copying_core_config_files {
    # The ConfigWriter module forget all comment inside the %Config hash
    # this is why I'm using it only for Databases.pm as it's only useful
    # for this module. Also the formatting is a bit crap
    system("mkdir -p $configdir/Bio/EnsEMBL/{Analysis,Pipeline}/Config");
    system("mkdir -p $configdir/Bio/EnsEMBL/Analysis/Config/GeneBuild");
    foreach my $example_file ('/Bio/EnsEMBL/Analysis/Config/General.pm',
                              '/Bio/EnsEMBL/Analysis/Config/Blast.pm',
                              '/Bio/EnsEMBL/Analysis/Config/GeneBuild/KillListFilter.pm') {
        print 'Copying ', $configdir, $example_file, "...\n" if ($verbose);
        copy($dir_enscode.'/ensembl-analysis/modules'.$example_file.'.example', $configdir.$example_file);
    }
    foreach my $example_file ('/Bio/EnsEMBL/Pipeline/Config/General.pm',
                              '/Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm') {
        print 'Copying ', $configdir, $example_file, "...\n" if ($verbose);
        copy($dir_enscode.'/ensembl-pipeline/modules'.$example_file.'.example', $configdir.$example_file);
    }
    my $databases_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
        -modulename => 'Bio::EnsEMBL::Analysis::Config::Databases',
        -is_example => 1,
        );
    $databases_cfg->moduledir($configdir);
    $databases_cfg->key_by_parent('REFERENCE_DB', '-dbname', $dbname);
    $databases_cfg->key_by_parent('REFERENCE_DB', '-dbhost', $dbhost);
    $databases_cfg->key_by_parent('REFERENCE_DB', '-dbuser', $dbuser);
    $databases_cfg->key_by_parent('REFERENCE_DB', '-dbport', $dbport);
    $databases_cfg->key_by_parent('REFERENCE_DB', '-dbpass', $dbpass);
    print 'Copying ', $databases_cfg->get_module_path, "...\n" if ($verbose);
    $databases_cfg->write_config();
    print "Config files copied.\n";
}


sub usage {
    print <<EOF

Usage:

$0 -ftp_species_path <ftp_species_path> -local_path <local_path> -dbhost <dbhost> [-dbport <dbport>] -dbname <dbname> -dbuser <dbuser> -dbpass <dbpass> -enscode <enscode_base_directory> -ensgbscripts <genebuild_scripts_directory>  -ass_name <ass_name> [-ena] -ena_tax <ena_tax> -ena_id <ena_id> [-skip <skip>] [-only <only>] [-verbose] [-help]

-ftp_species_path	Species path to assembly main directory at NCBI FTP site (it usually starts with "/genbank/genomes/").

-local_path		Local path where the 'ftp_species_path' contents will be downloaded and the dumped out tables after creating the reference db will be saved. The script will create it if necessary.

-dbhost    		host name where the reference database will be created

-dbport    		what port to connect (default 3306)

-dbname    		name for the new reference db that will be created 

-dbuser    		what username to connect as

-dbpass    		what password to use

-prod_dbhost  production db host name

-prod_dbport  production db port

-prod_dbname  production db name 

-prod_dbuser  production db username

-prod_dbpass  production db password

-enscode      path to the directories containing ensembl, ensembl-pipeline, ensembl-analysis. If not specified ENSCODE environment variable will be used

-ensgbscripts path to the Genebuilder scripts directory. If not specified ENSGBSCRIPT envrionment variable will be used.

-ass_name		  Unique name that identifies the assembly. Check the genebuilders' sap_assembly_loading document if you do not know how to choose it.

-ena			    Use this option if you want the get_ENA_sequences.pl script to fetch the contigs from the ENA ftp site. If this is not used, the contigs will be fetched from the NCBI ftp site.

-ena_tax		  the taxonomy code (mam, rod, mus, hum or vrt)

-ena_id			  the four letter code of the WGS project (http://www.ncbi.nlm.nih.gov/genome/assembly/, search your assembly name)

-genebuild_id Use this argument to override the lookup of your genebuilder id

-skip			    Number of steps to be skipped. 0 by default. For example: if skip=2, the script will begin by running the step 3 (download contigs from NCBI).
              Please note that some manually cleaning could be necessary before re-running some steps.

-only			    Run only the specified step. For example, it can be set like "-only 4" to run only step number 4. The 'skip' parameter will be ignored if set.
              Please note that some manually cleaning could be necessary before re-running some steps.

-until        Run up to and including this step number. For example until=4 will run steps 1, 2, 3, and 4

-verbose     	Use this option to get more print statements to follow
              the script. Set to 0 (not verbose) by default to get
              only the final summary.

-help		     	Show usage.

Example:

bsub -M 1000 -R 'select[mem>1000] rusage[mem=1000]' -o zebrafish_load_assembly.out -e zebrafish_load_assembly.err "perl load_assembly.pl -ftp_species_path /genomes/all/GCA/000/002/035/GCA_000002035.4_GRCz11/GCA_000002035.4_GRCz11_assembly_structure/ -local_path OUTPUT/GRCz11_assembly_loading -dbhost genebuild1 -dbname gb_zebrafish_11 -dbuser *** -dbpass *** -prod_dbhost prod_host -prod_dbname production_db -prod_dbuser *** -prod_dbpass *** -enscode ENSEMBL_CODE_BASE -ass_name GRCz11 -ena -ena_tax vrt -ena_id ahzz -verbose"

EOF
}

1;
