# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

HiveLoadPatches.pm

=head1 DESCRIPTION

This module automates the GRC patch loading into an Ensembl database by running the following steps:

Load patches

Checks

=head1 OPTIONS

-ftp_path       NCBI FTP path where the alt_scaffolds directory is located.

-output_path    Output path where output and log files will be written.

-dbhost         database host name

-dbport         database port

-dbname         database name

-dbuser         database username to connect as

-dbpass         database password to use

-cs_version     Coordinate system version.

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPatches -ftp_path https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.23_GRCh38.p8/GCA_000001405.23_GRCh38.p8_assembly_structure/PATCHES/alt_scaffolds -output_path $SCR9/ -dbhost genebuildX -dbname core_97 -dbuser *** -dbpass *** -dbport DBPORT -cs_version GRCh38

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPatches;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);

sub param_defaults {
    return {
      ftp_path => undef,
      output_path => undef,
      dbhost => undef,
      dbname => undef,
      dbuser => undef,
      dbpass => undef,
      dbport => undef,
      cs_version => undef,
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;

  $self->param_required('ftp_path');
  $self->param_required('output_path');
  $self->param_required('dbhost');
  $self->param_required('dbport');
  $self->param_required('dbname');
  $self->param_required('dbuser');
  $self->param_required('dbpass');
  $self->param_required('cs_version');
  
  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('output_path') =~ /\/$/)) {
    $self->param('output_path',$self->param('output_path')."/");
  }

  # create output dir if it does not exist
  if (not -e $self->param('output_path')) {
    run_command("mkdir -p ".$self->param('output_path'),"Create output path.");
  }

  download_patches($self->param('ftp_path'),
                   $self->param('output_path'));
  
  remove_old_rawcomputes_karyobands($self->param('dbhost'),
                                    $self->param('dbport'),
                                    $self->param('dbname'),
                                    $self->param('dbuser'),
                                    $self->param('dbpass')); 
  
  my $patches_removed = remove_deprecated_patches($self->param('dbhost'),
                                                  $self->param('dbport'),
                                                  $self->param('dbname'),
                                                  $self->param('dbuser'),
                                                  $self->param('dbpass'),
                                                  $self->param('output_path'));
  
  if ($patches_removed) {
    recalculate_markers_map_weights($self->param('dbhost'),
                                    $self->param('dbport'),
                                    $self->param('dbname'),
                                    $self->param('dbuser'),
                                    $self->param('dbpass')); 
  }
  
  my ($alt_scaffold_fixed_filepath,$alt_scaffold_fna_fixed_filepath) = remove_patches_on_patches($self->param('dbhost'),
                                                                                                 $self->param('dbport'),
                                                                                                 $self->param('dbname'),
                                                                                                 $self->param('dbuser'),
                                                                                                 $self->param('dbpass'),
                                                                                                 $self->param('output_path'));
  
  load_patches($self->param('dbhost'),
               $self->param('dbport'),
               $self->param('dbname'),
               $self->param('dbuser'),
               $self->param('dbpass'),
               $self->param('output_path'),
               $self->param('ftp_path'),
               $alt_scaffold_fixed_filepath);
  
  check_patches($self->param('output_path'),
                $self->param('dbhost'),
                $self->param('dbport'),
                $self->param('dbname'),
                $self->param('dbuser'),
                $self->param('dbpass'),
                $alt_scaffold_fixed_filepath,
                $alt_scaffold_fna_fixed_filepath,
                $self->param('cs_version'));

}

sub download_patches() {
# download the 5 patches files from the ftp_path (ie 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_.../GCA_...assembly_structure/PATCHES/alt_scaffolds'):
# alt.scaf.agp.gz, alt.scaf.fna.gz, alt_scaffold_placement.txt, patch_type and assembly_report.txt
# to the local directory 'local_dir'
 my ($self) = @_;

  my ($ftp_path,$local_dir) = @_;

  my $wget_verbose = "-nv";

  my @dirs = split(/\//,$ftp_path); # get array of names of subdirectories
  my @cleaned_dirs = map {$_ ? $_ : ()} @dirs; # delete empty elements in the array
  my $numDirs = @cleaned_dirs; # count the number of dirs we need to skip to be able to put the files in the local path specified as parameter

  my $cmd = "wget --no-proxy ".$wget_verbose." -r -nH --cut-dirs=".$numDirs." --reject *.rm.out.gz --reject *.asn --reject *.gff -P ".$local_dir." ".$ftp_path;
  my $return = system($cmd);
  if ($return) {
    $self->throw("Could not download the AGP, FASTA and info files. Commandline used:\n".$cmd);
  }

  # get the name of the directory where the assembly report file is
  my @ftp_path_array = split('/',$ftp_path);
  my $ass_report_dir = pop(@ftp_path_array); # Primary assembly dir
  $ass_report_dir = pop(@ftp_path_array);
  $ass_report_dir = pop(@ftp_path_array);
  $ass_report_dir = pop(@ftp_path_array);

  # my $link  = $ftp_path."/../../../".$ass_report_dir."_assembly_report.txt";
  # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.8_GRCm38.p6/GCA_000001635.8_GRCm38.p6_assembly_report.txt
  my $link = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.8_GRCm38.p6/GCA_000001635.8_GRCm38.p6_assembly_report.txt';

  if (system("wget ".$wget_verbose." -nH -P ".$local_dir." https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.8_GRCm38.p6/GCA_000001635.8_GRCm38.p6_assembly_report.txt -O ".$local_dir."/assembly_report.txt")) {
       $self->throw("Could not download *_assembly_report.txt file from ".$ftp_path."/../../../ to ".$local_dir.". Please, check that both paths are valid.");
  }
  else {
       print("Assembly report file was downloaded\n");
  }
 # if (system("wget ".$wget_verbose." -nH -P ".$local_dir." ".$ftp_path."/../../../".$ass_report_dir."_assembly_report.txt -O ".$local_dir."/assembly_report.txt")) {
 #   $self->throw("Could not download *_assembly_report.txt file from ".$ftp_path."/../../../ to ".$local_dir.". Please, check that both paths are valid.");
 # }
 # else {
 #   print("Assembly report file was downloaded\n");
 # }
  
  # check if the 5 files were downloaded
  my $num_files = int(`find $local_dir -type f | wc -l`);
  if ($num_files != 5) {
    $self->throw("Files were not downloaded successfully from ".$ftp_path." to ".$local_dir.". Please, check that both paths are valid and your output path has not been used before.");
  } else {
    print("$num_files files were downloaded\n");
  }

  $cmd = "gunzip $local_dir/*.gz";
  $return = system($cmd);
  if ($return) {
    $self->throw("Could not gunzip the .gz files. Commandline used:\n".$cmd);
  } else {
  	print("Files gunzipped successfully\n");
  }

  # print the total number of patches for your records
  my $num_patches = run_command("grep -v '^#' $local_dir/patch_type | wc -l","Getting the number of patches for your records...");
  print("Total number of patches: $num_patches\n");
}

sub remove_old_rawcomputes_karyobands() {
# Remove old patch raw compute and karyotype bands

  my ($host,$port,$name,$user,$pass) = @_;

  run_command("mysql -h $host -P $port -u $user -p$pass $name < ".'$ENSEMBL_ANALYSIS/scripts/assembly_patches/remove_patch_raw_compute.sql',"Removing raw computes on patches...");
  
  run_command("mysql -h $host -P $port -u $user -p$pass $name < ".'$ENSEMBL_ANALYSIS/scripts/assembly_patches/remove_patch_karyotype.sql',"Removing karyotype bands on patches..."); # we decided not to have karyotype bands on patches because projections do not make sense so there should not be any but just in case they come from the old days...
}

sub remove_deprecated_patches() {
# Remove altered and deprecated patches
# Return 1 if any patch was removed. Otherwise, return 0.

  my ($host,$port,$name,$user,$pass,$output_dir) = @_;
  
  my $alt_patch_file = $output_dir."/alt_patch_output.txt";
  my $delete_patch_file = $output_dir."/delete_patch.sql";
  
  run_command('perl $ENSEMBL_ANALYSIS/scripts/assembly_patches/remove_updated_and_deprecated_patches.pl '."-user $user -pass $pass -host $host -port $port -dbname $name -patchtype_file $output_dir/patch_type -central_cs scaffold -sqlfile $delete_patch_file > $alt_patch_file","Getting the deprecated patches to remove...");

  my $res;
  
  $res = run_command("grep modified $alt_patch_file | cat","Getting the list of modified patches...");
  print("Modified patches: $res\n");
  
  $res = run_command("grep removed $alt_patch_file | cat","Getting the list of removed patches...");
  print("Removed patches: $res\n");
  
  $res = run_command("grep 'is a new patch' $alt_patch_file | cat","Getting the list of new patches...");
  print("New patches: $res\n");
  
  if (-e $delete_patch_file and !(-z $delete_patch_file)) { # exists and not empty
    run_command("mysql -h $host -P $port -u $user -p$pass $name < $output_dir/delete_patch.sql","Deleting the modified and removed patches...");
    return 1;
  } else {
    print("There is not any patch to delete.\n");
    return 0;
  }
}

sub recalculate_markers_map_weights() {
# Since markers were added to patches that went through a full build or a markers update in the past, if the patches containing markers have been deleted, the weights may be wrong for remaining markers and the MarkerFeatures healthcheck will fail. To fix that, we need to recalculate the markers map weights by running the following script in the specified database.

  my ($host,$port,$name,$user,$pass) = @_;
  
  run_command('perl $ENSEMBL_ANALYSIS/scripts/markers/unmapped_markers.pl '."-host $host -port $port -user $user -pass $pass -dbname $name -port $port -logic_name marker -max_duplicates 3","Recalculating markers map weights...");
}

sub remove_patches_on_patches() {
# Since patches on patches are not supported, this sub removes the corresponding lines in the alt_scaffold_placement.txt file and in the alt.scaf.fna file
# output_dir is the directory where the alt_scaffold_placement.txt file is located and where the fixed file will be written.
# Return the filepath to the file containing the file after removing the patches on pathes.
  my ($host,$port,$name,$user,$pass,$output_dir) = @_;
  
  my $alt_scaffold_file = $output_dir."/alt_scaffold_placement.txt";
  my $alt_scaffold_file_fixed = $alt_scaffold_file.".fixed";
  my $alt_scaffold_file_seqs_to_remove = $alt_scaffold_file.".seqs_to_remove";
  my $alt_scaffold_fna_file = $output_dir."/alt.scaf.fna";
  my $alt_scaffold_fna_file_fixed = $output_dir."/alt.scaf.fna.fixed";

  run_command("rm $alt_scaffold_file_seqs_to_remove | cat","Removing old $alt_scaffold_file_seqs_to_remove if exists...");
  run_command("rm $alt_scaffold_fna_file_fixed | cat","Removing old $alt_scaffold_fna_file_fixed if exists...");

  #my $cmd = "awk '{if (\$2 == \"Primary\" || \$1 == \"#alt_asm_name\") print}' $alt_scaffold_file > $alt_scaffold_file_fixed";
  my $cmd = "IFS=\$'\\n';
             set -f;
             for line in `cat $alt_scaffold_file`; do
               PARENTACCorNAME=`echo \$line | awk '{if(\$2 == \"C57BL/6J\") {print \$7} else {if(\$7 == \"parent_acc\") {print \$7} else {print \$6}}}'`;
               ACC=`echo \$line | awk '{if(\$2 == \"C57BL/6J\") {print \$4} else {print \$4}}'`;
               RES=`mysql -u$user \\
                          -p$pass \\
                          -h$host \\
                          -P$port \\
                          -D$name -NB \\
                          -e\"SELECT COUNT(*)
                              FROM seq_region_synonym
                              WHERE synonym=\\\"\$PARENTACCorNAME\\\"
                                AND seq_region_id NOT IN (SELECT seq_region_id
                                                          FROM assembly_exception
                                                          WHERE exc_type IN (\\\"PATCH_NOVEL\\\",\\\"PATCH_FIX\\\",\\\"HAP\\\"));\"`;
               if [ \$RES == 1 ] || [ \$PARENTACCorNAME == \"parent_acc\" ]; then echo \$line >> $alt_scaffold_file_fixed ; else echo \$ACC >> $alt_scaffold_file_seqs_to_remove; fi;
             done";

  run_command($cmd,"Removing patches on patches from file $alt_scaffold_file...");

  #$cmd = "grep -v ^# $alt_scaffold_file | awk '{if (\$2 != \"Primary\") print \$4}' | cat > $alt_scaffold_file_seqs_to_remove";
  #run_command($cmd,"Getting the list of patch sequences to remove from the fasta file $alt_scaffold_fna_file...");

  run_command("fastaremove $alt_scaffold_fna_file $alt_scaffold_file_seqs_to_remove > $alt_scaffold_fna_file_fixed","Removing patch sequences from the fasta file $alt_scaffold_fna_file...");

  # if there was nothing to remove 
  if ((!(-e $alt_scaffold_fna_file_fixed)) or # not exists
      (-e $alt_scaffold_fna_file_fixed and -z $alt_scaffold_fna_file_fixed) # exists and empty 
     ) {
  	run_command("cat $alt_scaffold_fna_file > $alt_scaffold_fna_file_fixed","There was not any patch on a patch to remove. Using same fna file as fna fixed...");
  }

  my $total_num_patches = run_command("grep -vc ^# $alt_scaffold_file",
                                      "Getting the total number of patches...",
                                      run_command("grep -c '^>' $alt_scaffold_fna_file","Getting the total number of patches..."));
  print ("Total number of patches: $total_num_patches\n");
  
  $total_num_patches = run_command("grep -vc ^# $alt_scaffold_file_fixed",
                                   "Getting the total number of patches after removing the patches on patches...",
                                   run_command("grep -c '^>' $alt_scaffold_fna_file_fixed","Getting the total number of patches after removing the patches on patches..."));
  print ("Total number of patches after removing patches on patches: $total_num_patches\n");
  
  return ($alt_scaffold_file_fixed,$alt_scaffold_fna_file_fixed);
}

sub load_patches() {
# it creates the sql file that will be used to load the patches into the specified database
# 'ftp_dir' from NCBI will be used to extract the assembly name and the assembly GCA accession 

  my ($host,$port,$name,$user,$pass,$output_dir,$ftp_dir,$alt_scaffold_filepath) = @_;

#Find the GCA accession:
#GCA_000001405.20
  my $patches_filepath = $output_dir."/patches.sql";

  my $ass_acc = (split('_',(split('/',$ftp_dir))[9]))[0]."_".
                (split('_',(split('/',$ftp_dir))[9]))[1];
  my $ass_name = (split('_',(split('/',$ftp_dir))[9]))[2];

  run_command('perl $ENSEMBL_ANALYSIS/scripts/assembly_patches/assembly_exception_load.pl '."-user $user -pass $pass -host $host -port $port -dbname $name -mapping $output_dir/alt.scaf.agp -txt $alt_scaffold_filepath -patchtype $output_dir/patch_type -assembly_name $ass_name -assembly_acc $ass_acc -coord_system scaffold -sql_file $patches_filepath","Creating the file $patches_filepath that will be used to load the patches into the database...");

  run_command("mysql -h $host -P $port -u $user -p$pass $name < $patches_filepath","Loading the patches into the database $name:$host");
}

sub check_patches() {
# Run checks to ensure that the patch loading was completed successfully.

  my ($output_path,
      $host,$port,$name,$user,$pass,
      $alt_scaffold_fixed_filepath,
      $alt_scaffold_fna_fixed_filepath,
      $cs_version) = @_;

  compare_db_with_downloaded($output_path."/compare_db_with_downloaded/",$host,$port,$name,$user,$pass,$alt_scaffold_fixed_filepath,$alt_scaffold_fna_fixed_filepath);
  check_asm_exception($output_path,$host,$port,$name,$user,$pass,$cs_version,$alt_scaffold_fixed_filepath);


}

sub compare_db_with_downloaded() {
# This will run a script to write the patch sequences to a specified directory
# where it will write assorted files while doing comparisons between the patches
# in the db and the sequences from the GRC.

  my ($output_path,$host,$port,$name,$user,$pass,$alt_scaffold_fixed_filepath,$alt_scaffold_fna_fixed_filepath) = @_;

  # create output dir if it does not exist
  if (not -e $output_path) {
    run_command("mkdir -p ".$output_path,"Create output path.");
  }

  run_command('perl $ENSEMBL_ANALYSIS/scripts/assembly_patches/compare_db_patch_seq_with_downloaded.pl '."-dbport $port -dbname $name -dbhost $host -dbport $port -dbuser $user -dbpass $pass -write_dir $output_path -alt_scaf $alt_scaffold_fixed_filepath -downloaded_fasta $alt_scaffold_fna_fixed_filepath","Comparing database patch sequences with downloaded sequences...");

  run_command("wc -l $output_path/differences.out","Checking that the comparison output file does not contain any difference...",0);

}

sub check_asm_exception() {
# This compares the sections of the sequence at the beginning and end of the assembly exception (patch)
# with the reference equivalent.

  my ($output_path,$host,$port,$name,$user,$pass,$cs_version,$alt_scaffold_fixed_filepath) = @_;

## THESE SHOULD TAKE INTO ACCOUNT THE b ORIENTATION IN THE ALT_SCAFFOLD_PLACEMENT.TXT FILE SO
## IF THERE IS A b WE LOOK INTO THE gff FILE OF THE ALIGNMENT AND CHECK THAT THE
## STRAND FOR THE ALIGNMENT AT THE END OR START IS REVERSED WHEN THERE IS A MISMATCH OF starts have different
## seqs or ends have different seqs. This sould be done in the check_asm_exception.pl script so it does not
## throw the WARNING.

#  my $output_file = $output_path."/checker_exceptions_plus_warn.txt";
#
#  # run on current core database
#  run_command('perl $ENSEMBL_ANALYSIS/scripts/assembly_patches/check_asm_exception.pl '."-host $host -port $port -user ensro -dbname $name -file test_all -all -coord_system_version $cs_version >& $output_file","Checking boundaries of the assembly exception sequences in the current database...");
#
#  # check that the number of sequences on patches that have different starts and different ends matches
#  # the number of non-zero alt_start_tail and alt_stop_tail in the alt_scaffold_placement.txt file
#  run_command("grep -v REF $output_file | grep -A9 'type PATCH' | grep -c WARN",
#              "Checking that the number of sequences on patches that have different starts and different ends matches the number of non-zero alt_start_tail and alt_stop_tail in the $alt_scaffold_fixed_filepath file",
#              run_command("grep -v '^#' $alt_scaffold_fixed_filepath | awk '{print \$(NF-1)\",\"\$NF}' | grep -v '0,0' | wc -l","Getting the number of non-zero alt_start_tail and alt_stop_tail in the $alt_scaffold_fixed_filepath file..."));

}

sub write_output {
  my $self = shift;

  return 1;
}

1;
