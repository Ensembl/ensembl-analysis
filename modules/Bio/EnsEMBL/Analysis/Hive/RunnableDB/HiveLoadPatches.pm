# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

-output_path    Output path where output and log files will be written.

-dbhost         database host name

-dbport         database port (default 3306)

-dbname         database name

-dbuser         database username to connect as

-dbpass         database password to use

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPatches -output_path $SCR9/ -dbhost genebuildX -dbname core_97 -dbuser *** -dbpass *** -dbport 3306

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPatches;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use File::Basename;
use File::Find;
use List::Util qw(sum);

sub param_defaults {
    return {
      output_path => undef,
      dbhost => undef,
      dbname => undef,
      dbuser => undef,
      dbpass => undef,
      dbport => 3306,
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;

  $self->param_required('output_path');
  $self->param_required('dbhost');
  $self->param_required('dbname');
  $self->param_required('dbuser');
  $self->param_required('dbpass');

  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('output_path') =~ /\/$/)) {
    $self->param('output_path',$self->param('output_path')."/");
  }

  # create output dir if it does not exist
  if (not -e $self->param('output_path')) {
    run_command("mkdir -p ".$self->param('output_path'),"Create output path.",0);
  }

  download_patches();
  remove_old_rawcomputes_karyobands(); 
  remove_deprecated_patches();
  recalculate_markers_map_weights();
  remove_patches_on_patches();
  load_patches();
  check_patches();
  
  
sub download_patches() {

  my ($self,$ftp_path,$local_dir) = @_;

  my $wget_verbose = "-nv";

  my @dirs = split(/\//,$ftp_path); # get array of names of subdirectories
  my @cleaned_dirs = map {$_ ? $_ : ()} @dirs; # delete empty elements in the array
  my $numDirs = @cleaned_dirs; # count the number of dirs we need to skip to be able to put the files in the local path specified as parameter

  my $cmd = "wget --no-proxy ".$wget_verbose." -r -nH --cut-dirs=".$numDirs." --reject *.rm.out.gz -P ".$local_dir." ".$ftp_path;
  my $return = system($cmd);
  if($return) {
    $self->throw("Could not download the AGP, FASTA and info files. Commandline used:\n".$cmd);
  }

  # check if no file was downloaded
  my $num_files = int(`find $local_dir -type f | wc -l`);
  if ($num_files == 0) {
    $self->throw("No file was downloaded from ".$ftp_path." to ".$local_dir.". Please, check that both paths are valid.");
  } else {
    say "$num_files files were downloaded";
  }

  # get the name of the directory where the assembly report file is
  my @ftp_path_array = split('/',$ftp_path);
  my $ass_report_dir = pop(@ftp_path_array); # Primary assembly dir
  $ass_report_dir = pop(@ftp_path_array);
  $ass_report_dir = pop(@ftp_path_array);
  $ass_report_dir = pop(@ftp_path_array);

  my $link  = $ftp_path."/../../../".$ass_report_dir."_assembly_report.txt";

  if (system("wget ".$wget_verbose." -nH -P ".$local_dir." ".$ftp_path."/../../../".$ass_report_dir."_assembly_report.txt -O ".$local_dir."/assembly_report.txt")) {
    $self->throw("Could not download *_assembly_report.txt file from ".$ftp_path."/../../../ to ".$local_dir.". Please, check that both paths are valid.");
  }
  else {
    say "Assembly report file was downloaded";
  }
  

  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.20_GRCh38.p5/GCA_000001405.20_GRCh38.p5_assembly_structure/PATCHES/alt_scaffolds/FASTA/alt.scaf.f*a.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.20_GRCh38.p5/GCA_000001405.20_GRCh38.p5_assembly_structure/PATCHES/alt_scaffolds/AGP/alt.scaf.agp.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.20_GRCh38.p5/GCA_000001405.20_GRCh38.p5_assembly_structure/PATCHES/alt_scaffolds/alt_scaffold_placement.txt
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.20_GRCh38.p5/GCA_000001405.20_GRCh38.p5_assembly_structure/PATCHES/alt_scaffolds/patch_type
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.20_GRCh38.p5/GCA_000001405.20_GRCh38.p5_assembly_report.txt

farm3-head3 /nfs/ensembl/cgg/human_assembly_patches/patch_release_5_sep15> ll
total 5576
-rw-r--r-- 1 cgg ensembl    6123 Sep 28 16:33 alt.scaf.agp.gz
-r--r--r-- 1 cgg ensembl 5613973 Sep 25 17:40 alt.scaf.fna.gz
-rw-r--r-- 1 cgg ensembl    7487 Sep 28 16:33 alt_scaffold_placement.txt
-rw-r--r-- 1 cgg ensembl   54187 Sep 28 16:33 GCA_000001405.20_GRCh38.p5_assembly_report.txt
-rw-r--r-- 1 cgg ensembl    1897 Sep 28 16:33 patch_type

gunzip *.gz

Get the total number of patches for your records and further checks:

farm3-head3 /nfs/ensembl/cgg/human_assembly_patches/patch_release_5_sep15> grep -v '^#' patch_type | wc -l
62
farm3-head3 /nfs/ensembl/cgg/human_assembly_patches/patch_release_5_sep15> grep -v '^#' alt_scaffold_placement.txt | wc -l
62
NUM_PATCHES
}

sub remove_old_rawcomputes_karyobands() {
    
}

sub remove_deprecated_patches() {
    
}

sub recalculate_markers_map_weights() {
    
}

sub remove_patches_on_patches() {
    
}

sub load_patches() {
    
}

sub check_patches() {
    
}

  download_patches($self->param('dbhost'),
                                 $self->param('dbport'),
                                 $self->param('dbuser'),
                                 $self->param('dbpass'),
                                 $self->param('dbname')) if ($self->param('skip') < 1 or $self->param('only') == 1);

  methionine_to_stop_codon($self->param('dbhost'),
                           $self->param('dbport'),
                           $self->param('dbuser'),
                           $self->param('dbpass'),
                           $self->param('dbname'),
                           $self->param('dnadbhost'),
                           $self->param('dnadbport'),
                           $self->param('dnadbname'),
                           $self->param('check_vega_met_stop_dir'),
                           $self->param('output_path')) if ($self->param('skip') < 2 or $self->param('only') == 2);

  set_ncrna_host_gene_attribute($self->param('dbhost'),
                                $self->param('dbport'),
                                $self->param('dbuser'),
                                $self->param('dbpass'),
                                $self->param('dbname')) if ($self->param('skip') < 3 or $self->param('only') == 3);

  truncate_tsf_table($self->param('dbhost'),
                     $self->param('dbport'),
                     $self->param('dbuser'),
                     $self->param('dbpass'),
                     $self->param('dbname')) if ($self->param('skip') < 4 or $self->param('only') == 4);

  add_attribute_to_GAGE_cluster($self->param('dbhost'),
                                $self->param('dbport'),
                                $self->param('dbuser'),
                                $self->param('dbpass'),
                                $self->param('dbname')) if ($self->param('skip') < 5 or $self->param('only') == 5);
  return 1;
}



Remove old patch raw compute and karyotype bands 

mysql -h genebuild12 -u ensadmin -pensembl carlos_homo_sapiens_core_83 < $ENSEMBL_ANALYSIS/scripts/assembly_patches/remove_patch_raw_compute.sql
mysql -h genebuild12 -u ensadmin -pensembl carlos_homo_sapiens_core_83 < $ENSEMBL_ANALYSIS/scripts/assembly_patches/remove_patch_karyotype.sql
 
Remove altered and deprecated patches 

mv ~/enscode/ensembl-personal/cgg/env/GRCh38.p4_e83.sh ~/enscode/ensembl-personal/cgg/env/GRCh38.p5_e83.sh

vi $ENSEMBL_PERSONAL/cgg/env/GRCh38.p5_e83.sh
  #export BASE="/nfs/ensembl/cgg/homo_sapiens/e79"
  export BASE="/nfs/ensembl/cgg/human_assembly_patches/patch_release_5_sep15"

source $ENSEMBL_PERSONAL/cgg/env/GRCh38.p5_e83.sh

cd $BASE
perl $ENSEMBL_ANALYSIS/scripts/assembly_patches/remove_updated_and_deprecated_patches.pl -user ensro -host genebuild12 -dbname carlos_homo_sapiens_core_83 -patchtype_file $BASE/patch_type -central_cs scaffold > $SCR9/alt_patch_output.txt

Get the list of altered patches for your records:

 grep modified $SCR9/alt_patch_output.txt
 none

Get the list of deprecated patches for your records:

 grep removed $SCR9/alt_patch_output.txt
 none

Get the list of new patches for your records:

farm3-head3 /nfs/ensembl/cgg/human_assembly_patches/patch_release_5_sep15> grep 'is a new patch' $SCR9/alt_patch_output.txt
HSCHR1_4_CTG3 is a new patch
HSCHR1_3_CTG3 is a new patch
HSCHR1_5_CTG32_1 is a new patch
HSCHR4_2_CTG4 is a new patch
HSCHR4_8_CTG12 is a new patch
HSCHR4_9_CTG12 is a new patch
HSCHR6_1_CTG10 is a new patch
HG2072_PATCH is a new patch
HSCHR9_1_CTG6 is a new patch
HSCHR9_1_CTG7 is a new patch
HSCHR10_1_CTG6 is a new patch
HG2334_PATCH is a new patch
HG2116_PATCH is a new patch
HSCHR12_2_CTG1 is a new patch
HSCHR13_1_CTG7 is a new patch
HSCHR13_1_CTG8 is a new patch
HG2139_PATCH is a new patch
HSCHR16_5_CTG1 is a new patch
HSCHR16_4_CTG3_1 is a new patch
HSCHR18_5_CTG1_1 is a new patch
HG2213_PATCH is a new patch
HG26_PATCH is a new patch
HSCHR22_6_CTG1 is a new patch
HSCHR22_7_CTG1 is a new patch

Delete the altered and deprecated patches:

cat $BASE/delete_patch.sql
none

#mysql -h genebuild7 -u ensadmin -p*** CORE_DB < $BASE/delete_patch.sql
 
cd $ENSCODE/ensj-healthcheck

I have to change branch because the latest code does not work.

git checkout task/merge/human/GRCh38_e79

> vi database.properties
5 host1=genebuild12
  6 port1=3306
  7 user1=ensro
  8 password1=
  9 driver1=org.gjt.mm.mysql.Driver

bsub -q normal -M2000 -R 'select[mem>2000] rusage[mem=2000]' -o $SCR9/hc_before_loading_patch.log "./run-healthcheck.sh -d carlos_homo_sapiens_core_83 -output info -species human -type core post_genebuild >& $SCR9/hc_before_loading_patch.txt"

#git checkout task/merge/human/GRCh38_e81

less $SCR9/hc_before_loading_patch.txt

The following FAILED healthcheck needs to be fixed:

 org.ensembl.healthcheck.testcase.generic.MarkerFeatures

Common problems and fixes: 

Since markers were added to patches that went through a full build or a markers update in the past, if the patches containing markers have been deleted, the weights may be wrong for remaining markers and the MarkerFeatures healthcheck will fail. To fix that, we need to recalculate the markers map weights by running the following script:

#bsub -M900 -R 'select[mem>900] rusage[mem=900]' -o /lustre/scratch110/sanger/cgg/homo_sapiens79/marker.out -e /lustre/scratch110/sanger/cgg/homo_sapiens79/marker.err perl $ENSEMBL_ANALYSIS/scripts/markers/unmapped_markers.pl -host genebuild7 -user ensadmin -pass ensembl -dbname carlos_homo_sapiens_core_79_38 -port 3306 -logic_name marker -max_duplicates 3

This job can take up to 2 hours to run.

Check that it was successfully completed:

#less /lustre/scratch110/sanger/cgg/homo_sapiens79/marker.out

Check that the number of markers is similar to the one got for previous release:

# cat /lustre/scratch110/sanger/cgg/homo_sapiens79/marker.err
Found 296045 markers, checking / setting marker_feature map weights....
Found 0 duplicates and 0 unmapped markers
Creating unmapped entries...
Deleting duplicated and unmapped markers...
Done

#For example, for human, NUMBER_OF_MARKERS ~ 299818. 

#bsub -q normal -M2000 -R 'select[mem>2000] rusage[mem=2000]' -o /lustre/scratch110/sanger/cgg/homo_sapiens79/hc_before_loading_patch2.log "./run-healthcheck.sh -d carlos_homo_sapiens_core_79_38 -output info -species human -type core post_genebuild >& /lustre/scratch110/sanger/cgg/homo_sapiens79/hc_before_loading_patch2.txt"

carlos_homo_sapiens_core_79_38:  No marker features have priorities greater than the threshold (50)

#Again but that should be ok.

I didn't run the markers recalculation because I know there has not been any patch removed.


Find the GCA accession:
GCA_000001405.20

Generate the 'sql.txt' file that you will use to load patches (use bsub if you have a huge amount of patch):


There's a patch on a patch which we cannot process at the moment so I'll comment its line.

vi $BASE/alt_scaffold_placement.txt
#PATCHES  ALT_REF_LOCI_2    HG2139_PATCH  KN538374.1  SCAFFOLD  HSCHR15_4_CTG8  KI270905.1  REGION1A  + 1 4998962 1 5161414 0 0

perl $ENSEMBL_ANALYSIS/scripts/assembly_patches/assembly_exception_load.pl -user ensro -host genebuild12 -dbname carlos_homo_sapiens_core_83 -mapping $BASE/alt.scaf.agp -txt $BASE/alt_scaffold_placement.txt -patchtype $BASE/patch_type -assembly_name GRCh38.p5 -assembly_acc GCA_000001405.20 -coord_system scaffold -sql_file $BASE/sql.txt

Load the patches into the CORE_DB:

mysql -h genebuild12 -u ensadmin -pensembl carlos_homo_sapiens_core_83 < $BASE/sql.txt

compare_db_patch_seq_with_downloaded.pl 

This script will write the patch sequences to a specified directory, where it will write assorted files while doing comparisons between the patches in the db and the sequences from the GRC.

mkdir -p $SCR9/checks/download_seq/write
cp $BASE/alt_scaffold_placement.txt $SCR9/checks/download_seq/alt_scaffold_placement.txt
cp $BASE/alt.scaf.fna $SCR9/checks/download_seq/alt.scaf.fna

cd $SCR9/checks/download_seq

vi alt.scaf.fna
(replace header of KN538374.1 with ">DELETE")

vi to_remove.txt
DELETE 

fastaremove alt.scaf.fna to_remove.txt > alt.scaf.fna.fixed

grep -c ">" alt.scaf.fna.fixed
61
grep -c ">" alt.scaf.fna
62

[ensadmin@genebuild12:3306:carlos_homo_sapiens_core_83]> select count(*) from assembly_exception where exc_type='HAP' or exc_type like '%PATCH%';
+----------+
| count(*) |
+----------+
|      322 |
+----------+
1 row in set (0.00 sec)

[ensadmin@genebuild12:3306:carlos_homo_sapiens_core_83]> select count(*) from assembly_exception where exc_type like '%PATCH%';
+----------+
| count(*) |
+----------+
|       61 |
+----------+
1 row in set (0.00 sec)

perl $ENSEMBL_PERSONAL/genebuilders/sanity_scripts/compare_db_patch_seq_with_downloaded.pl -dbport 3306 -dbname carlos_homo_sapiens_core_83 -dbhost genebuild12 -dbuser ensro -write_dir write -alt_scaf alt_scaffold_placement.txt -download alt.scaf.fna.fixed

The last line in the output will contain the path to the file listing the differences (if any). Check the file, if there are differences, there is a problem:

wc -l write/differences.out

0 write/differences.out

    nt_transform_check.pl 

This script works around the three levels of the assembly checking that it is possible to map the whole way around the triangle. The script deals with the branching that occurs when the patches are in place (i.e. two mappings at one level of the assembly).

for CHR in `mysql -h genebuild12 -P 3306 -u ensro -Dcarlos_homo_sapiens_core_83 -e'select name from seq_region where coord_system_id=4 AND name like "CHR%"' -NB`
do
bsub -q long -M 1000 -R 'select[mem>1000] rusage[mem=1000]' -o $SCR9/checks/nt_transform_check_multi_$CHR.log "perl $ENSEMBL_PERSONAL/genebuilders/sanity_scripts/nt_transform_check.pl -h genebuild12 -u ensro -dbname carlos_homo_sapiens_core_83 -port 3306 -central_coordsystem scaffold -path GRCh38 -chromosome_name $CHR > $SCR9/checks/nt_transform_check_multi_$CHR.out"
done

1st Oct 2015

grep "Successfully completed" $SCR9/checks/nt_transform_check_*.log | wc -l
322

ls $SCR9/checks/nt_transform_check_*.log | wc -l
322

farm3-head2 /lustre/scratch109/sanger/cgg/homo_sapiens83> bsub -I "grep PROBLEM $SCR9/checks/nt_transform_check*.out"
Job <5052236> is submitted to default queue <normal>.
<<Waiting for dispatch ...>>
<<Starting on bc-31-2-05>>
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HG142_HG150_NOVEL_TEST.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HG151_NOVEL_TEST.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HG2116_PATCH.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HG2217_PATCH.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG1_1.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG1_2.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG2.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG3.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG5.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG6.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG6.out:PROBLEM: contig::AP006285.2:56928:190267:1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG7.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG8.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_2_CTG1_1.out:PROBLEM: contig::AC139749.4:68834:128565:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR11_3_CTG1.out:PROBLEM: contig::AC139749.4:12502:118739:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR17_5_CTG4.out:PROBLEM: contig::AC006070.1:80668:145532:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR17_6_CTG4.out:PROBLEM: contig::AC006070.1:99441:145532:-1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR19KIR_ABC08_AB_HAP_C_P_CTG3_1.out:PROBLEM: contig::AC245128.3:122720:213609:1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG2.out:PROBLEM: contig::AL022318.2:126968:127636:1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG5.out:PROBLEM: contig::Z82185.1:14244:35506:1 - this shouldn't happen
/lustre/scratch109/sanger/cgg/homo_sapiens83/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG6.out:PROBLEM: contig::AP000344.1:138998:184051:1 - this shouldn't happen

Same as in previous release. This is ok.

#All these are HAPs which were already checked when the assembly was loaded so are fine. I'll check the new one on a patch (CHR_HG2217_PATCH).

#From the first time the assembly was loaded:
"This has to do with the last grep ":1" below since they are all reports about contigs with no matches. It doesn't worry me since these contigs will likely match something on other regions. I'll check that below."

#farm3-head1 ~> bsub -I "grep 'MAPPED ONLY' /lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check*.out"
Job <2739893> is submitted to default queue <normal>.
<<Waiting for dispatch ...>>
<<Starting on bc-24-1-03>>

#ok

#farm3-head3 ~> bsub -I "grep 'contigs with no matches' /lustre/scratch110/sanger/cgg/homo_sapiens81/checks/nt_transform_check*.out | grep -v '^0'"
Job <2741164> is submitted to default queue <normal>.
<<Waiting for dispatch ...>>
<<Starting on bc-19-4-13>>
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG126_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG1362_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG142_HG150_NOVEL_TEST.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG151_NOVEL_TEST.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG1832_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2021_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2022_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2030_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2058_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2062_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2066_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2095_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2104_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2128_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2191_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2216_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2217_PATCH.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2232_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2233_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2241_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2242_HG2243_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2244_HG2245_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2247_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2249_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2288_HG2289_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2291_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG23_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG986_PATCH.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR10_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR10_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR10_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR10_1_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG1_1.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG1_2.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG2.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG3.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG5.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG6.out:2 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG7.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_1_CTG8.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_2_CTG1_1.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_2_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR11_3_CTG1.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_1_CTG11.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_1_CTG31.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_1_CTG32_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_1_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_2_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_2_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_3_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_3_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_4_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_4_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_5_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_5_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_6_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR12_7_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_2_CTG31.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_2_CTG32_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_2_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR13_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR13_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR13_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR13_1_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR13_1_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR13_1_CTG6.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_3_CTG31.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_3_CTG32_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR14_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR14_2_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR14_3_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR14_7_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_4_CTG31.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_4_CTG32_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_1_CTG8.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_2_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_2_CTG8.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_3_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_3_CTG8.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_4_CTG8.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR15_5_CTG8.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR16_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR16_1_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR16_2_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR16_3_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR16_4_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR16_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_10_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_1_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_1_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_1_CTG9.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_2_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_2_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_2_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_2_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_3_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_3_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_4_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_5_CTG4.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_6_CTG4.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_7_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_8_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR17_9_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_1_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_1_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_2_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_2_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_2_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_3_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_4_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_ALT21_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR18_ALT2_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_1_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_2_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_2_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_3_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_3_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_4_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_4_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19_5_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_ABC08_A1_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_ABC08_AB_HAP_C_P_CTG3_1.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_ABC08_AB_HAP_T_P_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH05_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH05_B_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH06_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH06_BA1_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH08_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH08_BAX_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH13_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH13_BA2_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH15_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_FH15_B_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_G085_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_G085_BA1_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_G248_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_G248_BA2_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_GRC212_AB_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_GRC212_BA1_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_LUCE_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_LUCE_BDEL_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_RP5_B_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_RSH_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_RSH_BA2_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_T7526_A_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19KIR_T7526_BDEL_HAP_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19LRC_COX1_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19LRC_COX2_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19LRC_LRC_I_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19LRC_LRC_J_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19LRC_LRC_S_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19LRC_LRC_T_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19LRC_PGF1_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR19LRC_PGF2_CTG3_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR1_ALT2_1_CTG32_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR20_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR20_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR20_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR20_1_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR21_1_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR21_2_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR21_3_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR21_4_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR21_5_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR21_6_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR21_8_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_1_CTG15.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_1_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_1_CTG7_2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_1_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG2.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG5.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG6.out:1 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_1_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_2_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_3_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_4_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR22_5_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_2_CTG15.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_2_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_2_CTG7_2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_2_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_3_CTG15.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_3_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_3_CTG7_2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_4_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_4_CTG7_2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR2_5_CTG7_2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_1_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_2_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_2_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_3_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_3_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_3_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_4_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_4_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_5_CTG2_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_5_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_6_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_7_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_8_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR3_9_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_1_CTG12.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_1_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_1_CTG6.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_1_CTG8_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_1_CTG9.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_2_CTG12.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_3_CTG12.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_4_CTG12.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_5_CTG12.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_6_CTG12.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR4_7_CTG12.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_1_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_1_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_2_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_2_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_2_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_3_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_3_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_3_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_4_CTG1_1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_4_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_5_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_6_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR5_7_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_1_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_1_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_1_CTG6.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_1_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_1_CTG8.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_1_CTG9.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_8_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_MHC_APD_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_MHC_COX_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_MHC_DBB_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_MHC_MANN_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_MHC_MCF_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_MHC_QBL_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR6_MHC_SSTO_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_1_CTG4_4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_1_CTG6.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_1_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_2_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_2_CTG4_4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_2_CTG6.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_2_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR7_3_CTG6.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_1_CTG6.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_1_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_2_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_2_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_3_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_3_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_4_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_4_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_5_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_5_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_6_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_6_CTG7.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_7_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_8_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR8_9_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR9_1_CTG1.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR9_1_CTG2.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR9_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR9_1_CTG4.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHR9_1_CTG5.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHRX_1_CTG3.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHRX_2_CTG12.out:0 contigs with no matches.
/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HSCHRX_2_CTG3.out:0 contigs with no matches.

I'll check the one on the patch as the rest were checked when the assembly was loaded and they were used somewhere else.

/lustre/scratch110/sanger/cgg/homo_sapiens79/checks/nt_transform_check_multi_CHR_HG2217_PATCH.out
AC139749.4

It turns out that this contig was already checked for e76 and it mapped somewhere else so it is ok.

#
#Check that these contigs have matches somewhere else:
#AC139749.4
#AP006285.2
#AC006070.1
#AC245128.3
#AL022318.2
#Z82185.1
#AP000344.1

#cgg@farm3-head2:cgg/homo_sapiens76> egrep 'AC139749.4|AP006285.2|AC006070.1|AC245128.3|AL022318.2|Z82185.1|AP000344.1' $SCR9/GRCh38_load_assembly/Primary_Assembly/*all.agp
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000101.2     952544  954933  16      F       AC139749.4      129287  131676  -
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000101.2     955225  955392  18      F       AC139749.4      128828  128995  -
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000101.2     955655  1015386 20      F       AC139749.4      68834   128565  -
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000101.2     1045972 1077810 22      F       AC139749.4      12502   44340   -
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000101.2     1078209 1090311 24      F       AC139749.4      1       12103   -
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000101.2     1463474 1634855 39      F       AP006285.2      18886   190267  +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000132.2     14225598        14371129        174     F       AC006070.1      1       145532  -
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000140.2     27457503        27671111        696     F       AC245128.3      1       213609  +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     4644479 4752444 148     O       AP000344.1      30916   138881  +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     4752561 4769880 150     O       AP000344.1      138998  156317  +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     18103763        18117593        564     F       Z82185.1        105     13935   +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     18117903        18137745        566     F       Z82185.1        14244   34086   +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     18138117        18139166        568     F       Z82185.1        34457   35506   +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     20174730        20200679        644     F       AL022318.2      101     26050   +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     20200821        20223461        646     F       AL022318.2      26192   48832   +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     20223749        20302265        648     F       AL022318.2      49120   127636  +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     20302556        20358665        650     F       AL022318.2      127927  184036  +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/placed_scaffolds_all.agp:GL000155.2     20358969        20375772        652     F       AL022318.2      184340  201143  +

#cgg@farm3-head2:cgg/homo_sapiens76> egrep 'GL000101.2|GL000132.2|GL000140.2|GL000155.2' $SCR9/GRCh38_load_assembly/Primary_Assembly/assembled_chromosomes/AGP/*.agp
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/assembled_chromosomes/AGP/chr11.agp:CM000673.2 60001   50821348        3       F       GL000101.2      1       50761348        +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/assembled_chromosomes/AGP/chr17.agp:CM000679.2 26935981        81742542        19      F       GL000132.2      1       54806562 +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/assembled_chromosomes/AGP/chr19.agp:CM000681.2 27240875        58607616        15      F       GL000140.2      1       31366742 +
/lustre/scratch109/sanger/cgg/homo_sapiens76/GRCh38_load_assembly/Primary_Assembly/assembled_chromosomes/AGP/chr22.agp:CM000684.2 18709565        49973865        71      F       GL000155.2      1       31264301 +

#Ok. All of them are there so they map to other places. Good.


    check_asm_exception.pl 

This compares the sections of the sequence at the beginning and end of the assembly exception (patch) with the reference equivalent. To cut down on the work, you should run this script against your database and a previously OK'd database. Then you can identify the new cases and only check them.

Run on previous core database on staging:

perl $ENSEMBL_PERSONAL/genebuilders/sanity_scripts/check_asm_exception.pl -host ens-staging1 -port 3306 -user ensro -dbname homo_sapiens_core_82_38 -file test_all -all -coord_system_version GRCh38 >& /lustre/scratch109/sanger/cgg/homo_sapiens83/checks/checker_exceptions_plus_warn_prev.txt

Run on current core database:

perl $ENSEMBL_PERSONAL/genebuilders/sanity_scripts/check_asm_exception.pl -host genebuild12 -port 3306 -user ensro -dbname carlos_homo_sapiens_core_83 -file test_all -all -coord_system_version GRCh38 >& /lustre/scratch109/sanger/cgg/homo_sapiens83/checks/checker_exceptions_plus_warn.txt

Check the jobs finished successfully:
less /lustre/scratch109/sanger/cgg/homo_sapiens83/checks/checker_exceptions_plus_warn_prev.txt
less /lustre/scratch109/sanger/cgg/homo_sapiens83/checks/checker_exceptions_plus_warn.txt

Compare the differences:

grep WARN /lustre/scratch109/sanger/cgg/homo_sapiens81/checks/checker_exceptions_plus_warn.txt | sort > /lustre/scratch109/sanger/cgg/homo_sapiens83/checks/checker_exceptions_plus_warn.grep
grep WARN /lustre/scratch109/sanger/cgg/homo_sapiens81/checks/checker_exceptions_plus_warn_prev.txt | sort > /lustre/scratch109/sanger/cgg/homo_sapiens83/checks/checker_exceptions_plus_warn_prev.grep

Get the new warnings

diff /lustre/scratch109/sanger/cgg/homo_sapiens83/checks/checker_exceptions_plus_warn.grep /lustre/scratch109/sanger/cgg/homo_sapiens83/checks/checker_exceptions_plus_warn_prev.grep

Go through the file

    Your patch has mismatched sequences at one end 

    Check the alt_start_tail and alt_stop_tail for the patch (and the strand) in alt_scaffold_placement.txt

        It's ok if the corresponding value is not 0
        Otherwise there might be a problem 

    Your patch has mismatched sequences at both ends, this should NOT happen 

    Check each of the sequence end, hopefully only one is problematic

        You might need to report to GRC 

ok!


Get the following numbers:

[ensadmin@genebuild12:3306:carlos_homo_sapiens_core_83]> select count(*) from seq_region;
+----------+
| count(*) |
+----------+
|   267729 |
+----------+
1 row in set (0.00 sec)

[ensadmin@genebuild12:3306:carlos_homo_sapiens_core_83]> select count(*) from assembly ;
+----------+
| count(*) |
+----------+
|  1186210 |
+----------+
1 row in set (0.00 sec)

[ensadmin@genebuild12:3306:carlos_homo_sapiens_core_83]> select count(*) from assembly_exception;
+----------+
| count(*) |
+----------+
|      324 |
+----------+
1 row in set (0.00 sec)

Update meta_coord:

perl $ENSEMBL/misc-scripts/meta_coord/update_meta_coord.pl -dbpattern ^carlos_homo_sapiens_core_83$ -host genebuild12 -port 3306 -user ensadmin -pass ensembl
Original meta_coord table backed up in carlos_homo_sapiens_core_83_1.meta_coord.backup
Updating assembly_exception table entries... done
Updating density_feature table entries... done
Updating ditag_feature table entries... done
Updating dna_align_feature table entries... done
Updating exon table entries... done
Updating gene table entries... done
Updating intron_supporting_evidence table entries... done
Updating karyotype table entries... done
Updating marker_feature table entries... done
Updating misc_feature table entries... done
Updating prediction_exon table entries... done
Updating prediction_transcript table entries... done
Updating protein_align_feature table entries... done
Updating repeat_feature table entries... done
Updating simple_feature table entries... done
Updating transcript table entries... done
==> Done with carlos_homo_sapiens_core_83/1
==> All done.


sub write_output {
  my $self = shift;

  return 1;
}

sub readthrough_transcripts_tagged {
  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;

  my $readthrough_transcripts_sql = "update transcript_attrib ta, transcript t set ta.attrib_type_id = (select attrib_type_id from attrib_type where code = 'readthrough_tra') where (ta.value like '%read%through%') and ta.transcript_id = t.transcript_id;";

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$readthrough_transcripts_sql\"",
              "\n===== STEP 1: Readthrough transcripts tagged =====\n");

  my $num_readthrough_transcripts_sql = "select count(*) from transcript_attrib ta,transcript t, attrib_type at where code = 'readthrough_tra' and (ta.value like '%read%through%') and ta.transcript_id = t.transcript_id;";

  my $num_updated = run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$num_readthrough_transcripts_sql\"");

  print("\nThere are ".int($num_updated)." readthrough transcripts tagged now.\n");
}

sub methionine_to_stop_codon {
  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname,$dnadbhost,$dnadbport,$dnadbname,$check_vega_met_stop_dir,$output_path) = @_;

  my $ids_file = "$output_path/havana_coding_transcript_ids.txt";
  my $full_ids_file = "$output_path/full_length_havana_coding_transcript_ids.txt";

  my $dump_ids_sql = "select distinct(tr.transcript_id) from translation tr, transcript t, transcript_attrib ta, gene g where t.gene_id = g.gene_id and t.transcript_id = ta.transcript_id and tr.transcript_id = t.transcript_id;";
  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$dump_ids_sql\" > $ids_file",
              "\n===== STEP 2: Methionine to stop codon =====\n");

  my $num_ids = int(run_command("wc -l $ids_file | cut -f1 -d' '"));

  run_command("perl $check_vega_met_stop_dir/check_vega_met_stop.pl -dbname $dbname -host $dbhost -dna_dbname $dnadbname -dna_host $dnadbhost -dna_port $dnadbport < $ids_file > $full_ids_file");

  my $num_wrong_ids = int(
     run_command(
       "awk '\$3 != \"ok\" { nok++ } END { print nok }' $full_ids_file")
  );
  my $num_correct_ids = $num_ids - $num_wrong_ids;

  # Check that $num_wrong_ids is similar to the previous release one (if previous release exists).
  warning("Methionine to stop codon check: Please check manually that the number of wrong transcripts $num_wrong_ids is similar to the one in previous release (if exists). Some of the previous releases numbers are: e71: 29260; e70: 29381 e69: 56594; e68: 53703");
}

sub set_ncrna_host_gene_attribute {
  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;

  my $ncrna_host_sql = "update gene_attrib ga, gene g set ga.attrib_type_id = (select attrib_type_id from attrib_type where code = 'ncrna_host') where ga.gene_id = g.gene_id and (value like 'transcribed%' or value = 'ncrna_host');";
  my $check_sql_1 = "select count(*) from gene_attrib where (value like 'transcribed%' or value = 'ncrna_host');";
  my $check_sql_2 = "select count(*) from gene g, gene_attrib ga where ga.attrib_type_id = (select attrib_type_id from attrib_type where code = 'ncrna_host') and ga.gene_id = g.gene_id;";

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$ncrna_host_sql\"",
              "\n===== STEP 3: Set ncrna_host gene attribute =====\n");

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$check_sql_2\"",
              "",
              run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$check_sql_1\"",
                          "\nChecking that the number of updated attrib_type_id matches the number of values that were found.\n")
             );
}

sub truncate_tsf_table {
  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;

  my $tsf_bak_name = "tsf_bak_".time();

  my $num_tsf = run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from transcript_supporting_feature;'",
                            "\n===== STEP 4: Truncate transcript_supporting_feature table =====\n");

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'rename table transcript_supporting_feature to $tsf_bak_name; create table transcript_supporting_feature like $tsf_bak_name;'");

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from $tsf_bak_name;'","",$num_tsf);

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from transcript_supporting_feature;'","",0);

}

sub add_attribute_to_GAGE_cluster {
  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;


  my $num_att = run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"select count(*) from transcript,gene_attrib ga where value = 'gene_cluster_GAGE' and transcript.gene_id=ga.gene_id;\"",
                            "\n===== STEP 5: Add attribute to GAGE cluster =====\n");

  my $insert_sql = <<END;
insert into transcript_attrib (transcript_id,attrib_type_id,value)
select t.transcript_id,(select attrib_type_id from attrib_type where code = 'gene_cluster'),'gene_cluster_GAGE'
from transcript t, gene_attrib ga, attrib_type at
where t.gene_id = ga.gene_id and ga.attrib_type_id = at.attrib_type_id and ga.value = 'gene_cluster_GAGE';
END

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$insert_sql\"");

  my $check_sql = "select count(*) from transcript_attrib where attrib_type_id=(select attrib_type_id from attrib_type where code = 'gene_cluster');";
  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$check_sql\"",
              "",
              $num_att);
}

1;
