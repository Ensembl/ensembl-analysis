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


use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use File::Basename;
use File::Spec::Functions qw(catfile);
use List::Util qw(sum);

my $output_path;
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
my $analysis_scripts;
my $uniprot_filename;
my $no_external_db = 0;
my $no_backup = 0;
my $clean = 0;
my $help = 0;
my $verbose = 0;
my $do_daf = 1;
my $do_paf = 1;
my $do_ise = 0;
my $collapse_ise = 0;

GetOptions('output_path:s' => \$output_path, 
           'dbhost|host|h=s'      => \$dbhost,
           'dbname|db|D=s'      => \$dbname,
           'dbuser|user|u=s'      => \$dbuser,
           'dbpass|pass|p=s'      => \$dbpass,
           'dbport|port|P=s'      => \$dbport,
           'prod_dbuser=s' => \$prod_dbuser,
           'prod_dbpass=s' => \$prod_dbpass,
           'prod_dbhost=s' => \$prod_dbhost,
           'prod_dbname=s' => \$prod_dbname,
           'prod_dbport=s' => \$prod_dbport,
           'analysis_scripts=s'     => \$analysis_scripts,
           'uniprot_filename=s' => \$uniprot_filename, 
           'no_external_db!'=> \$no_external_db, 
           'no_backup!'    => \$no_backup,       
           'clean!'        => \$clean,           
           'daf!'        => \$do_daf,           
           'paf!'        => \$do_paf,           
           'ise!'        => \$do_ise,
           'core!'      => \$collapse_ise,
           'help'          => \$help,
           'verbose'       => \$verbose);
print $0, "\n";

($analysis_scripts) = $0 =~ /(.*)\/[^\/]+$/ unless $analysis_scripts ;
if (!$output_path or !$dbport or !$dbhost or !$dbname or !$dbuser or !$dbpass or !$analysis_scripts or $help)
{
    &usage;
    exit(1);
}

if( !$no_external_db )
{
    if ($do_paf) {
      if (!$uniprot_filename) {
          warn("You haven't specified if you want the external db id to be populated...\nUniprot filename: $uniprot_filename\nno_external_db: $no_external_db");
          &usage;
          exit(1);
      }
      throw( "uniprot_filename needs to be a file (typically like /data/blastdb/Ensembl/uniprot_yyyy_mm/entry_loc)" ) unless ( -f $uniprot_filename);
    }
    if (!$prod_dbhost or !$prod_dbname or !$prod_dbuser)
    {
        warn('You did not specify all parameters for the production db: '.join(' ', $prod_dbhost, $prod_dbuser, $prod_dbname, $prod_dbpass));
        &usage;
        exit(1);
    }

}


#-------
# BEGIN
#-------
$output_path =~ s/\/$//;
$analysis_scripts =~ s/\/$//;

if (system("mkdir -p $output_path")) {
  throw("Cannot create output_path $output_path.");
}

my $num_fixed = 0; # will be calculated later
my $fix_supporting_evidence_script = catfile($analysis_scripts, 'fix_supporting_evidence_links.pl');
my $test_regex_script = catfile($analysis_scripts, 'test_regexes.pl');
my $assign_db_id_script = catfile($analysis_scripts, 'assign_external_db_ids.pl');
my $ise_script = catfile($analysis_scripts, 'fix_intron_supporting_evidence_links.pl');
foreach my $script ($fix_supporting_evidence_script, $test_regex_script, $assign_db_id_script, $ise_script) {
    throw("Cannot access $script") unless (-e $script);
}


if ($verbose) {
  print("\nhost: $dbhost\n");
  print("port: $dbport\n");
  print("name: $dbname\n");
  print("uniprot file: $uniprot_filename\n") if (!$no_external_db and $do_paf) ;
}

my @files_to_delete;

# ================
# BACK UP DATABASE
# ================

my $db_backup_file = $output_path.'/'.$dbname.'_backup.sql';

if ($no_backup) {
  print("\n===== SKIPPING: BACK UP DATABASE =====\n");
} else {
  print("\n===== BACK UP DATABASE =====\n") if ($verbose);

  print("\nBacking up database (dumping daf, paf, tsf and sf tables)...\n") if ($verbose);

  if(system("mysqldump -h$dbhost -P$dbport -u$dbuser -p$dbpass --quick --add-drop-table $dbname dna_align_feature protein_align_feature transcript_supporting_feature supporting_feature > $db_backup_file")) {
    throw("The dump was not completed. Please, check that you have enough disk space in the output path ($output_path) as well as writing permission.");
  }
  else {
    print("\nThe backup dump was completed successfully into file $db_backup_file. It contains daf, paf, tsf and sf tables dumps.\n") if ($verbose);
    print("\nRemember that you will be able to undo the changes made to your database by running the following command:\n");
    print("mysql -h$dbhost -P$dbport -u$dbuser -p*** -D$dbname < $db_backup_file\n");
    push(@files_to_delete, $db_backup_file);
  }
}

my @types;
push(@types, 'dna_align_feature') if ($do_daf);
push(@types, 'protein_align_feature') if ($do_paf);

my %synonyms = (
    protein_align_feature => 'paf',
    dna_align_feature     => 'daf',
    );
my %moltype = (
    protein_align_feature => 'protein',
    dna_align_feature     => 'dna',
    );
my %human_readable = (
    protein_align_feature => 'protein align features',
    dna_align_feature     => 'dna align features',
);

foreach my $type (@types) {
    my $dumped_file;
    my $num_ori = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'SELECT count(*) FROM $type'`);
    if (!$num_ori) {
        print "The table $type is empty, SKIPPING\n";
        next;
    }
    my $num_fixed;
    if ($no_external_db) {
        print("\n===== SKIPPING: ASSIGN EXTERNAL DB IDS =====\n") if ($verbose);
        $dumped_file = $output_path."/$synonyms{$type}.dump";
        # dump daf here since it was not dumped in the skipped assign external db step
        print("\nDumping '$type' table...\n") if ($verbose);
        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -q -e'SELECT * FROM $type' > $dumped_file")) {
            throw("Could not dump $type table\n");
        }
        # check
        $num_fixed = line_count($dumped_file);
        if ($num_ori != $num_fixed) {
          throw("The number of '$type' rows ($num_ori) does not match the number of dumped features ($num_fixed).");
        } else {
          print("The number of '$type' rows ($num_ori) matches the number of dumped features ($num_fixed). Great!\n") if ($verbose);
          print("Output file: $dumped_file`\n") if ($verbose);
        }
    }
    else {
# ======================
# ASSIGN EXTERNAL DB IDS
# ======================

        print("\n===== ASSIGN EXTERNAL DB IDS =====\n") if ($verbose);
        print "\nCreating the regular expression file for your ", uc($human_readable{$type}), " to assign external db ids...\n" if ($verbose);

        my $cfg_file = $output_path.'/'.$type.'.cfg';
        if (-e $cfg_file) {
            if(system("cp $cfg_file $cfg_file.$$")) {
                throw("Could not back up $cfg_file. You might want to delete this file manually or work in a different directory");
            }
            else {
                warning("$cfg_file exists!\nBackup is $cfg_file.$$");
                if (system("rm $cfg_file")) {
                    throw("Could not delete $cfg_file");
                }
            }
        }
        push(@files_to_delete, $cfg_file);
        my $cfg_log_file = $output_path.'/'.$type.'.cfg.log';
        push(@files_to_delete, $cfg_log_file);
        if (system("perl $test_regex_script -dbname $dbname -dbhost $dbhost -dbport $dbport -dbuser $dbuser -dbpass $dbpass -type $synonyms{$type} -main_regex_file $analysis_scripts/$synonyms{$type}_regexes.dat -output_config_file $cfg_file > $cfg_log_file")) {
            throw("Could not execute $test_regex_script\n");
        }

        print("Output files: $cfg_file and $cfg_log_file\n") if ($verbose);

        my $num_logic_names = line_count("grep Analysis $cfg_file | sort | uniq | ");
        my $num_all_matched = int(`grep -c "GOOD - All acc matched a regular expression" $cfg_log_file`);

        if ($num_logic_names <= 0) {
            warning("Could not find any logic name for $human_readable{$type} analyses.");
        } elsif ($num_logic_names != $num_all_matched) {
            throw("Some acc do not match a regular expression. The number of logic names found for $human_readable{$type} ($num_logic_names) does not match the number of 'GOOD - All acc matched a regular expression' ($num_all_matched).\nPLEASE, HAVE A LOOK AT $cfg_file and $cfg_log_file\n");
        } else {
            print("\nAll acc matched a regular expression for all analyses ($num_logic_names). Great!\n") if ($verbose);
        }

        print "\nAssigning external DB IDs to your ", uc($human_readable{$type}), "...\n" if ($verbose);
        my $assign_command = "perl $assign_db_id_script -masterhost $prod_dbhost -masterport $prod_dbport -masterdbname $prod_dbname -masteruser $prod_dbuser -host $dbhost -port $dbport -user $dbuser -pass $dbpass -dbname $dbname -conf $cfg_file -feature_type $moltype{$type} -dumpdir $output_path -update_only_null_rows";
        $assign_command .= " -uniprot_filename $uniprot_filename" if ($type eq 'protein_align_feature');
        if (system($assign_command)) {
            throw("Could not execute $assign_command\n");
        }

# check
        my $dump_fixed = "$output_path/table_dumps_$dbname.$moltype{$type}.fixed";
        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -q -e'SELECT * FROM $type where external_db_id IS NOT NULL and external_db_id != 0' >> $dump_fixed")) {
            throw("Could not dump $type rows where external_db_id is not null\n");
        }
        $num_fixed = line_count($dump_fixed);
        if ($num_ori != $num_fixed) {
            throw("The number of '$type' rows ($num_ori) does not match the number of external DB ID fixed features ($num_fixed).");
        } else {
            print("The number of '$type' rows ($num_ori) matches the number of external DB ID fixed features ($num_fixed). Great!\n") if ($verbose);
            print("Output file: $dump_fixed\n") if ($verbose);
        }
        $dumped_file = $dump_fixed;
    }
# ===============================
# DUMP FEATURE TABLES
# ===============================

        print("\n===== DUMP SUPPORTING FEATURE TABLES =====\n") if ($verbose);
        my $tsf_dumped_file = $output_path."/$synonyms{$type}_tsf.dump";
        my $sf_dumped_file = $output_path."/$synonyms{$type}_sf.dump";

        print "\nUsing file: $dumped_file\n" if ($verbose);
        push(@files_to_delete, $dumped_file);

# tsf
        print("\nDumping 'transcript_supporting_feature' table...\n") if ($verbose);
        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -N -q -e\"SELECT transcript_id,feature_type,feature_id FROM transcript_supporting_feature WHERE feature_type = '$type'\" > $tsf_dumped_file")) {
            throw("Could not dump TSF features\n");
        }
        push(@files_to_delete, $tsf_dumped_file);
# check
        my $num_tsf = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e"SELECT count(*) FROM transcript_supporting_feature WHERE feature_type = '$type'"`);
        my $num_dumped_tsf = line_count($tsf_dumped_file);
        if ($num_tsf != $num_dumped_tsf) {
          throw("The number of 'transcript_supporting_feature' rows ($num_tsf) does not match the number of dumped features ($num_dumped_tsf).");
        } else {
          print("The number of 'transcript_supporting_feature' rows ($num_tsf) matches the number of dumped features ($num_dumped_tsf). Great!\n") if ($verbose);
          print("Output file: $tsf_dumped_file`\n") if ($verbose);
        }

# sf
        print("\nDumping 'supporting_feature' table...\n") if ($verbose);
        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -q -e\"SELECT exon_id,feature_type,feature_id FROM supporting_feature WHERE feature_type = '$type'\" > $sf_dumped_file")) {
            throw("Could not dump SF features\n");
        }
        push(@files_to_delete, $sf_dumped_file);
# check
        my $num_sf = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e"SELECT count(*) FROM supporting_feature WHERE feature_type = '$type'"`);
        my $num_dumped_sf = line_count($sf_dumped_file);
        if ($num_sf != $num_dumped_sf) {
          throw("The number of 'supporting_feature' rows ($num_sf) does not match the number of dumped features ($num_dumped_sf).");
        } else {
          print("The number of 'supporting_feature' rows ($num_sf) matches the number of dumped features ($num_dumped_sf). Great!\n") if ($verbose);
          print("Output file: $sf_dumped_file`\n") if ($verbose);
        }
# =================
# SORT FILES
# =================

        print("\n===== SORT $synonyms{$type} FILES =====\n") if ($verbose);

        print("\nSorting $dumped_file...\n") if ($verbose);
        if (system("sort -u $dumped_file > $dumped_file.unique")) {
          throw("Failed sorting $dumped_file.unique");
        }
        push(@files_to_delete, $dumped_file.'.unique');

        my $num_sorted_dumped = line_count("$dumped_file.unique");
        if ($num_sorted_dumped != $num_fixed) {
          throw("The number of sorted dumped features ($num_sorted_dumped) does not match the number of dumped features ($num_fixed) for $dumped_file.unique and $dumped_file");
        } else {
          print("The number of sorted dumped features ($num_sorted_dumped) matches the number of dumped features ($num_fixed). Great!\n") if ($verbose);
          print("Output file: $dumped_file.unique`\n") if ($verbose);
        }

        print("\nSorting $dumped_file.unique by seq_region_id, seq_region_start and seq_region_end...\n") if ($verbose);
        if (system("sort -nk2 -nk3 -nk4 $dumped_file.unique > $dumped_file.unique.sorted")) {
          throw("Failed sorting $dumped_file.unique.sorted");
        }
        push(@files_to_delete, $dumped_file.'.unique.sorted');

        my $num_resorted_dumped = line_count("$dumped_file.unique.sorted");
        if ($num_resorted_dumped != $num_sorted_dumped) {
          throw("The number of resorted dumped features ($num_resorted_dumped) does not match the number of sorted dumped features ($num_sorted_dumped) for $dumped_file.unique.sorted and $dumped_file.unique");
        } else {
          print("The number of resorted dumped features ($num_resorted_dumped) matches the number of sorted dumped features ($num_sorted_dumped). Great!\n") if ($verbose);
          print("Output file: $dumped_file.unique.sorted`\n") if ($verbose);
        }

# ==========================
# FIX SUPPORT EVIDENCE LINKS
# ==========================

        my $fixed_tsf_file = $output_path."/$synonyms{$type}.fixed.tsf";
        my $fixed_sf_file = $output_path."/$synonyms{$type}.fixed.sf";
        my $tsf_fixed_tsf_file = $output_path."/$synonyms{$type}.tsf.fixed.tsf";
        my $sf_fixed_sf_file = $output_path."/$synonyms{$type}.sf.fixed.sf";

        print("\n===== FIX $type SUPPORT EVIDENCE LINKS =====\n") if ($verbose);

        print("\nFixing tsf supporting evidence links...\n") if ($verbose);
        if (system("perl $fix_supporting_evidence_script -outdaf $fixed_tsf_file -outsf $tsf_fixed_tsf_file -indaf $dumped_file.unique.sorted -insf $tsf_dumped_file")) {
            throw("Could not execute $fix_supporting_evidence_script\n");
        }
        push(@files_to_delete, $fixed_tsf_file, $fixed_sf_file);
# checks
        print("Checking tsf supporting evidence links...\n") if ($verbose);
        my $num_fixed_tsf = line_count($fixed_tsf_file);
        if ($num_fixed_tsf != $num_resorted_dumped) {
          throw("The number of fixed resorted dumped features ($num_fixed_tsf) does not match the number of resorted dumped features ($num_resorted_dumped) for $fixed_tsf_file and $dumped_file.unique.sorted");
        } else {
          print("The number of fixed resorted dumped features ($num_fixed_tsf) matches the number of resorted dumped features ($num_resorted_dumped). Great!\n") if ($verbose);
          print("Output file: $fixed_tsf_file\n") if ($verbose);
        }

        my $num_tsf_fixed_tsf = line_count($tsf_fixed_tsf_file);
        if ($num_tsf_fixed_tsf != $num_dumped_tsf) {
          throw("The number of fixed resorted dumped features ($num_tsf_fixed_tsf) does not match the number of resorted dumped features ($num_dumped_tsf) for $tsf_fixed_tsf_file and $tsf_dumped_file");
        } else {
          print("The number of fixed resorted dumped features ($num_tsf_fixed_tsf) matches the number of dumped features ($num_dumped_tsf). Great!\n") if ($verbose);
          print("Output file: $tsf_fixed_tsf_file\n") if ($verbose);
        }
        push(@files_to_delete, $tsf_fixed_tsf_file);

        print("\nFixing sf supporting evidence links...\n") if ($verbose);
        if (system("perl $fix_supporting_evidence_script -outdaf $fixed_sf_file -outsf $sf_fixed_sf_file -indaf $dumped_file.unique.sorted -insf $sf_dumped_file")) {
          throw("Could not execute $fix_supporting_evidence_script");
        }
# checks
        print("Checking sf supporting evidence links...\n") if ($verbose);
        my $num_fixed_sf = line_count($fixed_sf_file);
        if ($num_fixed_sf != $num_resorted_dumped) {
          throw("The number of fixed resorted dumped features ($num_fixed_sf) does not match the number of resorted dumped features ($num_resorted_dumped) for $fixed_sf_file and $dumped_file.unique.sorted");
        } else {
          print("The number of fixed resorted dumped features ($num_fixed_sf) matches the number of resorted dumped features ($num_resorted_dumped). Great!\n") if ($verbose);
          print("Output file: $fixed_sf_file\n") if ($verbose);
        }

        my $num_sf_fixed_sf = line_count($sf_fixed_sf_file);
        if ($num_sf_fixed_sf != $num_dumped_sf) {
          throw("The number of fixed resorted dumped features ($num_sf_fixed_sf) does not match the number of dumped features ($num_dumped_sf) for $sf_fixed_sf_file and $sf_dumped_file");
        } else {
          print("The number of fixed resorted dumped features ($num_sf_fixed_sf) matches the number of dumped features ($num_dumped_sf). Great!\n") if ($verbose);
          print("Output file: $sf_fixed_sf_file\n") if ($verbose);
        }
        push(@files_to_delete, $sf_fixed_sf_file);

# additional check
        print("\nChecking $fixed_tsf_file and $fixed_sf_file...\n") if ($verbose);
        my $comm_check = line_count("comm -12 $fixed_tsf_file $fixed_sf_file | ");
        if ($comm_check != $num_resorted_dumped) {
          throw("The files $fixed_tsf_file and $fixed_sf_file are not identical");
        } else {
          print("The files $fixed_tsf_file and $fixed_sf_file are identical. Great!\n") if ($verbose);
        }

        print("\n===== REPLACING ALL NULL WITH \\N =====\n") if ($verbose);
        # sf and tsf file for the protein_align_feature table are the same
        if (system("sed -i 's/NULL/\\\\N/g' $fixed_sf_file")) {
            throw("Failed to replace NULL values");
        }
# =====================================
# DUMP dna_align_feature RELATIONSHIPS
# =====================================

        print("\n===== DUMP $type RELATIONSHIPS =====\n") if ($verbose);

        print("\nDumping $type - SF - exon relationships...\n") if ($verbose);
        my $check_sf_old = "$output_path/check_$synonyms{$type}_sf_old.ls";
        if (system("mysql -NB -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e\"SELECT sf.exon_id, hit_name FROM supporting_feature sf, $type xaf WHERE sf.feature_type = '$type' AND sf.feature_id = xaf.".$type."_id ORDER BY sf.exon_id,xaf.hit_name\" > $check_sf_old")) {
            throw("Could not dump $type - SF relations");
        }
        push(@files_to_delete, $check_sf_old);

        my $num_sf_old = line_count($check_sf_old);
        if ($num_sf_old != $num_sf_fixed_sf) {
          throw("The $type - SF - exon relationships dump failed. Output file: $check_sf_old");
        } else {
          print("The $type - SF - exon relationships were dumped successfully. Output file: $check_sf_old\n") if ($verbose);
        }

        print("\nDumping $type - TSF - transcript relationships...\n") if ($verbose);
        my $check_tsf_old = "$output_path/check_$synonyms{$type}_tsf_old.ls";
        if (system("mysql -NB -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e\"SELECT tsf.transcript_id, hit_name FROM transcript_supporting_feature tsf,$type xaf WHERE tsf.feature_type = '$type' AND tsf.feature_id = xaf.".$type."_id ORDER BY tsf.transcript_id,xaf.hit_name\" > $check_tsf_old")) {
            throw("Could not dump $type - TSF relations");
        }
        push(@files_to_delete, $check_tsf_old);

        my $num_tsf_old = line_count($check_tsf_old);
        if ($num_tsf_old != $num_tsf_fixed_tsf) {
          throw("The $type - TSF - transcript relationships dump failed. Output file: $check_tsf_old");
        } else {
          print("The $type - TSF - transcript relationships were dumped successfully. Output file: $check_tsf_old\n") if ($verbose);
        }

# ================================
# TRUNCATE dna_align_feature table
# ================================

        print("\n===== TRUNCATE $type TABLE =====\n");

        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e 'TRUNCATE TABLE $type'")) {
          throw("The $type table truncation failed.");
        } else {
          print("The $type table was truncated successfully.\n") if ($verbose);
        }
        foreach my $table ('supporting_feature', 'transcript_supporting_feature') {
            print("\n===== DELETE $type FROM $table =====\n") if ($verbose);

            my $num_sf_before = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e"SELECT count(*) FROM $table WHERE feature_type = '$type'"`);
            print("\nDeleting from $table where feature_type = '$type'...\n") if ($verbose);
            if(system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e \"DELETE FROM $table WHERE feature_type = '$type'\"")) {
              throw("The supporting features deletion failed");
            } else {
              print("The supporting features deletion ran successfully.\n") if ($verbose);
            }
        }
# ===============================
# LOAD SORTED FILES INTO DATABASE
# ===============================

        print("\n===== LOAD SORTED FILES INTO DATABASE =====\n");
        print("\nLoading $fixed_sf_file...\n") if ($verbose);
        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e\"LOAD DATA LOCAL INFILE '$fixed_sf_file' INTO TABLE $type\"")) {
            throw("Failed to load $type table\n");
        }

        my $num_after = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e"SELECT count(*) FROM $type"`);
        if ($num_ori != $num_after) {
          throw("The number of '$type' rows ($num_after) does not match the previous number of rows ($num_ori).");
        } else {
          print("The number of '$type' rows ($num_after) matches the previous number of rows ($num_ori). Great!\n") if ($verbose);
        }

        if (!$no_external_db) {
          my $num_db_null = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'SELECT count(*) FROM $type WHERE external_db_id IS NULL'`);
          if ($num_db_null > 0) {
            throw("$num_db_null external DB IDs are NULL after loading the $type table data.");
          } else {
            print("\n$num_db_null external DB IDs are NULL after loading the $type table data. Great!\n") if ($verbose);
          }
        }

        print("Loading $tsf_fixed_tsf_file...\n") if ($verbose);
        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e\"LOAD DATA LOCAL INFILE '$tsf_fixed_tsf_file' INTO TABLE transcript_supporting_feature\"")) {
            throw("Could not load tsf data for $type\n");
        }
        print("Loading $sf_fixed_sf_file...\n") if ($verbose);
        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e\"LOAD DATA LOCAL INFILE '$sf_fixed_sf_file' INTO TABLE supporting_feature\"")) {
            throw("Could not load sf data for $type\n");
        }


# ===============================
# FINAL CHECKS
# ===============================

#DAF
        print("\n===== FINAL CHECKS $type =====\n") if ($verbose);

        print("\nChecking $type to sf, sf to exon...\n") if ($verbose);

        my $check_sf_new = $output_path."/check_$synonyms{$type}_sf_new.ls";
        if (system("mysql -NB -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e\"SELECT sf.exon_id, hit_name FROM supporting_feature sf, $type xaf WHERE sf.feature_type = '$type' AND sf.feature_id = xaf.".$type."_id ORDER BY sf.exon_id,xaf.hit_name\" > $check_sf_new")) {
            throw("Could not dump new SF");
        }
        push(@files_to_delete, $check_sf_new);

        if (system("diff -q $check_sf_old $check_sf_new")) {
          throw("The $type optimization failed since $check_sf_old is different from $check_sf_new");
        } else {
          print("$check_sf_old is the same as $check_sf_new. Great!\n") if ($verbose);
        }

        print("\nChecking $type to tsf, tsf to transcript...\n") if ($verbose);

        my $check_tsf_new = $output_path."/check_$synonyms{$type}_tsf_new.ls";
        if (system("mysql -NB -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e\"SELECT tsf.transcript_id, hit_name FROM transcript_supporting_feature tsf,$type xaf WHERE tsf.feature_type = '$type' AND tsf.feature_id = xaf.".$type."_id ORDER BY tsf.transcript_id,xaf.hit_name\" > $check_tsf_new")) {
            throw("Could not dump new TSF");
        }
        push(@files_to_delete, $check_tsf_new);

        if (system("diff $check_tsf_old $check_tsf_new")) {
          throw("The $type optimization failed since $check_tsf_old is different from $check_tsf_new");
        } else {
          print("$check_tsf_old is the same as $check_tsf_new. Great!\n") if ($verbose);
        }

        print("\nThe $type table optimization finished SUCCESSFULLY.\n");



    }

    if ($do_ise) {
      my $ise_name = catfile($output_path, 'intron_supporting_evidence.dat');
      my $tise_name = catfile($output_path, 'transcript_intron_supporting_evidence.dat');
      my $check_old_name = catfile($output_path, 'check_old.dat');
      my $check_new_name = catfile($output_path, 'check_new.dat');
      foreach my $file ($ise_name, $tise_name, $check_old_name) {
        throw("'$file' already exists, something went wrong and needs to be fixed") if (-e $file);
      }
      my $check_query = 'SELECT t.seq_region_id, t.seq_region_start, t.seq_region_end, t.seq_region_strand, i.seq_region_id, i.seq_region_start, i.seq_region_end, i.seq_region_strand FROM transcript t LEFT JOIN transcript_intron_supporting_evidence tise ON t.transcript_id = tise.transcript_id LEFT JOIN intron_supporting_evidence i ON tise.intron_supporting_evidence_id = i.intron_supporting_evidence_id';
      foreach my $query (
      "'SELECT * FROM intron_supporting_evidence' | sort -nk3 -nk4 -nk5 -nk6 > ".$ise_name,
      "'SELECT * FROM transcript_intron_supporting_evidence' > ".$tise_name,
      "'$check_query' | sort -nk1 -nk2 -nk3 -nk4 > ".$check_old_name,
      ) {
        if (system("mysql -NB -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e $query")) {
          throw("could not execute $query");
        }
      }
      if (-s $ise_name) {
        my $cmd = 'perl '.$ise_script
            .' -indaf '.$ise_name
            .' -insf '.$tise_name
            .' -outdaf '.$ise_name.'.fixed'
            .' -outsf '.$tise_name.'.fixed';
        if ($collapse_ise) {
          my $analysis_id = `mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e "SELECT analysis_id FROM analysis WHERE logic_name LIKE '%merged_rnaseq_ise'"`;
          if ($analysis_id and $analysis_id > 0) {
            $cmd .= "-analysis_id $analysis_id";
          }
          else {
            throw('Failed to get analysis_id for "%merged_rnaseq_ise"');
          }
        }
        if (system($cmd)) {
          throw('Could not prepare ISEs');
        }
        foreach my $table ('intron_supporting_evidence', 'transcript_intron_supporting_evidence') {
          foreach my $query ("TRUNCATE $table", "LOAD DATA LOCAL INFILE \"".catfile($output_path, $table.'.dat.fixed')."\" INTO TABLE $table") {
            if (system("mysql -NB -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e '$query'")) {
              throw("could not execute $query");
            }
          }
        }
        my $query = "'SELECT t.seq_region_id, t.seq_region_start, t.seq_region_end, t.seq_region_strand, i.seq_region_id, i.seq_region_start, i.seq_region_end, i.seq_region_strand FROM transcript t LEFT JOIN transcript_intron_supporting_evidence tise ON t.transcript_id = tise.transcript_id LEFT JOIN intron_supporting_evidence i ON tise.intron_supporting_evidence_id = i.intron_supporting_evidence_id' | sort -nk1 -nk2 -nk3 -nk4 > ".$check_new_name;
        if (system("mysql -NB -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e $query")) {
          throw("could not execute $query");
        }
        if (system("diff -q $check_old_name $check_new_name")) {
          throw("Optimising ISE failed, you may need to reload $ise_name and $tise_name before restarting");
        }
      }
      else {
        warning("The intron_supporting_evidence table is empty");
      }
    }




if (!$no_backup and !$clean) {
  print("\nRemember that you can undo the changes made to your database by running the following command:\n");
  print("mysql -h$dbhost -P$dbport -u$dbuser -p*** -D$dbname < $db_backup_file\n");
}
elsif ($clean) {
    foreach my $file (@files_to_delete) {
        if (system("rm $file")) {
            warning("Could not delete $file");
        }
    }
  print("\nThe backup (if any) and temporary files have been deleted because the 'clean' option was set.\n");
}


#-------
# END
#-------

sub line_count {
  my ($file) = @_;

  my $count = 0;
  open(RH, "$file") || throw("Could not open '$file'");
  while(<RH>) {
    ++$count;
  }
  close(RH) || throw("Could not close '$file'");
  return $count;
}

sub usage {
    print <<EOF

Usage:

$0 -output_path <output_path> -dbhost <dbhost> -dbport <dbport> -dbname <dbname> -dbuser <dbuser> -dbpass <dbpass> -prod_dbhost <prod_dbhost> -prod_dbname <prod_dbname> -prod_dbuser <prod_dbuser> -prod_dbport <prod_dbport> -analysis_scripts <analysis_scripts> -uniprot_filename <uniprot_filename> [-daf 0] [-paf 0] [-verbose] [-help]

-output_path	Path where the output files and backup files will be written. It will be created if it does not exist.

-dbhost    		host name where the database is located

-dbport    		port number

-dbname    		database name

-dbuser    		what username to connect as

-dbpass    		what password to use

-dbhost       production db host name

-dbport       production db port number

-dbname       production database name

-dbuser       what username to connect as for the production db

-dbpass       what password to use for the production db

-analysis_scripts	path to ensembl-analysis/scripts, needed to run fix_supporting_evidence_links.pl,  or deduced from this scripts path if not specified

-uniprot_filename	full path to the uniprot filename required by the script which assigns the external DB IDs

-no_external_db	skip the external db assignments, do the sorting only

-no_backup		skip the backup

-clean			  delete any backup and output file once the final checks step has finished successfully

-daf          By default both daf and paf tables are optimised. Use "-daf 0" to not do daf.

-paf          By default both daf and paf tables are optimised. Use "-paf 0" to not do paf.

-verbose      Use this option to get more print statements to follow the script.

-help			    Show usage.


Examples:
# assign external DB IDs and sort features tables
bsub -M 3700 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af.out -e optimize_af.err "perl load_external_db_ids_and_optimize_af.pl -output_path /lustre/scratch101/sanger/cgg/CanFam3.1/optimize -dbhost genebuild1 -dbport 3306 -dbname cgg_dog_ref_test -dbuser ensadmin -dbpass *** -analysis_scripts ~/enscode/ensembl-analysis/scripts/genebuild -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose"

# sort features tables only
bsub -M 3700 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af_after_alt_seq_mapping.out -e optimize_af_after_alt_seq_mapping.err "perl load_external_db_ids_and_optimize_af.pl -output_path optimize_core_af_alt_seq_mapping -dbhost ens-staging1 -dbport 3306 -dbname homo_sapiens_core_70_37 -dbuser ensadmin -dbpass *** -analysis_scripts ~/enscode/ensembl-analysis/scripts/genebuild -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose -no_external_db"

# assign external DB IDs and sort features tables, clean
bsub -M 3700 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af.out -e optimize_af.err "perl load_external_db_ids_and_optimize_af.pl -output_path /lustre/scratch101/sanger/cgg/optimize -dbhost ens-staging1 -dbport 3306 -dbname homo_sapiens_core_70_37 -dbuser ensadmin -dbpass *** -analysis_scripts ~/enscode/ensembl-analysis/scripts/genebuild -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose -clean"

EOF
}

1;
