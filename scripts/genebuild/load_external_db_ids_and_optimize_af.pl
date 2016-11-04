#!/usr/bin/env perl

# $Source: /cvsroot/ensembl/ensembl-personal/genebuilders/scripts/load_external_db_ids_and_optimize_af.pl,v $
# $Revision: 1.35 $

=pod

=head1 NAME 

load_external_db_ids_and_optimize_af.pl

=head1 DESCRIPTION

This script dumps the dna_align_feature, transcript_support_feature and support_feature tables for optimizing the dna_align_feature table for a given database. It also assigns the external DB IDs for both the dna_align_feature and the protein_align_feature tables. The scripts test_regexes.pl, assign_external_db_ids.pl and fix_supporting_evidence_links.pl in ensembl-personal/genebuilders/scripts are required.

=head1 OPTIONS

-output_path		Path where the output files and backup files will be written. It will be created if it does not exist.

-dbhost    		host name where the database is located

-dbport    		what port to connect (default 3306)

-dbname    		database name

-dbuser    		what username to connect as

-dbpass    		what password to use

-ensgbscripts	path to local ensembl-personal/genebuilders/scripts, read from ENSGBSCRIPTS environment variable or deduced from this scripts path if not specified

-uniprot_filename	full path to the uniprot filename required by the script which assigns the external DB IDs

-no_external_db		skip the external db assignments, only do the sorting

-no_backup		skip the backup

-clean			delete any backup and output file once the final checks step has finished successfully


=head2 Output options:

        -verbose        Use this option to get more print statements to follow
                        the script. Set to 0 (not verbose) by default to get
                        only the final summary.

	-help		Show usage.

=head1 EXAMPLE USAGE

=head1

# assign external DB IDs and sort features tables
bsub -M 3700000 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af.out -e optimize_af.err "perl load_external_db_ids_and_optimize_af.pl -prod_dbuser *** -prod_dbpass *** -prod_dbhost *** -prod_dbname *** -prod_dbport *** -output_path /lustre/scratch101/sanger/cgg/CanFam3.1/optimize -dbhost genebuild1 -dbname cgg_dog_ref_test -dbuser ensadmin -dbpass *** -ensgbscripts /nfs/users/nfs_c/cgg/ensembl-personal/genebuilders/scripts -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose"

# sort features tables only
bsub -M 3700000 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af_after_alt_seq_mapping.out -e optimize_af_after_alt_seq_mapping.err "perl load_external_db_ids_and_optimize_af.pl -prod_dbuser *** -prod_dbpass *** -prod_dbhost *** -prod_dbname *** -prod_dbport *** -output_path optimize_core_af_alt_seq_mapping -dbhost ens-staging1 -dbname homo_sapiens_core_70_37 -dbuser ensadmin -dbpass *** -ensgbscripts /nfs/users/nfs_c/cgg/ensembl-personal/genebuilders/scripts -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose -no_external_db"

# assign external DB IDs and sort features tables, clean
bsub -M 3700000 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af.out -e optimize_af.err "perl load_external_db_ids_and_optimize_af.pl -prod_dbuser *** -prod_dbpass *** -prod_dbhost *** -prod_dbname *** -prod_dbport *** -output_path /lustre/scratch101/sanger/cgg/optimize -dbhost ens-staging1 -dbname homo_sapiens_core_70_37 -dbuser ensadmin -dbpass *** -ensgbscripts /nfs/users/nfs_c/cgg/ensembl-personal/genebuilders/scripts -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose -clean"

=cut


use strict;
use warnings;

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use File::Basename;
use List::Util qw(sum);

my $output_path;
my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my $prod_dbhost;
my $prod_dbuser;
my $prod_dbpass;
my $prod_dbport = 3306;
my $prod_dbname;
my $dir_ensgbscripts;
my $uniprot_filename;
my $no_external_db = 0;
my $no_backup = 0;
my $clean = 0;
my $help = 0;
my $verbose = 0;
my $do_daf = 1;
my $do_paf = 1;

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
           'ensgbscripts=s'     => \$dir_ensgbscripts, 
           'uniprot_filename=s' => \$uniprot_filename, 
           'no_external_db!'=> \$no_external_db, 
           'no_backup!'    => \$no_backup,       
           'clean!'        => \$clean,           
           'daf!'        => \$do_daf,           
           'paf!'        => \$do_paf,           
           'help'          => \$help,
           'verbose'       => \$verbose);
print $0, "\n";

($dir_ensgbscripts) = $0 =~ /(.*)\/[^\/]+$/ unless $dir_ensgbscripts ;
if (!$output_path or !$dbhost or !$dbname or !$dbuser or !$dbpass or !$dir_ensgbscripts or $help) 
{
    &usage;
    exit(1);
}

if( !$no_external_db )
{
    if (!$uniprot_filename) {
        warn("You haven't specified if you want the external db id to be populated...\nUniprot filename: $uniprot_filename\nno_external_db: $no_external_db");
        &usage;
        exit(1);
    }
    throw( "uniprot_filename needs to be a file (typically like /data/blastdb/Ensembl/uniprot_yyyy_mm/entry_loc)" ) unless ( -f $uniprot_filename);
}


#-------
# BEGIN
#-------
$output_path =~ s/\/$//;
$dir_ensgbscripts =~ s/\/$//;

if (system("mkdir -p $output_path")) {
  throw("Cannot create output_path $output_path.");
}

my $num_fixed = 0; # will be calculated later
my $fix_supporting_evidence_script = $dir_ensgbscripts.'/fix_supporting_evidence_links.pl';
my $test_regex_script = $dir_ensgbscripts.'/test_regexes.pl';
my $assign_db_id_script = $dir_ensgbscripts.'/assign_external_db_ids.pl';
foreach my $script ($fix_supporting_evidence_script, $test_regex_script, $assign_db_id_script) {
    throw("Cannot access $script") unless (-e $script);
}


print("\nhost: $dbhost\n") if ($verbose);
print("port: $dbport\n") if ($verbose);
print("name: $dbname\n") if ($verbose);
print("uniprot file: $uniprot_filename\n") if ( $verbose && !$no_external_db ) ;

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
        $num_fixed = int(`wc -l $dumped_file | awk '{print \$1}'`);
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
        if (system("perl $test_regex_script -dbname $dbname -dbhost $dbhost -dbuser $dbuser -dbpass $dbpass -type $synonyms{$type} -main_regex_file $dir_ensgbscripts/$synonyms{$type}_regexes.dat -output_config_file $cfg_file > $cfg_log_file")) {
            throw("Could not execute $test_regex_script\n");
        }

        print("Output files: $cfg_file and $cfg_log_file\n") if ($verbose);

        my $num_logic_names = int(`grep Analysis $cfg_file | sort | uniq | wc -l | awk '{print \$1}'`);
        my $num_all_matched = int(`grep -c "GOOD - All acc matched a regular expression" $cfg_log_file`);

        if ($num_logic_names <= 0) {
            warning("Could not find any logic name for $human_readable{$type} analyses.");
        } elsif ($num_logic_names != $num_all_matched) {
            throw("Some acc do not match a regular expression. The number of logic names found for $human_readable{$type} ($num_logic_names) does not match the number of 'GOOD - All acc matched a regular expression' ($num_all_matched).\nPLEASE, HAVE A LOOK AT $cfg_file and $cfg_log_file\n");
        } else {
            print("\nAll acc matched a regular expression for all analyses ($num_logic_names). Great!\n") if ($verbose);
        }

        print "\nAssigning external DB IDs to your ", uc($human_readable{$type}), "...\n" if ($verbose);
        if (system("perl $assign_db_id_script -masterhost $prod_dbhost -masterport $prod_dbport -masterdbname $prod_dbname -masteruser $prod_dbuser -host $dbhost -user $dbuser -pass $dbpass -dbname $dbname -conf $cfg_file -feature_type $moltype{$type} -dumpdir $output_path -update_only_null_rows -uniprot_filename $uniprot_filename")) {
        #if (system("bsub  -M 3700000 -R 'select[mem>3700] rusage[mem=3700]' -I perl $assign_db_id_script -host $dbhost -pass $dbpass -dbname $dbname -conf $cfg_file -feature_type $moltype{$type} -dumpdir $output_path -update_only_null_rows -uniprot_filename $uniprot_filename")) {
            throw("Could not execute $assign_db_id_script\n");
        }

# check
        my $dump_fixed = "$output_path/table_dumps_$dbname.$moltype{$type}.fixed";
        if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -q -e'SELECT * FROM $type where external_db_id IS NOT NULL and external_db_id != 0' >> $dump_fixed")) {
            throw("Could not dump $type rows where external_db_id is not null\n");
        }
        $num_fixed = int(`wc -l $dump_fixed | awk '{print \$1}'`);
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
        my $num_dumped_tsf = int(`wc -l $tsf_dumped_file | awk '{print \$1}'`);
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
        my $num_dumped_sf = int(`wc -l $sf_dumped_file | awk '{print \$1}'`);
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
        `sort -u $dumped_file > $dumped_file.unique`;
        push(@files_to_delete, $dumped_file.'.unique');

        my $num_sorted_dumped = int(`wc -l $dumped_file.unique | awk '{print \$1}'`);
        if ($num_sorted_dumped != $num_fixed) {
          throw("The number of sorted dumped features ($num_sorted_dumped) does not match the number of dumped features ($num_fixed) for $dumped_file.unique and $dumped_file");
        } else {
          print("The number of sorted dumped features ($num_sorted_dumped) matches the number of dumped features ($num_fixed). Great!\n") if ($verbose);
          print("Output file: $dumped_file.unique`\n") if ($verbose);
        }

        print("\nSorting $dumped_file.unique by seq_region_id, seq_region_start and seq_region_end...\n") if ($verbose);
        `sort -nk2 -nk3 -nk4 $dumped_file.unique > $dumped_file.unique.sorted`;
        push(@files_to_delete, $dumped_file.'.unique.sorted');

        my $num_resorted_dumped = int(`wc -l $dumped_file.unique.sorted | awk '{print \$1}'`);
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
        my $num_fixed_tsf = int(`wc -l $fixed_tsf_file | awk '{print \$1}'`);
        if ($num_fixed_tsf != $num_resorted_dumped) {
          throw("The number of fixed resorted dumped features ($num_fixed_tsf) does not match the number of resorted dumped features ($num_resorted_dumped) for $fixed_tsf_file and $dumped_file.unique.sorted");
        } else {
          print("The number of fixed resorted dumped features ($num_fixed_tsf) matches the number of resorted dumped features ($num_resorted_dumped). Great!\n") if ($verbose);
          print("Output file: $fixed_tsf_file\n") if ($verbose);
        }

        my $num_tsf_fixed_tsf = int(`wc -l $tsf_fixed_tsf_file | awk '{print \$1}'`);
        if ($num_tsf_fixed_tsf != $num_dumped_tsf) {
          throw("The number of fixed resorted dumped features ($num_tsf_fixed_tsf) does not match the number of resorted dumped features ($num_dumped_tsf) for $tsf_fixed_tsf_file and $tsf_dumped_file");
        } else {
          print("The number of fixed resorted dumped features ($num_tsf_fixed_tsf) matches the number of dumped features ($num_dumped_tsf). Great!\n") if ($verbose);
          print("Output file: $tsf_fixed_tsf_file\n") if ($verbose);
        }
        push(@files_to_delete, $tsf_fixed_tsf_file);

        print("\nFixing sf supporting evidence links...\n") if ($verbose);
        `perl $fix_supporting_evidence_script -outdaf $fixed_sf_file -outsf $sf_fixed_sf_file -indaf $dumped_file.unique.sorted -insf $sf_dumped_file`;
# checks
        print("Checking sf supporting evidence links...\n") if ($verbose);
        my $num_fixed_sf = int(`wc -l $fixed_sf_file | awk '{print \$1}'`);
        if ($num_fixed_sf != $num_resorted_dumped) {
          throw("The number of fixed resorted dumped features ($num_fixed_sf) does not match the number of resorted dumped features ($num_resorted_dumped) for $fixed_sf_file and $dumped_file.unique.sorted");
        } else {
          print("The number of fixed resorted dumped features ($num_fixed_sf) matches the number of resorted dumped features ($num_resorted_dumped). Great!\n") if ($verbose);
          print("Output file: $fixed_sf_file\n") if ($verbose);
        }

        my $num_sf_fixed_sf = int(`wc -l $sf_fixed_sf_file | awk '{print \$1}'`);
        if ($num_sf_fixed_sf != $num_dumped_sf) {
          throw("The number of fixed resorted dumped features ($num_sf_fixed_sf) does not match the number of dumped features ($num_dumped_sf) for $sf_fixed_sf_file and $sf_dumped_file");
        } else {
          print("The number of fixed resorted dumped features ($num_sf_fixed_sf) matches the number of dumped features ($num_dumped_sf). Great!\n") if ($verbose);
          print("Output file: $sf_fixed_sf_file\n") if ($verbose);
        }
        push(@files_to_delete, $sf_fixed_sf_file);

# additional check
        print("\nChecking $fixed_tsf_file and $fixed_sf_file...\n") if ($verbose);
        my $comm_check = int(`comm -12 $fixed_tsf_file $fixed_sf_file | wc -l | awk '{print \$1}'`);
        if ($comm_check != $num_resorted_dumped) {
          throw("The files $fixed_tsf_file and $fixed_sf_file are not identical");
        } else {
          print("The files $fixed_tsf_file and $fixed_sf_file are identical. Great!\n") if ($verbose);
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

        my $num_sf_old = int(`wc -l $check_sf_old | awk '{print \$1}'`);
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

        my $num_tsf_old = int(`wc -l $check_tsf_old | awk '{print \$1}'`);
        if ($num_tsf_old != $num_tsf_fixed_tsf) {
          throw("The $type - TSF - transcript relationships dump failed. Output file: $check_tsf_old");
        } else {
          print("The $type - TSF - transcript relationships were dumped successfully. Output file: $check_tsf_old\n") if ($verbose);
        }
# ================================
# TRUNCATE dna_align_feature table
# ================================

        print("\n===== TRUNCATE $type TABLE =====\n") if ($verbose);

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

        print("\n===== REPLACING ALL NULL WITH \\N =====\n") if ($verbose);
        # sf and tsf file for the protein_align_feature table are the same
        if (system("sed -i 's/NULL/\\\\N/g' $fixed_sf_file")) {
            throw("Failed to replace NULL values");
        }
        print("\n===== LOAD SORTED FILES INTO DATABASE =====\n") if ($verbose);
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

        if (`diff $check_sf_old $check_sf_new` eq "") {
          print("$check_sf_old is the same as $check_sf_new. Great!\n") if ($verbose);
        } else {
          throw("The $type optimization failed since $check_sf_old is different from $check_sf_new");
        }

        print("\nChecking $type to tsf, tsf to transcript...\n") if ($verbose);

        my $check_tsf_new = $output_path."/check_$synonyms{$type}_tsf_new.ls";
        if (system("mysql -NB -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e\"SELECT tsf.transcript_id, hit_name FROM transcript_supporting_feature tsf,$type xaf WHERE tsf.feature_type = '$type' AND tsf.feature_id = xaf.".$type."_id ORDER BY tsf.transcript_id,xaf.hit_name\" > $check_tsf_new")) {
            throw("Could not dump new TSF");
        }
        push(@files_to_delete, $check_tsf_new);

        if (`diff $check_tsf_old $check_tsf_new` eq "") {
          print("$check_tsf_old is the same as $check_tsf_new. Great!\n") if ($verbose);
        } else {
          throw("The $type optimization failed since $check_tsf_old is different from $check_tsf_new");
        }

        print("\nThe $type table optimization finished SUCCESSFULLY.\n");



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

sub usage {
    print <<EOF

Usage:

$0 -output_path <output_path> -dbhost <dbhost> [-dbport <dbport>] -dbname <dbname> -dbuser <dbuser> -dbpass <dbpass> -ensgbscripts <ensgbscripts> -uniprot_filename <uniprot_filename> [-verbose] [-help]

-output_path	Path where the output files and backup files will be written. It will be created if it does not exist.

-dbhost    		host name where the database is located

-dbport    		port number (default 3306)

-dbname    		database name

-dbuser    		what username to connect as

-dbpass    		what password to use

-ensgbscripts	path to local ensembl-personal/genebuilders/scripts, needed to run fix_supporting_evidence_links.pl, read from ENSGBSCRIPTS environment variable or deduced from this scripts path if not specified

-uniprot_filename	full path to the uniprot filename required by the script which assigns the external DB IDs

-no_external_db	skip the external db assignments, do the sorting only

-no_backup		skip the backup

-clean			delete any backup and output file once the final checks step has finished successfully

-verbose       	Use this option to get more print statements to follow the script.

-help			Show usage.


Examples:
# assign external DB IDs and sort features tables
bsub -M 3700000 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af.out -e optimize_af.err "perl load_external_db_ids_and_optimize_af.pl -output_path /lustre/scratch101/sanger/cgg/CanFam3.1/optimize -dbhost genebuild1 -dbname cgg_dog_ref_test -dbuser ensadmin -dbpass *** -ensgbscripts /nfs/users/nfs_c/cgg/ensembl-personal/genebuilders/scripts -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose"

# sort features tables only
bsub -M 3700000 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af_after_alt_seq_mapping.out -e optimize_af_after_alt_seq_mapping.err "perl load_external_db_ids_and_optimize_af.pl -output_path optimize_core_af_alt_seq_mapping -dbhost ens-staging1 -dbname homo_sapiens_core_70_37 -dbuser ensadmin -dbpass *** -ensgbscripts /nfs/users/nfs_c/cgg/ensembl-personal/genebuilders/scripts -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose -no_external_db"

# assign external DB IDs and sort features tables, clean
bsub -M 3700000 -R 'select[mem>3700] rusage[mem=3700]' -o optimize_af.out -e optimize_af.err "perl load_external_db_ids_and_optimize_af.pl -output_path /lustre/scratch101/sanger/cgg/optimize -dbhost ens-staging1 -dbname homo_sapiens_core_70_37 -dbuser ensadmin -dbpass *** -ensgbscripts /nfs/users/nfs_c/cgg/ensembl-personal/genebuilders/scripts -uniprot_filename /data/blastdb/Ensembl/uniprot_2013_05/entry_loc -verbose -clean"

EOF
}

1;
