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

# This script creates the regular expression file for
# your dna/protein align features to assign external db
# ids. It uses dictonary of regular expressions in the 
# main_regex_file file which it tries to match.
 
# It will warn if there are acc in your table which do not match to
# any regular expression. In such a case you have to edit the regular
# expression. take care, if you are too greedy with your regular
# expression it will match different external db ids and we don't want
# that.
#
# example command line :
#    perl test_regexes.pl \
#         -dbname danio_rerio_core_54_8 -dbhost ens-staging -dbuser ensro -type paf \
#         -main_regex_file paf_regexes.dat
# Use -overwrite to overwrite the config file, useful when you try to run the script
# several times.

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Copy;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);


my %opt = ( dbhost   => '',
            dbname   => '',
            dbuser   => 'ensro',
            dbpass   => '',
            dbport   => '3306',
            filename => '', );

GetOptions(
  \%opt,
  'dbhost|host|h=s',
  'dbuser|user|u=s',
  'dbpass|pass|p=s',
  'dbport|port|P=i',
  'dbname|db|D=s',
  'type=s',
  'rerun!',
  'main_regex_file=s',
  'output_config_file=s',
  'tmp_dir=s',
  'overwrite!',

);

unless ( $opt{type} ) {
  throw(" missing option : -type <paf|daf|dna|protein> \n");
}

unless ( $opt{dbhost} && $opt{dbname} ) {
  throw(" missing option : -dbhost <HOST> -dbname <DBNAME>  \n");
}

my $db = $opt{dbname};

my $dba =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new( -dbname => $opt{dbname},
                                       -user   => $opt{dbuser},
                                       -pass   => $opt{dbpass},
                                       -host   => $opt{dbhost},
                                       -port   => $opt{dbport}, );

# do 2 raw sql commands here

my $feature_table_name;
if ( $opt{type} =~ m/paf/ || $opt{type} =~ m/protein/ ) {
  $feature_table_name = "protein_align_feature";
} elsif ( $opt{type} =~ m/daf/ || $opt{type} =~ m/dna/ ) {
  $feature_table_name = "dna_align_feature";
} else {
  throw(" unrecognized option : -type <$opt{type}\t-- allowed types : [dna|daf|protein|paf]  \n"
  );
}

print "\nLooking up logic_names from database $opt{dbname} "
    . "for features TYPE $opt{type}\n";

my @logic_names;
my $query = "select logic_name from analysis a, $feature_table_name p "
          . "where a.analysis_id = p.analysis_id group by p.analysis_id ";

my $sth = $dba->dbc->prepare($query);
$sth->execute();

print "\nLogic_names in $opt{dbname} for $opt{type} : \n";

for ( @{ $sth->fetchall_arrayref() } ) {
  print " - $$_[0]\n";
  push @logic_names, $$_[0];
}

print "\n";
$sth->finish;

#  create 2nd file with hitname and logic_name
my $tmp_dir;
if ( $opt{tmp_dir} ) {
  $tmp_dir = "/" . $opt{tmp_dir};
} else {
  $tmp_dir = "/tmp/";
}

unless ( -e $tmp_dir ) {
  throw(" your temp dir $tmp_dir does not exist\n");
}

my $feature_file = $tmp_dir . "/feature_$opt{dbname}_$opt{type}.lst";

unless ( $opt{rerun} ) {
  print "Building hitname->logic_name file \n";

  # group by both hit_name and logic_name if hit_names are shared
  # between several analyses. Shouldn't change the regexes but will
  # get the correct numbers.
  $query = "select hit_name, logic_name "
         . "from analysis a, $feature_table_name p "
         . "where a.analysis_id = p.analysis_id "
         . "group by hit_name, logic_name";

  $sth = $dba->dbc->prepare($query);
  $sth->execute();

  print "reading $opt{type} table \n";
  open FF, ">$feature_file" || die " can't write to file $feature_file\n";

  for my $line ( @{ $sth->fetchall_arrayref() } ) {
    my $hit_name   = $$line[0];
    my $logic_name = $$line[1];
    print FF $hit_name . "\t" . $logic_name . "\n";
  }

  print " done - features have been written to $feature_file\n";
  close(FF);
  $sth->finish;
} else {
  print "RE-RUNNING -- I am using old file file $feature_file \n";
  throw("file $feature_file does not exist ") if ( !-e $feature_file );
}

my $output_config_file;
if ( $opt{output_config_file} ) {
  $output_config_file = $opt{output_config_file};
} else {
  $output_config_file = $db . "_" . $opt{type} . "s.cfg";
}

if ($opt{overwrite} and -e $output_config_file) {
    if (copy($output_config_file, $output_config_file.'.bak')) {
        throw("Could not create a backup of $output_config_file");
    }
}

#open FPOUT,">$output_config_file" || die " can't write to $output_config_file\n" ;
for my $lg (@logic_names) {
  print "Analysis\t" . $lg . "\n";
  get_counts_for_logic_name( $lg, $db, $opt{type}, $feature_file );
}

print "output-config file has been written to : $output_config_file\n";
print "Main reg exp file used : $opt{main_regex_file}\n";

sub get_counts_for_logic_name {
  my ( $logic_name, $db, $type, $feature_file ) = @_;

  open FPOUT,
    ">>$output_config_file" || die " can't write to $output_config_file\n";
  print FPOUT "Analysis\t" . $logic_name . "\n";

  # reads  file :
  # A0A039.1        Uniprot
  # A0A182.1        SimilarityGenewise
  # A0A181          TargetedGenewise

  my $type_file = "/tmp/test_regexes.conv.$$.$logic_name.tmp";
  `grep $logic_name $feature_file >  $type_file`;

# creates type files fore each analysis : all uniprots in 1 file, all simgw in one file etc
  print "\n";
  my $nline = `wc -l $type_file | awk '{ print \$1 }'`;
  chomp $nline;
  print "\n\nHave $nline ids for $logic_name \n";

  open FPREGEX, "<$opt{main_regex_file}";

  my $tot_matched = 0;
  my @matched_regexes;
  while (<FPREGEX>) {
    chomp;
    my @words = split /\t/, $_;
    # just a check

    if ( $feature_table_name =~ m/protein_align_feature/ ) {
      if ( $words[1] =~ m/(dna|unigene)/i ) {
        throw("are you sure you are using the correct options? "
            . "Looks like you use dna-file to match protein-accession numbers\n"
        );
      }
    }
    if ( $feature_table_name =~ m/dna_align_feature/ ) {
      if ( $words[1] =~ m/(protein|uniprot)/i ) {
        throw("are you sure you are using the correct options? "
            . "Looks like you use protein-regex-file to match dna-accession numbers\n"
        );
      }
    }

    my $pat = $words[0];
    my $cnt = `awk '{ print \$1 }' $type_file | egrep -ce $pat`;
    if ( $cnt > 0 ) {
      $tot_matched += $cnt;
      print $pat . "\t" . $words[1] . "\t" . $cnt;
      print FPOUT " Regex\t" . $pat . "\t" . $words[1] . "\n";
      push @matched_regexes, $pat;
    }
  } ## end while (<FPREGEX>)
  print "Total matched lines = $tot_matched\n\n";
  my $all_lines_matched = 1;

  if ( $tot_matched < $nline ) {
    print "\n\nWARNING !!! You have unmatched accession-numbers for $logic_name !!!\n\n";
    my $egrep_string = "awk '{ print \$1 }' $type_file";
    foreach my $pat (@matched_regexes) {
      $egrep_string .= " | egrep -v -- $pat";
    }
    print $egrep_string . "\n";
    print "Unmatched lines:\n";
    system($egrep_string);
    print FPOUT " Default\t????\n";
    $all_lines_matched = 0;
  }
  close FPREGEX;
  if ( $all_lines_matched == 1 ) {
    print "GOOD - All acc matched a regular expression\n\n";
    unlink($type_file);
  }
  close FPOUT;
} ## end sub get_counts_for_logic_name
