=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

# NCBI provides SQL Server files, which are not compatible with MySQL.
# This script converts the SQL Server files provided by NCBI to MySQL files named "new_*.sql"
# and replaces empty values between tabs with \N in "Interpretations.txt".
# It also creates the extra table "EnsemblWithdrawals" which contains a list of actions related to the withdrawal comments from the annotators. If a CCDS has been withdrawn we may or may not want to remove the transcripts for our geneset altogether. This file tells us what to do for each withdrawal.

use warnings;
use strict;
use Getopt::Long;

# directory where the input files are located and the output files will be written
my $dir = undef;

&GetOptions(
  'dir:s'   => \$dir,
);

unless (defined $dir){
  die "need to define -dir \n";
}

#3 files to edit:

my $sql = "sql/";
my $data = "data/"; 

my $create_tables = "createTables.sql";
my $interpretations = "Interpretations.txt";
my $create_keys = "createKeys.sql";

#fix create_tables

open (IN, $dir."/".$sql.$create_tables) or die "can't open $dir/".$sql.$create_tables."\n";
open (OUT, ">".$dir."/".$sql."new_".$create_tables) or die "can't open $dir/".$sql."new_".$create_tables."\n";

$/ = "\n\n";
my $line = "";

while(<IN>){
  
  next if ($_=~/ALTER TABLE/);
  
  my $entry = $_;
  $entry =~s/\[dbo\]\.//g;
  $entry =~s/CREATE\s+TABLE\s+dbo\./CREATE TABLE /g;
  $entry =~s/\bGO\b/;/gi;
  $entry =~s/\[gap_count\] \[int\] NOT NULL ,/gap_count int NOT NULL/ ;
  $entry =~s/\[//g;
  $entry =~s/\]//g;
  $entry =~s/IDENTITY\s*\(1,\s*1\)/UNIQUE AUTO_INCREMENT/g;
  $entry =~s/General_//g;
  $entry =~s/NOT FOR REPLICATION //g;
  $entry =~s/ ON PRIMARY TEXTIMAGE_ON PRIMARY/ COLLATE=latin1_swedish_ci ENGINE=MyISAM/g;
  $entry =~s/ ON PRIMARY/ COLLATE=latin1_swedish_ci ENGINE=MyISAM/g;

  $line .=  "$entry\n";

}
$line =~ s/,\s*\n*\s*\)\s*ENGINE/\n) ENGINE/g;
print OUT $line;
$line = undef;

close IN;

#Add entry for extra table detailing Withdrawals:
print OUT "\n\nCREATE TABLE EnsemblWithdrawals (\n".
	"\tccds_uid int NOT NULL ,\n".
	"\taction ENUM( 'Keep', 'Remove transcript', 'Remove gene' ), \n".
	"\tcomment text COLLATE Latin1_BIN NULL\n". 
") COLLATE=latin1_swedish_ci ENGINE=MyISAM\n".
";\n";

close OUT;

#fix interpretations

open (IN, $dir."/".$data.$interpretations) or die "can't open $dir/".$data.$interpretations."\n";
open (OUT, ">".$dir."/".$data."new_".$interpretations) or die "can't open $dir/".$data."new_".$interpretations."\n";


$/ = "\n";

while(<IN>){
  
  my $entry = $_;
  $entry =~s/\t\t/\t\\N\t/g;
  $entry =~s/\t\t/\t\\N\t/g;
  $entry =~s/\t$/\t\\N/g;
  while ($entry=~/\r/){
    $entry =~s/\r//g;
  }
  print OUT "$entry";

}

close IN;
close OUT;


#fix create_keys

open (IN, $dir."/".$sql.$create_keys) or die "can't open $dir/".$sql.$create_keys."\n";
open (OUT, ">".$dir."/".$sql."new_".$create_keys) or die "can't open $dir/".$sql."new_".$create_keys."\n";


$/ = "\n\n";

while(<IN>){
  
   
  my $entry = $_;

  next if ($entry =~ /---------------------- Primary keys ------------------------/) ;
  next if ($entry =~ /---------------------- Foreign keys ------------------------/) ;


  $entry =~s/\[dbo\]\.//g;
  $entry =~s/\bGO\b/;/gi;
  $entry =~s/\[//g;
  $entry =~s/\]//g;
  $entry =~s/ WITH  FILLFACTOR = 90  ON PRIMARY//g;
 
  #separate out the alter_table clauses
  
  if ($entry =~/FOREIGN KEY/){
    my @rows = split/CONSTRAINT/, $entry;
    my $first_line = shift @rows;
    foreach my $constraint (@rows){
      $constraint = join " CONSTRAINT ", $first_line, $constraint;
      $constraint =~s/\),/\);\n/;
      print OUT $constraint."\n";
    }
  }else{
    print OUT "$entry";
  }

}

close IN;

#Add entry for extra table:
print OUT "\n\nALTER TABLE EnsemblWithdrawals ADD\n". 
	"\tCONSTRAINT PK_ensemblWithdrawals PRIMARY KEY  CLUSTERED\n". 
	"\t(\n".
	"\t\tccds_uid\n".
	"\t)\n". 
";\n";

print OUT "\n\nALTER TABLE EnsemblWithdrawals ADD\n". 
	 "\tCONSTRAINT  FK_EnsemblWithdrawals_CcdsUids FOREIGN KEY\n". 
	"\t(\n".
	"\t\tccds_uid\n".
	"\t) REFERENCES CcdsUids (\n".
	"\t\tccds_uid\n".
	");\n";

close OUT;
