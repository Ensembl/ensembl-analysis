#!/usr/bin/env perl
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

#
# script to add external_db_ids to the hit_names
# in protein- & dna-align_feature tables
#
# -example:
#  perl assign_external_db_ids.pl -conf fugu_dafs.cfg -feature_type [dna, protein]
#       -dumpdir /fugu_output
# -use option -nodump if the table was successfully dumped before and should be re-used
# -every analysis / regex must be defined in the config file
# -use LSF to run!
# -if successful, do a 'delete from table' & 'load data local infile "/absolutepath/table.fixed" into table dna_align_feature;'
#


# The full_path_to_uniprot file looks like this and holds all of uniprot
#P83011.1        STD     uniprot_sprot_vertebrates.dat
#P83010.1        STD     uniprot_sprot_vertebrates.dat
#Q5ZLQ6.1        STD     uniprot_sprot_vertebrates.dat
#Q5XGC8.1        STD     uniprot_sprot_vertebrates.dat
# grep -w PRE ../Ensembl/uniprot_2013_01.entry_loc |wc -l
#29266939
# grep -w STD ../Ensembl/uniprot_2013_01.entry_loc | wc -l
#538849



use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);


sub trim($);

my $host         = '';
my $user         = '';
my $dbname       = '';
my $port         = 3306;
my $pass         = '';
my $conf_file    = undef;
my $feature_type = undef;
my $dump_dir     = "";
my $dump_option  = 1; #do the fixing in the file in any case
my $nodump       = 0; #use existing dump
my $VERBOSE      = 0;
my $update_only_null_rows = 0 ; # updates only rows which have a null-value in the external db id field 
my $full_path_to_uniprot; # ../Ensembl/uniprot_2013_05/entry_loc 

# database containing the latest version of the external_db table
my $masterhost = "production_host" ;
my $masteruser = "read_only_user" ;
my $masterdbname = "production_dbname" ;
my $masterport = "3306" ;
my $masterpass = "" ;
$| = 1;

&GetOptions(
  'dbhost|host|h:s'        => \$host,
  'dbuser|user|u:s'        => \$user,
  'dbpass|pass|p:s'        => \$pass,
  'dbname|db|D:s'      => \$dbname,
  'dbport|port|P:n'        => \$port,
  'conf_file:s'   => \$conf_file, 
  'feature_type:s'=> \$feature_type,
  'uniprot_filename:s'=> \$full_path_to_uniprot,
  'dumpdir:s'     => \$dump_dir,
  'dump_option!'     => \$dump_dir,
  'nodump!'       => \$nodump, 
  'update_only_null_rows!' => \$update_only_null_rows, 
  'masterhost:s'       => \$masterhost,
  'masteruser:s'       => \$masteruser,
  'masterdbname:s'       => \$masterdbname,
  'masterport:n'       => \$masterport,
  'masterpass:s'       => \$masterpass,
);


#
# do some checks on the config before proceeding
#
if(!$host || !$dbname || !$conf_file || !$feature_type || !$dump_dir){
  die "Please define at least host, dbname, conf_file, feature_type, uniprot_filename and dumpdir.\n";
}
if ($feature_type ne "protein" && $feature_type ne "dna") {
  die "Feature type must be either 'protein' or 'dna' (it's $feature_type)";
} 
unless ( -e $dump_dir ) {
   die "HEY ! Your dump dir $dump_dir does not exist.\nTRY :    mkdir -p $dump_dir  \n" ;
}

if ($feature_type eq 'protein') {
  die('You need to provide a UniProtKB accession list like /some/path/uniprot_2018_01/entry_loc') unless ($full_path_to_uniprot);
  die("HEY ! Your uniprot_filename $full_path_to_uniprot does not exist.\n") unless ( -e $full_path_to_uniprot );
}

#
# connect to core database
#
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -dbname => $dbname,
  -pass   => $pass,
) or die "can not connect to $dbname@$host\n";


#Read config
# This is the output_config_file produced by running test_regexes.pl 
my $sections = read_config($conf_file);


#Get list of external dbs
# looks in the production database
my %extdb_names_hash = %{get_extdb_names_hash()};

# Looks in your core db
my %logic_names_hash = %{get_logic_names_hash($db)};



# First check if there's only one external_db_id for each analysis and whether or not 'Uniprot/Uncertain' is one of them
#   If only one analysis for each but not 'Uniprot/Uncertain' then we can just quickly set the external_db_ids with the table in place
#   else we need to dump it to avoid having to use the regexes on the db which are slow
# Also check that the external_db_ids and analyses exist here to stop problems later

my $only_single = 1;
my $have_uniprot_uncertain = 0;
foreach my $section (@$sections) {
  my @exdbs = @{$section->exdbs};

  if (!exists($logic_names_hash{$section->analysis})) {
    print STDERR  "Couldn't find analysis " . $section->analysis . " in db\n";
  }

  foreach my $exdb_name (@exdbs) {
    if (!exists($extdb_names_hash{$exdb_name}) && $exdb_name ne "Uniprot/Uncertain") {
      die "Couldn't find exdb name " . $exdb_name . " in db\n";
    }
    if ($exdb_name eq "Uniprot/Uncertain") {
      $have_uniprot_uncertain = 1;
    }
  }
  if ($section->default) {
    if (!exists($extdb_names_hash{$section->default})) {
      die "Couldn't find exdb name " . $section->default . " in db\n";
    }
  }

  if (scalar(@exdbs) > 1) {
    print "Have multiple exdbs\n";
    $only_single = 0;
  }
}


# a backup in case things go wrong:
my $dump_file = $dump_dir."/table_dumps_" .$dbname . ".". $feature_type . ".sql";
# and the file you hope to load in to your db:
my $load_file = $dump_dir."/table_dumps_" .$dbname . ".". $feature_type . ".fixed";



if (($only_single && !$have_uniprot_uncertain) and !$dump_option) {
  # For each section make the sql statement to change the 
  #"update " . $type . " _align_feature set external_db_id = XXX where analysis_id = XXX";
  print "Only single external db\n\n";
  open FPSQL,">$load_file";
  foreach my $section (@{$sections}) {
    my $exdb_name = $section->exdbs->[0];
    my $query = "update " . $feature_type . "_align_feature set external_db_id=" .
      $extdb_names_hash{$exdb_name} . " where analysis_id=" . $logic_names_hash{$section->analysis}->dbID . ";";
    print "$query\n";
    print FPSQL "$query\n";
    #my $sth = $db->dbc->prepare($query);
    #$sth->execute;
    #$sth->finish;
  }
  close FPSQL;
} else {
  print "More than one external_db found, and/or Uniprot/Uncertain.\n\n"; 
  #my ($uniprot_acc_to_class, $uniprot_acc_no_ver_to_class);
  my $uniprot_hash;

  if(!$nodump){
    if ( -e $dump_file ) {
      die "HEY ! Your dump_file $dump_file exists; please delete \n" ;
    }


    # dump table
    print "Dumping table.\n";  
    print " you can still stop ...................\n"  ; 
    sleep(5) ; 
    print " ..................dumping...............\n"; 

    # we can use less memory if we do this one analysis at a time
    foreach my $section (@$sections) {
      print "  doing ".$logic_names_hash{$section->analysis}->logic_name."\n";
      my $analysis_dbID = $logic_names_hash{$section->analysis}->dbID; 
      system("mysql -q -h$host -uensro -P$port -D$dbname -N -B -e 'select * from " . $feature_type . "_align_feature where analysis_id = $analysis_dbID' | sed -e 's/NULL/\\\\N/g' >> $dump_file");
    }
  }

  if ($have_uniprot_uncertain) { 
    print "Loading stuff for uniprot... this will take > 2 GB of memory\n" ; 
    #($uniprot_acc_to_class, $uniprot_acc_no_ver_to_class) = load_uniprot_mapping_hash($full_path_to_uniprot);
    $uniprot_hash = load_uniprot_mapping_hash($full_path_to_uniprot);
  }
  print "Now assigning external db ids...\n";

  my $anal_col;
  my $hit_name_col;
  my $exdb_col;
  if ($feature_type eq 'protein') {
    $anal_col = 8;
    $hit_name_col = 7;
    $exdb_col = 13;
  } else {
    # dna_align_feature table has a hit_strand column
    $anal_col = 9;
    $hit_name_col = 8;
    $exdb_col = 14;
  }

  my %anal_id_hash = %{get_anal_id_hash($db)};

  my %logic_name_to_section_hash;
  foreach my $section (@$sections) {
    $logic_name_to_section_hash{$section->analysis} = $section;
  }

  if ( $nodump && ! -e $dump_file ) {  
    die " You used option -nodump but there is no dump file i can read. create a dump first by not using -nodump option \n" ; 
  } 
  
  # go through assigning the ext db ids
  open FPDUMP,"<$dump_file" || die " cant open file $dump_file" ; 
  open FPLOAD,">$load_file" || die " cant read $load_file \n";  

  DUMP: while (<FPDUMP>) {
    my @words = split(/\t/, $_);
    my $logic_name = $anal_id_hash{$words[$anal_col]};

    my $section = $logic_name_to_section_hash{$logic_name}; 
    unless ($section) {  
      #print "skipping section $logic_name\n" ; 
      next DUMP ; 
    } 
    my $hit_name = $words[$hit_name_col];
    my $found = 0;
    my $exdb = undef; 


    if ( $update_only_null_rows ) { 
      unless (  $words[$exdb_col]=~m/\\N/ || int($words[$exdb_col])==0 || $words[$exdb_col] eq 'NULL') {
        # skipping as external db id is already set   
      # print "skipping $words[$exdb_col]\n" ;
       next DUMP;
      }
    } 

    #print "testing regexes .....\n" ; 
   REGEX: foreach my $regex (@{$section->get_regexes}) {
    #  print "Testing " . $regex->{regex} . " hit_name = " . $words[$hit_name_col] . " logic_name = " . $logic_name . "\n";
       #print $hit_name . "\n" ;  
      if ($hit_name =~ /$regex->{regex}/) {
        #print "Match to " . $regex->{regex} . "\n";

        $exdb = $regex->{exdb};

        if ($exdb eq "Uniprot/Uncertain") {
          # remove version from hit nme because the uniprot hash doesn't use versions in the key
          $hit_name =~s/\.\d+//;
          if (exists $uniprot_hash->{$hit_name}) {
            if ($uniprot_hash->{$hit_name} == 2) {
              $exdb = "Uniprot/SPTREMBL";
            } elsif ($uniprot_hash->{$hit_name} == 1) {
              $exdb = "Uniprot/SWISSPROT";
            }
            #print "$hit_name $uniprot_acc_to_class->{$hit_name}\n"
          } else {
            print "Couldn't find $hit_name for analysis $logic_name : setting to Uniprot/SPTREMBL\n";
            $exdb = "Uniprot/SPTREMBL";
          }
        }
        $found = 1;
        last REGEX;
      }else {  
        #print "hit name does not match regex \n" ; 
      } 
    }

    if ($found) { 
    #  print "found workd\n" ; 
      $words[$exdb_col] = $extdb_names_hash{$exdb};
    } elsif (!$found && $section->default) {
    #  print "found word section default\n" ; 
      $words[$exdb_col] = $extdb_names_hash{$section->default};
    } else {
    #  print "Didn't find matching regex for " .  $hit_name  . " in analysis " . $logic_name . "\n";
      if($VERBOSE){
	foreach my $regex (@{$section->get_regexes}) {
	  print "Tested " . $regex->{regex} . " hit_name = " . $words[$hit_name_col] . " logic_name = " . $logic_name . "\n";
	}
      }
    }
    print FPLOAD join("\t",@words);
  }
  print "closing file \n" ; 
  close FPDUMP;
  close FPLOAD;

  # drop table
  # load table 
  
  $uniprot_hash = undef;
  #$uniprot_acc_to_class = undef;
  #$uniprot_acc_no_ver_to_class = undef;
}

print "Done\n";
print "\n\n\nNow you should upload the data into your database using the following commands: \n" ; 

my $backup_table; 
if ( $feature_type =~m/protein/) {  
  $backup_table = "paf_bak" ; 
}else {  
  $backup_table = "daf_bak" ; 
} 
print "mysql -h $host -u $user -p$pass -P$port -D $dbname -e\'show tables like \"$backup_table\"\'\n"; 
print "mysql -h $host -u $user -p$pass -P$port -D $dbname -e\"rename table ".$feature_type."_align_feature to $backup_table \"; \n" ; 
print "mysql -h $host -u $user -p$pass -P$port -D $dbname -e\"create table ".$feature_type."_align_feature like $backup_table \"; \n" ; 
print "mysql -h $host -u $user -p$pass -P$port -D $dbname -e\'load data local infile \"$load_file\" into table ".$feature_type."_align_feature \'; \n" ; 
# -if successful, do a 'delete from table' & 'load data local infile "/absolutepath/table.fixed" into table dna_align_feature;'
exit;




sub get_extdb_names_hash {

  my $masterdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $masterhost,
    -user   => $masteruser,
    -port   => $masterport,
    -dbname => $masterdbname,
    -pass   => $masterpass,
  ) or die "can not connect to $masterdbname@$masterhost\n";

  my %namehash;

  my $query = "select db_name, external_db_id from external_db " ;
  my $sth = $masterdb->dbc->prepare($query) ;
  $sth->execute() ;

  while ( my $arrayref = $sth->fetchrow_arrayref() ) {
    my ( $name, $dbid ) = @$arrayref ;
    $namehash{$name} = $dbid ;
  }

  $sth->finish ;

  return \%namehash;
}

sub get_logic_names_hash {
  my $db = shift;

  my $aa = $db->get_AnalysisAdaptor;

  my @anals = @{$aa->fetch_all};

  my %analhash;
  foreach my $anal (@anals) {
    $analhash{$anal->logic_name} = $anal;
  }

  return \%analhash;
}

sub get_anal_id_hash {
  my $db = shift;

  my $aa = $db->get_AnalysisAdaptor;

  my @anals = @{$aa->fetch_all};

  my %analhash;
  foreach my $anal (@anals) {
    $analhash{$anal->dbID} = $anal->logic_name;
  }

  return \%analhash;
}

sub read_config {
  my $conf_file = shift;

  open (FP_CONF,"<$conf_file") || die " cant read config file ....\n" ;

  my @sections;
  my $cur_section;
  while (<FP_CONF>) {
    chomp;
    $_ = trim($_);
    next if (/^#/);
    next if (/^$/);

    my ($key,$value) = split /\t/,$_,2;

    print "Key = $key value = $value\n";
    if ($key eq 'Analysis') {
      $cur_section = ConfigSection->new(-ANALYSIS => $value);
      push @sections,$cur_section;
    } elsif ($key eq 'Regex') {
      my ($regex,$exdb) = split /\t/,$value;
      print $regex . " " . $exdb ."\n";
      $cur_section->add_regex($regex,$exdb);
    } elsif ($key eq 'Default') {
      $cur_section->default($value);
    } else {
      die "Error: unknown line type in config file. Line = " . $_;
    }
  }
  close FP_CONF;
  return \@sections;
}


sub trim($) {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}


sub load_uniprot_mapping_hash {
  my ($file_name) = @_;

  my %uniprot;
  my %dbs = ( PRE => 2, STD => 1);

  my $count = 0;
  print "Loading $file_name\n";
  open FP_UPDATE, "<$file_name" or die "can not open file $file_name.\n";
  while (<FP_UPDATE>) {
    chomp;
    my ($acc,$db) = split('\s+', $_, 2);
	next if ($db =~ /PRE/);
    # remove version
    $acc =~ s/\.\d+//;
    $uniprot{$acc} = $dbs{$db};
  }
  close FP_UPDATE;
  print "Finished loading ".(scalar(keys %uniprot))." uniprot accessions into memory\n";

  return (\%uniprot);
}


package ConfigSection;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);

sub new {
  my ($caller,@args) = @_;
  my $class = ref($caller) || $caller;


  my ($analysis) = rearrange([qw(ANALYSIS)], @args);

  return bless({'analysis' => $analysis}, $class);

}

sub analysis {
  my $self = shift;

  return $self->{analysis};
}

sub add_regex {
  my ($self,$regex,$exdb) = @_;

  if (!defined($exdb)) {
    die "Must pass both regex and exdb to add_regex";
  }
 
  if (!exists $self->{regexes}) {
    $self->{regexes} = [];
  }
  my %hash;
  $regex = trim($regex);
  $regex =~ s/^'//;
  $regex =~ s/'$//;
  print "After trimming ' regex is $regex \n" if $VERBOSE;
  $hash{regex} = $regex;
  $hash{exdb} = $exdb;

  push @{$self->{regexes}},\%hash;
}

sub get_regexes {
  my $self = shift;

  return $self->{regexes};
}

sub exdbs {
  my $self = shift;

  my %exdbs;
  foreach my $regex (@{$self->get_regexes}) {
    $exdbs{$regex->{exdb}} = 1;
  }

  return [ keys %exdbs ];
}

sub default {
  my $self = shift;

  if (@_ && $self->{default}) {
    die "Default already defined for " . $self->analysis;
  }

  $self->{default} = shift if (@_);

  return $self->{default};
}

sub trim($) {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}
1;
