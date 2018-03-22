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

=pod

=head1 NAME

load_Ditags.pl

=head1 DESCRIPTION

 Script to set up the Ditag analysis:

 - reads & checks config files Ditag.pm & Databases.pm
 - parses file with ditag FASTA entries
 - stores them to the ditag table
 - generates input ids for the DitagAlign analysis

 It writes the SQL inserts to a file first which is then
 executed with the mysql client to achieve a significant
 speedup.

=head1 SYNOPSIS

 -type   <> type of ditag to load
 -option <> (optional) run only <A> loading of ditags
                             or <B> creation of input-ids
 -write     (optional) write-flag
 -delete/nodelete (optional) (not) delete the sql files after storage
            if not given user will be prompted
 -help      (optional) this text

 If you are confident, just call it like this to get ready to run:
 perl load_Ditags.pl -type ZZ11 -write -delete
 Make sure the type matches the type in the Config::ExonerateTags.pm file.

=head1 CONTACT

http://lists.ensembl.org/mailman/listinfo/dev

=cut

use warnings ;
use strict;
use Bio::EnsEMBL::Map::Ditag;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::ExonerateTags;
use Bio::EnsEMBL::Analysis::Config::Databases;
use Getopt::Long;

my ($tagtype, $help);
my $write      = 0;
my $delete     = undef;
my $option     = "";
my $logic_name = "SubmitDitag";
my $analysis;
my %config;

GetOptions(
	   'type=s'    => \$tagtype,
	   'option:s'  => \$option,
	   'write!'    => \$write,
	   'delete!'   => \$delete,
           'logic_name=s'=> \$logic_name,
           'analysis=s'  => \$analysis,
	   'help!'     => \$help,
	  );

if($help or !$tagtype){
   exec('perldoc', $0);
   exit 1;
}

my $db = setup();
print STDERR "\nReading from Bio::EnsEMBL::Analysis::Config::Databases && ".
             "Bio::EnsEMBL::Analysis::Config::ExonerateTags.\nLooking at database " .
               $$DATABASES{'REFERENCE_DB'}{'-dbname'}.":".
               $$DATABASES{'REFERENCE_DB'}{'-host'}.":".
               $$DATABASES{'REFERENCE_DB'}{'-port'}."\n\n";

if(!$option or $option eq "A"){
  save_ditags($db, $tagtype, $write, $delete);
}
if((!$option and $write) or $option eq "B"){
  make_input_ids($db, $write, $delete);
}


 #################

sub setup {

  #READ CONFIG FILE VARIABLES INTO LOCAL HASH
  #default variables
  my $default_entry = $DITAG_CONFIG->{DEFAULT};
  foreach my $config_var (keys %{$default_entry}) {
    $config{$config_var} = $default_entry->{$config_var};
  }
  $default_entry = $DITAG_CONFIG->{$analysis};
  foreach my $config_var (keys %{$default_entry}) {
    $config{$config_var} = $default_entry->{$config_var};
  }
  #specific variables
  if (exists $DITAG_CONFIG->{$tagtype}) {
    my $entry = $DITAG_CONFIG->{$tagtype};
    foreach my $config_var (keys %{$entry}) {
      $config{$config_var} = $entry->{$config_var};
    }
  }
  if ( !$$DATABASES{'REFERENCE_DB'}{'-dbname'} or
       !$$DATABASES{'REFERENCE_DB'}{'-host'} or
       !$$DATABASES{'REFERENCE_DB'}{'-port'} or
       !$$DATABASES{'REFERENCE_DB'}{'-user'} or
       !$$DATABASES{'REFERENCE_DB'}{'-pass'} ){ 
    throw "\nDatabase parameters missing from file ".$$DATABASES{'REFERENCE_DB'}{'-host'}.
      "Bio::EnsEMBL:Analysis::Config::GeneBuild::Databases\n";
  }

  #CHECK ALL VARIABLES NEEDED
  my @configvars = qw( QUERYFILES GENOMICSEQS IIDREGEXP BATCHSIZE TMPDIR PROGRAM OPTIONS );
  foreach my $configvar (@configvars){
    my $configvarref = $config{$configvar};
    if(!$configvarref){
      throw "\nAnalysis parameters $configvar missing from file ".
	    "Bio::EnsEMBL:Analysis::Config::ExonerateTags\n";
    }
    if($configvar =~ m/DIR/){
      if(!-e $configvarref){
	throw "\nCouldn t find/access file or directory ".
	      "from file Bio::EnsEMBL:Analysis::Config::ExonerateTags\n";
      }
    }
  }

  #CONNECT TO PIPELINE DATABASE
  my $pipedb = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(
							    -host    => $$DATABASES{'REFERENCE_DB'}{'-host'},
							    -port    => $$DATABASES{'REFERENCE_DB'}{'-port'},
							    -user    => $$DATABASES{'REFERENCE_DB'}{'-user'},
							    -pass    => $$DATABASES{'REFERENCE_DB'}{'-pass'},
							    -dbname  => $$DATABASES{'REFERENCE_DB'}{'-dbname'}
							   )
					  or throw "cant connect to pipedatabase.";
  return($pipedb);
}


sub save_ditags {
  my ($db, $tagtype, $write, $delete) = @_;

  my ($probe_name, $sequence);

  #when importing mappings from other sources
  #there might be information about the tag number
  #within the header line
  my $expect_count = $config{TAGCOUNT};

  my $ditagcount   = 0;
  my $count        = 0;
  my $foo          = 0;
  my $ans          = "n";
  my %seqs;
  my $ditag_adaptor = $db->get_ditagAdaptor
       or throw("Couldn t get DitagAdaptor. Check analysis tables.");

  my $tmpfile = $config{TMPDIR}.'/ditag_inserts.'.$tagtype.'.sql';
  my $queryfile = $config{QUERYFILES}{$tagtype};

  open( DITAGS, "<" . $queryfile )
    or throw( "Couldn t open ditag file " . $queryfile );
  open( TOFILE, ">" . $tmpfile )
    or throw( "Couldn t open output file " . $tmpfile );

  while (<DITAGS>) {
    chomp;
    if ($_ =~ /^>/) {

      $_ =~ /$config{IIDREGEXP}/;

      if(defined $1) {
        $probe_name = $1;
        if(defined $2) {
	  if($expect_count){
	    $ditagcount = $2;
	  }
	  else{
	    $probe_name .= "-".$2;
	  }
        }
      }
      else {
        throw ("header ".$_." not parsed correctly.");
      }
    }
    else {
      $sequence = $_;
      if ( !$probe_name ) {
        warning("sequence without id? [$sequence] ");
      }
      elsif($sequence =~ /[^ACTGN]/ig){
        warning("strange characters in sequence! [$sequence] ");
      }
      else {
	if(!defined($seqs{$sequence})){
	  my %tag;
	  $tag{'probe_name'} = $probe_name;
	  $tag{'tagtype'}    = $tagtype;
	  $tag{'count'}      = ($ditagcount or 1);
	  $seqs{$sequence}   = \%tag;
	}
	else{
	  $seqs{$sequence}->{'count'}++;
	}
      }
      $probe_name = '';
    }

  }

  foreach my $tag_sequence (keys %seqs){
    print TOFILE $ditag_adaptor->print_creation(
						$seqs{$tag_sequence}->{'probe_name'},
						$seqs{$tag_sequence}->{'tagtype'},
						$seqs{$tag_sequence}->{'count'},
						$tag_sequence
					       );
    $count++;
  }

  close(TOFILE) or throw("cant close file ".$tmpfile);
  close(DITAGS) or throw("cant close file ".$queryfile);
  print STDERR "\nSQL INSERTS FOR DITAGS WRITTEN TO $tmpfile [$count].\n";

  if($write){
    #LOAD TO DATABASE 
    my $cmd = "mysql -h $$DATABASES{'REFERENCE_DB'}{'-host'} ".
               " -P$$DATABASES{'REFERENCE_DB'}{'-port'} ".
               " -u$$DATABASES{'REFERENCE_DB'}{'-user'} ".
               " -p$$DATABASES{'REFERENCE_DB'}{'-pass'} ".
               " -D$$DATABASES{'REFERENCE_DB'}{'-dbname'} ".
               " < $tmpfile " ;

    if(system($cmd)) { 
      throw("couldn t load file $tmpfile to database $$DATABASES{'REFERENCE_DB'}{'-dbname'}");
    }
    print STDERR "DITAGS STORED IN ".  $$DATABASES{'REFERENCE_DB'}{'-dbname'}.":".
               $$DATABASES{'REFERENCE_DB'}{'-host'}.":".
               $$DATABASES{'REFERENCE_DB'}{'-port'}."\n";

    if(!defined $delete){
      print STDERR "REMOVE FILE? ";
      $ans = <>;
      chomp $ans;
    }
    elsif($delete){
      $ans = "y";
    }
    if($ans =~ /^y/i){
      system("rm $tmpfile");
    }
  }
  else{
    print STDERR "DITAGS NOT STORED IN DATABASE YET.\n"
  }

}


sub make_input_ids {
  my ($db, $write, $delete) = @_;

  my $ans          = "n";
  my $analyisObj   = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
  my $ditagadaptor = $db->get_ditagAdaptor() or throw("no ditag adaptor.");
  my $ana_id       = $analyisObj->dbID;
  my $input_type   = $db->get_AnalysisAdaptor->fetch_analysis_input_id_type($analyisObj);

  my $number_of_ditags = scalar @{ $ditagadaptor->fetch_all_by_type($tagtype) };
  if(!$number_of_ditags){ throw "Couldnt find any ditags of type \"$tagtype\" in db."; }
     #desired 1 000 * 2 000 = 2 000 000
  my $number_of_chunks = ($number_of_ditags / $config{BATCHSIZE}) | 1;
  my $tmpfile = $config{TMPDIR}.'/ditag_inputids.'.$tagtype.'.sql';
  print STDERR "Number of Ditags: $number_of_ditags, Chunks: $number_of_chunks.\n";

  open(SQL, ">$tmpfile") or throw("cant create temp sql file $tmpfile.");

  for(my $i=1; $i<=$number_of_chunks; $i++){
    my $iid = "ditag.".$tagtype.".".$i;
    print SQL "INSERT into input_id_analysis ".
          "values ('$iid', '$input_type', $ana_id, now(), '', '', 0);\n";
  }
  close(SQL) or throw("cant close file ".$tmpfile);
  print STDERR "SQL INSERTS FOR INPUT IDS WRITTEN TO $tmpfile.\n";

  if($write){
    #LOAD TO DATABASE
    my $cmd = "mysql -h $$DATABASES{'REFERENCE_DB'}{'-host'} ".
               " -P$$DATABASES{'REFERENCE_DB'}{'-port'} ".
               " -u$$DATABASES{'REFERENCE_DB'}{'-user'} ".
               " -p$$DATABASES{'REFERENCE_DB'}{'-pass'} ".
               " -D$$DATABASES{'REFERENCE_DB'}{'-dbname'} ".
               " < $tmpfile " ; 
    if(system($cmd)) { 
      throw("couldn t load file $tmpfile to database $$DATABASES{'REFERENCE_DB'}{'-dbname'}");
    }
    print STDERR "INPUT IDS CREATED IN ".
               $$DATABASES{'REFERENCE_DB'}{'-dbname'}.":".
               $$DATABASES{'REFERENCE_DB'}{'-host'}.":".
               $$DATABASES{'REFERENCE_DB'}{'-port'}."\n";
    if(!defined $delete){
      print STDERR "REMOVE FILE? ";
      $ans = <>;
      chomp $ans;
    }
    elsif($delete){
      $ans = "y";
    }
    if($ans =~ /^y/i){
      system("rm $tmpfile");
    }
  }
  else{
    print STDERR "NOT STORED IN DATABASE YET.\n"
  }
}

1;
