#!/usr/local/ensembl/bin/perl -w

#script to fill out the analysis_description table
#expects text file with logic_name, description


use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

$! = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;

&GetOptions( 
            'dbhost=s'      => \$dbhost,
            'dbname=s'      => \$dbname,
            'dbuser=s'      => \$dbuser,
            'dbpass=s'      => \$dbpass,
            'dbport=s'      => \$dbport,
           );



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $dbhost,
   -user   => $dbuser,
   -dbname => $dbname,
   -pass   => $dbpass,
   -port   => $dbport,
  );

my $aa = $db->get_AnalysisAdaptor;
my $analyses = $aa->fetch_all;
my %descriptions;

while(<>){
  chomp;
  #print "Parsing ".$_."\n";
  my ($logic_name, $description) = split /\,/, $_;
  $description =~ s/^\s+//;
  $description =~ s/\s+$//;
  #print "Have logic_name ".$logic_name." and description ".$description."\n";
  $descriptions{$logic_name} = $description;
}

foreach my $analysis(@$analyses){
  my $logic_name = $analysis->logic_name;
  if(!defined($descriptions{$logic_name})){
    warn($logic_name." seems not to have an analysis in the table ");
    next;
  }
  $logic_name  =~ s/\_/ /;
  my $description = $descriptions{$logic_name};
  if(!$description){ $description  = "--"; }
  $analysis->description($description);
  $analysis->display_label($logic_name);
  print "setting ".$logic_name." / ".$description."\n";
  $aa->update($analysis) if($analysis->description);
}
