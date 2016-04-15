#!/usr/bin/env perl

# This scripts fixes the supporting evidence links
# as part of the align feature tables optimisation.

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

$! = 1;

my $daf_outfile;
my $daf_infile;
my $sf_infile;
my $sf_outfile;

GetOptions(
            'outdaf=s'       => \$daf_outfile,
            'outsf=s'        => \$sf_outfile,
            'indaf=s'        => \$daf_infile,
            'insf=s'         => \$sf_infile,
           );

if (!$daf_outfile || !$sf_outfile ) {
  throw("Need 2 outfiles: -outdaf and -outsf");
}
if (!$daf_infile || !$sf_infile ) {
  throw("Need 2 infiles: -indaf and -insf");
}

print STDERR "DAFIN\t$daf_infile\nSFIN\t$sf_infile\nDAFOUT\t$daf_outfile\nSFOUT\t$sf_outfile\n";

# open infile
open(DAFIN, "<$daf_infile") or die ("Can't read $daf_infile $! \n");
open(SFIN,  "<$sf_infile")  or die ("Can't read $sf_infile $! \n");

# open outfile
open(DAFOUT,">$daf_outfile") or die ("couldn't open file ".$daf_outfile." $!"); 
open(SFOUT, ">$sf_outfile")  or die ("couldn't open file ".$sf_outfile." $!");

# read in daf list;
my %daf_positions;
my $count = 0;

while (<DAFIN>) {
  my $line = $_;
  
  my ($old_daf_id) = $line =~ /^(\d+)/;
  $count++;
  $line =~ s/^\d+/$count/;
  print DAFOUT $line;
  # key = old daf id
  # value = new daf id
  $daf_positions{$old_daf_id} = $count; 
}
close (DAFIN);
close (DAFOUT);

# read in the tsf or sf file
# link to the daf files
# and write out to the new file

while (<SFIN>) {
  my $line = $_;
  chomp $line;
  my @fields = split (/\t/,$line);

  if (scalar(@fields) != 3) {
    throw("Must have 3 fields for -infile: $line");
  }

  # define variables
  my $supporting_feature_id = $fields[0];
  my $feature_type = $fields[1];
  my $feature_id = $fields[2];

  if (!exists $daf_positions{$feature_id}|| !defined $daf_positions{$feature_id}) {
    throw("Cannot find daf matching to supporting feature: $feature_id $feature_type $supporting_feature_id");
  }

  print SFOUT "$supporting_feature_id\t$feature_type\t".$daf_positions{$feature_id}."\n"; 

} # close infile
close (SFOUT);
close (SFIN);
