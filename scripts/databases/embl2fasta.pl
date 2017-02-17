use warnings;
use strict;
use feature 'say';


my $infile = $ARGV[0];
my $outfile = $infile;
$outfile =~ s/\.dat/\.fasta/;

my $state = 0;
my $current_accession = "";
my $current_version = "";
my $current_seq = "";
open(IN,$infile);
open(OUT,">$outfile");
while(<IN>) {
  my $line = $_;
  if($state == 0 && $line =~ /^AC +([a-zA-Z\d]+)\;/) {
    $current_accession = $1;
    $state = 1;
  }

  elsif($state == 1 && $line =~ /^DT.+sequence version (\d+)/) {
    $current_version = $1;
    $state = 2;
  }

  elsif($state == 1 && $line =~ /^SQ +SEQUENCE +/) {
    die "Failed to parse version for ".$current_accession."\nExiting";
  }

  elsif($state == 2 && $line =~ /^SQ +SEQUENCE +/) {
    $state = 3;
  }

  elsif($state == 3 && !($line =~ /^\/\//)) {
    $current_seq .= $line;
  }

  elsif($state == 3 && $line =~ /^\/\//) {
    $current_seq =~ s/ //g;
    print OUT '>'.$current_accession.".".$current_version."\n".$current_seq;
    $state = 0;
    $current_accession = '';
    $current_version = '';
    $current_seq = '';
  }

}
close OUT;
close IN;

exit;
