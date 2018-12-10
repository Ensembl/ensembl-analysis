#!usr/bin/perl

use strict;
use warnings;

open (T_ID_LIST, "<DataFiles/transcript_IDs_all.txt");
open (T_RS_MATCH_LIST, "<DataFiles/CARS_Refseq_TIDs.txt");

my $t_line = "";
my $c1 = 0;
my $c2 = 0;
my $c3 = 0;

while (<T_ID_LIST>){
    $t_line = $t_line." | ".$_;
    $c1++;
}

close (T_ID_LIST);

while (<T_RS_MATCH_LIST>){
    if (index($t_line, $_) != -1){
        $c2++;
    }
    $c3++;
}

close (T_RS_MATCH_LIST);

print "Total transcript IDs to be compared : $c1\n";
print "Total transcript IDs in RefSeq file : $c3\n";
print "Total transcript IDs matched : $c2\n";

exit;