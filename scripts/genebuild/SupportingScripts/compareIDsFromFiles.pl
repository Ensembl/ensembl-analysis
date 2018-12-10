#!usr/bin/perl

use strict;
use warnings;

open (F1, '<DataFiles/t_ids_w_intropolis');
open (F2, '<DataFiles/t_ids_wo_intropolis.txt');

open (O1, '>DataFiles/match.txt');
open (O2, '>DataFiles/dont_match.txt');

my $t1;
my $t2;

while (!eof(F1)&& !eof(F2)){
    $t1 = <F1>;
    $t2 = <F2>;
    chomp $t1;
    chomp $t2;
    
    
    if ($t1 eq $t2){
        print O1 $t1,"\n"; 
    }
    else{
        print O2 $t1,"\t",$t2,"\n"
    }
}

close (F1);
close (F2);
close (O1);
close (O2);

exit;