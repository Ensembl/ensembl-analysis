#!usr/bin/perl

package Utilities;
#use exporter;
use strict; 
use warnings;
use List::Util qw ( min max );

#@EXPORT

sub sumNums {
    my @array = @_;
    my $total = 0;
    foreach my $num ( @array ){
        #print ($num, "\n");
        $total += $num;
    }
    return $total;
}

sub normalizeArray {
    #print "\nNormalized values for above gene: \n";
    my @array = @_;
    my $max = max @array;
    my @normArray;
    for (@array){
        if ($_!=0){
            push (@normArray, $_/$max);
        }        
    }
    # @normArray = sort {$b<=>$a} @normArray;
#     my $i=0;
#     for ( @normArray ){
#         if ($_ == 1){
#             print ++$i,"\n";
#         }
#         #print "$_, ";
#}
    return @normArray;
}

1;