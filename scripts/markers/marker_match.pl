#!/usr/local/ensembl/bin/perl -w

# get a complete but non-redundant set of marker definitions

use strict;

my $infile1 = shift;
my $infile2 = shift;
my %marker = ();

#store file 1
open(IN,  "<$infile1") || die "cant open file $infile1\n" ; 
while(<IN>){ 
  chomp;
  my @line = split("\t", $_);
  push @{$marker{$line[0]}}, $_ ;   
}
close(IN);

#add file 2
open(IN,  "<$infile2")|| die "cant open file $infile1\n" ; ;
while(<IN>){
  chomp;
  my @line = split("\t", $_);
    push @{$marker{$line[0]}}, $_ ;  
}
close(IN);

# combine them 
foreach my $id (keys %marker){    
   my ( %names, %accs ) ; 
   my ($display_id, $lprim, $rprim, $dist, $name, $junk, $acc, $species) ;  

   for my $l (@{$marker{$id}}){    
     ($display_id, $lprim, $rprim, $dist, $name, $junk, $acc, $species) = split /\t/, $l;   

     # getting name unique 
     unless ($name=~m/-/){  
       if ($name=~m/;/) { 
         my @na = split/\;/,$name ;  
         @names{@na}=(); 
       } else { 
         $names{$name}=() ; 
       } 
     } 

     # getting acc unique 
     unless ($acc=~m/-/) {  
       if ($acc=~m/;/) { 
         my @ac = split/\;/,$acc ;  
         @accs{@ac}=(); 
       } else { 
         $accs{$acc}=() ; 
       } 
     }

   } 
   print "$display_id\t$lprim\t$rprim\t$dist\t";   
   unless (scalar(keys %names)==0) { 
     print join (";",keys %names) ; 
   } else { 
     print "\t-\t" ; 
   }
   print "\t$junk\t" ;  
   unless (scalar(keys %accs)==0) { 
     print join (";",keys %accs) ; 
   }else { 
     print "\t-\t" ; 
   }
   print "\t$species\n" ; 
}


__END__

87      AAAAACACAAGTTTCATACATCACA       AATGTAACTGTACCCTTCTGCATG        -       D9S1986 -       G07334;Z39132   Mus musculus
