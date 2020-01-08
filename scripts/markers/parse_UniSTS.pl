#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

use warnings ;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;


#  SETTING THE DEFAULT VALUES FOR DB-CONNECTION ETC
my %opt = (
           file => 'UniSTS.sts',
          );

&GetOptions( 
             \%opt, 
            'file=s',  
            'species=s',
            'outfile=s',  
           );


print "reading file $opt{file}\n" ;  

my (%unists,%id_org_cnt) ;   
my @species_entries ; 
open(UNI, $opt{file}) || die " Can't read file $opt{file}\n" ;  

FILE:
while (my $line = <UNI>) {   
    next if $line=~m/^\#/; 
    chomp($line) ;  

    my ( $id, $lprim, $rprim, $len, $name, $chr, $acc, $org )= split /\t/,$line ;    

    # consistency-check 
    if ( exists $unists{$id}{$org} ) {    
      if ( $unists{$id}{$org}{entry} eq $line ) {
        print "2 redundant entries in file $opt{file} for $id - $org - I skip this line \n" ; 
        next FILE ;  
      }
      print "ENTRY: $line\n" ;  
      print "ENTRY: $opt{$id}{$org}\n" ; 
      die("File has 2 entrie for the same marker in the same species : $id --> $org\n") ; 
    }  
    $unists{$id}{$org}{entry} = $line ;    
    $unists{$id}{$org}{len} = $len;    
    $unists{$id}{$org}{name} = $name ; 
    $unists{$id}{$org}{rprim} = $rprim; 
    $unists{$id}{$org}{lprim} = $lprim ; 
    $unists{$id}{$org}{chr} = $chr ; 
    $unists{$id}{$org}{acc} = $acc; 
    $unists{$id}{$org}{org} = $org;   

    $id_org_cnt{$id}{cnt}++; 
    push @species_entries, $id if ( $org =~m/$opt{species}/) ; 
}
close(UNI) ;  

open(O, "> $opt{outfile}") || die("Can't write to outfile $opt{outfile}\n"); 

print scalar(@species_entries) . " entries for $opt{species} found\n" ; 

MARKER: 
for my $id ( @species_entries ) {   
  my ( $len, $name, $len_org , $name_org ) ; 

  if ( $id_org_cnt{$id}{cnt} == 1 ) {  
    # unique marker only for dog   
    print O "$unists{$id}{$opt{species}}{entry}\n" ; 
    next MARKER; 
  } 
  # try to stick info out of file together ...   
  #
  #  get ACC and LEN attribute and check if rprim and lprim are the same 
   my %uni_orgs = %{$unists{$id}};  

 
   for my $org ( keys %uni_orgs ) {  

     if ( ($uni_orgs{$opt{species}}{len}=~m/\d/) 
          && ($uni_orgs{$opt{species}}{name} !~ m/-/) ) {  
          # entry is complete for $opt{species} - we don't need to do anything  
          print O "$unists{$id}{$opt{species}}{entry}\n" ;  
          next MARKER ;  
     } 
   
     # 
     # get length and name from same species 
     #   
     if ( $uni_orgs{$org}{len}=~m/\d/  &&  $uni_orgs{$org}{name} ne '-' ) {   
        $len_org = $org ; 
        $len = $uni_orgs{$org}{len};
        $name_org = $org ; 
        $name = $uni_orgs{$org}{name}; 
        #print " name $uni_orgs{$org}{name}\n" ;   
     }

     
     if ( $uni_orgs{$org}{len}=~m/\d/  && !$len ) {  
        #print " len $uni_orgs{$org}{len}\n" ; 
        $len_org = $org ; 
        $len = $uni_orgs{$org}{len};
     }  
    
     #
     # get length if length is not defiend for $opt{species}  
     #   
     
     if ( $uni_orgs{$org}{name} ne '-' && !$name) {  
        #print " name $uni_orgs{$org}{name}\n" ;   
        $name_org = $org ; 
        $name = $uni_orgs{$org}{name}; 
     } 

   }  
   if ( $len && $name ) {  
      unless ( $len_org eq $name_org ){  
         print "name and length come from different organisms\n";   
         print "$uni_orgs{$len_org}{entry}\n"; 
         print "$uni_orgs{$name_org}{entry}\n"; 
     } else {   
       # check if lprim and rprim are the same 
       if ( $uni_orgs{$len_org}{rprim} eq $uni_orgs{$opt{species}}{rprim} 
            &&  $uni_orgs{$len_org}{lprim} eq $uni_orgs{$opt{species}}{lprim}) {    

            my @ent = qw ( rprim lprim name chr acc ) ; 
     
            print O "$id\t$uni_orgs{$opt{species}}{rprim}\t$uni_orgs{$opt{species}}{lprim}\t" . 
                  "$len\t$name\t$uni_orgs{$opt{species}}{chr}\t$uni_orgs{$opt{species}}{acc}\t". 
                   "$opt{species}\n" ; 
       }
     } 
   } 
} 
close(O) ;    
print "**done**\n" ; 




