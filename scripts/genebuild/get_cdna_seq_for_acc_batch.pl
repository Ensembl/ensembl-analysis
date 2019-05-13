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

use warnings ;
use strict;
use Bio::EnsEMBL::DBSQL::DBConnection; 
use Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor; 
use Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch;   
 
use Bio::EnsEMBL::Pipeline::Analysis; 
use Bio::EnsEMBL::Utils::Exception qw ( throw warning ) ;  
use Bio::EnsEMBL::Utils::Exception qw ( throw warning ) ;  
use Bio::EnsEMBL::Analysis::Tools::Utilities qw ( get_input_arg ) ; 

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor; 
use Getopt::Long;
$|=1; 
#  SETTING THE DEFAULT VALUES FOR DB-CONNECTION ETC
my %opt = ( 
           #infile  => "acc.all" , 
           infile  => "tmp.acc" , 
           outfile => "fasta.seq" , 
          );

# 
&GetOptions( 
             \%opt , 
            'infile=s', 
            'outfile=s', 
           ) ;   

# read file with acc's to fetch 
my @all_acc ; 
open (PAF,$opt{infile}) || die " cant read file $opt{infile}\n" ;    
  for ( <PAF> ) { 
    chomp ; 
    push @all_acc , $_ ;  
  } 
close(PAF);  

# create Mfetch seqfetcher 
my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch->new();    
$obj->verbose(0); 
$obj->options("-d emblrelease");   
  
# chunk acc ...
my @hit_id_chunks;   
print " have " .scalar(@all_acc) . " ids to investigate ...\n" ;   

while ( @all_acc) {  
  my @tmp = splice @all_acc, 0, 5000 ; 
  push @hit_id_chunks, \@tmp  ; 
}   

system("rm $opt{outfile}") ;  

my (  $entries, $not_found ) ;  
my @all_not_found ;    
my $count = 0 ; 
for my $acc_ref  ( @hit_id_chunks )  {  
     ($entries, $not_found ) =  @{$obj->get_Seq_BatchFetch($acc_ref) } ;   

    for ( @$entries ) {
     open (FA,">>$opt{outfile}") || die "cant write to file $opt{outfile}\n" ;
     print FA "$_\n" ;
     close(FA) ;  
     $count++ if /^>/; 
    }

    print "  $count  entries found \n" ;  
    push @all_not_found , @$not_found ;   
    print scalar( @all_not_found) . " entries not found \n" ;  
}     


  #
  # single fetchs  with wildcards - second round  
  #

  print "\nNo description found by mfetch for " . scalar(@all_not_found) . " entries !!!! - we try with wildcards now ...  \n" ;    

  my @not_found_2 ;    
  my @wildcard_found ;  

  for my $acc ( @all_not_found ) {    
     print "trying $acc\n" ; 
     my @entry =  @{$obj->get_Seq_by_acc_wildcard($acc)} ;  

     if ( scalar(@entry > 0 )) { 
       push @wildcard_found, @entry ; 
     }else {  
        warning( "No entry found for $acc\% with wildcard-search ") ;  
        push @not_found_2, $acc ; 
     }
  } 

    for ( @wildcard_found) {
     open (FA,">>$opt{outfile}") || die "cant write to file $opt{outfile}\n" ;
     print FA "$_\n" ;
     close(FA) ;  
     $count++ if /^>/; 
    }

open (A,">no_acc_found.txt") || die "Can't write to file \n" ; 
 for ( @not_found_2) {  
    print; 
    print A "entry_really_not_found_by_mfetch_for  $_\n" ; 
 }  
close(A); 
