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
           outfile => "protein_classification.sql" , 
           filename_acc_no_desc => "acc_no_desc_found.txt" ,
           paf_file => 'paf_file.txt' ,  
           all_entries_hash_file => 'all_entries_hash_file.HASH',
          );

# 
my $default_paf_logic_name = 'Uniprot' ; 
my $paf_logic_name = '';
&GetOptions( 
             \%opt , 
            'config=s',             
            'delete' , # avoid interactive mode and overwrite files 
            'test!', 
            'output_dir=s', 
            'outfile=s',
            'paf_logic_name=s',  
            'filename_acc_no_desc=s',  
            'verbose!', 
            'paf_file=s',  
            'rerun!', 
           ) ;   


my %paf;  
my %paf_hash; 
     open (PAF,$opt{paf_file}) || die " cant read file $opt{paf_file}\n" ;  
     for ( <PAF> ) {  
       my ( $logic_name, $hit_name ) = split ; 
       ${$paf_hash{$logic_name}}{$hit_name} = 1 ;   
     } 
    close(PAF);  


print "fetching *done*\n" ; 

my @protein_analysis; 
my @tmp_all_ids; 
for my $lg( keys %paf_hash ) {  
  print "Analysis : $lg   hits_found:  " . scalar(keys %{$paf_hash{$lg}}) . "\n";   
  push @tmp_all_ids, keys %{$paf_hash{$lg}}; 
}  
my %hslice ; 
@hslice{@tmp_all_ids}=1 ; 
my @all_ids = keys %hslice ;  
print "\n\nHave " . scalar(@all_ids) . " ID's to fetch ...\n" ; 
#sleep(3);  



# fetch analysis for input_id types  


  print "running mfetch to fetch description etc ..\n" ; 
  my @hit_ids_investigate = @all_ids ; 

  my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch->new();    
  $obj->verbose(1) if $opt{verbose}; 
  my @fields_to_fetch  = qw ( Taxon acc ) ;  
  
  # split into small chunks to save data every 5000 entries ...    
  
  my @hit_id_chunks;   

  print " have " .scalar(@hit_ids_investigate) . " ids to investigate ...\n" ;   

  # add routine which 
  #   - checks if $opt{all_entries_hash_file} exists  
  #   - reads  $opt{all_entries_hash_file}  
  #   comparse it with @hit_ids_investigate
  #   if entry is found in @hit_ids_investigate and exists in all_entries which has been read remove it from hit_ids_to_investigate ... 
  #   
  #   also check array if it's in the list of no-match freatues 
  #   remove those as well such we don't re-fetch them 
  # thank ask use if you like to proceeed or delete file and restart ....    
  my %all_entries ;  
  if ( 0 ) { 
  if ( -e $opt{all_entries_hash_file} ) {   
    print "hash file found ... try to read it ... \n " ;  
    #dbmopen(%all_entries, "$opt{all_entries_hash_file}", 0444) or die "Cant open  $opt{all_entries_hash_file} file\n"; 
    #dbmclose(%all_entries );
    print scalar ( keys %all_entries  ) . " keys found - comparing ...\n" ;   
    my %acc_not_in_mfetch ; 
    my @not_fetched_now ;   
    my $alr_fetched = 0; 
    my $no_match = 0 ; 
    for my $acc  ( @hit_ids_investigate ) {    
      unless ( $all_entries{$acc} ) {  
        unless ( $acc_not_in_mfetch{$acc} ) { 
         push @not_fetched_now, $acc ;  
        }else {  
         $no_match++ ; 
        }  
      } else {  
        $alr_fetched++; 
      } 
    }   
    print " ... hash read ... \n $no_match entries found in hash which don't have entry in mfetch\n" ; 
    print " ... hash read ... \n $alr_fetched entries found which i've already fetched\n" ; 
    print " i will carry on with " . scalar(@not_fetched_now) . " entries which i haven't fetched now ....\n" ;   
   @hit_ids_investigate = @not_fetched_now ;  
   sleep(5); 
 } 
 }  

  while ( @hit_ids_investigate) {  
    my @tmp = splice @hit_ids_investigate, 0, 5000 ; 
    push @hit_id_chunks, \@tmp  ; 
  }  

  my (  $entries, $not_found ) ; 
  my @all_not_found ;   
  for my $acc_ref  ( @hit_id_chunks )  {  
     ($entries, $not_found ) =  @{$obj->get_Entry_Fields($acc_ref ,\@fields_to_fetch) } ;   
     %all_entries = ( %all_entries, %$entries ) ;     

     #dbmopen(%all_entries, $opt{all_entries_hash_file}, 0666) or die "Cant open testdb file\n";
     #dbmclose(%all_entries);

       print keys ( %all_entries ) . " entries in all_entries\n" ;  
       print scalar( @all_not_found) . " entries not found \n" ;  
       push @all_not_found , @$not_found ;   
  }    

  my @not_found = @all_not_found ; 

  print "\nNO description found by mfetch for " . scalar(@not_found) . " entries !!!! - we try with wildcards now ...  \n" ;    
     for ( @not_found ) {  
       print "$_ not_found\n" ; 
     }  
  # print now try to  fetch with wildcards     - second round 
  my @not_found_2 ;  
  for my $acc ( @not_found ) {   
     #print scalar(keys %all_entries)  . " in all entries \n" ; 
     my $nacc = $acc ; 
     $nacc=~s/\..*//g; 
     #print " fetching $acc ---->  $nacc\n" ;   

     my $search_acc = " -i \"acc:$nacc\%\"";  
     my ($entries, $not_found ) =  @{$obj->get_Entry_Fields($search_acc, \@fields_to_fetch) }   ;   

     if ( scalar(@$not_found) > 0 ) {  
        warning( "No entry found for $nacc\% with wildcard-search ") ; 
        push @not_found_2, $acc ;  
     }else {  

       #print "try adding re-fetched values to all_entries{$acc}\n" ; 
       my %e = %$entries ;    

       if ( exists $all_entries{$acc} ) { 
             print "\n there's already an entry for $acc :\n" ; 
             my %stored_entry = %{$all_entries{$acc}} ; 
             for my $se ( keys %stored_entry ) {   
               print "entry_stored: $acc  $se $stored_entry{$se}\n" ; 
             }
             foreach my $k ( keys %e ) {   
                print "k : $k\n" ; 
                my %tmp = %{$e{$k}};  
                for my $field ( keys %tmp  ) {   
                   print " fetched now : acc $acc $field  $tmp{$field} \n" ;  
               }  
            } 
            print ( "entry already exists: $acc\n") ; 
        }

       if ( 1) { 
       foreach my $k ( keys %e ) {   
         print "k : $k\n" ; 
         my %tmp = %{$e{$k}};  
         for my $field ( keys %tmp  ) {   
            print " fetched now : acc $acc $field  $tmp{$field} \n" ;  
              $all_entries{$acc}{$field}  = $tmp{$field} ;     
           }  
         } 
      }     
     }
    }
  #print scalar(keys %all_entries)  . " in all entries \n" ; 

  for my $acc( keys  %all_entries ) {  
     print "STORED   $acc  " ; 
     for my $k ( keys %{$all_entries{$acc}}) {  
          print   ${$all_entries{$acc}}{$k}." ";
     }  
    print "\n" ; 
  } 


  open (A,">$opt{filename_acc_no_desc}") || die "Can't write to file \n" ; 
  for ( @not_found_2) { 
     print A "entry_really_not_found_by_mfetch_for  $_\n" ; 
  }  
  close(A); 
  print " descriptions fetched ........\n" ;  



sub get_and_exclude_protein_hits_by_field { 
  my ($entries,$field,$expr_to_match,$expr_to_exclude) = @_ ;   
  my %result ;     

  unless ( $expr_to_exclude ) {  
    $expr_to_exclude = [];
  } 
  my $expression_to_match = 0 ; 
  if ( scalar(@$expr_to_match) > 0 ) { 
    $expression_to_match = 1
  }  

  ACC: for my $acc ( keys %$entries ) {       
     my $attr_to_test = ${$entries}{$acc}{$field} ;   
   
     # this is the routine if you like to match AND exclude    
     #print "processing $acc\n" ;   

     if ( $expression_to_match ) {  
       for my $string_match ( @$expr_to_match ) {   
         # print "testing : $string_match  - $attr_to_test " ; 
         if ($attr_to_test=~m/$string_match/ ) {   
           print "MATCH  $attr_to_test $string_match $acc \n" if $opt{verbose} ; 
            # we have a match 
            if ( scalar(@$expr_to_exclude) > 0 ) {   
              for my $exclude ( @$expr_to_exclude ) {  
                if ( $attr_to_test =~m/$exclude/ ) {   
                  print "excluding : $exclude - $attr_to_test \n" if $opt{verbose} ;  
                  # we need to exclude this as it's listed to be excluded .......
                } else { 
                   $result{$acc}= ${$entries}{$acc};  
                } 
              }             
            } else { 
              $result{$acc}= ${$entries}{$acc};  
            } 
            next ACC; 
         } else {  
           if ( $opt{verbose} ) {  
             print " sorry - no match. the ACC $acc does not match the ".
             " $field-field \"$string_match\" criteria [ $attr_to_test=~m/$string_match/ :-( ] \n" ;
           }  

         }
       }
     } else {    
        #print "we only have to exclude stuff, not to match stuff - for field $field .... \n" ; 
         for my $exclude ( @$expr_to_exclude ) {  
            if ( $attr_to_test =~m/$exclude/ ) {   
                #print "exclude matches : $exclude $attr_to_test \n" ; 
                # we need to exclude this as it's listed to be excluded .......
             } else { 
                  #print "no match ----storing \n" ;  
                  $result{$acc}= ${$entries}{$acc};  
             } 
        }             
     } 
  }   # end of big loop ..... 
  return \%result;
}



 
 
