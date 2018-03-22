# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 NAME

prune_keys_in_cdna2genome_targetted_proteinseq.pl

=head1 SYNOPSIS

prune_keys_in_cdna2genome_targetted_proteinseq.pl

=head1 DESCRIPTION

This script cleans up a mixture of headers with either 1 accession
number or two accession numbers in fasta protein files.  The single
accession numbers are standard ones from e.g. downloading a protein
file from SwissProt or RefSeq. The duo accession numbers are often 
the result of running the ensembl-analysis/scripts/get_protein_seq_from_cdnas2.pl
script, which can be used to retrieve protein sequences for *cDNAs* 
use warnings ;
used in cdna2genome targetted genebuild step.

Example duo accession number headers:

       >NM_001105760.2  NP_001099230.2
 
       >AF227200.1  AAF98119.1

In the same fasta file, we might have single-accessioned headers with 
these protein accession numbers (from SwissProt/RefSeq download):
        
       >NP_001099230.2
     
       >AAF98119.1

As protein indexing won't work with headers with two accessions 
(>NM_001105760.2  NP_001099230.2), we need to "factorise" the headers, 
even if it means duplicating sequences in the proteome file.  Also, headers
in the output file should be unique, so accessions NP_001099230.2 and 
AAF98119.1 in the above example will each appear once only. 

Following the above example, the output will look like this:

>NM_001105760.2
 MAETNEEVAVLVQRVVKDITNAFRRNPHIDEIGLIPCPEARYNRSPIVLVENKLGVESWC
.....

>NP_001099230.2
 MAETNEEVAVLVQRVVKDITNAFRRNPHIDEIGLIPCPEARYNRSPIVLVENKLGVESWC
.....

>AF227200.1
 MYKPVDPHSRMQSTYSYGMRGGAYPPRYFYPFPVPPLLYQVELSVGGQQFNGKGKMRPPV
..... 

>AAF98119.1
 MYKPVDPHSRMQSTYSYGMRGGAYPPRYFYPFPVPPLLYQVELSVGGQQFNGKGKMRPPV
.....


No. of redundant keys (accessions) and which one they are will be
written to STDOUT.


=head1 OPTIONS

-combined_infile   The redundant fasta file of protein sequences with a
                   mixture of single- and duo-accessioned headers.

-out_file          Output file with cleaned, non-redundant, single headers


=head1 EXAMPLE

perl prune_keys_in_cdna2genome_targetted_proteinseq.pl -combined_infile All_targetted_proteome_redun.fa -out_file All_targetted_proteome_final.fa > pruning_error.log

=cut



#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::SeqIO; 

#  SETTING THE DEFAULT VALUES FOR DB-CONNECTION ETC
my %opt; 
&GetOptions( 
             \%opt , 
            'combined_infile=s', # contans seq from cdna2genome and proteins used for pmatch / targettedExonerate  
            'out_file=s', # contans seq from cdna2genome and proteins used for pmatch / targettedExonerate  
           ) ; 


unless ( $opt{combined_infile} ) {  
     print STDERR " \n\n use -combined_infile <FILE> \n" ; 
     exit(0); 
} 
my %written_seqs; 
my $seq_in  = Bio::SeqIO->new(   -file => "<$opt{combined_infile}", -format => "fasta" ); 
my $seq_out = Bio::SeqIO->new( -file => ">$opt{out_file}", -format => 'fasta' );  


while( my $seq = $seq_in->next_seq() ) {

          my $prim_id = $seq->display_id ; 
          my $sec_id = $seq->desc;   

          if ( exists $written_seqs{$prim_id} ) {  
            print STDOUT "Primary ID $prim_id has already been written previously, skipped.\n" ;  
          }  

          if ( exists $written_seqs{$sec_id} ) { 
            print STDOUT "Secondary ID $sec_id has already been written previously, skipped.\n" ;  
          }  
            
           if ( $prim_id ) {  
              unless ( $written_seqs{$prim_id} ) { 
                $seq->display_id($prim_id) ;  
                $seq->desc("");
                $seq_out->write_seq($seq); 
                $written_seqs{$prim_id} = 1 ; 
              }
           } 

           if ( $sec_id ) {  
              unless ( $written_seqs{$sec_id} ) { 
                $seq->display_id($sec_id) ;  
                $seq->desc("");
                $seq_out->write_seq($seq); 
                $written_seqs{$sec_id} = 1 ; 
              }
           }
        }

exit(0);
my %acc_keys ; 

my %redundant_keys;  
open(C,"$opt{combined_infile}") || die "can't read file \n" ;   
my  @file = <C> ;  
my %np2nm;   

 #
 # we need to chop off the versions as the indexing is not done on versions laters os 
 #  xxxx.1 will be redundante with xxxx.2 
 # 
my @fa_seq; 
for ( my $i=0 ; $i<=@file; $i++) {  
   my $line = $file[$i]; 
   if ( $line =~m/^>/ ) {   
      $line=~s/^>//; 
      my @acc = split /\s+/,$line ;     
      map { s/\..*// } @acc;
      
      $np2nm{$acc[1]}=$acc[0];  
     
      for my $a ( @acc ) { 
        unless ( $acc_keys{$a}) {  
           $acc_keys{$a}=1;
        } else {  
          #print "key exists twice ...$a \n" ;  
          $redundant_keys{$a}=1; 
        }  
      }
   } 
} 

print STDOUT "found " .scalar(keys %redundant_keys ) . " redundant keys \n" ;  

open(O,">$opt{out_file}") || die " cant write to file $opt{out_file}\n" ; 

for my $line ( @file ) {  
   if ( $line =~m/^>/ ) {     
      $line=~s/>//; 
      my @acc = split /\s+/,$line ;    
      map { s/\..*// } @acc;
      if ( scalar(@acc) > 1 ) {  
         # we have 2 acc in header 
         if ( exists $redundant_keys{$acc[1]} ) {   
           $line=~s/$acc[1]\..*//g; 
         } 
      }   
      $line = ">".$line; 
   }
   print O $line ; 
} 
close(O); 
