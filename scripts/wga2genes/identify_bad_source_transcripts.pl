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

# identify_bad_source_transcripts.pl
# 
# Scans the given database for protein_coding genes and
# transcripts that are unsuitable for projecting onto
# a low-coverage genome via a whole-genome alignment
#
# The resulting set can act as the starting point for a "kill list"  
# of genes and transcripts that will be ignored during the
# WGA2Genes process. 

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass, 
    $look_for_dodgy,
    $print_gsi , 
);

$dbuser = 'ensro';
$dbport = 3306;

GetOptions(
	'dbname|db|D=s' => \$dbname,
	'dbuser|user|u=s' => \$dbuser,
	'dbhost|host|h=s' => \$dbhost,
	'dbport|port|P=s' => \$dbport,
	'dbpass|pass|p=s' => \$dbpass,
        'dodgy'    => \$look_for_dodgy,
        'print_gsi' => \$print_gsi ,   # prints gene stable id of rejected genes
);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);
my @slices;
if (@ARGV) {
  foreach my $sid (@ARGV) {
    push @slices, $db->get_SliceAdaptor->fetch_by_name($sid);
  }
} else {
  @slices = @{$db->get_SliceAdaptor->fetch_all('toplevel')};
}


my %reject;
my @rejected_g; 

foreach my $sl (sort { $a->length <=> $b->length } @slices) {
  my @orig_genes = @{$sl->get_all_Genes_by_type('protein_coding')};
  my @usable_genes;

  #
  # Reject Mitochondrial genes
  #
  if ($sl->seq_region_name eq 'MT') {
    map { $reject{$_->stable_id} =  "Mitochondrial protein" } @orig_genes;   

    for my $g ( @orig_genes ) { 
      push @rejected_g, $g->stable_id;  
    } 

    next;
  }

  my @genes;

  #
  # Reject protein_coding genes which contain not even a single protein-coding transcript. Sometimes
  # this happens in the human gene set when inconsistencies arise during the Ensembl-Havana merge.
  #

  foreach my $gene (@orig_genes) {
    my $coding_transcript_cnt = 0;
    my @transcripts = @{$gene->get_all_Transcripts};
    foreach my $tr(@transcripts) {
      if ($tr->biotype eq 'protein_coding') {
         $coding_transcript_cnt ++;
      }
    }
    if ($coding_transcript_cnt == 0) {
       # print STDERR "gene ".$gene->stable_id." contains ".scalar(@transcripts)." transcripts but none of them are protein_coding.\n";
       push @rejected_g, $gene->stable_id;
    } else {
       push (@usable_genes, $gene);
    }
  }


  foreach my $g (@usable_genes) {
    my ($gst, $gen);

    my @trans;
    foreach my $t (@{$g->get_all_Transcripts}) {

      if ( $t->translation ) {  
        my $x_tr = $t->translateable_seq ; 
        my $x_tl = $t->translation->seq ;  
  
        my $ratio = length($x_tr) / length($x_tl ) ;    
          if ( $x_tl=~m/^X/ ) { 
             $reject{$t->stable_id} = "Partial start-codon ( Translation starts with X )" ;  
             push @rejected_g, $g->stable_id;  
          } 
        }
        
      if (length($t->translateable_seq) % 3 != 0) {
        $reject{$t->stable_id} = "Non modulo-3 coding length";
        push @rejected_g, $g->stable_id;  
      } else {
        if ($look_for_dodgy) {
          my ($cst, $cen);
          
          foreach my $e (@{$t->get_all_translateable_Exons}) {
            $cst = $e->start if not defined $cst or $cst > $e->start;
            $cen = $e->end   if not defined $cen or $cen < $e->end;
          }
          
          push @trans, {
            start => $cst,
            end   => $cen,
            tran => $t,
          }
        }
      }
    }

    if (@trans) {
      foreach my $t (@trans) {
        $gst = $t->{start} if not defined $gst or $gst > $t->{start};
        $gen = $t->{end}   if not defined $gen or $gen < $t->{end};
      }
      
      push @genes, {
        start => $gst,
        end   => $gen,
        gene  => $g,
        transcripts => \@trans,
      }
    }
  }



  @genes = sort { $a->{start} <=> $b->{start} } @genes;

  my @clusters;
  foreach my $g (@genes) {
    if (not @clusters or $clusters[-1]->{end} < $g->{start} - 1) {
      push @clusters, {
        start => $g->{start},
        end   => $g->{end},
        genes => [$g],
      };
    } else {
      if ($g->{end} > $clusters[-1]->{end}) {
        $clusters[-1]->{end} = $g->{end};
      }
      push @{$clusters[-1]->{genes}}, $g;
    }
  }

  # we're interested in clusters that become separated 
  # when you remove on (apparently problemattic) gene

  foreach my $c (@clusters) {
    if (@{$c->{genes}} > 1) {      
      # indentify gene in the cluster where one transcript overlap
      # genomically with other genes, but other transcript(s) do
      # not

      foreach my $g (@{$c->{genes}}) {
        my @inside;

        foreach my $og (@{$c->{genes}}) {
          next if $og eq $g;
          
          if ($og->{start} > $g->{start} and 
              $og->{end}   < $g->{end}) {
            push @inside, $og;
          }
        }

        if (@inside) {
          # other genes fit inside this one. But is it just a dogdy
          # transcript that is causing this overlap? 
          my (@overlaps, @doesnt);
          foreach my $t (@{$g->{transcripts}}) {
            my $all_inside = 1;
            foreach my $in (@inside) {
              if ($in->{start} < $t->{start} or
                  $in->{end} >   $t->{end}) {
                $all_inside = 0;
                last;
              }
            }
            if ($all_inside) {
              push @overlaps, $t;
            } else {
              push @doesnt, $t;
            }
          }

          if (scalar(@overlaps) == 1 and scalar(@doesnt) > 1) {
            map { $reject{$_->{tran}->stable_id} = "Dodgy transcript" } @overlaps;
          }
        }
      }
    }
  }
} 

   my %tmp ; 
   @tmp{ @rejected_g }=1; 
  if ($print_gsi) {   
    print join ("\n",  keys %tmp ) ;  
    exit(0); 
  } 

my %reasons;
foreach my $id (keys %reject) {
  push @{$reasons{$reject{$id}}}, $id;
}

foreach my $reason (keys %reasons) {
  foreach my $id (@{$reasons{$reason}}) {
    print "$id\t$reason\n";
  }
}
