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
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::DB::Flat::BinarySearch;
use Bio::SeqIO;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $index_dir,
    $index_file,
    $guess_gene_names,
    $all_gene_names,
);

$dbuser = 'ensro';
$dbport = 3306;

GetOptions(
            'dbname|db|D=s' => \$dbname,
            'dbuser|user|u=s' => \$dbuser,
            'dbhost|host|h=s' => \$dbhost,
            'dbport|port|P=s' => \$dbport,
            'dbpass|pass|p=s' => \$dbpass,
            'indexdir=s' => \$index_dir,
            'indexfile=s' => \$index_file,
            'guessgenenames' => \$guess_gene_names,
            'allgenenames' => \$all_gene_names,
            );

my $qy_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
        '-dbname' => $dbname,
        '-host' => $dbhost,
        '-user' => $dbuser,
        '-port' => $dbport,
        '-pass' => $dbpass
        );

die "You must supply an index directory and file\n"
    if not defined $index_dir or not defined $index_file;

my %classes = ('Ig-Heavy'        => "Immunoglobulin heavy chain",,
               'Ig-Light-Lambda' => "Immunoglobulin Lambda light chain" ,
               'Ig-Light-Kappa'  => "Immunoglobulin Kappa light chain",
               'TcR-Alpha'       => "T-cell receptor alpha",
               'TcR-Beta'        => "T-cell receptor beta",
               'TcR-Delta'       => "T-cell receptor delta",
               'TcR-Gamma'       => "T-cell receptor gamma",
               );
               

my %biotypes = (
  'V_segment' => 'V gene segment',
  'D_segment' => 'D gene segment',
  'J_segment' => 'J gene segment',
  'C_segment' => 'C gene segment',
                );


               

my $index = new Bio::DB::Flat::BinarySearch(-index_dir => $index_dir,
                                         -dbname    => $index_file,
                                         -format    => 'IMGT_embl',
                                            );

my %genes_by_slice;

if (@ARGV) {
  foreach my $slid (@ARGV) {
    my $sl = $qy_db->get_SliceAdaptor->fetch_by_name($slid);
    my $tl_sl = $qy_db->get_SliceAdaptor->fetch_by_region('toplevel',
                                                       $sl->seq_region_name);
    foreach my $tp (keys %biotypes) {
      my @genes = @{$sl->get_all_Genes_by_type($tp)};
      @genes = map { $_->transfer($tl_sl) } @genes;
      push @{$genes_by_slice{$sl->seq_region_name}}, @genes;
    }
  }

} else {
  foreach my $tp (keys %biotypes) {
    foreach my $g (@{$qy_db->get_GeneAdaptor->fetch_all_by_type($tp)}) {
      push @{$genes_by_slice{$g->slice->seq_region_name}}, $g;
    }
  }
}

my (%gene_to_name, %name_to_gene, %gene_to_desc, %tran_to_gene, %all_genes);

foreach my $sr_id (keys %genes_by_slice) {
  my @genes = sort { $a->start <=> $b->start } @{$genes_by_slice{$sr_id}};

  my (%seq_info_cache);

  foreach my $g (@genes) {

    $all_genes{$g->dbID} = $g;
    my (%support, %tran_support);

    foreach my $e (@{$g->get_all_Exons}) {
      foreach my $f (@{$e->get_all_supporting_features}) {
        $support{$f->hseqname} = $f->percent_id;
      }
    }
    foreach my $t (@{$g->get_all_Transcripts}) {
      my ($f) = @{$t->get_all_supporting_features};
      $tran_support{$f->hseqname} = $f->percent_id;

      my ($id, $reg_st, $reg_en) = $f->hseqname =~ /^(\S+)\.(\d+)\-(\d+)(\-\d)?$/;
      
      if (not exists $seq_info_cache{$id}) {
        $seq_info_cache{$id} = &extract_seq_info($id);
      }
      my $info = $seq_info_cache{$id};
  
      printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
             $g->stable_id ? $g->stable_id : $g->dbID,
             $t->stable_id ? $t->stable_id : $t->dbID,
             "IMGT/LIGM_DB",
             $f->hseqname,
             $f->hseqname,
             "", 
             "XREF");
    }

    my (%gnames, %descs, %supported_descs);


    foreach my $hseqname (keys %support) {
      my ($id, $reg_st, $reg_en) = $hseqname =~ /^(\S+)\.(\d+)\-(\d+)(\-\d)?$/;

      if (not exists $seq_info_cache{$id}) {
        $seq_info_cache{$id} = &extract_seq_info($id);
      }

      my $info = $seq_info_cache{$id};
      my @greg = &find_overlapping_genes($reg_st, $reg_en, $info->{gene_list});

      # store 2 bits of info for each symbol:
      # - the number of support-sequences that refer to that symbol;
      # - the average percent-id of the evidences that refer to that symbol

      foreach my $reg (@greg) {
        if (not exists $gnames{$reg->{gene}}) {
          $gnames{$reg->{gene}} = {
            name    => $reg->{gene},
            count   => 0,
            percent => 0,
          };
        } 
        $gnames{$reg->{gene}}->{count}++;
        $gnames{$reg->{gene}}->{percent} += $support{$hseqname};
        if (exists $tran_support{$hseqname}) {
          $gnames{$reg->{gene}}->{supported} = 1;
        }
      }

      my @desc_words = grep { exists($classes{$_}) } @{$info->{keywords}};
      if (scalar(@desc_words) == 1) {
        $descs{$desc_words[0]}++;
        if (exists $tran_support{$hseqname}) {
          $supported_descs{$desc_words[0]} = 1;
        }
      }
    }


    $gene_to_name{$g->dbID} = {};
    $gene_to_desc{$g->dbID} = {};

    if (scalar(keys %descs) > 1 and scalar(keys(%supported_descs)) == 1) {
      %descs = %supported_descs;
    }

    map { $gene_to_name{$g->dbID}->{$_} = $gnames{$_} } keys %gnames;
    map { $name_to_gene{$_}->{$g->dbID} = 1 } keys %gnames;
    map { $gene_to_desc{$g->dbID}->{$_} = 1 } keys %descs;
  }
}

&find_best_name_and_desc(\%gene_to_name, \%gene_to_desc, \%name_to_gene);

foreach my $gid (sort { 
  $all_genes{$a}->slice->seq_region_name cmp $all_genes{$b}->seq_region_name or 
      $all_genes{$a}->start <=> $all_genes{$b}->start
    } keys %gene_to_name) {
  my $g = $all_genes{$gid};

  my $desc = "";
  if (defined $gene_to_desc{$gid}) {
    $desc = $gene_to_desc{$gid};
  }

  my @names; 
  if (@{$gene_to_name{$gid}}) {
    if ($guess_gene_names or
        scalar(@{$gene_to_name{$gid}}) == 1) {
      push @names, $gene_to_name{$gid}->[0]->{name};
    } elsif ($all_gene_names) {
      foreach my $nh (@{$gene_to_name{$gid}}) {
        push @names, sprintf("%s(%d,%.1f%s)", 
                             $nh->{name}, 
                             $nh->{count}, 
                             $nh->{av_percent}, 
                             $nh->{supported} ? "*" : "");
      }
    }
  }
  
  foreach my $t (@{$g->get_all_Transcripts}) {
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
           $g->stable_id ? $g->stable_id : $g->dbID,
           $t->stable_id ? $t->stable_id : $t->dbID,
           "IMGT/GENE_DB",
           @names ? join(":", @names) : "NONAME", 
           @names ? join(":", @names) : "NONAME", 
           $desc ? join(" ", $classes{$desc}, $biotypes{$t->biotype}) : "",
           "KNOWN");
  }
}


###############################################################

sub extract_seq_info {
  my ($id) = @_;

  my $seq = $index->get_Seq_by_id($id);

  my @gene_regs;
  foreach my $val ($seq->get_SeqFeatures) {
    if ($val->primary_tag =~ /^(gene)|(\S+\-GENE)|(\S\-REGION)$/ and
        $val->annotation->get_Annotations("gene")) {
      my ($tag_val) = $val->annotation->get_Annotations("gene");

      my $name = $tag_val->value;
      $name =~ s/^\s*(.+\S+)\s*/$1/;

      # ignore the symbol as "complex" if it contains whitespace
      if ($tag_val->value !~ /\s+/) {
        push @gene_regs, {
          start => $val->location->start,
          end   => $val->location->end,
          gene  => $tag_val->value,
        };
      }
    }
  }

  @gene_regs = sort { $a->{start} <=> $b->{start} or $a->{end} <=> $b->{end} } @gene_regs;

  my $max_end_upstream;
  my $idx_of_max_end_upstream;

  for (my $i=0; $i < @gene_regs; $i++) {
    my $this_gene = $gene_regs[$i];
    
    if (not defined $max_end_upstream or
        $this_gene->{end} > $max_end_upstream) {
      $max_end_upstream = $this_gene->{end};
      $idx_of_max_end_upstream = $i;
    }

    $this_gene->{max_end_up} = $max_end_upstream;
    $this_gene->{max_end_up_idx} = $idx_of_max_end_upstream;
  }

  return {
    description => $seq->description,
    keywords    => [$seq->get_keywords],
    gene_list   => \@gene_regs,
  };
}


sub find_overlapping_genes {
  my ($st, $en, $list) = @_;

  # find min j s.t. segs[j].pos.s > end_pos
  #  strategy: binary search 

  my $left = 0;
  my $right = scalar(@$list);

  while ($left < $right) {
    my $mid = int (($left + $right) / 2);

    if ($list->[$mid]->{start} < $en) {
      $left = $mid + 1;
    } else {
      $right = $mid;
    }
  }

  my @overlapping;
  for(my $j = $left -1; $j >= 0; $j--) {

    if ($list->[$j]->{max_end_up} < $st) {
      last;
    } elsif ($list->[$j]->{start} <= $en and
             $list->[$j]->{end} >= $st) {
      push @overlapping, $list->[$j];
    }
  }

  return @overlapping;
}

sub find_best_name_and_desc {
  my ($g2name, $g2desc, $name2g) = @_;

  my @gids = keys %$g2name;

  foreach my $gid (sort {$a <=> $b} @gids) {
    my ($gname, $gdesc);
  
    my @gnames = keys  %{$g2name->{$gid}};
    my @descs = keys %{$g2desc->{$gid}};

    my @nh = values %{$g2name->{$gid}};
    foreach my $nh (@nh) {
      $nh->{av_percent} = $nh->{percent} / $nh->{count};
    } 
    @nh = sort { $b->{av_percent} <=> $a->{av_percent} or $b->{count} <=> $a->{count} } @nh;
    $g2name->{$gid} = \@nh;
    
    if (scalar(@descs) == 1) {
      ($gdesc) = @descs;
    } elsif (scalar(@descs) > 1) {
      # see if there is a clear winner
      @descs = sort { $g2desc->{$gid}->{$b} <=> $g2desc->{$gid}->{$a} } @descs;
      if ($g2desc->{$gid}->{$descs[0]} > $g2desc->{$gid}->{$descs[1]}) {
        ($gdesc) = @descs;
      }
    }

    $g2desc->{$gid} = $gdesc;
  }
}
