#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;
use buildchecks::ScriptUtils;
use Getopt::Long;

use Bio::Tools::CodonTable;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;

my $translate  = 1;
my $use_intron = 0;
my $host;
my $user;
my $dbname;
my $port;
my $comphost ;
my $compuser ;
my $compdbname;
my $compport;
my $outhost = undef;
my $outuser = undef;
my $outpass = undef;
my $outdbname  = undef;
my $outport = undef;
my $dnahost;
my $dnauser ;
my $dnadbname  = undef;
my $dnaport;
my @chromosomes;
my $file          = 'stdout';
my $usephase      = 1;
my $keepredundant = 0;
my $schema        = 19;
my @defcompgenetypes = ('ncbigene');
my @defgenetypes = ( 'ensembl', 'Known', 'Novel_CDS' );
my @genetypes;
my @compgenetypes;
my @analyses;
my $coordsystem    = 'chromosome';
my $set1_name      = 'hinxton';
my $set2_name      = 'ncbi';
my $trimstop       = 0;
my $display_name   = 0;
my $show_gene_type = 1;
my $start;
my $end;

my $selcys_file;

my $path;# = 'NCBI34';

$| = 1;

GetOptions( 'host:s'          => \$host,
            'user:s'          => \$user,
            'dbname:s'        => \$dbname,
            'path:s'          => \$path,
            'port:n'          => \$port,
            'comphost:s'      => \$comphost,
            'compuser:s'      => \$compuser,
            'compdbname:s'    => \$compdbname,
            'compport:n'      => \$compport,
            'outhost:s'       => \$outhost,
            'outuser:s'       => \$outuser,
            'outpass:s'       => \$outpass,
            'outdbname:s'     => \$outdbname,
            'outport:n'       => \$outport,
            'dnahost:s'       => \$dnahost,
            'dnauser:s'       => \$dnauser,
            'dnadbname:s'     => \$dnadbname,
            'dnaport:n'       => \$dnaport,
            'schema:n'        => \$schema,
            'chromosomes:s'   => \@chromosomes,
            'compgenetypes:s' => \@compgenetypes,
            'genetypes:s'     => \@genetypes,
            'analyses:s'      => \@analyses,
            'translate!'      => \$translate,
            'intron!'         => \$use_intron,
            'phasecheck!'     => \$usephase,
            'redundant!'      => \$keepredundant,
            'selcys_file:s'   => \$selcys_file,
            'file:s'          => \$file,
            'coordsystem:s'   => \$coordsystem,
            'set1_name:s'     => \$set1_name,
            'set2_name:s'     => \$set2_name,
            'trimstop'        => \$trimstop,
            'display_name'    => \$display_name,
            'show_gene_type'  => \$show_gene_type,
            'start:n'         => \$start,
            'end:n'           => \$end );

if ( scalar(@chromosomes) ) {
  @chromosomes = split( /,/, join( ',', @chromosomes ) );
}

if ( scalar(@compgenetypes) ) {
  @compgenetypes = split( /,/, join( ',', @compgenetypes ) );
} else {
  @compgenetypes = @defcompgenetypes;
}

if ( scalar(@genetypes) ) {
  @genetypes = split( /,/, join( ',', @genetypes ) );
} else {
  @genetypes = @defgenetypes;
}

if ( scalar(@analyses) ) {
  @analyses = split( /,/, join( ',', @analyses ) );
}

my $db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -user   => $user,
                                      -port   => $port,
                                      -dbname => $dbname );

my $dnadb;
if ($dnadbname) {
  $dnadb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                        -user   => $dnauser,
                                        -port   => $dnaport,
                                        -dbname => $dnadbname );

  $db->dnadb($dnadb);
}

my $sa   = $db->get_SliceAdaptor();
my $ga   = $db->get_GeneAdaptor();
my $dafa = $db->get_DnaAlignFeatureAdaptor();

my $compdb =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $comphost,
                                      -user   => $compuser,
                                      -port   => $compport,
                                      -dbname => $compdbname );

if ($dnadb) {
  $compdb->dnadb($dnadb);
}

my $compsa   = $compdb->get_SliceAdaptor();
my $compga   = $compdb->get_GeneAdaptor();
my $compdafa = $compdb->get_DnaAlignFeatureAdaptor();

my $outsa;
my $outga;
my $outdafa;
my $outdb;

if ($outdbname) {
  $outdb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $outhost,
                                        -user   => $outuser,
                                        -port   => $outport,
                                        -pass   => $outpass,
                                        -dbname => $outdbname, );

  if ($dnadb) {
    $outdb->dnadb($dnadb);
  }

  $outsa   = $outdb->get_SliceAdaptor();
  $outga   = $outdb->get_GeneAdaptor();
  $outdafa = $outdb->get_DnaAlignFeatureAdaptor();
}

my %selcys_hash;
if ($selcys_file) {
  open FPSC,"<$selcys_file" or die "Failed opening $selcys_file for read";
  while (<FPSC>) {
    chomp;
    $selcys_hash{$_} = 1;
  }
  close FPSC;
}

my %types_hash;
$types_hash{$set1_name} = \@genetypes;
$types_hash{$set2_name} = \@compgenetypes;

my $chrhash;
if ( $schema == 20 ) {
  $chrhash = get_chrlengths_v20( $db, $path, $coordsystem );
} else {
  $chrhash = get_chrlengths_v19( $db, $path );
}

filter_to_chr_list( \@chromosomes, $chrhash, $db->dbc->dbname );

# why disable buffering here...
$| = 1;
if ( $file ne "stdout" ) {
  open FP, ">$file";
} else {
  open FP, ">-";
}

my $old_fh = select(FP);

# ... and again here...
$| = 1;
select($old_fh);

my $n_match_trans                      = 0;
my $n_complete_translating_match_trans = 0;
my $n_twoway_cluster                   = 0;
my $n_set2_only_cluster                = 0;
my $n_set1_only_cluster                = 0;
my $n_set2_only_unclustered            = 0;
my $n_set1_only_unclustered            = 0;
my $n_cluster                          = 0;
my $n_unclustered                      = 0;
my $n_match_cluster                    = 0;
my $n_gene                             = 0;
my $n_comp_gene                        = 0;
my $n_trans                            = 0;
my $n_comp_trans                       = 0;

foreach my $chr (@{sort_chr_names($chrhash)}) {
  my $chrstart = (defined $start ? $start : 1);
  my $chrend   = (defined $end ? $end : $chrhash->{$chr});

  print STDERR "Chr $chr from $chrstart to " . $chrend. "\n";

  clear_coding_exons_cache();

#  $chrstart = 1;
#  $chrend = 5000000;
#  $chrstart = 31000000;
#  $chrend   = 32500000;

  my $n_chr_match_trans = 0;
  my $n_chr_complete_translating_match_trans = 0;
  my $n_chr_match_cluster = 0;

  my %outgenes;

  my $slice;
  my $comp_slice;
  my $slicename;
  if ($schema == 20) {
    $slicename = "$coordsystem:$path:$chr:$chrstart:$chrend:1";
    # print "Slice = $slicename\n";
    $slice = $sa->fetch_by_name($slicename);
    $comp_slice = $compsa->fetch_by_name($slicename);
  } else {
    #$slicename = "$chr:$chrstart-$chrend";
    #print "Slice = $slicename\n";
    $slice = $sa->fetch_by_chr_start_end($chr,$chrstart,$chrend);
    $comp_slice = $compsa->fetch_by_chr_start_end($chr,$chrstart,$chrend);
  }


  my @allgenes;

  print "Fetching genes\n";
  my %genes_hash;
  my @genes;
 # Hack to run 17 manual/ensembl comparison
  my @munged_genetypes;

  foreach my $genetype (@genetypes) {
    $genes_hash{$genetype} = $slice->get_all_Genes_by_type($genetype);
    print FP "Got " . scalar(@{$genes_hash{$genetype}}) . " $genetype genes\n";
 # Hack to run 17 manual/ensembl comparison
    foreach my $gene (@{$genes_hash{$genetype}}){
      $gene->biotype($gene->biotype . "_set1");
    }
    push @munged_genetypes,$genetype . "_set1";
    push @genes,@{$genes_hash{$genetype}};
  }

  if (scalar(@analyses)) {
    my @filtered_genes;
    foreach my $gene (@genes) {
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        my $found = 0;
        foreach my $logicname (@analyses) {
          if ($trans->analysis->logic_name eq $logicname) {
            $found = 1;
          }
        }
        if (!$found) {
          remove_transcript_from_gene($gene,$trans);
        }
      }
      if (scalar(@{$gene->get_all_Transcripts})) {
        push @filtered_genes,$gene;
      }
    }
    @genes = @filtered_genes;
  }
  print "Done fetching genes (fetched " . scalar(@genes) .")\n";

 # Hack to run 17 manual/ensembl comparison
  $types_hash{$set1_name} = \@munged_genetypes;

  $n_gene += scalar(@genes);

  push @allgenes,@genes;

  print "Fetching comp genes\n";
  my %comp_genes_hash;
  my @comp_genes;
  foreach my $genetype (@compgenetypes) {
    $genes_hash{$genetype} = $comp_slice->get_all_Genes_by_type($genetype);
    print FP "Got " . scalar(@{$genes_hash{$genetype}}) . " $genetype genes\n";
    push @comp_genes,@{$genes_hash{$genetype}};
  }
  print "Done fetching comp genes (fetched " . scalar(@comp_genes) .")\n";

  if ($trimstop) {
    print "Trimming stops\n";
    foreach my $gene ((@genes,@comp_genes)) {
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        if ($trans->translation) {
          my $coding_end   = $trans->cdna_coding_end;
          my $cdna_seq     = uc($trans->spliced_seq);
          my $endseq       = substr($cdna_seq,$coding_end-3,3);

          if ($endseq eq "TAG" || $endseq eq "TGA" || $endseq eq "TAA") {
            if ($trans->translation->end > 3) {
              $trans->translation->end($trans->translation->end()-3);
            } else {
              print "NOTE Moving stop would move it into previous exon\n";
            }
          }
        }
      }
    }
  }

  my $ntrans = count_trans(\@genes);
  my $ncomptrans = count_trans(\@comp_genes);

  $n_trans += $ntrans;
  $n_comp_trans += $ncomptrans;

  print FP "GENE INFO: Chr $chr " . scalar(@genes) .      " genes " .      $ntrans . " trans " .
                                    scalar(@comp_genes) . " comp_genes " . $ncomptrans . " comp_trans\n";

  $n_comp_gene += scalar(@comp_genes);

  push @allgenes,@comp_genes;

  my ($clusters, $unclustered) = cluster_Genes(\@allgenes, \%types_hash);

  print "Got " . scalar(@$clusters) . " clusters\n";

  $n_cluster += scalar(@$clusters);

  my $n_chr_set1_only_unclustered = 0;
  my $n_chr_set2_only_unclustered = 0;
  foreach my $uncl (@$unclustered) {
    my @inc_sets = @{$uncl->get_sets_included};
    print FP "UNCLUSTERED: Chr $chr " . $inc_sets[0];
    foreach my $gene (@{$uncl->get_Genes}) {
      my $max_exons_in_trans = 0;
      my $n_trans_with_stop = 0;
      my $total_stops_in_gene = 0;
      my $n_trans_with_fshift = 0;
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        my $nexon_in_trans = scalar(@{$trans->get_all_Exons});

        if ($nexon_in_trans > $max_exons_in_trans) {
          $max_exons_in_trans = $nexon_in_trans;
        }
        if ($trans->translation) {
          my $tlnseq = $trans->translate->seq;
          if ($tlnseq =~ /\*/) {
            $n_trans_with_stop++;
            $total_stops_in_gene += ($tlnseq =~ tr/*/*/);
          }

          my $prev_exon=undef;
          FEX:
          foreach my $exon (@{$trans->get_all_translateable_Exons}) {
            if (defined($prev_exon)) {
              my $intron_len;
              if ($trans->strand == 1) {
                $intron_len = $exon->start - $prev_exon->end -1;
              } else {
                $intron_len = $prev_exon->start - $exon->end -1;
              }
              if ($intron_len <= 7) {
                $n_trans_with_fshift++;
                last FEX;
              }
            }
            $prev_exon = $exon;
          }
        }
      }
      print FP " " . get_gene_id($gene) . " " . scalar(@{$gene->get_all_Transcripts}) . " " . scalar(@{$gene->get_all_Exons}) . " " . $max_exons_in_trans . " " . $n_trans_with_stop . " " . $n_trans_with_fshift . " " . $total_stops_in_gene . " " . $gene->biotype . " " . $gene->start . " " . $gene->end;
    }
    if ($inc_sets[0] eq $set1_name) {
      $n_set1_only_unclustered++;
      $n_chr_set1_only_unclustered++;

      # Temporary to see if in CCDS
      my $in_ccds = 0;
      foreach my $gene (@{$uncl->get_Genes}) {
        foreach my $trans (@{$gene->get_all_Transcripts}) {
          foreach my $xref (@{$trans->get_all_DBLinks}) {
            if ($xref->dbname eq 'CCDS') {
              $in_ccds++;
            }
          }
        }
      }
      print FP " in ccds = $in_ccds";
    } elsif ($inc_sets[0] eq $set2_name) {
      $n_set2_only_unclustered++;
      $n_chr_set2_only_unclustered++;

      # Temporary to see if in CCDS
      my $in_ccds = 0;
      foreach my $gene (@{$uncl->get_Genes}) {
        foreach my $trans (@{$gene->get_all_Transcripts}) {
          foreach my $xref (@{$trans->get_all_DBLinks}) {
            if ($xref->dbname eq 'CCDS') {
              $in_ccds++;
            }
          }
        }
      }
      print FP " in ccds = $in_ccds";
    } else {
      #print "Shouldn't happen not $set1_name or $set2_name in unclustered!!\n"
      print "Shouldn't happen, neither $set1_name nor $set2_name in unclustered!!\n"
    }
    print FP "\n";
  }

  $n_unclustered += scalar(@$unclustered);

  my @twoways;
  my $n_set1_redund = 0;
  my $n_set2_redund = 0;

  foreach my $cluster (@$clusters) {
    my @inc_sets = @{$cluster->get_sets_included};
    if (scalar(@inc_sets) == 2) {

      if (!$keepredundant) {
        $n_set1_redund += remove_redundant_cds_from_cluster($cluster,$set1_name);
        $n_set2_redund += remove_redundant_cds_from_cluster($cluster,$set2_name);
      }

      push @twoways, $cluster;
    } elsif ($inc_sets[0] eq $set1_name) {
      $n_set1_only_cluster++;
      print "SINGLE SET CLUSTER: $set1_name ";
      foreach my $gene (@{$cluster->get_Genes_by_Set($set1_name)}) {
        print " " . get_gene_id($gene);
      }
      print "\n";
    } elsif ($inc_sets[0] eq $set2_name) {
      print "SINGLE SET CLUSTER: $set2_name ";
      foreach my $gene (@{$cluster->get_Genes_by_Set($set2_name)}) {
        print " " . get_gene_id($gene);
      }
      print "\n";
      $n_set2_only_cluster++;
    } else {
      print "Shouldn't happen not $set1_name or $set2_name in single set cluster!!\n"
    }
  }
  print "Got " . scalar(@twoways) . " twoways\n";

  $n_twoway_cluster += scalar(@twoways);

  printf FP
      "CLUSTER INFO: Chr %10.10s %4d clusters %4d unclustered %4d $set1_name"
    . "_only_uncl %4d $set2_name"
    . "_only_uncl %4d $set1_name"
    . "_redund_cds %4d $set2_name"
    . "_redund_cds %4d twoways\n", $chr, scalar(@$clusters),
    scalar(@$unclustered), $n_chr_set1_only_unclustered,
    $n_chr_set2_only_unclustered,
    $n_set1_redund, $n_set2_redund, scalar(@twoways);

  CLUSTER: foreach my $cluster (@twoways) {
    my @genes = @{$cluster->get_Genes_by_Set($set1_name)};
    my @comp_genes = @{$cluster->get_Genes_by_Set($set2_name)};


    my $cluster_match = 0;
    foreach my $gene (@genes) {
      my %hadmatch;
      foreach my $comp_gene (@comp_genes) {
        TRANS: foreach my $trans (@{$gene->get_all_Transcripts}) {
          my @exons;
          if (!$use_intron) {
            if ($translate) {
              if (!defined($trans->translation)) {
                next TRANS;
              }
              @exons = sort {$a->start <=> $b->start} @{$trans->get_all_translateable_Exons};
            } else {
              @exons = sort {$a->start <=> $b->start} @{$trans->get_all_Exons};
            }
          } else {
            @exons = sort {$a->start <=> $b->start} @{$trans->get_all_Introns};
          }


          COMPTRANS: foreach my $comp_trans (@{$comp_gene->get_all_Transcripts}) {
            my @comp_exons;
            if (!$use_intron) {
              if ($translate) {
                if (!defined($comp_trans->translation)) {
                  print FP "No translation\n";
                  next COMPTRANS;
                }
                @comp_exons = sort {$a->start <=> $b->start} @{$comp_trans->get_all_translateable_Exons};
              } else {
                @comp_exons = sort {$a->start <=> $b->start} @{$comp_trans->get_all_Exons};
              }
            } else {
              @comp_exons = sort {$a->start <=> $b->start} @{$comp_trans->get_all_Introns};
            }

            if (scalar(@comp_exons) != scalar(@exons)) {
                print FP "Number of exons differs: ", scalar(@comp_exons), " ", scalar(@exons), "\n";
                next COMPTRANS;
            }
#            next COMPTRANS if (scalar(@comp_exons) != scalar(@exons));

            my $e_count = 0;
            my $n_exon  = scalar(@exons);

            foreach my $exon (@exons) {
              my $comp_exon = shift @comp_exons;

              my $first_e;
              my $last_e;
              if ($exon->strand == 1) {
                $first_e = ($e_count == 0);
                $last_e  = ($e_count == ($n_exon-1));
              } else {
                $first_e = ($e_count == ($n_exon-1));
                $last_e  = ($e_count == 0);
              }

              # print "Comparing " . get_exon_id($exon) . " "  . $exon->start . " " . $exon->end . " with " . get_exon_id($comp_exon) .  "  " .$comp_exon->start . "  " . $comp_exon->end . "\n";
              if ($exon->start != $comp_exon->start ||
                  $exon->end != $comp_exon->end ||
                  (!$use_intron && $exon->strand != $comp_exon->strand) ) {
                # print "Failed \n";
               print "Comparing " . get_exon_id($exon) . " "  . $exon->start . " " . $exon->end . " with " . get_exon_id($comp_exon) .  "  " .$comp_exon->start . "  " . $comp_exon->end . "\n";
                next COMPTRANS;
              }
              if ($usephase) {
                if (($exon->phase != $comp_exon->phase) && !$first_e){
                  print FP "Phase diff (e_count = $e_count) for " . get_exon_id($exon) . " (" . $exon->phase . ") and exon " . get_exon_id($comp_exon) . "(" . $comp_exon->phase . ")\n";
                  next COMPTRANS;
                }
                if (($exon->end_phase != $comp_exon->end_phase) && !$first_e && !$last_e) {
                  print FP "End phase diff for " . get_exon_id($exon) . " (" . $exon->end_phase . ") and exon " . get_exon_id($comp_exon) . "(" . $comp_exon->end_phase . ")\n";
                  next COMPTRANS;
                }
              }
              $e_count++;
            }
            my $is_complete = is_complete_cds($trans);
            my $has_stops   = contains_stops($trans);
            my $has_frameshifts = has_frameshift_exons($trans);
            my $has_noncons_splice = (Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils::are_splice_sites_canonical($trans,$translate) == 1) ? 0 : 1;

            my $is_selcys = 0;
            if (exists($selcys_hash{get_transcript_id($comp_trans)})) {
              $is_selcys = 1;
              $has_stops = 0;
            }

            if (exists($selcys_hash{get_transcript_id($trans)})) {
              $is_selcys = 1;
              $has_stops = 0;
            }

            printf FP "MATCH:\t%s\t%-40.40s\t%-40.40s\t%-40.40s\t%-40.40s\t%1d\t%1d\t%1d\t%1d\t%3d\t%3d\t%3d",
                  $chr, get_gene_id($gene), get_transcript_id($trans),
                  get_gene_id($comp_gene), get_transcript_id($comp_trans),
                  $is_complete,  $has_stops,
                  $has_frameshifts, $has_noncons_splice,
                  scalar(@{$trans->get_all_Exons}),
                  scalar(@{$comp_trans->get_all_Exons}),
                  count_coding_exons($trans);

            $hadmatch{$trans} = 1;

            if ($is_complete && !$has_stops && !$has_frameshifts) {
              if (!$has_noncons_splice) {
                print FP " REAL";
              }
              print FP " GOOD";
              if ($is_selcys) {
                print FP " SELCYS";
              }
              $n_chr_complete_translating_match_trans++;
              $n_complete_translating_match_trans++;
              if ($outdb) {
                if (!exists($outgenes{$gene})) {
                  my $outgene = copy_empty_gene($gene);
                  $outgenes{$gene} = $outgene;
                }
                my $outgene = $outgenes{$gene};

                my $outtrans = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils::clone_Transcript($trans);
                $outgene->add_Transcript($outtrans);
              }
            }
            print FP "\n";

            $n_match_trans++;
            $n_chr_match_trans++;
            $cluster_match=1;
          }
        }
      }
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        if (!$hadmatch{$trans}) {
          my $is_complete = is_complete_cds($trans);
          my $has_stops   = contains_stops($trans);
          my $has_frameshifts = has_frameshift_exons($trans);
          my $has_noncons_splice = (Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils::are_splice_sites_canonical($trans,$translate) == 1) ? 0 : 1;

          my $is_ccds = 0;
          foreach my $xref (@{$trans->get_all_DBLinks}) {
            if ($xref->dbname eq 'CCDS') {
              $is_ccds=1;
              last;
            }
          }

          print FP "UNMATCHED TRANSCRIPT: ", get_transcript_id($trans) , " " , $trans->slice->seq_region_name , " " ,
                $trans->start, " ", $trans->end, " ",
                $is_complete,  " ",
                $has_stops, " ",
                $has_frameshifts, " ",
                $has_noncons_splice, " ",
                $is_ccds, " ",
                get_gene_id($gene), "\n";
        }
      }
    }
    # Temporary to see if in CCDS
    my $in_ccds = 0;
    # Temporary to see if in CCDS
    ##foreach my $gene (@comp_genes) {
    foreach my $gene (@genes) {
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        foreach my $xref (@{$trans->get_all_DBLinks}) {
          if ($xref->dbname eq 'CCDS') {
            $in_ccds++;
          }
        }
      }
    }

    if ($cluster_match) {
      $n_match_cluster++;
      $n_chr_match_cluster++;
      print FP "MATCHED ";
    } else {
      print FP "UNMATCHED ";
    }
    print FP "CLUSTER: Chr $chr";

    my $start = $genes[0]->start;
    my $end   = $genes[0]->end;
    my $nclusttrans = 0;
    my $nclustcodingtrans = 0;
    my $nclustfulllentrans = 0;
    foreach my $gene (@genes) {
      if ($gene->start < $start) {
        $start = $gene->start;
      }
      if ($gene->end > $end) {
        $end = $gene->end;
      }
      $nclusttrans += scalar(@{$gene->get_all_Transcripts});
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        if ($trans->translation) {
          $nclustcodingtrans++;
        }
        if (is_complete_cds($trans)) {
          $nclustfulllentrans++;
        }
      }
    }
    print FP " $start $end";

    $start = $comp_genes[0]->start;
    $end   = $comp_genes[0]->end;
    my $nclustcomptrans = 0;
    my $nclustcompcodingtrans = 0;
    my $nclustcompfulllentrans = 0;
    foreach my $gene (@comp_genes) {
      if ($gene->start < $start) {
        $start = $gene->start;
      }
      if ($gene->end > $end) {
        $end = $gene->end;
      }
      $nclustcomptrans += scalar(@{$gene->get_all_Transcripts});
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        if ($trans->translation) {
          $nclustcompcodingtrans++;
        }
        if (is_complete_cds($trans)) {
          $nclustcompfulllentrans++;
        }
      }
    }
    print FP " $start $end";
    print FP " $in_ccds " . scalar(@genes) . " $nclusttrans $nclustcodingtrans $nclustfulllentrans " . scalar(@comp_genes) . " $nclustcomptrans $nclustcompcodingtrans $nclustcompfulllentrans";

    my $firstone = 1;
    foreach my $gene (@genes) {
      if ($firstone) {
        print FP " ";
      } else {
        print FP ","
      }
      $firstone = 0;
      print FP get_gene_id($gene);
    }
    print FP " with ";
    $firstone = 1;
    foreach my $comp_gene (@comp_genes) {
      if ($firstone) {
        print FP " ";
      } else {
        print FP ","
      }
      $firstone = 0;
      print FP get_gene_id($comp_gene);
    }
    print FP "\n";
  }

  printf FP "SUMMARY: Chr %10.10s %4d genes %4d trans %4d comp_genes %4d comp_trans %4d clusters %4d unclustered %4d twoways %4d matched %4d good matches %4d matched_clusters\n",
                 $chr,scalar(@genes),$ntrans,scalar(@comp_genes),$ncomptrans,
                 scalar(@$clusters),scalar(@$unclustered),scalar(@twoways),
                 $n_chr_match_trans, $n_chr_complete_translating_match_trans,
                 $n_chr_match_cluster;
  if ($outdb) {
    my $out_slice;
    if ($schema == 20) {
      $slicename = "$coordsystem:$path:$chr:$chrstart:$chrend:1";
      $out_slice = $outsa->fetch_by_name($slicename);
    } else {
      $out_slice = $outsa->fetch_by_chr_start_end($chr,$chrstart,$chrend);
    }
    write_genes($outga, $out_slice, \%outgenes);
  }
}
print FP "OVERALL: Total matched trans $n_match_trans\n";
print FP "OVERALL: Total good trans matches $n_complete_translating_match_trans\n";
print FP "OVERALL: Total $set1_name genes $n_gene\n";
print FP "OVERALL: Total $set2_name genes $n_comp_gene\n";
print FP "OVERALL: Total $set1_name transcripts $n_trans\n";
print FP "OVERALL: Total $set2_name transcripts $n_comp_trans\n";
print FP "OVERALL: Total clusters $n_cluster\n";
print FP "OVERALL: Total matched clusters $n_match_cluster\n";
print FP "OVERALL: Total twoway clusters $n_twoway_cluster\n";
print FP "OVERALL: Total $set1_name only    clusters $n_set1_only_cluster\n";
print FP "OVERALL: Total $set2_name only    clusters $n_set2_only_cluster\n";
print FP "OVERALL: Total unclustered $n_unclustered\n";
print FP "OVERALL: Total $set1_name only unclustered $n_set1_only_unclustered\n";
print FP "OVERALL: Total $set2_name only unclustered $n_set2_only_unclustered\n";
print "Done\n";
close(FP);

sub get_gene_id {
  my ($gene) = @_;

  my $name = undef;

  if ($display_name) {
    if ($gene->external_name) {
      $name = $gene->external_name . " (" . $gene->stable_id . ")";
#      $name = $gene->external_name;
    }
  }

  $name = $gene->stable_id if ($gene->stable_id && !$name);

  if (!$name) {
    if (scalar(@{$gene->get_all_Transcripts}) == 1) {
      my $transcript = $gene->get_all_Transcripts->[0];

      my @tsfs = @{$transcript->get_all_supporting_features};
      if (scalar(@tsfs)) {
       $name = $tsfs[0]->hseqname;
      }
    } else {
     #print "FAILED To get tsf from gene with " . scalar(@{$gene->get_all_Transcripts}) . " transcripts, geneid = " . $gene->dbID . " ";
    }
  }
  if (!$name) {
    $name = $gene->dbID;
  }

  if ($show_gene_type) {
    $name .= ' '.$gene->biotype;
    $name .= '_'.$gene->status if ($gene->status);
    $name .= '_'.$gene->source;
  }
  return $name;
}

sub get_transcript_id {
  my ($transcript) = @_;

  if ($display_name) {
    if ($transcript->external_name) {
  #    return $transcript->external_name;
    }
  }

  return $transcript->stable_id if ($transcript->stable_id);

  my @tsfs = @{$transcript->get_all_supporting_features};
  if (scalar(@tsfs)) {
    return $tsfs[0]->hseqname;
  }

  return $transcript->dbID;
}

sub get_exon_id {
  my ($exon) = @_;

  return $exon->stable_id if ($exon->stable_id);

  return $exon->dbID;
}

sub write_genes {
  my ($ga, $slice, $geneshash)=@_;

  foreach my $gene ( values %$geneshash ) {
    my $tcount = 1;
    foreach my $trans ( @{$gene->get_all_Transcripts} ) {
      my $get2 = $trans->translation->stable_id;
      trim_to_CDS($trans,$tcount);

      my @exons= @{$trans->get_all_Exons};
      my $get  = $trans->translation;
      $trans->_translation_id(undef);

      foreach my $exon (@exons) {
        $exon->stable_id;
        $exon->contig($slice);
        $exon->get_all_supporting_features;
      }
      $tcount++;
    }
    prune_Exons($gene);
    $gene->transform;
    $ga->store($gene);
  }
}

sub trim_to_CDS {
  my ($trans,$tnum) = @_;

  my @exons=tweak_translateable($trans,$tnum);

  $trans->flush_Exons;

  foreach my $exon (@exons) {
    $trans->add_Exon($exon);
  }

  $trans->translation->start_Exon($exons[0]);
  $trans->translation->end_Exon($exons[$#exons]);

  $trans->translation->start(1);
  $trans->translation->end($exons[$#exons]->length);
}

sub tweak_translateable {
  my ($trans,$tnum) = @_;

  my @exons = @{$trans->get_all_translateable_Exons};

# The start and end exons may have been clipped
# This hack gives them the stable_id of the unclipped exon

  if ($trans->translation->start_Exon != $exons[0]) {
    $exons[0]->stable_id($trans->translation->start_Exon->stable_id . "S$tnum");
    $exons[0]->version($trans->translation->start_Exon->version);
    $exons[0]->created($trans->translation->start_Exon->created);
    $exons[0]->modified($trans->translation->start_Exon->modified);

    $exons[0]->add_supporting_features(@{$trans->translation->start_Exon->get_all_supporting_features});
  }

  if ($trans->translation->start_Exon != $exons[$#exons]) {
    $exons[$#exons]->stable_id($trans->translation->end_Exon->stable_id . "E$tnum");
    $exons[$#exons]->version($trans->translation->end_Exon->version);
    $exons[$#exons]->created($trans->translation->end_Exon->created);
    $exons[$#exons]->modified($trans->translation->end_Exon->modified);

    # If not same as start exon
    if ($#exons) {
      $exons[$#exons]->add_supporting_features(@{$trans->translation->end_Exon->get_all_supporting_features});
    }
  }

  return @exons;
}


sub copy_empty_gene {

  print ">>> copy_empty_gene\n";
  my ($gene) = @_;
  my $newgene = new Bio::EnsEMBL::Gene;
  if ($gene->type){
    $newgene->type( $gene->biotype);
  }
  if ( defined $gene->analysis ){
    $newgene->analysis($gene->analysis);
  }
  if ( defined $gene->stable_id ){
    $newgene->stable_id( $gene->stable_id );
    $newgene->version( $gene->version );
    $newgene->modified( $gene->modified );
    $newgene->created( $gene->created );
  }

  return $newgene;
}

sub remove_redundant_cds_from_gene {
  my ($gene) = @_;
  my @trans_array = @{$gene->get_all_Transcripts};
  my @trans_to_remove;

  for (my $i=0; $i<scalar(@trans_array); $i++) {
    my $trans = $trans_array[$i];

    next if (!defined($trans));
    next if (!$trans->translation);

    my @exons = sort {$a->start <=> $b->start} @{$trans->get_all_translateable_Exons};

    COMPTRANS:
    for (my $j=$i+1; $j<scalar(@trans_array); $j++) {
      my $comp_trans = $trans_array[$j];

      next if (!defined($comp_trans));
      next if (!$comp_trans->translation);
      next if ($trans == $comp_trans);

      my @comp_exons = sort {$a->start <=> $b->start} @{$comp_trans->get_all_translateable_Exons};

      next if (scalar(@comp_exons) != scalar(@exons));

      foreach my $exon (@exons) {
        my $comp_exon = shift @comp_exons;
        if ($exon->start != $comp_exon->start ||
            $exon->end != $comp_exon->end ||
            $exon->strand != $comp_exon->strand) {
          next COMPTRANS;
        }
      }
      $trans_array[$j] = undef;
      push @trans_to_remove, $comp_trans;
    }
  }

  foreach my $trans (@trans_to_remove) {
    remove_transcript_from_gene($gene,$trans);
  }

  return scalar(@trans_to_remove);
}

sub remove_redundant_cds_from_cluster {
  my ($cluster, $set) = @_;
  my @genes = @{$cluster->get_Genes_by_Set($set)};
  my $n_removed = 0;

  foreach my $gene (@genes) {
    $n_removed += remove_redundant_cds_from_gene($gene);
  }


# So now we have a set of genes which have no duplication within themselves

  my @unique_trans_array;

  foreach my $gene (@genes) {
    my @trans_from_gene = @{$gene->get_all_Transcripts};

    TRANS:
    for (my $i=0; $i<scalar(@trans_from_gene); $i++) {
      my $trans = $trans_from_gene[$i];

      next if (!$trans->translation);

      my @exons = sort {$a->start <=> $b->start} @{$trans->get_all_translateable_Exons};

      COMPTRANS:
      for (my $j=0; $j<scalar(@unique_trans_array); $j++) {
        my $comp_trans = $unique_trans_array[$j];

        next if (!$comp_trans->translation);

        my @comp_exons = sort {$a->start <=> $b->start} @{$comp_trans->get_all_translateable_Exons};

        next if (scalar(@comp_exons) != scalar(@exons));

        foreach my $exon (@exons) {
          my $comp_exon = shift @comp_exons;
          if ($exon->start != $comp_exon->start ||
              $exon->end != $comp_exon->end ||
              $exon->strand != $comp_exon->strand) {
            next COMPTRANS;
          }
        }
        print "NOTE: Removing extra trans " . get_transcript_id($trans) . " from cluster\n";
        $n_removed++;
        remove_transcript_from_gene($gene,$trans);
        next TRANS;
      }
    }
    push @unique_trans_array,@{$gene->get_all_Transcripts};
  }

  return $n_removed;
}

sub remove_transcript_from_gene {
  my ($gene, $trans_to_del)  = @_;

  my @newtrans;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans != $trans_to_del) {
      push @newtrans,$trans;
    }
  }

# The naughty bit!
  $gene->{_transcript_array} = [];

  print "REM removed " . get_transcript_id($trans_to_del) . "\n";

  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }

  return scalar(@newtrans);
}


sub count_trans {
  my ($genes) = @_;
  my $count = 0;

  foreach my $gene (@$genes) {
    $count += scalar(@{$gene->get_all_Transcripts});
  }

  return $count;
}

sub has_frameshift_exons {
  my ($t) = @_;

  my @exons = @{$t->get_all_Exons};

  if ( scalar (@exons ) == 1 ){
    return 0;

  } elsif (scalar(@exons) > 1) {
    @exons = sort{ $a->start <=> $b->start } @exons;

    for(my $i=0; $i<$#exons; $i++){
      my $intlen = $exons[$i+1]->start - $exons[$i]->end - 1;
      if ( $intlen < 10 ){
        return 1;
      }
    }
    return 0;

  } else{
    return 0;
  }
# Ensure the world sees this the right way round
  #$t->sort;
}

sub count_coding_exons {
  my ($t) = @_;

  return 0 if (!$t->translation);

  return (scalar(@{$t->get_all_translateable_Exons}));
}

sub gene_is_all_partial {
  my ($gene) = @_;

  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if (is_complete_cds($trans)) {
      return 0;
    }
  }
  return 1;
}


sub is_significant_overlap {
  my ($gene1,$gene2,$sigamount) = @_;

  my @genes1 = ($gene1);
  my @genes2 = ($gene2);
  my $tc1 = genes_to_Transcript_Cluster(\@genes1);
  my $tc2 = genes_to_Transcript_Cluster(\@genes2);

  my @exonclusters1 = cluster_exons_in_transcript_cluster($tc1);
  my @exonclusters2 = cluster_exons_in_transcript_cluster($tc2);

  my $overlap = 0;
  foreach my $cluster1 (@exonclusters1) {
    foreach my $cluster2 (@exonclusters2) {
      my ($a,$over,$b) = $cluster1->overlap_extent($cluster2);
      $overlap += $over;
    }
  }
  print "Overlap = $overlap\n";
  if ($overlap >= $sigamount) {
    return 1;
  } else {
    return 0;
  }
}

sub is_complete_cds {
  my ($trans) = @_;

  return 0 if (!defined($trans->translation));

  my $tln = $trans->translation;

  #$trans->sort;

  my $coding_start = $trans->cdna_coding_start;
  my $coding_end   = $trans->cdna_coding_end;
  my $cdna_seq     = uc($trans->spliced_seq);

  my $startseq     = substr($cdna_seq,$coding_start-1,3);
  my $endseq       = substr($cdna_seq,$coding_end-3,3);

#  print "Codons: " . $startseq . " " . $endseq . "\n";

  return 0 if ($startseq ne "ATG");
  return 0 if ($endseq ne "TAG" && $endseq ne "TGA" && $endseq ne "TAA");

  return 1;
}

sub contains_stops {
  my ($trans) = @_;

  if (!defined($trans->translation)) {
    print "Warning: Called contains_stops with non translating transcript\n";
    return 0;
  }

  my $tln = $trans->translation;

  #$trans->sort;

  my $pepseq = $trans->translate->seq;

  return 1 if ($pepseq =~ /\*/);

  return 0;
}

sub check_for_gene_duplication {
  my ($genes) = @_;

  my $tc = genes_to_Transcript_Cluster($genes);

  my @exonclusters = cluster_exons_in_transcript_cluster($tc);

  my $most_exon_transcript = find_Transcript_with_most_Exons_in_Transcript_Cluster($tc);

  my $max_num_exons_transcript = scalar(@{$most_exon_transcript->get_all_Exons()});

  print "Number of exon clusters = " . scalar(@exonclusters) .
        " max number of exons in transcript = " . $max_num_exons_transcript . "\n";

# Matches against multiple wrongly combined genes (duplicates or gene cluster)
  if (scalar(@exonclusters) >= 2.0 * $max_num_exons_transcript) {
    print "NOTE: Possible gene duplication\n";
  }
}

sub find_Longest_CDS_of_Gene {
  my $gene = shift;

  my $maxlen = 0;
  my $nexonmax = 0;
  my $longest;
  print "For gene " . get_gene_id($gene);
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    my $len = 0;
    if ($trans->translation) {
      my @trans_exons = @{$trans->get_all_translateable_Exons()};
      foreach my $exon (@trans_exons) {
  #      print "exon start = " . $exon->start . " end = " . $exon->end . " length " . $exon->length . "\n";
        $len+= $exon->length;
      }
  #    print "Len = $len\n";
      if ($len > $maxlen) {
  #      print "New max len = $len\n";
        $maxlen = $len;
        $longest = $trans;
        $nexonmax = scalar($trans->get_all_Exons)
      }
    }
  }
  if (!defined($longest)) {
    print "Didn't find a longest transcript for gene "
      . get_gene_id($gene) . " type " . $gene->biotype . "\n";
  } else {
    print " longest transcript is " . get_transcript_id($longest) . "\n";
  }
  return ($longest,$maxlen,$nexonmax);
}



sub cluster_Genes {
  my ($genes, $types_hash) = @_;

  return ([],[]) if (!scalar(@$genes));

  my @sorted_genes = sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end } @$genes;

  print "Clustering ".scalar( @sorted_genes )." genes...\n";

  my @clusters;

  foreach my $gene (@sorted_genes) {
    my @matching_clusters;
  CLUSTER:
    foreach my $cluster (@clusters) {
      if ($gene->end  >= $cluster->start &&
          $gene->start <= $cluster->end) {
        foreach my $cluster_gene (@{$cluster->get_Genes}){
          if ($gene->end  >= $cluster_gene->start &&
              $gene->start <= $cluster_gene->end) {

            if (_compare_Genes( $gene, $cluster_gene)) {
              push (@matching_clusters, $cluster);
              next CLUSTER;
            }
          }
        }
      }
    }

    if (scalar(@matching_clusters) == 0) {
      my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster->new();
      foreach my $set_name (keys %$types_hash) {
        $newcluster->gene_Types($set_name,$types_hash->{$set_name});
      }
      $newcluster->put_Genes([$gene]);

      push(@clusters,$newcluster);

    } elsif (scalar(@matching_clusters) == 1) {
      $matching_clusters[0]->put_Genes([$gene]);

    } else {
      # Merge the matching clusters into a single cluster
      my @new_clusters;

      my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster->new();

      foreach my $set_name (keys %$types_hash) {
        $merged_cluster->gene_Types($set_name,$types_hash->{$set_name});
      }

      my %match_cluster_hash;
      foreach my $clust (@matching_clusters) {
        $merged_cluster->put_Genes($clust->get_Genes);

        $match_cluster_hash{$clust} = $clust;
      }
      $merged_cluster->put_Genes([$gene]);

      push @new_clusters,$merged_cluster;

      # Add back non matching clusters
      foreach my $clust (@clusters) {
        if (!exists($match_cluster_hash{$clust})) {
          push @new_clusters,$clust;
        }
      }
      @clusters =  @new_clusters;
    }
  }

  my @new_clusters;
  my @unclustered;
  foreach my $cl (@clusters){
    if ( $cl->get_Gene_Count == 1 ){
      push @unclustered, $cl;
    }
    else{
      push( @new_clusters, $cl );
    }
  }

  return (\@new_clusters, \@unclustered);
}

=head2 _compare_Genes()

 Title: _compare_Genes
 Usage: this internal function compares the exons of two genes on overlap

=cut

sub _compare_Genes {
  my ($gene1,$gene2) = @_;

#  print "Comparing " . $gene1->stable_id . " to " . $gene2->stable_id . "\n";
#
  if ($gene1->end < $gene2->start || $gene1->start > $gene2->end) {
#    print "Gene 1  " . $gene1->start . " " . $gene1->end . " \n";
#    print "Gene 2  " . $gene1->start . " " . $gene1->end . " \n";
#    print "Failed extents check - returning 0\n";
    return 0;
  }

  if (!$translate) {
    foreach my $exon1 (@{$gene1->get_all_Exons}){
      foreach my $exon2 (@{$gene2->get_all_Exons}){
        if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
          #print "Passed all exon overlap check - returning 1\n";
  	  return 1;
        }
      }
    }
  } else {
    my $exons1 = get_coding_exons_for_gene($gene1);
    my $exons2 = get_coding_exons_for_gene($gene2);
    foreach my $exon1 (@$exons1) {
      foreach my $exon2 (@$exons2) {
        if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
          #print "Passed CDS overlap check - returning 1\n";
  	  return 1;
        }
      }
    }
  }
  #print "Failed overlap check (translate = $translate) - returning 0\n";
  return 0;
}
#########################################################################

{
  my %coding_exon_cache;

  sub clear_coding_exons_cache {
    %coding_exon_cache = ();
  }

  sub get_coding_exons_for_gene {
    my ($gene) = @_;

    if (exists($coding_exon_cache{$gene})) {
      return $coding_exon_cache{$gene};
    } else {
      my %coding_hash;
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        next if (!$trans->translation);
        foreach my $exon (@{$trans->get_all_translateable_Exons}) {
          $coding_hash{$exon} = $exon;
        }
      }

      # my @coding = sort { $a->start <=> $b->start } values %coding_hash;
      my @coding = values %coding_hash;

      $coding_exon_cache{$gene} = \@coding;
      return $coding_exon_cache{$gene};
    }
  }
}

sub prune_Exons {
  my ($gene) = @_;

  my @unique_Exons;

  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_exon and translation->end_Exon

  #print STDERR "Pruning exons\n";

  my %exonhash;

  foreach my $tran (@{ $gene->get_all_Transcripts }) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty

      UNI: foreach my $uni (@unique_Exons) {
        if ($uni->start  == $exon->start  &&
            $uni->end    == $exon->end    &&
            $uni->strand == $exon->strand &&
            $uni->phase  == $exon->phase   &&
            $uni->end_phase == $exon->end_phase)
        {
          $found = $uni;
          last UNI;
        }
      }
        #print STDERR " Exon " . $exon->stable_id . "\n";
        #print STDERR " Phase " . $exon->phase . " EndPhase " . $exon->end_phase . "\n";
        #print STDERR " Strand " . $exon->strand . " Start " . $exon->start . " End ". $exon->end ."\n";

      if (defined($found)) {
        #print STDERR " Duplicate\n";
        push (@newexons, $found);
        if ($tran->translation) {
          if ($exon == $tran->translation->start_Exon) {
            $tran->translation->start_Exon($found);
          }

          if ($exon == $tran->translation->end_Exon) {
            $tran->translation->end_Exon($found);
          }
        }
      } else {
        #print STDERR "New = " . $exon->stable_id . "\n";

        ### This is nasty for the phases - sometimes exons come back with
        ### the same stable id and different phases - we need to strip off
        ### the stable id if we think we have a new exon but we've
        ### already seen the stable_id

        if (defined($exon->stable_id) && defined($exonhash{$exon->stable_id})) {
           #print STDERR "Already seen stable id " . $exon->stable_id . " - removing stable_id\n";
           $exon->{_stable_id} = undef;
           #print STDERR "Exon id " .$exon->stable_id . "\n";
        }
        push (@newexons,     $exon);
        push (@unique_Exons, $exon);
      }
      if (my $stable = $exon->stable_id) {
        $exonhash{$stable} = 1;
      }
    }
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }

  my @exons = @{$gene->get_all_Exons};

  %exonhash = ();

  foreach my $ex (@exons) {
    if (my $stable = $ex->stable_id) {
      $exonhash{$stable}++;
    }
  }

  while (my ($id, $count) = each %exonhash) {
    print STDERR "Exon id $id seen $count times\n" if $count > 1;
  }
}

