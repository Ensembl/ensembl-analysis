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
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DnaPepAlignFeature;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $test,
    $verbose,
    $gene_logic_name,
    $sf_logic_name,
    %ug_feats, 
    %transcripts,
    %attributes,
    @out_genes,
);

$dbport = 3306;

# Load the attributes by default
my $load_attributes = 1;

GetOptions(
            'dbname|db|D=s'     => \$dbname,
            'dbuser|user|u=s'     => \$dbuser,
            'dbhost|host|h=s'     => \$dbhost,
            'dbport|port|P=s'     => \$dbport,
            'dbpass|pass|p=s'     => \$dbpass,
            'genelogic=s'  => \$gene_logic_name,
            'sflogic=s'    => \$sf_logic_name,
            'test'         => \$test,
            'verbose'      => \$verbose,
            # Use -noattributes if not loading the attributes to the database
            'attributes!'   => \$load_attributes,
            );

die "You must supply a gene/transcript logic name with -genelogic\n" 
    if not defined $gene_logic_name;
die "You must supply a supporting_feature logic name with -sflogic\n" 
    if not defined $sf_logic_name;


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                             '-dbname' => $dbname,
                                             '-host' => $dbhost,
                                             '-user' => $dbuser,
                                             '-port' => $dbport,
                                             '-pass' => $dbpass);

my $gene_analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($gene_logic_name);

if (not defined $gene_analysis) {
  $verbose and print STDERR "Storing new entry in analysis for WGA2Genes gene\n";

  $gene_analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $gene_logic_name,
                                               -module => 'WGA2Genes',
                                               -gff_source => 'ensembl',
                                               -gff_feature => 'gene');
  $db->get_AnalysisAdaptor->store($gene_analysis);
}

my $sf_analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($sf_logic_name);

if (not defined $sf_analysis) {
  $verbose and print STDERR "Storing new entry in analysis for WGA2Genes sup.feats\n";

  $sf_analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $sf_logic_name,
                                             -module => 'WGA2Genes',
                                             -gff_source => 'ensembl',
                                             -gff_feature => 'feature');
  $db->get_AnalysisAdaptor->store($sf_analysis);
}


my ($gs_cs, $prev_tl_cs) = 
    sort {$a->rank <=> $b->rank} @{$db->get_CoordSystemAdaptor->fetch_all};

$verbose and print STDERR "Gene scaffold coord-system = '" . $gs_cs->name . "'\n";
$verbose and print STDERR "Previous top level coord-system = '" . $prev_tl_cs->name . "'\n";

$verbose and print STDERR "Reading GENE file...\n";

while(<>) {
  #
  # Pseudo GFF. Fields:
  #
  # seq_name
  # source
  # type
  # start
  # end
  # strand
  # 
  # Then for exons we have
  #
  # phase
  # end_phase
  # group
  # 
  # And for supporting features we have
  # 
  # group
  #

  if (/^\#\#\-ATTRIBUTE\s+(.+)/) {
    my $line = $1;
    my ($tid, $code, $name, $val, $desc);
    if ($line =~ /transcript=(\S+)/) {
      $tid = $1;
    }
    if ($line =~ /code=(\S+)/) {
      $code = $1;
    }
    if ($line =~ /value=(\S+)/) {
      $val = $1;
    }
    if ($line =~ /name=\"([^\"]+)\"/) {
      $name = $1;
    }
    if ($line =~ /desc=\"([^\"]+)\"/) {
      $desc = $1;
    }

    if (defined($tid) and defined($code) and defined($val)) {
      $name = $code if not defined $name;
      $desc = $code if not defined $desc;

      push @{$attributes{$tid}}, {
        code => $code,
        value => $val,
        name => $name,
        desc => $desc,
      };
    }
    next;
  } elsif (/^\#/) {
    next;
  }

  my @l = split(/\t/, $_);

  my ($seq_name,
      $source,
      $type,
      $start,
      $end,
      $strand) = ($l[0], $l[1], $l[2], $l[3], $l[4], $l[5]);


  if ($l[2] eq 'Exon') {
    my ($phase, $end_phase) = ($l[6], $l[7]);

    my %group = split /[\s|;=]+/, $l[8];

    my ($exon_id)       = $group{exon};
    my ($transcript_id) = $group{transcript};

    push @{$transcripts{$seq_name}->{$transcript_id}}, {
      exon_id => $exon_id,
      start   => $start,
      end     => $end,
      strand  => $strand,
      phase   => $phase,
      end_phase => $end_phase,
    };

  } elsif ($l[2] eq 'Supporting') {
    my %group = split /[\s|;=]+/, $l[6];

    my ($hit_name) = $group{hname};
    my ($hit_start) = $group{hstart};
    my ($hit_end) = $group{hend};
    my ($exon_id) = $group{exon};
    my ($external_db_id) = $group{external_db_id};
    my ($hcoverage) = $group{hcoverage};

    push @{$ug_feats{$exon_id}}, {
      start => $start,
      end   => $end,
      strand => $strand,
      hseqname => $hit_name,
      hstart => $hit_start,
      hend  => $hit_end,
      external_db_id => $external_db_id,
      hcoverage => $hcoverage,
    };
  }
}



$verbose and printf(STDERR "Constructing transcripts for %d seq_regions...\n", 
                    scalar(keys %transcripts));

my $sr_count = 0;

foreach my $target_id (keys %transcripts) {
  my ($slice, @transcripts);

  if ($target_id =~ /^GeneScaffold/) {
    $slice = $db->get_SliceAdaptor->fetch_by_region($gs_cs->name,
                                                    $target_id);
  } else {
    $slice = $db->get_SliceAdaptor->fetch_by_region($prev_tl_cs->name,
                                                    $target_id);
  }

  foreach my $transcript_id (keys %{$transcripts{$target_id}}) {
    my $tran = Bio::EnsEMBL::Transcript->new(-analysis => $gene_analysis,
                                             -biotype  => 'protein_coding');

    my (@tran_ug_feats);

    foreach my $ex_h (@{$transcripts{$target_id}->{$transcript_id}}) {
      my $exon = Bio::EnsEMBL::Exon->new(-start     => $ex_h->{start},
                                         -end       => $ex_h->{end},
                                         -strand    => $ex_h->{strand},
                                         -phase     => $ex_h->{phase},
                                         -end_phase => $ex_h->{end_phase},
                                         -slice     => $slice);
      if (exists($ug_feats{$ex_h->{exon_id}})) {
        my @exon_ug_feats;
        foreach my $f_hash (@{$ug_feats{$ex_h->{exon_id}}}) {
          my $f = Bio::EnsEMBL::FeaturePair->new(-slice    => $slice,
                                                 -start    => $f_hash->{start},
                                                 -end      => $f_hash->{end},
                                                 -strand   => $f_hash->{strand},
                                                 -hseqname => $f_hash->{hseqname},
                                                 -hstart   => $f_hash->{hstart},
                                                 -hend     => $f_hash->{hend},
                                                 -external_db_id => $f_hash->{external_db_id},
                                                 -hcoverage => $f_hash->{hcoverage},
                                                 -analysis => $sf_analysis);

          push @exon_ug_feats, $f;
          push @tran_ug_feats, $f;
        }

        my $sf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@exon_ug_feats);
        $exon->add_supporting_features($sf);
      }

      $tran->add_Exon($exon);
    }

    if (@tran_ug_feats) {
      my $sf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tran_ug_feats);
      $tran->add_supporting_features($sf);
    }

    # Tranlation

    my $translation = Bio::EnsEMBL::Translation->new;
    $translation->start_Exon($tran->get_all_Exons->[0]);
    $translation->end_Exon($tran->get_all_Exons->[-1]);
    $translation->start(1);
    $translation->end($tran->get_all_Exons->[-1]->length);
    $tran->translation($translation);

    # Attributes
    if ($load_attributes) {
      foreach my $attr_hash (@{$attributes{$transcript_id}}) {
        my $attr = Bio::EnsEMBL::Attribute->new(-code => $attr_hash->{code},
                                                -name => $attr_hash->{name},
                                                -desc => $attr_hash->{desc},
                                                -value => $attr_hash->{value});
        $tran->add_Attributes($attr);
      }
    }

    push @transcripts, $tran;
  }

  foreach my $clust (&cluster_by_exon_overlap(@transcripts)) {
    # check for exon-overlap between result genes;
    my $gene = Bio::EnsEMBL::Gene->new(-analysis => $gene_analysis,
                                       -biotype  => 'protein_coding');
    foreach my $tran (@$clust) {
      $gene->add_Transcript($tran);
    }

    if ($verbose and scalar(@$clust) > 1) {
      print STDERR "Found a transcript cluster on $target_id\n";
    }
  
    push @out_genes, $gene;
  }

  $sr_count++;
  if ($verbose and $sr_count % 1000 == 0) {
    print STDERR "Made transcripts for $sr_count seq_regions...\n";
  }
}

if ($test) {
  foreach my $g (@out_genes) {
    print &gene_string($g);
  }
} else {
  $verbose and print STDERR "Writing genes...\n";

  my $current_gene;

  eval {
    my $gene_count = 0;

    foreach my $g (@out_genes) {
      $current_gene = $g;
      $db->get_GeneAdaptor->store($g);

      $gene_count++;
      if ($verbose and $gene_count % 1000 == 0) {
        print STDERR "Written $gene_count genes...\n";
      }
    }
  };
  if ($@) {
    die "Failure during writing of gene:\n" . &gene_string($current_gene) . "\n";
  } else {
    print "Successfully wrote ", scalar(@out_genes), " transcripts\n";
  }
}


sub cluster_by_exon_overlap {
  my (@trans) = @_;

  my @clusters;
  foreach my $t (@trans) {
    my @exons = @{$t->get_all_Exons};

    my @overlapping;

    foreach my $c (@clusters) {
      my $overlap = 0;

      CT: foreach my $ct (@$c) {
        my @c_exons = @{$ct->get_all_Exons};

        for(my $i=0; $i<@exons; $i++) {
          for(my $j=0; $j < @c_exons; $j++) {
            if ($exons[$i]->strand == $c_exons[$j]->strand and
                $exons[$i]->start <= $c_exons[$j]->end and
                $exons[$i]->end >= $c_exons[$j]->start) {
              $overlap = 1;
              last CT;
            }
          }
        }
      }
      if ($overlap) {
        push @overlapping, $c;
      }
    }

    if (@overlapping) {
      my ($first, @others) = @overlapping;
      foreach my $o (@others) {
        push @$first, @$o;
        @$o = ();
      }
      push @$first, $t;
    } else {
      push @clusters, [$t];
    }
  }

  my @final_clusts;
  foreach my $c (@clusters) {
    if (@$c) {
      push @final_clusts, $c;
    }
  }

  return @final_clusts;
}


sub gene_string {
  my $g = shift;

  my $str = "";

  foreach my $t (@{$g->get_all_Transcripts}) {
    $str .= "TRANSCRIPT\n";
    foreach my $e (@{$t->get_all_Exons}) {
      $str .= sprintf("%s %d %d %d %d %d\n", 
                      $e->slice->seq_region_name, 
                      $e->start, 
                      $e->end, 
                      $e->strand, 
                      $e->phase, 
                      $e->end_phase);
      foreach my $sf (@{$e->get_all_supporting_features}) {
        $str .=  sprintf(" SUPPORT: %d %d %s %d %d\n", 
                         $sf->start, 
                         $sf->end, 
                         $sf->hseqname, 
                         $sf->hstart, 
                         $sf->hend);
      }
    }
  }

  return $str;
}
