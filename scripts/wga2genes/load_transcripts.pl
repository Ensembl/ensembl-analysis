#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;

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
    $logic_name,
    %ug_feats, 
    %transcripts,
    %attributes,
    @out_genes,
);

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'logic=s'  => \$logic_name,
            'test'     => \$test,
            'verbose'  => \$verbose,
            );

die "You must supply a logic name with -logic_name\n" 
    if not defined $logic_name;

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                             '-dbname' => $dbname,
                                             '-host' => $dbhost,
                                             '-user' => $dbuser,
                                             '-port' => $dbport,
                                             '-pass' => $dbpass);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);

die "Analysis for logic name '$logic_name' does not exist in the database\n"
    if not defined $analysis;

$verbose and print STDERR "Reading GFF file...\n";

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

  next if $source ne $logic_name;

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
    
    push @{$ug_feats{$exon_id}}, {
      start => $start,
      end   => $end,
      strand => $strand,
      hseqname => $hit_name,
      hstart => $hit_start,
      hend  => $hit_end,
    };
  }
}

$verbose and print STDERR "Constructing transcripts...\n";

foreach my $target_id (keys %transcripts) {
  my $slice;
  if ($target_id =~ /^SCAFFOLD/) {
    $slice = $db->get_SliceAdaptor->fetch_by_region('scaffold',
                                                    $target_id);
  } else {
    $slice = $db->get_SliceAdaptor->fetch_by_region('genescaffold',
                                                    $target_id);
  }

  foreach my $transcript_id (keys %{$transcripts{$target_id}}) {
    my $tran = Bio::EnsEMBL::Transcript->new;

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
                                                 -analysis => $analysis);

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
    foreach my $attr_hash (@{$attributes{$transcript_id}}) {
      my $attr = Bio::EnsEMBL::Attribute->new(-code => $attr_hash->{code},
                                              -name => $attr_hash->{name},
                                              -desc => $attr_hash->{desc},
                                              -value => $attr_hash->{value});
      $tran->add_Attributes($attr);
    }

    my $gene = Bio::EnsEMBL::Gene->new;
    $gene->type($logic_name);
    $gene->analysis($analysis);

    $gene->add_Transcript($tran);

    push @out_genes, $gene;
  }
}

if ($test) {
  foreach my $g (@out_genes) {
    foreach my $t (@{$g->get_all_Transcripts}) {
      print "TRANSCRIPT\n";
      foreach my $e (@{$t->get_all_Exons}) {
        printf "%s %d %d %d %d %d\n", $e->slice->seq_region_name, $e->start, $e->end, $e->strand, $e->phase, $e->end_phase;
        foreach my $sf (@{$e->get_all_supporting_features}) {
          printf " SUPPORT: %d %d %s %d %d\n", $sf->start, $sf->end, $sf->hseqname, $sf->hstart, $sf->hend;
        }
      }
    }
  }
} else {
  $verbose and print STDERR "Writing genes...\n";

  eval {
    foreach my $g (@out_genes) {
      $db->get_GeneAdaptor->store($g);
    }
  };
  if ($@) {
    die "Something went wrong while writing genes: $@\n";
  } else {
    print "Successfully wrote ", scalar(@out_genes), " transcripts\n";
  }
}
