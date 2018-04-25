use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;

my $source_dbname = '';
my $source_user   = '';
my $source_host   = '';
my $source_port   = '';

my $dna_dbname = '';
my $dna_user   = '';
my $dna_host   = '';
my $dna_port   = '';

my $out_dbname = '';
my $out_user   = '';
my $out_host   = '';
my $out_port   = '';
my $out_pass   = '';


my $slice_name = ''; # e.g. 'primary_assembly:Astyanax_mexicanus-2.0:23:1:20000000:1'
my $select_most_supported_transcripts = 1;
my $source_type = 'protein';

my $options = GetOptions ("source_user=s"   => \$source_user,
                          "source_host=s"   => \$source_host,
                          "source_port=i"   => \$source_port,
                          "source_dbname=s" => \$source_dbname,
                          "out_host=s"      => \$out_host,
                          "out_port=i"      => \$out_port,
                          "out_dbname=s"    => \$out_dbname,
                          "out_user=s"      => \$out_user,
                          "out_pass=s"      => \$out_pass,
                          "dna_dbname=s"    => \$dna_dbname,
                          "dna_host=s"      => \$dna_host,
                          "dna_port=i"      => \$dna_port,
                          "dna_user=s"      => \$dna_user,
                          "slice_name=s"    => \$slice_name);

my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dna_port,
  -user    => $dna_user,
  -host    => $dna_host,
  -dbname  => $dna_dbname);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $source_port,
  -user    => $source_user,
  -host    => $source_host,
  -dbname  => $source_dbname);

my $output_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $out_port,
  -user    => $out_user,
  -pass    => $out_pass,
  -host    => $out_host,
  -dbname  => $out_dbname);

$db->dnadb($dna_db);
$output_db->dnadb($dna_db);

my $initial_genes;

if($slice_name) {
  my $slice_adaptor = $db->get_SliceAdaptor();
  my $slice = $slice_adaptor->fetch_by_name($slice_name);
  say "Processing slice: ".$slice->name;
  $initial_genes = $slice->get_all_Genes;
} else {
  say "No slice name provided, so fetching all genes";
  $initial_genes = $db->get_GeneAdaptor->fetch_all();
}

unless(scalar(@$initial_genes)) {
  throw("Found no genes in input db");
}

my $genes = [];
foreach my $initial_gene (@{$initial_genes}) {
  my $transcripts = $initial_gene->get_all_Transcripts();
  if(scalar(@$transcripts > 1)) {
    throw("Found a multi-transcript gene with dbID ".$initial_gene->dbID.". The script is built for single transcript genes");
  }

  my $initial_transcript = shift($transcripts);
  # Skip transcripts that look dodgy because of frameshift or small introns
  if($initial_transcript->biotype =~ /\_pe3\_/ || $initial_transcript->biotype =~ /\_del$/ || $initial_transcript->biotype =~ /old\_/ ||
     $initial_transcript->biotype =~ /\_sub\_/ || $initial_transcript->biotype =~ /\_int\_/) {
    next;
  } elsif(avg_intron_size($initial_transcript)) {
    next;
  } else {
    push(@$genes,$initial_gene);
  }
}

my $biotypes_hash = get_all_biotypes($genes);
my $biotypes_array = [keys(%$biotypes_hash)];

my $types_hash;
$types_hash->{genes} = $biotypes_array;

say "  Clustering genes from input_dbs...";
my ($clusters, $unclustered) = cluster_Genes($genes,$types_hash);
my $output_genes = process_clusters($clusters,$unclustered);

if($select_most_supported_transcripts) {
  my ($clusters, $unclustered) = cluster_Genes($output_genes,$types_hash);
  my $most_supported_output_genes = select_most_supported_transcripts($clusters,$unclustered);
  write_output($most_supported_output_genes,$output_db);
} else {
  write_output($output_genes,$output_db);
}
exit;


sub write_output {
  my ($output_genes,$output_db) = @_;

  my $gene_adaptor = $output_db->get_GeneAdaptor();
  foreach my $output_gene (@{$output_genes}) {
    empty_Gene($output_gene);
    $gene_adaptor->store($output_gene);
  }
}


sub get_all_biotypes {
  my ($genes) = @_;

  my $master_biotypes_hash = {};
  foreach my $gene (@{$genes}) {
    unless($master_biotypes_hash->{$gene->biotype}) {
      $master_biotypes_hash->{$gene->biotype} = 1;
    }
  }
  return($master_biotypes_hash);
}


sub process_clusters {
  my ($clustered,$unclustered) = @_;

  # This is expecting single transcript gene clusters
  my $all_clusters = [@{$clustered},@{$unclustered}];

  my $output_genes = [];
  foreach my $single_cluster (@{$all_clusters}) {

    my $cluster_genes = $single_cluster->get_Genes();
    my $max_cluster_exons = 0;
    my $max_cluster_seq = 0;
    my $transcript_keys = {};
    my $coding_transcript_keys = {};

    say "Analysing cluster of size: ".scalar(@{$cluster_genes});

    foreach my $single_gene (@{$cluster_genes}) {
      my $transcripts = $single_gene->get_all_Transcripts();
      my $transcript = shift(@{$transcripts});

      my $seq_length = $transcript->length();
      my $exon_count = scalar(@{$transcript->get_all_Exons()});
      if($exon_count > $max_cluster_exons) {
        $max_cluster_exons = $exon_count;
      }

      if($seq_length > $max_cluster_seq) {
        $max_cluster_seq = $seq_length;
      }
    }

    # Get 80 percent of the max exon count. If the remainder is > 0.5 round the count up
    my $raw_max_cluster_exon_threshold  = (($max_cluster_exons * 4) / 5);
    my $remainder = ($raw_max_cluster_exon_threshold - int($raw_max_cluster_exon_threshold));
    if($remainder > 0.5) {
      $raw_max_cluster_exon_threshold++;
    }

    my $max_cluster_exon_threshold = int($raw_max_cluster_exon_threshold);
    foreach my $single_gene (@{$cluster_genes}) {
      my $transcripts = $single_gene->get_all_Transcripts();
      my $transcript = shift(@{$transcripts});

      my $key = generate_transcript_key($transcript);
      if($transcript_keys->{$key}) {
        say "Skipping redundant transcript: ".$transcript->dbID();
        next;
      } else {
        $transcript_keys->{$key} = 1;
        # Not really sure what I was planning with this code for the coding_key. Maybe I'll remember some day
        # my $coding_key = generate_transcript_key($transcript,1);
        # if($coding_transcript_keys->{$coding_key}) {
        #
        # } else {
        #   $coding_transcript_keys->{$coding_key} = [$transcript->length,scalar(@{$transcript->get_all_Exons})];
        # }
      }

      my $seq_length = $transcript->length();
      my $exon_count = scalar(@{$transcript->get_all_Exons()});
      # This is only an initial implementation. But bascially the idea is to only take single/two exons models
      # If there aren't overlapping models with more exons, if they're all that is present then only
      # keep them if there's some seq (will filter out a lot of RNA-seq fragement). Then keep everything with
      # 5 exons or more, in the case of less than > 2 exons and < 5 exons it will take anything that has the
      # the same amount of exons as the transcript with the most exons in the cluster or anything with more
      # sequence than the transcript with the most amount of exons. These are all just based on things I've
      # observed in the stats for various gene sets, but should do a reasonable job of cleaning a set


      if($max_cluster_exons < 5) {
        my $potentially_keep = assess_transcript_stats($transcript,$source_type);
        unless($potentially_keep) {
          next;
        }
      }

      say "EC/MCT: ".$exon_count.":".$max_cluster_exon_threshold;
      if($exon_count >= $max_cluster_exon_threshold) {
        say "1: ".$transcript->dbID;
        push(@{$output_genes},$single_gene);
      } elsif($seq_length >= $max_cluster_seq) {
        say "2: ".$transcript->dbID;
        push(@{$output_genes},$single_gene);
      } else {
        say "3: ".$transcript->dbID;
        say "Skipping transcript: ".$transcript->dbID();
      }
    }
    say "Output gene array size: ".scalar(@{$output_genes});
  }
  return($output_genes);
}


sub assess_transcript_stats {
  my ($transcript,$source_type) = @_;

  my $min_translatable_seq_length = 150;
  if($source_type eq 'rnaseq') {
    $min_translatable_seq_length = 500;
  }

  my $translateable_seq = $transcript->translateable_seq();
  if($translateable_seq) {
    my $start_codon  = uc( substr( $translateable_seq, 0, 3 ) );
    my $end_codon  = uc( substr( $translateable_seq, -3 ) );
#    unless($start_codon eq "ATG" && (($end_codon eq 'TAG') || ($end_codon eq 'TAA') || ($end_codon eq 'TGA')) && length($translateable_seq) >= 300) {
    unless($start_codon eq "ATG" && length($translateable_seq) >= $min_translatable_seq_length) {
      return(0);
    } else {
      return(1);
    }
  } else {
    return(0);
  }
}


sub generate_transcript_key {
  my ($transcript,$cds_only) = @_;


  my $key = $transcript->seq_region_name.":".
            $transcript->strand.":".
            $transcript->translateable_seq.":";

  my $exons;

  if($cds_only) {
    $exons = $transcript->get_all_translateable_Exons;
  } else {
    $exons = $transcript->get_all_Exons;
  }

  foreach my $exon (@{$exons}) {
    $key .= $exon->seq_region_start.":".$exon->seq_region_end.":";
  }

  return $key;
}


sub avg_intron_size {
  my ($transcript) = @_;

  my $short_length = 75;
  my $max_short = 2;
  my $introns = $transcript->get_all_Introns();
  my $total_introns = scalar(@$introns);
  my $total_length = 0;
  my $total_short = 0;
  foreach my $intron (@$introns) {
    my $length = $intron->length();
    say "FERGAL INTRON LEN: ".$length;
    $total_length += $length;
    if($length < $short_length) {
      $total_short++;
    }
  }

  if($total_introns) {
    my $avg_length = $total_length / $total_introns;
    if($avg_length < $short_length) {
      say "Transcript avg intron length too short: ".$transcript->dbID." ".$avg_length." ".$total_short;
      return 1;
    } elsif($total_short > $max_short) {
      say "Too many short introns: ".$transcript->dbID." ".$avg_length." ".$total_short;
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}


sub select_most_supported_transcripts {
  my ($clustered,$unclustered) = @_;

  # This is expecting single transcript gene clusters
  my $all_clusters = [@{$clustered},@{$unclustered}];

  my $most_supported_output_genes = [];
  foreach my $single_cluster (@{$all_clusters}) {

    my $cluster_genes = $single_cluster->get_Genes();
    my $max_cluster_exons = 0;
    my $max_cluster_seq = 0;
    my $transcript_keys = {};
    my $coding_transcript_keys = {};
    my $observed_exon_counts = {};
    say "Analysing cluster of size: ".scalar(@{$cluster_genes});

    # Generating the observer exon counts
    foreach my $single_gene (@{$cluster_genes}) {
      my $exons = $single_gene->get_all_Exons;
      foreach my $exon (@$exons) {
        my $exon_key = $exon->seq_region_start.":".$exon->seq_region_end.":".$exon->strand;
        if($observed_exon_counts->{$exon_key}) {
          $observed_exon_counts->{$exon_key}++;
        } else {
          $observed_exon_counts->{$exon_key} = 1;
	}
      }
    } # end foreach my $single_gene

    # Generating calculating observed exon count scores per gene
    my $top_gene;
    my $current_top_score = 0;
    foreach my $single_gene (@{$cluster_genes}) {
      my $total_observed = 0;
      my $exons = $single_gene->get_all_Exons;
      foreach my $exon (@$exons) {
        my $exon_key = $exon->seq_region_start.":".$exon->seq_region_end.":".$exon->strand;
        my $observed_count = $observed_exon_counts->{$exon_key};
        $total_observed += $observed_count;
      }
      my $avg = $total_observed / scalar(@$exons);
      if($avg > $current_top_score) {
        $top_gene = $single_gene;
        $current_top_score = $avg;
      } elsif($avg == $current_top_score) {
        if(scalar(@$exons) > scalar(@{$top_gene->get_all_Exons})) {
          $top_gene = $single_gene;
        }
      }
    } # foreach my $single_gene
    push(@$most_supported_output_genes,$top_gene);
  }

  return($most_supported_output_genes);
}
