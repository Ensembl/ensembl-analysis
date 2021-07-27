=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=pod

=head1 NAME

report_projection_stats.pl

=head1 SYNOPSIS

perl report_projection_stats.pl --reference_user READONLYUSER --reference_host REFHOST --reference_port REFPORT --reference_dbname REFDBNAME --query_user READONLYUSER --query_host QHOST --query_port QPORT --query_dbname QDBNAME --cluster_by_transcript_sid 1 --reference_cs_version GRCh37 --query_cs_version GRCh38 --output_path /my/output/path

=head1 DESCRIPTION

This script compares two gene sets from a reference and a query databases based on clustering the genes by cluster_Genes algorithm or
by making the link between the reference and the query genes using the reference transcript id stored in the query transcripts (--cluster_by_transcript_sid 1). The latter can be run if HiveProjectionMinimap has been used to project the genes from the reference to the query database.

=head1 CONTACT

Post general queries to B<http://lists.ensembl.org/mailman/listinfo/dev>

=head1 APPENDIX
=cut

use warnings;
use strict;
use feature 'say';
use Digest::MD5 qw(md5);
use Data::Dumper;

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Getopt::Long qw(:config no_ignore_case);
use List::MoreUtils qw(uniq);
use File::Spec::Functions qw(catfile);

my $reference_cs_version;
my $query_cs_version;

my $reference_dbname;
my $reference_user;
my $reference_host;
my $reference_port;

my $query_dbname;
my $query_user;
my $query_host;
my $query_port;

my $toplevel_region_name;
my $reference_canonical_only = 0;
my $query_canonical_only = 0;
my $cluster_by_transcript_sid = 1;

my $gene_biotypes;
my $output_path;

my $options = GetOptions ("reference_user=s"          => \$reference_user,
                          "reference_host=s"          => \$reference_host,
                          "reference_port=i"          => \$reference_port,
                          "reference_dbname=s"        => \$reference_dbname,
                          "query_user=s"              => \$query_user,
                          "query_host=s"              => \$query_host,
                          "query_port=i"              => \$query_port,
                          "query_dbname=s"            => \$query_dbname,
                          "toplevel_region_name=s"    => \$toplevel_region_name,
                          "reference_canonical_only!" => \$reference_canonical_only,
                          "query_canonical_only!"     => \$query_canonical_only,
                          "gene_biotypes=s"           => \$gene_biotypes,
                          "reference_cs_version=s"    => \$reference_cs_version,
                          "query_cs_version=s"        => \$query_cs_version,
                          "cluster_by_transcript_sid" => \$cluster_by_transcript_sid,
                          "output_path=s"             => \$output_path);

my $reference_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $reference_port,
  -user    => $reference_user,
  -host    => $reference_host,
  -dbname  => $reference_dbname);
$reference_db->dnadb($reference_db);

my $query_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $query_port,
  -user    => $query_user,
  -host    => $query_host,
  -dbname  => $query_dbname);
$query_db->dnadb($query_db);

my $reference_slice_adaptor = $reference_db->get_SliceAdaptor();
my $query_slice_adaptor = $query_db->get_SliceAdaptor();
my $slices;

my $num_not_found_reference_regions = 0;
my $num_not_found_reference_genes = 0;
my $num_not_found_reference_transcripts = 0;
my $num_split_reference_genes = 0;
my $num_split_query_genes = 0;

my $file_prefix = time().'_'.$reference_cs_version."_".$query_cs_version;

open(my $not_found_reference_regions_fh,'>',catfile($output_path,$file_prefix.'_not_found_reference_regions.tsv')) || die "Couldn't open file *not_found_reference_regions.tsv, $!";
open(my $not_found_reference_genes_fh,'>',catfile($output_path,$file_prefix.'_not_found_reference_genes.tsv')) || die "Couldn't open file *not_found_reference_genes.tsv, $!";
open(my $not_found_reference_transcripts_fh,'>',catfile($output_path,$file_prefix.'_not_found_reference_transcripts.tsv')) || die "Couldn't open file *not_found_reference_transcripts.tsv, $!";
open(my $split_reference_genes_fh,'>',catfile($output_path,$file_prefix.'_split_reference_genes.tsv')) || die "Couldn't open file *split_reference_genes.tsv, $!";
open(my $split_query_genes_fh,'>',catfile($output_path,$file_prefix.'_split_query_genes.tsv')) || die "Couldn't open file *split_query_genes.tsv, $!";

my $reference_gene_adaptor = $reference_db->get_GeneAdaptor();
my $reference_transcript_adaptor = $reference_db->get_TranscriptAdaptor();
my $query_gene_adaptor = $query_db->get_GeneAdaptor();
my $query_transcript_adaptor = $query_db->get_TranscriptAdaptor();

if ($toplevel_region_name) {
  
  my $reference_slice = $reference_slice_adaptor->fetch_by_region('toplevel',
                                                                  $toplevel_region_name,
                                                                  undef,
                                                                  undef,
                                                                  undef,
                                                                  $reference_cs_version);
  if (!$reference_slice) {
    print STDERR "Could not fetch reference toplevel region name slice ".$toplevel_region_name."\n";
    exit;
  }

  my $query_slice = $query_slice_adaptor->fetch_by_region('toplevel',
                                                           $toplevel_region_name,
                                                           undef,
                                                           undef,
                                                           undef,
                                                           $query_cs_version);

  if (!$query_slice) {
    print STDERR "Could not fetch query toplevel region name slice ".$toplevel_region_name."\n";
    exit;
  }

  push(@$slices,[$reference_slice,$query_slice]);

} else {
  my $reference_slices = $reference_slice_adaptor->fetch_all('toplevel',$reference_cs_version);
  foreach my $reference_slice (@$reference_slices) {
    #my $query_slice = $query_slice_adaptor->fetch_by_name($reference_slice->name());
    
    # it is required that the seq region names of the reference and query assemblies match and their regions are equivalent
    # if there are multiple coord systems in the query database, the latest one will be used (fetched as 'toplevel')

    my $reference_slice_sr_name = $reference_slice->seq_region_name();

    # fetch_by_region(coord_system_name,seq_region_name,start,end,strand,version,no_fuzz)
    my $query_slice = $query_slice_adaptor->fetch_by_region('toplevel',
                                                            $reference_slice_sr_name,
                                                            undef,
                                                            undef,
                                                            undef,
                                                            $query_cs_version);

    if ($query_slice) {
      push(@$slices,[$reference_slice,$query_slice]);
    } else {
      print STDERR "Reference region ".$reference_slice_sr_name." not found in query database.\n";
      $num_not_found_reference_regions++;
      say $not_found_reference_regions_fh join("\t",
                                               $reference_cs_version,
                                               $reference_slice_sr_name,
                                               $query_cs_version,
                                               $reference_slice_sr_name);

      # We want to include the missing genes and transcripts in the stats.
      # Since they will not be fetched later on, we loop through them now.
      my $missing_reference_genes = ();
      if ($gene_biotypes) {
        foreach my $biotype (split(',',$gene_biotypes)) {
          push(@$missing_reference_genes,@{$reference_slice->get_all_Genes_by_type($biotype)});
        }
      } else {
        $missing_reference_genes = $reference_slice->get_all_Genes();
      }

      foreach my $missing_reference_gene (@$missing_reference_genes) {
        foreach my $missing_reference_transcript (@{$missing_reference_gene->get_all_Transcripts()}) {
          print STDERR "Reference transcript ".$missing_reference_transcript->dbID()." not found in the query database.\n";
          $num_not_found_reference_transcripts++;
          say $not_found_reference_transcripts_fh ensembl_obj2str($missing_reference_transcript);
        }
        print STDERR "No transcript id from reference gene ".$missing_reference_gene->dbID()." found in any query gene stable id.\n";
        $num_not_found_reference_genes++;
        say $not_found_reference_genes_fh ensembl_obj2str($missing_reference_gene);
      }
    }
  }
}

my $all_comparison_results = {};
foreach my $slice_pair (@$slices) {

  # md5 hash is used to store features where a exact match is considered
  # The string hash is for when partial matches are needed (e.g an intron string being a substring of another intron string)
  my $feature_md5s_hash = {};
  my $feature_string_hash = {};

  my $reference_slice = ${$slice_pair}[0];
  my $query_slice = ${$slice_pair}[1];

  say "Reference slice: ".$reference_slice->name();
  say "Query slice: ".$query_slice->name();

  my $reference_genes = ();
  my $query_genes = ();
  if ($gene_biotypes) {
    foreach my $biotype (split(',',$gene_biotypes)) {
      push(@$reference_genes,@{$reference_slice->get_all_Genes_by_type($biotype)});
      push(@$query_genes,@{$query_slice->get_all_Genes_by_type($biotype)});
    }
  } else {
    $reference_genes = $reference_slice->get_all_Genes();
    $query_genes = $query_slice->get_all_Genes();
  }

  say "Reference gene count: ".scalar(@$reference_genes);
  say "Query gene count: ".scalar(@$query_genes);

  my $comparison_results;

  if ($cluster_by_transcript_sid) {

    foreach my $reference_gene (@$reference_genes) {

      my @query_gene_ids = ();
      my @query_genes = ();
      foreach my $reference_transcript (@{$reference_gene->get_all_Transcripts()}) {
        my $reference_transcript_id = $reference_transcript->dbID();
        # note that the query (or projected or target) transcripts stable IDs store the reference transcript IDs
        my $query_transcript = $query_transcript_adaptor->fetch_by_stable_id($reference_transcript_id);
        if ($query_transcript) {
          my $query_gene = $query_transcript->get_Gene();
          push(@query_gene_ids,$query_gene->dbID());
          push(@query_genes,$query_gene);
        } else {
          print STDERR "Reference transcript ".$reference_transcript_id." not found in the query database.\n";
          $num_not_found_reference_transcripts++;
          say $not_found_reference_transcripts_fh ensembl_obj2str($reference_transcript);
        }
      }

      my @uniq_query_gene_ids = uniq(@query_gene_ids);
      my $num_uniq_query_gene_ids = scalar(@uniq_query_gene_ids);
      if ($num_uniq_query_gene_ids == 0) {
        my $reference_gene_id = $reference_gene->dbID();
        print STDERR "No transcript id from reference gene ".$reference_gene_id." found in any query gene stable id.\n";
        $num_not_found_reference_genes++;
        say $not_found_reference_genes_fh ensembl_obj2str($reference_gene);

      } elsif ($num_uniq_query_gene_ids > 1) {
        print STDERR "Transcripts in reference gene ".$reference_gene->dbID()." found in multiple query genes ".join(',',@uniq_query_gene_ids)."\n";
        $num_split_reference_genes++;
        $num_split_query_genes += $num_uniq_query_gene_ids;
        say $split_reference_genes_fh ensembl_obj2str($reference_gene);
        foreach my $uniq_query_gene_id (@uniq_query_gene_ids) {
          my $split_query_gene = $query_gene_adaptor->fetch_by_dbID($uniq_query_gene_id);
          say $split_query_genes_fh ensembl_obj2str($split_query_gene);
        }

        my @reference_genes = ($reference_gene);
        my @uniq_query_genes = ();
        foreach my $uniq_query_gene_id (@uniq_query_gene_ids) {
          push(@uniq_query_genes,$query_gene_adaptor->fetch_by_dbID($uniq_query_gene_id));
        }

        my $comparison_results = compare_genes(\@reference_genes,\@uniq_query_genes);
        @$all_comparison_results{ keys %$comparison_results } = values %$comparison_results;

      } else { # one2one
        my @reference_genes = ($reference_gene);
        my @uniq_query_genes = ($query_genes[0]); # all elements in this array are the same gene, choose any
        my $comparison_results = compare_genes(\@reference_genes,\@uniq_query_genes);
        @$all_comparison_results{ keys %$comparison_results } = values %$comparison_results;
      }

    } # end foreach my reference_gene

  } else { # if cluster_by_transcript_id
    $comparison_results = process_genes($reference_genes,$query_genes);
  }
  @$all_comparison_results{ keys %$comparison_results } = values %$comparison_results;
}

if (%$all_comparison_results) {
  print_results($all_comparison_results);
  print $num_not_found_reference_regions." reference regions not found in the query database.\n";
  print $num_not_found_reference_genes." reference genes not found (no single transcript of these genes found in the query database).\n";  
  print $num_not_found_reference_transcripts." reference transcripts not found in the query database.\n";
  print $num_split_reference_genes." reference genes split into ".$num_split_query_genes." query genes.\n";

  close($not_found_reference_regions_fh);
  close($not_found_reference_genes_fh);
  close($not_found_reference_transcripts_fh);
  close($split_reference_genes_fh);
  close($split_query_genes_fh);

} else {
  print "No comparison results.\n";
}

=head2 process_genes
 Arg [1]    : Arrayref of Bio::EnsEMBL::Gene
 Arg [2]    : Arrayref of Bio::EnsEMBL::Gene
 Description: This returns what type of matches there are for each reference gene
              based on exon and intron strings using the clustering algorithm first.
 Returntype : Hashref with the dbID of the gene as the key.
 Exceptions : None
=cut

sub process_genes {
  my ($reference_genes,$query_genes) = @_;

  my $all_comparison_results = {};
  my $mixed_loci = 0;
  my $reference_loci = 0;
  my $query_loci = 0;

  my $reference_biotypes_hash = get_all_biotypes([@$reference_genes]);
  my $reference_biotypes_array = [keys(%$reference_biotypes_hash)];

  # To separate out the genes more easily modify the biotypes
  # for the query genes to be unique
  foreach my $query_gene (@$query_genes) {
    $query_gene->biotype('query_',$query_gene->biotype());
  }
  my $query_biotypes_hash = get_all_biotypes([@$query_genes]);
  my $query_biotypes_array = [keys(%$query_biotypes_hash)];

  my $types_hash;
  my $reference_set_name = 'reference_genes';
  my $query_set_name = 'query_genes';
  $types_hash->{$reference_set_name} = $reference_biotypes_array;
  $types_hash->{$query_set_name} = $query_biotypes_array;

  my $all_genes = [@$reference_genes,@$query_genes];
  say "Clustering genes from input_dbs...";
  my ($clusters, $unclustered) = cluster_Genes($all_genes,$types_hash);
  foreach my $cluster (@$clusters) {
    my $reference_cluster_genes = $cluster->get_Genes_by_Set($reference_set_name);
    my $query_cluster_genes = $cluster->get_Genes_by_Set($query_set_name);

    if(scalar(@$reference_cluster_genes) and scalar(@$query_cluster_genes)) {
      $mixed_loci++;
      my $comparison_results = compare_genes($reference_cluster_genes,$query_cluster_genes);
      @$all_comparison_results{ keys %$comparison_results } = values %$comparison_results;
    } elsif(scalar(@$reference_cluster_genes)) {
      $reference_loci++;
    } else {
      $query_loci++;
    }
  } # End foreach my $cluster

  foreach my $singleton (@$unclustered) {
    # Note will loop even if this should just be single genes
    my $reference_cluster_genes = $singleton->get_Genes_by_Set($reference_set_name);
    my $query_cluster_genes = $singleton->get_Genes_by_Set($query_set_name);

    if(scalar(@$reference_cluster_genes) and scalar(@$query_cluster_genes)) {
      $mixed_loci++;
      my $comparison_results = compare_genes($reference_cluster_genes,$query_cluster_genes);
      @$all_comparison_results{ keys %$comparison_results } = values %$comparison_results;
    } elsif(scalar(@$reference_cluster_genes)) {
      $reference_loci++;
    } else {
      $query_loci++;
    }
  }

  say "Mixed loci: ".$mixed_loci;
  say "Reference loci: ".$reference_loci;
  say "Query loci: ".$query_loci;

  return $all_comparison_results;
}

=head2 compare_genes
 Arg [1]    : Arrayref of Bio::EnsEMBL::Gene
 Arg [2]    : Arrayref of Bio::EnsEMBL::Gene
 Description: This returns what type of matches there are for each reference gene
              based on exon and intron strings.
 Returntype : Hashref with the dbID of the gene as the key.
 Exceptions : None
=cut

sub compare_genes {
  my ($reference_genes,$query_genes) = @_;

  my $query_gene_exons = {};
  my $query_gene_introns = {};
  my $query_gene_cds_exons = {};
  my $query_gene_cds_introns = {};
  my $query_gene_md5s = {};
  my $query_gene_cds_md5s = {};

  # This is what will be passed back, comparison numbers for each reference gene passed in
  my $comparison_results = {};

  # This loops through the query genes and records features in terms of introns, exons, cds
  # and the actual structures of the transcripts.
  foreach my $query_gene (@$query_genes) {
    my $query_gene_id = $query_gene->dbID();
    my $query_transcripts;
    if ($query_canonical_only) {
      $query_transcripts = [$query_gene->canonical_transcript()];
    } else {
      $query_transcripts = $query_gene->get_all_Transcripts();
    }

    foreach my $query_transcript (@$query_transcripts) {
      my $exons = $query_transcript->get_all_Exons();
      my $introns = $query_transcript->get_all_Introns();
      my $exon_string = generate_exon_string($exons);
      my $intron_string = generate_intron_string($introns);
      my $transcript_string = $exon_string;

      foreach my $exon (@$exons) {
        my $single_exon_string = generate_exon_string([$exon]);
        $query_gene_exons->{$single_exon_string} = 1;
      }

      foreach my $intron (@$introns) {
        my $single_intron_string = generate_intron_string([$intron]);
        $query_gene_introns->{$single_intron_string} = 1;
      }

      if ($query_transcript->translation()) {
        my $cds_exons = $query_transcript->get_all_CDS();
        my $cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
        my $cds_introns = $cds_transcript->get_all_Introns();
        my $cds_exon_string = generate_exon_string($cds_exons);
        my $cds_intron_string = generate_intron_string($cds_introns);
        $transcript_string .= ":".$cds_exon_string;
        $query_gene_cds_md5s->{$cds_exon_string} = 1;

        foreach my $cds_exon (@$cds_exons) {
          my $single_exon_string = generate_exon_string([$cds_exon]);
          $query_gene_cds_exons->{$single_exon_string} = 1;
        }

        foreach my $cds_intron (@$cds_introns) {
          my $single_intron_string = generate_intron_string([$cds_intron]);
          $query_gene_cds_introns->{$single_intron_string} = 1;
        }
      }

      $query_gene_md5s->{$transcript_string} = 1;
    } # End foreach my $query_transcript
  } # End foreach my $query_gene

  foreach my $reference_gene (@$reference_genes) {
    my $reference_biotype = $reference_gene->biotype();
    my $total_exons = 0;
    my $total_introns = 0;
    my $total_cds_exons = 0;
    my $total_cds_introns = 0;
    my $total_transcripts = 0;
    my $total_cds_transcripts = 0;
    my $matched_exons = 0;
    my $matched_introns = 0;
    my $matched_cds_exons = 0;
    my $matched_cds_introns = 0;
    my $matched_transcripts = 0;
    my $matched_cds_transcripts = 0;
    my $matched_gene = 0;
    my $matched_cds_gene = 0;
    my $is_cds_gene = 0;

    my $reference_gene_id = $reference_gene->dbID();
    my $reference_transcripts;
    if($reference_canonical_only eq "yes") {
      $reference_transcripts = [$reference_gene->canonical_transcript()];
    } else {
      $reference_transcripts = $reference_gene->get_all_Transcripts();
    }
    foreach my $reference_transcript (@$reference_transcripts) {
      $total_transcripts++;
      my $exons = $reference_transcript->get_all_Exons();
      my $introns = $reference_transcript->get_all_Introns();
      my $exon_string = generate_exon_string($exons);
      my $intron_string = generate_intron_string($introns);
      my $transcript_string = $exon_string;

      foreach my $exon (@$exons) {
        $total_exons++;
        my $single_exon_string = generate_exon_string([$exon]);
        if($query_gene_exons->{$single_exon_string}) {
          $matched_exons++;
        }
      }

      foreach my $intron (@$introns) {
        $total_introns++;
        my $single_intron_string = generate_intron_string([$intron]);
        if($query_gene_introns->{$single_intron_string}) {
          $matched_introns++;
        }
      }

      if ($reference_transcript->translation()) {
        $total_cds_transcripts++;
        my $cds_exons = $reference_transcript->get_all_CDS();
        my $cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
        my $cds_introns = $cds_transcript->get_all_Introns();
        my $cds_exon_string = generate_exon_string($cds_exons);
        my $cds_intron_string = generate_intron_string($cds_introns);
        $transcript_string .= ":".$cds_exon_string;#$reference_transcript->translate->seq();

        if($query_gene_cds_md5s->{$cds_exon_string}) {
          $matched_cds_transcripts++;
        }

        foreach my $cds_exon (@$cds_exons) {
          $total_cds_exons++;
          my $single_exon_string = generate_exon_string([$cds_exon]);
          if($query_gene_cds_exons->{$single_exon_string}) {
            $matched_cds_exons++;
          }
        }

        foreach my $cds_intron (@$cds_introns) {
          $total_cds_introns++;
          my $single_intron_string = generate_intron_string([$cds_intron]);
          if($query_gene_cds_introns->{$single_intron_string}) {
            $matched_cds_introns++;
          }
        }
      }

      if($query_gene_md5s->{$transcript_string}) {
        $matched_transcripts++;
      }

    } # End foreach my $reference_transcript

    if ($total_cds_transcripts) {
      $is_cds_gene = 1;
    }

    if ($matched_transcripts eq $total_transcripts) {
      $matched_gene++;
    }

    if ($matched_cds_transcripts and ($matched_cds_transcripts == $total_cds_transcripts)) {
      $matched_cds_gene++;
    }

    $comparison_results->{$reference_gene_id} = {'biotype' => $reference_biotype,
                                                 'exons' => [$matched_exons,$total_exons],
                                                 'introns' => [$matched_introns,$total_introns],
                                                 'cds_exons' => [$matched_cds_exons,$total_cds_exons],
                                                 'cds_introns' => [$matched_cds_introns,$total_cds_introns],
                                                 'transcript_structure' => [$matched_transcripts,$total_transcripts],
                                                 'transcript_cds_structure' => [$matched_cds_transcripts,$total_cds_transcripts],
                                                 'gene_structure' => [$matched_gene,1],
                                                 'gene_cds_structure' => [$matched_cds_gene,$is_cds_gene]};
  } # End foreach my $reference_gene

  return $comparison_results;

  # By this point we have processed all the reference genes so that means we have the required info in various hashes
  # Each hash has the reference gene id as the initial key and this points to whatever else needs to be tracked
  # That way if there are multiple reference genes in a cluster it will be easy to separate them back out
  # The info being tracked:
  # All the individual exons in the reference genes
  # All the individual cds exons in the reference genes
  # All the individual introns in the reference genes
  # All the individual cds introns in the reference genes
  # An md5 of each transcript structure in the gene
  # An md4 of each CDS structure in the gene
  # Note that the above use of all can mean all, or just all for the canonical transcript, depending on what has been set
  # Everything that is recorded in the hashes is set to 0. When we find a match in the query set we'll just increment
  # the count. That way any feature with 1 or more as a value will have been matched exactly in the query set

}

=head2 print_results
 Arg [1]    : Hashref containing the stats to print.
 Description: It prints to STDOUT the stats found in Arg [1].
 Returntype : None
 Exceptions : None
=cut

sub print_results {
  my ($results) = @_;


  my $all_gene_structures_match = 0;
  my $all_gene_structures = 0;
  my $all_gene_cds_structures_match = 0;
  my $all_gene_cds_structures = 0;
  my $all_transcript_structures_match = 0;
  my $all_transcript_structures = 0;
  my $all_transcript_cds_structures_match = 0;
  my $all_transcript_cds_structures = 0;
  my $all_exons_match = 0;
  my $all_exons = 0;
  my $all_introns_match = 0;
  my $all_introns = 0;
  my $all_cds_exons_match = 0;
  my $all_cds_exons = 0;
  my $all_cds_introns_match = 0;
  my $all_cds_introns = 0;


  foreach my $gene_id (keys(%$results)) {
    my $gene_results = $results->{$gene_id};
    my $gs_match = $gene_results->{'gene_structure'};
    $all_gene_structures_match += ${$gs_match}[0];
    $all_gene_structures += ${$gs_match}[1];
    my $gs_cds_match = $gene_results->{'gene_cds_structure'};
    $all_gene_cds_structures_match += ${$gs_cds_match}[0];
    $all_gene_cds_structures += ${$gs_cds_match}[1];
    my $ts_match = $gene_results->{'transcript_structure'};
    $all_transcript_structures_match += ${$ts_match}[0];
    $all_transcript_structures += ${$ts_match}[1];
    my $ts_cds_match = $gene_results->{'transcript_cds_structure'};
    $all_transcript_cds_structures_match += ${$ts_cds_match}[0];
    $all_transcript_cds_structures += ${$ts_cds_match}[1];
    my $exons_match = $gene_results->{'exons'};
    $all_exons_match += ${$exons_match}[0];
    $all_exons += ${$exons_match}[1];
    my $introns_match = $gene_results->{'introns'};
    $all_introns_match += ${$introns_match}[0];
    $all_introns += ${$introns_match}[1];
    my $cds_exons_match = $gene_results->{'cds_exons'};
    $all_cds_exons_match += ${$cds_exons_match}[0];
    $all_cds_exons += ${$cds_exons_match}[1];
    my $cds_introns_match = $gene_results->{'cds_introns'};
    $all_cds_introns_match += ${$cds_introns_match}[0];
    $all_cds_introns += ${$cds_introns_match}[1];
  }


  my $match_gs_perc = sprintf("%.2f", ($all_gene_structures_match/$all_gene_structures) * 100);
  my $match_gs_cds_perc = sprintf("%.2f", ($all_gene_cds_structures_match/$all_gene_cds_structures) * 100);
  my $match_ts_perc = sprintf("%.2f", ($all_transcript_structures_match/$all_transcript_structures) * 100);
  my $match_ts_cds_perc = sprintf("%.2f", ($all_transcript_cds_structures_match/$all_transcript_cds_structures) * 100);
  my $match_exons_perc = sprintf("%.2f", ($all_exons_match/$all_exons) * 100);
  my $match_introns_perc = sprintf("%.2f", ($all_introns_match/$all_introns) * 100);
  my $match_cds_exons_perc = sprintf("%.2f", ($all_cds_exons_match/$all_cds_exons) * 100);
  my $match_cds_introns_perc = sprintf("%.2f", ($all_cds_introns_match/$all_cds_introns) * 100);

  say "Matched/Total gene structures: ".$all_gene_structures_match."/".$all_gene_structures.", ".$match_gs_perc."%";
  say "Matched/Total CDS gene structures: ".$all_gene_cds_structures_match."/".$all_gene_cds_structures.", ".$match_gs_cds_perc."%";

  say "Matched/Total transcript structures: ".$all_transcript_structures_match."/".$all_transcript_structures.", ".$match_ts_perc."%";
  say "Matched/Total CDS transcript structures: ".$all_transcript_cds_structures_match."/".$all_transcript_cds_structures.", ".$match_ts_cds_perc."%";
  say "Matched/Total exon structures: ".$all_exons_match."/".$all_exons.", ".$match_exons_perc."%";
  say "Matched/Total intron structures: ".$all_introns_match."/".$all_introns.", ".$match_introns_perc."%";
  say "Matched/Total CDS exon structures: ".$all_cds_exons_match."/".$all_cds_exons.", ".$match_cds_exons_perc."%";
  say "Matched/Total CDS intron structures: ".$all_cds_introns_match."/".$all_cds_introns.", ".$match_cds_introns_perc."%";
}

=head2 generate_intron_string
 Arg [1]    : Arrayref of Bio::EnsEMBL::Intron
 Description: It concatenates a few attributes from the Ensembl object in Arg[1] into a string.
 Returntype : String
 Exceptions : None
=cut

sub generate_intron_string {
  my ($intron_array) = @_;

  my $intron_string = "";
  foreach my $intron (@{$intron_array}) {
    $intron_string .= $intron->length().":".$intron->strand().":";
  }
  return $intron_string;
}

=head2 get_all_biotypes
 Arg [1]    : Arrayref of Bio::EnsEMBL::Gene
 Description: It stores all biotypes from the genes in Arg[1] in a hash.
 Returntype : Hashref
 Exceptions : None
=cut

sub get_all_biotypes {
  my ($master_genes_array) = @_;

  my $master_biotypes_hash = {};

  foreach my $gene (@{$master_genes_array}) {
    unless($master_biotypes_hash->{$gene->biotype}) {
      $master_biotypes_hash->{$gene->biotype} = 1;
    }
  }
  return $master_biotypes_hash;
}

=head2 generate_exon_string
 Arg [1]    : Arrayref of Bio::EnsEMBL::Exon
 Description: It concatenates a few attributes from the Ensembl object in Arg[1] into a string.
 Returntype : String
 Exceptions : None
=cut

sub generate_exon_string {
  my ($exons) = @_;

  my $exon_string = "";
  foreach my $exon (@$exons) {
    $exon_string .= $exon->length().":".$exon->strand.":";
  }
  return $exon_string;
}

=head2 ensembl_obj2str
 Arg [1]    : Bio::EnsEMBL::Gene or Bio::EnsEMBL::Transcript
 Description: It concatenates a few attributes from the Ensembl object in Arg[1] into a string.
 Returntype : String
 Exceptions : None
=cut

sub ensembl_obj2str {
  my $ensembl_obj = shift();

  my $display_id = "";
  my $description = "";
  my $display_xref = $ensembl_obj->display_xref();
  if ($display_xref) {
    $display_id = $display_xref->display_id();
    $description = $display_xref->description();
  }

  return join("\t",
              $ensembl_obj->dbID(),
              $ensembl_obj->stable_id(),
              $ensembl_obj->biotype(),
              $ensembl_obj->seq_region_name(),
              $ensembl_obj->seq_region_strand(),
              $ensembl_obj->seq_region_start(),
              $ensembl_obj->seq_region_end(),
              $display_id,
              $description);
}

1;
