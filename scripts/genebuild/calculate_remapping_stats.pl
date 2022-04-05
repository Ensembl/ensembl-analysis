=head1 LICENSE

 Copyright [2022] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 DESCRIPTION

 A script to compare genes and transcripts in a reference database to those in a target database and generate mapping stats. Works under the assumption that the
 mapping process used the Minimap2Remap module and also that pipeline in general as it bases the stats off values in the gene and transcript descriptions, which
 are produced by that module and the mapping pipeline in general. Assuming the data are compatible it will provide stats in terms of mapping percentages for genes
 and transcripts across various biotypes. It also utilises the information about transcript coverage and percent id to provide stats on those

=cut

use warnings;
use strict;
use feature 'say';
use Digest::MD5 qw(md5);

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $coord_system = 'toplevel';

my $reference_dbname = '';
my $reference_user   = '';
my $reference_host = '';
my $reference_port;

my $query_dbname = '';
my $query_user   = '';
my $query_host = '';
my $query_port;

my $ids_to_map_list;
my $output_dir = '';
my $output_file_prefix = "mapping_stats";

my $xy_scanner;

my $options = GetOptions ("reference_user=s"          => \$reference_user,
                          "reference_host=s"          => \$reference_host,
                          "reference_port=i"          => \$reference_port,
                          "reference_dbname=s"        => \$reference_dbname,
                          "query_user=s"              => \$query_user,
                          "query_host=s"              => \$query_host,
                          "query_port=i"              => \$query_port,
                          "query_dbname=s"            => \$query_dbname,
                          "output_dir=s"              => \$output_dir,
                          "output_file_prefix=s"      => \$output_file_prefix,
                          "xy_scanner=s"              => \$xy_scanner);


my $input_genes_file = $output_file_prefix.".input_gene_ids.txt";
my $missing_genes_file = $output_file_prefix.".missing_gene_ids.txt";
my $translation_seqs_file = $output_file_prefix.".pep.fa";
my $mapping_stats_file = $output_file_prefix.".mapping_stats.txt";
open(INPUT,">".$output_dir."/".$input_genes_file);
open(MISSING,">".$output_dir."/".$missing_genes_file);
open(TRANSLATION,">".$output_dir."/".$translation_seqs_file);
open(MAPPING,">".$output_dir."/".$mapping_stats_file);

my $reference_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $reference_port,
  -user    => $reference_user,
  -host    => $reference_host,
  -dbname  => $reference_dbname);

my $query_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $query_port,
  -user    => $query_user,
  -host    => $query_host,
  -dbname  => $query_dbname);


my $source_gene_ids_hash = {};
if($ids_to_map_list) {
  open(IN,$ids_to_map_list);
  while(<IN>) {
    my $line = $_;
    chomp($line);

    unless($line =~ /^ENS/) {
      next;
    }
    $source_gene_ids_hash->{$line} = 1;
  }
}


my $source_gene_adaptor = $reference_db->get_GeneAdaptor();
my $target_gene_adaptor = $query_db->get_GeneAdaptor();
my $source_slice_adaptor = $reference_db->get_SliceAdaptor();
my $target_slice_adaptor = $query_db->get_SliceAdaptor();

my $source_gene_info = {};
my $source_transcript_info = {};
my $source_slices = $source_slice_adaptor->fetch_all('toplevel');


# Decided to go by slice instead of the id list here, I think fetching genes individually by id would be too slow
# Maybe just sql to pull the info from the gene/transcript table would be better, but will go with the API for now
foreach my $slice (@$source_slices) {
  my $region_name = $slice->seq_region_name();

#  unless($region_name eq '1' || $region_name eq '10') {
#    next;
#  }

  say "Processing source slice: ".$region_name;
  my $genes = $slice->get_all_Genes();
  unless($ids_to_map_list) {
    say "Filtering ".scalar(@$genes)." initial genes";
    $genes = filter_genes($genes,$source_gene_ids_hash,$xy_scanner);
    say "Have ".scalar(@$genes)." genes post filtering";
  }
  foreach my $gene (@$genes) {
    # Just process genes in the original mapping list
    if($source_gene_ids_hash->{$gene->stable_id()}) {
      my $gene_biotype = $gene->biotype();
      my $gene_biotype_group = $gene->get_Biotype->biotype_group();
      my $gene_region_name = $gene->seq_region_name();
      $source_gene_info->{$gene->stable_id()}->{'gene_biotype'} = $gene_biotype;
      $source_gene_info->{$gene->stable_id()}->{'gene_biotype_group'} = $gene_biotype_group;
      $source_gene_info->{$gene->stable_id()}->{'gene_region_name'} = $gene_region_name;
      my $transcripts = $gene->get_all_Transcripts();
      foreach my $transcript (@$transcripts) {
        my $transcript_stable_id = $transcript->stable_id();
        my $transcript_biotype = $transcript->biotype();
        my $transcript_biotype_group = $transcript->get_Biotype->biotype_group();
        $source_transcript_info->{$transcript_stable_id}->{'transcript_biotype'} = $transcript_biotype;
        $source_transcript_info->{$transcript_stable_id}->{'transcript_biotype_group'} = $transcript_biotype_group;
      }
    } # End if(source_gene_ids_hash
  }
}

say "Finished processing source slices for genes and transcripts";

# This is slightly different for the target genes/transcripts as we want to parse the descriptions since they
# should have the most accurate record of the source stable ids and also cov/percent id for transcripts
my $target_gene_info = {};
my $target_transcript_info = {};
my $target_slices = $target_slice_adaptor->fetch_all('toplevel');
# NOTE forgot about the fact that there are potentially multiple mappings, so should put something in about that
foreach my $slice (@$target_slices) {
  my $region_name = $slice->seq_region_name();

  say "Processing target slice: ".$region_name;
  my $genes = $slice->get_all_Genes();
  foreach my $gene (@$genes) {
    #my $gene_description = $gene->description();
    #$gene_description =~ /;parent_gene=(.+);mapping_type=(.+)$/;
    #my $gene_versioned_stable_id = $1;
    #my $gene_type = $2;
    #my $gene_stable_id = $gene_versioned_stable_id;
    #$gene_stable_id =~ s/\.\d+//;
    my ($gene_stable_id) = @{$gene->get_all_Attributes('proj_parent_g')};

    # Just process genes in the original mapping list
    if($source_gene_ids_hash->{$gene_stable_id}) {
      my $gene_biotype = $gene->biotype();
      my $gene_biotype_group = $gene->get_Biotype->biotype_group();
      $target_gene_info->{$gene_stable_id}->{'gene_biotype'} = $gene_biotype;
      $target_gene_info->{$gene_stable_id}->{'gene_biotype_group'} = $gene_biotype_group;
      #$target_gene_info->{$gene_stable_id}->{'gene_type'} = $gene_type;

      my $transcripts = $gene->get_all_Transcripts();
      foreach my $transcript (@$transcripts) {
        if ($transcript->translation()) {
          say TRANSLATION ">".$transcript->stable_id();
          say TRANSLATION $transcript->translation->seq();
        }

       	my $transcript_biotype = $transcript->biotype();
        my $transcript_biotype_group = $transcript->get_Biotype->biotype_group();
        my $transcript_description = $transcript->description();

        # DUE TO BUG IN FINDPARALOGUES NOT ADDING DESCRIPTION
        unless($transcript_description) {
          next;
        }

        $transcript_description =~ /;parent_transcript=(.+);mapping_coverage=(.+);mapping_identity=(.+)$/;
        my $transcript_versioned_stable_id = $1;
        my $transcript_coverage = $2;
        my $transcript_perc_id = $3;

        # DUE TO BUG IN FINDPARALOGUES NOT ADDING DESCRIPTION
        unless($transcript_versioned_stable_id and defined($transcript_coverage) and defined($transcript_perc_id)) {
          next;
        }

        my $transcript_stable_id = $transcript_versioned_stable_id;
        $transcript_stable_id =~ s/\.\d+//;
        $target_transcript_info->{$transcript_stable_id}->{'transcript_biotype'} = $transcript_biotype;
        $target_transcript_info->{$transcript_stable_id}->{'transcript_biotype_group'} = $transcript_biotype_group;
        $target_transcript_info->{$transcript_stable_id}->{'coverage'} = $transcript_coverage;
        $target_transcript_info->{$transcript_stable_id}->{'perc_id'} = $transcript_perc_id;
      }
    }
  }
}


say "Finished processing target slices for genes and transcripts";

say "Counting source and mapped biotypes for genes and transcripts";
my $total_gene_count_by_biotype = {};
my $mapped_gene_count_by_biotype = {};
my $total_transcript_count_by_biotype = {};
my $mapped_transcript_count_by_biotype = {};
my $mapped_transcripts_count = 0;
my $mapped_transcripts_coverage = 0;
my $mapped_transcripts_percent_id = 0;

# At this point a straightforward comparison is enough
my @missing_gene_info = ();
foreach my $gene_id (keys(%{$source_gene_info})) {
  my $gene_biotype = $source_gene_info->{$gene_id}->{'gene_biotype'};
  my $gene_region_name = $source_gene_info->{$gene_id}->{'gene_region_name'};
  if($total_gene_count_by_biotype->{$gene_biotype}) {
    $total_gene_count_by_biotype->{$gene_biotype}++;
  } else {
    $total_gene_count_by_biotype->{$gene_biotype}  = 1;
  }

  if($target_gene_info->{$gene_id}) {
    if($mapped_gene_count_by_biotype->{$gene_biotype}) {
      $mapped_gene_count_by_biotype->{$gene_biotype}->{'count'}++;
    } else {
      $mapped_gene_count_by_biotype->{$gene_biotype}->{'count'} = 1;
    }
  } else {
    my $missing_string = $gene_id." ".$gene_region_name." ".$gene_biotype;
    push(@missing_gene_info,$missing_string);
  }
}


foreach my $transcript_id (keys(%{$source_transcript_info})) {
  my $transcript_biotype = $source_transcript_info->{$transcript_id}->{'transcript_biotype'};
  if($total_transcript_count_by_biotype->{$transcript_biotype}) {
    $total_transcript_count_by_biotype->{$transcript_biotype}++;
  } else {
    $total_transcript_count_by_biotype->{$transcript_biotype} = 1;
  }

  if($target_transcript_info->{$transcript_id}) {

    $mapped_transcripts_count++;
    $mapped_transcripts_coverage += $target_transcript_info->{$transcript_id}->{'coverage'};
    $mapped_transcripts_percent_id += $target_transcript_info->{$transcript_id}->{'perc_id'};
  
    if($mapped_transcript_count_by_biotype->{$transcript_biotype}) {
      $mapped_transcript_count_by_biotype->{$transcript_biotype}->{'count'}++;
      $mapped_transcript_count_by_biotype->{$transcript_biotype}->{'coverage'} += $target_transcript_info->{$transcript_id}->{'coverage'};
      $mapped_transcript_count_by_biotype->{$transcript_biotype}->{'perc_id'} += $target_transcript_info->{$transcript_id}->{'perc_id'};
    } else {
      $mapped_transcript_count_by_biotype->{$transcript_biotype}->{'count'} = 1;
      $mapped_transcript_count_by_biotype->{$transcript_biotype}->{'coverage'} = $target_transcript_info->{$transcript_id}->{'coverage'};
      $mapped_transcript_count_by_biotype->{$transcript_biotype}->{'perc_id'} = $target_transcript_info->{$transcript_id}->{'perc_id'};
    }

    if($target_transcript_info->{$transcript_id}->{'coverage'} < 90 or $target_transcript_info->{$transcript_id}->{'perc_id'} < 95) {
      if($mapped_transcript_count_by_biotype->{$transcript_biotype}->{'problematic_count'}) {
        $mapped_transcript_count_by_biotype->{$transcript_biotype}->{'problematic_count'}++;
      } else {
        $mapped_transcript_count_by_biotype->{$transcript_biotype}->{'problematic_count'} = 1;
      }
    }
  }
}

say "Finished counting source and mapped biotypes for genes and transcripts";


my @gene_biotyes = sort { $total_gene_count_by_biotype->{$b} <=> $total_gene_count_by_biotype->{$a} }keys %$total_gene_count_by_biotype;
my @transcript_biotyes = sort { $total_transcript_count_by_biotype->{$b} <=> $total_transcript_count_by_biotype->{$a} }keys %$total_transcript_count_by_biotype;

say "Missing gene list:";
for my $missing_gene_string (@missing_gene_info) {
  say MISSING $missing_gene_string;
}

say MAPPING "\nGene mapping stats by biotype:";
print_mapping_stats(\@gene_biotyes,$total_gene_count_by_biotype,$mapped_gene_count_by_biotype);

say MAPPING "\nTranscript mapping stats by biotype:";
print_mapping_stats(\@transcript_biotyes,$total_transcript_count_by_biotype,$mapped_transcript_count_by_biotype);

close INPUT;
close MISSING;
close TRANSLATION;
close MAPPING;


sub print_mapping_stats {
  my ($biotypes,$total_count_by_biotype,$mapped_count_by_biotype) = @_;

  my $overall_source = 0;
  my $overall_mapped = 0;
  my $overall_coverage = 0;
  my $overall_perc_id = 0;
  my $overall_problematic_count = 0;
  my $coverage_biotype_count = 0;
  foreach my $biotype (@$biotypes) {
    my $total = $total_count_by_biotype->{$biotype};
    my $mapped = $mapped_count_by_biotype->{$biotype}->{'count'};
    unless($mapped) {
      $mapped = 0;
    }

    $overall_source += $total;
    $overall_mapped += $mapped;
    my $mapping_percentage = sprintf("%.2f", (($mapped/$total) * 100));
    my $mapping_string =  "  ".$biotype.": ".$mapped."/".$total." (".$mapping_percentage.")";
    if($mapped_count_by_biotype->{$biotype}->{'coverage'}) {
      my $coverage_percent = sprintf("%.2f", ($mapped_count_by_biotype->{$biotype}->{'coverage'}/$mapped));
      my $id_percent = sprintf("%.2f", ($mapped_count_by_biotype->{$biotype}->{'perc_id'}/$mapped));
      $mapping_string .= ", Coverage: ".$coverage_percent."%, Percent id: ".$id_percent."%";
      $overall_coverage += $coverage_percent;
      $overall_perc_id += $id_percent;
      $coverage_biotype_count++;
      my $problematic_mapping_count = 0;
      if($mapped_count_by_biotype->{$biotype}->{'problematic_count'}) {
        $problematic_mapping_count = $mapped_count_by_biotype->{$biotype}->{'problematic_count'};
      }
      $mapping_string .= ", Problematic: ".$problematic_mapping_count;
      $overall_problematic_count += $problematic_mapping_count;
    }
    say MAPPING $mapping_string;
  }

  my $overall_mapping_percent = sprintf("%.2f", (($overall_mapped/$overall_source) * 100));
  my $overall_string = "  Overall mapping: ".$overall_mapped."/".$overall_source." (".$overall_mapping_percent."), ".($overall_source - $overall_mapped)." missing";
  if($overall_coverage) {
    $overall_coverage = sprintf("%.2f", ($mapped_transcripts_coverage/$mapped_transcripts_count));
    $overall_perc_id = sprintf("%.2f", ($mapped_transcripts_percent_id/$mapped_transcripts_count));
    $overall_string .= ", Coverage: ".$overall_coverage."%, Percent id: ".$overall_perc_id."%, Problematic: ".$overall_problematic_count;
  }
  say MAPPING $overall_string;
}


sub filter_genes {
  my ($genes,$source_gene_ids_hash,$xy_scanner) = @_;

  # This should mirror the code in the Minimap2Remap filtering to give comparable numbers
  my $filtered_genes = [];
  foreach my $gene (@$genes) {

    # If xy_scanner has something in it at this point then it means one or both of the chromosomes are missing
    # from the target gene set so we want to skip the genes from the source gene set
    if ($xy_scanner and ($gene->seq_region_name() eq 'X' or $gene->seq_region_name() eq 'Y')) {
      if ($xy_scanner eq 'None') {
        say "Skipping gene  ".$gene->stable_id()." (".$gene->biotype().", ".$gene->seq_region_name."), reason: xy_scanner None";
        next;
      } elsif($xy_scanner eq 'X' and $gene->seq_region_name() eq 'Y') {
        say "Skipping gene  ".$gene->stable_id()." (".$gene->biotype().", ".$gene->seq_region_name."), reason: xy_scanner X and gene seq region name Y";
        next;
      } elsif($xy_scanner eq 'Y' and $gene->seq_region_name() eq 'X') {
        say "Skipping gene  ".$gene->stable_id()." (".$gene->biotype().", ".$gene->seq_region_name."), reason: xy_scanner Y and gene seq region name X";
        next;
      }
    }

    if($gene->seq_region_name() eq 'MT') {
      say "Skipping gene  ".$gene->stable_id()." (".$gene->biotype().", ".$gene->seq_region_name."), reason: MT gene";
      next;
    }

    my $is_readthrough = 0;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      if($is_readthrough) {
         last;
      }

      my $attributes = $transcript->get_all_Attributes();
      foreach my $attribute (@{$attributes}) {
        if($attribute->value eq 'readthrough') {
          $is_readthrough = 1;
          last;
        }
      } # foreach my $attribute
    } # End foreach my $transcript

    unless($is_readthrough) {
      $source_gene_ids_hash->{$gene->stable_id()} = 1;
      say INPUT $gene->stable_id." ".$gene->seq_region_name." ".$gene->biotype;
      push(@$filtered_genes,$gene);
    } else {
      say "Skipping gene  ".$gene->stable_id()." (".$gene->biotype().", ".$gene->seq_region_name."), reason: gene has readthrough transcript";
    }
  }

  return($filtered_genes);
}
