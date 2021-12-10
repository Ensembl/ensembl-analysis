#!/usr/bin/env perl
# Copyright [2021] EMBL-European Bioinformatics Institute
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

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw (attach_Slice_to_Transcript attach_Analysis_to_Transcript calculate_exon_phases);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Getopt::Long qw(:config no_ignore_case);

my $analysis_name = "ensembl";
my $module_name = "Anno";
my $specify_strand;
my $output_path;
my $genome_file;
my $input_gtf_file;
my $output_gtf_file;
my $region_details;
my $compute_translations;

GetOptions(
            'output_path=s'         => \$output_path,
            'input_gtf_file=s'      => \$input_gtf_file,
            'output_gtf_file=s'     => \$output_gtf_file,
            'genome_file=s'         => \$genome_file,
            'region_details=s'      => \$region_details,
            'compute_translations!' => \$compute_translations,
            'analysis_name=s'       => \$analysis_name,
            );

if (!(-e $input_gtf_file)) {
  die "Could not open the GTF file, path used: ".$input_gtf_file;
}

if ($genome_file && -e $genome_file) {
  setup_fasta(-FASTA => $genome_file);
}

unless($region_details =~ /(.+)\.rs(\d+)\.re(\d+)/) {
  die "Issues with the seq region details"
}

my $region_name = $1;
my $region_start = $2;
my $region_end = $3;
my $region_length = $region_end - $region_start + 1;

my $slice = Bio::EnsEMBL::Slice->new(-start             => $region_start,
                                     -end               => $region_end,
                                     -strand            => 1,
                                     -seq_region_name   => $region_name,
                                     -seq_region_length => $region_length);

my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);
my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => $analysis_name,
                                          -module     => $module_name);

my ($genes_to_process,$genes_to_copy) = load_genes($input_gtf_file,$slice,$analysis,%strand_conversion);

say "Genes loaded, have ".scalar(@$genes_to_process)." coding genes to process. A further ".scalar(@$genes_to_copy)." non-coding genes will be copied";

my $cleaned_genes = clean_genes($genes_to_process);
say "Finished cleaning genes. Have ".scalar(@$cleaned_genes)." cleaned genes";

my $all_genes = [@$genes_to_process,@$genes_to_copy,@$cleaned_genes];
my $filtered_genes = filter_lncRNAs($all_genes);

write_to_gtf_file($filtered_genes,$output_gtf_file,$analysis_name);

exit;


sub filter_lncRNAs {
  my ($all_genes) = @_;

  my $coding_genes = [];
  my $lncRNA_genes = [];
  my $filtered_lncRNAs = [];
  my $other_genes = [];
  foreach my $gene (@$all_genes) {
    if($gene->biotype =~ /protein_coding/) {
      push(@$coding_genes,$gene);
    } elsif($gene->biotype =~ /lncRNA/) {
      push(@$lncRNA_genes,$gene);
    } else {
      push(@$other_genes,$gene);
    }
  }

  my $coding_ranges = get_ranges($coding_genes);
  foreach my $gene (@$lncRNA_genes) {
    my $start = $gene->seq_region_start;
    my $end = $gene->seq_region_end;
    my $overlaps = 0;
    foreach my $coord_pair (@$coding_ranges) {
      my $coding_start = ${$coord_pair}[0];
      my $coding_end = ${$coord_pair}[1];
      if (($start <= $coding_end) and ($end >= $coding_start)) {
        $overlaps = 1;
        last;
      }
    }

    unless($overlaps) {
      push(@$filtered_lncRNAs,$gene);
    }
  }

  return([@$coding_genes,@$filtered_lncRNAs,@$other_genes]);
}


sub get_ranges {
  my ($genes) = @_;

  my $ranges = [];
  foreach my $gene (@$genes) {
    push(@$ranges,[$gene->seq_region_start,$gene->seq_region_end]);
  }
  return($ranges);
}


sub clean_genes {
  my ($initial_genes) = @_;

  my $params = { min_size_utr_exon => 30,
                 ratio_5prime_utr => .3,
                 ratio_3prime_utr => .6,
                 ratio_same_transcript => .02,
                 ratio_max_allowed_difference => .05,
                 ratio_expansion => 3,
                 minimum_expanding_number_for_single_transcript => 2,
                 ratio_transcript_fragment => 3,
                 ratio_exon_expansion => 2,
                 ratio_utrs => 2,
                 store_rejected => 0,
                 copy_biotypes_to_ignore => {
                   low_coverage => 1,
                   CRISPR => 1,
                   broken_gene => 1
                 }
               };

  my ($extra_genes, $rejected) = clean_utrs($initial_genes,
                                            $params->{'min_size_utr_exon'},
                                            $params->{'ratio_5prime_utr'},
                                            $params->{'ratio_3prime_utr'},
                                            $params->{'ratio_same_transcript'},
                                            $params->{'ratio_max_allowed_difference'},
                                            $params->{'ratio_expansion'},
                                            $params->{'minimum_expanding_number_for_single_transcript'},
                                            $params->{'ratio_transcript_fragment'},
                                            $params->{'ratio_exon_expansion'},
                                            $params->{'ratio_utrs'},
                                            $params->{'store_rejected'});

  say "Initial: ".scalar(@$initial_genes);
  say "Extra: ".scalar(@$extra_genes);
  say "Rejected: ".scalar(@$rejected);
  return($extra_genes);
}


sub write_to_gtf_file {
  my ($genes_to_write,$output_gtf_file,$analysis_name) = @_;

  my $output_transcript_seq_file = $output_gtf_file.".cdna";
  my $output_transcript_prot_file = $output_gtf_file.".prot";
  open(OUT1,">".$output_gtf_file);
  open(OUT2,">".$output_transcript_seq_file);
  open(OUT3,">".$output_transcript_prot_file);
  my $gene_count = 1;
  my $transcript_count = 1;
  foreach my $gene (@$genes_to_write) {
    my $gene_id = $analysis_name."_".$gene_count;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      my $transcript_id = $analysis_name."_".$transcript_count;
      my $exons = $transcript->get_all_Exons();
      my $translation = $transcript->translation();
      my $record = build_gtf_record($transcript,$exons,$analysis_name,$gene_id,$transcript_id,$translation);
      foreach my $line (@$record) {
        say OUT1 $line;
      }

      say OUT2 ">".$transcript_id."\n".$transcript->seq->seq();

      if($translation) {
        say OUT3 ">".$transcript_id."\n".$translation->seq();
      }
      $transcript_count++;
    }
    $gene_count++;
  }
  close OUT1;
  close OUT2;
  close OUT3;
}


sub build_gtf_record{
  my ($transcript,$exons,$analysis_name,$gene_id,$transcript_id,$translation) = @_;

  my $record = [];
  my $strand = "+";
  if($transcript->strand() == -1) {
    $strand = "-";
  }

  my $transcript_attribs = 'gene_id "'.$gene_id.'"; transcript_id "'.$transcript_id.'"; biotype "'.$transcript->biotype().'";';
  if($translation) {
    # This is a simple encoding for the translation in the output file
    my $start_exon = $translation->start_Exon();
    my $end_exon = $translation->end_Exon();
    my $translation_coords = ' translation_coords "'.$start_exon->start().':'.$start_exon->end().':'.$translation->start().':'.$end_exon->start().':'.$end_exon->end().':'.$translation->end().'";';
    $transcript_attribs .= $translation_coords;
  }

  if($transcript->is_canonical()) {
    $transcript_attribs .= " canonical_transcript;"
  }

  my @transcript_cols = ($transcript->slice->seq_region_name(),$analysis_name,'transcript',$transcript->start(),$transcript->end(),'.',$strand,'.',$transcript_attribs);
  my $transcript_line = join("\t",@transcript_cols);
  push(@$record,$transcript_line);


  my $exon_attribs_generic = 'gene_id "'.$gene_id.'"; transcript_id "'.$transcript_id.'";';
  my $exon_rank = 1;
  foreach my $exon (@$exons) {
    my $exon_attribs = $exon_attribs_generic.' exon_number "'.$exon_rank.'";';
    my @exon_cols = ($transcript->slice->seq_region_name(),$analysis_name,'exon',$exon->start(),$exon->end(),'.',$strand,'.',$exon_attribs);
    my $exon_line = join("\t",@exon_cols);
    push(@$record,$exon_line);
    $exon_rank++;
  }

  return($record);
}

sub load_genes {
  my ($input_gtf_file,$slice,$analysis,%strand_conversion) = @_;

  my $genes_to_process = [];
  my $genes_to_copy = [];
  my $current_gene_id = "";
  my $current_transcript_id = "";
  my $current_translation_coords = "";
  my $current_biotype = "";
  my $exons = [];
  my $transcripts = [];
  my $transcripts_by_gene_id = {};

  say "Reading GTF";
  open(IN,$input_gtf_file);
  while(<IN>) {
    my $line = $_;
    my @eles = split("\t",$line);

    unless(scalar(@eles) == 9) {
      next;
    }

    my $gtf_region = $eles[0];
    unless($gtf_region eq $slice->seq_region_name()) {
      next;
    }

    my $gtf_type = $eles[2];
    my $gtf_start = $eles[3];
    my $gtf_end = $eles[4];
    my $gtf_strand = $strand_conversion{$eles[6]};
    unless($gtf_strand) {
      throw("Issue with parsing strand")
    }

    my $attributes = set_attributes($eles[8]);
    my $gtf_gene_id = $attributes->{'gene_id'};
    my $gtf_transcript_id = $attributes->{'transcript_id'};

    if($gtf_type eq 'transcript' && !$current_transcript_id) {
      # This if for the very first transcript
      $current_transcript_id = $gtf_transcript_id;
      $current_gene_id = $gtf_gene_id;
      my $gtf_biotype = $attributes->{'biotype'};
      my $gtf_translation_coords = $attributes->{'translation_coords'};
      if($gtf_translation_coords) {
        unless($gtf_biotype) {
          $gtf_biotype = 'protein_coding';
        }
        $current_translation_coords = $gtf_translation_coords;
      }

      unless($gtf_biotype) {
        $gtf_biotype = 'not_set';
      }
      $current_biotype = $gtf_biotype;
      next;
    } elsif($gtf_type eq 'transcript') {
      # This is for moving onto a new transcript, process the old one
      if($$exons[0]->strand() == 1) {
        $exons = [sort { $a->start <=> $b->start } @{$exons}];
      } else {
        $exons = [sort { $b->start <=> $a->start } @{$exons}];
      }

      my $transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $exons);
      $transcript->stable_id($current_transcript_id);
      $transcript->biotype($current_biotype);
      $transcript->analysis($analysis);
      $transcript->{'translation_coords'} = $current_translation_coords;

      if($transcripts_by_gene_id->{$current_gene_id}) {
        push(@{$transcripts_by_gene_id->{$current_gene_id}},$transcript);
      } else {
        $transcripts_by_gene_id->{$current_gene_id} = [$transcript];
      }

      # Now set the current variables to the ones from the new transcript
      $current_transcript_id = $gtf_transcript_id;
      $current_gene_id = $gtf_gene_id;
      my $gtf_biotype = $attributes->{'biotype'};
      my $gtf_translation_coords = $attributes->{'translation_coords'};
      if($gtf_translation_coords) {
        unless($gtf_biotype) {
          $gtf_biotype = 'protein_coding';
        }
        $current_translation_coords = $gtf_translation_coords;
      } else {
        $current_translation_coords = "";
      }

      unless($gtf_biotype) {
        $gtf_biotype = 'not_set';
      }
      $current_biotype = $gtf_biotype;
      $exons = [];
    } elsif($gtf_type eq 'exon') {
      # Build an exon and add to the current array of exons
      my $exon = Bio::EnsEMBL::Exon->new(
                                          -START     => $gtf_start,
                                          -END       => $gtf_end,
                                          -STRAND    => $gtf_strand,
                                          -PHASE     => -1,
                                          -END_PHASE => -1);
      $exon->{'region_name'} = $gtf_region;
      push(@$exons,$exon);
    } else {
      throw("Found an unexpected type in the GTF file, expected transcript or exon, found: ".$gtf_type);
    }

  }

  say "Finished reading GTF";
  close IN;

  if(scalar(@$exons)) {
    if($$exons[0]->strand() == 1) {
      $exons = [sort { $a->start <=> $b->start } @{$exons}];
    } else {
      $exons = [sort { $b->start <=> $a->start } @{$exons}];
    }

    my $transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $exons);
    $transcript->stable_id($current_transcript_id);
    $transcript->biotype($current_biotype);
    $transcript->analysis($analysis);
    $transcript->{'translation_coords'} = $current_translation_coords;

    if($transcripts_by_gene_id->{$current_gene_id}) {
      push(@{$transcripts_by_gene_id->{$current_gene_id}},$transcript);
    } else {
      $transcripts_by_gene_id->{$current_gene_id} = [$transcript];
    }
  }


  my $genes_by_region = {};
  foreach my $gene_id (keys(%$transcripts_by_gene_id)) {
    my $transcripts = $transcripts_by_gene_id->{$gene_id};
    my $exons = ${$transcripts}[0]->get_all_Exons();
    my $region_name = ${$exons}[0]->{'region_name'};
    $genes_by_region->{$region_name}->{$gene_id} = $transcripts;
  }

  foreach my $gene_id (keys(%{$genes_by_region->{$region_name}})) {
    my $gene_is_coding = 0;
    my $transcripts = $genes_by_region->{$region_name}->{$gene_id};
    foreach my $transcript (@$transcripts) {
      attach_Slice_to_Transcript($transcript,$slice);
      attach_Analysis_to_Transcript($transcript,$analysis);
      if($transcript->{'translation_coords'}) {
        build_translation($transcript,$transcript->{'translation_coords'});
        $gene_is_coding = 1;
      } elsif($compute_translations) {
        compute_translation($transcript);
        $gene_is_coding = 1;
      }
    }

    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($analysis);
    $gene->stable_id($gene_id);
    foreach my $transcript (@$transcripts) {
      $gene->add_Transcript($transcript);
    }

    if($gene_is_coding) {
      push(@$genes_to_process,$gene);
    } else {
      push(@$genes_to_copy,$gene);
    }
  } # End foreach my $gene_id

  return($genes_to_process,$genes_to_copy);
}



sub sort_genes_by_region {
  my ($genes) = @_;

  my $genes_by_region = {};
  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts();
    my $exons = ${$transcripts}[0]->get_all_Exons();
    my $region_name = ${$exons}[0]->{'region_name'};
    if($genes_by_region->{$region_name}) {
      push(@{$genes_by_region->{$region_name}},$gene);
    } else {
      $genes_by_region->{$region_name} = [$gene];
    }
  }

  return($genes_by_region);
}


sub set_attributes {
  my ($attribute_string) = @_;
  my $attribute_pairs = {};

  my @attribute_array = split(";",$attribute_string);
  foreach my $attribute (@attribute_array) {
    my @pairs = split(" ",$attribute);
    if(scalar(@pairs) == 2) {
      $pairs[1] =~ s/\"//g;
      $attribute_pairs->{$pairs[0]} = $pairs[1];
    } else {
      $attribute =~ s/ //g;
      $attribute_pairs->{$attribute} = 1;
    }
  }

  return($attribute_pairs);
}


sub build_translation {
  my ($transcript,$translation_coords) = @_;

  unless($translation_coords =~ /(\d+)\:(\d+)\:(\d+)\:(\d+)\:(\d+)\:(\d+)/) {
    throw("Issue parsing translation coords, coords string:\n".$translation_coords);
  }

  my $start_exon_start = $1;
  my $start_exon_end = $2;
  my $start_exon_offset = $3;
  my $end_exon_start = $4;
  my $end_exon_end = $5;
  my $end_exon_offset = $6;

  my $exons = $transcript->get_all_Exons();
  my $start_exon;
  my $end_exon;
  foreach my $exon (@$exons) {
    if($exon->start() == $start_exon_start && $exon->end() == $start_exon_end) {
      $start_exon = $exon;
    }

    if($exon->start() == $end_exon_start && $exon->end() == $end_exon_end) {
      $end_exon = $exon;
    }
  }

  unless($start_exon && $end_exon) {
    throw("Could not find matching start/end exon from the translation coords, translation coords:\n".$translation_coords);
  }

  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($start_exon);
  $translation->start($start_exon_offset);
  $translation->end_Exon($end_exon);
  $translation->end($end_exon_offset);
  $transcript->translation($translation);
  calculate_exon_phases($transcript,0);
}
