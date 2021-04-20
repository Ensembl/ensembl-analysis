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
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases features_overlap clone_Transcript);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(clone_Exon);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Getopt::Long qw(:config no_ignore_case);

my $analysis_name = "ensembl";
my $module_name = "Anno";

my $specify_strand;
my $set_canonical = 1;
my $genome_file;
my $gtf_file;
my $slice_name;
my $compute_translations;
my $output_gtf_file;
GetOptions(
            'specify_strand=s'  => \$specify_strand,
            'gtf_file=s' => \$gtf_file,
            'genome_file=s' => \$genome_file,
            'slice_name=s' => \$slice_name,
            'compute_translations!' => \$compute_translations,
            'output_gtf_file=s' => \$output_gtf_file);


if (!(-e $gtf_file)) {
  die "Could not open the GTF file, path used: ".$gtf_file;
}

if ($genome_file && -e $genome_file) {
  setup_fasta(-FASTA => $genome_file);
}

my $slice_hash = {};
my $slices = fetch_slices_from_fasta($genome_file);
foreach my $slice (@{$slices}) {
  my $seq_region_name = $slice->seq_region_name;
  if($slice_name && $seq_region_name ne $slice_name) {
    next;
  }
  $slice_hash->{$seq_region_name} = $slice;
}


my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);
my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => $analysis_name,
                                          -module     => $module_name);

my $current_gene_id = "";
my $current_transcript_id = "";
my $current_translation_coords = "";
my $current_biotype = "";
my $exons = [];
my $transcripts = [];
my $transcripts_by_gene_id = {};

say "Reading GTF";
open(IN,$gtf_file);
while(<IN>) {
  my $line = $_;
  my @eles = split("\t",$line);

  unless(scalar(@eles) == 9) {
    next;
  }

  my $gtf_region = $eles[0];
  unless($slice_hash->{$gtf_region}) {
    next;
  }

  my $gtf_slice = $slice_hash->{$gtf_region};
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
  unless($transcripts_by_gene_id->{$gtf_gene_id}) {
    $transcripts_by_gene_id->{$gtf_gene_id} = [];
  }

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
    $transcript->slice($$exons[0]->slice());
    $transcript->analysis($analysis);

    if($current_translation_coords) {
      build_translation($transcript,$current_translation_coords);
    } elsif($compute_translations) {
      compute_translation($transcript);
    }

#    push(@$transcripts,$transcript);
    push(@{$transcripts_by_gene_id->{$current_gene_id}},$transcript);

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
                                        -SLICE     => $gtf_slice,
                                        -PHASE     => -1,
                                        -END_PHASE => -1);
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
  $transcript->slice($$exons[0]->slice());
  $transcript->analysis($analysis);

  if($current_translation_coords) {
    build_translation($transcript,$current_translation_coords);
  } elsif($compute_translations) {
    compute_translation($transcript);
  }

  push(@{$transcripts_by_gene_id->{$current_gene_id}},$transcript);
}

my $genes = [];
foreach my $gene_id (keys(%$transcripts_by_gene_id)) {
  my $transcripts = $transcripts_by_gene_id->{$gene_id};
  my $new_gene = Bio::EnsEMBL::Gene->new();
  $new_gene->strand(${$transcripts}[0]->strand());
  $new_gene->biotype(${$transcripts}[0]->biotype());
  $new_gene->slice(${$transcripts}[0]->slice());
  $new_gene->analysis(${$transcripts}[0]->analysis());
  foreach my $transcript (@$transcripts) {
    $new_gene->add_Transcript($transcript);
  }
  push(@$genes,$new_gene);
}

set_canonical_transcripts($genes);

my $transcript_seq_file = $gtf_file.".cdna.fa";
my $transcript_canonical_seq_file = $gtf_file.".canonical.cdna.fa";
my $protein_seq_file = $gtf_file.".prot.fa";
my $protein_canonical_seq_file = $gtf_file.".canonical.prot.fa";

open(OUT1,">".$transcript_seq_file);
open(OUT2,">".$transcript_canonical_seq_file);
open(OUT3,">".$protein_seq_file);
open(OUT4,">".$protein_canonical_seq_file);
foreach my $gene (@$genes) {
  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@$transcripts) {
    my $transcript_seq = $transcript->seq->seq();
    say OUT1 ">".$transcript->stable_id."\n".$transcript_seq;
    if($transcript->is_canonical()) {
      say OUT2 ">".$transcript->stable_id."\n".$transcript_seq;
    }

    if($transcript->translation()) {
      my $translation_seq = $transcript->translation->seq();
      say OUT3 ">".$transcript->stable_id."\n".$translation_seq;
      if($transcript->is_canonical()) {
        say OUT4 ">".$transcript->stable_id."\n".$translation_seq;
      }
    }
  } # End foreach my $transcript
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;

exit;


sub sort_genes_by_strand {
  my ($genes) = @_;

  my $forward_genes = [];
  my $reverse_genes = [];
  foreach my $gene (@$genes) {
    if($gene->strand == 1) {
      push(@$forward_genes,$gene);
    } else {
      push(@$reverse_genes,$gene);
    }
  }

  my $sorted_forward_genes = [(sort { $a->start() <=> $b->start() } @$forward_genes)];
  my $sorted_reverse_genes = [(sort { $a->start() <=> $b->start() } @$reverse_genes)];

  return($sorted_forward_genes,$sorted_reverse_genes);
}



sub fetch_slices_from_fasta {
  my ($genome_file) = @_;

  my @slice_info = ();
  my $slices = [];
  my $current_slice_name = "";
  my $current_seq = "";
  open(IN,$genome_file);
  while(<IN>) {
    my $line = $_;
    chomp($line);
    if($line =~ /\>(.+)/ && $current_slice_name) {
      my $new_slice_name = $1;
      my $slice_length = length($current_seq);
      push(@slice_info,[$current_slice_name,$slice_length]);
      $current_seq = "";
      $current_slice_name = $new_slice_name;
    } elsif($line =~ /\>(.+)/) {
      $current_slice_name = $1;
    } else {
      $current_seq .= $line;
    }
  }
  close IN;

  if($current_seq) {
    my $slice_length = length($current_seq);
    push(@slice_info,[$current_slice_name,$slice_length]);
  }

  foreach my $slice_details (@slice_info) {
    my $slice = Bio::EnsEMBL::Slice->new(-seq_region_name   => ${$slice_details}[0],
                                         -start             => 1,
                                         -end               => ${$slice_details}[1],
                                         -strand            => 1,
                                         -seq_region_length => ${$slice_details}[1]);
    push(@$slices,$slice);
  }
  return($slices);
}


sub create_single_transcript_genes {
  my ($transcripts) = @_;

  # Creates single transcript genes for use with things like genebuilder
  my $single_transcript_genes = [];
  foreach my $transcript (@$transcripts) {
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->biotype($transcript->biotype());
    $gene->slice($transcript->slice());
    $gene->analysis($transcript->analysis());
    $gene->add_Transcript($transcript);
    push(@$single_transcript_genes,$gene);
  }

  return($single_transcript_genes);
}


sub sort_features_by_slice {
  my ($features) = @_;

  my $features_by_slice = {};
  foreach my $feature (@$features) {
    my $slice_name = $feature->slice->name();
    unless(exists($features_by_slice->{$slice_name})) {
      $features_by_slice->{$slice_name} = [];
    }
    push(@{$features_by_slice->{$slice_name}},$feature);
  }

  return($features_by_slice);
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


sub set_canonical_transcripts {
  my ($genes) = @_;

  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts();
    my $current_canonical = pop(@$transcripts);

    foreach my $transcript (@{$gene->get_all_Transcripts()}) {
      my $translation = $transcript->translation();
      if($translation) {
        if(!($current_canonical->translation())) {
          $current_canonical = $transcript;
        } elsif($current_canonical->translation->length() < $translation->length()) {
          $current_canonical = $transcript;
        } elsif($current_canonical->translation->length() == $translation->length() && $current_canonical->length() < $transcript->length()) {
          $current_canonical = $transcript;
        }
      } elsif(!($current_canonical->translation()) && $current_canonical->length() < $transcript->length()) {
        $current_canonical = $transcript;
      }
    }
    $gene->canonical_transcript($current_canonical);
  }
  return($genes);
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
