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
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw (attach_Slice_to_Transcript attach_Analysis_to_Transcript calculate_exon_phases);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Getopt::Long qw(:config no_ignore_case);

my $analysis_name = "ensembl";
my $module_name = "Anno";

my $host;
my $port;
my $user;
my $pass;
my $dbname;
my $dna_host;
my $dna_port;
my $dna_user;
my $dna_dbname;
my $region_details;
my $specify_strand;
my $set_canonical = 1;
my $output_path;
my $genome_file;
my $gtf_file;
my $slice_name;
my $protein_coding_biotype;
my $non_coding_biotype;
my $compute_translations;
my $load_type;
my $make_single_transcript_genes = 0;

GetOptions( 'host|dbhost|h=s'       => \$host,
            'port|dbport|P=s'       => \$port,
            'user|dbuser|u=s'       => \$user,
            'pass|dbpass|p=s'       => \$pass,
            'dbname|db|D=s'     => \$dbname,
            'dna_host=s'   => \$dna_host,
            'dna_port=s'   => \$dna_port,
            'dna_user=s'   => \$dna_user,
            'dna_dbname=s' => \$dna_dbname,
            'specify_strand=s'  => \$specify_strand,
            'output_path=s' => \$output_path,
            'gtf_file=s' => \$gtf_file,
            'genome_file=s' => \$genome_file,
            'slice_name=s' => \$slice_name,
            'compute_translations!' => \$compute_translations,
            'protein_coding_biotype=s' => \$protein_coding_biotype,
            'non_coding_biotype=s' => \$non_coding_biotype,
            'load_type=s' => \$load_type,
            'analysis_name=s' => \$analysis_name,
            'make_single_transcript_genes!' => \$make_single_transcript_genes);

if (!(-e $gtf_file)) {
  die "Could not open the GTF file, path used: ".$gtf_file;
}

if ($genome_file && -e $genome_file) {
  setup_fasta(-FASTA => $genome_file);
}

my $dna_dba;
if ($dna_host and $dna_user and $dna_port and $dna_dbname) {
  $dna_dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      '-host'   => $dna_host,
      '-user'   => $dna_user,
      '-port'   => $dna_port,
      '-dbname' => $dna_dbname,
    );
}

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      '-host'   => $host,
      '-user'   => $user,
      '-pass'   => $pass,
      '-port'   => $port,
      '-dbname' => $dbname,
  );

if ($dna_dba) {
  $dba->dnadb($dna_dba);
}

my $slice_adaptor = $dba->get_SliceAdaptor();
my $slice_hash = {};

my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);
my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => $analysis_name,
                                          -module     => $module_name);

if($load_type eq 'gene') {
  load_genes($dba,$gtf_file,$slice_hash,$analysis,%strand_conversion);
} elsif($load_type eq 'single_line_feature' and ($analysis_name eq 'dust' or $analysis_name eq 'repeatmask_red' or $analysis_name eq 'trf')) {
  load_repeats($dba,$gtf_file,$slice_hash,$analysis,%strand_conversion);
} elsif($load_type eq 'single_line_feature' and ($analysis_name eq 'cpg' or $analysis_name eq 'eponine')) {
  load_simple_features($dba,$gtf_file,$slice_hash,$analysis,%strand_conversion);
}

exit;


sub load_genes {
  my ($dba,$gtf_file,$slice_hash,$analysis,%strand_conversion) = @_;

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


  my $gene_adaptor = $dba->get_GeneAdaptor();
  foreach my $region_name (keys(%$genes_by_region)) {
    my $slice = $slice_adaptor->fetch_by_region('toplevel',$region_name);
    foreach my $gene_id (keys(%{$genes_by_region->{$region_name}})) {
      my $transcripts = $genes_by_region->{$region_name}->{$gene_id};
      foreach my $transcript (@$transcripts) {
        attach_Slice_to_Transcript($transcript,$slice);
        attach_Analysis_to_Transcript($transcript,$analysis);
        if($transcript->{'translation_coords'}) {
          build_translation($transcript,$transcript->{'translation_coords'});
        } elsif($compute_translations) {
          compute_translation($transcript);
        }
      }

      if($make_single_transcript_genes) {
        foreach my $transcript (@$transcripts) {
          if($transcript->biotype() eq 'tRNA_pseudogene') {
            next;
          }
          my $gene = Bio::EnsEMBL::Gene->new();
          $gene->analysis($analysis);
          $gene->stable_id($gene_id);
          $gene->add_Transcript($transcript);
          $gene->biotype($transcript->biotype());
          if($set_canonical) {
            $gene->canonical_transcript($transcript);
          }
          if($protein_coding_biotype || $non_coding_biotype) {
            override_biotypes([$gene],$protein_coding_biotype,$non_coding_biotype);
          }
          empty_Gene($gene);
          $gene_adaptor->store($gene);
          undef($gene);
        }
      } else {
        my $gene = Bio::EnsEMBL::Gene->new();
        $gene->analysis($analysis);
        $gene->stable_id($gene_id);
        foreach my $transcript (@$transcripts) {
          $gene->add_Transcript($transcript);
        }
        if($set_canonical) {
          set_canonical_transcript($gene);
          $gene->biotype($gene->canonical_transcript->biotype());
        }
        if($protein_coding_biotype || $non_coding_biotype) {
          override_biotypes([$gene],$protein_coding_biotype,$non_coding_biotype);
        }
        if($gene->biotype() eq 'tRNA_pseudogene') {
          next;
        }

        empty_Gene($gene);
        $gene_adaptor->store($gene);
        undef($gene);
      }
    } # End foreach my $gene_id
    undef($slice);
  }
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


sub load_repeats {
  my ($dba,$gtf_file,$slice_hash,$analysis,%strand_conversion) = @_;

  my $repeat_types = { trf => 'Tandem repeats',
                       dust => 'Dust',
                     };
  my $repeat_features = [];
  my $features_by_region = {};
  say "Reading GTF";
  open(IN,$gtf_file);
  while(<IN>) {
    my $line = $_;
    my @eles = split("\t",$line);

    unless(scalar(@eles) == 9) {
      next;
    }

    my $gtf_region = $eles[0];
    my $gtf_type = $eles[2];
    my $gtf_start = $eles[3];
    my $gtf_end = $eles[4];
    my $gtf_strand = $strand_conversion{$eles[6]};
    unless($gtf_strand) {
      throw("Issue with parsing strand")
    }

    my $attributes = set_attributes($eles[8]);
    my $repeat_type = $repeat_types->{$analysis->logic_name};
    my $repeat_consensus_seq = 'N';
    if($attributes->{'repeat_consensus'}) {
      $repeat_consensus_seq = $attributes->{'repeat_consensus'};
    }
    my $repeat_score = 0;
    if($attributes->{'score'}) {
      $repeat_score = $attributes->{'score'};
    }
    my $feature_factory = Bio::EnsEMBL::Analysis::Tools::FeatureFactory->new();
    my $repeat_consensus = $feature_factory->create_repeat_consensus($analysis->logic_name,$analysis->logic_name, $repeat_type, $repeat_consensus_seq);
    my $repeat_feature = $feature_factory->create_repeat_feature($gtf_start, $gtf_end, 1, $repeat_score, 1,($gtf_end - $gtf_start + 1),$repeat_consensus);

    if($features_by_region->{$gtf_region}) {
      push(@{$features_by_region->{$gtf_region}},$repeat_feature);
    } else {
      $features_by_region->{$gtf_region} = [$repeat_feature];
    }
  }
  say "Finished reading GTF file";

  my $repeat_feature_adaptor = $dba->get_RepeatFeatureAdaptor();
  foreach my $region_name (keys(%$features_by_region)) {
    my $slice = $slice_adaptor->fetch_by_region('toplevel',$region_name);
    my $features = $features_by_region->{$region_name};
    foreach my $feature (@$features) {
      $feature->analysis($analysis);
      $feature->slice($slice);
      $repeat_feature_adaptor->store($feature);
      undef($feature)
    }
    undef($slice);
  }
}


sub load_simple_features {
  my ($dba,$gtf_file,$slice_hash,$analysis,%strand_conversion) = @_;

  my $simple_features = [];

  my $features_by_region = {};
  say "Reading GTF";
  open(IN,$gtf_file);
  while(<IN>) {
    my $line = $_;
    my @eles = split("\t",$line);

    unless(scalar(@eles) == 9) {
      next;
    }

    my $gtf_region = $eles[0];
    my $gtf_type = $eles[2];
    my $gtf_start = $eles[3];
    my $gtf_end = $eles[4];
    my $gtf_strand = $strand_conversion{$eles[6]};
    unless($gtf_strand) {
      throw("Issue with parsing strand")
    }

    my $attributes = set_attributes($eles[8]);
    my $simple_feature = Bio::EnsEMBL::SimpleFeature->new
    (
     -start => $gtf_start,
     -end => $gtf_end,
     -strand => $gtf_strand,
     -score => 0,
     -display_label => '',
     -seqname => '',
    );

    if($features_by_region->{$gtf_region}) {
      push(@{$features_by_region->{$gtf_region}},$simple_feature);
    } else {
      $features_by_region->{$gtf_region} = [$simple_feature];
    }
  }
  say "Finished reading GTF file";

  my $simple_feature_adaptor = $dba->get_SimpleFeatureAdaptor();
  foreach my $region_name (keys(%$features_by_region)) {
    my $slice = $slice_adaptor->fetch_by_region('toplevel',$region_name);
    my $features = $features_by_region->{$region_name};
    foreach my $feature (@$features) {
      $feature->analysis($analysis);
      $feature->slice($slice);
      $simple_feature_adaptor->store($feature);
      undef($feature)
    }
  }
}


sub override_biotypes {
  my ($final_genes,$protein_coding_biotype,$non_coding_biotype) = @_;

  foreach my $gene (@$final_genes) {
    my $canonical_transcript = $gene->canonical_transcript();
    if($protein_coding_biotype && $canonical_transcript->translation()) {
      $gene->biotype($protein_coding_biotype)
    } elsif($non_coding_biotype && !$canonical_transcript->translation()) {
      $gene->biotype($non_coding_biotype);
    }

    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript(@$transcripts) {
      if($protein_coding_biotype && $transcript->translation()) {
        $transcript->biotype($protein_coding_biotype);
      } elsif($non_coding_biotype && !$transcript->translation()) {
        $transcript->biotype($non_coding_biotype);
      }
    }
  }
}


sub create_single_transcript_genes {
  my ($transcripts) = @_;

  # Creates single transcript genes for use with things like genebuilder
  my $single_transcript_genes = [];
  foreach my $transcript (@$transcripts) {
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->stable_id($transcript->stable_id);
    $gene->biotype($transcript->biotype);
    $gene->slice($transcript->slice());
    $gene->analysis($transcript->analysis);
    $gene->add_Transcript($transcript);
    $gene->canonical_transcript($transcript);
    push(@$single_transcript_genes,$gene);
  }

  return($single_transcript_genes);
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


sub set_canonical_transcript {
  my ($gene) = @_;

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
