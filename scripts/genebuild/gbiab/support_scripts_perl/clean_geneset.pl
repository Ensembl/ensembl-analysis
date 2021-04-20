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

my $host;
my $port;
my $user;
my $pass;
my $dbname;
my $dna_host;
my $dna_port;
my $dna_user;
my $dna_dbname;
my $specify_strand;
my $set_canonical = 1;
my $genome_file;
my $gtf_file;
my $slice_name;
my $compute_translations;
my $output_gtf_file;
GetOptions(
            'host|dbhost|h=s'       => \$host,
            'port|dbport|P=s'       => \$port,
            'user|dbuser|u=s'       => \$user,
            'pass|dbpass|p=s'       => \$pass,
            'dbname|db|D=s'     => \$dbname,
            'dna_host=s'   => \$dna_host,
            'dna_port=s'   => \$dna_port,
            'dna_user=s'   => \$dna_user,
            'dna_dbname=s' => \$dna_dbname,
            'specify_strand=s'  => \$specify_strand,
            'gtf_file=s' => \$gtf_file,
            'genome_file=s' => \$genome_file,
            'slice_name=s' => \$slice_name,
            'compute_translations!' => \$compute_translations,
            'output_gtf_file=s' => \$output_gtf_file);

#my $dna_dba;
#if ($dna_host and $dna_user and $dna_port and $dna_dbname) {
#  $dna_dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
#      '-host'   => $dna_host,
#      '-user'   => $dna_user,
#      '-port'   => $dna_port,
#      '-dbname' => $dna_dbname,
#    );
#}

#my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
#      '-host'   => $host,
#      '-user'   => $user,
#      '-pass'   => $pass,
#      '-port'   => $port,
#      '-dbname' => $dbname,
#  );

#if ($dna_dba) {
#  $dba->dnadb($dna_dba);
#}

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

#my $slice_adaptor = $dba->get_SliceAdaptor();
#if($slice_name) {
#  my $slice = $slice_adaptor->fetch_by_region('toplevel',$slice_name);
#  say "Slice name: ".$slice->name;
#  $slice_hash->{$slice_name} = $slice;
#} else {
#  my $slice_array = $slice_adaptor->fetch_all('toplevel');
#  foreach my $slice (@{$slice_array}) {
#    my $seq_region_name = $slice->seq_region_name;
#    $slice_hash->{$seq_region_name} = $slice;
#  }
#}

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

    push(@$transcripts,$transcript);

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

  push(@$transcripts,$transcript);
}

my $final_genes = [];
my $single_transcript_genes = create_single_transcript_genes($transcripts);
my $single_transcript_genes_by_slice = sort_features_by_slice($single_transcript_genes);
foreach my $slice_name (keys(%$single_transcript_genes_by_slice)) {
  my $slice_genes = $single_transcript_genes_by_slice->{$slice_name};
  my $slice = $slice_hash->{$slice_name};
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
                     -query => $slice,
                     -analysis => $analysis,
                     -genes => $slice_genes,
                     -output_biotype => 'gbuild',
                     -max_transcripts_per_cluster => 100,
                     -min_short_intron_len => 7,
                     -max_short_intron_len => 15,
                     -blessed_biotypes => {},
                     -skip_readthrough_check => 1,
                     -max_exon_length => 50000,
                     -coding_only => 1,
                   );
  $runnable->run();

  my $collapsed_slice_genes = $runnable->output();
  my $genes_with_unique_ids = set_ids($collapsed_slice_genes);
  my $genes_with_cannonicals = set_canonical_transcripts($genes_with_unique_ids);
  my $genes_with_biotypes = set_gene_biotypes($genes_with_cannonicals);
  my $post_lncrna_filtered_genes = remove_overlapping_lncrnas($genes_with_biotypes);
  my $collapsed_lncrna_genes = collapse_lncrna_genes($post_lncrna_filtered_genes,$analysis);
  my $genes_without_readthroughs = remove_potential_readthroughs($collapsed_lncrna_genes,$analysis);
  push(@$final_genes,@$genes_without_readthroughs);
#  my $genes_cleanded_utrs = clean_utrs($genes_without_readthroughs);
#  push(@$final_genes,@$genes_cleanded_utrs);
}

#my $gene_adaptor = $dba->get_GeneAdaptor();
#foreach my $gene (@$final_genes) {
#  $gene->analysis($analysis);
#  empty_Gene($gene);
#  $gene_adaptor->store($gene);
#}

$final_genes = set_ids($final_genes);
write_to_gtf_file($final_genes,$output_gtf_file,$analysis_name);

exit;


sub write_to_gtf_file {
  my ($genes_to_write,$output_gtf_file,$analysis_name) = @_;

  open(OUT1,">".$output_gtf_file);
  my $gene_count = 1;
  my $transcript_count = 1;
  foreach my $gene (@$genes_to_write) {
    my $gene_id = "gene_".$gene_count;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      my $transcript_id = "transcript_".$transcript_count;
      my $exons = $transcript->get_all_Exons();
      my $translation = $transcript->translation();
      my $record = build_gtf_record($transcript,$exons,$analysis_name,$gene_id,$transcript_id,$translation);
      foreach my $line (@$record) {
        say OUT1 $line;
      }
      $transcript_count++;
    }
    $gene_count++;
  }
  close OUT1;
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


#sub set_attributes {
#  my ($attribute_string) = @_;
#  my $attribute_pairs = {};

#  my @attribute_array = split(";",$attribute_string);
#  foreach my $attribute (@attribute_array) {
#    my @pairs = split(" ",$attribute);
#    if(scalar(@pairs) == 2) {
#      $pairs[1] =~ s/\"//g;
#      $attribute_pairs->{$pairs[0]} = $pairs[1];
#    }
#  }

#  return($attribute_pairs);
#}


sub set_ids {
  my ($genes) = @_;
  # Note: This sets unique gene/transcript stable ids. This is used initially on the slices to help process genes
  # and then re-applied on the final gene set to make completely unique ids
  my $gene_id_index = 1;
  my $transcript_id_index = 1;
  foreach my $gene (@$genes) {
    $gene->stable_id("gene_".$gene_id_index);
    $gene_id_index++;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      $transcript->stable_id("transcript_".$transcript_id_index);
      $transcript_id_index++;
    }
  }
  return($genes);
}


sub set_gene_biotypes {
  my ($genes) = @_;

  # Note: This runs off the assumption that the input is just lncRNAs and protein coding genes
  foreach my $gene (@$genes) {
    $gene->biotype('lncRNA');
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      if($transcript->translation()) {
        $gene->biotype('protein_coding');
        $transcript->biotype('protein_coding');
      } else {
        $transcript->biotype('lncRNA');
      }
    }
  }
  return($genes);
}


sub collapse_lncrna_genes {
  my ($genes,$analysis) = @_;
  my $collapsed_genes = [];
  my $lncrna_genes = [];
  foreach my $gene (@$genes) {
    unless($gene->biotype() eq 'lncRNA') {
      push(@$collapsed_genes,$gene);
      next;
    }
    push(@$lncrna_genes,$gene);
  }

  unless(scalar(@$lncrna_genes)) {
    return($collapsed_genes);
  }

  my $slice = ${$lncrna_genes}[0]->slice();
  my $biotype = ${$lncrna_genes}[0]->biotype();
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
                   -query => $slice,
                   -analysis => $analysis,
                   -genes => $lncrna_genes,
                   -output_biotype => $biotype,
                   -max_transcripts_per_cluster => 100,
                   -min_short_intron_len => 7,
                   -max_short_intron_len => 15,
                   -blessed_biotypes => {},
                   -skip_readthrough_check => 1,
                   -max_exon_length => 50000,
                 );
  $runnable->run();
  my $collapsed_lncrna_genes = $runnable->output();

  push(@$collapsed_genes,@$collapsed_lncrna_genes);
  return($collapsed_genes);
}


sub remove_potential_readthroughs {
  my ($genes,$analysis) = @_;

  my $max_coding_exon_skip = 3;
  my $filtered_genes = [];
  foreach my $gene (@$genes) {
    unless($gene->biotype() eq 'protein_coding') {
      push(@$filtered_genes,$gene);
      next;
    }

    my $canonical_transcript = $gene->canonical_transcript();
    my $filtered_transcripts = [];
    my $transcripts = $gene->get_all_Transcripts();
    my $coding_exons = get_unique_coding_exons($transcripts);
    foreach my $transcript (@$transcripts) {
      my $remove_transcript = 0;
      my $introns = $transcript->get_all_Introns();
      foreach my $intron (@$introns) {
        my $skipped_coding_exon_count = 0;
        foreach my $coding_exon (@$coding_exons) {
          if(features_overlap($intron,$coding_exon)) {
            $skipped_coding_exon_count++;
            if($skipped_coding_exon_count > $max_coding_exon_skip) {
              $remove_transcript = 1;
              last;
            }
          }
        }
        if($skipped_coding_exon_count > $max_coding_exon_skip) {
          last;
        }
      } # End foreach my $intron
      unless($remove_transcript) {
        push(@$filtered_transcripts,$transcript);
      }
    } # End foreach my $transcript

    if(scalar(@$filtered_transcripts) == 0) {
      push(@$filtered_transcripts,$canonical_transcript);
    }

    unless(scalar(@$filtered_transcripts) == scalar(@$transcripts)) {
      my $single_transcript_genes = create_single_transcript_genes($filtered_transcripts);
      my $slice = ${$single_transcript_genes}[0]->slice();
      my $biotype = ${$single_transcript_genes}[0]->biotype();
      my $analysis = ${$single_transcript_genes}[0]->analysis();
      my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
                     -query => $slice,
                     -analysis => $analysis,
                     -genes => $single_transcript_genes,
                     -output_biotype => $biotype,
                     -max_transcripts_per_cluster => 100,
                     -min_short_intron_len => 7,
                     -max_short_intron_len => 15,
                     -blessed_biotypes => {},
                     -skip_readthrough_check => 1,
                     -max_exon_length => 50000,
                     -coding_only => 1,
                   );
      $runnable->run();
      my $collapsed_genes = $runnable->output();
      push(@$filtered_genes,@$collapsed_genes);
    } else {
      push(@$filtered_genes,$gene);
    }
  }

  return($filtered_genes);
}


sub get_unique_coding_exons {
  my ($transcripts) = @_;

  my $seen_coding_exons = {};
  my $coding_exons = [];
  foreach my $transcript (@$transcripts) {
    my $cds_exons = $transcript->get_all_CDS();
    foreach my $cds_exon (@$cds_exons) {
      my $cds_exon_string = $cds_exon->start.":".$cds_exon->end;
      unless($seen_coding_exons->{$cds_exon_string}) {
        push(@$coding_exons,$cds_exon);
        $seen_coding_exons->{$cds_exon_string} = 1;
      }
    }
  }
  return($coding_exons);
}


sub remove_overlapping_lncrnas {
  my ($genes) = @_;

  my ($forward_genes,$reverse_genes) = sort_genes_by_strand($genes);

  my $filtered_forward_genes = remove_overlapping_lncrnas_by_strand($forward_genes);
  my $filtered_reverse_genes = remove_overlapping_lncrnas_by_strand($reverse_genes);

  return([@$filtered_forward_genes,@$filtered_reverse_genes]);
}


sub remove_overlapping_lncrnas_by_strand {
  my ($stranded_genes) = @_;

  my $gene_ids_to_remove = {};
  my $filtered_genes = [];
  for(my $i=0; $i<scalar(@$stranded_genes)-1; $i++) {
    my $gene_left = ${$stranded_genes}[$i];
    if($gene_ids_to_remove->{$gene_left->stable_id}) {
      next;
    }

    my $biotype_left = $gene_left->biotype();
    for(my $j = $i+1; $j<scalar(@$stranded_genes); $j++) {
      my $gene_right = ${$stranded_genes}[$j];
      if($gene_ids_to_remove->{$gene_right->stable_id}) {
        next;
      }

      my $biotype_right = $gene_right->biotype();
      unless(features_overlap($gene_left,$gene_right)) {
        last;
      }

      if(features_overlap($gene_left,$gene_right)) {
        if($biotype_left eq $biotype_right) {
          next;
        } if($biotype_left eq 'lncRNA') {
          $gene_ids_to_remove->{$gene_left->stable_id()} = 1;
        } else {
          $gene_ids_to_remove->{$gene_right->stable_id()} = 1;
        }
      }
    } # for(my $j = $i+1;
  } # End for(my $i=0;

  foreach my $gene (@$stranded_genes) {
    unless($gene_ids_to_remove->{$gene->stable_id()}) {
      push(@$filtered_genes,$gene);
    }
  }
  return($filtered_genes);
}


sub clean_utrs {
  my ($genes) = @_;

  my $cleaned_utr_genes = [];
  my $genes_to_clean = [];
  foreach my $gene (@$genes) {
    if($gene->biotype() eq 'protein_coding') {
      push(@$genes_to_clean,$gene);
    } else {
      push(@$cleaned_utr_genes,$gene);
    }
  }

  my ($forward_genes,$reverse_genes) = sort_genes_by_strand($genes_to_clean);

  my $filtered_forward_genes = clean_utrs_by_strand($forward_genes);
  my $filtered_reverse_genes = clean_utrs_by_strand($reverse_genes);

  return([@$filtered_forward_genes,@$filtered_reverse_genes]);
}


sub clean_utrs_by_strand {
  my ($stranded_genes) = @_;

  my $updated_genes = [];
  for(my $i=0; $i<scalar(@$stranded_genes)-1; $i++) {
    my $gene_left = ${$stranded_genes}[$i];
    my $transcripts_left = $gene_left->get_all_Transcripts();
    my $cds_feature_left = generate_cds_feature($transcripts_left);

    for(my $j = $i+1; $j<scalar(@$stranded_genes); $j++) {
      my $gene_right = ${$stranded_genes}[$j];
      unless(features_overlap($gene_left,$gene_right)) {
        last;
      }

      my $transcripts_right = $gene_right->get_all_Transcripts();
      my $cds_feature_right = generate_cds_feature($transcripts_right);

      my $gene_left_updated = 0;
      my $gene_right_updated = 0;

      my $new_gene_left = remove_overlapping_utr_exons($transcripts_left,$cds_feature_right,$gene_left);
      my $new_gene_right = remove_overlapping_utr_exons($transcripts_right,$cds_feature_left,$gene_right);

      if($new_gene_left) {
        ${$stranded_genes}[$i] = $new_gene_left;
        $transcripts_left = $new_gene_left->get_all_Transcripts();
      }

      if($new_gene_right) {
        ${$stranded_genes}[$j] = $new_gene_right;
      }
    } # for(my $j = $i+1;
  } # End for(my $i=0;

  return($stranded_genes);
}


sub remove_overlapping_utr_exons {
  my ($transcripts,$cds_feature,$gene) = @_;

  # This builds a new gene and builds transcripts removing and exons that are UTR that overlap the cds region of the other gene
  my $new_gene = Bio::EnsEMBL::Gene->new();
  $new_gene->strand(${$transcripts}[0]->strand());
  $new_gene->biotype(${$transcripts}[0]->biotype());
  $new_gene->slice(${$transcripts}[0]->slice());
  $new_gene->analysis(${$transcripts}[0]->analysis());

  my $updated_gene = 0;
  foreach my $transcript (@$transcripts) {
    my $updated_transcript = 0;
    unless(features_overlap($transcript,$cds_feature)) {
      my $cloned_transcript = clone_Transcript($transcript);
      $new_gene->add_Transcript($transcript);
    }

    my $exons = $transcript->get_all_Exons();
    my $updated_exons = [];
    say "F1!!!!!!!!!!!!!!1";
    foreach my $exon (@$exons) {
      # Push and coding or partially coding exons onto the updated pile
      say "Exon: ".$exon->phase." ".$exon->end_phase;
      unless($exon->phase == -1 && $exon->end_phase == -1) {
        my $cloned_exon = clone_Exon($exon);
        push(@$updated_exons,$cloned_exon);
      }

      # If the exon doesn't overlap with the cds then push it onto the updated_pipe
      # Otherwise skip it and mark both the gene and the transcript as having been updated
      unless(features_overlap($exon,$cds_feature)) {
        my $cloned_exon = clone_Exon($exon);
        push(@$updated_exons,$cloned_exon);
      } else {
        $updated_gene = 1;
        $updated_transcript = 1;
      }
    }

    if($updated_transcript) {
      unless(scalar(@$updated_exons)) {
        throw("ISSUE WITH EXONS");
      }
      my $new_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $updated_exons);
      $new_transcript->biotype(${$transcripts}[0]->strand());
      $new_transcript->biotype(${$transcripts}[0]->biotype());
      $new_transcript->slice(${$transcripts}[0]->slice());
      $new_transcript->analysis(${$transcripts}[0]->analysis());
      $new_transcript->translation(${$transcripts}[0]->translation());
#      my $cloned_start_exon = clone_Exon($transcript->translation->start_Exon());
#      my $cloned_end_exon = clone_Exon($transcript->translation->end_Exon());
#      my $translation = Bio::EnsEMBL::Translation->new();
#      $translation->start_Exon($cloned_start_exon);
#      $translation->start($transcript->translation->start());
#      $translation->end_Exon($cloned_end_exon);
#      $translation->end($transcript->translation->end());
#      $new_transcript->translation($translation);
#      calculate_exon_phases($new_transcript,0);

 #     calculate_exon_phases($new_transcript,0);
      $new_gene->add_Transcript($new_transcript);
      say "TRANS: ".$new_transcript->translation->seq;
    } else {
      $new_gene->add_Transcript($transcript);
    }
  }

  if($updated_gene) {
    say "ORIGINAL GENE: ".$gene->start." ".$gene->end." ".$gene->strand;
    say "NEW GENE: ".$new_gene->start." ".$new_gene->end." ".$new_gene->strand;
    my $new_transcripts = $new_gene->get_all_Transcripts();
    foreach my $transcript (@$new_transcripts) {
      say "TRANSCRIPT STRAND: ".$transcript->strand." ".$transcript->start()." ".$transcript->end();
    }
#    throw("TEST");
    return($new_gene);
  } else {
    return(0);
  }
}


sub generate_cds_feature {
  my ($transcripts) = @_;

  # This generate an exon that spans the entire cds region covered by the transcripts
  # in a gene
  my $min_start;
  my $max_end;
  foreach my $transcript (@$transcripts) {
    my $cds_exons = $transcript->get_all_CDS();
    foreach my $cds_exon (@$cds_exons) {
      my $start = $cds_exon->start();
      my $end = $cds_exon->end();
      if(!$min_start || $start < $min_start) {
        $min_start = $start;
      }

      if(!$max_end || $end > $max_end) {
        $max_end = $end;
      }
    }
  }

  my $cds_feature = Bio::EnsEMBL::Exon->new(
                                            -START     => $min_start,
                                            -END       => $max_end,
                                            -STRAND    => ${$transcripts}[0]->strand(),
                                            -SLICE     => ${$transcripts}[0]->slice(),
                                            -PHASE     => -1,
                                            -END_PHASE => -1);
  return($cds_feature);

}


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
