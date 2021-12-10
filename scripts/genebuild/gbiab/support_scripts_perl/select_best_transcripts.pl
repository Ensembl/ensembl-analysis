use warnings;
use strict;
use feature 'say';
use List::Util qw( min max );

use File::Spec::Functions;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(clone_Exon);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases clone_Transcript features_overlap);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(contains_internal_stops compute_translation);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;

my $host;
my $port;
my $user;
my $pass;
my $dbname;
my $dna_host;
my $dna_port;
my $dna_user;
my $dna_dbname;

my $analysis_name = "anno";
my $module_name = "anno";

my $region_details = '';
my $specify_strand;
my $good_biotype = 'transcriptomic';
my $bad_biotype = 'transcriptomic_flagged';
my $all_cds_exons = 0;
my $join_transcripts = 0;
my $clean_transcripts = 0;
my $cds_search = 0;
my $write_protein_seqs = 1;
my $write_transcript_seqs = 0;
my $write_single_transcript_genes = 1;
my $set_canonical = 1;
my $skip_db_write = 0;
my $output_path;

my $genome_file = '';
my $input_gtf_file = '';
my $output_gtf_file = '';
my $final_biotype = 'not_set';

GetOptions( 'gtf_file=s'             => \$input_gtf_file,
            'region_details=s'       => \$region_details,
            'specify_strand=s'       => \$specify_strand,
            'cds_search!'            => \$cds_search,
            'write_protein_seqs!'    => \$write_protein_seqs,
            'write_transcript_seqs!' => \$write_transcript_seqs,
            'skip_db_write!'         => \$skip_db_write,
            'output_path=s'          => \$output_path,
            'input_gtf_file=s'       => \$input_gtf_file,
            'output_gtf_file=s'      => \$output_gtf_file,
            'genome_file=s'          => \$genome_file,
            'final_biotype=s'        => \$final_biotype,
            'join_transcripts!'      => \$join_transcripts,
            'clean_transcripts!'     => \$clean_transcripts,
            'all_cds_exons!'         => \$all_cds_exons);

print "FROM GBIAB\n";
print "$input_gtf_file\n$output_gtf_file\n$region_details\n";
if (!(-e $input_gtf_file)) {
  die "Could not open the GTF file, path used: ".$input_gtf_file;
}

setup_fasta(-FASTA => $genome_file);

unless($region_details =~ /(.+)\.rs(\d+)\.re(\d+)/) {
  die "Issues with the seq region details"
}

my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => $analysis_name,
                                          -module     => $module_name);

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

open(IN,$input_gtf_file);
my $current_gene_id = "";
my $current_transcript_id = "";
my $current_gtf_region_name = "";
my $exons = [];
my $transcripts = [];
my $genes = [];
my $exon;
my $gtf_region_name;
my $new_record = 0;
my $biotype = 'not_set';

say "Reading in GTF";
while (<IN>) {
  my $current_line = $_;

  if($current_line =~ /^\#/) {
    next;
  }

  my @columns = split("\t",$current_line);
  $gtf_region_name = $columns[0];

  if($gtf_region_name ne $region_name) {
    next;
  }

  my $type = $columns[2];
  if($type eq 'transcript' && $current_transcript_id) {
    $new_record = 1;
    next;
  }

  if($type ne 'exon') {
    next;
  }

  my $start = $columns[3];
  my $end = $columns[4];
  my $strand = $strand_conversion{$columns[6]};

  unless($strand) {
    $strand = 1;
  }


  if($specify_strand && $specify_strand != $strand) {
    next;
  }

  my $phase = ($columns[7] =~ /\./) ? undef : $columns[7];
  my $attributes = set_attributes($columns[8]);

  my $gene_id = $attributes->{'gene_id'};
  my $transcript_id = $attributes->{'transcript_id'};

  if($attributes->{'biotype'}) {
    $biotype = $attributes->{'biotype'}
  }

  $exon = Bio::EnsEMBL::Exon->new(
                                      -START     => $start,
                                      -END       => $end,
                                      -STRAND    => $strand,
                                      -SLICE     => $slice,
                                      -PHASE     => -1,
                                      -END_PHASE => -1);

  # This is weak if the transcript id is not unique
  if($new_record) {
    $new_record = 0;
    if($$exons[0]->strand == 1) {
      $exons = [sort { $a->start <=> $b->start } @{$exons}];
    } else {
      $exons = [sort { $b->start <=> $a->start } @{$exons}];
    }

    my $transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $exons);
    $transcript->stable_id($current_transcript_id);
    $transcript->biotype($biotype);
    $transcript->slice($slice);
    compute_translation($transcript);
    push(@$transcripts,$transcript);
    $exons = [];
    push(@$exons,$exon);
  } else {
    push(@$exons,$exon);
  }
  $current_transcript_id = $transcript_id;
  $current_gtf_region_name = $gtf_region_name;
}


if(scalar(@$exons)) {
  if($$exons[0]->strand == 1) {
    $exons = [sort { $a->start <=> $b->start } @{$exons}];
  } else {
    $exons = [sort { $b->start <=> $a->start } @{$exons}];
  }

  my $transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $exons);
  $transcript->stable_id($current_transcript_id);
  $transcript->biotype($biotype);
  $transcript->slice($slice);
#  compute_translation($transcript);
  push(@$transcripts,$transcript);
}

say "Finished reading GTF";

#foreach my $transcript (@$transcripts) {
#  if($transcript->translation) {
#    say ">".$transcript->stable_id;
#    say $transcript->translation->seq();
#    say "Transcript seq:\n".$transcript->seq->seq();
#    say "Translation seq:\n".$transcript->translation->seq();
#  }
#}


my $transcripts_by_slice = sort_transcripts_by_slice($transcripts);
foreach my $slice_name (keys(%{$transcripts_by_slice})) {
  say "Processing region: ".$slice_name;

  my $initial_transcripts = $transcripts_by_slice->{$slice_name};
  my $sorted_transcripts = remove_overlapping_exons($initial_transcripts);

  # We want to keep the original set of transcripts for processing, but make a copy
  # for modifying/merging
  my $cloned_transcripts = [];
  foreach my $transcript(@$sorted_transcripts) {
    $transcript->biotype('orig');
    push(@$cloned_transcripts,clone_Transcript($transcript));
  }

  my $joined_transcripts;
  if($join_transcripts) {
    $joined_transcripts = join_transcripts($cloned_transcripts);
    push(@$joined_transcripts,@$sorted_transcripts);
    say "Transcript count after joining: ".scalar(@{$joined_transcripts});
  } else {
    $joined_transcripts = $sorted_transcripts;
  }


  say "Computing translations";
  foreach my $transcript (@$joined_transcripts) {
    if($all_cds_exons) {
      my $exons = $transcript->get_all_Exons();
      my $start_exon = $$exons[0];
      my $end_exon = $$exons[scalar(@$exons)-1];
      my $translation = Bio::EnsEMBL::Translation->new();
      $translation->start_Exon($start_exon);
      $translation->start(1);
      $translation->end_Exon($end_exon);
      $translation->end($end_exon->length());
      $transcript->translation($translation);
      calculate_exon_phases($transcript,0);
    } else {
      compute_translation($transcript);
    }
  }


  my $cleaned_transcripts;
  if($clean_transcripts) {
    say "Cleaning transcripts";
    $cleaned_transcripts = clean_initial_transcripts($joined_transcripts);
  } else {
    $cleaned_transcripts = $joined_transcripts;
  }

  say "Transcripts post initial cleaning: ".scalar(@$cleaned_transcripts);
  say "Processing transcripts to look for additional ORFs in UTR regions";

  my $processed_transcripts = [];
  if($cds_search) {
    my $processed_exon_strings = {};
    # At this point the transcripts array has all the transcripts and each has a translation
    # Now we want to reprocess them to find additional ORFs
    # Iteratively, UTRs are made into new transcripts and added back to the pile
    # This reprocessing of UTRs continues until and of the following stop conditions are met:
    # 1: The UTR length is < 350bp
    # 2: The longest ORF from the computed translation is < 100aa
    # 3: The ORF does not start with a met or end with a stop
    while(scalar(@$cleaned_transcripts)) {
      # Note the original transcript should be kept by process_transcript, as each unique exon string it encounters
      # is added to $processed_transcripts. So nothing should be lost as such
      my $transcript = pop(@$cleaned_transcripts);
      process_transcript($transcript,$cleaned_transcripts,$processed_transcripts,$processed_exon_strings);
    }
    say "Transcripts after scanning UTRs for potential ORFs: ".scalar(@$processed_transcripts);
    # Should consider code here to add back UTR by substiting transcripts from the original set
    # that had UTR and the same CDS with the transcripts in processed transcript
    # A UTR cleaner later could cut UTR off overlapping genes
  } else {
    $processed_transcripts = $cleaned_transcripts;
  }

  my $initial_genes = create_single_transcript_genes($processed_transcripts);
  say "Initial genes: ".scalar(@$initial_genes);
  my $select_genes = select_geneset($initial_genes);
  say "Select genes: ".scalar(@$select_genes);
  my $final_genes = build_final_geneset($select_genes);
  say "Final genes: ".scalar(@$final_genes);

  my $genes_to_write = prep_genes_for_writing($final_genes,$write_single_transcript_genes,$analysis);

  say "Number of prepped genes: ".scalar(@$genes_to_write);

  if($set_canonical) {
    foreach my $gene (@$genes_to_write) {
      set_canonical_transcript($gene);
    }
  }

  write_to_gtf_file($genes_to_write,$output_gtf_file,$analysis_name,$final_biotype);
}

exit;


sub write_to_gtf_file {
  my ($genes_to_write,$output_gtf_file,$analysis_name,$final_biotype) = @_;

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
      my $record = build_gtf_record($transcript,$exons,$analysis_name,$gene_id,$transcript_id,$translation,$final_biotype);
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
  my ($transcript,$exons,$analysis_name,$gene_id,$transcript_id,$translation,$final_biotype) = @_;

  my $record = [];
  my $strand = "+";
  if($transcript->strand() == -1) {
    $strand = "-";
  }

  my $transcript_attribs = 'gene_id "'.$gene_id.'"; transcript_id "'.$transcript_id.'"; biotype "'.$final_biotype.'";';
  if($translation) {
    # This is a simple encoding for the translation in the output file
    my $start_exon = $translation->start_Exon();
    my $end_exon = $translation->end_Exon();
    my $translation_coords = ' translation_coords "'.$start_exon->start().':'.$start_exon->end().':'.$translation->start().':'.$end_exon->start().':'.$end_exon->end().':'.$translation->end().'";';
    $transcript_attribs .= $translation_coords;
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


sub prep_genes_for_writing {
  my ($final_genes,$write_single_transcript_genes,$analysis) = @_;

  my $genes_to_write = [];
  foreach my $gene (@$final_genes) {
    # This array will store exactly what to write. By default we would write the gene as constructed to this point
    # However if the run is just to process a particular type of transcripts into a cleaned set, we may want to
    # pull those back out into single transcript genes for futher processing later. In that case the array will
    # be replace with single transcript genes

    my $transcripts = $gene->get_all_Transcripts();
    if($write_single_transcript_genes && scalar(@$transcripts) > 1) {
      foreach my $transcript (@$transcripts) {
        my $new_gene = Bio::EnsEMBL::Gene->new(
                           -START     => $transcript->start(),
                           -END       => $transcript->end(),
                           -STRAND    => $transcript->strand,
                           -SLICE     => $transcript->slice());
        $new_gene->add_Transcript($transcript);
        $new_gene->analysis($analysis);
        push(@$genes_to_write,$new_gene);
      }
    } else {
      push(@$genes_to_write,$gene);
    }
  }
  return($genes_to_write);
}


sub remove_overlapping_exons {
  my ($transcripts) = @_;
  my $filtered_transcripts = [];

  foreach my $transcript (@$transcripts) {
    my $exons = $transcript->get_all_Exons();
    my $flagged = 0;
    for(my $i=0; $i<scalar(@$exons)-1; $i++) {
      my $exon1 = ${$exons}[$i];
      my $exon2 = ${$exons}[$i+1];
      if(features_overlap($exon1,$exon2)) {
        $flagged = 1;
        last;
      }
    }

    unless($flagged) {
      push(@$filtered_transcripts,$transcript);
    } else {
      say "Flagged transcript with start/end ".$transcript->seq_region_start." ".$transcript->seq_region_end.", will remove";
    }
  }
  return($filtered_transcripts);
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
  } # End foreach my $transcript

  $gene->canonical_transcript($current_canonical);
}


sub join_transcripts {
  my ($transcripts_to_join) = @_;

  my $extension_length = 250;
  # The plan
  # First loop through all the transcripts and extend the terminal exon by the extension length
  # Then cluster on genomic overlap
  # Go through the clusters, foreach transcript figure out if there is an extension possible
  # An extension is possible if:
  # 1) There is an overlap between the two transcript on the terminal intron
  # 2) There is an overlap between two transcript on the terminal exons
  #    For this one it doesn't matter if one of the exons overlaps an intron in the order (or both)
  # Join the transcript to:
  # 1) The candidate with the most exons
  # 2) The candidate that will extend the transcript the most
  # Tricks:
  # 1) Put all the exons on the forward strand temporarily, that way only a simple sort is needed instead of strand logic
  # 2) Resort to handle 5' and 3' in the same way if possible

  adjust_transcript_ends($transcripts_to_join,$extension_length);

  # Now after extending make them into genes to cluster on genomic overlap
  my $initial_genes = create_single_transcript_genes($transcripts_to_join);
  my $joined_transcripts = process_clusters_for_joining($initial_genes);
  $joined_transcripts =  restore_strand($joined_transcripts);

  return($joined_transcripts);
}


sub process_clusters_for_joining {
  my ($genes) = @_;

  my $all_joined_transcripts = [];
  my $biotypes_hash = get_all_biotypes([@$genes]);
  my $biotypes_array = [keys(%$biotypes_hash)];

  my $types_hash;
  $types_hash->{genes} = [@$biotypes_array];

  say "Clustering genes from input_dbs...";
  my ($clusters, $unclustered) = cluster_Genes($genes,$types_hash);
  foreach my $cluster (@$clusters) {
    my $genes = $cluster->get_Genes();
    my $extracted_transcripts = extract_transcripts_from_genes($genes);
    my $joined_transcripts = calculate_joins($extracted_transcripts);
    foreach my $transcript (@$joined_transcripts) {
      push(@$all_joined_transcripts,$transcript);
    }
  } # End foreach my $cluster

  return($all_joined_transcripts);
}


sub place_on_forward_strand {
  my ($transcripts_to_change) = @_;
  my $forward_transcripts = [];
  foreach my $transcript (@$transcripts_to_change) {
    my $original_strand = $transcript->strand;
    $transcript->{'original_strand'} = 1;
    if($transcript->strand == 1) {
      push(@$forward_transcripts,$transcript);
      next;
    }

    my $forward_exons = [];
    my $exons = $transcript->get_all_Exons();
    foreach my $exon (@$exons) {
      my $forward_exon = Bio::EnsEMBL::Exon->new(
                           -START     => $exon->start(),
                           -END       => $exon->end(),
                           -STRAND    => 1,
                           -SLICE     => $exon->slice(),
                           -PHASE     => -1,
                           -END_PHASE => -1);
      push(@$forward_exons,$forward_exon);
    }

    my $sorted_forward_exons = [sort { $a->start <=> $b->start } @{$forward_exons}];
    my $forward_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $sorted_forward_exons);
    $forward_transcript->strand(1);
    $forward_transcript->stable_id($transcript->stable_id());
    $forward_transcript->biotype($good_biotype);
    $forward_transcript->slice($transcript->slice());
    $forward_transcript->analysis($analysis);
    $forward_transcript->{'original_strand'} = $original_strand;
    push(@$forward_transcripts,$forward_transcript);
  } # foreach my $transcript
  return($forward_transcripts);
}


sub restore_strand {
  my ($transcripts_to_change) = @_;
  my $reverse_transcripts = [];
  foreach my $transcript (@$transcripts_to_change) {
    my $original_strand = $transcript->{'original_strand'};
    if($original_strand == 1) {
      push(@$reverse_transcripts,$transcript);
      next;
    }

    my $reverse_exons = [];
    my $exons = $transcript->get_all_Exons();
    foreach my $exon (@$exons) {
      my $reverse_exon = Bio::EnsEMBL::Exon->new(
                           -START     => $exon->start(),
                           -END       => $exon->end(),
                           -STRAND    => -1,
                           -SLICE     => $exon->slice(),
                           -PHASE     => -1,
                           -END_PHASE => -1);
      push(@$reverse_exons,$reverse_exon);
    }

    my $sorted_reverse_exons =  [sort { $b->end <=> $a->end } @{$reverse_exons}];
    my $reverse_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $sorted_reverse_exons);
    $reverse_transcript->strand(-1);
    $reverse_transcript->stable_id($transcript->stable_id());
    $reverse_transcript->biotype($good_biotype);
    $reverse_transcript->slice($transcript->slice());
    $reverse_transcript->analysis($analysis);
    copy_transcript_attribs($transcript,$reverse_transcript);
    push(@$reverse_transcripts,$reverse_transcript);
  } # foreach my $transcript
  return($reverse_transcripts);
}


sub calculate_joins {
  my ($extracted_transcripts) = @_;

  my $joined_transcripts = [];
  # Put all the transcripts on the forward strand
  $extracted_transcripts = place_on_forward_strand($extracted_transcripts);


  # Assign ids and index the introns
#  my $transcript_index = 1;
  my $intron_index = {};
#  $transcript_index = index_transcripts($extracted_transcripts,$transcript_index,$intron_index);

  # If we've already seen an exon string, there's no point in processing it again
  my $exon_string_record = {};

  # Loop through all the transcripts and calculate the two best candidate matches for both the
  for(my $i=0; $i<scalar(@$extracted_transcripts); $i++) {
    # Foreach of these transcripts you want to make up to two new transcripts, then push them on the end of the array
    # This could be disasterous, but it might just work. To begin with I'm not going to clean out the old transcripts
    # but this will be assessed in testing. Also could consider making very light versions of the new transcripts,
    # i.e. constructing exons from scratch and not adding any slices. That could keep the footprint very low until some
    # cleaning of the cluster was carried out. For example you might want to process each cluster fully to a final set
    # including all the cleaning and collapsing, that way you could keep memory spikes in check
    # Anyway, regardless, foreach transcript:
    # 1: Find the join that will give the most length
    # 2: Find the join that will give the most introns
    # 3: Construct both of these. If they're identical remove one
    # 4: Pass any constructed transcripts onto the end of the array
    # This way there should

    my $longest_seq_length = 0;
    my $most_exons_count = 0;
    my $best_joined_transcripts = [];
    my $transcript_i = ${$extracted_transcripts}[$i];
    my $transcript_i_exons = $transcript_i->get_all_Exons();
    my $transcript_i_exon_string = generate_exon_string($transcript_i_exons);

    if($exon_string_record->{$transcript_i_exon_string}) {
      next;
    } else {
      $exon_string_record->{$transcript_i_exon_string} = 1;
    }

    my $single_exon = 0;
    if(scalar(@{$transcript_i_exons}) == 1) {
      $single_exon = 1;
    }

    for(my $j=0; $j<scalar(@$extracted_transcripts); $j++) {
      # No point in comparing to itself
      if($i == $j) {
        next;
      }

      my $transcript_j = ${$extracted_transcripts}[$j];

      # First check is there any overlap between i and j, if not then move on
      unless(features_overlap($transcript_i,$transcript_j)) {
        next;
      }

      my $joined_transcript;
      #      if($single_exon) {
#        $joined_transcript = single_exon_join($transcript_i,$transcript_j);
      #        if($joined_transcript && $joined_transcript->{'original_strand'} == -1) {
#        }
#      } else {
       $joined_transcript = multi_exon_join($transcript_i,$transcript_j);
#      }

      # Check if the joined transcrpit has the longest seq or most exons
      if($joined_transcript) {
        my $joined_transcript_exon_count = scalar(@{$joined_transcript->get_all_Exons()});
        my $joined_transcript_length = $joined_transcript->length();

        # If it's longer, then make the new longest
        if($joined_transcript_length > $longest_seq_length) {
          $longest_seq_length = $joined_transcript_length;
          ${$best_joined_transcripts}[0] = $joined_transcript;
        }

        # If the exon count is the same as the current best, check the length of it versus the current best
        if($joined_transcript_exon_count == $most_exons_count) {
          if($joined_transcript_length > ${$best_joined_transcripts}[1]->length()) {
            ${$best_joined_transcripts}[1] = $joined_transcript;
          }
        }

        # If the exon count is higher than the current best then make the current best
        if($joined_transcript_exon_count > $most_exons_count) {
          $most_exons_count = $joined_transcript_exon_count;
          ${$best_joined_transcripts}[1] = $joined_transcript;
        }
      } # End if($joined_transcript)
    } # End for(my $j=0;

    # You'll want to determine if the two best transcripts are the same or not
    # if they are just add one, otherwise add both to the end of the transcript array
    # Make sure to assign an id, might also need to index the introns, and recond strand of original
    if(${$best_joined_transcripts}[0]) {
      ${$best_joined_transcripts}[0]->{'original_strand'} = $transcript_i->{'original_strand'};
    }

    if(${$best_joined_transcripts}[1]) {
      ${$best_joined_transcripts}[1]->{'original_strand'} = $transcript_i->{'original_strand'};
    }

    # Check if the transcripts are being identical and only keep one if they are
    my $final_joined_transcripts = [];
    my $exon_string1;
    my $exon_string2;

    if(${$best_joined_transcripts}[0]) {
      $exon_string1 = generate_exon_string(${$best_joined_transcripts}[0]->get_all_Exons());
    }

    if(${$best_joined_transcripts}[1]) {
      $exon_string2 = generate_exon_string(${$best_joined_transcripts}[1]->get_all_Exons());
    }

    if($exon_string1 && $exon_string1 eq $exon_string2) {
       push(@$final_joined_transcripts,${$best_joined_transcripts}[0]);
     } else {
       if($exon_string1) {
        push(@$final_joined_transcripts,${$best_joined_transcripts}[0]);
      }
       if($exon_string2) {
        push(@$final_joined_transcripts,${$best_joined_transcripts}[1]);
      }
     } # End else

    # Now index them and push them
    if(scalar(@$final_joined_transcripts)) {
#      $transcript_index = index_transcripts($final_joined_transcripts,$transcript_index,$intron_index);
      foreach my $final_joined_transcript (@$final_joined_transcripts) {
        my $final_transcript_exons = $final_joined_transcript->get_all_Exons();
        my $final_transcript_exon_string = generate_exon_string($final_transcript_exons);
        unless($exon_string_record->{$final_transcript_exon_string}) {
          push(@$extracted_transcripts,$final_joined_transcript);
          push(@$joined_transcripts,$final_joined_transcript);
          $exon_string_record->{$final_transcript_exon_string} = 1;
        }
      }
    } # End if(scalar
  } # End for(my $i=0;

  return($joined_transcripts);
}


sub single_exon_join {
  my ($transcript1,$transcript2) = @_;

  my $exons1 = $transcript1->get_all_Exons();
  my $exons2 = $transcript2->get_all_Exons();
  unless(scalar(@$exons1) == 1) {
    throw("In the single exon join subroutine, but had a mutli-exon first transcript. First transcript should be single exon");
  }

  my $exon1 = clone_Exon(${$exons1}[0]);

  # A few scenarios here:
  # 1: exon1 is completely contained in an exon2
  #    -> nothing should be done
  # 2: An exon2 is completely contained in exon1 and no other overlap
  #    -> nothing should be done
  # 3: An exon2 partially overlaps exon1 across the boundary
  #    -> A join should be calculated and exon1 should get overwritten
  #    -> This should be iterative, e.g if there are different 5' and 3' exons overlapping
  #    -> then they should be progressively picked up
  #    -> all exons2 before the first overlap and after the last overlap should be added to the  new exon set
  my $five_prime_exon2s = [];
  my $three_prime_exon2s = [];

  # These two variable check if there has been a modificiation to either the 5' or 3' boundary
  my $modified_five_prime = 0;
  my $modified_three_prime = 0;
  for(my $i=0; $i<scalar(@$exons2); $i++) {
    my $exon2 = ${$exons2}[$i];
    if(features_overlap($exon1,$exon2)) {
      # Check the different conditions
      # Note that originally I had thought it was clever to change seq_region_start/end. But this actually
      # does not update the start/end in general. It's still better to generally use seq_region_start/end
      # but not for adjusting boundaries
      if($exon2->seq_region_start() <= $exon1->seq_region_start() && $exon2->seq_region_end() >= $exon1->seq_region_end()) {
        # This means the first transcript is completely contained in an exon of the second
        # Just return in this case
        return;
      }

      if($exon2->seq_region_start() > $exon1->seq_region_start() && $exon2->seq_region_end() < $exon1->seq_region_start()) {
        # This means the exon2 is completely in exon1 and can be ignored
        next;
      } elsif($exon2->seq_region_end() >= $exon1->seq_region_start() && $exon2->seq_region_start() < $exon1->seq_region_start()) {
        # This means there is an overlap on the 5' boundary and exon1 should be updated to match exon2
        $exon1->start($exon2->start());
        $modified_five_prime = 1;
      } elsif($exon2->seq_region_start() <= $exon1->seq_region_end() && $exon2->seq_region_end() > $exon1->seq_region_end()) {
        # This means there is an overlap on the 3' boundary and exon1 should be updated to match exon2
        $exon1->end($exon2->end());
        $modified_three_prime = 1;
      }
    } else {
      if($exon2->seq_region_end() < $exon1->seq_region_start()) {
        push(@{$five_prime_exon2s},$exon2);
      } else {
        push(@{$three_prime_exon2s},$exon2);
      }
    } # End else
  } # End for(my $i=0;

  unless($modified_five_prime || $modified_three_prime) {
    return;
  }

  my $merged_exons = [];
  if($modified_five_prime) {
    push(@$merged_exons,@{$five_prime_exon2s});
  }

  push(@$merged_exons,$exon1);

  if($modified_three_prime) {
    push(@$merged_exons,@{$three_prime_exon2s});
  }

  $merged_exons = [sort { $a->start <=> $b->start } @{$merged_exons}];

  # Since these exons might be modified in other merges, clone them before making the new transcript
  $merged_exons = clone_exon_array($merged_exons);

  say "Creating joined transcript from single exon transcript. Joined transcript has ".scalar(@$merged_exons)." exons";
  my $joined_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $merged_exons);
  $joined_transcript->stable_id($transcript1->stable_id());
  $joined_transcript->biotype($good_biotype);
  $joined_transcript->slice($transcript1->slice());
  $joined_transcript->analysis($analysis);
  copy_transcript_attribs($transcript1,$joined_transcript);

  return($joined_transcript);
}


sub clone_exon_array {
  my ($exons) = @_;
  my $cloned_exons = [];
  foreach my $exon (@$exons) {
    push(@$cloned_exons,clone_Exon($exon));
  }
  return($cloned_exons);
}


#sub index_transcripts {
#  my ($transcripts_to_index,$transcript_index,$intron_index) = @_;
#  foreach my $transcript (@$transcripts_to_index) {
#    $transcript->{'internal_id'} = $transcript_index;
#    $transcript_index++;
#    my $introns = $transcript->get_all_Introns();
#    for(my $i=0; $i<scalar(@$introns); $i++) {
      # The intron string is the main key, which point to the transcript id and that has a
      # value of the index of the intron, I've made it 1-based since sometimes we will just
      # want use the value in conditionals to check presence
#      my $intron = ${$introns}[$i];
#      my $intron_string = generate_intron_string([$intron]);
#      $intron_index->{$intron_string}->{$transcript->{'internal_id'}} = $i+1;
#    }
#  }
#  return($transcript_index);
#}


sub multi_exon_join {
  my ($transcript1,$transcript2) = @_;

  my $cloned_transcript1 = $transcript1;#clone_Transcript($transcript1);
  my $exons1 = $cloned_transcript1->get_all_Exons();
  my $exons2 = $transcript2->get_all_Exons();

  my $cloned_transcipt1_exon_string = generate_exon_string($exons1);
  my $transcript2_exon_string = generate_exon_string($exons2);

  if($cloned_transcipt1_exon_string eq $transcript2_exon_string) {
    return;
  }

  # Do this in two parts
  # First examine for matches on the terminal introns, if there are then:
  # 1: If t1 has more exons past the boundary, replace the end exon of t1 with the corresponding exon from t2
  #    and add extra exons past the boundary from t2
  # 2: If there are no more exons then pick the largest exon between t1 and t2
  # Consider just looking at the intron/exon boundary
  # Keep the transcript that adds the most exons and the one that adds the most length
  # If there are no intron matches then look for exon overlap with the terminal exon
  # 1:

  # Note this is probably not something that needs to be different for single/multi exon
  # As all comparisons are being done and judged, just look at it from the perspective
  # of what exon from exons2 is overlaps the 5'/3' boundary exon and is closest to or
  # crosses the boundary. When you figure out what exon that is for the 5' and 3' sides
  # you adjust the start/end to match
  # Just need the exons that are on or over the start/end to decide to join. Do not need
  # introns to make the decision

  my $modified_five_prime = 0;
  my $modified_three_prime = 0;

  my $five_prime_exon2s = [];
  my $three_prime_exon2s = [];

  # Work out the internal exons for transcript1
  my $internal_exon1s = [];
  for(my $i=1; $i<scalar(@$exons1)-1; $i++) {
    push(@$internal_exon1s,${$exons1}[$i]);
  }

  my $terminal_exon1s = [${$exons1}[0]];
  if(scalar(@$exons1) > 1) {
    push(@$terminal_exon1s,${$exons1}[$#$exons1]);
  }

  for(my $i=0; $i<scalar(@$terminal_exon1s); $i++) {
    my $exon1 = ${$terminal_exon1s}[$i];
    for(my $j=0; $j<scalar(@$exons2); $j++) {
      my $exon2 = ${$exons2}[$j];
      if(features_overlap($exon1,$exon2)) {
        # Check the different conditions
        if($exon2->seq_region_start() > $exon1->seq_region_start() && $exon2->seq_region_end() < $exon1->seq_region_start()) {
          # This means the exon2 is completely in exon1 and can be ignored
          next;
        } elsif($exon2->seq_region_end() >= $exon1->seq_region_start() && $exon2->seq_region_start() <= $exon1->seq_region_start() && $i == 0) {
          # This means there is an overlap on the 5' boundary and exon1 should be updated to match exon2
          $exon1->start($exon2->start());
          $modified_five_prime = 1;
        } elsif($exon2->start() <= $exon1->end() && $exon2->end() >= $exon1->end() && $i == $#$terminal_exon1s) {
         # This means there is an overlap on the 3' boundary and exon1 should be updated to match exon2
          $exon1->end($exon2->end());
          $modified_three_prime = 1;
        }
      } else {
        if($exon2->end() < $cloned_transcript1->start() && $i == 0) {
          push(@{$five_prime_exon2s},$exon2);
        } elsif($exon2->start > $cloned_transcript1->end() && $i == $#$terminal_exon1s) {
          push(@{$three_prime_exon2s},$exon2);
        }
      } # End else
    } # End for(my $j=0;
  } # End for(my $i=0;

  unless($modified_five_prime || $modified_three_prime) {
    return;
  }

  my $merged_exons = [];
  if($modified_five_prime && scalar(@{$five_prime_exon2s})) {
    push(@$merged_exons,@{$five_prime_exon2s});
  }

  push(@$merged_exons,${$terminal_exon1s}[0]);

  if(scalar(@$internal_exon1s)) {
    push(@$merged_exons,@{$internal_exon1s});
  }

  if(scalar(@$terminal_exon1s) > 1) {
    push(@$merged_exons,${$terminal_exon1s}[$#$terminal_exon1s]);
  }

  if($modified_three_prime && scalar(@{$three_prime_exon2s})) {
    push(@$merged_exons,@{$three_prime_exon2s});
  }

  $merged_exons = [sort { $a->start <=> $b->start } @{$merged_exons}];

  # Since these exons might be modified in other merges, clone them before making the new transcript
  $merged_exons = clone_exon_array($merged_exons);

  my $merged_exon_string = generate_exon_string($merged_exons);

  if($merged_exon_string eq $cloned_transcipt1_exon_string || $merged_exon_string eq $transcript2_exon_string){
    return;
  }


  say "Creating joined transcript from single exon transcript. Joined transcript has ".scalar(@$merged_exons)." exons";
  my $joined_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $merged_exons);
  $joined_transcript->stable_id($transcript1->stable_id());
  $joined_transcript->biotype($good_biotype);
  $joined_transcript->slice($transcript1->slice());
  $joined_transcript->analysis($analysis);
  copy_transcript_attribs($transcript1,$joined_transcript);

  return($joined_transcript);
}


sub adjust_transcript_ends {
  my ($transcripts_to_extend,$extension_length) = @_;

  # You should put a cloning step in here if it's better to keep the unmodified transcript
  # For now just modify them directly
  foreach my $transcript (@$transcripts_to_extend) {
    my $exons = $transcript->get_all_Exons();
    if($transcript->strand() == 1) {
      my $new_start = ${$exons}[0]->start() - $extension_length;
      my $new_end = ${$exons}[$#{$exons}]->end() + $extension_length;
      if($new_start < 1) {
        $new_start = 1;
      }

      if($new_end > $transcript->slice->length()) {
        $new_end = $transcript->slice->length();
      }
    } else {
      my $new_start = ${$exons}[$#{$exons}]->start() - $extension_length;
      my $new_end = ${$exons}[0]->end() + $extension_length;
      if($new_start < 1) {
        $new_start = 1;
      }

      if($new_end > $transcript->slice->length()) {
        $new_end = $transcript->slice->length();
      }
    }
  } # foreach my $transcript
}


sub clean_initial_transcripts {
  my ($transcripts) = @_;

  my $min_multi_exon_cds_length = 300;
  my $min_single_exon_cds_length = 450;
  my $min_intron_length = 75;

  # Removes initial transcripts that don't pass some basic criteria
  # For example if the cds does not have a start/stop, or is short, or non canonical introns
  # It also removes short introns and recalculates the cds to see if it makes a difference
  my $cleaned_transcripts = [];
  my $progress_index = 1;
  foreach my $transcript (@$transcripts) {
    if($progress_index % 10 == 0) {
      say "Processing transcript ".$progress_index." of ".scalar(@$transcripts);
    }
    $progress_index++;

    # We are not going to skip transcripts without a cds initially, as that might be fixed by later
    # removing dodgy introns and recomputing
    my $cds_seq = $transcript->translateable_seq();
    my $cds_length = 0;
    if($cds_seq) {
      $cds_length = length($cds_seq);
    }

    # First thing to check are the introns, as this might lead to a significant change to the cds
    my $introns = $transcript->get_all_Introns();
    if(scalar(@$introns) >= 1) {
      # Check the introns and if there are small introns and removing them increases the cds length
      # then add the modified transcript back to the pile. Note that this is only going to do one
      # pass on each transcript for speed, i.e. remove all dodgy introns and recalulcate, not do
      # any sort of iterative removal/recompute
      my $modified_transcript = check_introns($transcript,$introns,$min_intron_length);
      if($modified_transcript) {
        compute_translation($modified_transcript);
        my $modified_cds_seq = $modified_transcript->translateable_seq();
        if(length($modified_cds_seq) >= $cds_length) {
          push(@$transcripts,$modified_transcript);
        }
      }
    } # End if(scalar(@$introns) >= 1)

    # Now validate the cds and push to the clean set if it passes
    my $valid_cds = 0;
    if(scalar(@$introns) >= 1) {
      $valid_cds = validate_cds($cds_seq,$min_multi_exon_cds_length);
    } else {
      $valid_cds = validate_cds($cds_seq,$min_single_exon_cds_length);
    }

    if($valid_cds) {
      push(@$cleaned_transcripts,$transcript);
    } else {
      say "CDS not valid!!!";
    }
  } # foreach my $transcript (@$transcripts)

  return($cleaned_transcripts);
}


sub check_introns {
  my ($transcript,$introns,$min_intron_length) = @_;

  my $exons = $transcript->get_all_Exons();
  my $was_modified = 0;
  my $modified_exons = [];
  my $previous_exon = $$exons[0];
  push(@$modified_exons,clone_Exon($previous_exon));
  for(my $intron_index = 0; $intron_index < scalar(@$introns); $intron_index++) {
    my $intron = ${$introns}[$intron_index];
    my $next_exon = ${$exons}[$intron_index+1];

    if($intron->length() < $min_intron_length) {
      my $merged_exon = merge_exons($previous_exon,$next_exon);
      $was_modified = 1;
      ${$modified_exons}[$#{$modified_exons}] = $merged_exon;
      $previous_exon = $merged_exon;
    } else {
      push(@$modified_exons,clone_Exon($next_exon));
      $previous_exon = $next_exon;
    }
  }

  my $modified_transcript;
  if($was_modified) {
    if($$modified_exons[0]->strand == 1) {
      $modified_exons = [sort { $a->start <=> $b->start } @{$modified_exons}];
    } else {
      $modified_exons = [sort { $b->end <=> $a->end } @{$modified_exons}];
    }

    $modified_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $modified_exons);
    $modified_transcript->stable_id($transcript->stable_id());
    $modified_transcript->biotype($good_biotype);
    $modified_transcript->slice($transcript->slice());
    $modified_transcript->analysis($analysis);
    compute_translation($modified_transcript);

#    say "Org:";
#    print_exons($exons);
#    say "Mod:";
#    print_exons($modified_exons);
    return($modified_transcript);
  }
}


sub merge_exons {
  my ($exon1,$exon2) = @_;


  my @coords = ($exon1->seq_region_start(),$exon1->seq_region_end(),$exon2->seq_region_start(),$exon2->seq_region_end());
  my $start = min(@coords);
  my $end = max(@coords);
  #  if($exon1->strand == -1) {
#    $start = $exon2->seq_region_start();
#    $end = $exon1->seq_region_end();
#  }

  my $merged_exon = Bio::EnsEMBL::Exon->new(
                      -START     => $start,
                      -END       => $end,
                      -STRAND    => $exon1->strand(),
                      -SLICE     => $exon1->slice(),
                      -PHASE     => -1,
                      -END_PHASE => -1);

  return($merged_exon);
}


sub validate_cds {
  my ($cds_seq,$min_cds_length) = @_;
  # Checks the cds has a proper start/stop and passes the length cutoff
#  unless($cds_seq =~ /^ATG/ && $cds_seq =~ /(TAA|TAG|TGA)$/ && length($cds_seq) >= $min_cds_length) {
  unless(length($cds_seq) >= $min_cds_length) {
    return(0);
  }

  return(1);
}


sub create_single_transcript_genes {
  my ($transcripts) = @_;

  # Creates single transcript genes for use with things like genebuilder
  my $single_transcript_genes = [];
  foreach my $transcript (@$transcripts) {
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->stable_id($transcript->stable_id());
    $gene->biotype($transcript->biotype);
    $gene->slice($transcript->slice());
    $gene->analysis($analysis);
    $gene->add_Transcript($transcript);
    push(@$single_transcript_genes,$gene);
  }

  return($single_transcript_genes);
}


sub extract_transcripts_from_genes {
  my ($genes) = @_;
  my $extracted_transcripts = [];
  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      push(@$extracted_transcripts,$transcript);
    }
  }
  return($extracted_transcripts);
}


sub process_transcript {
  my ($transcript,$transcripts,$processed_transcripts,$processed_exon_strings) = @_;

  my $consider_cds_only = 0; # NOTE: This is disabled/enabled directly from here, when we want to actually add this it needs to be a real param
  if($consider_cds_only) {
    # This is for pulling out cds seqs without UTR (for example in single cell stuff where it's mostly tightly packed single exon genes)
    # So the transcripts coming out of this are all then without UTR
    my $cds_exons = $transcript->get_all_translateable_Exons();
    my $cds_exon_string = generate_exon_string($cds_exons);
    if($processed_exon_strings->{$cds_exon_string}) {
      return;
    }

    $processed_exon_strings->{$cds_exon_string} = 1;
    my $cds_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $cds_exons);
#    $cds_transcript->{'internal_id'} = $transcript->{'internal_id'};
    $cds_transcript->{'original_strand'} = $transcript->{'original_strand'};
    $cds_transcript->stable_id($transcript->stable_id());
    $cds_transcript->biotype($good_biotype);
    $cds_transcript->slice($transcript->slice());
    $cds_transcript->analysis($transcript->analysis());
    compute_translation($cds_transcript);
    push(@$processed_transcripts,$cds_transcript);
  } else {
    # This checks if the exon string has been seen and if not it's added to the processed_transcripts array and marked
    # If it has been seen already, then this just returns
    my $exons = $transcript->get_all_Exons();
    my $exon_string = generate_exon_string($exons);
    if($processed_exon_strings->{$exon_string}) {
      return;
    }
    $processed_exon_strings->{$exon_string} = 1;
  }

  my $utr_5p_features = $transcript->get_all_five_prime_UTRs();
  my $utr_3p_features = $transcript->get_all_three_prime_UTRs();

  my $utr_5p_exons = create_exons_from_utr($utr_5p_features);
  my $utr_3p_exons = create_exons_from_utr($utr_3p_features);

  my $min_5p_exon_count = 2;
  my $min_3p_exon_count = 2;
  my $min_5p_length = 500;
  my $min_3p_length = 500;
  my $cds_buffer_length = 100;
  my $utr_5p_length = cumulative_feature_length($utr_5p_exons,$cds_buffer_length);
  my $utr_3p_length = cumulative_feature_length($utr_3p_exons,$cds_buffer_length);

  # If these are exon strings that have not been seen before, we want to build a transcript out of them and then assuming it
  # passed out critera they will get added back to the transcripts array, which means they'll come back into this method later
  # and be further broke down if possible
  unless($processed_exon_strings->{generate_exon_string($utr_5p_exons)} || (scalar(@$utr_5p_exons) < $min_5p_exon_count) || ($utr_5p_length < $min_5p_length)) {
    generate_new_transcript($utr_5p_exons,$transcript,$transcripts);
  }

  unless($processed_exon_strings->{generate_exon_string($utr_3p_exons)} || (scalar(@$utr_3p_exons) < $min_3p_exon_count) || ($utr_3p_length < $min_3p_length)) {
    generate_new_transcript($utr_3p_exons,$transcript,$transcripts);
  }

  # Now for the original transcript, trim the 3' UTR as needed
  # It is very unusual to have real introns in the 3' UTR. Usually an intron in the 3' UTR is incorrect, because of the reconstructor
  # or because of thing like transposons
  # First cut off any additional introns (these will be searched for ORFs part of the recursive nature of this process)
  # Then trim the remaining UTR based on looking for the PAS signal
  if(scalar(@$utr_3p_exons)) {
    $transcript = remove_three_prime_exons($transcript,$utr_3p_features);
    $transcript = trim_3prime_utr($transcript);
  }

  if(scalar(@$utr_5p_exons)) {
    $transcript = remove_five_prime_exons($transcript);
  }

  push(@$processed_transcripts,$transcript);
}


sub cumulative_feature_length {
  my ($features,$buffer_length) = @_;

  # If a buffer has been specificed take that into account, for example it would be unusual for one another gene to end closer than 100bp of the next
  # so remove the buffer when considering the availble length for a potential new CDS to fit into
  unless($buffer_length) {
    $buffer_length = 0;
  }

  my $cumulative_length = 0;
  foreach my $feature (@$features) {
    $cumulative_length += $feature->length();
  }

   $cumulative_length -= $buffer_length;
   return($cumulative_length);
}


sub remove_three_prime_exons {
  my ($transcript,$utr_3p_features) = @_;

  # Look at exons in pairs
  my $exons = $transcript->get_all_Exons();
  my $kept_exons = [];
  my $skip = 1;

  my $closest_utr_exon;
  for(my $i=scalar(@$exons)-1; $i>=0; $i--) {
    my $exon = ${$exons}[$i];
    if($exon->is_coding($transcript) and $skip) {
      if($i<scalar(@$exons)-1) {
        $closest_utr_exon = ${$exons}[$i+1];
      }
      $skip = 0;
    }

    unless($skip) {
      push(@$kept_exons,$exon);
    }
  }

  # Treat the first pure UTR exon (if it exists) as a special case. If it's close to the end of the CDS and it's not abnormally short, then allow it
  my $max_first_intron_dist = 100;
  my $min_first_intron_length = 250;
  my $min_closest_exon_length = 100;
  my $max_utr_length = 1000;
  if($closest_utr_exon and ${$utr_3p_features}[0]->length() < $max_utr_length and $closest_utr_exon->length() >= $min_closest_exon_length) {
    my $intron_size = 0;
    if($closest_utr_exon->strand == 1) {
      $intron_size = $closest_utr_exon->seq_region_start() - ${$kept_exons}[0]->seq_region_end() + 1;
    } else {
      $intron_size = ${$kept_exons}[0]->seq_region_start() - $closest_utr_exon->seq_region_end() + 1;
    }

    if($intron_size >= $min_first_intron_length) {
      push(@$kept_exons,$closest_utr_exon);
    }
  }


  my @unsorted_exons = @{$kept_exons};
  my @sorted_exons = ();
  if (${$kept_exons}[0]->strand == 1) {
    @sorted_exons = sort {$a->start <=> $b->start} @unsorted_exons;
  } else {
    @sorted_exons = sort {$b->start <=> $a->start} @unsorted_exons;
  }

  $kept_exons = \@sorted_exons;

  my $last_exon = ${$kept_exons}[-1];
  my $coding_offset = 0;
  # First figure out if there's a coding offset to take into account
  if($last_exon->is_coding($transcript)) {
    my $cds_end_coord = $last_exon->coding_region_end($transcript);
    if($last_exon->strand() == 1) {
      $coding_offset = $cds_end_coord - $last_exon->seq_region_start();
    } else {
      $coding_offset = $last_exon->seq_region_end() - $cds_end_coord;
    }
  }

  # If the length of the UTR of the final exon is > the max length then adjust
  my $last_exon_utr_length = $last_exon->length - $coding_offset;
  if($last_exon_utr_length > $max_utr_length) {
    my $diff = $last_exon_utr_length - $max_utr_length;
    if($last_exon->strand() == 1 and $diff > 0) {
      $last_exon->end($last_exon->end()-$diff);
    } else {
      $last_exon->start($last_exon->start()+$diff);
    }
  }


  my $new_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $kept_exons);
  $new_transcript->stable_id($transcript->stable_id);
  $new_transcript->biotype($transcript->biotype());
  $new_transcript->slice($transcript->slice());
  $new_transcript->analysis($transcript->analysis());
  copy_transcript_attribs($transcript,$new_transcript);

  return($new_transcript);
}


sub remove_five_prime_exons {
  my ($transcript) = @_;

  unless($transcript->translation()) {
    compute_translation($transcript);
  }

  # Want to go through the 5p UTR exons
  # Drop most 5' exons if UTR only and length < min_terminal_length
  # For remaining exons, go through and work at what exon the UTR is already greater than the allowed length
  # Stop at that exon and cut back to max allowed length

  # Look at exons in pairs
  my $exons = $transcript->get_all_Exons();
  my $kept_exons = [];
  my $skip = 1;
  my $min_utr_terminal_exon_length = 100;
  for(my $i=0; $i<scalar(@$exons); $i++) {
    my $exon = ${$exons}[$i];
    if($exon->is_coding($transcript)) {
      $skip = 0;
    }

    # If we're hit a coding exon we'll push it as skip will have switched to 0 at that point
    # and everything else will be pushed
    # If we're still in the 5' UTR exons push the exon if it meats the min size criteria
    unless($skip) {
      push(@$kept_exons,$exon);
    } elsif($skip and $exon->length > $min_utr_terminal_exon_length) {
      push(@$kept_exons,$exon);
    }
  }

  my @unsorted_exons = @{$kept_exons};
  my @sorted_exons = ();
  if (${$kept_exons}[0]->strand == 1) {
    @sorted_exons = sort {$a->start <=> $b->start} @unsorted_exons;
  } else {
    @sorted_exons = sort {$b->start <=> $a->start} @unsorted_exons;
  }

  $kept_exons = \@sorted_exons;

  my $max_5p_utr_length = 300;
  # Now loop through the sorted exons and find the index if the first coding exon, then work back in the 5' direction
  my $cds_exon_index = 0;
  foreach my $exon (@$kept_exons) {
    if($exon->is_coding($transcript)) {
      last;
    }
    $cds_exon_index++;
  }

  my $cumulative_exon_length = 0;
  my $exons_5p_to_keep = [];
  for(my $i=$cds_exon_index; $i>=0; $i--) {
    my $exon = ${$kept_exons}[$i];
    my $coding_offset = 0;
    if($exon->is_coding($transcript)) {
      my $cds_start_coord = $exon->coding_region_start($transcript);
      if($exon->strand() == 1) {
        $coding_offset = $exon->seq_region_end - $cds_start_coord;
      } else {
        $coding_offset = $cds_start_coord - $exon->seq_region_start;
      }
    }

    my $exon_utr_length = $exon->length() - $coding_offset;
    my $min_remaining_length = 50;
    # At this point we know the length of the current UTR and the cumulative length
    # If the length of the of the current utr exon added to the cumulative total > max_5p_utr_length (minus any coding offset)
    # then this exon needs to either be dropped or it's the last one (technically the first one). First determine if once trimmed
    # the length is >= min_remaining_length, if it is then trim and keep. If it's not but it's a coding exon, then you don't need
    # do anything and just keep as is. Otherwise just drop it. Once any of these conditions are hit, stop and add all the exons
    # after the coding exon index
    my $current_cumulative_length = $cumulative_exon_length + $exon_utr_length;
    if($current_cumulative_length > $max_5p_utr_length) {
      my $overlimit_length = $current_cumulative_length - $max_5p_utr_length;
      my $trimmed_length = $exon_utr_length - $overlimit_length;
      if($trimmed_length < $min_remaining_length and $exon->is_coding($transcript)) {
        push(@$exons_5p_to_keep,$exon);
        last;
      } elsif($trimmed_length < $min_remaining_length) {
        last;
      } else {
        # In this case we want to trim the exon
        if($exon->strand() == 1) {
          $exon->start($exon->start()+$overlimit_length);
        } else {
          $exon->end($exon->end()-$overlimit_length);
        }
        push(@$exons_5p_to_keep,$exon);
        last;
      }
    }

    push(@$exons_5p_to_keep,$exon);
    $cumulative_exon_length = $current_cumulative_length;
  }


  my $final_exons = $exons_5p_to_keep;
  for(my $i=$cds_exon_index+1; $i<scalar(@$kept_exons); $i++) {
    my $exon = ${$kept_exons}[$i];
    push(@$final_exons,$exon);
  }

#  unless(scalar()) {

#  }
  my $new_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $final_exons);
  $new_transcript->stable_id($transcript->stable_id);
  $new_transcript->biotype($transcript->biotype());
  $new_transcript->slice($transcript->slice());
  $new_transcript->analysis($transcript->analysis());
  copy_transcript_attribs($transcript,$new_transcript);

  return($new_transcript);
}



sub create_exons_from_utr {
  my ($utr_features) = @_;

  my $exons = [];
  foreach my $utr (@$utr_features) {
    my $exon = Bio::EnsEMBL::Exon->new(
                                        -START     => $utr->start(),
                                        -END       => $utr->end(),
                                        -STRAND    => $utr->strand(),
                                        -SLICE     => $utr->slice(),
                                        -PHASE     => -1,
                                        -END_PHASE => -1);
    push(@$exons,$exon);
  }

  return($exons);
}


sub copy_transcript_attribs {
  my ($transcript1,$transcript2) = @_;

  # Copy attribs from transcript1 to transcript2
#  $transcript2->{'internal_id'} = $transcript1->{'internal_id'};
  $transcript2->{'original_strand'} = $transcript1->{'original_strand'};
}


sub generate_new_transcript {
  my ($utr_exons,$transcript,$transcripts) = @_;

  my $min_multi_exon_cds_length = 300;
  my $min_single_exon_cds_length = 450;

  my $utr_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $utr_exons);
  $utr_transcript->stable_id($transcript->stable_id);
  $utr_transcript->biotype($good_biotype);
  $utr_transcript->slice($transcript->slice());
  $utr_transcript->analysis($transcript->analysis());
  copy_transcript_attribs($transcript,$utr_transcript);

  # If it's not at least a little bigger than the min cds length (i.e. to have some UTR), return
  # To be honest this doesn't make too much sense, 50bp is a bit random
  unless($utr_transcript->length() >= ($min_single_exon_cds_length + 50)) {
    return;
  }

  compute_translation($utr_transcript);
  my $cds_seq = $utr_transcript->translateable_seq();

  my $valid_cds = 0;
  if(scalar(@$utr_exons) == 1) {
    $valid_cds = validate_cds($cds_seq,$min_single_exon_cds_length);
  } elsif(scalar(@$utr_exons) >= 1) {
    $valid_cds = validate_cds($cds_seq,$min_multi_exon_cds_length);
  }

  unless($valid_cds) {
    return;
  } else {
    push(@$transcripts,$utr_transcript);
  }
}


sub trim_3prime_utr {
  my ($transcript) = @_;

  # Only going to use the most common nuclear PAS signal
  my $pas_signal = 'AATAAA';
  my $cleavage_signal = 'CA';
  my $max_no_cleavage = 1000;

  unless($transcript->three_prime_utr) {
    warning("The trim_3prime_utr_short_read was called on a transcript with no 3 prime UTR. Nothing to trim");
    return($transcript);
  }

  my $exons = $transcript->get_all_Exons;
  my $final_exon = ${$exons}[$#{$exons}];
  my $translation_end_exon = $transcript->translation->end_Exon;

  my $coding_offset = 0;
  my $final_exon_seq = $final_exon->seq->seq;
  if($final_exon->start == $translation_end_exon->start) {
    $coding_offset = $transcript->translation->end;
  }

  my $found_pas = 0;
  my $pas_start = 0;
  my $pas_end = 0;
  my $cleavage_site = 0;
  while($final_exon_seq =~ /$pas_signal/g && !$found_pas) {
    # Set to 1 base offset for ease
    $pas_start = $-[0] + 1;
    $pas_end =  $+[0];

    if($pas_start <= $coding_offset) {
      next;
    }

    say "Found PAS signal in final exon seq at the following coords: ".$pas_start."..".$pas_end;
    $found_pas = 1;
    last;
  }

  if($found_pas) {
    # There are 15-30bp between the end of the pas signal and the cleavage site
    # as pas_end is already shifted to 1bp offset, just add 14
    my $post_pas_seq = substr($final_exon_seq,$pas_end + 14,15);
    say $post_pas_seq;
    if($post_pas_seq =~ /CA/) {
      $cleavage_site = $pas_end + 14 + $+[0];
      say "Cleavage site found within 15-30bp range of PAS signal";
      say $cleavage_site;
    } else {
      $cleavage_site = $pas_end + 30;
      say "Cleavage site not found within 15-30bp range of PAS signal, setting to 30bp downstream:";
      say $cleavage_site;
    }
  } else {
    if((length($final_exon_seq) - $coding_offset) > $max_no_cleavage) {
      $cleavage_site = $coding_offset + $max_no_cleavage;
    }
    say "Could not find PAS signal in 3' UTR, will use max_no_cleavage as a cut-off:";
    say $cleavage_site;
  }

  if($cleavage_site >= length($final_exon_seq)) {
    say "Not cleaving as proposed cleavage site is at or over the end of the final exon";
    return($transcript);
  }

  if($cleavage_site) {
    if($final_exon->strand == 1 && $cleavage_site > $coding_offset) {
      $final_exon->end($final_exon->start + $cleavage_site - 1);
      $transcript->end($final_exon->end);
    } elsif($final_exon->strand == -1 && $cleavage_site > $coding_offset) {
      $final_exon->start($final_exon->end - $cleavage_site + 1);
      $transcript->start($final_exon->start);
    }
  }

  return($transcript);
}


sub collapse_transcripts {
  my ($transcripts) = @_;

  # This remove collapse transcripts that are contained in other transcripts
  # If a transcript has the same intron pattern as another transcript then they are merged
  # keeping the longest terminal exon on each end
}


sub build_final_geneset {
  my ($genes) = @_;

  my $final_genes = [];
  my $good_genes_by_slice = sort_features_by_slice($genes);

  foreach my $slice_name (keys(%$good_genes_by_slice)) {
    my $genes = $good_genes_by_slice->{$slice_name};
    my $slice = ${$genes}[0]->slice();
    my $biotype = ${$genes}[0]->biotype();
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
                     -query => $slice,
                     -analysis => $analysis,
                     -genes => $genes,
                     -output_biotype => $biotype,
                     -max_transcripts_per_cluster => 100,
                     -min_short_intron_len => 7,
                     -max_short_intron_len => 15,
                     -blessed_biotypes => {},
                     -skip_readthrough_check => 1,
                     -max_exon_length => 50000,
                     -coding_only => 1,
                   );
    $runnable->run;
#    say "Created ".scalar(@{$runnable->output})." output genes";
    push(@$final_genes,@{$runnable->output});
  }

  return($final_genes);
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


sub select_geneset {
  my ($genes) = @_;

  my $min_score = 1;
  my $stats = {
    'cds_exon_length' => 150,
    'genomic_span' => 20000,
    'cds_exons' => 5,
    'cds_length' => 300,
    'cds_ratio' => 0.85,
    'small_orf_length' => 300,
    'small_genomic_span' => 1000,
    'single_exon_cds_length' => 450,
  };

  my @cds_ratio = ();
  my @cds_length = ();
  my @cds_exon_count = ();
  my @cds_exon_lengths = ();
  my @exon_count = ();
  my @intron_sizes = ();
  my @genomic_span = ();
  my @scores = ();

  my $final_genes = [];
  my $initial_good_genes = [];
  my $initial_bad_genes = [];
  my $revised_good_genes = [];
  my $exon_strings = {};
  foreach my $gene (@$genes) {
    my $transcript = ${$gene->get_all_Transcripts}[0];
    my $cds_exons =  $transcript->get_all_CDS();
    my $cds_exon_string = generate_exon_string($cds_exons);
    if($exon_strings->{$cds_exon_string}) {
      $exon_strings->{$cds_exon_string}++;
    } else {
      $exon_strings->{$cds_exon_string} = 1;
    }
  }


  foreach my $gene (@$genes) {
    my $biotype = $gene->biotype;

    my $transcript = ${$gene->get_all_Transcripts()}[0];

    my $exons = $transcript->get_all_Exons();
    my $cds_exons = $transcript->get_all_CDS();
    my $cds_exon_string = generate_exon_string($cds_exons);
    my $cds_exon_length = 0;

    unless(scalar(@$cds_exons)) {
      push(@$initial_bad_genes,$gene);
      next;
    }

    foreach my $cds_exon (@$cds_exons) {
      $cds_exon_length += $cds_exon->length;
    }

    $cds_exon_length = $cds_exon_length / scalar(@$cds_exons);
    my $cds_ratio = scalar(@$cds_exons) / scalar(@$exons);
    my $cds_seq = $transcript->translateable_seq;
    my $cds_length = length($cds_seq);
    my $genomic_span = $transcript->seq_region_end - $transcript->seq_region_start + 1;
    my $introns = $transcript->get_all_Introns();
    foreach my $intron (@$introns) {
      push(@intron_sizes,$intron->length);
    }


    my $score = 0;
    if(scalar(@{$cds_exons}) == 1 && $cds_length >= $stats->{'single_exon_cds_length'}) {
      $score++;
    } elsif($cds_length >= $stats->{'cds_length'}) {
      $score++;
    }

# The higher euk script had a lot of stats, but not really applicable to something like plasmodium
# Consider adding things back in

    #    if($cds_length >= $stats->{'cds_length'}) {
#      $score++;
#    }

    #    if($exon_strings->{$cds_exon_string}) {
#      $score += 2 * log($exon_strings->{$cds_exon_string});
#    }

    #    if($cds_exon_length >= $stats->{'cds_exon_length'}) {
#      $score++;
#    }

    #    if($cds_ratio >= $stats->{'cds_ratio'}) {
#      $score++;
#    }

    #    if(scalar(@$cds_exons) >= $stats->{'cds_exons'}) {
#      $score++;
#      my $diff = scalar(@$cds_exons) - $stats->{'cds_exons'};
    #      if($diff) {
#        $score += log($diff);
#      }
#    }

    #    if($cds_seq =~ /^ATG/) {
#      $score++;
#    }

    #    if($cds_seq =~ /(TAA|TAG|TGA)$/) {
#      $score++;
#    }

    push(@scores,$score);

#    say "Score: ".$score;
    if($score < $min_score) {
      push(@$initial_bad_genes,$gene);
    } else {
      push(@$initial_good_genes,$gene);
    }
  }

  $revised_good_genes = process_initial_good_genes($initial_good_genes,$initial_bad_genes);
#  $revised_good_genes = flag_small_models($revised_good_genes,$initial_bad_genes,$stats);

  foreach my $gene (@$revised_good_genes) {
    my $transcript = ${$gene->get_all_Transcripts()}[0];
    $transcript->biotype($good_biotype);
    $gene->biotype($good_biotype);
    push(@$final_genes,$gene);
    #$gene_adaptor->store($gene);
  }


  my $final_bad_genes = final_classification($revised_good_genes,$initial_bad_genes);
  foreach my $gene (@$final_bad_genes) {
    my $transcript = ${$gene->get_all_Transcripts()}[0];
    $transcript->biotype($bad_biotype);
    $gene->biotype($bad_biotype);
    push(@$final_genes,$gene);
    #$gene_adaptor->store($gene);
  }

  return($final_genes);
}


sub final_classification {
  my ($good_genes,$bad_genes) = @_;

  my $final_bad_genes = [];

  my $bad_biotype = 'transcriptomic_check';
  foreach my $gene (@$bad_genes) {
    $gene->biotype($bad_biotype);
  }

  my $good_biotypes_hash = get_all_biotypes([@$good_genes]);
  my $good_biotypes_array = [keys(%$good_biotypes_hash)];
  my $bad_biotypes_array = [$bad_biotype];

  my $types_hash;
  $types_hash->{genes} = [@$good_biotypes_array,$bad_biotype];

  my $all_genes = [@$good_genes,@$bad_genes];
  say "Clustering genes from input_dbs...";
  my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap($all_genes,$types_hash);
  foreach my $cluster (@$clusters) {
    my $good_cluster_genes = $cluster->get_Genes_by_Type($good_biotypes_array);
    my $bad_cluster_genes = $cluster->get_Genes_by_Type($bad_biotypes_array);

    unless(scalar(@$bad_cluster_genes)) {
       next;
     }

    if (!scalar(@$good_cluster_genes) and scalar(@$bad_cluster_genes)) {
      push(@$final_bad_genes,@$bad_cluster_genes);
      next;
    }

    # We have good and bad, so update the bad to ingore
    foreach my $gene (@$bad_cluster_genes) {
      $gene->biotype('transcriptomic_ignore');
      push(@$final_bad_genes,$gene);
    }
  } # End foreach my $cluster

  foreach my $singleton (@$unclustered) {
    # Note will loop even if this should just be single genes
    my $bad_single_genes = $singleton->get_Genes_by_Type($bad_biotypes_array);
    foreach my $gene (@$bad_single_genes) {
      push(@$final_bad_genes,$gene);
    }
  }
  return($final_bad_genes);
}


sub process_initial_good_genes {
  my ($good_genes,$bad_genes) = @_;

  my $revised_good_genes = [];
  my $biotypes_hash = get_all_biotypes([@$good_genes,@$bad_genes]);
  my $biotypes_array = [keys(%$biotypes_hash)];

  my $types_hash;
  $types_hash->{genes} = $biotypes_array;

  say "Clustering genes from input_dbs...";
  my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap($good_genes,$types_hash);

  foreach my $cluster (@$clusters) {
    my $intron_strings = {};
    my $cluster_genes = $cluster->get_Genes();
    foreach my $gene (@$cluster_genes) {
      my $transcript = ${$gene->get_all_Transcripts()}[0];
      my $cds_exons = $transcript->get_all_CDS();
      # This is faster than calling to get all CDS introns from the API
      my $cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
      copy_transcript_attribs($transcript,$cds_transcript);
      my $introns = $cds_transcript->get_all_Introns();
      my $intron_string = generate_intron_string($introns);
      $intron_strings->{$intron_string} = 1;
    }

    # This could obviously be sped up, but the cluster should be small
    my $cluster_good_genes = [];
    my $max_cluster_cds_exons = 0;
    my $max_cluster_orf_length = 0;
    my @unique_strings = keys(%$intron_strings);
    foreach my $gene (@$cluster_genes) {
      my $transcript = ${$gene->get_all_Transcripts()}[0];
      my $cds_exons = $transcript->get_all_CDS();
      # This is faster than calling to get all CDS introns from the API
      my $cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
      copy_transcript_attribs($transcript,$cds_transcript);
      my $introns = $cds_transcript->get_all_Introns();
      my $intron_string = generate_intron_string($introns);
      my $is_substring = 0;
      my $cds_exon_count = scalar(@{$transcript->get_all_CDS()});
      my $cds_length = length($transcript->translateable_seq);
      foreach my $unique_string (@unique_strings) {
        if($unique_string eq $intron_string) {
          next;
        }

        if($unique_string =~ /$intron_string/) {
          $is_substring = 1;
          last;
        }
      }

      if($is_substring) {
        push(@$bad_genes,$gene);
      } else {
        push(@$cluster_good_genes,$gene);
        if($cds_exon_count > $max_cluster_cds_exons) {
          $max_cluster_cds_exons = $cds_exon_count;
        }

        if($cds_length > $max_cluster_orf_length) {
          $max_cluster_orf_length = $cds_length;
        }
      }
    } # End foreach my $gene (@$cluster_genes)

    foreach my $gene (@$cluster_good_genes) {
      my $transcript = ${$gene->get_all_Transcripts()}[0];
      my $cds_exon_count = scalar(@{$transcript->get_all_CDS()});
      my $cds_length = length($transcript->translateable_seq);
      if($cds_exon_count < ($max_cluster_cds_exons * 0.75) && $cds_length < ($max_cluster_orf_length * 0.75)) {
        push(@$bad_genes,$gene);
      } else {
        push(@$revised_good_genes,$gene);
      }
    }
  } # End foreach my $cluster (@$clusters)

  foreach my $singleton (@$unclustered) {
    # Note will loop even if this should just be single genes
    my $single_genes = $singleton->get_Genes();
    foreach my $gene (@$single_genes) {
      push(@$revised_good_genes,$gene);
    }
  }

  return($revised_good_genes);
}


sub flag_small_models {
  my ($good_genes,$bad_genes,$stats) = @_;

  my $revised_good_genes = [];
  my $small_orf_length = $stats->{'small_orf_length'};
  my $small_genomic_span = $stats->{'small_genomic_span'};
  my $two_exon_scaling = 2;

  foreach my $gene (@$good_genes) {
    my $transcript = ${$gene->get_all_Transcripts()}[0];
    my $cds_seq = $transcript->translateable_seq;

    unless($cds_seq) {
      push(@$bad_genes,$gene);
      next;
    }

    my $cds_length = length($cds_seq);
    my $cds_exon_count = scalar(@{$transcript->get_all_CDS});
    my $genomic_span = $transcript->seq_region_end - $transcript->seq_region_start + 1;

    if($cds_exon_count == 1) {
      push(@$bad_genes,$gene);
      next;
    }

    if($cds_exon_count == 2) {
      unless($cds_seq =~ /^ATG/ && $cds_seq =~ /(TAA|TAG|TGA)$/) {
        push(@$bad_genes,$gene);
        next;
      }

      unless(length($cds_seq) >= $small_orf_length * $two_exon_scaling) {
        push(@$bad_genes,$gene);
        next;
      }
    }

    if($genomic_span < $small_genomic_span || $cds_length < $small_orf_length) {
      push(@$bad_genes,$gene);
      next;
    } else {
      push(@$revised_good_genes,$gene);
    }
  }
  return($revised_good_genes);
}


sub sort_transcripts_by_slice {
  my ($transcripts) = @_;

  my $sorted_transcript_hash = {};
  foreach my $transcript (@$transcripts) {
    my $seq_region_name = $transcript->seq_region_name();
    unless($sorted_transcript_hash->{$seq_region_name}) {
      $sorted_transcript_hash->{$seq_region_name} = [];
    }
    push(@{$sorted_transcript_hash->{$seq_region_name}},$transcript);
  }
  return($sorted_transcript_hash);
}


sub generate_intron_string {
  # Assumes stranded data
  my ($intron_array) = @_;

  my $intron_string = "";
  my $count = 0;
  foreach my $intron (@{$intron_array}) {
    my $start = $intron->start();
    my $end = $intron->end();
    $intron_string .= $start."..".$end.":";
  }
  return($intron_string);
}


sub get_all_biotypes {
  my ($master_genes_array) = @_;

  my $master_biotypes_hash = {};

  foreach my $gene (@{$master_genes_array}) {
    unless($master_biotypes_hash->{$gene->biotype}) {
      $master_biotypes_hash->{$gene->biotype} = 1;
    }
  }
  return($master_biotypes_hash);
}


sub extract_genes {
  my ($clusters,$type) = @_;
  my $genes = [];
  foreach my $cluster (@$clusters) {
    if($type) {
      push(@$genes,@{$cluster->get_Genes_of_Type($type)});
    } else {
      push(@$genes,@{$cluster->get_Genes});
    }
  }
  return($genes);
}


sub generate_exon_string {
  my ($exons) = @_;

  my $exon_string = "";
  foreach my $exon (@$exons) {
    $exon_string .= $exon->start.":".$exon->end.":".$exon->strand;
  }

  return($exon_string);
}


sub log_array {
  my (@a) = @_;

  for(my $i=0; $i < scalar(@a); $i++) {
    $a[$i] = log($a[$i]);
  }

  return(@a);
}


sub print_exons {
  my ($exons) = @_;

  my $exon_string = "";
  foreach my $exon (@$exons) {
    $exon_string .= "(".$exon->seq_region_start."..".$exon->seq_region_end."):";
  }
  say $exon_string.$$exons[0]->strand();
}
