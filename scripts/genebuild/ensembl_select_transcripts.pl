#!/usr/bin/env perl

# NB:
# We had a method to take the transcript as canonical if it contains all coding exons in the gene
# however, this is not a good idea as transcripts with transposable elements that may actually prevent
# them from being transcribed or translated, or are completely on existent, will automatically 
# get chosen over what could be much better candidates. It is available in previous git commits

use strict;
use warnings;
use feature 'say';
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(empty_Transcript);

my $ensembl_version = 94;
my $outdir = '/path/to/output_dir/';

# If you are making an output db
my $clone_db_script = '/path/to/clone_database.ksh';


my $dbname = 'homo_sapiens_core_'.$ensembl_version.'_38';
my $dbuser = 'WRITE_USER';
my $dbpass = 'WRITE_PASS';
my $dbhost = '';
my $dbport = '';

my $otherfdbname = 'homo_sapiens_otherfeatures_'.$ensembl_version.'_38';
my $otherfdbuser = 'READ_USER';
my $otherfdbhost = '';
my $otherfdbport = '';

my $vardbname = 'homo_sapiens_variation_'.$ensembl_version.'_38';
my $vardbuser = 'READ_USER';
my $vardbhost = '';
my $vardbport = '';

# The following options need to be specified - 0 for no 1 for yes
my $download_data = 1;
my $reward_refseq_match = 1;

# The weighting of the sources depends on whether or not you want to give extra reward to transcripts with a RefSeq match
my $appris_weight;
my $tsl_weight;
my $uniprot_weight;
my $refseq_canonical_weight;
my $refseq_match_weight;

if ($reward_refseq_match) {
  $appris_weight = 13;
  $tsl_weight = 12;
  $uniprot_weight = 2;
  $refseq_canonical_weight = 9;
  $refseq_match_weight = 4;
} else {
  $appris_weight = 11;
  $tsl_weight = 10;
  $uniprot_weight = 2;
  $refseq_canonical_weight = 4;
  $refseq_match_weight = 0;
}

my $length_fraction = 0.75;

my $log_file = $outdir.'/transcript_selection.log';
my $gene_file = $outdir.'/all_genes_transcripts.txt';
my $canonicals_file = $outdir.'/canonicals.txt';
my $alert_file1 = $outdir.'/multiple_transcripts_exons_'.$ensembl_version.'.txt';
my $alert_file2 = $outdir.'/multiple_transcripts_variants_'.$ensembl_version.'.txt';


# The following are biotypes that we'll consider for the canonical transcripts. This list can be updated as needed.
my @allowed_biotypes = ('protein_coding','TR_V_gene','TR_J_gene','TR_D_gene','TR_C_gene','IG_D_gene','IG_J_gene','IG_V_gene','IG_C_gene');

GetOptions(
  "dbuser|user|u=s" => \$dbuser,          # string
  "dbhost|host|h=s" => \$dbhost,          # string
  "dbport|port|P=i" => \$dbport,          # numeric
  "dbpass|pass|p=s" => \$dbpass,          # string
  "dbname|db|D=s" => \$dbname,            # string
  "otherfuser=s" => \$otherfdbuser,       # string
  "otherfhost=s" => \$otherfdbhost,       # string
  "otherfport=i" => \$otherfdbport,       # numeric
  "otherfdbname=s" => \$otherfdbname,     # string
  "varuser=s" => \$vardbuser,             # string
  "varhost=s" => \$vardbhost,             # string
  "varport=i" => \$vardbport,             # numeric
  "vardbname=s" => \$vardbname,           # string
);


my $coredb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $dbhost,
  -port => $dbport,
  -user => $dbuser,
  -dbname => $dbname,
  -pass => $dbpass,
);

my $otherfdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $otherfdbhost,
  -port => $otherfdbport,
  -user => $otherfdbuser,
  -dbname => $otherfdbname,
);

my $vardb = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
  -host => $vardbhost,
  -port => $vardbport,
  -user => $vardbuser,
  -dbname => $vardbname,
);

my $vfa = $vardb->get_VariationFeatureAdaptor();
say "DUMPER: ".Dumper($vfa);
my $sa = $coredb->get_SliceAdaptor();
my $ta = $coredb->get_TranscriptAdaptor();
my $ga = $coredb->get_GeneAdaptor();
my $aa = $coredb->get_AttributeAdaptor();

# make the output directory
system ("mkdir -p $outdir");

# get the gene names
say "Retrieving gene names";
my %xref;
my $xref_select = $coredb->dbc->prepare("SELECT DISTINCT g.stable_id, display_label FROM gene g LEFT JOIN object_xref ox ON g.gene_id = ox.ensembl_id JOIN xref x \
                                         ON x.xref_id = ox.xref_id LEFT JOIN transcript t ON t.gene_id = ox.ensembl_id LEFT JOIN translation tn ON t.transcript_id = tn.transcript_id \
                                         WHERE t.source IN ('ensembl','havana','ensembl_havana') AND tn.stable_id IS NOT NULL AND external_db_id = 1100");
$xref_select->execute();
while (my $xref_row = $xref_select->fetchrow_arrayref()){
  my ($gene_id, $label) = @$xref_row;
  $xref{$gene_id} = $label;
}
$xref_select->finish;

# get the appris P1s, P2s and P3s from the core db
say "Retrieving APPRIS data";
my %appris;
my $appris_select = $coredb->dbc->prepare("SELECT DISTINCT t.stable_id, ta.value FROM gene g JOIN transcript t ON g.gene_id = t.gene_id LEFT JOIN transcript_attrib \
                                           ta ON t.transcript_id = ta.transcript_id LEFT JOIN translation tn ON ta.transcript_id = tn.transcript_id WHERE ta.value IN \
                                           ('principal1','principal2','principal3') AND ta.attrib_type_id = 427 AND tn.translation_id IS NOT NULL");
$appris_select->execute();

while (my $appris_row = $appris_select->fetchrow_arrayref()){
  my ($appris_transcript, $value) = @$appris_row;
  if ($value eq 'principal1') {
    $appris{$appris_transcript} = 1;
  } elsif ($value eq 'principal2') {
    $appris{$appris_transcript} = 2;
  } elsif ($value eq 'principal3') {
    $appris{$appris_transcript} = 3;
  }
}
$appris_select->finish;


# get the TSL1s from the core db
say "Retrieving TSL data";
my %tsl;
my $tsl_select = $coredb->dbc->prepare("SELECT DISTINCT t.stable_id AS transcript_stable_id FROM gene g JOIN transcript t ON g.gene_id = t.gene_id LEFT JOIN transcript_attrib \
                                        ta ON t.transcript_id = ta.transcript_id LEFT JOIN translation tn ON ta.transcript_id = tn.transcript_id WHERE ta.value = 'tsl1' AND \
                                        ta.attrib_type_id = 428 AND tn.translation_id is NOT NULL");
$tsl_select->execute();
while (my $tsl_transcript = $tsl_select->fetchrow()){
  $tsl{$tsl_transcript} = 1;
}
$tsl_select->finish;


# get the Ensembl canonicals
say "Getting existing canonical data";
my %ensembl_canonical;
my $ensembl_select = $coredb->dbc->prepare("SELECT DISTINCT t.stable_id AS transcript_stable_id FROM gene g JOIN transcript t ON g.canonical_transcript_id = t.transcript_id \
                                           LEFT JOIN translation tn ON t.transcript_id = tn.transcript_id WHERE tn.stable_id IS NOT NULL AND g.source IN \
                                           ('ensembl','ensembl_havana','havana')");
$ensembl_select->execute();
while (my $ensembl_transcript = $ensembl_select->fetchrow()){
  $ensembl_canonical{$ensembl_transcript} = 1;
}
$ensembl_select->finish;

if ($download_data) {
  download_data($outdir);
}

# parse HGNC info
say "Parsing HGNC info";
my %hgnc;
open (HGNC, $outdir.'/hgnc_complete_set.txt');

while (my $line = <HGNC>) {
  my @line = split (/\t/, $line);
  if ($line[19]) {
    $hgnc{$line[19]} = $line[0];
  }
}
close (HGNC);


say "Parsing RefSeq info";
my %refseq_canonical;
my %refseq_match;
my %refseq_canonical_match;

open (REFSEQ, $outdir.'/refseq_variant1.out');

# go through the refseq isoform 1 file and add the accession to a hash
while (my $line = <REFSEQ>) {
  my ($a, $b) = split (/>/, $line);
  my @refseq_line = split (/\s+/, $b);
  $refseq_canonical {$refseq_line[0]} = 1;
}

close (REFSEQ);

say "Retrieving RefSeq attribute data";
my $refseq_select = $otherfdb->dbc->prepare("SELECT DISTINCT stable_id, value FROM transcript t, transcript_attrib ta, attrib_type at WHERE t.transcript_id = ta.transcript_id AND \
                                            ta.attrib_type_id = at.attrib_type_id AND ta.attrib_type_id = 510");
$refseq_select->execute();
while (my $refseq_row = $refseq_select->fetchrow_arrayref()){
  my ($refseq_acc, $ensembl_ids) = @$refseq_row;

  my @ens_ids = split (/:/, $ensembl_ids);
  pop @ens_ids;
  foreach my $stable_id (@ens_ids) {
    unless (exists $refseq_canonical_match {$stable_id}) {
      if (exists $refseq_canonical {$refseq_acc}) {
        $refseq_canonical_match {$stable_id} = $refseq_acc;
      } else {
        unless (exists $refseq_match {$stable_id}) {
          $refseq_match {$stable_id} = $refseq_acc;
        }
      }
    }
  }
}
$refseq_select->finish;


say "Processing UniProt data";
my %uniprot;

open (UNI_MAPPING, $outdir.'/ENST_uniprot_mapping.out');

while (my $line = <UNI_MAPPING>) {
  my ($uniprot,$desc,$transcript) = split (/\s+/, $line);
  chomp $transcript;

  if ($uniprot =~ /-/) {
    if ($uniprot =~ /-1$/) {
      $uniprot{$transcript} = 1;
    }
  } else {
    $uniprot{$transcript} = 1;
  }
}
close (UNI_MAPPING);

# open a log file for output messgaes
open (LOGFILE, ">$log_file");


say "Processing genes";

# Now assign the scores
# Print a file containing all the genes and transcripts
open (GENEFILE, ">$gene_file");

print GENEFILE "#gene_stable_id\tgene_coding_exons\tgene_pathogenic_variants\tgene_biotype\ttranscript_stable_id\ttrans_coding_exons\ttrans_pathogenic_variants\ttranscript_length\tchromosome\tgene_start\tgene_end\tstrand\ttranscript_score\tregion\n";

my $genes_select = $coredb->dbc->prepare("SELECT DISTINCT g.gene_id, g.stable_id, g.biotype, sr.name, g.seq_region_start, g.seq_region_end, g.seq_region_strand, exc_type FROM gene g JOIN transcript t ON g.gene_id = t.gene_id \
                                         JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id LEFT JOIN assembly_exception ae ON g.seq_region_id = ae.seq_region_id LEFT \
                                         JOIN translation tn ON t.transcript_id = tn.transcript_id WHERE t.source IN ('ensembl','havana','ensembl_havana') AND tn.stable_id IS NOT NULL");

$genes_select->execute();

open (ALERT_FILE1, ">$alert_file1");
print ALERT_FILE1 "List of genes for which more than 1 transcript is required to cover all the coding exons\n";
open (ALERT_FILE2, ">$alert_file2");
print ALERT_FILE2 "List of genes for which more than 1 transcript is required to cover all the coding pathogenic variants\n";

my %gene_hash_array;

my %number_variants;
my %coding_exons;

my %pathogenic_variants;
my %transcript_set;

my %trans_ordered_by_variants;
my %trans_ordered_by_exons;

my %trans_score;
my %biotype;
my %chromosome;
my %start;
my %end;
my %strand;
my %region;
my %length;
my %gene_id;

while (my $gene_row= $genes_select->fetchrow_arrayref()){
  my ($gene_id, $gene, $g_bt, $g_chr, $g_start, $g_end, $g_strand, $exc) = @$gene_row;
  $gene_id {$gene} = $gene_id;
  $biotype{$gene} = $g_bt;
  $chromosome{$gene} = $g_chr;
  $start{$gene} = $g_start;
  $end{$gene} = $g_end;
  $strand{$gene} = $g_strand;
  if ($exc) {
    $region{$gene} = $exc;
  } else {
    $region{$gene} = "NULL";
  }
}
$genes_select->finish;

my $slices = $sa->fetch_all('toplevel');
my $slice_count = scalar(@$slices);
my $processing_count = 0;
foreach my $slice (@$slices) {
  $processing_count++;
  say "Processing slice: ".$slice->name." (".$processing_count."/".$slice_count.")";
  my $genes = $slice->get_all_Genes();
foreach my $generef (@$genes) {
  my $gene = $generef->stable_id();
  my @transcripts = @{$generef->get_all_Transcripts()};
  my @coding_transcripts;
  my @coding_trans_ids;

  foreach my $element (@transcripts) {
    if (grep $_ eq ($element->biotype), @allowed_biotypes) { # check the transcript has a biotype in the allowed list of biotypes
      push(@coding_transcripts, $element);
    }
  }
  # Get all the transcripts for the gene and for each of them add to the coding exons hash
  my @gene_exons;
  my @gene_variants;
  unless ($coding_exons{$gene}) {
    foreach my $transcript (@coding_transcripts) {
      my $trans = $transcript->stable_id;
      my $score = 0;
      push (@coding_trans_ids, $trans);

      if (exists $appris{$trans}) {
        if ($appris{$trans} == 1) {
          $score = $score + $appris_weight;
        }
      }
      if ($tsl{$trans}) {
        $score = $score + $tsl_weight;
      }
      if ($refseq_canonical_match{$trans}) {
        $score = $score + $refseq_canonical_weight;
      } elsif ($refseq_match{$trans}) {
        $score = $score + $refseq_match_weight;  
      }
      if ($uniprot{$trans}) {
        $score = $score + $uniprot_weight;
      }

      $trans_score{$trans} = $score;

      my $trans_length = 0;

      # Get the coding exons for transcript and gene
      # This could perhaps be modified to include all exons
      my @trans_coding_exons = @{$transcript->get_all_translateable_Exons()};
      my @trans_coding_ids;
      foreach my $exon_ref (@trans_coding_exons) {
        push (@trans_coding_ids, $exon_ref->stable_id);
      }
      $coding_exons{$trans} = [@trans_coding_ids];
      push (@gene_exons, @trans_coding_ids);

      my @trans_variants;
      # Get the pathogenic variants and the length - note that we're only including coding exons
      foreach my $t_exon (@trans_coding_exons) {
        $trans_length+= $t_exon->length;
        my $t_exon_slice = $sa->fetch_by_region('toplevel', $t_exon->seq_region_name, $t_exon->start, $t_exon->end);
        foreach my $exon_vf ( @{ $vfa->fetch_all_by_Slice($t_exon_slice) } ) {
          my @csstates = @{$exon_vf->get_all_clinical_significance_states()};
          if ( grep( /pathogenic/, @csstates ) ) {
            push (@trans_variants, $exon_vf->name());
            push (@gene_variants, $exon_vf->name());
          }
        }
      }
      $length{$trans} = $trans_length;
      $pathogenic_variants{$trans} = [uniq(@trans_variants)];
      $number_variants{$trans} = scalar (@{$pathogenic_variants{$trans}});
    }
  }
  $coding_exons{$gene} = [uniq(@gene_exons)];
  $pathogenic_variants{$gene} = [uniq(@gene_variants)];
  $number_variants{$gene} = scalar (@{$pathogenic_variants{$gene}});

  $gene_hash_array{$gene} = [@coding_trans_ids];
  if (scalar (@coding_transcripts) == 0) {
    delete $gene_hash_array{$gene};
  } elsif (scalar (@coding_transcripts) > 1) {
    my @pathogenic_transcripts = coverage_sorter(\@coding_trans_ids, \%pathogenic_variants, \%trans_score, scalar(@{$pathogenic_variants{$gene}}), "array");
    $trans_ordered_by_variants{$gene} = [@pathogenic_transcripts];
    my @exon_transcripts = coverage_sorter(\@coding_trans_ids, \%coding_exons, \%trans_score, scalar (@{$coding_exons{$gene}}), "array");
    $trans_ordered_by_exons{$gene} = [@exon_transcripts];
    # print into the ALERT files if more than 1 transcript is needed for coverage
    if (scalar (@pathogenic_transcripts) > 1) {
      print ALERT_FILE2 $gene, ":\t@pathogenic_transcripts\n";
    }
    if (scalar (@exon_transcripts) > 1) {
      print ALERT_FILE1 $gene, ":\t@exon_transcripts\n";
    }
  }
  foreach my $coding_transcript(@{$gene_hash_array{$gene}}) {
    print GENEFILE $gene, "\t", scalar (@{$coding_exons{$gene}}), "\t", scalar (@{$pathogenic_variants{$gene}}), "\t", $biotype{$gene}, "\t", $coding_transcript, "\t", scalar (@{$coding_exons{$coding_transcript}}), "\t", scalar (@{$pathogenic_variants{$coding_transcript}}), "\t", $length{$coding_transcript}, "\t", $chromosome{$gene}, "\t", $start{$gene}, "\t", $end{$gene}, "\t", $strand{$gene}, "\t", $trans_score{$coding_transcript}, "\t", $region{$gene}, "\n";
  }
}
} # @$slices

say "Finished processing genes";

close (GENEFILE);
close (ALERT_FILE1);
close (ALERT_FILE2);

# now for each gene create an array ordered according to length and another according to score and compare
# then assign the canonical based on the outcome
my %canonical;
my %highest_scoring_trans;
my %longest_trans;
my %highest_no_variants;
my %description;
my %highest_no_coding_exons;

say "Assigning canonicals";
open (CAN_FILE, ">$canonicals_file");
print CAN_FILE "#gene_id\tnumber_gene_variants\tgene_name\thgnc_acc\tgene_biotype\treference_canonical\tnumber_trans_variants\ttrans_length\tchr\tstart\tend\tstrand\tregion\tscore\thighest_scoring\tlongest\tmost_variants\tdescription\n";

for my $gene ( keys %gene_hash_array ) {
  my @trans_ordered_by_variants;
  my @trans_ordered_by_exons;
  my @longest_keys;
  my @high_score_keys;
  my @high_variant_keys;
  my @high_exon_keys;
  my $high_variants;
  my $high_score;
  my $longest;
  my %length_hash;
  my %score_hash;
  my %variants_hash;
  my %coding_exons_hash;
  my @coding_exons_array;

  # If there's no hgnc_id then set this to be a hyphen
  unless (exists $hgnc{$gene}){
    $hgnc{$gene} = "-";
  }

  # If there's only 1 transcript then just set that as the canonical and move on to the next
  if (scalar @{$gene_hash_array{$gene}} == 1) {
    $canonical{$gene} = "@{$gene_hash_array{$gene}}";
    chomp $canonical{$gene};
    $description{$gene} = 'Only 1 transcript for gene';
    $high_score = $longest = $high_variants = 'N/A';
  } elsif (scalar @{$gene_hash_array{$gene}} > 1) {
    foreach my $element (@{$gene_hash_array{$gene}}) {
      $length_hash{$element} = $length{$element};
      $score_hash{$element} = $trans_score{$element};
      $variants_hash{$element} = scalar (@{$pathogenic_variants{$element}});
      $coding_exons_hash{$element} = scalar (@{$coding_exons{$element}});
    }

    # transcript with the highest number of coding exons
    @trans_ordered_by_exons = @{$trans_ordered_by_exons{$gene}};    
    $highest_no_coding_exons{$gene} = $trans_ordered_by_exons[0];
    @high_exon_keys = grep { $coding_exons_hash{$_} eq $coding_exons_hash{$highest_no_coding_exons{$gene}} } keys %coding_exons_hash;

    # transcript with the highest number of variants
    @trans_ordered_by_variants = @{$trans_ordered_by_variants{$gene}};
    $highest_no_variants{$gene} = $trans_ordered_by_variants[0];
    @high_variant_keys = grep { $variants_hash{$_} eq $variants_hash{$highest_no_variants{$gene}} } keys %variants_hash;

    # a list of the transcripts sorted by their corresponding lengths
    my @transcript_lengths = sort {$length{$b} <=> $length{$a}} values $gene_hash_array{$gene};
    $longest_trans{$gene} = $transcript_lengths[0];
    @longest_keys = grep { $length_hash{$_} eq $length_hash{$longest_trans{$gene}} } keys %length_hash;

    # a list of the transcripts sorted by their corresponding scores
    my @scores = sort {$trans_score{$b} <=> $trans_score{$a}} values $gene_hash_array{$gene};
    $highest_scoring_trans{$gene} = $scores[0];
    @high_score_keys = grep { $score_hash{$_} eq $score_hash{$highest_scoring_trans{$gene}} } keys %score_hash;

    # if the highest scoring transcript is also the longest (or it's at least 75% the length of the longest) then just choose as canonical
    if ($transcript_lengths[0] eq $scores[0]) {
      $canonical{$gene} = $scores[0];
      $description{$gene} = 'Selected transcript is the longest and highest scoring transcript';
    } elsif ($length_hash{$scores[0]} > $length_hash{$transcript_lengths[0]}*$length_fraction) {
      $canonical{$gene} = $scores[0];
      $description{$gene} = 'Selected transcript is the highest scoring transcript';
    } else {
      # go through the transcripts in order of descending score until we find one that's at least 75% of the longest one
      for (my $i = 1; $i < scalar (@scores); $i++) {
        if ($scores[$i] eq $transcript_lengths[0]) {
          $canonical{$gene} = $scores[$i];
          $description{$gene} = 'Selected transcript is the longest transcript';
          last;
        } elsif ($length_hash{$scores[$i]} > $length_hash{$transcript_lengths[0]}*$length_fraction) {
          $canonical{$gene} = $scores[$i];
          $description{$gene} = 'Selected transcript is at least '.$length_fraction.' times the length of the longest transcript';
          last;
        }
      }
    }
    # throw a warning if one hasn't been assigned at this stage
    unless ($canonical{$gene}) {
      print "ERROR: canonical for $gene not yet defined\n";
    }

    # Now check if there are other transcripts belonging to the gene that have the same score
    # If they do then check if they are at least 75% of the length of the longest transcript
    # if there's still more than 1 then check if any are APPRIS P2s:
    #  - if only 1 is then set it as canonical
    #  - if more than 1 is then check Ensembl canonicals (explained below)
    #  - If none are then go to APPRIS P3s
    # APPRIS P3s:
    #  - If just 1 set it as canonical
    #  - if more than 1 then check Ensembl canonicals
    # if a canonical hasn't been decided yet then take whichever was chosen as the Ensembl canonical
    # if there isn't a clear choice then pick whichever was chosen as the HAVANA canonical
    # if there isn't a clear choice by this stage take whichever is the longest
    # and if this doesn't make a clear choice then randomly select one
    my $index = 0;
    my @new_scores;
    foreach my $t (@scores) {
      if ($trans_score{$t} eq $trans_score{$canonical{$gene}} && $length{$t} > $length_hash{$transcript_lengths[0]}*$length_fraction) {
        push (@new_scores, $t);
      }
      $index++;
    }
    if (scalar(@new_scores > 1)) {
      # I need to put in a check if the @new_scores array contains either 1 hit or no hits
      my @appris2;
      my @appris3;
      my @appris_hits;
      my @undecided;

      foreach my $hit (@new_scores) {
        if (exists $appris{$hit}) {
          if ($appris{$hit} eq 2) {
            push (@appris2, $hit);
          } elsif ($appris{$hit} eq 3) {
            push (@appris3, $hit);
          }
        }
      }
      # First check if there are APPRIS P2s
      # if none then check for P3s
      # if 1 set it as canonical
      # if more than 1 check for Ensembl canonical
      if (scalar(@appris2 == 1)) {
        $canonical{$gene} = "@appris2";
        chomp $canonical{$gene};
        $description{$gene} = $description{$gene} . ' and an APPRIS P2';
      } elsif (scalar(@appris2 == 0)) {
        if (scalar(@appris3 == 1)) {
          $canonical{$gene} = "@appris3";
          chomp $canonical{$gene};
          $description{$gene} = $description{$gene} . ' and an APPRIS P3';
        } elsif (scalar(@appris3 == 0)) {
          # the array to be checked for Ensembl hits will consist of all the @new_scores
          @undecided = @new_scores;
        } elsif (scalar(@appris3 > 1)) {
          # the array to be checked for Ensembl hits will consist of just these hits
          @undecided = @appris3;
        }
      } elsif (scalar(@appris2 > 1)) {
        # the array to be checked for Ensembl hits will consist of just these hits
        @undecided = @appris2;
      }

      # now go through the undecided array and create an ensembl array of hits from it
      # only if the undecided array is defined do we want to do this
      if (scalar(@undecided == 1 )) {
        print "ERROR: There's only 1 undecided canonical for $gene at this stage, there should be more!\n";
      } elsif (scalar(@undecided > 1)) {
        my @ensembl_hits;
        foreach my $hit (@undecided) {
          if ($ensembl_canonical{$hit}) {
            push (@ensembl_hits, $hit);
          }
        }

        if (scalar(@ensembl_hits == 1)) {
          $canonical{$gene} = "@ensembl_hits";
          chomp $canonical{$gene};
          $description{$gene} = $description{$gene} . ' and the Ensembl canonical';
        } elsif (scalar(@ensembl_hits > 1)) {
          my @sorted = sort { $length_hash{$b} <=> $length_hash{$a} } @ensembl_hits;
          $canonical{$gene} = $sorted[0];
          $description{$gene} = $description{$gene} . ' and the longest Ensembl canonical';
        } elsif (scalar(@ensembl_hits == 0)) {
          # if there is no ensembl hit at this stage then just take the longest
          my @sorted = sort { $length_hash{$b} <=> $length_hash{$a} } @new_scores;
          $canonical{$gene} = $sorted[0];
        }
      }
    }
  }
  # retrieve the xref
  my $name;
  if (exists $xref{$gene}) {
    $name = $xref{$gene};
  } else {
    $name = $gene;
  }

  if (exists $canonical{$gene}) {
    unless ($longest) {
      if (scalar @longest_keys == 1) {
        $longest = $longest_trans{$gene};
      } else {
        $longest = join(",", @longest_keys);
      }
    }
    unless ($high_score) {
      if (scalar @high_score_keys == 1) {
        $high_score = $highest_scoring_trans{$gene};
      } else {
        $high_score = join(",",@high_score_keys);
      }
    }
    unless ($high_variants) {    
      $high_variants = join(",", @high_variant_keys);
    }
    if ($number_variants{$gene} == 0) {
      $high_variants = "N/A";
    } elsif ($number_variants{$gene} == $number_variants{$canonical{$gene}}) {
      $description{$gene} .= ". Canonical transcript encompasses all of gene's pathogenic variants.";
    } elsif ($number_variants{$gene} == $number_variants{$highest_scoring_trans{$gene}}) {
      # only take the high scoring transcripts that also have the correct number of variants
      $description{$gene} .= ". Highest scoring transcript encompasses all of gene's pathogenic variants.";
    } elsif ($number_variants{$gene} == $number_variants{$highest_no_variants{$gene}}) {
      $description{$gene} .= ". $highest_no_variants{$gene} contains all of gene's pathogenic variants.";
    }

    # print the canonicals into a file
    print CAN_FILE $gene, "\t", $number_variants{$gene}, "\t", $name, "\t", $hgnc{$gene}, "\t", $biotype{$gene}, "\t", $canonical{$gene}, "\t", $number_variants{$canonical{$gene}}, "\t", $length{$canonical{$gene}}, "\t", $chromosome{$gene}, "\t", $start{$gene}, "\t", $end{$gene}, "\t", $strand{$gene}, "\t", $region{$gene}, "\t", $trans_score{$canonical{$gene}}, "\t", $high_score, "\t", $longest, "\t", $high_variants, "\t", $description{$gene}, "\n";
    
    # store a gene attribute called "select_transcript" to store the "select" transcript
    my $gene_object = $ga->fetch_by_stable_id($gene);
    my $select_attrib = Bio::EnsEMBL::Attribute->new(-CODE => 'select_transcript',
                                                     -NAME => 'Select transcript',
                                                     -DESCRIPTION => '"select" transcript stable ID.',
                                                     -VALUE => $canonical{$gene});
    my @attributes = ();
    push(@attributes,$select_attrib);
    $aa->store_on_Gene($gene_object,\@attributes);
    print LOGFILE "Transcript ".$canonical{$gene}." stored as select transcript for gene ".$gene." in the core database.\n";
  } else {
    print LOGFILE "No protein-coding canonical transcript defined for $gene\n"; 
  }
}
close (CAN_FILE);
close (LOGFILE);

exit;

#################
#  SUBROUTINES  #
#################

sub download_data {
  my $outdir = shift();
 
  # get the HGNC accessions
  say "Downloading HGNC accessions";
  system ("bsub -J HGNC_download wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt -P $outdir");

  # get RefSeq canonicals
  say "Downloading RefSeq canonicals";
  system ("bsub -J refseq_download wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/*rna.fna* -P $outdir");

  # get the UniProt canonicals
  say "Downloading UniProt canonicals";
  system ("bsub -J uniprot_download wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz -P $outdir");

  # unzip all the gz files
  say "Unzipping downloaded files";
  system ("bsub -w 'done(HGNC_download) && done(refseq_download) && done(uniprot_download)' -K gunzip $outdir/*gz");

  # parse refseq info
  say "Parsing RefSeq variant info";
  system ("grep 'variant 1,' $outdir/*fna > $outdir/refseq_variant1.out");

  # parse the uniprot info
  say "Parsing UniProt files";
  system ("grep 'ENST' $outdir/HUMAN_9606_idmapping.dat > $outdir/ENST_uniprot_mapping.out");
}

# return only unique entries in an array
sub uniq {
  my %temp_hash = map { $_, 0 } @_;
  return keys %temp_hash;
} 

# The following subroutine takes as input the transcript list, 2 hashes you wish to sort by
# a value which should be the number of expected exons or variants
# and a string that indicates whether the value of the first hash should be an array or a value
# It will then return an array of the transcripts that are needed for full coverage for the first parameter
# This parameter will usually be pathogenic variants or coding exons

sub coverage_sorter {
  my ($transcript_list, $primary_hash_ref, $secondary_hash_ref, $value, $type) = @_;
  my %primary_hash = %$primary_hash_ref;
  my %secondary_hash = %$secondary_hash_ref;
  my @sorted_transcripts;
  my @coverage_list;
  my @final_transcript_list;
  # sort the transcript list by the first then the second hash
  # If the type is value then it means the first hash value is a value,
  # array means the first hash value is an array
  # we assume the second hash value is a number but this can be updated if you're ordering according to
  # two arrays
  if ($type eq 'value') {
    @sorted_transcripts = sort {$primary_hash{$b} <=> $primary_hash{$a} or $secondary_hash{$b} <=> $secondary_hash{$a}} @$transcript_list;
  } elsif ($type eq 'array') {
    @sorted_transcripts = sort {scalar(@{$primary_hash{$b}}) <=> scalar(@{$primary_hash{$a}}) or $secondary_hash{$b} <=> $secondary_hash{$a}} @$transcript_list;
  }
  @final_transcript_list = $sorted_transcripts[0];
  my $count = scalar(@{$primary_hash{$sorted_transcripts[0]}}); #scalar (uniq(@coverage_list));;
  @coverage_list = @{$primary_hash{$sorted_transcripts[0]}};
  unless ($count == $value) {
    my $i = 1;
    while ($i < scalar (@$transcript_list) && scalar (uniq(@coverage_list)) < $value) {
      uniq (@coverage_list);
      my @temp_array = @coverage_list;
      push (@temp_array, @{$primary_hash{$sorted_transcripts[$i]}});
      uniq (@temp_array);
      if (scalar (@temp_array) > scalar (@coverage_list)) {
        push (@final_transcript_list, $sorted_transcripts[$i]);
        push (@coverage_list, @{$primary_hash{$sorted_transcripts[$i]}});
        $i++;
      }
    }
  }
  return @final_transcript_list;
}

sub transcript_sorter {
  my ($output_array_ref, $input_hash_ref) = @_;
  my @output_array = @$output_array_ref;
  my %input_hash = %$input_hash_ref;
  @output_array = sort { $input_hash{$b} <=> $input_hash{$a} } keys %input_hash;
  return $output_array[0];
}

sub get_top_keys {
  my ($query, $input_hash_ref) = @_;
  my %input_hash = %$input_hash_ref;
  my @output_array = grep { $input_hash{$_} eq $input_hash{$query} } keys %input_hash;
  return @output_array;
}
