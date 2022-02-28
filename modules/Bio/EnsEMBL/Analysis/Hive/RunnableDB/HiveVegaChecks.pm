# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

# Usage examples:
# Generic:
# standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaChecks -dbname my_vega_after_merge -dbhost genebuild1 -dnadbname my_dna_db -dnadbhost genebuild1 -dnadbport 3306 -port 3306 -user ensro -coord_system toplevel -path Zv9 [-chromosome 13] [-user *** -pass *** -write 1] [-sql_output $SCR9/vega_check.sql] [-affix 1] [-biotypes_extension 1] [-dbtype vega]"

# Vega check before the merge:
# standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaChecks -dbname my_vega_fixed_before_merge -dbhost genebuild1 -dnadbname my_dna_db -dnadbhost genebuild1 -dnadbport 3306 -port 3306 -user ensro -coord_system toplevel -path Zv9 -sql_output $SCR9/vega_check_before.sql -stdout_file REPORTS_DIR/vega_checks_before.out -stderr_file REPORTS_DIR/vega_checks_before.err"

# Vega check after the merge:
# standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaChecks -dbname my_vega_after_merge -dbhost genebuild1 -dnadbname my_dna_db -dnadbhost genebuild1 -dnadbport 3306 -port 3306 -user ensro -coord_system toplevel -path Zv9 -sql_output $SCR9/vega_check_after.sql] -affix 1 -biotypes_extension 1 -stdout_file REPORTS_DIR/vega_checks_after.out -stderr_file REPORTS_DIR/vega_checks_after.err"

# Vega check after the merge and after cleaning biotypes affixes:
# standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaChecks -dbname my_vega_after_merge -dbhost genebuild1 -dnadbname my_dna_db -dnadbhost genebuild1 -dnadbport 3306 -port 3306 -user ensro -coord_system toplevel -path Zv9 -sql_output $SCR9/vega_check_after.sql] -biotypes_extension 1 -stdout_file REPORTS_DIR/vega_checks_after.out -stderr_file REPORTS_DIR/vega_checks_after.err"

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaChecks;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Getopt::Long;
use integer; # to get integers instead of floats

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Merge::vega_check qw(get_combos get_biotype_groups get_actions);

sub param_defaults {
    return {
      dbname => '',
      dbhost => '',
      dnadbname => '',
      dnadbhost => '',
      dnadbport => '3306',
      coord_system => '',
      path => '',
      sql_output => '',
      dbtype => '', # can be 'vega' or '' (empty string)
      port => '3306',
      user => 'ensro',
      pass => '',
      chromosome => '',
      write => '',
      affix => 0, # perform the checks by using the biotypes with or without the prefixes and suffixes like weird_, _Ens, _hav, ... ; without affixes by default
      biotypes_extension => 0, # use additional biotypes in 'gene_trans_type_combination_extension'; default=no
      stdout_file => '', # full path to the file where the standard output will be printed
      stderr_file => '' # full path to the file where the standard error will be printed
    }
}

sub fetch_input {
  my $self = shift;

  $self->param('gene_trans_type_combination',get_combos('ensembl'));
  $self->param('gene_trans_type_combination_extension',get_combos('ensembl_extension'));
  $self->param('ensembl_biotype_groups',get_biotype_groups('ensembl'));
  $self->param('gene_coding_biotypes',${$self->param('ensembl_biotype_groups')}{'gene_coding'});
  $self->param('gene_non_coding_biotypes',${$self->param('ensembl_biotype_groups')}{'gene_non_coding'});
  $self->param('transcript_coding_biotypes',${$self->param('ensembl_biotype_groups')}{'transcript_coding'});
  $self->param('transcript_non_coding_biotypes',${$self->param('ensembl_biotype_groups')}{'transcript_non_coding'});
  
  return 1;
}

sub run {
  my $self = shift;

  my $TXT_NO_TRANSLATION_TRANSCRIPT_TO_PROCESSED_TRANSCRIPT = "Translation not found for given transcript. Transcript biotype should be";
  my $TXT_TRANSLATION_STOP_TRANSCRIPT_TO_POLY = "Translation has a stop codon. Transcript biotype should be";
  my $TXT_TRANSLATION_NO_STOP_TRANSCRIPT_TO_CODING = "Translation without stop codon found. Transcript biotype should be";

  my $TXT_GENE_TO_POLY = "Gene contains at least one '\%polymorphic_pseudogene\%' transcript. Gene biotype should be";
  my $TXT_GENE_TO_CODING = "Gene contains at least one '\%protein_coding\%' transcript and no '\%polymorphic_pseudogene\%' transcripts. Gene biotype should be";
  my $TXT_GENE_TO_PROCESSED_TRANSCRIPT = "Gene contains no '\%polymorphic_pseudogene\%' transcripts, no '\%protein_coding\%' transcripts and at least one '\%processed_transcript\%' transcript. Gene biotype should be";

  my $TXT_TRANS_BIOTYPE_NOT_FOUND = "Unknown transcript biotype. Please, correct this problem and run this script again.";
  my $TXT_GENE_BIOTYPE_NOT_FOUND = "Unknown gene biotype. Please, correct this problem and run this script again.";
  my $TXT_GENE_TRANS_BIOTYPE_NOT_ALLOWED = "Gene and transcript biotypes combination is not allowed. Please, correct this problem and run this script again.";

  # allowed affixes for gene biotypes:
  my @GENE_PREFIXES = ('','weird_');
  my @GENE_SUFFIXES = ('','_Ens','_Hav','_Ens_Hav');

  # allowed affixes for transcript biotypes:
  my @TRANSCRIPT_PREFIXES = ('','weird_');
  my @TRANSCRIPT_SUFFIXES = ('','_ens','_hav','_hav_mrgd','_hav_m','_hav_mrg','_mrgd');

  my $num_poly_genes = 0;
  my $num_protein_coding_genes = 0;
  my $num_processed_transcript_genes = 0;
  my $num_polymorphic_transcript = 0;
  my $num_protein_coding_transcript = 0;
  my $num_processed_transcript = 0;
  my $num_unknown_genes = 0;
  my $num_unknown_transcripts = 0;
  my $num_gt_comb_not_allowed = 0;

  open STDOUT, '>', $self->param('stdout_file') if ($self->param('stdout_file'));
  open STDERR, '>', $self->param('stderr_file') if ($self->param('stderr_file'));

  if (not $self->param('user') or not $self->param('dbhost') or not $self->param('dbname')) {
    $self->throw("DB connection parameters missing. Can't connect.");
  }

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbname => $self->param('dbname'),
                                              -host   => $self->param('dbhost'),
                                              -user   => $self->param('user'),
                                              -pass   => $self->param('pass'),
                                              -port   => $self->param('port'));

  my $coord_system_adaptor = $db->get_CoordSystemAdaptor();
  $coord_system_adaptor->fetch_all->[0]->version($self->param('path'));

  # The DNA is in the ref db
  if ($self->param('dnadbname')) {
    if (not $self->param('dnadbhost')) {
      $self->throw("Can't connect to DNA DB. Please specify the host the DNA DB is on.");
    }
    my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $self->param('dnadbname'),
                                                    -host   => $self->param('dnadbhost'),
                                                    -user   => $self->param('user'),
                                                    -pass   => $self->param('pass'),
                                                    -port   => $self->param('dnadbport'));
    $db->dnadb($dnadb);
  }

  # we want that the biotypes regexp match the 'base' biotype and also
  # the 'base' biotype plus extensions (prefix and/or suffix) like 'weird_','_hav','_hav_mrgd','_Ens',...
  # depending on the parameter 'affix' value
  my $gene_prefixes_str  = ""; $gene_prefixes_str = join("|",@GENE_PREFIXES) if ($self->param('affix'));
  my $gene_suffixes_str = ""; $gene_suffixes_str = join("|",@GENE_SUFFIXES) if ($self->param('affix'));
  my $transcript_prefixes_str = ""; $transcript_prefixes_str = join("|",@TRANSCRIPT_PREFIXES) if ($self->param('affix'));
  my $transcript_suffixes_str = ""; $transcript_suffixes_str = join("|",@TRANSCRIPT_SUFFIXES) if ($self->param('affix'));
  my $gene_coding_biotypes_str = join("|",@{$self->param('gene_coding_biotypes')});
  my $gene_non_coding_biotypes_str = join("|",@{$self->param('gene_non_coding_biotypes')});
  my $gene_biotypes_str = $gene_coding_biotypes_str."|".$gene_non_coding_biotypes_str;

  my $transcript_coding_biotypes_str = join("|",@{$self->param('transcript_coding_biotypes')});
  my $transcript_non_coding_biotypes_str = join("|",@{$self->param('transcript_non_coding_biotypes')});
  my $transcript_biotypes_str = $transcript_coding_biotypes_str."|".$transcript_non_coding_biotypes_str;

  my $ga = $db->get_GeneAdaptor();
  my $ta = $db->get_TranscriptAdaptor();
  my $sa = $db->get_SliceAdaptor();

  my @slices;
  my $all_slices = $sa->fetch_all( 'toplevel', $self->param('path'), 1, undef );
  if ($self->param('chromosome')) { # if chr was defined as parameter, choose one slice
    foreach my $sl ( @{$all_slices} ) {
      if ( $sl->seq_region_name eq $self->param('chromosome') ) {
        $slices[0] = $sl;
      }
    }
  } else { # if chr was not defined, choose all toplevel slices
    @slices = @{$all_slices};
  }

  if ($self->param('sql_output')) {
    if (-e $self->param('sql_output')) {
      open(SQL_FILE,">".$self->param('sql_output'));
      close(SQL_FILE);
      print("\nThe file ".$self->param('sql_output')." has been truncated.\n");
    }
  }

  print("\nGene coding biotypes:\n$gene_coding_biotypes_str\n\n");
  print("Gene non-coding biotypes:\n$gene_non_coding_biotypes_str\n\n");
  print("Transcript coding biotypes:\n$transcript_coding_biotypes_str\n\n");
  print("Transcript non-coding biotypes:\n$transcript_non_coding_biotypes_str\n\n");

  print("Gene biotype prefixes:\n$gene_prefixes_str\n\n") if ($gene_prefixes_str);
  print("Gene biotype suffixes:\n$gene_suffixes_str\n\n") if ($gene_suffixes_str);
  print("Transcript biotype prefixes:\n$transcript_prefixes_str\n\n") if ($transcript_prefixes_str);
  print("Transcript biotype suffixes:\n$transcript_suffixes_str\n\n") if ($transcript_suffixes_str);

  foreach my $slice (@slices)
  {
    # Fetch all genes on the slice and load their transcripts at the same time
    my @genes = @{$slice->get_all_Genes(undef, undef, 1, undef, undef)};

    print("--------------- Processing ".scalar(@genes)." genes from slice: ".$slice->seq_region_name."\n");

    foreach my $gene (@genes) {

      my $polymorphic_transcript = 0;
      my $protein_coding_transcript = 0;
      my $processed_transcript = 0;

      my $g_biotype = $gene->biotype();
      my $g_biotype_basename = get_biotype_basename($g_biotype,$gene_biotypes_str,$gene_prefixes_str,$gene_suffixes_str); # remove any prefix and suffix

      if ((exists(${$self->param('gene_trans_type_combination')}{$g_biotype_basename})) or
          ($self->param('biotypes_extension') and exists(${$self->param('gene_trans_type_combination_extension')}{$g_biotype_basename}))) {

        my $g_allowed_transcript_biotypes_str = "";
        if (exists(${$self->param('gene_trans_type_combination')}{$g_biotype_basename})) {
          $g_allowed_transcript_biotypes_str .= join("|",@{${$self->param('gene_trans_type_combination')}{$g_biotype_basename}});
        }

        if ($self->param('biotypes_extension') and exists(${$self->param('gene_trans_type_combination_extension')}{$g_biotype_basename})) {
          $g_allowed_transcript_biotypes_str .= "|" if (!($g_allowed_transcript_biotypes_str eq ""));
          $g_allowed_transcript_biotypes_str .= join("|",@{${$self->param('gene_trans_type_combination_extension')}{$g_biotype_basename}});
        }

        my @transcripts = @{$gene->get_all_Transcripts()};
        my $num_transcripts = scalar(@transcripts);

        foreach my $trans (@transcripts) {
          my $t_biotype = $trans->biotype();
          # check if transcript biotype is known
          if ($t_biotype =~ /^($transcript_prefixes_str)($transcript_biotypes_str)($transcript_suffixes_str)$/) {

            # check if gene-transcript biotype combination is allowed
            if ($t_biotype !~ /^($transcript_prefixes_str)($g_allowed_transcript_biotypes_str)($transcript_suffixes_str)$/) {
              print_gene_transcript_info($gene,$trans,$TXT_GENE_TRANS_BIOTYPE_NOT_ALLOWED) unless $self->param('dbtype') =~ /vega/;
              $num_gt_comb_not_allowed++ unless $self->param('dbtype') =~ /vega/;
            }

            if (!$trans->translation) { # transcript has no translation
              if ($t_biotype =~ /^($transcript_prefixes_str)($transcript_coding_biotypes_str)($transcript_suffixes_str)$/) { # but its biotype is protein coding (including polymorphic_pseudogene)
                # transcript biotype should be processed_transcript, which means it is non-coding
                my ($t_biotype_prefix,$t_biotype_suffix) = get_biotype_affix($t_biotype,$transcript_coding_biotypes_str,$transcript_prefixes_str,$transcript_suffixes_str);
                my $new_suggested_biotype = $t_biotype_prefix.'processed_transcript'.$t_biotype_suffix;

                print_gene_transcript_info($gene,$trans,$TXT_NO_TRANSLATION_TRANSCRIPT_TO_PROCESSED_TRANSCRIPT,$new_suggested_biotype);
                update_gt("transcript",$ta,$trans,$new_suggested_biotype,$self->param('write'),$self->param('sql_output'));

                $num_processed_transcript++; # keep track of the total number of updates to this type of transcripts
              }
              # we want to keep track of the presence of processed_transcript/non-coding transcripts to decide about the gene biotype later on
              $processed_transcript = 1;
            } elsif ($trans->translate->seq =~ /\*/ and $t_biotype !~ /TR_.*gene/ and $t_biotype !~ /IG_.*gene/) { # transcript has translation with stop codon
              if ($t_biotype !~ /.*polymorphic.*/)  { # but its biotype is not polymorphic or polymorphic_pseudogene
                # transcript biotype should be polymorphic_pseudogene
                my ($t_biotype_prefix,$t_biotype_suffix) = get_biotype_affix($t_biotype,$transcript_biotypes_str,$transcript_prefixes_str,$transcript_suffixes_str);
                my $new_suggested_biotype = $t_biotype_prefix.'polymorphic_pseudogene'.$t_biotype_suffix;

                print_gene_transcript_info($gene,$trans,$TXT_TRANSLATION_STOP_TRANSCRIPT_TO_POLY,$new_suggested_biotype);
                update_gt("transcript",$ta,$trans,$new_suggested_biotype,$self->param('write'),$self->param('sql_output'));

                $num_polymorphic_transcript++; # keep track of the total number of updates to this type of transcripts
              }
              # we want to keep track of the presence of polymorphic_pseudogene transcripts to decide about the gene biotype later on
              # regardless it was a polymorphic transcript before or after this check
              $polymorphic_transcript = 1;
            } else { # transcript has translation without stop codon
              if (($t_biotype =~ /^($transcript_prefixes_str)($transcript_non_coding_biotypes_str)($transcript_suffixes_str)$/)) { # but its biotype is non-coding
                # transcript biotype should be protein_coding
                my ($t_biotype_prefix,$t_biotype_suffix) = get_biotype_affix($t_biotype,$transcript_biotypes_str,$transcript_prefixes_str,$transcript_suffixes_str);
                my $new_suggested_biotype = $t_biotype_prefix.'protein_coding'.$t_biotype_suffix;

                print_gene_transcript_info($gene,$trans,$TXT_TRANSLATION_NO_STOP_TRANSCRIPT_TO_CODING,$new_suggested_biotype);
                update_gt("transcript",$ta,$trans,$new_suggested_biotype,$self->param('write'),$self->param('sql_output'));
                $num_protein_coding_transcript++; # keep track of the total number of updates to this type of transcripts
              }
              # we want to keep track of the presence of protein_coding transcripts to decide about the gene biotype later on
              # regardless it was a protein_coding transcript before or after this check
              if ($t_biotype !~ /.*polymorphic.*/) {
                $protein_coding_transcript = 1;
              } else {
                $polymorphic_transcript = 1;
              }
            }

          } else { # transcript biotype unknown
            print_gene_transcript_info($gene,$trans,$TXT_TRANS_BIOTYPE_NOT_FOUND) unless $self->param('dbtype') =~ /vega/;
            $num_unknown_transcripts++;
          }
        } # foreach transcripts

        # decide the gene biotype
        # precedence: polymorphic_pseudogene > protein_coding > processed_transcript
        # note that we want to keep the gene biotype if no change has been suggested
        if ($polymorphic_transcript > 0) {
          if ($g_biotype !~ /.*polymorphic.*/) {
            my ($g_biotype_prefix,$g_biotype_suffix) = get_biotype_affix($g_biotype,$gene_biotypes_str,$gene_prefixes_str,$gene_suffixes_str);
            my $new_suggested_biotype = $g_biotype_prefix.'polymorphic_pseudogene'.$g_biotype_suffix;
            print_gene_info($gene,$TXT_GENE_TO_POLY,$new_suggested_biotype);
            update_gt("gene",$ga,$gene,$new_suggested_biotype,$self->param('write'),$self->param('sql_output'));
            $num_poly_genes++;
          }
        } elsif ($protein_coding_transcript > 0) {
          if (($g_biotype =~ /^($gene_prefixes_str)($gene_non_coding_biotypes_str)($gene_suffixes_str)$/) or ($g_biotype =~ /.*polymorphic.*/)) { # if gene biotype is protein coding excluding polymorphic (or "is non-coding or polymorphic")
            my ($g_biotype_prefix,$g_biotype_suffix) = get_biotype_affix($g_biotype,$gene_biotypes_str,$gene_prefixes_str,$gene_suffixes_str);
            my $new_suggested_biotype = $g_biotype_prefix.'protein_coding'.$g_biotype_suffix;
            print_gene_info($gene,$TXT_GENE_TO_CODING,$new_suggested_biotype);
            update_gt("gene",$ga,$gene,$new_suggested_biotype,$self->param('write'),$self->param('sql_output'));
            $num_protein_coding_genes++;
          }
        } elsif ($processed_transcript > 0 and !gene_has_assembly_error_attribute($gene)) {
          # note we want to skip the protein_coding gene biotypes set by Havana for genes with NoTransRefError attrib
          if ($g_biotype !~ /^($gene_prefixes_str)($gene_non_coding_biotypes_str)($gene_suffixes_str)$/) {
            my ($g_biotype_prefix,$g_biotype_suffix) = get_biotype_affix($g_biotype,$gene_biotypes_str,$gene_prefixes_str,$gene_suffixes_str);
            my $new_suggested_biotype = $g_biotype_prefix.'processed_transcript'.$g_biotype_suffix;
            print_gene_info($gene,$TXT_GENE_TO_PROCESSED_TRANSCRIPT,$new_suggested_biotype);
            update_gt("gene",$ga,$gene,$new_suggested_biotype,$self->param('write'),$self->param('sql_output'));
            $num_processed_transcript_genes++;
          }
        } #else {} # we want to keep the gene biotype if no change has been suggested
      } else { # if does not exist gene_trans_type_combination 
        print_gene_info($gene,$TXT_GENE_BIOTYPE_NOT_FOUND) unless $self->param('dbtype') =~ /vega/;
        $num_unknown_genes++;
      }
    } # foreach gene
  } # foreach slice

  my $should_be_or_has_been = "should be";
  $should_be_or_has_been = "has been" if ($self->param('write'));

  print("\n--- SUMMARY ---\n\n");

  print("The number of unknown gene biotypes is: $num_unknown_genes\n") unless $self->param('dbtype') =~ /vega/;
  print("The number of unknown transcript biotypes is: $num_unknown_transcripts\n\n") unless $self->param('dbtype') =~ /vega/;
  print("The number of gene-transcript biotype combinations not allowed is: $num_gt_comb_not_allowed\n\n") unless $self->param('dbtype') =~ /vega/;

  print("The number of genes whose biotype $should_be_or_has_been updated to '\%protein_coding\%' is: $num_protein_coding_genes\n");
  print("The number of genes whose biotype $should_be_or_has_been updated to '\%polymorphic_pseudogene\%' is: $num_poly_genes\n");
  print("The number of genes whose biotype $should_be_or_has_been updated to '\%processed_transcript\%' is: $num_processed_transcript_genes\n\n");

  print("The number of transcripts whose biotype $should_be_or_has_been updated to '\%protein_coding\%' is: $num_protein_coding_transcript\n");
  print("The number of transcripts whose biotype $should_be_or_has_been updated to '\%polymorphic_pseudogene\%' is: $num_polymorphic_transcript\n");
  print("The number of transcripts whose biotype $should_be_or_has_been updated to '\%processed_transcript\%' is: $num_processed_transcript\n\n");

  if ($self->param('sql_output')) {
    my $not_sql_output = '';
    $not_sql_output = " NOT" if (!(-e $self->param('sql_output')));

    print("The SQL file ".$self->param('sql_output')." that $should_be_or_has_been run to apply the suggested updates has$not_sql_output been written.\n");
  }
    return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub print_gene_transcript_info {
  my $gene = shift;
  my $trans = shift;
  my $message = shift;
  my $new_biotype = shift;

  my $g_id = $gene->dbID;
  my $g_stable_id = "undefined"; $g_stable_id = $gene->stable_id if ($gene->stable_id);
  my $g_biotype = $gene->biotype;
  my $t_id = $trans->dbID;
  my $t_stable_id = "undefined"; $t_stable_id = $trans->stable_id if ($trans->stable_id);
  my $t_biotype = $trans->biotype;

  my $str_to_print = "Gene ID: $g_id, Gene stable_id: $g_stable_id, Gene biotype: $g_biotype, Transcript ID: $t_id, Transcript stable_id: $t_stable_id, Transcript biotype: $t_biotype. $message";
  $str_to_print .= " \'$new_biotype\'" if (defined($new_biotype));
  $str_to_print .= ".\n";
  print($str_to_print);
}

sub print_gene_info {
  my $gene = shift;
  my $message = shift;
  my $new_biotype = shift;

  my $g_id = $gene->dbID;
  my $g_stable_id = "undefined"; $g_stable_id = $gene->stable_id if ($gene->stable_id);
  my $g_biotype = $gene->biotype;

  my $str_to_print = "Gene ID: $g_id, Gene stable_id: $g_stable_id, Gene biotype: $g_biotype. $message";
  $str_to_print .= " \'$new_biotype\'" if (defined($new_biotype));
  $str_to_print .= ".\n";
  print($str_to_print);
}

sub get_biotype_affix {
  my $biotype = shift;
  my $biotypes_str = shift; # biotypes regex string with no affixes
  my $prefixes_str = shift;
  my $suffixes_str = shift;
  my $prefix = '';
  my $suffix = '';

  if ($biotype =~ /^($prefixes_str)($biotypes_str)($suffixes_str)$/) {
    $prefix = $1;
    $suffix = $3;
  }
  return $prefix,$suffix;
}

sub get_biotype_basename {
  my $biotype = shift;
  my $biotypes_str = shift; # biotypes regex string with no affixes
  my $prefixes_str = shift;
  my $suffixes_str = shift;
  my $biotype_basename = '';

  if ($biotype =~ /^($prefixes_str)($biotypes_str)($suffixes_str)$/) {
    $biotype_basename = $2;
  }
  return $biotype_basename;
}

sub update_gt {
  my $gt_str = shift; 			# it will be string "gene" or "transcript"
  my $gta = shift; 			# gene/transcript adaptor
  my $gt = shift;			# gene/transcript object
  my $new_suggested_biotype = shift;	# new biotype
  my $write = shift;			# update db or not
  my $sql_output = shift;		# output file for SQL commands

  my $gt_id = $gt->dbID;
  my $gt_stable_id = $gt->stable_id;

  if ($write) {
    $gt->biotype($new_suggested_biotype);
    $gta->update($gt);

    print("$gt_str $gt_id biotype has been updated to \'$new_suggested_biotype\'.\n");
  }

  if ($sql_output) {
    open(SQL_FILE,">>$sql_output");
    if ($gt_stable_id) {
      print SQL_FILE "update $gt_str set biotype='$new_suggested_biotype' where stable_id='$gt_stable_id';\n";
    }
    else {
      print SQL_FILE "update $gt_str set biotype='$new_suggested_biotype' where $gt_str"."_id=$gt_id;\n";
    }
    close(SQL_FILE);
  }
}

sub gene_has_assembly_error_attribute {

  my $gene = shift;

  my @attribs = @{$gene->get_all_Attributes('NoTransRefError')};
  return (scalar(@attribs) > 0);
}

1;
