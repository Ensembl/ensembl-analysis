#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2021] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Digest::MD5 qw(md5);

use Bio::EnsEMBL::Utils::Exception qw(warning throw);

use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(cluster_Genes_by_coding_exon_overlap make_types_hash);
use Bio::EnsEMBL::DBSQL::DBAdaptor;


# Connection to the target DB
my $old_host;
my $old_port = '3306';
my $old_user;
my $old_pass;
my $old_dbname;
my $havana_host;
my $havana_port = '3306';
my $havana_user;
my $havana_pass;
my $havana_dbname;
my $new_host;
my $new_port = '3306';
my $new_user;
my $new_pass;
my $new_dbname;
my $iid;
my $genome_file;
my $write;

&GetOptions (
            'oh|old_host=s'   => \$old_host,
            'oP|old_port=s'   => \$old_port,
            'ou|old_user=s'   => \$old_user,
            'op|old_pass=s'   => \$old_pass,
            'od|old_dbname=s' => \$old_dbname,
            'hh|havana_host=s'   => \$havana_host,
            'hP|havana_port=s'   => \$havana_port,
            'hu|havana_user=s'   => \$havana_user,
            'hp|havana_pass=s'   => \$havana_pass,
            'hd|havana_dbname=s' => \$havana_dbname,
            'nh|new_host=s'   => \$new_host,
            'nP|new_port=s'   => \$new_port,
            'nu|new_user=s'   => \$new_user,
            'np|new_pass=s'   => \$new_pass,
            'nd|new_dbname=s' => \$new_dbname,
            'gf|genome_file=s' => \$genome_file,
            'iid=s'           => \$iid,
            'write!'          => \$write,
        );

my $old_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
            '-host'   => $old_host,
            '-port'   => $old_port,
            '-user'   => $old_user,
            '-pass'   => $old_pass,
            '-dbname' => $old_dbname,
);

my $new_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
            '-host'   => $new_host,
            '-port'   => $new_port,
            '-user'   => $new_user,
            '-pass'   => $new_pass,
            '-dbname' => $new_dbname,
);

#throw("You need a region name not '$iid'") unless ($iid);

print scalar(localtime), "\n";
if ($genome_file) {
  setup_fasta(-FASTA => $genome_file);
}
# We need to get the last id for gene, exon, transcript and translation
my $highest_stable_ids = undef;
if (!$havana_dbname) {
  $highest_stable_ids = {
    gene        => get_highest_stable_id('gene', $old_db, $new_db),
    transcript  => get_highest_stable_id('transcript', $old_db, $new_db),
    exon        => get_highest_stable_id('exon', $old_db, $new_db),
    translation => get_highest_stable_id('translation', $old_db, $new_db),
  };
}

my $stable_id_prefix = $old_db->get_MetaContainer->single_value_by_key('species.stable_id_prefix');
my $unique_old_genes = {};
my $unique_old_transcripts = {};
my $unique_old_translations = {};
my $unique_old_exons = {};
my $unique_new_genes = {};
my $unique_new_transcripts = {};
my $unique_new_translations = {};
my $unique_new_exons = {};

#foreach my $old_slice (@{$old_db->get_SliceAdaptor->fetch_all('toplevel')}) {
my $old_slice = $old_db->get_SliceAdaptor->fetch_by_name($iid);
print scalar(localtime), ' start ', $old_slice->seq_region_name, "\n";
my $new_slice = $new_db->get_SliceAdaptor->fetch_by_name($old_slice->name);

my $old_genes = $old_slice->get_all_Genes(undef, undef, 1);
my $new_genes = $new_slice->get_all_Genes(undef, undef, 1);

my ($old_gene_keys, $old_transcript_keys) = generate_hash_keys($old_genes, $stable_id_prefix, $highest_stable_ids, $unique_old_genes, $unique_old_transcripts, $unique_old_translations, $unique_old_exons);
my ($new_gene_keys, $new_transcript_keys) = generate_hash_keys($new_genes, $stable_id_prefix, $highest_stable_ids, $unique_new_genes, $unique_new_transcripts, $unique_new_translations, $unique_new_exons);

my $stable_id_events;
my $objects_to_update;

# I will first loop through all the old genes using the unique keys, start:end:num_transcripts
# Because most of the genes should be the same, I should only have few genes to try to find if
# I am looking at the correct gene
# If I can find the same key in the set of new genes, it means than nothing or almost nothing changed
# If I cannot find the same key, I look at all the transcript keys and hopefully most of the transcripts
# will be unchanged and I can then check the stable id or the name to be sure I am looking at the same
# gene.
# For the hopefully very few remaining genes, I will need to check thoroughly each objects
my @coordinate_check;
foreach my $gene_key (keys %$old_gene_keys) {
  my $old_gene = $old_gene_keys->{$gene_key};
  my $old_transcripts = $old_gene->get_all_Transcripts;
  my $gene_needs_update = 0;
  my $gene_version_change = 0;
  if (exists $new_gene_keys->{$gene_key}) {
#   If the gene_key exists in the new database I would expect to have exactly the same gene
#   or maybe some of the transcripts will be a bit different
#   So I will loop through each transcripts to find its match in the new database
    print $new_gene_keys->{$gene_key}->stable_id_version, " FOUND\n";
    my $new_gene = $new_gene_keys->{$gene_key};
    $new_gene->{_gb_found} = 0;
    my $new_transcripts = $new_gene->get_all_Transcripts;
    my @old_missed_transcripts;
    OLD_TRANSCRIPT: foreach my $old_transcript (@$old_transcripts) {
      foreach my $new_transcript (@$new_transcripts) {
#       If the transcript_key exists in the new database, I will have the same exon-intron structure
#       with the following possible changes
#         - different stable id, if the transcript was deleted and recreated
#         - different translation if the start has been updated (stop is less likely)
#           or with a non standard amino acid
#         - different cDNA if a new attribute has been added
        if (${$old_transcript->{_gb_key}} eq ${$new_transcript->{_gb_key}}) {
          my $old_exons = $old_transcript->get_all_Exons;
          my $new_exons = $new_transcript->get_all_Exons;
          for (my $index = 0; $index < @$old_exons; $index++) {
            my $exon_needs_update = 0;
            if ($old_exons->[$index]->stable_id ne $new_exons->[$index]->stable_id) {
              if (exists $unique_old_exons->{$new_exons->[$index]->stable_id}) {
                $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exons->[$index]->stable_id}, $new_exons->[$index]));
              }
              else {
                $new_exons->[$index]->created_date('');
                $new_exons->[$index]->modified_date('');
                $new_exons->[$index]->version(1);
                $exon_needs_update = 1;
              }
              print '   EXON KEY ID ', $old_exons->[$index]->stable_id_version, ' -> ', $new_exons->[$index]->stable_id_version, "\n";
            }
            elsif ($old_exons->[$index]->version != $new_exons->[$index]->version) {
              $exon_needs_update = 1 if (object_needs_update($old_exons->[$index], $new_exons->[$index]));
              print '   EXON KEY VERSION ', $old_exons->[$index]->stable_id_version, ' -> ', $new_exons->[$index]->stable_id_version, "\n";
            }
            push(@$objects_to_update, $new_exons->[$index]) if ($exon_needs_update);
          }
          process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations);
          print ' TRANSCRIPT GENE KEY TRANSCRIPT KEY ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
          next OLD_TRANSCRIPT;
        }
#       Otherwise I'm checking the stable id, we know we should be in the same genes so the transcript stable ids should match
#       If the stable ids match, it is more likely that all stable will match however at least one the following would have happened:
#         - exon with different phase whether a correction or a real change
#         - a change in the peptide, usually an extension via the start codon position
#         - an extension of the transcript
        elsif ($old_transcript->stable_id eq $new_transcript->stable_id) {
          my $old_exons = $old_transcript->get_all_Exons;
          my $new_exons = $new_transcript->get_all_Exons;
          if ($old_transcript->spliced_seq ne $new_transcript->spliced_seq) {
            my $new_index = 0;
            foreach my $old_exon (@$old_exons) {
              for (my $index = $new_index; $index < @$new_exons; $index++) {
                my $exon_needs_update = 0;
                my $new_exon = $new_exons->[$index];
                if ($old_exon->start == $new_exon->start and $old_exon->end == $new_exon->end) {
                  if ($old_exon->stable_id ne $new_exon->stable_id) {
                    if (exists $unique_old_exons->{$new_exon->stable_id}) {
                      $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon));
                    }
                    else {
                      $new_exon->created_date('');
                      $new_exon->modified_date('');
                      $new_exon->version(1);
                      $exon_needs_update = 1;
                    }
                    print '   EXON DEEP ID ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                  }
                  elsif ($old_exon->version != $new_exon->version) {
                    $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon));
                    print '   EXON DEEP VERSION ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                  }
                  $new_index = $index+1;
                  push(@$objects_to_update, $new_exon) if ($exon_needs_update);
                  last;
                }
                elsif ($old_exon->start <= $new_exon->end and $old_exon->end >= $new_exon->start) {
                  if ($old_exon->stable_id ne $new_exon->stable_id) {
                    if (exists $unique_old_exons->{$new_exon->stable_id}) {
                      $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon, 1));
                    }
                    else {
                      $new_exon->created_date('');
                      $new_exon->modified_date('');
                      $new_exon->version(1);
                      $exon_needs_update = 1;
                    }
                    print '   EXON DEEP INC ID ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                  }
                  elsif ($old_exon->version+1 != $new_exon->version) {
                    $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon, 1));
                    print '   EXON DEEP INC VERSION ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                  }
                  $new_index = $index+1;
                  push(@$objects_to_update, $new_exon) if ($exon_needs_update);
                  last;
                }
              }
            }
          }
          else {
            for (my $index = 0; $index < @$old_exons; $index++) {
              my $exon_needs_update = 0;
              my $old_exon = $old_exons->[$index];
              my $new_exon = $new_exons->[$index];
              if ($old_exon->stable_id ne $new_exon->stable_id) {
                if (exists $unique_old_exons->{$new_exon->stable_id}) {
                  $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon));
                }
                else {
                  $new_exon->created_date('');
                  $new_exon->modified_date('');
                  $new_exon->version(1);
                  $exon_needs_update = 1;
                }
                print '   EXON ID ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
              }
              elsif ($old_exon->version != $new_exon->version) {
                $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon));
                print '   EXON VERSION ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
              }
              push(@$objects_to_update, $new_exon) if ($exon_needs_update);
            }
          }
          process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations);
          print ' TRANSCRIPT GENE KEY TRANSCRIPT STABLEID ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
          next OLD_TRANSCRIPT;
        }
      }
      push(@old_missed_transcripts, $old_transcript);
    }
#   If @old_missed_transcripts has one element or more, I could not find the match of some of my old transcripts
#   in the new database so I will look at all the transcripts left in the new gene and calculate a score for each
#   of them:
#     - exon coordinate match -> 1
#     - exon stable id match -> 1
#     - exon overlap
#       - if small difference -> .4
#       - phase update -> .2
#       - end_phase update -> .2
#   Then I'm comparing at the number of exons of the old transcript times .6, if the score is above, it is considered
#   for a possible match.
#   Finally, I take the transcript with the best score. If none if found, the old transcript is labelled as deleted.
#   All the new transcripts without a match are considered new
    if (@old_missed_transcripts) {
      find_transcript_by_exon_overlap(\@old_missed_transcripts, $new_transcripts, $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations);
      print ' TRANSCRIPT GENE KEY TRANSCRIPT MISSED ', "\n";
      foreach my $new_transcript (@$new_transcripts) {
        if (!exists $new_transcript->{_gb_found}) {
          print '  TRANSCRIPT NEW ', $new_transcript->stable_id_version, "\n";
          push(@$stable_id_events, [undef, $new_transcript]);
          if ($new_transcript->translation and exists $unique_old_translations->{$new_transcript->translation->stable_id}) {
            print 'ERROR translation part of a new transcript ', $new_transcript->translation->stable_id, "\n";
          }
        }
      }
    }
#   If the gene stable id is different, I look at the gene name if present. This is for information
    if ($old_gene->stable_id ne $new_gene->stable_id) {
      my $old_name_attrib = $old_gene->get_all_Attributes('name');
      if ($old_name_attrib and @$old_name_attrib) {
        my $new_name_attrib = $new_gene->get_all_Attributes('name');
        if ($new_name_attrib and @$new_name_attrib) {
          if ($old_name_attrib->[0]->value eq $new_name_attrib->[0]->value) {
            print '  TRANSCRIPT DIFF GENEID SAME NAME RECOVER ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
          }
          else {
            print '  TRANSCRIPT DIFF GENEID DIFF NAME RECOVER ', $old_gene->stable_id_version, ' ', $old_name_attrib->[0]->value, ' -> ', $new_gene->stable_id_version, ' ', $new_name_attrib->[0]->value, "\n";
          }
        }
        else {
          print '  TRANSCRIPT DIFF GENEID NO NEW NAME RECOVER ', $old_gene->stable_id_version, ' ', $old_name_attrib->[0]->value, ' -> ', $new_gene->stable_id_version, "\n";
        }
      }
      else {
        print '  TRANSCRIPT DIFF GENEID RECOVER ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
      }
      if (exists $unique_old_genes->{$new_gene->stable_id}) {
        $gene_needs_update = 1 if (object_needs_update($unique_old_genes->{$new_gene->stable_id}, $new_gene));
      }
      else {
        $new_gene->created_date('');
        $new_gene->modified_date('');
        $new_gene->version(1);
        $gene_needs_update = 1;
      }
      print ' GENE ID ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
      push(@$stable_id_events, [$old_gene, $new_gene]);
    }
    else {
      if ($gene_version_change) {
        if ($old_gene->version+1 != $new_gene->version) {
          print '  GENE VERSION INC ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
        }
        push(@$stable_id_events, [$old_gene, $new_gene]);
        object_needs_update($old_gene, $new_gene, 1);
      }
      else {
        if ($old_gene->version != $new_gene->version) {
          print '  GENE VERSION ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
        }
        $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene));
      }
    }
    push(@$objects_to_update, $new_gene) if ($gene_needs_update);
  }
  else {
#   If the gene_key is not found, I will first look at the old transcript and try to find them with their transcript_key
#   or their stable id in the new database
#   If none of them are found, I will need to look at all the new gene left and us the coordinates to make sure that they
#   have been deleted
    print $old_gene->stable_id_version, " NOMATCH\n";
    my $gene_version_change = 0;
    my %possible_new_genes;
    my @missed_old_transcripts;
    foreach my $old_transcript (@$old_transcripts) {
      if (exists $new_transcript_keys->{${$old_transcript->{_gb_key}}}) {
        my $new_transcript = $new_transcript_keys->{${$old_transcript->{_gb_key}}};
        $possible_new_genes{$new_transcript->{_gb_gene}} = $new_transcript->{_gb_gene};
        if ($old_gene->stable_id ne $new_transcript->{_gb_gene}->stable_id) {
          my $old_name_attrib = $old_gene->get_all_Attributes('name');
          if ($old_name_attrib and @$old_name_attrib) {
            my $new_name_attrib = $new_transcript->{_gb_gene}->get_all_Attributes('name');
            if ($new_name_attrib and @$new_name_attrib) {
              if ($old_name_attrib->[0]->value eq $new_name_attrib->[0]->value) {
                print '  TRANSCRIPT DIFF GENEID SAME NAME RECOVER ', $old_gene->stable_id_version, ' -> ', $new_transcript->{_gb_gene}->stable_id_version, "\n";
              }
              else {
                print '  TRANSCRIPT DIFF GENEID DIFF NAME RECOVER ', $old_gene->stable_id_version, ' ', $old_name_attrib->[0]->value, ' -> ', $new_transcript->{_gb_gene}->stable_id_version, ' ', $new_name_attrib->[0]->value, "\n";
              }
            }
            else {
              print '  TRANSCRIPT DIFF GENEID NO NEW NAME RECOVER ', $old_gene->stable_id_version, ' ', $old_name_attrib->[0]->value, ' -> ', $new_transcript->{_gb_gene}->stable_id_version, "\n";
            }
          }
          else {
            print '  TRANSCRIPT DIFF GENEID RECOVER ', $old_gene->stable_id_version, ' -> ', $new_transcript->{_gb_gene}->stable_id_version, "\n";
          }
        }
        my $old_exons = $old_transcript->get_all_Exons;
        my $new_exons = $new_transcript->get_all_Exons;
        for (my $index = 0; $index < @$old_exons; $index++) {
          my $exon_needs_update = 0;
          if ($old_exons->[$index]->stable_id ne $new_exons->[$index]->stable_id) {
            if (exists $unique_old_exons->{$new_exons->[$index]->stable_id}) {
              $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exons->[$index]->stable_id}, $new_exons->[$index]));
            }
            else {
              $new_exons->[$index]->created_date('');
              $new_exons->[$index]->modified_date('');
              $new_exons->[$index]->version(1);
              $exon_needs_update = 1;
            }
            print '   EXON KEY ID ', $old_exons->[$index]->stable_id_version, ' -> ', $new_exons->[$index]->stable_id_version, "\n";
          }
          elsif ($old_exons->[$index]->version != $new_exons->[$index]->version) {
            $exon_needs_update = 1 if (object_needs_update($old_exons->[$index], $new_exons->[$index]));
            print '   EXON KEY VERSION ', $old_exons->[$index]->stable_id_version, ' -> ', $new_exons->[$index]->stable_id_version, "\n";
          }
        }
        process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations);
        print ' TRANSCRIPT GENE NOKEY TRANSCRIPT KEY ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
      }
      elsif (exists $unique_new_transcripts->{$old_transcript->stable_id}) {
        my $new_transcript = $unique_new_transcripts->{$old_transcript->stable_id};
        $possible_new_genes{$new_transcript->{_gb_gene}} = $new_transcript->{_gb_gene};
        if ($old_transcript->slice->name eq $new_transcript->slice->name) {
          if (!find_transcript_by_exon_overlap([$old_transcript], [$new_transcript], $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations)) {
            print 'ERROR stable_id on a different transcript ', $old_transcript->stable_id_version, , ' ', $unique_new_transcripts->{$old_transcript->stable_id}->{_gb_gene}->stable_id, "\n";
          }
        }
        else {
          print 'ERROR stable_id on a different region ', $old_transcript->stable_id_version, ' ', $old_transcript->slice->name, ' => '. $new_transcript->slice->name, "\n";
        }
      }
      else {
        push(@missed_old_transcripts, $old_transcript);
        print '  TRANSCRIPT MISSED RECOVER ', $old_transcript->stable_id_version, "\n";
      }
    }
#   If none of the transcripts have been found, I will need to check with coordinates later before saying it has been deleted
    if (@missed_old_transcripts == @$old_transcripts) {
      print '  TRANSCRIPT FOUND NONE RECOVER ', $old_gene->stable_id_version, "\n";
      push(@coordinate_check, $old_gene);
    }
    elsif (@missed_old_transcripts == 0) {
#     If all the transcripts have been found,
#       - the gene has been merged into another if they are all in the same new gene
#         or an already existing one
#       - the gene has been broken into two genes if they are in two different genes
#       - I don't expect to have my transcripts in three different genes so die
      if (keys %possible_new_genes == 1) {
        my ($new_gene) = values %possible_new_genes;
        if ($old_gene->stable_id ne $new_gene->stable_id) {
          if (exists $unique_old_genes->{$new_gene->stable_id} and exists $unique_new_genes->{$old_gene->stable_id}) {
            print 'ERROR The old stable id is present in the new database but all its transcripts have been assigned to a different gene ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
          }
          else {
            if (exists $unique_old_genes->{$new_gene->stable_id} and !exists $unique_new_genes->{$old_gene->stable_id}) {
              print '  TRANSCRIPT MERGED RECOVER ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
            }
            elsif (@$old_transcripts != @{$new_gene->get_all_Transcripts}) {
              print '  TRANSCRIPT HAVANA DELETE RECOVER ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
            }
            else {
              print '  TRANSCRIPT ERROR RECOVER ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
            }
            $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene, 1));
            push(@$stable_id_events, [$old_gene, $new_gene]);
          }
        }
        else {
          if (@$old_transcripts != @{$new_gene->get_all_Transcripts}) {
            $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene, 1));
          }
          else {
            $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene, 1));
          }
          push(@$stable_id_events, [$old_gene, $new_gene]);
        }
      }
      elsif (keys %possible_new_genes == 2) {
        if (exists $unique_new_genes->{$old_gene->stable_id}) {
          foreach my $possible_new_gene (values %possible_new_genes) {
            if ($old_gene->stable_id ne $possible_new_gene->stable_id) {
              if (exists $unique_old_genes->{$possible_new_gene->stable_id} and !exists $unique_new_genes->{$old_gene->stable_id}) {
                print '  TRANSCRIPT MERGED RECOVER ', $old_gene->stable_id_version, ' -> ', $possible_new_gene->stable_id_version, "\n";
              }
              elsif (@$old_transcripts != @{$possible_new_gene->get_all_Transcripts}) {
                print '  TRANSCRIPT HAVANA DELETE RECOVER ', $old_gene->stable_id_version, ' -> ', $possible_new_gene->stable_id_version, "\n";
              }
              else {
                print '  TRANSCRIPT ERROR RECOVER ', $old_gene->stable_id_version, ' -> ', $possible_new_gene->stable_id_version, "\n";
              }
              $gene_needs_update = 1 if (object_needs_update($old_gene, $possible_new_gene, 1));
              push(@$stable_id_events, [$old_gene, $possible_new_gene]);
            }
            else {
              if (@$old_transcripts != @{$possible_new_gene->get_all_Transcripts}) {
                $gene_needs_update = 1 if (object_needs_update($old_gene, $possible_new_gene, 1));
              }
              else {
                $gene_needs_update = 1 if (object_needs_update($old_gene, $possible_new_gene, 1));
              }
              push(@$stable_id_events, [$old_gene, $possible_new_gene]);
            }
          }
        }
        else {
          print 'ERROR gene stable id not found in new database, unlikely event ', $old_gene->stable_id_version, "\n";
        }
      }
      else {
        print 'ERROR transcripts in more than 2 different genes ', $old_gene->stable_id, ' -> ', join(', ', grep {$_->stable_id_version} values %possible_new_genes), "\n";
      }
    }
    else {
      print '  TRANSCRIPT FOUND RECOVER ', $old_gene->stable_id_version, ' ', scalar(@$old_transcripts), ' ', scalar(@missed_old_transcripts), "\n";
      my $found_transcript = 0;
      foreach my $missed_old_transcript (@missed_old_transcripts) {
        my @possible_new_transcripts;
        foreach my $possible_new_gene (values %possible_new_genes) {
          push(@possible_new_transcripts, @{$possible_new_gene->get_all_Transcripts});
        }
        if (!find_transcript_by_exon_overlap([$missed_old_transcript], \@possible_new_transcripts, $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations)) {
          print '  TRANSCRIPT MISSING DELETED RECOVER ', $missed_old_transcript->stable_id_version, ' within ', $old_gene->stable_id_version, "\n";
          push(@$stable_id_events, [$missed_old_transcript, undef]);
        }
      }
    }
  }
}
# This is the last chance, I will first check the stable id because at one point you have to trust them
# Then I will check the coordinates of the genes which don't have _gb_found set has it means no transcript
# of theirs match an old transcript
if (@coordinate_check) {
  foreach my $old_gene (@coordinate_check) {
    my $old_transcripts = $old_gene->get_all_Transcripts;
    if (exists $unique_new_genes->{$old_gene->stable_id}) {
      my $new_gene = $unique_new_genes->{$old_gene->stable_id};
      if (exists $new_gene->{_gb_found}) {
        print 'ERROR new_gene has _gb_found but for a different gene ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, ' ', $new_gene->{_gb_found}, "\n";
      }
      else {
        $new_gene->{_gb_found} = 0;
        my $new_transcripts = $new_gene->get_all_Transcripts;
        my $gene_version_change = 0;
        my $gene_needs_update = 0;
        CHECK_TRANSCRIPT: foreach my $old_transcript (@$old_transcripts) {
          foreach my $new_transcript (@$new_transcripts) {
            if (!exists $new_transcript->{_gb_found}) {
              if ($old_transcript->stable_id eq $new_transcript->stable_id) {
                my $old_exons = $old_transcript->get_all_Exons;
                my $new_exons = $new_transcript->get_all_Exons;
                if ($old_transcript->spliced_seq ne $new_transcript->spliced_seq) {
                  my $new_index = 0;
                  foreach my $old_exon (@$old_exons) {
                    for (my $index = $new_index; $index < @$new_exons; $index++) {
                      my $exon_needs_update = 0;
                      my $new_exon = $new_exons->[$index];
                      if ($old_exon->start == $new_exon->start and $old_exon->end == $new_exon->end) {
                        if ($old_exon->stable_id ne $new_exon->stable_id) {
                          if (exists $unique_old_exons->{$new_exon->stable_id}) {
                            $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon));
                          }
                          else {
                            $new_exon->created_date('');
                            $new_exon->modified_date('');
                            $new_exon->version(1);
                            $exon_needs_update = 1;
                          }
                          print '   EXON DEEP ID ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                        }
                        elsif ($old_exon->version != $new_exon->version) {
                          $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon));
                          print '   EXON DEEP VERSION ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                        }
                        $new_index = $index+1;
                        push(@$objects_to_update, $new_exon) if ($exon_needs_update);
                        last;
                      }
                      elsif ($old_exon->start <= $new_exon->end and $old_exon->end >= $new_exon->start) {
                        if ($old_exon->stable_id ne $new_exon->stable_id) {
                          if (exists $unique_old_exons->{$new_exon->stable_id}) {
                            $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon, 1));
                          }
                          else {
                            $new_exon->created_date('');
                            $new_exon->modified_date('');
                            $new_exon->version(1);
                            $exon_needs_update = 1;
                          }
                          print '   EXON DEEP INC ID ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                        }
                        elsif ($old_exon->version+1 != $new_exon->version) {
                          $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon, 1));
                          print '   EXON DEEP INC VERSION ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                        }
                        $new_index = $index+1;
                        push(@$objects_to_update, $new_exon) if ($exon_needs_update);
                        last;
                      }
                    }
                  }
                }
                else {
                  for (my $index = 0; $index < @$old_exons; $index++) {
                    my $exon_needs_update = 0;
                    my $old_exon = $old_exons->[$index];
                    my $new_exon = $new_exons->[$index];
                    if ($old_exon->stable_id ne $new_exon->stable_id) {
                      if (exists $unique_old_exons->{$new_exon->stable_id}) {
                        $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon));
                      }
                      else {
                        $new_exon->created_date('');
                        $new_exon->modified_date('');
                        $new_exon->version(1);
                        $exon_needs_update = 1;
                      }
                      print '   EXON ID ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                    }
                    elsif ($old_exon->version != $new_exon->version) {
                      $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon));
                      print '   EXON VERSION ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
                    }
                    push(@$objects_to_update, $new_exon) if ($exon_needs_update);
                  }
                }
                process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations);
                print ' TRANSCRIPT GENE STABLEID TRANSCRIPT COORD ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
                next CHECK_TRANSCRIPT;
              }
            }
          }
          print '  TRANSCRIPT NOT FOUND CHECK ', $old_transcript->stable_id_version, "\n";
          push(@$stable_id_events, [$old_transcript, undef]);
        }
        if (@$old_transcripts != @$new_transcripts) {
          $gene_version_change = 1;
        }
        if ($gene_version_change) {
          $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene, 1));
          push(@$stable_id_events, [$old_gene, $new_gene]);
        }
        else {
          $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene));
        }
        push(@$objects_to_update, $new_gene) if ($gene_needs_update);
      }
    }
    else {
      my $found_gene;
      foreach my $new_gene (@$new_genes) {
        if (!exists $new_gene->{_gb_found}) {
          if ($old_gene->strand == $new_gene->strand and $old_gene->start <= $new_gene->end and $old_gene->end >= $new_gene->start) {
            my $old_name_attrib = $old_gene->get_all_Attributes('name');
            if ($old_name_attrib and @$old_name_attrib) {
              my $new_name_attrib = $new_gene->get_all_Attributes('name');
              if ($new_name_attrib and @$new_name_attrib) {
                if ($old_name_attrib->[0]->value eq $new_name_attrib->[0]->value) {
                  print ' GENE DIFF GENEID SAME NAME CHECK ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
                }
                else {
                  print ' GENE DIFF GENEID DIFF NAME CHECK ', $old_gene->stable_id_version, ' ', $old_name_attrib->[0]->value, ' -> ', $new_gene->stable_id_version, ' ', $new_name_attrib->[0]->value, "\n";
                }
              }
              else {
                print ' GENE DIFF GENEID NO NEW NAME CHECK ', $old_gene->stable_id_version, ' ', $old_name_attrib->[0]->value, ' -> ', $new_gene->stable_id_version, "\n";
              }
            }
            else {
              print ' GENE DIFF GENEID CHECK ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
            }
            find_transcript_by_exon_overlap($old_transcripts, $new_gene->get_all_Transcripts, $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations);
            print ' TRANSCRIPT GENE COORD TRANSCRIPT COORD ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
          }
        }
      }
      if ($found_gene) {
        push(@$stable_id_events, [$old_gene, $found_gene]);
      }
      else {
        push(@$stable_id_events, [$old_gene, undef]);
        foreach my $old_transcript (@$old_transcripts) {
          push(@$stable_id_events, [$old_transcript, undef]);
        }
      }
    }
  }
}
foreach my $gene (values %$new_gene_keys) {
  my $transcripts = $gene->get_all_Transcripts;
  if (!exists $gene->{_gb_found}) {
    if ($gene->version != 1) {
      $gene->version(1);
    }
    $gene->created_date('');
    $gene->modified_date('');
    push(@$objects_to_update, $gene);
    push(@$stable_id_events, [undef, $gene]);
    $gene->{_gb_found} = 0;
  }
  if ($gene->{_gb_found} != @$transcripts) {
    foreach my $transcript (@$transcripts) {
      if (!exists $transcript->{_gb_found}) {
        if ($transcript->version != 1) {
          $transcript->version(1);
        }
        $transcript->created_date('');
        $transcript->modified_date('');
        push(@$objects_to_update, $transcript);
        push(@$stable_id_events, [undef, $transcript]);
        if (exists $unique_old_transcripts->{$transcript->stable_id}) {
          print 'ERROR old transcript not processed ', $transcript->stable_id_version, ' ', $gene->stable_id_version, "\n";
        }
      }
    }
  }
}
print scalar(localtime), ' start ', $old_slice->seq_region_name, "\n";
#}
print scalar(localtime), "\n";
print "To update: ", scalar(@$objects_to_update), "\n";
print "stable id events: ", scalar(@$stable_id_events), "\n";

if ($write) {
# Create the mapping session
  my $old_db_schema = $old_db->get_MetaContainer->get_schema_version;
  my $old_assembly = $old_db->get_CoordSystemAdaptor->get_default_version;
  my $new_db_schema = $new_db->get_MetaContainer->get_schema_version;
  my $new_assembly = $new_db->get_CoordSystemAdaptor->get_default_version;

  my @old_db_name = split('_', $old_db->dbname);
  my @new_db_name = @old_db_name;
  if ($old_db_schema+1 < $new_db_schema) {
    $old_db_name[-2] = $new_db_schema-1;
  }
  else {
    $new_db_name[-2] = $old_db_schema+1;
  }
  if ($old_db_schema != $old_db_name[-2] or $new_db_schema != $new_db_name[-2]) {
    warning('Schema mismatch, using ', $old_db_name[-2], ' -> ', $new_db_name[-2], ' instead of ', $old_db_name[-2], ' -> ', $new_db_name[-2], "\n");
    $old_db_schema = $old_db_name[-2];
    $new_db_schema = $new_db_name[-2];
  }
  $new_db->dbc->do("INSERT IGNORE INTO mapping_session (old_db_name, old_release, old_assembly, new_db_name, new_release, new_assembly, created) VALUES('".join('_', @old_db_name)."', $old_db_schema, '$old_assembly', '".join('_', @new_db_name)."', $new_db_schema, '$new_assembly', NOW())");
  my $sth = $new_db->dbc->prepare("SELECT mapping_session_id, created FROM mapping_session WHERE old_db_name = '".join('_', @old_db_name)."' AND old_release = $old_db_schema AND old_assembly = '$old_assembly' AND new_db_name = '".join('_', @new_db_name)."' AND new_release = $new_db_schema AND new_assembly = '$new_assembly'");
  $sth->execute;
  my ($mapping_session_id, $mapping_session_date) = $sth->fetchrow_array;

# Deleting data from a possible failed run
  foreach my $table (qw(gene_archive peptide_archive stable_id_event)) {
    $new_db->dbc->do("DELETE FROM $table WHERE mapping_session_id = $mapping_session_id");
    $new_db->dbc->do("ALTER TABLE $table AUTO_INCREMENT = 1");
  }
  foreach my $object_to_update (@$objects_to_update) {
    $object_to_update->created_date($mapping_session_date) unless ($object_to_update->created_date);
    $object_to_update->modified_date($mapping_session_date) unless ($object_to_update->modified_date);
    $object_to_update->adaptor->update;
  }
  my $sth_stable_id_event = $new_db->dbc->prepare("INSERT INTO stable_id_mapping (mapping_session_id, old_stable_id, old_version, new_stable_id, new_version, type, score) VALUES($mapping_session_id, ?, ?, ?, ?, ?, ?)");
  my $sth_stable_id_peptide_event = $new_db->dbc->prepare("INSERT INTO stable_id_mapping (mapping_session_id, old_stable_id, old_version, new_stable_id, new_version, type, score) VALUES($mapping_session_id, ?, ?, ?, ?, ?, ?)");
  my $sth_gene_archive = $new_db->dbc->prepare("INSERT INTO gene_archive (mapping_session_id, gene_stable_id, gene_version, transcript_stable_id, transcript_version, translation_stable_id, translation_version, peptide_archive_id) VALUES($mapping_session_id, ?, ?, ?, ?, ?, ?, ?)");
  my $sth_peptide_archive = $new_db->dbc->prepare("INSERT INTO peptide_archive (md5_checksum, peptide_seq) VALUES(?, ?)");
  foreach my $event (@$stable_id_events) {
    my $old_object = $event->[0];
    my $new_object = $event->[1];
    my $old_translation;
    my $new_translation;
    if ($old_object) {
      $sth_stable_id_event->bind_param(1, $old_object->stable_id);
      $sth_stable_id_event->bind_param(2, $old_object->version);
      my @type = split('::', ref($old_object));
      $sth_stable_id_event->bind_param(5, lc($type[-1]));
      if ($type[-1] eq 'Transcript') {
        $sth_gene_archive->bind_param(1, $old_object->{_gb_gene}->stable_id);
        $sth_gene_archive->bind_param(2, $old_object->{_gb_gene}->version);
        $sth_gene_archive->bind_param(3, $old_object->stable_id);
        $sth_gene_archive->bind_param(4, $old_object->version);
        $old_translation = $old_object->translation;
        if ($old_translation) {
          my $seq = $old_translation->seq;
          $sth_peptide_archive->bind_param(1, md5($seq));
          $sth_peptide_archive->bind_param(2, $seq);
          $sth->execute;
          $sth_gene_archive->bind_param(5, $old_translation->stable_id);
          $sth_gene_archive->bind_param(6, $old_translation->version);
          $sth_gene_archive->bind_param(7, $sth->last_insert_id);
        }
        else {
          $sth_gene_archive->bind_param(5, undef);
          $sth_gene_archive->bind_param(6, undef);
          $sth_gene_archive->bind_param(7, undef);
        }
        $sth_gene_archive->execute;
      }
    }
    else {
      $sth_stable_id_event->bind_param(1, undef);
      $sth_stable_id_event->bind_param(2, undef);
    }
    if ($new_object) {
      $sth_stable_id_event->bind_param(3, $new_object->stable_id);
      $sth_stable_id_event->bind_param(4, $new_object->version);
      my @type = split('::', ref($new_object));
      if ($type[-1] eq 'Transcript') {
        $new_translation = $new_object->translation;
      }
      $sth_stable_id_event->bind_param(5, lc($type[-1]));
      if (exists $new_object->{_gb_score}) {
        $sth_stable_id_event->bind_param(6, $new_object->{_gb_score});
      }
      else {
        $sth_stable_id_event->bind_param(6, .99);
      }
      $sth_stable_id_event->execute;
    }
    else {
      $sth_stable_id_event->bind_param(3, undef);
      $sth_stable_id_event->bind_param(4, undef);
    }
    $sth_stable_id_event->execute;
    if ($old_translation and $new_translation) {
      if ($old_translation->stable_id_version ne $new_translation->stable_id_version) {
        $sth_stable_id_event->bind_param(1, $old_translation->stable_id);
        $sth_stable_id_event->bind_param(2, $old_translation->version);
        $sth_stable_id_event->bind_param(3, $new_translation->stable_id);
        $sth_stable_id_event->bind_param(4, $new_translation->version);
        $sth_stable_id_event->bind_param(5, 'translation');
        $sth_stable_id_event->bind_param(6, .99);
        $sth_stable_id_event->execute;
      }
    }
    elsif ($old_translation and !$new_translation) {
        $sth_stable_id_event->bind_param(1, $old_translation->stable_id);
        $sth_stable_id_event->bind_param(2, $old_translation->version);
        $sth_stable_id_event->bind_param(3, undef);
        $sth_stable_id_event->bind_param(4, undef);
        $sth_stable_id_event->bind_param(5, 'translation');
        $sth_stable_id_event->bind_param(6, 0);
        $sth_stable_id_event->execute;
    }
    elsif (!$old_translation and $new_translation) {
        $sth_stable_id_event->bind_param(1, undef);
        $sth_stable_id_event->bind_param(2, undef);
        $sth_stable_id_event->bind_param(3, $new_translation->stable_id);
        $sth_stable_id_event->bind_param(4, $new_translation->version);
        $sth_stable_id_event->bind_param(5, 'translation');
        $sth_stable_id_event->bind_param(6, 0);
        $sth_stable_id_event->execute;
    }
  }
}


sub find_transcript_by_exon_overlap {
  my ($old_transcripts, $new_transcripts, $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations) = @_;

#   Calculate a score for each new_transcripts against old_transcript where old_gene and new_gene are considered the same gene:
#     - exon coordinate match -> 1
#     - exon stable id match -> 1
#     - exon overlap
#       - if small difference -> .4
#       - phase update -> .2
#       - end_phase update -> .2
#   Then I'm comparing at the number of exons of the old transcript times .6, if the score is above, it is considered
#   for a possible match.
#   Finally, I take the transcript with the best score. If none if found, the old transcript is labelled as deleted.
  my $old_transcripts_found = 0;
  foreach my $old_transcript (@$old_transcripts) {
    my $old_exons = $old_transcript->get_all_Exons;
    my @possible_new_transcripts;
    foreach my $new_transcript (@$new_transcripts) {
      if (!exists $new_transcript->{_gb_found}) {
        if ($old_transcript->start <= $new_transcript->end and $old_transcript->end >= $new_transcript->start) {
          my $new_exons = $new_transcript->get_all_Exons;
          my $new_index = 0;
          my $score = 0;
          foreach my $old_exon (@$old_exons) {
            for (my $index = $new_index; $index < @$new_exons; $index++) {
              my $new_exon = $new_exons->[$index];
              if ($old_exon->stable_id eq $new_exon->stable_id
                  or $old_exon->start == $new_exon->start and $old_exon->end == $new_exon->end) {
                $score += 1;
                $new_index = $index+1;
                last;
              }
              elsif ($old_exon->start <= $new_exon->end and $old_exon->end >= $new_exon->start) {
                if (abs($old_exon->length-$new_exon->length) < 10) {
                  $score += .4;
                }
                if ($old_exon->phase == $new_exon->phase) {
                  $score += .2;
                }
                if ($old_exon->end_phase == $new_exon->end_phase) {
                  $score += .2;
                }
                $new_index = $index+1;
                last;
              }
            }
          }
          $score /= scalar(@$old_exons);
          if ($score >= .6
              or ($old_transcript->stable_id eq $new_transcript->stable_id and $old_transcript->{_gb_gene}->stable_id eq $new_transcript->{_gb_gene}->stable_id)) {
            $new_transcript->{_gb_score} = $score;
            push(@possible_new_transcripts, $new_transcript);
          }
        }
      }
    }
    if (@possible_new_transcripts) {
      foreach my $possible_transcript (sort {$b->{_gb_score} <=> $a->{_gb_score}} @possible_new_transcripts) {
        $possible_transcript->{_gb_found} = 1;
        if (exists $possible_transcript->{_gb_gene}->{_gb_found}) {
          $possible_transcript->{_gb_gene}->{_gb_found}++;
        }
        else {
          $possible_transcript->{_gb_gene}->{_gb_found} = 1;
        }
        print '  TRANSCRIPT RECOVER ', $old_transcript->stable_id_version, ' -> ', $possible_transcript->stable_id_version, ' ', $possible_transcript->{_gb_score}, "\n";
        $old_transcripts_found++;
        foreach my $old_exon (@$old_exons) {
          my $new_exons = $possible_transcript->get_all_Exons;
          my $new_index = 0;
          my $score = 0;
          foreach my $old_exon (@$old_exons) {
            my $exon_needs_update = 0;
            for (my $index = $new_index; $index < @$new_exons; $index++) {
              my $new_exon = $new_exons->[$index];
              if ($old_exon->stable_id eq $new_exon->stable_id) {
                if ($old_exon->start == $new_exon->start and $old_exon->end == $new_exon->end) {
                  $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon));
                }
                else {
                  $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon, 1));
                }
              }
              elsif ($old_exon->start == $new_exon->start and $old_exon->end == $new_exon->end) {
                if (exists $unique_old_exons->{$new_exons->[$index]->stable_id}) {
                  $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exons->[$index]->stable_id}, $new_exons->[$index]));
                }
                else {
                  $new_exons->[$index]->created_date('');
                  $new_exons->[$index]->modified_date('');
                  $new_exons->[$index]->version(1);
                  $exon_needs_update = 1;
                }
              }
              push(@$objects_to_update, $new_exons->[$index]) if ($exon_needs_update);
              $new_index = $index+1;
              last;
            }
          }
        }
        process_transcript($old_transcript, $possible_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations);
        last;
      }
      if (@possible_new_transcripts > 1) {
        foreach my $possible_transcript (sort {$b->{_gb_score} <=> $a->{_gb_score}} @possible_new_transcripts) {
          print '  TRANSCRIPT RECOVER LOOP ', $old_transcript->stable_id_version, ' -> ', $possible_transcript->stable_id_version, ' ', $possible_transcript->{_gb_score}, "\n";
        }
      }
    }
    else {
      print '  TRANSCRIPT RECOVER DELETED ', $old_transcript->stable_id_version, "\n";
      push(@$stable_id_events, [$old_transcript, undef]);
      if ($old_transcript->translation and exists $unique_new_translations->{$old_transcript->translation->stable_id}) {
        print 'ERROR translation part of a new transcript ', $old_transcript->translation->stable_id, "\n";
      }
    }
  }
  print "WARNING $old_transcripts_found\n";
  return $old_transcripts_found;
}

sub generate_exon_keys {
  my ($exons) = @_;

  my %exon_keys;
  foreach my $exon (@$exons) {
    $exon_keys{$exon->start.':'.$exon->end} = $exon;
  }
  return \%exon_keys;
}


sub generate_hash_keys {
  my ($genes, $stable_id_prefix, $high_stable_ids, $unique_genes, $unique_transcripts, $unique_translations, $unique_exons) = @_;

  my %gene_keys;
  my %transcript_keys;
  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts;
    foreach my $transcript (@$transcripts) {
      my $transcript_id;
      foreach my $exon (@{$transcript->get_all_Exons}) {
        $transcript_id .= join(':', $exon->start, $exon->end, $exon->phase, $exon->end_phase, '');
        if (!$exon->stable_id and $high_stable_ids) {
          $exon->stable_id(sprintf("%sE%011d", $stable_id_prefix, $high_stable_ids->{exon}++));
          $exon->version(1);
        }
        if (exists $unique_exons->{$exon->stable_id} and $unique_exons->{$exon->stable_id} != $exon) {
          warning($exon->stable_id.' is not unique');
        }
        $unique_exons->{$exon->stable_id} = $exon;
      }
      my $translation = $transcript->translation;
      if ($translation) {
        if (!$translation->stable_id and $high_stable_ids) {
          $translation->stable_id(sprintf("%sP%011d", $stable_id_prefix, $high_stable_ids->{translation}++));
          $translation->version(1);
        }
        if (exists $unique_translations->{$translation->stable_id} and $unique_translations->{$translation->stable_id} != $translation) {
          warning($translation->stable_id.' is not unique');
        }
        $unique_translations->{$translation->stable_id} = $translation;
        $transcript_id .= join(':', $transcript->cdna_coding_start, $transcript->cdna_coding_end, '');
      }
      $transcript_id .= $transcript->strand;
      $transcript_keys{$transcript_id} = $transcript;
      $transcript->{_gb_key} = \$transcript_id;
      $transcript->{_gb_gene} = $gene;
      if (!$transcript->stable_id and $high_stable_ids) {
        $transcript->stable_id(sprintf("%sT%011d", $stable_id_prefix, $high_stable_ids->{transcript}++));
        $transcript->version(1);
      }
      if (exists $unique_transcripts->{$transcript->stable_id} and $unique_transcripts->{$transcript->stable_id} != $transcript) {
        warning($transcript->stable_id.' is not unique');
      }
      $unique_transcripts->{$transcript->stable_id} = $transcript;
    }
    my $gene_key = join(':', $gene->start, $gene->end, $gene->strand, scalar(@$transcripts));
    $gene_keys{$gene_key} = $gene;
    if (!$gene->stable_id and $high_stable_ids) {
      $gene->stable_id(sprintf("%sG%011d", $stable_id_prefix, $high_stable_ids->{gene}++));
      $gene->version(1);
    }
    if (exists $unique_genes->{$gene->stable_id} and $unique_genes->{$gene->stable_id} != $gene) {
      warning($gene->stable_id.' is not unique');
    }
    $unique_genes->{$gene->stable_id} = $gene;
  }
  return \%gene_keys, \%transcript_keys;
}


sub get_highest_stable_id {
  my ($table, $old_db, $new_db) = @_;

  my $old_sth_max = $old_db->prepare("SELECT MAX(stable_id) FROM $table");
  $old_sth_max->execute;
  my ($old_high_stable_id) = $old_sth_max->fetchrow_array;
  my ($highest_id) = $old_high_stable_id =~ /(\d+)/;
  my $new_sth_max = $new_db->prepare("SELECT MAX(stable_id) FROM $table");
  $new_sth_max->execute;
  my ($new_high_stable_id) = $old_sth_max->fetchrow_array;
  $new_high_stable_id =~ /(\d+)/;
  $highest_id = $1 if ($1 > $old_high_stable_id);
  return $highest_id+1;
}


sub object_needs_update {
  my ($old_object, $new_object, $increment) = @_;

  my $object_needs_update = 0;
  my $expected_version = $old_object->version;
  $expected_version += 1 if ($increment);
  if ($new_object->created_date ne $old_object->created_date) {
    $new_object->created_date($old_object->created_date);
    $object_needs_update = 1;
  }
  if ($new_object->modified_date ne $old_object->modified_date) {
    $new_object->modified_date($old_object->modified_date);
    $object_needs_update = 1;
  }
  if ($new_object->version ne $expected_version) {
    $new_object->version($expected_version);
    $object_needs_update = 1;
  }
  return $object_needs_update;
}


sub process_translations {
  my ($old_translation, $new_translation, $objects_to_update, $unique_old_translations) = @_;

  my $transcript_stable_id_event = 0;
  if ($old_translation and $new_translation) {
    my $translation_needs_update = 0;
    if ($old_translation->seq eq $new_translation->seq) {
      if ($old_translation->stable_id ne $new_translation->stable_id) {
        if (exists $unique_old_translations->{$new_translation->stable_id}) {
          $translation_needs_update = 1 if (object_needs_update($unique_old_translations->{$new_translation->stable_id}, $new_translation));
        }
        else {
          $new_translation->created_date('');
          $new_translation->modified_date('');
          $new_translation->version(1);
          $translation_needs_update = 1;
        }
        $transcript_stable_id_event = 1;
        print '   TRANSLATION KEY ID ', $old_translation->stable_id_version, ' -> ', $new_translation->stable_id_version, "\n";
      }
      elsif ($old_translation->version != $new_translation->version) {
        $translation_needs_update = 1 if (object_needs_update($old_translation, $new_translation));
        print '   TRANSLATION KEY VERSION ', $old_translation->stable_id_version, ' -> ', $new_translation->stable_id_version, "\n";
      }
    }
    else {
      $transcript_stable_id_event = 1;
      if ($old_translation->stable_id ne $new_translation->stable_id) {
        if (exists $unique_old_translations->{$new_translation->stable_id}) {
          $translation_needs_update = 1 if (object_needs_update($unique_old_translations->{$new_translation->stable_id}, $new_translation));
        }
        else {
          $new_translation->created_date('');
          $new_translation->modified_date('');
          $new_translation->version(1);
          $translation_needs_update = 1;
        }
        print '   TRANSLATION KEY FULL ID ', $old_translation->stable_id_version, ' -> ', $new_translation->stable_id_version, "\n";
      }
      elsif ($old_translation->version+1 != $new_translation->version) {
        $translation_needs_update = 1 if (object_needs_update($old_translation, $new_translation, 1));
        print '   TRANSLATION KEY INC VERSION ', $old_translation->stable_id_version, ' -> ', $new_translation->stable_id_version, "\n";
      }
    }
    push(@$objects_to_update, $new_translation) if ($translation_needs_update);
  }
  return $transcript_stable_id_event;
}


sub process_transcript {
  my ($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations, $stable_id_events) = @_;

  my $transcript_version_change = 0;
  my $transcript_stable_id_event = 0;
  my $transcript_needs_update = 0;
  if ($old_transcript->spliced_seq ne $new_transcript->spliced_seq) {
    $transcript_version_change = 1;
  }
  $transcript_stable_id_event = 1 if (process_translations($old_transcript->translation, $new_transcript->translation, $objects_to_update, $unique_old_translations));
  if ($old_transcript->stable_id ne $new_transcript->stable_id) {
    if (exists $unique_old_transcripts->{$new_transcript->stable_id}) {
      $transcript_needs_update = 1 if (object_needs_update($unique_old_transcripts->{$new_transcript->stable_id}, $new_transcript));
    }
    else {
      $new_transcript->created_date('');
      $new_transcript->modified_date('');
      $new_transcript->version(1);
      $transcript_needs_update = 1;
    }
    $transcript_stable_id_event = 1;
    print '  TRANSCRIPT KEY ID ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
  }
  else {
    if ($transcript_version_change) {
      if ($old_transcript->version+1 != $new_transcript->version) {
        print '  TRANSCRIPT KEY INC VERSION ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
      }
      $transcript_needs_update = 1 if (object_needs_update($old_transcript, $new_transcript, 1));
    }
    else {
      if ($old_transcript->version != $new_transcript->version) {
        print '  TRANSCRIPT KEY VERSION ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
      }
      $transcript_needs_update = 1 if (object_needs_update($old_transcript, $new_transcript));
    }
  }
  push(@$objects_to_update, $new_transcript) if ($transcript_needs_update);
  push(@$stable_id_events, [$old_transcript, $new_transcript]) if ($transcript_stable_id_event);
  $new_transcript->{_gb_found} = 1;
  if (exists $new_transcript->{_gb_gene}->{_gb_found}) {
    $new_transcript->{_gb_gene}->{_gb_found}++;
  }
  else {
    $new_transcript->{_gb_gene}->{_gb_found} = 1;
  }
}
