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
use POSIX qw(strftime);

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
my $genome_file;

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

print scalar(localtime), "\n";
if ($genome_file) {
  setup_fasta(-FASTA => $genome_file);
}
my $stable_id_prefix = $old_db->get_MetaContainer->single_value_by_key('species.stable_id_prefix');

# We need to get the last id for gene, exon, transcript and translation
my $highest_stable_ids = undef;
if ($havana_dbname) {
  my $havana_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
              '-host'   => $havana_host,
              '-port'   => $havana_port,
              '-user'   => $havana_user,
              '-pass'   => $havana_pass,
              '-dbname' => $havana_dbname,
  );

  $highest_stable_ids = {
    gene        => get_highest_stable_id_from_havana('gene', $havana_db, $new_db),
    transcript  => get_highest_stable_id_from_havana('transcript', $havana_db, $new_db),
    exon        => get_highest_stable_id_from_havana('exon', $havana_db, $new_db),
    translation => get_highest_stable_id_from_havana('translation', $havana_db, $new_db),
  };
  foreach my $table (qw(gene transcript translation exon)) {
    my $sth_wrong_id = $new_db->dbc->prepare("SELECT ${table}_id, stable_id FROM $table WHERE stable_id LIKE 'OTT%'");
    $sth_wrong_id->execute;
    my $sth_update_row = $new_db->dbc->prepare("UPDATE $table SET stable_id = CONCAT('$stable_id_prefix".uc(substr($table, 0, 1))."', LPAD(?, 11, 0)), version = 1 WHERE ${table}_id = ?");
    my $sth_new_id = $new_db->dbc->prepare("SELECT stable_id FROM $table WHERE ${table}_id = ?");
    my $extra = '';
    if ($table eq 'gene' or $table eq 'transcript') {
      $extra = 'is_current = 1 AND'
    }
    while (my $row = $sth_wrong_id->fetchrow_arrayref) {
      $sth_update_row->bind_param(1, $highest_stable_ids->{$table}++);
      $sth_update_row->bind_param(2, $row->[0]);
      $sth_update_row->execute;
      $sth_new_id->bind_param(1, $row->[0]);
      $sth_new_id->execute;
      my ($new_id) = $sth_new_id->fetchrow_array;

      print "UPDATE $table SET stable_id = '$new_id', version = 1 WHERE $extra stable_id = '", $row->[1], "';\n";
    }
  }
  my $duplicate_error = 0;
  foreach my $table (qw(gene transcript translation exon)) {
    my $sth_duplicated = $new_db->dbc->prepare("SELECT stable_id FROM $table WHERE stable_id LIKE '$stable_id_prefix%' GROUP BY stable_id HAVING COUNT(stable_id) > 1");
    $sth_duplicated->execute;
    my $sth_old_genes = $old_db->dbc->prepare("SELECT DISTINCT(g.stable_id) FROM gene g LEFT JOIN transcript t ON g.gene_id = t.gene_id LEFT JOIN exon_transcript et ON t.transcript_id = et.transcript_id LEFT JOIN exon e ON et.exon_id = e.exon_id WHERE e.stable_id = ?");
    my $sth_update_row = $new_db->dbc->prepare("UPDATE gene g, transcript t, exon_transcript et, exon e SET e.stable_id = CONCAT('${stable_id_prefix}E', LPAD(?, 11, 0)), e.version = 1 WHERE g.gene_id = t.gene_id AND t.transcript_id = et.transcript_id AND et.exon_id = e.exon_id AND e.stable_id = ? AND g.stable_id != ?");
    my $sth_select_dup = $new_db->dbc->prepare("SELECT exon_id FROM exon WHERE stable_id = ?");
    my $sth_update_randomrow = $new_db->dbc->prepare("UPDATE exon SET stable_id = CONCAT('${stable_id_prefix}E', LPAD(?, 11, 0)), version = 1 WHERE exon_id = ?");
    while (my $row = $sth_duplicated->fetchrow_arrayref) {
      if ($table eq 'exon') {
        my $stable_id = $row->[0];
        $sth_old_genes->bind_param(1, $stable_id);
        $sth_old_genes->execute;
        my $genes = $sth_old_genes->fetchall_arrayref;
        if (@$genes and @$genes == 1) {
          $sth_update_row->bind_param(1, $highest_stable_ids->{exon}++);
          $sth_update_row->bind_param(2, $stable_id);
          $sth_update_row->bind_param(3, $genes->[0]->[0]);
          $sth_update_row->execute;
        }
        elsif (@$genes and @$genes == 0) {
          $sth_select_dup->bind_param(1, $stable_id);
          $sth_select_dup->execute;
          my $count = 0;
          while (my $row_r = $sth_select_dup->fetchrow_arrayref) {
            if ($count) {
              $sth_update_randomrow->bind_param(1, $row_r->[0]);
              $sth_update_randomrow->execute;
            }
            $count++;
          }
        }
        else {
          warning("Something went wrong for $stable_id and ".join(' ', @$genes));
        }
      }
      else {
        $duplicate_error = 1;
        last;
      }
    }
  }
  if ($duplicate_error) {
    die("You have duplication in table which needs manual intervention");
  }
}
else {
  $highest_stable_ids = {
    gene        => get_highest_stable_id_from_core('gene', $old_db, $new_db),
    transcript  => get_highest_stable_id_from_core('transcript', $old_db, $new_db),
    exon        => get_highest_stable_id_from_core('exon', $old_db, $new_db),
    translation => get_highest_stable_id_from_core('translation', $old_db, $new_db),
  };
}

my $unique_old_genes = {};
my $unique_old_transcripts = {};
my $unique_old_translations = {};
my $unique_old_exons = {};
my $unique_new_genes = {};
my $unique_new_transcripts = {};
my $unique_new_translations = {};
my $unique_new_exons = {};

my $stable_id_events = [];
my $objects_to_update = [];

foreach my $old_slice (@{$old_db->get_SliceAdaptor->fetch_all('toplevel', undef, 1)}) {
#  next unless ($old_slice->seq_region_name eq '1');
  print scalar(localtime), ' start ', $old_slice->seq_region_name, "\n";
  my $new_slice = $new_db->get_SliceAdaptor->fetch_by_name($old_slice->name);

  my $old_genes = $old_slice->get_all_Genes(undef, undef, 1);
  my $new_genes = $new_slice->get_all_Genes(undef, undef, 1);

  my ($old_gene_keys, $old_transcript_keys) = generate_hash_keys($old_genes, $stable_id_prefix, $highest_stable_ids, $unique_old_genes, $unique_old_transcripts, $unique_old_translations, $unique_old_exons);
  my ($new_gene_keys, $new_transcript_keys) = generate_hash_keys($new_genes, $stable_id_prefix, $highest_stable_ids, $unique_new_genes, $unique_new_transcripts, $unique_new_translations, $unique_new_exons);

# I will first loop through all the old genes using the unique keys, start:end:num_transcripts
# Because most of the genes should be the same, I should only have few genes to try to find if
# I am looking at the correct gene
# If I can find the same key in the set of new genes, it means than nothing or almost nothing changed
# If I cannot find the same key, I look at all the transcript keys and hopefully most of the transcripts
# will be unchanged and I can then check the stable id or the name to be sure I am looking at the same
# gene.
# For the hopefully very few remaining genes, I will need to check thoroughly each objects
  my @coordinate_check;
  my @missing_old_genes;
  foreach my $old_gene (@$old_genes) {
    my $old_transcripts = $old_gene->get_all_Transcripts;
    my $gene_key = join(':', $old_gene->start, $old_gene->end, $old_gene->strand, scalar(@$old_transcripts));
    my $gene_needs_update = 0;
    my $gene_stable_id_event = 0;
    my @missing_old_transcripts;
    if (exists $unique_new_genes->{$old_gene->stable_id}) {
      my $new_gene = $unique_new_genes->{$old_gene->stable_id};
      print 'FOUND STABLE ID ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
      my $new_transcripts = $new_gene->get_all_Transcripts;
      foreach my $old_transcript (@$old_transcripts) {
        my $old_exons = $old_transcript->get_all_Exons;
        if (exists $unique_new_transcripts->{$old_transcript->stable_id}) {
          my $new_transcript = $unique_new_transcripts->{$old_transcript->stable_id};
          print ' FOUND STABLE ID ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
          if ($old_gene->stable_id ne $new_transcript->{_gb_gene}->stable_id) {
            print "WARNING STABLE ID gene switch for ", $old_transcript->stable_id_version, ' ', $old_gene->stable_id_version, ' -> ', $new_transcript->{_gb_gene}->stable_id, "\n";
          }
          process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations, $stable_id_events);
          process_exons($old_transcript, $new_transcript, $old_exons, $unique_old_exons, $objects_to_update);
        }
        elsif (exists $new_transcript_keys->{${$old_transcript->{_gb_key}}}) {
          my $new_transcript = $new_transcript_keys->{${$old_transcript->{_gb_key}}};
          print ' FOUND KEY ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
          if ($old_gene->stable_id ne $new_transcript->{_gb_gene}->stable_id) {
            print "WARNING KEY gene switch for ", $old_transcript->stable_id_version, ' ', $old_gene->stable_id_version, ' -> ', $new_transcript->{_gb_gene}->stable_id, "\n";
          }
          process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations, $stable_id_events);
          process_exons($old_transcript, $new_transcript, $old_exons, $unique_old_exons, $objects_to_update);
        }
        else {
          push(@missing_old_transcripts, $old_transcript);
        }
      }
      if (@missing_old_transcripts) {
        foreach my $missing_old_transcript (@missing_old_transcripts) {
          print ' MISSING ', $missing_old_transcript->stable_id_version, "\n";
          find_transcript_by_exon_overlap([$missing_old_transcript], $new_transcripts, $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations);
        }
      }
      if ($old_gene->start != $new_gene->start
          or $old_gene->end != $new_gene->end
          or @$old_transcripts != @$new_transcripts) {
        if (object_needs_update($old_gene, $new_gene, 1)) {
          push(@$objects_to_update, $new_gene);
        }
        push(@$stable_id_events, [$old_gene, $new_gene]);
      }
      else {
        if (object_needs_update($old_gene, $new_gene)) {
          push(@$objects_to_update, $new_gene);
        }
      }
    }
    elsif (exists $new_gene_keys->{$gene_key}) {
#   If the gene_key exists in the new database I would expect to have exactly the same gene
#   or maybe some of the transcripts will be a bit different
#   So I will loop through each transcripts to find its match in the new database
      my $new_gene = $new_gene_keys->{$gene_key};
      print 'FOUND KEY ', $old_gene->stable_id_version, ' ', $new_gene->stable_id_version, "\n";
      my $new_transcripts = $new_gene->get_all_Transcripts;
      foreach my $old_transcript (@$old_transcripts) {
        my $old_exons = $old_transcript->get_all_Exons;
        if (exists $unique_new_transcripts->{$old_transcript->stable_id}) {
          my $new_transcript = $unique_new_transcripts->{$old_transcript->stable_id};
          print ' FOUND STABLE ID ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
          if ($old_gene->stable_id ne $new_transcript->{_gb_gene}->stable_id) {
            print "WARNING STABLE ID gene switch for ", $old_transcript->stable_id_version, ' ', $old_gene->stable_id_version, ' -> ', $new_transcript->{_gb_gene}->stable_id, "\n";
          }
          process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations, $stable_id_events);
          process_exons($old_transcript, $new_transcript, $old_exons, $unique_old_exons, $objects_to_update);
        }
        elsif (exists $new_transcript_keys->{${$old_transcript->{_gb_key}}}) {
          my $new_transcript = $new_transcript_keys->{${$old_transcript->{_gb_key}}};
          print ' FOUND KEY ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
          if ($old_gene->stable_id ne $new_transcript->{_gb_gene}->stable_id) {
            print "WARNING KEY gene switch for ", $old_transcript->stable_id_version, ' ', $old_gene->stable_id_version, ' -> ', $new_transcript->{_gb_gene}->stable_id, "\n";
          }
          process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations, $stable_id_events);
          process_exons($old_transcript, $new_transcript, $old_exons, $unique_old_exons, $objects_to_update);
        }
        else {
          push(@missing_old_transcripts, $old_transcript);
        }
      }
      if (@missing_old_transcripts) {
        foreach my $missing_old_transcript (@missing_old_transcripts) {
          print ' MISSING ', $missing_old_transcript->stable_id_version, "\n";
          find_transcript_by_exon_overlap([$missing_old_transcript], $new_transcripts, $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations);
        }
      }
      if (exists $unique_old_genes->{$new_gene->stable_id}) {
        if (!exists $new_gene->{_gb_found}) {
          if ($unique_old_genes->{$new_gene->stable_id}->start != $new_gene->start
              or $unique_old_genes->{$new_gene->stable_id}->end != $new_gene->end
              or @{$unique_old_genes->{$new_gene->stable_id}->get_all_Transcripts} != @$new_transcripts) {
            if (object_needs_update($unique_old_genes->{$new_gene->stable_id}, $new_gene, 1)) {
              push(@$objects_to_update, $new_gene);
            }
          }
          else {
            if (object_needs_update($unique_old_genes->{$new_gene->stable_id}, $new_gene)) {
              push(@$objects_to_update, $new_gene);
            }
          }
        }
      }
      else {
        $new_gene->version(1);
        $new_gene->created_date('');
        $new_gene->modified_date('');
      }
      push(@$stable_id_events, [$old_gene, $new_gene]);
    }
    else {
#   If the gene_key is not found, I will first look at the old transcript and try to find them with their transcript_key
#   or their stable id in the new database
#   If none of them are found, I will need to look at all the new gene left and us the coordinates to make sure that they
#   have been deleted
      print 'NO MATCH ', $old_gene->stable_id_version, "\n";
      my %possible_new_genes;
      my @missed_old_transcripts;
      foreach my $old_transcript (@$old_transcripts) {
        my $old_exons = $old_transcript->get_all_Exons;
        if (exists $unique_new_transcripts->{$old_transcript->stable_id}) {
          my $new_transcript = $unique_new_transcripts->{$old_transcript->stable_id};
          print ' FOUND STABLE ID ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
          $possible_new_genes{$new_transcript->{_gb_gene}} = $new_transcript->{_gb_gene};
          if ($old_transcript->slice->name eq $new_transcript->slice->name) {
            if (!find_transcript_by_exon_overlap([$old_transcript], [$new_transcript], $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations)) {
              print 'WARNING stable_id on a different transcript ', $old_transcript->stable_id_version, , ' ', $unique_new_transcripts->{$old_transcript->stable_id}->{_gb_gene}->stable_id, "\n";
            }
          }
          else {
            print 'ERROR stable_id on a different region ', $old_transcript->stable_id_version, ' ', $old_transcript->slice->name, ' => '. $new_transcript->slice->name, "\n";
          }
        }
        elsif (exists $new_transcript_keys->{${$old_transcript->{_gb_key}}}) {
          my $new_transcript = $new_transcript_keys->{${$old_transcript->{_gb_key}}};
          print ' FOUND KEY ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
          $possible_new_genes{$new_transcript->{_gb_gene}} = $new_transcript->{_gb_gene};
          process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations, $stable_id_events);
          process_exons($old_transcript, $new_transcript, $old_exons, $unique_old_exons, $objects_to_update);
          print ' TRANSCRIPT GENE NOKEY TRANSCRIPT KEY ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
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
          if ($old_gene->stable_id eq $new_gene->stable_id) {
            $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene, 1));
          }
          else {
            if (exists $unique_old_genes->{$new_gene->stable_id} and !exists $unique_new_genes->{$old_gene->stable_id}) {
              print '  TRANSCRIPT MERGED RECOVER ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
              if ($unique_old_genes->{$new_gene->stable_id}->start == $new_gene->start
                  and $unique_old_genes->{$new_gene->stable_id}->end == $new_gene->end
                  and $unique_old_genes->{$new_gene->stable_id}->start == $new_gene->start) {
                $gene_needs_update = 1 if (object_needs_update($unique_old_genes->{$new_gene->stable_id}, $new_gene));
              }
              else {
                $gene_needs_update = 1 if (object_needs_update($unique_old_genes->{$new_gene->stable_id}, $new_gene, 1));
              }
            }
            elsif (@$old_transcripts != @{$new_gene->get_all_Transcripts}) {
              print '  TRANSCRIPT HAVANA DELETE RECOVER ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
              $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene, 1));
            }
            else {
              print '  TRANSCRIPT ERROR RECOVER ', $old_gene->stable_id_version, ' -> ', $new_gene->stable_id_version, "\n";
              $gene_needs_update = 1 if (object_needs_update($old_gene, $new_gene, 1));
            }
          }
          push(@$objects_to_update, $new_gene) if ($gene_needs_update);
          push(@$stable_id_events, [$old_gene, $new_gene]);
        }
        elsif (keys %possible_new_genes == 2) {
          if (exists $unique_new_genes->{$old_gene->stable_id}) {
            foreach my $possible_new_gene (values %possible_new_genes) {
              if ($old_gene->stable_id eq $possible_new_gene->stable_id) {
                $gene_needs_update = 1 if (object_needs_update($old_gene, $possible_new_gene, 1));
              }
              else {
                if (exists $unique_old_genes->{$possible_new_gene->stable_id} and !exists $unique_new_genes->{$old_gene->stable_id}) {
                  print '  TRANSCRIPT MERGED RECOVER ', $old_gene->stable_id_version, ' -> ', $possible_new_gene->stable_id_version, "\n";
                  if ($unique_old_genes->{$possible_new_gene->stable_id}->start == $possible_new_gene->start
                      and $unique_old_genes->{$possible_new_gene->stable_id}->end == $possible_new_gene->end
                      and $unique_old_genes->{$possible_new_gene->stable_id}->start == $possible_new_gene->start) {
                    $gene_needs_update = 1 if (object_needs_update($unique_old_genes->{$possible_new_gene->stable_id}, $possible_new_gene));
                  }
                  else {
                    $gene_needs_update = 1 if (object_needs_update($unique_old_genes->{$possible_new_gene->stable_id}, $possible_new_gene, 1));
                  }
                }
                elsif (@$old_transcripts != @{$possible_new_gene->get_all_Transcripts}) {
                  print '  TRANSCRIPT HAVANA DELETE RECOVER ', $old_gene->stable_id_version, ' -> ', $possible_new_gene->stable_id_version, "\n";
                  $gene_needs_update = 1 if (object_needs_update($old_gene, $possible_new_gene, 1));
                }
                else {
                  print '  TRANSCRIPT ERROR RECOVER ', $old_gene->stable_id_version, ' -> ', $possible_new_gene->stable_id_version, "\n";
                  $gene_needs_update = 1 if (object_needs_update($old_gene, $possible_new_gene, 1));
                }
              }
              push(@$objects_to_update, $possible_new_gene) if ($gene_needs_update);
              push(@$stable_id_events, [$old_gene, $possible_new_gene]);
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
          my $gene_needs_update = 0;
          CHECK_TRANSCRIPT: foreach my $old_transcript (@$old_transcripts) {
            my $old_exons = $old_transcript->get_all_Exons;
            foreach my $new_transcript (@$new_transcripts) {
              if (!exists $new_transcript->{_gb_found}) {
                if ($old_transcript->stable_id eq $new_transcript->stable_id) {
                  process_exons($old_transcript, $new_transcript, $old_exons, $unique_old_exons, $objects_to_update);
                  process_transcript($old_transcript, $new_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations, $stable_id_events);
                  print ' TRANSCRIPT GENE STABLEID TRANSCRIPT COORD ', $old_transcript->stable_id_version, ' -> ', $new_transcript->stable_id_version, "\n";
                  next CHECK_TRANSCRIPT;
                }
              }
            }
            print '  TRANSCRIPT NOT FOUND CHECK ', $old_transcript->stable_id_version, "\n";
            push(@$stable_id_events, [$old_transcript, undef]);
          }
          if (@$old_transcripts != @$new_transcripts) {
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
              if (find_transcript_by_exon_overlap($old_transcripts, $new_gene->get_all_Transcripts, $stable_id_events, $objects_to_update, $unique_new_transcripts, $unique_old_translations)) {
                $found_gene = $new_gene;
              }
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
  foreach my $gene (@$new_genes) {
    my $transcripts = $gene->get_all_Transcripts;
    if (!exists $gene->{_gb_found}) {
      $gene->version(1);
      $gene->created_date('');
      $gene->modified_date('');
      push(@$objects_to_update, $gene);
      push(@$stable_id_events, [undef, $gene]);
      $gene->{_gb_found} = 0;
    }
    if ($gene->{_gb_found} != @$transcripts) {
      foreach my $transcript (@$transcripts) {
        if (!exists $transcript->{_gb_found}) {
          $transcript->version(1);
          $transcript->created_date('');
          $transcript->modified_date('');
          push(@$objects_to_update, $transcript);
          push(@$stable_id_events, [undef, $transcript]);
          if (exists $unique_old_transcripts->{$transcript->stable_id}) {
            print 'WARNING old transcript not processed ', $transcript->stable_id_version, ' ', $gene->stable_id_version, "\n";
          }
          my $translation = $transcript->translation;
          if ($translation) {
            $translation->version(1);
            $translation->created_date('');
            $translation->modified_date('');
            push(@$objects_to_update, $translation);
          }
        }
      }
    }
  }
  print scalar(localtime), ' start ', $old_slice->seq_region_name, "\n";
}
foreach my $exon (values %$unique_new_exons) {
  my $exon_needs_update = 0;
  if (!exists $exon->{_gb_found}) {
    if (exists $unique_old_exons->{$exon->stable_id}) {
      print 'WARNING old exon not processed ', $unique_old_exons->{$exon->stable_id}->stable_id_version, ' ', $exon->stable_id_version, "\n";
      if ($unique_old_exons->{$exon->stable_id}->start == $exon->start and $unique_old_exons->{$exon->stable_id}->end == $exon->end) {
        $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$exon->stable_id}, $exon));
      }
      else {
        $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$exon->stable_id}, $exon, 1));
      }
    }
    else {
      $exon->version(1);
      $exon->created_date('');
      $exon->modified_date('');
      $exon_needs_update = 1;
    }
    push(@$objects_to_update, $exon) if ($exon_needs_update);
  }
}
print scalar(localtime), "\n";
print "To update: ", scalar(@$objects_to_update), "\n";
print "stable id events: ", scalar(@$stable_id_events), "\n";

print scalar(localtime), ' start mapping session', "\n";
# Create the mapping session
my $old_db_schema = $old_db->get_MetaContainer->get_schema_version;
my $old_assembly = $old_db->get_CoordSystemAdaptor->get_default_version;
my $new_db_schema = $new_db->get_MetaContainer->get_schema_version;
my $new_assembly = $new_db->get_CoordSystemAdaptor->get_default_version;

my $new_dbc = $new_db->dbc;

my @old_db_name = split('_', $old_db->dbc->dbname);
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
my $sth_tables = $new_dbc->prepare("SHOW TABLES");
$sth_tables->execute;
my %tables = (
  gene => 1,
  transcript => 1,
  translation => 1,
  exon => 1,
);
while ( my $row = $sth_tables->fetchrow_arrayref) {
  if ($row->[0] =~ /(\w+)_si_bak/) {
    $tables{$1} = 0;
  }
}
foreach my $table (keys %tables) {
  if ($tables{$table}) {
    $new_dbc->do("CREATE TABLE ${table}_si_bak SELECT * FROM $table");
  }
}
my $sth_select_mapping_session = $new_dbc->prepare("SELECT mapping_session_id, UNIX_TIMESTAMP(created) FROM mapping_session WHERE old_db_name = '".join('_', @old_db_name)."' AND old_release = $old_db_schema AND old_assembly = '$old_assembly' AND new_db_name = '".join('_', @new_db_name)."' AND new_release = $new_db_schema AND new_assembly = '$new_assembly'");
$sth_select_mapping_session->execute;
my ($mapping_session_id, $mapping_session_date) = $sth_select_mapping_session->fetchrow_array;
if (!$mapping_session_id or !$mapping_session_date) {
  $new_dbc->do("INSERT INTO mapping_session (old_db_name, old_release, old_assembly, new_db_name, new_release, new_assembly, created) VALUES('".join('_', @old_db_name)."', $old_db_schema, '$old_assembly', '".join('_', @new_db_name)."', $new_db_schema, '$new_assembly', NOW())");
  $sth_select_mapping_session->execute;
  ($mapping_session_id, $mapping_session_date) = $sth_select_mapping_session->fetchrow_array;
}

# Deleting data from a possible failed run
foreach my $table (qw(gene_archive stable_id_event)) {
  $new_dbc->do("DELETE FROM $table WHERE mapping_session_id = $mapping_session_id");
  $new_dbc->do("ALTER TABLE $table AUTO_INCREMENT = 1");
}
$new_dbc->do("DELETE FROM pa USING peptide_archive pa LEFT JOIN gene_archive ga ON pa.peptide_archive_id = ga.peptide_archive_id WHERE ga.peptide_archive_id IS NULL");

print scalar(localtime), ' update objects', "\n";
my $sth_exon_update = $new_dbc->prepare("UPDATE exon SET version = ?, created_date = ?, modified_date = ? WHERE exon_id = ?");
my $sth_translation_update = $new_dbc->prepare("UPDATE translation SET version = ?, created_date = ?, modified_date = ? WHERE translation_id = ?");
my $sth_transcript_update = $new_dbc->prepare("UPDATE transcript SET version = ?, created_date = ?, modified_date = ? WHERE transcript_id = ?");
my $sth_gene_update = $new_dbc->prepare("UPDATE gene SET version = ?, created_date = ?, modified_date = ? WHERE gene_id = ?");
my $sth;
foreach my $object_to_update (@$objects_to_update) {
  if ($object_to_update->isa('Bio::EnsEMBL::Exon')) {
    $sth = $sth_exon_update;
  }
  elsif ($object_to_update->isa('Bio::EnsEMBL::Transcript')) {
    $sth = $sth_transcript_update;
  }
  elsif ($object_to_update->isa('Bio::EnsEMBL::Translation')) {
    $sth = $sth_translation_update;
  }
  elsif ($object_to_update->isa('Bio::EnsEMBL::Gene')) {
    $sth = $sth_gene_update;
  }
  print $object_to_update->stable_id_version, ' ', $object_to_update->created_date, ' ', $object_to_update->modified_date, ' ==> ';
  $object_to_update->created_date($mapping_session_date) unless ($object_to_update->created_date);
  $object_to_update->modified_date($mapping_session_date) unless ($object_to_update->modified_date);
  $sth->bind_param(1, $object_to_update->version);
  $sth->bind_param(2, strftime("%F %T", localtime($object_to_update->created_date)));
  $sth->bind_param(3, strftime("%F %T", localtime($object_to_update->modified_date)));
  $sth->bind_param(4, $object_to_update->dbID);
  $sth->execute;
  print $object_to_update->stable_id_version, ' ', $object_to_update->created_date, ' ', $object_to_update->modified_date, "\n";
}
print scalar(localtime), ' disabling keys', "\n";
foreach my $table (qw(gene_archive peptide_archive stable_id_event)) {
  $new_dbc->do("ALTER TABLE $table DISABLE KEYS");
}
print scalar(localtime), ' stable id mapping session', "\n";
my $sth_stable_id_event = $new_dbc->prepare("INSERT INTO stable_id_event (mapping_session_id, old_stable_id, old_version, new_stable_id, new_version, type, score) VALUES($mapping_session_id, ?, ?, ?, ?, ?, ?)");
my $sth_gene_archive = $new_dbc->prepare("INSERT INTO gene_archive (mapping_session_id, gene_stable_id, gene_version, transcript_stable_id, transcript_version, translation_stable_id, translation_version, peptide_archive_id) VALUES($mapping_session_id, ?, ?, ?, ?, ?, ?, ?)");
my $sth_peptide_archive = $new_dbc->prepare("INSERT INTO peptide_archive (md5_checksum, peptide_seq) VALUES(?, ?)");
foreach my $event (@$stable_id_events) {
  my $old_object = $event->[0];
  my $new_object = $event->[1];
  my $old_translation;
  my $new_translation;
  my $score = 0;
  if ($old_object) {
    print 'OLD ', $old_object->stable_id_version, ' ';
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
        $sth_peptide_archive->execute;
        $sth_gene_archive->bind_param(5, $old_translation->stable_id);
        $sth_gene_archive->bind_param(6, $old_translation->version);
        $sth_gene_archive->bind_param(7, $sth_peptide_archive->last_insert_id);
      }
      else {
        $sth_gene_archive->bind_param(5, undef);
        $sth_gene_archive->bind_param(6, 0);
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
    print 'NEW ', $new_object->stable_id_version;
    $sth_stable_id_event->bind_param(3, $new_object->stable_id);
    $sth_stable_id_event->bind_param(4, $new_object->version);
    my @type = split('::', ref($new_object));
    if ($type[-1] eq 'Transcript') {
      $new_translation = $new_object->translation;
    }
    $sth_stable_id_event->bind_param(5, lc($type[-1]));
    if (exists $new_object->{_gb_score}) {
      $score = $new_object->{_gb_score};
    }
    else {
      $score = .99;
    }
  }
  else {
    $sth_stable_id_event->bind_param(3, undef);
    $sth_stable_id_event->bind_param(4, undef);
  }
  print "\n";
  $sth_stable_id_event->bind_param(6, $score);
  $sth_stable_id_event->execute unless ($old_object and $new_object and $old_object->stable_id_version eq $new_object->stable_id_version);
  if ($old_translation and $new_translation) {
    if ($old_translation->stable_id_version ne $new_translation->stable_id_version) {
      $sth_stable_id_event->bind_param(1, $old_translation->stable_id);
      $sth_stable_id_event->bind_param(2, $old_translation->version);
      $sth_stable_id_event->bind_param(3, $new_translation->stable_id);
      $sth_stable_id_event->bind_param(4, $new_translation->version);
      $sth_stable_id_event->bind_param(5, 'translation');
      $sth_stable_id_event->bind_param(6, .99);
      $sth_stable_id_event->execute;
      print 'OLD ', $old_translation->stable_id_version, ' ==> ', 'NEW ', $new_translation->stable_id_version, "\n";
    }
    else {
      print "NOT STORING FOR ", $old_translation->stable_id_version, ' and ', $new_translation->stable_id_version, "\n";
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
      print 'OLD ', $old_translation->stable_id_version, "\n";
  }
  elsif (!$old_translation and $new_translation) {
      $sth_stable_id_event->bind_param(1, undef);
      $sth_stable_id_event->bind_param(2, undef);
      $sth_stable_id_event->bind_param(3, $new_translation->stable_id);
      $sth_stable_id_event->bind_param(4, $new_translation->version);
      $sth_stable_id_event->bind_param(5, 'translation');
      $sth_stable_id_event->bind_param(6, 0);
      $sth_stable_id_event->execute;
      print 'NEW ', $new_translation->stable_id_version, "\n";
  }
}
print scalar(localtime), ' stable id mapping session done', "\n";
foreach my $table (qw(gene_archive peptide_archive stable_id_event)) {
  $new_dbc->do("ALTER TABLE $table ENABLE KEYS");
}
print scalar(localtime), ' keys enabled', "\n";


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
          if ($score >= .6 or $old_transcript->stable_id eq $new_transcript->stable_id) {
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
        process_transcript($old_transcript, $possible_transcript, $objects_to_update, $unique_old_transcripts, $unique_old_translations, $stable_id_events);
        process_exons($old_transcript, $possible_transcript, $old_exons, $unique_old_exons, $objects_to_update);
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
        print 'WARNING translation part of a new transcript ', $old_transcript->translation->stable_id, "\n";
      }
    }
  }
  print "WARNING $old_transcripts_found\n";
  return $old_transcripts_found;
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
        if (!$exon->stable_id) {
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
        if (!$translation->stable_id) {
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
      if (!$transcript->stable_id) {
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
    if (!$gene->stable_id) {
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


sub get_highest_stable_id_from_havana {
  my ($table, $havana_db) = @_;

  my $havana_sth_max = $havana_db->dbc->prepare("SELECT MAX(${table}_pool_id) FROM ${table}_stable_id_pool");
  $havana_sth_max->execute;
  my ($high_stable_id) = $havana_sth_max->fetchrow_array;
  return $high_stable_id+1;
}


sub get_highest_stable_id_from_core {
  my ($table, $old_db, $new_db) = @_;

  my $old_sth_max = $old_db->dbc->prepare("SELECT MAX(stable_id) FROM $table");
  $old_sth_max->execute;
  my ($old_high_stable_id) = $old_sth_max->fetchrow_array;
  my ($highest_id) = $old_high_stable_id =~ /(\d+)/;
  my $new_sth_max = $new_db->dbc->prepare("SELECT MAX(stable_id) FROM $table");
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
  print $old_object->stable_id_version, ': ', $new_object->stable_id, ' ', $new_object->version, ' ', $new_object->created_date, ' ', $new_object->modified_date, ' => ', "\t";
  if ($new_object->created_date ne $old_object->created_date) {
    $new_object->created_date($old_object->created_date);
    $object_needs_update = 1;
  }
  if ($increment) {
    $new_object->modified_date('');
    $object_needs_update = 1;
  }
  elsif ($new_object->modified_date ne $old_object->modified_date) {
    $new_object->modified_date($old_object->modified_date);
    $object_needs_update = 1;
  }
  if ($new_object->version ne $expected_version) {
    $new_object->version($expected_version);
    $object_needs_update = 1;
  }
  print $new_object->version, ' ', $new_object->created_date, ' ', $new_object->modified_date, ' => ', $object_needs_update, "\n";
  return $object_needs_update;
}


sub process_translations {
  my ($old_translation, $new_translation, $objects_to_update, $unique_old_translations) = @_;

  my $transcript_stable_id_event = 0;
  my $translation_needs_update = 0;
  if ($old_translation and $new_translation) {
    if ($old_translation->seq eq $new_translation->seq) {
      if ($old_translation->stable_id eq $new_translation->stable_id) {
        if ($old_translation->version != $new_translation->version) {
          print '   TRANSLATION KEY VERSION ', $old_translation->stable_id_version, ' -> ', $new_translation->stable_id_version, "\n";
        }
        $translation_needs_update = 1 if (object_needs_update($old_translation, $new_translation));
      }
      else {
        if (exists $unique_old_translations->{$new_translation->stable_id}) {
          if ($unique_old_translations->{$new_translation->stable_id}->seq eq $new_translation->seq) {
            $translation_needs_update = 1 if (object_needs_update($unique_old_translations->{$new_translation->stable_id}, $new_translation));
          }
          else {
            $translation_needs_update = 1 if (object_needs_update($unique_old_translations->{$new_translation->stable_id}, $new_translation, 1));
          }
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
    }
    else {
      $transcript_stable_id_event = 1;
      if ($old_translation->stable_id ne $new_translation->stable_id) {
        if (exists $unique_old_translations->{$new_translation->stable_id}) {
          if ($unique_old_translations->{$new_translation->stable_id}->seq eq $new_translation->seq) {
            $translation_needs_update = 1 if (object_needs_update($unique_old_translations->{$new_translation->stable_id}, $new_translation));
          }
          else {
            $translation_needs_update = 1 if (object_needs_update($unique_old_translations->{$new_translation->stable_id}, $new_translation, 1));
          }
        }
        else {
          $new_translation->created_date('');
          $new_translation->modified_date('');
          $new_translation->version(1);
          $translation_needs_update = 1;
        }
        print '   TRANSLATION KEY FULL ID ', $old_translation->stable_id_version, ' -> ', $new_translation->stable_id_version, "\n";
      }
      else {
        if ($old_translation->version+1 != $new_translation->version) {
          print '   TRANSLATION KEY INC VERSION ', $old_translation->stable_id_version, ' -> ', $new_translation->stable_id_version, "\n";
        }
        $translation_needs_update = 1 if (object_needs_update($old_translation, $new_translation, 1));
      }
    }
  }
  elsif (($old_translation and !$new_translation) or (!$old_translation and $new_translation)) {
    $transcript_stable_id_event = 1;
    if (!$old_translation) {
      $new_translation->created_date('');
      $new_translation->modified_date('');
      $new_translation->version(1);
      $translation_needs_update = 1;
    }
  }
  push(@$objects_to_update, $new_translation) if ($translation_needs_update);
  if ($transcript_stable_id_event) {
    print $old_translation ? $old_translation->stable_id_version : 'NEW', ' -> ', $new_translation ? $new_translation->stable_id_version : 'DELETED', "\n";
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
    $transcript_stable_id_event = 1;
  }
  $transcript_stable_id_event = 1 if (process_translations($old_transcript->translation, $new_transcript->translation, $objects_to_update, $unique_old_translations));
  if ($old_transcript->stable_id eq $new_transcript->stable_id) {
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
  else {
    if (exists $unique_old_transcripts->{$new_transcript->stable_id}) {
      if (!exists $new_transcript->{_gb_found}) {
        if ($unique_old_transcripts->{$new_transcript->stable_id}->spliced_seq eq $new_transcript->spliced_seq) {
          $transcript_needs_update = 1 if (object_needs_update($unique_old_transcripts->{$new_transcript->stable_id}, $new_transcript));
        }
        else {
          $transcript_needs_update = 1 if (object_needs_update($unique_old_transcripts->{$new_transcript->stable_id}, $new_transcript, 1));
        }
      }
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

sub process_exons {
  my ($old_transcript, $new_transcript, $old_exons, $unique_old_exons, $objects_to_update) = @_;

  my $new_exons = $new_transcript->get_all_Exons;
  if ($old_transcript->{_gb_key} eq $new_transcript->{_gb_key} or $old_transcript->spliced_seq eq $new_transcript->spliced_seq) {
    for (my $index = 0; $index < @$old_exons; $index++) {
      my $exon_needs_update = 0;
      if ($old_exons->[$index]->stable_id eq $new_exons->[$index]->stable_id) {
        if (object_needs_update($old_exons->[$index], $new_exons->[$index])) {
          $exon_needs_update = 1;
          print '   EXON KEY VERSION ', $old_exons->[$index]->stable_id_version, ' -> ', $new_exons->[$index]->stable_id_version, "\n"
            if ($old_exons->[$index]->version != $new_exons->[$index]->version);
        }
      }
      else {
        if (exists $unique_old_exons->{$new_exons->[$index]->stable_id}) {
          if ($unique_old_exons->{$new_exons->[$index]->stable_id}->start == $new_exons->[$index]->start and $unique_old_exons->{$new_exons->[$index]->stable_id}->end == $new_exons->[$index]->end) {
            $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exons->[$index]->stable_id}, $new_exons->[$index]));
          }
          else {
            $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exons->[$index]->stable_id}, $new_exons->[$index], 1));
          }
          print '   EXON KEY ID ', $old_exons->[$index]->stable_id_version, ' -> ', $new_exons->[$index]->stable_id_version, "\n";
        }
        else {
          $new_exons->[$index]->created_date('');
          $new_exons->[$index]->modified_date('');
          $new_exons->[$index]->version(1);
          $exon_needs_update = 1;
          print '   EXON KEY NEW ID ', $old_exons->[$index]->stable_id_version, ' -> ', $new_exons->[$index]->stable_id_version, "\n";
        }
      }
      $new_exons->[$index]->{_gb_found} = 1;
      push(@$objects_to_update, $new_exons->[$index]) if ($exon_needs_update);
    }
  }
  else {
    my $new_index = 0;
    OLD_EXON: foreach my $old_exon (@$old_exons) {
      my @possible_exons;
      my $exon_needs_update = 0;
      print 'STARTING index ', $old_exon->stable_id_version, ' ', $new_index, "\n";
      for (my $index = $new_index; $index < @$new_exons; $index++) {
        my $new_exon = $new_exons->[$index];
        if ($old_exon->stable_id eq $new_exon->stable_id) {
          if ($old_exon->start == $new_exon->start and $old_exon->end == $new_exon->end) {
            $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon));
          }
          else {
            $exon_needs_update = 1 if (object_needs_update($old_exon, $new_exon, 1));
          }
          $new_index = $index+1;
          push(@$objects_to_update, $new_exon) if ($exon_needs_update);
          $new_exon->{_gb_found} = 1;
          print '   EXON DEEP ID ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
          print 'NEW ID index ', $new_index, "\n";
          next OLD_EXON;
        }
        elsif ($old_exon->start == $new_exon->start and $old_exon->end == $new_exon->end and $old_exon->phase == $new_exon->phase and $old_exon->end_phase == $new_exon->end_phase) {
          if (exists $unique_old_exons->{$new_exon->stable_id}) {
            if ($unique_old_exons->{$new_exon->stable_id}->start == $new_exon->start and $unique_old_exons->{$new_exon->stable_id}->end == $new_exon->end) {
              $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon));
            }
            else {
              $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon, 1));
            }
          }
          else {
            $new_exon->created_date('');
            $new_exon->modified_date('');
            $new_exon->version(1);
            $exon_needs_update = 1;
          }
          print '   EXON DEEP MATCH ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
          $new_index = $index+1;
          push(@$objects_to_update, $new_exon) if ($exon_needs_update);
          $new_exon->{_gb_found} = 1;
          print 'NEW DEEP index ', $new_index, "\n";
          next OLD_EXON;
        }
        elsif ($old_exon->start <= $new_exon->end and $old_exon->end >= $new_exon->start) {
          push(@possible_exons, $index);
        }
      }
      if (@possible_exons == 1) {
        my $new_exon = $new_exons->[$possible_exons[0]];
        if (exists $unique_old_exons->{$new_exon->stable_id}) {
          print 'WARNING EXON DEEP POSSIBLE ', $unique_old_exons->{$new_exon->stable_id}->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
          if ($unique_old_exons->{$new_exon->stable_id}->start == $new_exon->start and $unique_old_exons->{$new_exon->stable_id}->end == $new_exon->end) {
            $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon));
          }
          else {
            $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon, 1));
          }
        }
        else {
          $new_exon->created_date('');
          $new_exon->modified_date('');
          $new_exon->version(1);
          $exon_needs_update = 1;
        }
        print '   EXON DEEP POSSIBLE ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
        $new_index = $possible_exons[0];
        push(@$objects_to_update, $new_exon) if ($exon_needs_update);
        $new_exon->{_gb_found} = 1;
      }
      else {
        foreach my $possible_index ( sort {abs($new_exons->[$a]->length-$old_exon->length) >= abs($new_exons->[$b]->length-$old_exon->length)} @possible_exons) {
          my $new_exon = $new_exons->[$possible_index];
          if (exists $unique_old_exons->{$new_exon->stable_id}) {
            print 'WARNING EXON DEEP POSSIBLE MULTI ', $unique_old_exons->{$new_exon->stable_id}->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
            $exon_needs_update = 1 if (object_needs_update($unique_old_exons->{$new_exon->stable_id}, $new_exon));
          }
          else {
            $new_exon->created_date('');
            $new_exon->modified_date('');
            $new_exon->version(1);
            $exon_needs_update = 1;
          }
          print '   EXON DEEP POSSIBLE MULTI ', $old_exon->stable_id_version, ' -> ', $new_exon->stable_id_version, "\n";
          $new_index = $possible_index;
          push(@$objects_to_update, $new_exon) if ($exon_needs_update);
          $new_exon->{_gb_found} = 1;
        }
      }
      print 'NEW POSSSIBLE index ', $new_index, "\n";
    }
  }
}
