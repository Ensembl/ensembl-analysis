=head1 LICENSE

 Copyright [2019] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::SpliceSiteCorrection

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::SpliceSiteCorrection->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module takes a set of input genes and a file containing intron frequency counts and then
edits exons where the corresponding intron bounrdy is not frequenctly observed. Everything,
whether edited or unedited is by default written to an output db

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::SpliceSiteCorrection;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name is_canonical_splice);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    small_intron_size       => 75,
    low_frequency           => 5,
    modified_biotype        => 'consensus',
    write_unmodified        => 1,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for intron dumping
 Returntype : None

=cut

sub fetch_input {
  my ($self) = @_;

  say "Fetching input";

  my $intron_file = $self->param_required('combined_intron_file');
  unless(open(IN,$intron_file)) {
    $self->throw("Could not open intron file on specified path. Path: ".$intron_file);
  }

  # This bit will process the hash of introns
  my $intron_hash = {};
  my $combined_intron_hash = {};
  while(<IN>) {
    my $line = $_;
    chomp $line;
    $self->process_introns($line,$intron_hash,$combined_intron_hash);
  }
  close IN;

  $self->param('intron_hash',$intron_hash);
  $self->param('combined_intron_hash',$combined_intron_hash);

  my $input_dba = $self->hrdb_get_dba($self->param('input_db'));
  my $dna_dba;
  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $input_dba->dnadb($dna_dba);
  }

  # If we are dumping introns from genes the get the genes
  my $gene_adaptor = $input_dba->get_GeneAdaptor();
  my $gene_ids = $self->input_id();
  my $genes = [];

  # Note that the original script in Ensembl common has some flagging here to reduce the set. I have
  # left it out for the the moment so that every gene gets porcessed
  foreach my $gene_id (@$gene_ids) {
    push(@$genes,$gene_adaptor->fetch_by_dbID($gene_id));
  }

  say "Fetched ".scalar(@$genes)." input genes";
  $self->param('input_genes',$genes);

  my $target_dba = $self->hrdb_get_dba($self->param('target_db'));
  if($dna_dba) {
    $target_dba->dnadb($dna_dba);
  }
  $self->hrdb_set_con($target_dba,'target_db');

  my $target_slice_adaptor = $target_dba->get_SliceAdaptor;
  $self->param('target_slice_adaptor',$target_slice_adaptor);
}


sub run {
  my ($self) = @_;

  my $genes = $self->param('input_genes');
  my $intron_hash = $self->param('intron_hash');
  my $combined_intron_hash = $self->param('combined_intron_hash');
  my $target_slice_adaptor = $self->param('target_slice_adaptor');
  my $processed_output_genes = [];
  foreach my $gene (@$genes) {
    say "Processing gene: ".$gene->dbID;
    my $modified_gene = $self->fix_introns($gene,$intron_hash,$combined_intron_hash,$target_slice_adaptor,$self->param('low_frequency'));
    say "Modified gene start/end: ".$modified_gene->start."..".$modified_gene->end;
    if($modified_gene) {
      push(@$processed_output_genes,$modified_gene);
    } elsif($self->param('write_unmodified')) {
      push(@$processed_output_genes,$gene);
    }
  }
  $self->output($processed_output_genes);
}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $target_db = $self->hrdb_get_con('target_db');
  my $gene_adaptor = $target_db->get_GeneAdaptor();

  say "Finished processing. Preparing to write genes to db";
  foreach my $output_gene (@{$self->output}) {
    empty_Gene($output_gene);
    $gene_adaptor->store($output_gene);
  }
}


# Currently not in use
sub check_input_gene {
  my ($self,$gene,$combined_intron_hash,$slice_adaptor,$low_frequency) = @_;

  my $seq_region = $gene->seq_region_name;
  my $flag = 0;
  my $transcripts = $gene->get_all_Transcripts();
  if(scalar(@$transcripts > 1)) {
    $self->throw("Gene has more than one transcript, not designed for multi transcript genes. Gene dbID: ".$gene->dbID);
  }

  foreach my $transcript (@$transcripts) {
    my $introns = $transcript->get_all_Introns();
    foreach my $intron (@$introns) {
      my ($is_canonical,$donor,$acceptor) = is_canonical_splice($intron,$slice_adaptor,$gene->slice);
      my $intron_string = $intron->start.':'.$intron->end.':'.$intron->strand;
      my $frequency = 1;
      if($combined_intron_hash->{$seq_region}->{$intron_string}) {
        $frequency = $combined_intron_hash->{$seq_region}->{$intron_string};
      }

      unless($is_canonical) {
#        say "Non-canonical intron: ".$intron_string.":".$frequency." ".$donor."..".$acceptor;
        $flag = 1;
      } elsif($frequency <= $low_frequency) {
#        say "Low frequency intron: ".$intron_string.":".$frequency." ".$donor."..".$acceptor;
        $flag = 1;
      } else {
#        say "Canonical introns: ".$intron_string.":".$frequency." ".$donor."..".$acceptor;
      }
    }
  }

  if($flag) {
    return(1);
  }

  return(0);
}


sub process_introns {
  my ($self,$line,$intron_hash,$combined_intron_hash) = @_;

  my @line = split(',',$line);
  my $header = $line[0];
  my ($seq_region,$type) = split(':',$header);

  foreach(my $i=1; $i<scalar(@line); $i++) {
    my $intron_details = $line[$i];
    my ($start,$end,$strand,$count) = split(':',$intron_details);
    my $intron_string = $start.':'.$end.':'.$strand;

    if($combined_intron_hash->{$seq_region}->{$intron_string}) {
      $combined_intron_hash->{$seq_region}->{$intron_string} += $count;
    } else {
      $combined_intron_hash->{$seq_region}->{$intron_string} = $count;
    }

    if($intron_hash->{$seq_region}->{$type}->{$intron_string}) {
      $self->warning("An entry for the intron already existed, the file is expected to have unique entries at this point. Skipping this intron entry");
      next;
    } else {
      $intron_hash->{$seq_region}->{$type}->{$intron_string} = $count;
    }
  }
}


sub fix_introns {
  my ($self,$gene,$introns_hash,$combined_intron_hash,$slice_adaptor,$low_frequency) = @_;

  say "Attempting fix for gene";
  my $small_intron_size = $self->param_required('small_intron_size');
  my $max_edit_distance = 30;
  my $min_intron_ratio = 0.1;
  my $seq_region = $gene->seq_region_name;
  my $hash_intron_types = $introns_hash->{$seq_region};

  my $priority = {'cdna'      => 5,
                  'rnaseq'    => 4,
                  'long_read' => 3,
                  'projecton' => 2,
                  'protein'   => 1};

  my $transcripts = $gene->get_all_Transcripts();

  my $new_gene = Bio::EnsEMBL::Gene->new();
  $new_gene->analysis($gene->analysis);
  $new_gene->slice($gene->slice);
  $new_gene->biotype($gene->biotype);
  $new_gene->stable_id($gene->stable_id);
  $new_gene->version($gene->version);

  my $new_transcripts = [];
  foreach my $transcript (@$transcripts) {
    my $modified = 0;
    my $introns = $transcript->get_all_Introns();
    my $intron_index = 0;
    my $exons = $transcript->get_all_Exons();
    foreach my $intron (@$introns) {
      $intron_index++;
      my $intron_start = $intron->start();
      my $intron_end = $intron->end();
      my $intron_strand = $intron->strand;
      my $intron_string = $intron_start.":".$intron_end.":".$intron_strand;
      my $intron_count = 1;
      if($combined_intron_hash->{$seq_region}->{$intron_string}) {
        $intron_count = $combined_intron_hash->{$seq_region}->{$intron_string};
      }

      my ($is_canonical,$donor,$acceptor) = is_canonical_splice($intron,$slice_adaptor,$gene->slice);

#
#      if($is_canonical && $intron_count > $low_frequency) {
#        say "SKIPPING: ".$intron_string.":".$intron_count;
#        next;
#      }

      my $exon_left;
      my $exon_right;
      if($intron_strand == 1) {
        $exon_left = $$exons[$intron_index-1];
        $exon_right = $$exons[$intron_index];
      } else {
        $exon_left = $$exons[$intron_index];
        $exon_right = $$exons[$intron_index-1];
      }

      my $best_introns = {};
      foreach my $type_key (keys(%$hash_intron_types)) {
        say "TK: ".$type_key;
        my $hash_introns = $hash_intron_types->{$type_key};
        foreach my $hash_intron (keys(%$hash_introns)) {
          if($intron_string eq $hash_intron) {
            next;
          }

          my $hash_intron_count = $hash_introns->{$hash_intron};

          my ($hash_intron_start,$hash_intron_end,$hash_intron_strand) = split(':',$hash_intron);
          unless($intron_strand == $hash_intron_strand) {
            next;
          }
          unless((abs($hash_intron_start - $intron_start) <= $max_edit_distance) && (abs($hash_intron_end - $intron_end) <= $max_edit_distance)) {
            next;
          }

          if (!(exists $best_introns->{$type_key}) || $hash_intron_count >= ${$best_introns->{$type_key}}[1]) {
            $best_introns->{$type_key} = [$hash_intron,$hash_intron_count];
          }
        } # foreach my $hash_intron (keys(%$hash_introns))
      } # foreach my $type_key

      say "Checking for best potential intron";
      my $final_type;
      foreach my $best_type (keys(%$best_introns)) {
        say "Best type: ".$best_type;
        my ($hash_intron_start,$hash_intron_end,$hash_intron_strand) = split(':',${$best_introns->{$best_type}}[0]);
        my $hash_intron_count = ${$best_introns->{$best_type}}[1];
        say "Original intron: ".$intron_start."..".$intron_end.":".$intron_strand.":".$intron_count;
        say "Potential edit:  ".$hash_intron_start."..".$hash_intron_end.":".$hash_intron_strand.":".$hash_intron_count;
        unless($final_type) {
          $final_type=$best_type;
        } else {
          if($priority->{$best_type} > $priority->{$final_type}) {
            $final_type = $best_type;
          }
        }
      }

      unless($final_type) {
        say "No suitable replacement introns found";
        next;
      }

      my $final_intron = ${$best_introns->{$final_type}}[0];
      my $final_intron_count = ${$best_introns->{$final_type}}[1];
      say "Final selected intron: ".$final_intron;
      # These are currently the scenarios we might want to edit on
      # The frequency ratio should be the most useful in general as a lot of incorrect introns pass the other criteria
      if(($intron_count/$final_intron_count) < $min_intron_ratio || $intron_count <= $low_frequency || !$is_canonical || $intron->length <= $small_intron_size) {
        if($self->update_exon_boundaries($exon_left,$exon_right,$final_intron)) {
          $modified = 1;
        }
      } else {
        say "Intron passed checks, so will not modify";
      }
    } # foreach my $intron (@$introns)

    unless($modified) {
      say "The transcript was not modified";
      push(@$new_transcripts,$transcript);
    } else {
      say "The transcript was modified";
      foreach my $exon (@$exons) {
        $exon->phase(-1);
        $exon->end_phase(-1);
      }

      my $new_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $exons);
      $new_transcript->slice($transcript->slice);
      $new_transcript->analysis($transcript->analysis);
      $new_transcript->biotype($self->param_required('modified_biotype'));
      $new_transcript->stable_id($transcript->stable_id);
      $new_transcript->version($transcript->version);
      compute_translation($new_transcript);
      push(@$new_transcripts,$new_transcript);
    }
  } # end foreach my $transcript (@$transcripts)

  foreach my $new_transcript (@$new_transcripts) {
    $new_gene->add_Transcript($new_transcript);

    # Note this will modify the gene biotype if the original transcript biotype matches the modified biotype value
    # This could cause some confusion if the gene biotype differed from the transcript biotype, but if this were to
    # happen the confusion begins before this
    if($new_transcript->biotype eq $self->param_required('modified_biotype')) {
      $new_gene->biotype($self->param_required('modified_biotype'));
    }
  }

  return($new_gene);
}


sub update_exon_boundaries {
  my ($self,$exon_left,$exon_right,$final_intron) = @_;

  say "Updating exon boundaries";
 # print "Original exon left: ";
 # $self->print_exon($exon_left);
 # print "Original exon right: ";
 # $self->print_exon($exon_right);

  # Now update the exon boundaries
  my ($intron_start,$intron_end,$intron_strand) = split(':',$final_intron);
#  say "IL: ".($intron_start-1);
#  say "IR: ".($intron_end+1);
  my $diff_left = ($intron_start-1) - $exon_left->seq_region_end;
  my $diff_right = ($intron_end+1) - $exon_right->seq_region_start;
#  say "DL: ".$diff_left;
#  say "DR: ".$diff_right;

  my $slice_shift_left = $exon_left->end + $diff_left;
  my $slice_shift_right = $exon_right->start + $diff_right;

  unless($exon_left->start < $slice_shift_left && $slice_shift_right < $exon_right->end) {
    say "Issue with adjusted boundaries, will not modify";
    $self->warning("New boundaries would produce an exon with start >= end (or vice versa). Will not modify");
    return 0;
  }
  $exon_left->end($slice_shift_left);
  $exon_right->start($slice_shift_right);

  return 1;
#  print "Modified exon left: ";
#  $self->print_exon($exon_left);
#  print "Modified exon right: ";
#  $self->print_exon($exon_right);

}


sub print_exon {
  my ($self,$exon) = @_;

  my $exon_string = "(".$exon->seq_region_start."..".$exon->seq_region_end."):".$exon->strand;
  say $exon_string;
}


1;
