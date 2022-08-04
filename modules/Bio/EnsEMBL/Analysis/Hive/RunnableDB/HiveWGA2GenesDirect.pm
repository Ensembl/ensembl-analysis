#!/usr/bin/env perl

=head1 LICENSE

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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::HiveWGA2GenesDirect -

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $genscan = Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect
    ->new (-db         => $pipelinedb,
           -input_id   => $input_id
           -analysis   => $analysis );
  $genscan->fetch_input();
  $genscan->run();
  $genscan->write_output(); #writes to stdout


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect;

use warnings ;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(attach_Analysis_to_Gene attach_Slice_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(empty_Transcript);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(parse_timer create_file_name);

use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(replace_stops_with_introns);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my($self) = @_;
  $self->create_analysis;
  my $input_ids = $self->param('iid');

  my $max_internal_stops = $self->param('max_internal_stops');
  unless(defined($max_internal_stops)) {
    $self->warning("No max_internal_stops param found, defaulting to 1");
    $max_internal_stops = 1;
  }
  $self->max_internal_stops($max_internal_stops);

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');

  # Define the source transcript and target transcript dbs
  my $source_transcript_dba = $self->hrdb_get_dba($self->param('source_transcript_db'));
  my $target_transcript_dba = $self->hrdb_get_dba($self->param('target_transcript_db'));
  $target_transcript_dba->dnadb($target_dna_dba);
  $self->hrdb_set_con($source_transcript_dba,'source_transcript_db');
  $self->hrdb_set_con($target_transcript_dba,'target_transcript_db');

  # Define the compara db
  my $compara_dba = $self->hrdb_get_dba($self->param('compara_db'),undef,'Compara');
  $self->hrdb_set_con($compara_dba,'compara_db');

  # Get the genome db adpator
  my $genome_dba = $compara_dba->get_GenomeDBAdaptor;

  # Retrieve the production names for the query and target species
  my $source_species = $source_transcript_dba->get_MetaContainerAdaptor->get_production_name();
  my $target_species = $target_transcript_dba->get_MetaContainerAdaptor->get_production_name();

  my $source_genome_db = $genome_dba->fetch_by_core_DBAdaptor($source_transcript_dba);
  my $target_genome_db = $genome_dba->fetch_by_core_DBAdaptor($target_transcript_dba);

  ########
  # check that the default assembly for the query and target agrees
  # with that for the method_link_species_set GenomeDBs
  ########

  my $source_assembly = $source_genome_db->assembly;
  my $target_assembly = $target_genome_db->assembly;

  my ($source_assembly_version, $target_assembly_version);
  eval {
    $source_assembly_version = $source_transcript_dba->get_CoordSystemAdaptor->fetch_by_name('toplevel',$source_genome_db->assembly);
    $target_assembly_version = $target_transcript_dba->get_CoordSystemAdaptor->fetch_by_name('toplevel',$target_genome_db->assembly);
  };
  if($@) {
    $self->throw("Had trouble fetching coord systems for ".$source_genome_db->assembly . " and " .$target_genome_db->assembly . " from core dbs:\n".$@);
  }

  my $transcript_align_slices = {};
  my $mlss = $compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_by_method_link_type_GenomeDBs($self->param('method_link_type'),
                                                                                                [$source_genome_db,
                                                                                                $target_genome_db]);
  unless($mlss) {
    $self->throw("No MethodLinkSpeciesSet for :\n" .$self->param('method_link_type') . "\n" .$source_species . "\n" .$target_species);
  }

  foreach my $input_id (@$input_ids) {
    my $transcript = $source_transcript_dba->get_TranscriptAdaptor->fetch_by_dbID($input_id);
    
    if ($self->param('canonical') and !($transcript->is_canonical())) {
      say "Skipping transcript because it is not canonical and the parameter 'canonical' is 1.";
      next;
    } else {
      say "Processing transcript: ".$transcript->dbID;
    }

    #########
    # get the compara data: MethodLinkSpeciesSet, reference DnaFrag,
    # and all GenomicAlignBlocks
    #########
    my $dnafrag = $compara_dba->get_DnaFragAdaptor->fetch_by_GenomeDB_and_name($source_genome_db,$transcript->slice->seq_region_name);
    my $genomic_align_block_adaptor = $compara_dba->get_GenomicAlignBlockAdaptor;
    my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,$dnafrag,$transcript->start,$transcript->end);

    my (%chains, @chains);
    foreach my $genomic_align_block (@$genomic_align_blocks) {
      my $source_genomic_align = $genomic_align_block->reference_genomic_align;
      my ($target_genomic_align) = @{$genomic_align_block->get_all_non_reference_genomic_aligns};

      # fetch the target slice for later reference
      if (not exists $transcript->{'_target_slices'}->{$target_genomic_align->dnafrag->name}) {
        $transcript->{'_target_slices'}->{$target_genomic_align->dnafrag->name} = $target_transcript_dba->get_SliceAdaptor->fetch_by_region('toplevel',$target_genomic_align->dnafrag->name);
      }

      if ($genomic_align_block->reference_genomic_align->dnafrag_strand < 0) {
        $genomic_align_block->reverse_complement;
      }

      push @{$chains{$genomic_align_block->group_id}}, $genomic_align_block;
    }

    foreach my $chain_id (keys %chains) {
      push @chains, [
                    sort {
                       $a->reference_genomic_align->dnafrag_start <=> $b->reference_genomic_align->dnafrag_start;
                     } @{$chains{$chain_id}}
                    ];
    }

    unless (scalar@chains) {
      say "Transcript ".$transcript->dbID." has no alignment chains, skipping.";
      next;
    }
    $transcript->{'_genomic_align_block_chains'} = \@chains;
    $self->source_transcripts($transcript);
 }

}

sub run {
  my ($self) = @_;

  my @result_transcripts;
  my $timer = $self->param('timer');
  $timer = parse_timer($timer);
  my $failed_transcripts = [];

  foreach my $source_transcript (@{$self->source_transcripts}) {
     my $successfully_projected = 0;
     my $preliminary_transcripts = [];
     eval {
      local $SIG{ALRM} = sub { die "alarm clock restart\n" };
      alarm $timer; #schedule alarm in '$timer' seconds

      foreach my $chain (@{$source_transcript->{_genomic_align_block_chains}}) {
        my $gene_scaffold = Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold->new(
                                                                                     -genomic_align_blocks => $chain,
                                                                                     -from_slice    => $source_transcript->slice,
                                                                                     -to_slices     => $source_transcript->{'_target_slices'},
                                                                                     -transcripts   => [$source_transcript],
                                                                                     -max_readthrough_dist => $self->param('max_exon_readthrough_dist'),
                                                                                     -direct_target_coords => 1,
                                                                                   );
        my $proj_transcript;
        eval{
          $proj_transcript = $gene_scaffold->place_transcript($source_transcript);
        };
        if($@ || !$proj_transcript) {
          say "Projection failed on block chain for transcript ".$source_transcript->dbID;
          next;
        }

        $proj_transcript = $self->process_transcript($proj_transcript);
        unless($proj_transcript) {
          next;
        }

        if($self->param('calculate_coverage_and_pid')) {
          if($proj_transcript->translation->length < 20000 && $source_transcript->translation->length < 20000) {
            $self->realign_translation($source_transcript,$proj_transcript);
          }
          else {
            $self->warning('Not realigning translation as translation length >= 20000');
          }
        }

        my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_t', -VALUE => $source_transcript->stable_id.".".$source_transcript->version);
        $proj_transcript->add_Attributes($parent_attribute);
        $proj_transcript->stable_id($source_transcript->stable_id);
        push(@{$preliminary_transcripts}, $proj_transcript);
        $successfully_projected = 1;
      } # end foreach my $chain
      alarm 0; #reset alarm
    }; # end eval

    if($@ && $@ eq "alarm clock restart\n") {
      say "Projection failed for transcript ".$source_transcript->dbID." because of time limit on timer param";
    }

    unless($successfully_projected) {
      push(@{$failed_transcripts}, $source_transcript->dbID);
      next;
    }

    # Pick the best transcript and any others that are within 5 percent of the combined coverage and identity of the best
    my $selected_transcripts;
    if(scalar(@{$preliminary_transcripts}) == 1) {
      $selected_transcripts = $preliminary_transcripts;
    } else {
      $selected_transcripts = $self->select_best_transcripts($preliminary_transcripts);
    }

    foreach my $selected_transcript (@{$selected_transcripts}) {
      push(@result_transcripts, $selected_transcript);
    }

  } # end foreach my $source_transcript

  if ($self->TRANSCRIPT_FILTER){
    $self->output($self->filter->filter_results(\@result_transcripts));
  }
  else {
    $self->output(\@result_transcripts);
  }

  $self->param('_failed', $failed_transcripts);
  say "At the end of RUN, we had ", scalar(@{$self->output}), " transcripts";
}

sub write_output {
  my ($self) = @_;

  my $transcript_count = 0;
  my $target_transcript_dbc = $self->hrdb_get_con('target_transcript_db');
  my $target_gene_adaptor = $target_transcript_dbc->get_GeneAdaptor();
  my $failure_branch_code = $self->param('_branch_to_flow_to_on_fail');

  foreach my $transcript (@{$self->output}) {
    $transcript->analysis($self->analysis);
    $transcript->biotype('projection');
    empty_Transcript($transcript);

    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->biotype('projection');

    $gene->add_Transcript($transcript);
    attach_Slice_to_Gene($gene, $transcript->slice);
    attach_Analysis_to_Gene($gene, $self->analysis);
    $target_gene_adaptor->store($gene);
    $transcript_count++;
  }

  my $output_hash = {};
  my $failed_transcript_ids = $self->param('_failed');
  foreach my $failed_transcript_id (@$failed_transcript_ids) {
    $output_hash->{'iid'} = [$failed_transcript_id];
    $self->dataflow_output_id($output_hash, $failure_branch_code);
  }

  return 1;

}


sub runnable_failed {
  my ($self,$runnable_failed) = @_;
  unless ($self->param_is_defined('_runnable_failed')) {
    $self->param('_runnable_failed',[]);
  }
  if ($runnable_failed) {
    push (@{$self->param('_runnable_failed')},$runnable_failed);
  }
  return ($self->param('_runnable_failed'));
}

sub source_transcripts {
  my ($self,$transcript) = @_;
  unless ($self->param_is_defined('_source_transcripts')) {
    $self->param('_source_transcripts',[]);
  }

  if ($transcript) {
    push (@{$self->param('_source_transcripts')},$transcript);
  }
  return ($self->param('_source_transcripts'));
}

sub realign_translation {
  my ($self,$source_transcript,$projected_transcript) = @_;

  my $query_seq = $source_transcript->translate->seq();
  my $projected_seq = $projected_transcript->translate->seq();

  my $align_input_file = create_file_name('projected_align_', 'fa');
  my $align_output_file = create_file_name('projected_align_', 'aln');

  open(INPUT,">".$align_input_file);
  say INPUT ">query";
  say INPUT $query_seq;
  say INPUT ">target";
  say INPUT $projected_seq;
  close INPUT;

  my $align_program_path = 'muscle';

  my $cmd = $align_program_path." -in ".$align_input_file." -out ".$align_output_file;
  my $result = system($cmd);

  if($result) {
    throw("Got a non-zero exit code from alignment. Commandline used:\n".$cmd);
  }

  my $file = "";
  open(ALIGN,$align_output_file);
  while(<ALIGN>) {
    $file .= $_;
  }
  close ALIGN;

  unless($file =~ /\>.+\n(([^>]+\n)+)\>.+\n(([^>]+\n)+)/) {
    throw("Could not parse the alignment file for the alignment sequences. Alignment file: ".$align_output_file);
  }

  my $aligned_query_seq = $1;
  my $aligned_projected_seq = $3;

  $aligned_query_seq =~ s/\n//g;
  $aligned_projected_seq =~ s/\n//g;

  `rm $align_input_file`;
  `rm $align_output_file`;

  # Work out coverage
  my $coverage;
  my $temp = $aligned_projected_seq;
  my $projected_gap_count = $temp =~ s/\-//g;
  my $ungapped_query_seq = $aligned_query_seq;
  $ungapped_query_seq  =~ s/\-//g;

  if(length($ungapped_query_seq) == 0) {
    $coverage = 0;
  } else {
    $coverage = 100 - (($projected_gap_count/length($ungapped_query_seq)) * 100);
  }

  # Work out precent identity
  my $match_count = 0;
  my $aligned_positions = 0;
  for(my $j=0; $j<length($aligned_query_seq); $j++) {
    my $char_query = substr($aligned_query_seq,$j,1);
    my $char_target = substr($aligned_projected_seq,$j,1);
    if($char_query eq '-' || $char_target  eq '-') {
      next;
    }
    if($char_query eq $char_target) {
      $match_count++;
    }
    $aligned_positions++;
  }

  unless($aligned_positions) {
    throw("Pairwise alignment between the query sequence and the translation shows zero aligned positions. Something has gone wrong");
  }

  my $percent_id = ($match_count / $aligned_positions) * 100;

  # Get all exons and transcript supporting features
  my $transcript_supporting_features = $projected_transcript->get_all_supporting_features();
  my $exons = $projected_transcript->get_all_Exons();

  # Now clean these out
  $projected_transcript->flush_Exons();
  $projected_transcript->flush_supporting_features();

  # Loop through the TSFs and add the coverage and pid, then add back into transcript
  foreach my $transcript_supporting_feature (@{$transcript_supporting_features}) {
    $transcript_supporting_feature->hcoverage($coverage);
    $transcript_supporting_feature->percent_id($percent_id);
    $transcript_supporting_feature->hseqname($source_transcript->stable_id.".".$source_transcript->version);
    $projected_transcript->add_supporting_features($transcript_supporting_feature);
  }


  # Loop through exons, get supporting features for each, flush existing SF, add coverage and pid, add back to exon, add exon to transcript
  foreach my $exon (@{$exons}) {
    my $exon_supporting_features = $exon->get_all_supporting_features();
    $exon->flush_supporting_features();
    foreach my $exon_supporting_feature (@{$exon_supporting_features}) {
      $exon_supporting_feature->hcoverage($coverage);
      $exon_supporting_feature->percent_id($percent_id);
      $exon_supporting_feature->hseqname($source_transcript->stable_id.".".$source_transcript->version);
      $exon->add_supporting_features($exon_supporting_feature);
    }
    $projected_transcript->add_Exon($exon);
  }

}

sub select_best_transcripts {
  my ($self,$preliminary_transcripts) = @_;

  my $selected_transcripts = [];
  my $best_score = 0;
  foreach my $preliminary_transcript (@{$preliminary_transcripts}) {
    my @tsfs = @{$preliminary_transcript->get_all_supporting_features};
    my $cov = $tsfs[0]->hcoverage;
    my $pid = $tsfs[0]->percent_id;
    my $combined_score = $cov + $pid;
    if($combined_score > $best_score) {
      $best_score = $combined_score;
    }
    $preliminary_transcript->{'combined_score'} = $combined_score;
  }

  unless($best_score) {
    $self->throw('Issue with calculating the best score from projection result set. This should not happen.');
  }

  # At this point we know the highest value in terms of combined percent identity and coverage. Now
  # we want all transcripts that fall within 5 percent of this value
  foreach my $preliminary_transcript (@{$preliminary_transcripts}) {
    my $combined_score = $preliminary_transcript->{'combined_score'};
    #if($combined_score/$best_score >= 0.95) {
    if ($combined_score/$best_score >= 1 and scalar(@{$selected_transcripts}) == 0) { 
      push(@{$selected_transcripts},$preliminary_transcript);
    }
  }

  unless(scalar(@{$selected_transcripts})) {
    $self->throw("No transcripts selected in select_best_transcripts, something went wrong");
  }

  return($selected_transcripts);

}

#################################################################
# FUNCTION  : process_transcript
#
# Description:
#    Subjects the given transcript to a few tests, returning 
#    the transcript if they succeed, undef if not. If the
#    transcript contains less than $max_stops stops, these
#    are "spliced out"; otherwise the transcripts is rejected
#################################################################

sub process_transcript {
  my ($self,$tran) = @_;

  return 0 if not $tran;

  my (@processed_transcripts);

  my @exons = @{$tran->get_all_Exons};
  my ($tsf) = @{$tran->get_all_supporting_features};
  my $pep = $tran->translate->seq;


  if (CORE::length($pep) == 0) {
    logger_info("Rejecting proj of ".$tran->dbID()." because was just a stop codon");
    return 0;
  }

  my $num_stops = $pep =~ s/\*/\*/g;

  if($num_stops > $self->max_internal_stops()) {
    $self->warning("The internal stop count (".$num_stops.") is greater than the allowed value: ".$self->max_internal_stops());
    return 0;
  }

  ##################
  # number of stops is non-zero but acceptable. Need to 
  # operate on the transcript to jump over the stops
  if($num_stops) {
    $tran = replace_stops_with_introns($tran,$self->max_internal_stops());
  }

  return $tran;
}

###########################
# gets/sets
###########################

sub genomic_align_block_chains {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_gen_al_chains',$val);
  }

  return $self->param('_gen_al_chains');
}

sub target_slices {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_target_slices',$val);
  }

  unless($self->param_is_defined('_target_slices')) {
    $self->param('_target_slices',{});
  }

  return $self->param('_target_slices');
}

sub max_internal_stops {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_max_internal_stops',$val);
  }

  return $self->param('_max_internal_stops');
}

####################################
# config variable holders
####################################
#
# transcript editing and filtering
#


sub TRANSCRIPT_FILTER {
   my ($self, $val) = @_;

  if (defined $val) {
    $self->param('TRANSCRIPT_FILTER',$val);
  }

  if ($self->param_is_defined('TRANSCRIPT_FILTER')) {
    return $self->param('TRANSCRIPT_FILTER');
  }
  else {
    return;
  }
}


sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->param('_runnable_filter',$val);
  }

 # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
  if ($self->TRANSCRIPT_FILTER and !$self->param_is_defined('_runnable_filter')) {
    if (not ref($self->TRANSCRIPT_FILTER) eq "HASH" or
        not exists($self->TRANSCRIPT_FILTER->{OBJECT}) or
        not exists($self->TRANSCRIPT_FILTER->{PARAMETERS})) {

      $self->throw("FILTER in config for '".$self->analysis->logic_name."' must be a hash ref with elements:\n" .
            "  OBJECT : qualified name of the filter module;\n" .
            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
    } else {
      $self->require_module($self->TRANSCRIPT_FILTER->{OBJECT});
      $self->filter($self->TRANSCRIPT_FILTER->{OBJECT}->new(%{$self->TRANSCRIPT_FILTER->{PARAMETERS}}));
    }
  }
  if ($self->param_is_defined('_runnable_filter')) {
    return $self->param('_runnable_filter');
  }
  else {
    return;
  }
}


sub MIN_COVERAGE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('MIN_COVERAGE',$val);
  }

  return $self->param('MIN_COVERAGE');
}

1;
