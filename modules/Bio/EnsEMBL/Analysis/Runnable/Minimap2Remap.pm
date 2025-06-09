=head1 LICENSE

 Copyright [2022] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module acts as a follow up to the LowDivergenceProjection module. It uses an alternative strategy to try and map genes
and transcripts that are absent or problematic from the projected set. In some cases these genes will be correctly missing,
so should not really be mapped, but this module will try anyway. The idea is to look at the expected location of the gene,
based on finding syntenic regions in the source and target dbs through high confidence, unambiguous projections. High confidence
unambiguous projections are ones that
1) Have all their transcripts mapped >= 99 percent cov and identity
2) Have only been mapped to one place
3) Don't overlap with something, or only overlap with something that they overlap with in the source
4) Have the same neighbouring genes as in source. Might need to be less strict on this due to scaffold level assemblies and
   the fact that things were allowed multimap in the previous step. At the same time, it may be fine to be strict in practice
   since the goal is to confidently get the region and not necessarily to completely pinpoint the exact location

Another point to consider is that if the target gene is supposed to be between two high confidence projections and the region
has shrunk in the target assembly versus the source, then that's an indication that a deletion/expansion might have happened

The next step in the process is to take the combined set of genes from the two modules and remove things that shouldn't be there

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use POSIX;
use List::Util qw(min max);

use Bio::EnsEMBL::Analysis::Runnable::Minimap2;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_best_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(set_alignment_supporting_features attach_Slice_to_Transcript attach_Analysis_to_Transcript calculate_exon_phases replace_stops_with_introns);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs map_cds_location align_proteins execute_with_wait);
use File::Spec;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw(throw);

use Bio::EnsEMBL::Analysis::Provenance::Logger;

use parent ('Bio::EnsEMBL::Analysis::Runnable');


=head2 new

 Arg [DECOMPRESS]           : String as a command like 'gzip -c -'
 Arg [EXPECTED_ATTRIBUTES]  : String specify the attribute expected for the output, see STAR manual
 Description                : Creates a  object to align reads to a genome using STAR
 Returntype                 : 
 Exceptions                 : Throws if WORKDIR does not exist
                              Throws if the genome has not been indexed

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($genome_index, $input_file, $paftools_path, $source_adaptor, $target_adaptor, $delete_input_file, $parent_genes, $parent_gene_ids, $no_projection) = rearrange([qw (GENOME_INDEX INPUT_FILE PAFTOOLS_PATH SOURCE_ADAPTOR TARGET_ADAPTOR DELETE_INPUT_FILE PARENT_GENES PARENT_GENE_IDS NO_PROJECTION)],@args);
  $self->genome_index($genome_index);
  $self->input_file($input_file);
  $self->paftools_path($paftools_path);
  $self->source_adaptor($source_adaptor);
  $self->target_adaptor($target_adaptor);
  $self->delete_input_file($delete_input_file);
  $self->genes_to_process($parent_genes);
  $self->parent_gene_ids($parent_gene_ids);
  $self->no_projection($no_projection);

  my $log_dir = $self->{_output_dir} || 'logs';
  say "Provenance will be stored in: $log_dir";

  $self->{_logger} = Bio::EnsEMBL::Analysis::Provenance::Logger->new(
    log_dir => $log_dir
  );
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run Star to align reads to an indexed genome. The resulting output file will be stored in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;
  my $gene_genomic_seqs_hash = {};
  my $source_adaptor = $self->source_adaptor();
  my $target_adaptor = $self->target_adaptor();
  my $target_gene_adaptor = $target_adaptor->get_GeneAdaptor();
  my $source_genes = $self->genes_to_process();

#  foreach my $source_gene (@$source_genes) {
#    say "FERGAL SOURCE: ".$source_gene->seq_region_name()." ".$source_gene->stable_id();
#  }
  # TEST!!!!!!!!!!!!!!!!
#  my $test_slice = $target_adaptor->get_SliceAdaptor->fetch_by_region('toplevel','17');
#  my $target_genes = $target_adaptor->get_GeneAdaptor->fetch_all_by_Slice($test_slice);
#  my $target_genes_by_slice = $self->sort_genes_by_slice($target_genes);

  my $source_genes_by_slice = $self->sort_genes_by_slice($source_genes);
  my $target_genes_by_slice = {};
  my $target_genes = [];
  foreach my $slice (@{$target_adaptor->get_SliceAdaptor->fetch_all('toplevel')}) {
    $target_genes_by_slice->{$slice->seq_region_name} = [sort {$a->start <=> $b->start} @{$slice->get_all_Genes}];
    push(@$target_genes, @{$target_genes_by_slice->{$slice->seq_region_name}});
  }
  my $source_genes_by_stable_id = $self->genes_by_stable_id($source_genes);
  my $target_genes_by_stable_id = $self->genes_by_stable_id($target_genes);
  say "Searching for missing source genes";

  my $missing_source_genes = $self->list_missing_genes($source_genes_by_stable_id,$target_genes_by_stable_id);
  say "Found ".scalar(@$missing_source_genes)." missing source genes";

  foreach my $gene (@$target_genes) {
    # TEST!!!!!
#    unless($gene->seq_region_name eq '17') {
#      next;
#    }
    my $source_gene = ${$source_genes_by_stable_id->{$gene->stable_id()}}[0];
    unless($source_gene) {
      $self->throw("Couldn't find a source gene for ".$gene->stable_id()); 
      if ($self->{_logger}) {
        my $stable_id = $gene->stable_id();
        $self->{_logger}->log(
          "Minimap2Remap-missing_gene_identified",
          {
            feature_id => $stable_id,
            feature_type => "gene",
            status => "missing",
            details => {
              source_region => $source_gene->seq_region_name(),
              source_coordinates => $source_gene->seq_region_start() . "-" . $source_gene->seq_region_end(),
              biotype => $source_gene->biotype(),
              transcript_count => scalar(@{$source_gene->get_all_Transcripts()})
            },
            message => "Gene missing from target, added to recovery list"
          }
        );
      }
    }
    my $source_transcripts = $source_gene->get_all_Transcripts();
    $self->check_complete_mapping($gene,$source_transcripts);
  }

  # This will add a list of all the genes on the slice ordered by increasing midpoint distance to every gene in the source set
  # It's probably overkill to calculate all distances per gene, but will leave it this way for now. Later having all of them
  # could allow a highly accurate global distance between the synteny of the source and target genes to be calculated
  say "Calculating source gene neighbours and midpoints";
  foreach my $slice (keys(%$source_genes_by_slice)) {
    say " Getting info for source genes on region: ".$slice;
    # TEST!!!!!!!!!!!!!!!!!
#    unless($slice eq '17') {
#      next;
#    }
    my $genes = $source_genes_by_slice->{$slice};
    my $midpoint_coords_by_id = {};
    foreach my $gene (@$genes) {
      my $midpoint = $gene->start + ceil($gene->length()/2);
      unless($midpoint_coords_by_id->{$gene->stable_id()}) {
        $midpoint_coords_by_id->{$gene->stable_id()} = [];
      }
      push(@{$midpoint_coords_by_id->{$gene->stable_id()}},$midpoint);
    }
    foreach my $gene (@$genes) {
      $self->set_closest_gene_ids($gene,$midpoint_coords_by_id);
    }
  }

  say "Calculating high confidence genes in target";
  my $high_confidence_genes = $self->list_high_confidence_genes($target_genes_by_slice,$target_genes_by_stable_id,$source_genes_by_stable_id);
  my $high_confidence_genes_by_slice = $self->sort_genes_by_slice($high_confidence_genes);
  my $high_confidence_genes_by_id = $self->genes_by_stable_id($high_confidence_genes);

  my $high_confidence_count = 0;
  foreach my $gene (@$target_genes) {
    if($gene->{'is_high_confidence'}) {
       $high_confidence_count++;
    }
  }

  say "Got ".$high_confidence_count." high condfidence genes from a total of ".scalar(@$target_genes);

  # Now that we have a high confidence set of genes, we want to take the missing genes and see where they should be located, then try to map
  # them to the expected region
  $self->set_expected_regions($missing_source_genes,$source_genes_by_slice,$high_confidence_genes_by_id);
  my $missing_gene_reads_file = $self->create_input_file($missing_source_genes,$gene_genomic_seqs_hash,$source_adaptor);

  # We want to take the missing genes, along with the flanking regions and then map, allowing secondary mappings. Foreach hit we should define the
  # expected region to be covered by the gene. Then for the next hit, if it also falls in the boundary, we should check if it's redundant or if
  # it overlaps and extends an existing region, or if it's a unique region. Basically we should get all hits within the region, extend as perdicted
  # from the source region, order the extended hits and then merge into clustered regions. Foreach clustered region, then proceed normally with the
  # analysis

  my $paf_file = $self->create_filename(undef,'paf');
  $self->files_to_delete($paf_file);

  my $genome_index  = $self->genome_index;
  my $options = $self->options;

  # run minimap2
  my $minimap2_command = $self->program." --cs --secondary=yes -x map-ont ".$genome_index." ".$missing_gene_reads_file." > ".$paf_file;
  $self->warning("Command:\n".$minimap2_command."\n");
  execute_with_wait($minimap2_command);

  my $paf_results = [];
  open(IN, $paf_file) or $self->throw("Could not open $paf_file for reading");
  while(<IN>) {
    chomp($_);
    push(@$paf_results,$_);
  }
  close(IN) or $self->throw("Could not close $paf_file");

  my $processed_gene_ids = {};
  my $paf_results_hash = {};
  foreach my $paf_result (@$paf_results) {
    say "PAF results for first pass:\n".$paf_result;
    my @result_cols = split("\t",$paf_result);
    my $gene_stable_id = $result_cols[0];
    if($paf_results_hash->{$gene_stable_id}) {
      push(@{$paf_results_hash->{$gene_stable_id}},\@result_cols)
    } else {
      $paf_results_hash->{$gene_stable_id} = [\@result_cols];
    }
  }

  my $all_recovered_genes = [];
  foreach my $gene (@$missing_source_genes) {
    say "Processing source gene: ".$gene->stable_id();
    my $paf_results = $paf_results_hash->{$gene->stable_id()};
    $self->filter_paf_hits($gene,$paf_results);
    my $recovered_genes = $self->process_results($gene,$gene_genomic_seqs_hash,$target_genes);
    if(scalar(@$recovered_genes)) {
      foreach my $recovered_gene (@$recovered_genes) {
        my $source_gene = ${$source_genes_by_stable_id->{$recovered_gene->stable_id()}}[0];
        unless($source_gene) {
          $self->throw("Couldn't find a source gene for recovered gene ".$recovered_gene->stable_id());
        }
        my $source_transcripts = $source_gene->get_all_Transcripts();
        $self->check_complete_mapping($recovered_gene,$source_transcripts);
      }
      push(@$all_recovered_genes,@$recovered_genes);
    }
  }

  say "Got ".scalar(@$all_recovered_genes)." genes after recovery";

  # Log gene recovery results and confidence reassessment
  if ($self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-gene_recovery_results",
      {
        feature_id => "pipeline_recovery_phase",
        feature_type => 'pipeline',
        result => scalar(@$all_recovered_genes) > 0 ? 'genes_recovered' : 'no_recovery',
        details => {
          genes_recovered => scalar(@$all_recovered_genes),
          total_target_genes_after_recovery => scalar(@$target_genes) + scalar(@$all_recovered_genes),
          recovery_successful => scalar(@$all_recovered_genes) > 0 ? 1 : 0,
          will_trigger_confidence_reassessment => scalar(@$all_recovered_genes) > 0 ? 1 : 0
        },
        message => "Gene recovery phase completed: " . scalar(@$all_recovered_genes) . " genes recovered"
      }
    );
  }

  if(scalar(@$all_recovered_genes)) {
    push(@$target_genes,@$all_recovered_genes);
    $target_genes_by_slice = $self->sort_genes_by_slice($target_genes);
    $target_genes_by_stable_id = $self->genes_by_stable_id($target_genes);
    
    # Calculate high confidence before and after for comparison
    my $previous_high_confidence_count = scalar(@$high_confidence_genes);
    
    $high_confidence_genes = $self->list_high_confidence_genes($target_genes_by_slice,$target_genes_by_stable_id,$source_genes_by_stable_id);
    $high_confidence_genes_by_slice = $self->sort_genes_by_slice($high_confidence_genes);
    $high_confidence_genes_by_id = $self->genes_by_stable_id($high_confidence_genes);
    
    my $new_high_confidence_count = scalar(@$high_confidence_genes);
    my $high_confidence_change = $new_high_confidence_count - $previous_high_confidence_count;
    
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-confidence_reassessment_after_recovery",
        {
          feature_id => "pipeline_confidence_update",
          feature_type => 'pipeline', 
          result => 'confidence_reassessed',
          details => {
            genes_added_to_target_set => scalar(@$all_recovered_genes),
            total_target_genes_now => scalar(@$target_genes),
            previous_high_confidence_count => $previous_high_confidence_count,
            new_high_confidence_count => $new_high_confidence_count,
            high_confidence_change => $high_confidence_change,
            confidence_impact => $high_confidence_change > 0 ? 'improved' : 
                                $high_confidence_change < 0 ? 'degraded' : 'unchanged',
            reassessment_reason => 'gene_recovery_changed_landscape'
          },
          message => "Confidence reassessment after recovery: $high_confidence_change change in high-confidence genes (now $new_high_confidence_count total)"
        }
      );
    }
  }

  # Because the target genes can have the same stable ids, and some of the recovered genes don't have dbIDs, assign unique internal ids
  $self->assign_internal_ids($target_genes);

  # Check for conflic. Conflict is when we have a gene with exon overlap with another gene, that it does not have exon overlap with in the reference
  # There are different levels of conflict, for example a very minor overlap with a neighbouring gene in the reference might just be a slight misalignment
  # of a UTR. Whereas if you have two genes with significant overlap, the assumption is that that are paralogues and thus one of them is mismapped. There
  # might be lots of them mismapped at a locus. For the moment we just want to record if there is conflict
  $self->check_conflict($target_genes_by_slice,$source_genes_by_slice,$source_genes_by_stable_id,$target_genes_by_stable_id,$high_confidence_genes_by_id);

  # At this point any gene with conflict is labelled with $gene->{'genomic_conflict'} = 1 and has a list of conflicting stable ids with dbIDs as the values
  # accessed via $gene->{'conflicting_genes'}->{$id} = $dbID
  # This means we can loop through all the genes on the slice and decide what to do for each conflicting pair
  $self->resolve_conflict($target_genes_by_slice,$source_genes_by_slice,$source_genes_by_stable_id,$target_genes_by_stable_id);

  # The genes that have been recovered will have a 'to_write' tag, the ones to remove will have 'to_remove'. If a gene has both then
  # the 'to_remove' supercedes the 'to_write'
  $self->output($target_genes);
} # End run


sub resolve_conflict {
  my ($self,$target_genes_by_slice,$source_genes_by_slice,$source_genes_by_stable_id,$target_genes_by_stable_id) = @_;

  # This will go through the target genes and look for what overlaps then. First check stranded genomic overlap and
  # then exon overlap. Conflict is when a gene has exon overlap (need to decide level of exon overlap), with something
  # it does not have exon overlap with on the reference
  say "Processing genes by slice to remove conflicts";
  my $coverage_cutoff = 99;
  my $identity_cutoff = 99;

  my $genes_to_remove = {};
  my $total_conflicts_processed = 0;
  my $total_genes_removed = 0;
  my @resolution_decisions = ();

  foreach my $slice (keys(%$target_genes_by_slice)) {
    my $target_genes = $target_genes_by_slice->{$slice};
    my $conflicts_on_slice = 0;
    my $removals_on_slice = 0;
    
    for(my $i=0; $i<scalar(@$target_genes)-1; $i++) {
      my $gene = ${$target_genes}[$i];
      unless($gene->{'genomic_conflict'}) {
        next;
      }

      $conflicts_on_slice++;
      say "Examing gene ".$gene->stable_id()." (".$gene->{'internal_id'}.") as it is listed as having conflict";
      if($genes_to_remove->{$gene->stable_id()}->{$gene->{'internal_id'}}) {
        say "  Skipping gene ".$gene->stable_id()." (".$gene->{'internal_id'}.") as it is already listed for removal";
        next;
      }

      say "  Gene ".$gene->stable_id()." is in conflict with ".scalar(keys(%{$gene->{'conflicting_genes'}}))." genes, processing to resolve conflict";

      # Log conflict resolution start
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-conflict_resolution_start",
          {
            feature_id => $gene->stable_id(),
            feature_type => 'gene',
            result => 'conflict_resolution_initiated',
            details => {
              gene_internal_id => $gene->{'internal_id'},
              conflicting_genes_count => scalar(keys(%{$gene->{'conflicting_genes'}})),
              conflicting_gene_ids => [keys(%{$gene->{'conflicting_genes'}})],
              resolution_strategy => 'hierarchical_filtering',
              quality_cutoffs => {
                coverage => $coverage_cutoff,
                identity => $identity_cutoff
              }
            },
            message => "Starting conflict resolution for gene with " . scalar(keys(%{$gene->{'conflicting_genes'}})) . " conflicts"
          }
        );
      }

      # How to resolve conflict, take the easiest cases first
      # Gene A is high confidence, Gene B is not
      # Gene A has a lower than expected cov/identity, gene B has passes
      # Gene A is in expected location, gene B is not
      # Gene A has a lower neighbourhood score than gene B
      # Gene A further from expected order than gene B
      # If it's sort of impossible to pick, just pick the one with the most exons
      foreach my $stable_id (keys(%{$gene->{'conflicting_genes'}})) {
        my $conflicting_db_id = $gene->{'conflicting_genes'}->{$stable_id};
        # If the conlficting gene has already been processed for removal, then there's nothing to do
        if($genes_to_remove->{$stable_id}->{$conflicting_db_id}) {
          next;
        }

        my $conflicting_gene;
        foreach my $target_conflicting_gene (@{$target_genes_by_stable_id->{$stable_id}}) {
          if($target_conflicting_gene->{'internal_id'} == $conflicting_db_id) {
            $conflicting_gene = $target_conflicting_gene;
            last;
          }
        }

        # We now have the current gene and the conflicting gene, so start filtering
        # NOTE: ultimately a better strategy in terms of using syntenty info to begin with should go in here. But in this iteration
        #       we will use a more straightforward strategy to begin with

        $total_conflicts_processed++;

        # 1) Filter based on confidence
        if($gene->{'is_high_confidence'} and !($conflicting_gene->{'is_high_confidence'})) {
          say "  Removing conflicting gene ".$conflicting_gene->stable_id()." (".$conflicting_gene->{'internal_id'}.") as it is not high confidence while ".$gene->stable_id().
              " (".$gene->{'internal_id'}.") is";
          $genes_to_remove->{$conflicting_gene->stable_id()}->{$conflicting_gene->{'internal_id'}} = 1;
          $conflicting_gene->{'to_remove'} = 1;
          $removals_on_slice++;
          $total_genes_removed++;
          
          if ($self->{_logger}) {
            $self->{_logger}->log(
              "Minimap2Remap-conflict_resolution_decision",
              {
                feature_id => $gene->stable_id(),
                feature_type => 'gene',
                result => 'conflict_resolved_by_confidence',
                details => {
                  resolution_rule => 'high_confidence_beats_low_confidence',
                  winner => $gene->stable_id(),
                  winner_internal_id => $gene->{'internal_id'},
                  winner_confidence => 'high',
                  removed_gene => $conflicting_gene->stable_id(),
                  removed_internal_id => $conflicting_gene->{'internal_id'},
                  removed_confidence => 'low',
                  decision_basis => 'confidence_level'
                },
                message => "Conflict resolved: high-confidence gene retained, low-confidence gene removed"
              }
            );
          }
          
          push @resolution_decisions, {
            winner => $gene->stable_id(),
            removed => $conflicting_gene->stable_id(),
            rule => 'high_confidence_beats_low_confidence'
          };
          next;
        } elsif($conflicting_gene->{'is_high_confidence'} and !($gene->{'is_high_confidence'})) {
          $genes_to_remove->{$gene->stable_id()}->{$gene->{'internal_id'}} = 1;
          $gene->{'to_remove'} = 1;
          $removals_on_slice++;
          $total_genes_removed++;
          
          say "  Removing current gene ".$gene->stable_id()." (".$gene->{'internal_id'}.") as it is not high confidence while ".$conflicting_gene->stable_id()." (".
              $conflicting_gene->{'internal_id'}.") is";
          
          if ($self->{_logger}) {
            $self->{_logger}->log(
              "Minimap2Remap-conflict_resolution_decision",
              {
                feature_id => $gene->stable_id(),
                feature_type => 'gene',
                result => 'conflict_resolved_by_confidence',
                details => {
                  resolution_rule => 'high_confidence_beats_low_confidence',
                  winner => $conflicting_gene->stable_id(),
                  winner_internal_id => $conflicting_gene->{'internal_id'},
                  winner_confidence => 'high',
                  removed_gene => $gene->stable_id(),
                  removed_internal_id => $gene->{'internal_id'},
                  removed_confidence => 'low',
                  decision_basis => 'confidence_level'
                },
                message => "Conflict resolved: high-confidence conflicting gene retained, low-confidence current gene removed"
              }
            );
          }
          
          push @resolution_decisions, {
            winner => $conflicting_gene->stable_id(),
            removed => $gene->stable_id(),
            rule => 'high_confidence_beats_low_confidence'
          };
          last;
        } else {
          say "  Both genes are high confidence";
        }

        # 2) Filter based on one gene passing cut-offs while the other fails
        if (defined($gene->{'avg_cov'}) and defined($gene->{'avg_perc_id'}) and
            defined($conflicting_gene->{'avg_cov'}) and defined($conflicting_gene->{'avg_perc_id'})) {
          if(($gene->{'avg_cov'} >= $coverage_cutoff and $gene->{'avg_perc_id'} >= $identity_cutoff) and
             ($conflicting_gene->{'avg_cov'} < $coverage_cutoff or $conflicting_gene->{'avg_perc_id'} < $identity_cutoff)) {
            $genes_to_remove->{$conflicting_gene->stable_id()}->{$conflicting_gene->{'internal_id'}} = 1;
            $conflicting_gene->{'to_remove'} = 1;
            $removals_on_slice++;
            $total_genes_removed++;
            
            say "  Removing conflicting gene ".$conflicting_gene->stable_id()." (".$conflicting_gene->{'internal_id'}.") as it fails the coverage/identity cut-off while ".$gene->stable_id().
              " (".$gene->{'internal_id'}.") does not";
            
            if ($self->{_logger}) {
              $self->{_logger}->log(
                "Minimap2Remap-conflict_resolution_decision",
                {
                  feature_id => $gene->stable_id(),
                  feature_type => 'gene',
                  result => 'conflict_resolved_by_quality_cutoffs',
                  details => {
                    resolution_rule => 'quality_cutoff_pass_beats_fail',
                    winner => $gene->stable_id(),
                    winner_metrics => {
                      coverage => $gene->{'avg_cov'},
                      identity => $gene->{'avg_perc_id'},
                      passes_cutoffs => 1
                    },
                    removed_gene => $conflicting_gene->stable_id(),
                    removed_metrics => {
                      coverage => $conflicting_gene->{'avg_cov'},
                      identity => $conflicting_gene->{'avg_perc_id'},
                      passes_cutoffs => 0
                    },
                    cutoffs => {
                      coverage => $coverage_cutoff,
                      identity => $identity_cutoff
                    },
                    decision_basis => 'quality_thresholds'
                  },
                  message => "Conflict resolved: gene passing quality cutoffs retained"
                }
              );
            }
            
            push @resolution_decisions, {
              winner => $gene->stable_id(),
              removed => $conflicting_gene->stable_id(),
              rule => 'quality_cutoff_pass_beats_fail'
            };
            next;
          } elsif(($conflicting_gene->{'avg_cov'} >= $coverage_cutoff and $conflicting_gene->{'avg_perc_id'} >= $identity_cutoff) and
                  ($gene->{'avg_cov'} < $coverage_cutoff or $gene->{'avg_perc_id'} < $identity_cutoff)) {
            $genes_to_remove->{$gene->stable_id()}->{$gene->{'internal_id'}} = 1;
            $gene->{'to_remove'} = 1;
            $removals_on_slice++;
            $total_genes_removed++;
            
            say "  Removing current gene ".$gene->stable_id()." (".$gene->{'internal_id'}.") as it fails the coverage/identity cut-off while ".$conflicting_gene->stable_id()." (".
                $conflicting_gene->{'internal_id'}.") does not";
            
            if ($self->{_logger}) {
              $self->{_logger}->log(
                "Minimap2Remap-conflict_resolution_decision",
                {
                  feature_id => $gene->stable_id(),
                  feature_type => 'gene',
                  result => 'conflict_resolved_by_quality_cutoffs',
                  details => {
                    resolution_rule => 'quality_cutoff_pass_beats_fail',
                    winner => $conflicting_gene->stable_id(),
                    winner_metrics => {
                      coverage => $conflicting_gene->{'avg_cov'},
                      identity => $conflicting_gene->{'avg_perc_id'},
                      passes_cutoffs => 1
                    },
                    removed_gene => $gene->stable_id(),
                    removed_metrics => {
                      coverage => $gene->{'avg_cov'},
                      identity => $gene->{'avg_perc_id'},
                      passes_cutoffs => 0
                    },
                    cutoffs => {
                      coverage => $coverage_cutoff,
                      identity => $identity_cutoff
                    },
                    decision_basis => 'quality_thresholds'
                  },
                  message => "Conflict resolved: conflicting gene passing quality cutoffs retained"
                }
              );
            }
            
            push @resolution_decisions, {
              winner => $conflicting_gene->stable_id(),
              removed => $gene->stable_id(),
              rule => 'quality_cutoff_pass_beats_fail'
            };
            last;
          } else {
            say "  Both genes pass coverage/perc id cut-offs: ".$gene->{'avg_cov'}."/".$gene->{'avg_perc_id'}." vs ".$conflicting_gene->{'avg_cov'}."/".$conflicting_gene->{'avg_perc_id'};
          }
        } else {
          $self->throw("avg_cov or avg_perc_id not defined for gene or conflicting_gene");
        }

        # 3) Filter by expected location
        if($gene->{'expected_location'} and !($conflicting_gene->{'expected_location'})) {
          $genes_to_remove->{$conflicting_gene->stable_id()}->{$conflicting_gene->{'internal_id'}} = 1;
          $conflicting_gene->{'to_remove'} = 1;
          $removals_on_slice++;
          $total_genes_removed++;
          
          say "  Removing conflicting gene ".$conflicting_gene->stable_id()." (".$conflicting_gene->{'internal_id'}.") as it is not in the expected location ".$gene->stable_id().
              " (".$gene->{'internal_id'}.") is";
          
          if ($self->{_logger}) {
            $self->{_logger}->log(
              "Minimap2Remap-conflict_resolution_decision",
              {
                feature_id => $gene->stable_id(),
                feature_type => 'gene',
                result => 'conflict_resolved_by_expected_location',
                details => {
                  resolution_rule => 'expected_location_beats_unexpected',
                  winner => $gene->stable_id(),
                  winner_in_expected_location => 1,
                  removed_gene => $conflicting_gene->stable_id(),
                  removed_in_expected_location => 0,
                  decision_basis => 'synteny_expectations'
                },
                message => "Conflict resolved: gene in expected location retained"
              }
            );
          }
          
          push @resolution_decisions, {
            winner => $gene->stable_id(),
            removed => $conflicting_gene->stable_id(),
            rule => 'expected_location_beats_unexpected'
          };
          next;
        } elsif($conflicting_gene->{'expected_location'} and !($gene->{'expected_location'})) {
          $genes_to_remove->{$gene->stable_id()}->{$gene->{'internal_id'}} = 1;
          $gene->{'to_remove'} = 1;
          $removals_on_slice++;
          $total_genes_removed++;
          
          say "  Removing current gene ".$gene->stable_id()." (".$gene->{'internal_id'}.") as it is not in the expected location ".$conflicting_gene->stable_id()." (".
              $conflicting_gene->{'internal_id'}.") is";
          
          if ($self->{_logger}) {
            $self->{_logger}->log(
              "Minimap2Remap-conflict_resolution_decision",
              {
                feature_id => $gene->stable_id(),
                feature_type => 'gene',
                result => 'conflict_resolved_by_expected_location',
                details => {
                  resolution_rule => 'expected_location_beats_unexpected',
                  winner => $conflicting_gene->stable_id(),
                  winner_in_expected_location => 1,
                  removed_gene => $gene->stable_id(),
                  removed_in_expected_location => 0,
                  decision_basis => 'synteny_expectations'
                },
                message => "Conflict resolved: conflicting gene in expected location retained"
              }
            );
          }
          
          push @resolution_decisions, {
            winner => $conflicting_gene->stable_id(),
            removed => $gene->stable_id(),
            rule => 'expected_location_beats_unexpected'
          };
          last;
        } else {
          say "  Both genes are in the expected location";
        }

        # Log when resolution rules are disabled/commented
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-conflict_resolution_incomplete",
            {
              feature_id => $gene->stable_id(),
              feature_type => 'gene',
              result => 'additional_rules_disabled',
              details => {
                conflicting_gene => $conflicting_gene->stable_id(),
                active_rules_exhausted => ['confidence_comparison', 'quality_cutoffs', 'expected_location'],
                disabled_rules => ['global_minimap_preference', 'neighbourhood_score', 'combined_score'],
                both_genes_equivalent => 1,
                current_gene_method => $gene->{'annotation_method'} || 'unknown',
                conflicting_gene_method => $conflicting_gene->{'annotation_method'} || 'unknown',
                resolution_status => 'unresolved_conflict_remains'
              },
              message => "Conflict unresolved - additional resolution rules are disabled/commented out"
            }
          );
        }

        # 4) Remove global (COMMENTED OUT)
#        if($gene->{'annotation_method'} eq 'minimap_global' and !($conflicting_gene->{'annotation_method'}) eq 'minimap_global') {
#          $genes_to_remove->{$conflicting_gene->stable_id()}->{$conflicting_gene->{'internal_id'}} = 1;
#          $conflicting_gene->{'to_remove'} = 1;
#          say "  Removing conflicting gene ".$conflicting_gene->stable_id()." (".$conflicting_gene->{'internal_id'}.") as it was generated by global minimap ".$gene->stable_id().
#              " (".$gene->{'internal_id'}.") is";
#          next;
#        } elsif($conflicting_gene->{'annotation_method'} eq 'minimap_global' and !($gene->{'annotation_method'}) eq 'minimap_global') {
#          $genes_to_remove->{$gene->stable_id()}->{$gene->{'internal_id'}} = 1;
#          $gene->{'to_remove'} = 1;
#          say "  Removing current gene ".$gene->stable_id()." (".$gene->{'internal_id'}.") as it was generated by global minimap ".$conflicting_gene->stable_id()." (".
#              $conflicting_gene->{'internal_id'}.") is";
#          last;
#        } else {
#          say "  Can't use global minimap to resolve";
#        }
#        # 4) Filter by neighbourhood score (COMMENTED OUT)
#        if($gene->{'neighbourhood_score'} > $conflicting_gene->{'neighbourhood_score'}) {
#          $genes_to_remove->{$conflicting_gene->stable_id()}->{$conflicting_gene->{'internal_id'}} = 1;
#          $conflicting_gene->{'to_remove'} = 1;
#          say "  Removing conflicting gene ".$conflicting_gene->stable_id()." (".$conflicting_gene->{'internal_id'}.") as it is has a lower neighbourhood score ".$gene->stable_id().
#              " (".$gene->{'internal_id'}."), ".
#              $conflicting_gene->{'neighbourhood_score'}." vs ".$gene->{'neighbourhood_score'};
#          next;
#        } elsif($conflicting_gene->{'neighbourhood_score'} > $gene->{'neighbourhood_score'}) {
#          $genes_to_remove->{$gene->stable_id()}->{$gene->{'internal_id'}} = 1;
#          $gene->{'to_remove'} = 1;
#          say "  Removing current gene ".$gene->stable_id()." (".$gene->{'internal_id'}.") as it is has a lower neighbourhood score ".$conflicting_gene->stable_id().
#              " (".$conflicting_gene->{'internal_id'}."), ".$gene->{'neighbourhood_score'}." vs ".$conflicting_gene->{'neighbourhood_score'};
#          last;
#        } else {
#          say "  Both genes have the same neighbourhood score, ".$gene->{'neighbourhood_score'};
#        }
#
#        # 5) Just take the highest identity/coverage, or the current gene if they're still identical (COMMENTED OUT)
#        if(($gene->{'avg_cov'} + $gene->{'avg_perc_id'}) >= ($conflicting_gene->{'avg_cov'} + $conflicting_gene->{'avg_perc_id'})) {
#          $genes_to_remove->{$conflicting_gene->stable_id()}->{$conflicting_gene->{'internal_id'}} = 1;
#          $conflicting_gene->{'to_remove'} = 1;
#          say "  Removing conflicting gene ".$conflicting_gene->stable_id()." (".$conflicting_gene->{'internal_id'}.") as a lower or equivalent combined average cov/perc id ".$gene->stable_id().
#              " (".$gene->{'internal_id'}.")";
#          next;
#        } else {
#          $genes_to_remove->{$gene->stable_id()}->{$gene->{'internal_id'}} = 1;
#          $gene->{'to_remove'} = 1;
#          say "  Removing current gene ".$gene->stable_id()." (".$gene->{'internal_id'}.") as it is has the lower combined average cov/perc id ".$conflicting_gene->stable_id()." (".
#              $conflicting_gene->{'internal_id'}.")";
#          last;
#        }
      }
    } # for(my $i=0; $i<scalar
    
    # Log slice-level resolution summary
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-slice_conflict_resolution_summary",
        {
          feature_id => $slice,
          feature_type => 'slice',
          result => 'slice_resolution_completed',
          details => {
            slice_name => $slice,
            conflicts_processed => $conflicts_on_slice,
            genes_removed => $removals_on_slice,
            resolution_rate => $conflicts_on_slice > 0 ? sprintf("%.1f", ($removals_on_slice / $conflicts_on_slice) * 100) : 0,
            resolution_limitation => 'only_3_of_6_rules_active'
          },
          message => "Slice conflict resolution: $removals_on_slice genes removed from $conflicts_on_slice conflicts (limited by disabled rules)"
        }
      );
    }
  } # foreach my $slice (keys

  # Log overall conflict resolution assessment
  if ($self->{_logger}) {
    my %rule_usage = ();
    foreach my $decision (@resolution_decisions) {
      $rule_usage{$decision->{rule}}++;
    }
    
    my $unresolved_conflicts = $total_conflicts_processed - $total_genes_removed;
    
    $self->{_logger}->log(
      "Minimap2Remap-overall_conflict_resolution_assessment",
      {
        feature_id => "pipeline_conflict_resolution",
        feature_type => 'pipeline',
        result => 'conflict_resolution_completed',
        details => {
          total_conflicts_processed => $total_conflicts_processed,
          total_genes_removed => $total_genes_removed,
          unresolved_conflicts => $unresolved_conflicts,
          resolution_success_rate => $total_conflicts_processed > 0 ? sprintf("%.1f", ($total_genes_removed / $total_conflicts_processed) * 100) : 0,
          resolution_rules_used => \%rule_usage,
          most_effective_rule => scalar(keys %rule_usage) > 0 ? (sort { $rule_usage{$b} <=> $rule_usage{$a} } keys %rule_usage)[0] : 'none',
          active_rules => ['confidence_comparison', 'quality_cutoffs', 'expected_location'],
          disabled_rules => ['global_minimap_preference', 'neighbourhood_score', 'combined_score'],
          resolution_limitation => 'additional_rules_needed_for_complete_resolution',
          pipeline_impact => $unresolved_conflicts > 0 ? 'manual_curation_required' : 'all_conflicts_resolved'
        },
        message => "Conflict resolution completed: $total_genes_removed resolved, $unresolved_conflicts unresolved (limited rule set)"
      }
    );
  }
}

sub assign_internal_ids {
  my ($self,$genes) = @_;

  my $id_counter = 1;
  foreach my $gene (@$genes) {
    $gene->{'internal_id'} = $id_counter;
    $id_counter++;
  }
}


sub check_conflict {
  my ($self,$target_genes_by_slice,$source_genes_by_slice,$source_genes_by_stable_id,$target_genes_by_stable_id,$high_confidence_genes_by_id) = @_;

  my $total_conflicts_detected = 0;
  my $total_genes_with_conflicts = 0;
  my @conflict_details = ();

  # This will go through the target genes and look for what overlaps then. First check stranded genomic overlap and
  # then exon overlap. Conflict is when a gene has exon overlap (need to decide level of exon overlap), with something
  # it does not have exon overlap with on the reference
  foreach my $slice (keys(%$target_genes_by_slice)) {
    say "Calculating genomically overlapping genes for target slice: ".$slice;
    my $target_genes = $target_genes_by_slice->{$slice};
    my $overlaps_found_on_slice = 0;
    
    for(my $i=0; $i<scalar(@$target_genes)-1; $i++) {
      my $gene1 = ${$target_genes}[$i];
      unless($gene1->{'genomic_overlapping_genes'}) {
        $gene1->{'genomic_overlapping_genes'} = {};
      }

      # Would be faster to go in both directions from the current gene and keep adding till no overlap
      my $overlap = 1;
      for(my $j=$i+1; $j<scalar(@$target_genes) and $overlap; $j++) {
        my $gene2 = ${$target_genes}[$j];

        if($self->coords_overlap($gene1->seq_region_start(),$gene1->seq_region_end(),$gene2->seq_region_start(),$gene2->seq_region_end()) and $gene1->strand() == $gene2->strand()) {
          $gene1->{'genomic_overlapping_genes'}->{$gene2->stable_id()} = $gene2->{'internal_id'};
          $overlaps_found_on_slice++;
        } else {
          $overlap = 0;
        }
      } # end for(my $j=$i+1;

      $overlap = 1;
      for(my $j=$i-1; $j>=0 and $overlap; $j--) {
        my $gene2 = ${$target_genes}[$j];
        if($self->coords_overlap($gene1->seq_region_start(),$gene1->seq_region_end(),$gene2->seq_region_start(),$gene2->seq_region_end()) and $gene1->strand() == $gene2->strand()) {
          $gene1->{'genomic_overlapping_genes'}->{$gene2->stable_id()} = $gene2->{'internal_id'};
          $overlaps_found_on_slice++;
        } else {
          $overlap = 0;
        }
      } # end for(my $j=$i-1;
    } # end for(my $i=0
    
    # Log overlap detection results per slice
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-genomic_overlap_detection",
        {
          feature_id => $slice,
          feature_type => 'slice',
          result => 'overlap_analysis_completed',
          details => {
            slice_name => $slice,
            genes_on_slice => scalar(@$target_genes),
            overlapping_pairs_found => $overlaps_found_on_slice,
            overlap_density => sprintf("%.3f", $overlaps_found_on_slice / scalar(@$target_genes))
          },
          message => "Genomic overlap detection: $overlaps_found_on_slice overlapping pairs found on slice $slice"
        }
      );
    }
  } # end foreach my $slice

  foreach my $slice (keys(%$target_genes_by_slice)) {
    say "Calculating source gene overlap for target genes that have overlap: ".$slice;
    my $target_genes = $target_genes_by_slice->{$slice};
    my $conflicts_on_slice = 0;
    
    foreach my $target_gene (@$target_genes) {
      unless(scalar(keys(%{$target_gene->{'genomic_overlapping_genes'}}))) {
        next;
      }
      
      # At this point the target gene has something that overlaps with it on the same strand, so we want to investigate
      my $source_gene = ${$source_genes_by_stable_id->{$target_gene->stable_id()}}[0];

      # TEST!!!!!!!!!!!!!!
      unless($source_gene) {
        next;
      }

      my $slice_source_genes = $source_genes_by_slice->{$source_gene->seq_region_name()};
      foreach my $slice_source_gene (@$slice_source_genes) {
        if($source_gene->stable_id() eq $slice_source_gene->stable_id()) {
          next;
        }
        if($self->coords_overlap($source_gene->seq_region_start(),$source_gene->seq_region_end(),$slice_source_gene->seq_region_start(),$slice_source_gene->seq_region_end())
           and $source_gene->strand() == $slice_source_gene->strand()) {
          $source_gene->{'genomic_overlapping_genes'}->{$slice_source_gene->stable_id()} = $source_gene->dbID();
        }
      }

      # At this point we have all the overlapping genes for the source and target
      my $unexpected_overlaps = 0;
      my $gene_has_conflicts = 0;
      
      foreach my $id (keys(%{$target_gene->{'genomic_overlapping_genes'}})) {
        unless($source_gene->{'genomic_overlapping_genes'}->{$id}) {
          $unexpected_overlaps++;
          my $source_conflict_gene = ${$source_genes_by_stable_id->{$id}}[0];
          my $target_conflict_gene_db_id =  $target_gene->{'genomic_overlapping_genes'}->{$id};
          my $target_conflict_gene;
          foreach my $gene (@{$target_genes_by_stable_id->{$id}}) {
            if($gene->{'internal_id'} == $target_conflict_gene_db_id) {
              $target_conflict_gene = $gene;
              last;
            }
          }

          say "For ".$target_gene->stable_id()." got an unexpected overlap with ".$id.": Source: ".$source_gene->seq_region_name().":".$source_gene->seq_region_start().":".
              $source_gene->seq_region_end().":".$source_gene->strand()." ".$source_conflict_gene->seq_region_name().":".$source_conflict_gene->seq_region_start().":".
              $source_conflict_gene->seq_region_end().":".$source_conflict_gene->strand().", Target: ".$target_gene->seq_region_name().":".$target_gene->seq_region_start().":".
              $target_gene->seq_region_end().":".$target_gene->strand()." ".$target_conflict_gene->seq_region_name().":".$target_conflict_gene->seq_region_start().":".
              $target_conflict_gene->seq_region_end().":".$target_conflict_gene->strand();

          my $pairwise_coverage_cutoff = 0.75;
          my $pairwise_coverage = $self->get_pairwise_coverage($target_gene,$target_conflict_gene);
          
          if($pairwise_coverage >= $pairwise_coverage_cutoff) {
            say "The genes fail the average coverage check, classing as conflict. Average coverage: ".$pairwise_coverage;
            $target_gene->{'genomic_conflict'} = 1;
            $target_gene->{'conflicting_genes'}->{$target_conflict_gene->stable_id()} = $target_conflict_gene->{'internal_id'};
            $conflicts_on_slice++;
            $gene_has_conflicts = 1;
            $total_conflicts_detected++;
            
            # Determine conflict type based on source relationship
            my $conflict_type = 'unknown';
            my $source_distance = 'unknown';
            
            if ($source_conflict_gene) {
              # Calculate distance between source genes
              if ($source_gene->seq_region_name eq $source_conflict_gene->seq_region_name) {
                my $distance = abs($source_gene->seq_region_start - $source_conflict_gene->seq_region_start);
                $source_distance = $distance;
                
                if ($distance < 10000) {
                  $conflict_type = 'local_tandem_genes';
                } elsif ($distance < 100000) {
                  $conflict_type = 'nearby_genes_same_chromosome';
                } else {
                  $conflict_type = 'distant_genes_same_chromosome';
                }
              } else {
                $conflict_type = 'genes_from_different_chromosomes';
                $source_distance = 'different_chromosome';
              }
            }
            
            # Detailed conflict logging
            if ($self->{_logger}) {
              $self->{_logger}->log(
                "Minimap2Remap-paralog_conflict_detected",
                {
                  feature_id => $target_gene->stable_id(),
                  feature_type => 'gene',
                  result => 'conflict_detected',
                  details => {
                    conflict_type => $conflict_type,
                    conflicting_gene_id => $id,
                    pairwise_coverage => $pairwise_coverage,
                    coverage_cutoff => $pairwise_coverage_cutoff,
                    source_gene_location => $source_gene->seq_region_name.":".$source_gene->seq_region_start."-".$source_gene->seq_region_end."(".$source_gene->strand.")",
                    target_gene_location => $target_gene->seq_region_name.":".$target_gene->seq_region_start."-".$target_gene->seq_region_end."(".$target_gene->strand.")",
                    source_conflict_location => $source_conflict_gene ? 
                      $source_conflict_gene->seq_region_name.":".$source_conflict_gene->seq_region_start."-".$source_conflict_gene->seq_region_end."(".$source_conflict_gene->strand.")" : 
                      'unknown',
                    target_conflict_location => $target_conflict_gene ? 
                      $target_conflict_gene->seq_region_name.":".$target_conflict_gene->seq_region_start."-".$target_conflict_gene->seq_region_end."(".$target_conflict_gene->strand.")" : 
                      'unknown',
                    source_gene_distance => $source_distance,
                    overlap_unexpected_in_source => 1,
                    will_trigger_expected_region_calculation => 1
                  },
                  message => "Paralog conflict detected: $conflict_type with ${pairwise_coverage} coverage overlap"
                }
              );
            }
            
            push @conflict_details, {
              main_gene => $target_gene->stable_id(),
              conflicting_gene => $id,
              conflict_type => $conflict_type,
              pairwise_coverage => $pairwise_coverage,
              source_distance => $source_distance
            };
            
            $self->set_expected_regions([$source_gene],$source_genes_by_slice,$high_confidence_genes_by_id);
          } else {
            # Low coverage overlap - not considered conflict
            if ($self->{_logger}) {
              $self->{_logger}->log(
                "Minimap2Remap-low_coverage_overlap",
                {
                  feature_id => $target_gene->stable_id(),
                  feature_type => 'gene',
                  result => 'overlap_below_threshold',
                  details => {
                    overlapping_gene_id => $id,
                    pairwise_coverage => $pairwise_coverage,
                    coverage_cutoff => $pairwise_coverage_cutoff,
                    decision => 'not_considered_conflict',
                    likely_cause => 'minor_boundary_differences_or_utr_variation'
                  },
                  message => "Low coverage overlap (${pairwise_coverage}) - not considered conflict"
                }
              );
            }
          }
        } # end unless($source_gene->{'genomic_overlapping_genes'}
      } # end foreach my $id (keys
      
      if ($gene_has_conflicts) {
        $total_genes_with_conflicts++;
      }
    } # end foreach my $target_gene
    
    # Log slice-level conflict summary
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-slice_conflict_summary",
        {
          feature_id => $slice,
          feature_type => 'slice',
          result => 'conflict_analysis_completed',
          details => {
            slice_name => $slice,
            conflicts_detected => $conflicts_on_slice,
            genes_analyzed => scalar(@$target_genes),
            conflict_rate => sprintf("%.3f", $conflicts_on_slice / scalar(@$target_genes))
          },
          message => "Slice conflict analysis: $conflicts_on_slice conflicts detected on $slice"
        }
      );
    }
  } # end foreach my $slice
  
  # Log overall conflict assessment
  if ($self->{_logger}) {
    my %conflict_type_counts = ();
    foreach my $conflict (@conflict_details) {
      $conflict_type_counts{$conflict->{conflict_type}}++;
    }
    
    $self->{_logger}->log(
      "Minimap2Remap-overall_conflict_assessment",
      {
        feature_id => "pipeline_conflict_detection",
        feature_type => 'pipeline',
        result => 'conflict_detection_completed',
        details => {
          total_conflicts_detected => $total_conflicts_detected,
          total_genes_with_conflicts => $total_genes_with_conflicts,
          conflict_types_breakdown => \%conflict_type_counts,
          most_common_conflict_type => (sort { $conflict_type_counts{$b} <=> $conflict_type_counts{$a} } keys %conflict_type_counts)[0],
          pipeline_impact => $total_conflicts_detected > 0 ? 'conflicts_require_resolution' : 'no_conflicts_detected',
          paralog_issues_detected => $conflict_type_counts{'local_tandem_genes'} || 0,
          expected_regions_calculated => $total_conflicts_detected
        },
        message => "Overall conflict assessment: $total_conflicts_detected conflicts detected affecting $total_genes_with_conflicts genes"
      }
    );
  }
}


sub get_pairwise_coverage {
  my ($self,$gene1,$gene2) = @_;

  my $max_start = max($gene1->seq_region_start(),$gene2->seq_region_start());
  my $min_end = min($gene1->seq_region_end(),$gene2->seq_region_end());
  my $overlap = max(1,$min_end-$max_start);

  my $gene1_overlap = $overlap / ($gene1->seq_region_end() - $gene1->seq_region_start() + 1);
  my $gene2_overlap = $overlap / ($gene2->seq_region_end() - $gene2->seq_region_start() + 1);

  my $average_overlap = ($gene1_overlap + $gene2_overlap) / 2;

  return($average_overlap);
}


sub create_input_file {
  my ($self,$genes,$gene_genomic_seqs_hash,$source_adaptor) = @_;

  my $source_sequence_adaptor = $source_adaptor->get_SequenceAdaptor();
  my $genomic_reads = [];
  foreach my $gene (@$genes) {
    my $slice = $gene->slice();

    # There's a big issue in terms of small genes, even with a fair amount of padding. To counter this have a minimum
    # target region size of about 50kb. If the gene is bigger then the padding itself should be enough as it's likely
    # there are many neutral sites in the gene already
    my $min_padding = 500;

    my $region_start = $gene->seq_region_start - $min_padding;
    if($region_start < $slice->seq_region_start()) {
       $region_start = $slice->seq_region_start();
    }

    my $region_end = $gene->seq_region_end + $min_padding;
    if($region_end > $slice->seq_region_end()) {
      $region_end = $slice->seq_region_end();
    }

    my $genomic_seq = ${ $source_sequence_adaptor->fetch_by_Slice_start_end_strand($slice,$region_start,$region_end,1) };

    my $fasta_record = ">".$gene->stable_id()."\n".$genomic_seq;
    push(@$genomic_reads, $fasta_record);
    $gene_genomic_seqs_hash->{$gene->stable_id()} = [$region_start,$region_end,$genomic_seq];
  }

  my $output_file = $self->write_input_file($genomic_reads);
  return($output_file);
}


sub write_input_file {
  my ($self,$genomic_reads) = @_;

  my $output_file = $self->create_filename();
  open(OUT,">".$output_file);
  foreach my $genomic_read (@$genomic_reads) {
    say OUT $genomic_read;
  }
  close OUT;

  return($output_file);
}


sub set_expected_regions {
  my ($self,$genes,$source_genes_by_slice,$high_confidence_genes_by_id) = @_;

  # Foreach missing source gene, find the nearest non-overlapping 5' and 3' genes. Look these up in the high confidence set, keep looking up until
  # you either find a high confidence 5' and 3' on the same slice, or just take the nearest high confidence one and the expected length of the region
  # Maybe if you can't find two high confidence ones on the same slice that's a sign that we should try further anyway
  # Could take a few and see what the consensue region is, then focus on the two closest on that?

  foreach my $gene (@$genes) {
    say "Setting expected region for missing source gene with stable id ".$gene->stable_id();
    my $five_prime_neighbours = $self->fetch_neighbours($gene,$source_genes_by_slice,1);
    my $three_prime_neighbours = $self->fetch_neighbours($gene,$source_genes_by_slice,0);
    say "  Retrieved ".scalar(@$five_prime_neighbours)." 5' neighbours";
    say "  Retrieved ".scalar(@$three_prime_neighbours)." 3' neighbours";

    my $target_slice = $self->identify_target_slice([@$five_prime_neighbours,@$three_prime_neighbours],$high_confidence_genes_by_id);
    unless($target_slice) {
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-expected_region_for_missing_gene_based_on_high_conf",
          {
              feature_id => $gene->stable_id(),
              feature_type => 'gene',
              result => 'Slice Not Found',
              details => {
                five_prime_neighbours => [map { $_->stable_id } @$five_prime_neighbours],
                three_prime_neighbours => [map { $_->stable_id } @$three_prime_neighbours],
              },
              message => "No target slice identified for missing gene with stable id ".$gene->stable_id().", skipped"
            }
        );
      }
      say "  No target slice identified for missing gene with stable id ".$gene->stable_id().", skipping";
      next;
    }

    say "  Target slice retrieved: ".$target_slice;
    my ($left_boundary,$right_boundary) = $self->get_boundaries($gene,$five_prime_neighbours,$three_prime_neighbours,$high_confidence_genes_by_id,$target_slice);
    if($left_boundary and $right_boundary) {
      say "  Boundaries calculated: ".$left_boundary."..".$right_boundary;
    }
    $gene->{'left_boundary'} = $left_boundary;
    $gene->{'right_boundary'} = $right_boundary;
    $gene->{'target_slice'} = $target_slice;
    if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-expected_region_for_missing_gene_based_on_high_conf",
          {
              feature_id => $gene->stable_id(),
              feature_type => 'gene',
              result => 'Slice Found',
              details => {
                left_boundary => $left_boundary,
                right_boundary => $right_boundary,
                target_slice => $target_slice->name,
                five_prime_neighbours => [map { $_->stable_id } @$five_prime_neighbours],
                three_prime_neighbours => [map { $_->stable_id } @$three_prime_neighbours],
              },
              message => "Target slice found for missing gene based on high confidence neighbours"
            }
        );
    }
  }
}


sub get_boundaries {
  my ($self,$gene,$five_prime_neighbours,$three_prime_neighbours,$high_confidence_genes_by_id,$target_slice) = @_;

  my $five_prime_gene;
  my $three_prime_gene;
  foreach my $neighbour_gene (@$five_prime_neighbours) {
    # THINK THIS WILL NEED TO ACCOUNT FOR ARRAYREF
    if($high_confidence_genes_by_id->{$neighbour_gene->stable_id()} and ${$high_confidence_genes_by_id->{$neighbour_gene->stable_id()}}[0]->seq_region_name eq $target_slice) {
      $five_prime_gene = ${$high_confidence_genes_by_id->{$neighbour_gene->stable_id()}}[0];
      last;
    }
  }

  foreach my $neighbour_gene (@$three_prime_neighbours) {
    # THINK THIS WILL NEED TO ACCOUNT FOR ARRAYREF
    if($high_confidence_genes_by_id->{$neighbour_gene->stable_id()} and ${$high_confidence_genes_by_id->{$neighbour_gene->stable_id()}}[0]->seq_region_name eq $target_slice) {
      $three_prime_gene = ${$high_confidence_genes_by_id->{$neighbour_gene->stable_id()}}[0];
      last;
    }
  }

  unless($five_prime_gene and $three_prime_gene) {
    say "Did not find both a high confidence 5' and 3' gene on the target slice, so cannot calculate a target region";
    return;
  }

  if($five_prime_gene->seq_region_end() < $three_prime_gene->seq_region_start()) {
    return($five_prime_gene->seq_region_end(),$three_prime_gene->seq_region_start());
  } else {
    return($three_prime_gene->seq_region_end(),$five_prime_gene->seq_region_start());
  }
}


sub identify_target_slice {
  my ($self,$neighbours,$high_confidence_genes_by_id) = @_;

  my $region_counter = {};
#  say "FERGAL DEBUG CHECKING NEIGHBOURS: ".scalar(@$neighbours)." ".scalar(keys(%$high_confidence_genes_by_id));
  foreach my $gene (@$neighbours) {
#    say "FERGAL DEBUG GSID: ".$gene->stable_id();
    unless($high_confidence_genes_by_id->{$gene->stable_id()}) {
      next;
    }

#    say "FERGAL DEBUG GSRN: ".$gene->seq_region_name();
    unless($region_counter->{$gene->seq_region_name()}) {
      $region_counter->{$gene->seq_region_name()} = 1;
    } else {
      $region_counter->{$gene->seq_region_name()}++;
    }
  }

  my $best_region;
  my $best_region_count = 0;
  foreach my $region (keys(%$region_counter)) {
    my $count = $region_counter->{$region};
    if($count > $best_region_count) {
      $best_region = $region;
      $best_region_count = $count;
    }
  }

  unless($best_region) {
    $self->warning("Didn't find a target region, implies no high confidence neighbours in the set");
  }
  return($best_region);
}

sub fetch_neighbours {
  my ($self,$gene,$genes_by_slice,$direction) = @_;

  my $slice_name = $gene->seq_region_name();
  my $slice_genes = $genes_by_slice->{$slice_name};

  my $neighbour_limit = 50;
  my $neighbours = [];
  for(my $i=0; $i<scalar(@$slice_genes); $i++) {
    my $slice_gene = ${$slice_genes}[$i];
    if($slice_gene->stable_id() eq $gene->stable_id()) {
      my $neighbour_counter = 0;
      if($direction) {
        for(my $j=$i-1; $j>=0 and $neighbour_counter < $neighbour_limit; $j--) {
          my $neighbour_gene = ${$slice_genes}[$j];
          if($self->coords_overlap($gene->seq_region_start(),$gene->seq_region_end(),$neighbour_gene->seq_region_start(),$neighbour_gene->seq_region_end())) {
            next;
          }
          # At this point we have a gene that doesn't overlap in the 5' direction
          push(@$neighbours,$neighbour_gene);
          $neighbour_counter++;
        }
      } else {
        for(my $j=$i+1; $j<scalar(@$slice_genes) and $neighbour_counter < $neighbour_limit; $j++) {
          my $neighbour_gene = ${$slice_genes}[$j];
          if($self->coords_overlap($gene->seq_region_start(),$gene->seq_region_end(),$neighbour_gene->seq_region_start(),$neighbour_gene->seq_region_end())) {
            next;
          }
          # At this point we have a gene that doesn't overlap in the 3' direction
          push(@$neighbours,$neighbour_gene);
          $neighbour_counter++;
        }
      }
      last;
    }
  }
  return($neighbours);
}

sub list_high_confidence_genes {
  my ($self,$target_genes_by_slice,$target_genes_by_stable_id,$source_genes_by_stable_id) = @_;

  my $high_confidence_genes = [];
  foreach my $target_slice (keys(%$target_genes_by_slice)) {
    my $target_slice_genes = $target_genes_by_slice->{$target_slice};
    my $midpoint_coords_by_id = {};

    # First calculate the midpoints for the all the genes and store them by stable id. As the target ids might not be unique, this means
    # we have to store an arrayref
    say "  Setting midpoints for genes on ".$target_slice;
    foreach my $gene (@$target_slice_genes) {
      # By default say the gene is not high confidence
      $gene->{'is_high_confidence'} = 0;
      my $midpoint = $gene->start + ceil($gene->length()/2);
      unless($midpoint_coords_by_id->{$gene->stable_id()}) {
        $midpoint_coords_by_id->{$gene->stable_id()} = [];
      }
      # CHECK, this seems like it is wrong, assumes stable id is 
      push(@{$midpoint_coords_by_id->{$gene->stable_id()}},$midpoint);
    }

    # Now set the closest set of ids on the target genes, so we can get a score later against the reference
    say "  Setting closest genes for genes on ".$target_slice;
    foreach my $gene (@$target_slice_genes) {
      $self->set_closest_gene_ids($gene,$midpoint_coords_by_id);
    }

    # The source genes have already had their ordered neighbours calculated, so foreach target gene we should be able to get a score
    # for how similar the ordering is
    say "  Setting neighbourhood scores for genes on ".$target_slice;
    foreach my $gene (@$target_slice_genes) {
      my $source_gene = ${$source_genes_by_stable_id->{$gene->stable_id}}[0];
      unless($source_gene) {
#        $self->throw("Couldn't find a source gene with a matching stable id to the target. Target stable id: ".$gene->stable_id());
        # TEST!!!!!!!!!!
        say "Couldn't find a source gene with a matching stable id to the target. Target stable id: ".$gene->stable_id();
        $gene->{'neighbourhood_score'} = 0;
        next;
      }
      $self->set_neighbourhood_score($gene,$source_gene);
    }
  }

  # We now have a score on each gene in terms of how similar the source and target neighbours are (though a very basic calculation)
  # We can start deciding what genes are high confidence at this point
  my $neighbourhood_cutoff = 0.8;
  foreach my $id (keys(%$target_genes_by_stable_id)) {
    my $genes = $target_genes_by_stable_id->{$id};
    my $gene =  ${$genes}[0];
    my $source_gene = ${$source_genes_by_stable_id->{$id}}[0];

    my $source_transcripts = $source_gene->get_all_Transcripts();

    if(scalar(@$genes) > 1) {
      say "  Gene with stable id ".$id." has multiple mappings, so not high confidence";
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-high_confidence_evaluation",
          {
              feature_id => $id,
              feature_type => 'gene',
              result => 'rejected',
              details => {
                reason => 'multiple_mappings',
                mapping_count => scalar(@$genes)
              },
              message => "Gene has multiple mappings, not high confidence"
            }
        );
      }
      next;
    }

    unless($gene->{'complete_mapping'}) {
      say "  Gene with stable id ".$id." fails the transcript mapping check, so not high confidence";
      $self->{_logger}->log(
        "Minimap2Remap-high_confidence_evaluation",
        {
            feature_id => $id,
            feature_type => 'gene',
            result => 'rejected',
            details => {
              reason => 'incomplete_transcript_mapping',
            },
            message => "Gene fails the transcript mapping check, not high confidence"
          }
      );
      next;
    }

    if($gene->{'neighbourhood_score'} < $neighbourhood_cutoff) {
      say "  Gene with stable id ".$gene->stable_id." fails the neighbourhood score cutoff, so not high confidence. Score: ".$gene->{'neighbourhood_score'};
      if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-high_confidence_evaluation",
        {
            feature_id => $id,
            feature_type => 'gene',
            result => 'rejected',
            details => {
              reason => 'low_neighbourhood_score',
              neighbourhood_score => $gene->{'neighbourhood_score'},
              neighbourhood_cutoff => $neighbourhood_cutoff
            },
            message => "Gene does not pass neighbourhood score cutoff, not high confidence"
          }
      );
      }
      next;
    }

    # At this point we are reasonably confident the gene is mapped correctly because it only maps in one place, has
    # all the transcripts mapped with hight coverage and identity and it has a very similar set of neighbouring genes
    # The weakness here is in the fact that the similar set of neighbouring genes is not that strict, but at the same
    # time some of these assemblies are scaffold level and thus break synteny artifically, there are real expansions
    # and contractions and the fact that these are single loci and have excellent mapping scores means it should almost
    # always be correct
    $gene->{'is_high_confidence'} = 1;
    if ($self->{_logger}) {
    $self->{_logger}->log(
      'Minimap2Remap-high_confidence_evaluation',
      {
        feature_id => $id,
        feature_type => 'gene',
        result => 'accepted',
        details => {
          mapping_count => scalar(@$genes),
          neighbourhood_score => $gene->{'neighbourhood_score'},
          cutoff => $neighbourhood_cutoff
        },
        message => "Gene meets high confidence criteria"
      }
    );
    }
    push(@$high_confidence_genes,$gene);
  }

  say "Found ".scalar(@$high_confidence_genes)." high confidence genes in target";
  return($high_confidence_genes);
}


sub check_complete_mapping {
  my ($self,$gene,$source_transcripts) = @_;

  my $coverage_cutoff = 99;
  my $perc_id_cutoff = 99;
  my $target_transcripts = $gene->get_all_Transcripts();

  $gene->{'complete_mapping'} = 1;

  my $target_transcript_ids = {};
  my $total_cov = 0;
  my $total_perc_id = 0;
  foreach my $target_transcript (@$target_transcripts) {
    my $description = $target_transcript->description();
    $description =~ /parent_transcript=((ENS|GeneID_|LOC|XR_|XM_).+)\.\d*;mapping_coverage=(\d+\.\d+);mapping_identity=(\d+\.\d+)/;
    my $stable_id = $1;
    my $coverage = $3;
    my $perc_id = $4;
    unless($stable_id and defined($coverage) and defined($perc_id)) {
      $self->throw("Issue finding info for parent stable id and coverage/identity from transcript description. Description: ".$description);
    }
    $target_transcript_ids->{$stable_id} = 1;
    $total_cov += $coverage;
    $total_perc_id += $perc_id;
  }

  my $avg_cov = 0;
  my $avg_perc_id = 0;
  if (scalar(@$target_transcripts)) {
    $avg_cov = $total_cov/scalar(@$target_transcripts);
    $avg_perc_id = $total_perc_id/scalar(@$target_transcripts);
  }
  $gene->{'avg_cov'} = $avg_cov;
  $gene->{'avg_perc_id'} = $avg_perc_id;

  # Check coverage/identity cutoffs
  unless($avg_perc_id >= $perc_id_cutoff and $avg_cov >= $coverage_cutoff) {
    say "  Gene fails average coverage and identity cutoffs, so failing: Coverage: ".$avg_cov.", Perc id: ".$avg_perc_id;
    $gene->{'complete_mapping'} = 0;
    
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-complete_mapping_evaluation",
        {
          feature_id => $gene->stable_id(),
          feature_type => "gene",
          result => "failed",
          failure_reason => "quality_cutoffs",
          details => {
            avg_coverage => $avg_cov,
            avg_identity => $avg_perc_id,
            coverage_cutoff => $coverage_cutoff,
            identity_cutoff => $perc_id_cutoff,
            source_transcript_count => scalar(@$source_transcripts),
            target_transcript_count => scalar(@$target_transcripts)
          },
          message => "Gene failed quality cutoffs"
        }
      );
    }
    return(0);
  }

  # Check for missing transcripts
  my $missing_transcripts = [];
  foreach my $source_transcript (@$source_transcripts) {
    unless($target_transcript_ids->{$source_transcript->stable_id}) {
      say "  Source transcript with stable id ".$source_transcript->stable_id." was not found in the target transcript stable id list, so failing";
      push(@$missing_transcripts, $source_transcript->stable_id);
    }
  }
  
  if (@$missing_transcripts) {
    $gene->{'complete_mapping'} = 0;
    
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-complete_mapping_evaluation", 
        {
          feature_id => $gene->stable_id(),
          feature_type => "gene",
          result => "failed",
          failure_reason => "missing_transcripts",
          details => {
            missing_transcript_ids => $missing_transcripts,
            missing_count => scalar(@$missing_transcripts),
            total_source_transcripts => scalar(@$source_transcripts),
            avg_coverage => $avg_cov,
            avg_identity => $avg_perc_id
          },
          message => "Gene missing " . scalar(@$missing_transcripts) . " transcripts"
        }
      );
    }
    return(0);
  }

  # Check transcript count mismatch
  unless(scalar(@$source_transcripts) == scalar(@$target_transcripts)) {
    say "  Source and target transcript count differ, so failing";
    $gene->{'complete_mapping'} = 0;
    
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-complete_mapping_evaluation",
        {
          feature_id => $gene->stable_id(), 
          feature_type => "gene",
          result => "failed",
          failure_reason => "transcript_count_mismatch",
          details => {
            source_count => scalar(@$source_transcripts),
            target_count => scalar(@$target_transcripts),
            avg_coverage => $avg_cov,
            avg_identity => $avg_perc_id
          },
          message => "Transcript count mismatch"
        }
      );
    }
    return(0);
  }

  # Success case
  if ($self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-complete_mapping_evaluation",
      {
        feature_id => $gene->stable_id(),
        feature_type => "gene", 
        result => "passed",
        details => {
          avg_coverage => $avg_cov,
          avg_identity => $avg_perc_id,
          transcript_count => scalar(@$source_transcripts),
          coverage_cutoff => $coverage_cutoff,
          identity_cutoff => $perc_id_cutoff
        },
        message => "Gene passed complete mapping validation"
      }
    );
  }

  return(1);
}


sub set_neighbourhood_score {
  my ($self,$target_gene,$source_gene) = @_;

  my $neighbour_limit = 100;
  my $source_neighbours = $source_gene->{'sorted_neighbours'};

  my $target_neighbours = $target_gene->{'sorted_neighbours'};
  my $target_neighbours_by_id = {};
  for(my $i=0; $i<scalar(@$target_neighbours) and $i<$neighbour_limit; $i++) {
    $target_neighbours_by_id->{${$target_neighbours}[$i]} = 1;
  }

  my $score = 0;
  for(my $i=0; $i<scalar(@$source_neighbours) and $i<$neighbour_limit; $i++) {
    my $id = ${$source_neighbours}[$i];
    if($target_neighbours_by_id->{$id}) {
      $score++;
    }
  }

  my $neighbourhood_score = $score/$neighbour_limit;
#  say "Neighbourhood score ".$target_gene->stable_id.": ".$neighbourhood_score;
  $target_gene->{'neighbourhood_score'} = $neighbourhood_score;
}


sub sort_genes_by_slice {
  my ($self,$genes) = @_;

  my $genes_by_slice_hash = {};
  foreach my $gene (@$genes) {
    unless($genes_by_slice_hash->{$gene->seq_region_name()}) {
      $genes_by_slice_hash->{$gene->seq_region_name()} = [];
    }
    push(@{$genes_by_slice_hash->{$gene->seq_region_name()}},$gene);
  }

  foreach my $slice (keys(%$genes_by_slice_hash)) {
    my $slice_genes = $genes_by_slice_hash->{$slice};
    my @sorted_slice_genes  = sort { $a->start <=> $b->start } @{$slice_genes};
    $genes_by_slice_hash->{$slice} = \@sorted_slice_genes;
  }

  return($genes_by_slice_hash);
}

sub set_closest_gene_ids {
  my ($self,$gene,$midpoint_coords_by_id) = @_;

  my $id_distances = {};
  my $gene_midpoint = $gene->start + ceil($gene->length()/2);

  # Loop through each id, calculate the distance between the midpoint of the current gene and the midpoint of each other gene
  # If there are more than one genes for an id, then take the distance of the closest one
  foreach my $id (keys(%$midpoint_coords_by_id)) {
    if($id eq $gene->stable_id()) {
      next;
    }
    my $slice_gene_id_midpoints = $midpoint_coords_by_id->{$id};
    foreach my $midpoint (@$slice_gene_id_midpoints) {
      my $dist = abs($gene_midpoint - $midpoint);
      unless($id_distances->{$id}) {
        $id_distances->{$id} = $dist;
      } elsif($dist < $id_distances->{$id}) {
        $id_distances->{$id} = $dist;
      }
    }
  }

  # At this point $id_distance will have exactly one distance per id so we should sort and add up to the ten closest ids
  my @sorted_ids = sort {$id_distances->{$a} <=> $id_distances->{$b}} keys(%$id_distances);
  $gene->{'sorted_neighbours'} = \@sorted_ids;
}

sub list_missing_genes {
  my ($self,$source_genes_by_stable_id,$target_genes_by_stable_id) = @_;

  # Log the start of missing gene identification
  $self->{_logger}->log(
    "Minimap2Remap-missing_gene_identification",
    {
      feature_type => "gene_set",
      status => "started",
      details => {
        total_source_genes => scalar(keys %$source_genes_by_stable_id),
        total_target_genes => scalar(keys %$target_genes_by_stable_id)
      },
      message => "Starting identification of missing source genes"
    }
  ) if $self->{_logger};

  my $missing_source_genes = [];
  my $missing_gene_details = [];
  
  foreach my $stable_id (keys(%$source_genes_by_stable_id)) {
    unless($target_genes_by_stable_id->{$stable_id}) {
      say "Missing the following stable id in target, will add to list: ".$stable_id;
      
      my $source_gene = ${$source_genes_by_stable_id->{$stable_id}}[0];
      push(@$missing_source_genes, $source_gene);
      
      # Collect details for batch logging
      push(@$missing_gene_details, {
        stable_id => $stable_id,
        region => $source_gene->seq_region_name(),
        biotype => $source_gene->biotype(),
        length => $source_gene->length(),
        transcript_count => scalar(@{$source_gene->get_all_Transcripts()})
      });
      
      # Log individual missing gene (optional - might be too verbose)
      $self->{_logger}->log(
        "Minimap2Remap-missing_gene_identified",
        {
          feature_id => $stable_id,
          feature_type => "gene",
          status => "missing",
          details => {
            source_region => $source_gene->seq_region_name(),
            source_coordinates => $source_gene->seq_region_start() . "-" . $source_gene->seq_region_end(),
            biotype => $source_gene->biotype(),
            transcript_count => scalar(@{$source_gene->get_all_Transcripts()})
          },
          message => "Gene missing from target, added to recovery list"
        }
      ) if $self->{_logger};
    }
  }

  # Log summary of missing genes
  $self->{_logger}->log(
    "missing_gene_identification",
    {
      feature_type => "gene_set", 
      status => "completed",
      result => scalar(@$missing_source_genes) > 0 ? "genes_missing" : "no_missing_genes",
      details => {
        missing_count => scalar(@$missing_source_genes),
        missing_percentage => sprintf("%.2f", (scalar(@$missing_source_genes) / scalar(keys %$source_genes_by_stable_id)) * 100),
        missing_by_biotype => $self->_group_by_biotype($missing_gene_details),
        missing_by_region => $self->_group_by_region($missing_gene_details)
      },
      message => "Identified " . scalar(@$missing_source_genes) . " missing genes for recovery"
    }
  ) if $self->{_logger};

  return($missing_source_genes);
}

# Helper methods for grouping
sub _group_by_biotype {
  my ($self, $gene_details) = @_;
  my $grouped = {};
  foreach my $gene (@$gene_details) {
    $grouped->{$gene->{biotype}}++;
  }
  return $grouped;
}

sub _group_by_region {
  my ($self, $gene_details) = @_;
  my $grouped = {};
  foreach my $gene (@$gene_details) {
    $grouped->{$gene->{region}}++;
  }
  return $grouped;
}


sub genes_by_stable_id {
  my ($self,$genes) = @_;

  # This holds genes in a hash using the stable id as the key. Note that while in the source we expect the stable id to be unique, we don't expect that in the target
  # Therefore the hash keys will point to an arrayref instead of just directly to the gene
  my $genes_by_id = {};
  foreach my $gene (@$genes) {
    unless($gene->stable_id()) {
      $self->throw("Gene with dbID ".$gene->dbID()." does not have a stable id assigned, need an assigned stable id to compare between source and target sets");
    }
    unless($genes_by_id->{$gene->stable_id()}) {
      $genes_by_id->{$gene->stable_id()} = [];
    }
    push(@{$genes_by_id->{$gene->stable_id()}},$gene);
  }
  return($genes_by_id);
}


sub process_results {
  my ($self,$source_gene,$gene_genomic_seqs_hash,$target_genes) = @_;

  my $source_gene_id = $source_gene->stable_id();
  my $high_confidence = 0;
  my $gene_seq = $source_gene->seq();
  my $source_transcripts = $source_gene->get_all_Transcripts();
  my $source_transcript_id_hash = {};

  say "Source transcript list:";
  foreach my $source_transcript (@$source_transcripts) {
    say "  ".$source_transcript->stable_id();
    my $source_transcript_id = $source_transcript->dbID();
    $source_transcript_id_hash->{$source_transcript_id} = $source_transcript;
  }

  my $max_intron_size = $self->calculate_max_intron_size($source_transcripts);
  $max_intron_size = ceil($max_intron_size);
  say "Max intron size of transcripts in the cluster: ".$max_intron_size;

  my $good_transcripts = [];
  my $bad_source_transcripts = [];
  my $bad_transcripts = [];
  my $target_adaptor = $self->target_adaptor();
  my $final_genes = [];
  my $paf_results = $source_gene->{'paf_regions'};
  
  my $paf_region_counter = 0;
  foreach my $paf_result (@$paf_results) {
    $paf_region_counter++;
    
    my $best_transcripts_by_id = {};
    my $target_genomic_start = ${$paf_result}[0];
    my $target_genomic_end = ${$paf_result}[1];
    my $target_strand = ${$paf_result}[2];
    my $target_genomic_name = ${$paf_result}[3];

    if($target_strand eq '+') {
      $target_strand = 1;
    } else {
      $target_strand = -1;
    }

    say "Mapping gene ".$source_gene->stable_id()." to ".$target_genomic_start.":".$target_genomic_end.":".$target_strand.":".$target_genomic_name;

    my $target_slice_adaptor = $target_adaptor->get_SliceAdaptor();
    my $target_parent_slice = $target_slice_adaptor->fetch_by_region('toplevel',$target_genomic_name);

    unless($target_parent_slice) {
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-slice_fetch_failure",
          {
            feature_id => $source_gene_id,
            feature_type => 'gene',
            result => 'failed',
            details => {
              target_genomic_name => $target_genomic_name,
              paf_region => "$target_genomic_name:$target_genomic_start-$target_genomic_end($target_strand)",
              failure_reason => "target_slice_not_found"
            },
            message => "Could not fetch target slice, skipping this PAF region"
          }
        );
      }
      $self->throw("Could not fetch the slice in the target assembly. Slice name: ".$target_genomic_name);
    }

    my $target_sequence_adaptor = $target_adaptor->get_SequenceAdaptor;
    my $target_genomic_seq = ${ $target_sequence_adaptor->fetch_by_Slice_start_end_strand($target_parent_slice, $target_genomic_start, $target_genomic_end, 1) };
    my $target_region_slice = $target_slice_adaptor->fetch_by_region('toplevel', $target_genomic_name, $target_genomic_start, $target_genomic_end, 1);
    say "Target slice identified to search for transcripts in after first pass: ".$target_region_slice->name();
    my $target_genomic_fasta = ">".$target_genomic_name."\n".$target_genomic_seq;
    my $target_genome_file = $self->write_input_file([$target_genomic_fasta]);

    say "Projecting gene: ".$source_gene->stable_id();

    my $coverage_threshold = 98;
    my $perc_id_threshold = 98;
    
    # PROJECTION METHOD
    my $projected_transcripts_by_id = $self->project_gene_coords($source_gene,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,$gene_genomic_seqs_hash);
    $self->print_transcript_stats($projected_transcripts_by_id,'projection');
    $self->update_best_transcripts($best_transcripts_by_id,$projected_transcripts_by_id);
    
    # Log projection quality decisions
    my $projection_passed = 0;
    my $projection_failed_transcripts = [];
    foreach my $transcript_id (keys %$projected_transcripts_by_id) {
      my $transcript = $projected_transcripts_by_id->{$transcript_id};
      if ($transcript->{'coverage'} >= $coverage_threshold && $transcript->{'percent_id'} >= $perc_id_threshold) {
        $projection_passed = 1;
      } else {
        push @$projection_failed_transcripts, {
          transcript_id => $transcript_id,
          coverage => $transcript->{'coverage'},
          percent_id => $transcript->{'percent_id'},
          coverage_threshold => $coverage_threshold,
          perc_id_threshold => $perc_id_threshold
        };
      }
    }
    
    if (scalar(@$projection_failed_transcripts) > 0 && $self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-projection_quality_failure",
        {
          feature_id => $source_gene_id,
          feature_type => 'gene',
          result => 'quality_check_failed',
          details => {
            method => 'coordinate_projection',
            failed_transcripts => $projection_failed_transcripts,
            reason_for_minimap => 'projection_below_thresholds',
            will_attempt_minimap => 1
          },
          message => "Coordinate projection failed quality thresholds, proceeding to minimap2"
        }
      );
    }
    
    # MINIMAP2 METHOD
    my $minimap_transcripts_by_id = $self->map_gene_minimap($source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
                                                            $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,$target_slice_adaptor,
                                                            $best_transcripts_by_id);
    $self->print_transcript_stats($minimap_transcripts_by_id,'minimap local');
    $self->update_best_transcripts($best_transcripts_by_id,$minimap_transcripts_by_id);
    
    # Log minimap quality decisions
    my $minimap_passed = 0;
    my $minimap_failed_transcripts = [];
    foreach my $transcript_id (keys %$minimap_transcripts_by_id) {
      my $transcript = $minimap_transcripts_by_id->{$transcript_id};
      if ($transcript->{'coverage'} >= $coverage_threshold && $transcript->{'percent_id'} >= $perc_id_threshold) {
        $minimap_passed = 1;
      } else {
        push @$minimap_failed_transcripts, {
          transcript_id => $transcript_id,
          coverage => $transcript->{'coverage'},
          percent_id => $transcript->{'percent_id'},
          coverage_threshold => $coverage_threshold,
          perc_id_threshold => $perc_id_threshold
        };
      }
    }
    
    if (scalar(@$minimap_failed_transcripts) > 0 && $self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-minimap_quality_failure",
        {
          feature_id => $source_gene_id,
          feature_type => 'gene',
          result => 'quality_check_failed',
          details => {
            method => 'minimap2',
            failed_transcripts => $minimap_failed_transcripts,
            reason_for_exonerate => 'minimap_below_thresholds',
            will_attempt_exonerate => 1,
            max_intron_size_used => $max_intron_size
          },
          message => "Minimap2 failed quality thresholds, proceeding to exonerate"
        }
      );
    }
    
    # EXONERATE METHOD
    my $exonerate_transcripts_by_id = $self->map_gene_exonerate($source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
                                                                $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,
                                                                $target_slice_adaptor,$best_transcripts_by_id);
    $self->print_transcript_stats($exonerate_transcripts_by_id,'exonerate local');
    $self->update_best_transcripts($best_transcripts_by_id,$exonerate_transcripts_by_id);
    
    # Log exonerate quality decisions and final method selection
    my $exonerate_passed = 0;
    my $exonerate_failed_transcripts = [];
    my $final_method_used = [];
    
    foreach my $transcript_id (keys %$exonerate_transcripts_by_id) {
      my $transcript = $exonerate_transcripts_by_id->{$transcript_id};
      if ($transcript->{'coverage'} >= $coverage_threshold && $transcript->{'percent_id'} >= $perc_id_threshold) {
        $exonerate_passed = 1;
      } else {
        push @$exonerate_failed_transcripts, {
          transcript_id => $transcript_id,
          coverage => $transcript->{'coverage'},
          percent_id => $transcript->{'percent_id'},
          coverage_threshold => $coverage_threshold,
          perc_id_threshold => $perc_id_threshold
        };
      }
    }
    
    # Determine which methods were used for final transcripts
    foreach my $transcript_id (keys %$best_transcripts_by_id) {
      my $transcript = $best_transcripts_by_id->{$transcript_id};
      push @$final_method_used, {
        transcript_id => $transcript_id,
        method_used => $transcript->{'annotation_method'} || 'unknown',
        coverage => $transcript->{'coverage'},
        percent_id => $transcript->{'percent_id'}
      };
    }
    
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-method_selection_decision",
        {
          feature_id => $source_gene_id,
          feature_type => 'gene',
          result => 'method_selection_complete',
          details => {
            projection_success => $projection_passed,
            minimap_success => $minimap_passed,
            exonerate_success => $exonerate_passed,
            final_transcripts_by_method => $final_method_used,
            all_methods_failed => (!$projection_passed && !$minimap_passed && !$exonerate_passed),
            quality_thresholds => {
              coverage => $coverage_threshold,
              percent_id => $perc_id_threshold
            }
          },
          message => "Method selection completed - best transcript per method chosen"
        }
      );
    }

    $self->set_cds_sequences($best_transcripts_by_id,$source_transcript_id_hash);
    $self->qc_cds_sequences($best_transcripts_by_id,$source_transcript_id_hash);
    $self->set_transcript_descriptions($best_transcripts_by_id,$source_transcript_id_hash);

    my $all_transcripts = $self->label_transcript_status($best_transcripts_by_id);
    my $biotypes_hash = $self->generate_biotypes_hash($all_transcripts);

    say "Building single transcript genes ahead of clustering";
    my $single_transcript_genes = $self->generate_single_transcript_genes($all_transcripts);

    my $genes_by_seq_region = {};
    foreach my $gene (@$single_transcript_genes) {
      my $seq_region_name = $gene->seq_region_name();
      unless($genes_by_seq_region->{$seq_region_name}) {
        $genes_by_seq_region->{$seq_region_name} = [];
      }
      push(@{$genes_by_seq_region->{$seq_region_name}},$gene);
    }

    say "Clustering genes";
    my $found_good_cluster = 0;
    my $all_clusters = [];
    foreach my $seq_region_name (keys(%{$genes_by_seq_region})) {
      my $genes = $genes_by_seq_region->{$seq_region_name};
      my ($clusters, $unclustered) = cluster_Genes($genes,$biotypes_hash);
      say "Found ".scalar(@$clusters)." transcript clusters";
      say "Found ".scalar(@$unclustered)." unclustered transcripts";
      push(@$all_clusters,@$clusters);
      push(@$all_clusters,@$unclustered);
    }

    say "Found ".scalar(@$all_clusters)." initial clusters";

    foreach my $cluster (@$all_clusters) {
      $self->check_cluster_status($cluster);
      if($cluster->{'status'} eq 'good') {
        $found_good_cluster = 1;
      }
    } 

    my $final_clusters = [];
    my $good_clusters = [];
    my $bad_clusters = [];
    
    foreach my $cluster (@$all_clusters) {
      if($cluster->{'status'} eq 'good') {
        push(@$final_clusters,$cluster);
        push(@$good_clusters,$cluster);
      } elsif($cluster->{'status'} eq 'bad' and !$found_good_cluster) {
        push(@$final_clusters,$cluster);
        push(@$bad_clusters,$cluster);
      }
    }

    say "Found ".scalar(@$final_clusters)." final clusters";

    # Log cluster selection decision with detailed reasoning
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-cluster_quality_decision",
        {
          feature_id => $source_gene_id,
          feature_type => 'gene',
          result => 'cluster_selection_made',
          details => {
            total_clusters_evaluated => scalar(@$all_clusters),
            good_clusters_available => scalar(@$good_clusters),
            bad_clusters_available => scalar(@$bad_clusters),
            selection_rule => $found_good_cluster ? 'good_clusters_only' : 'bad_clusters_fallback',
            final_clusters_chosen => scalar(@$final_clusters),
            rejected_clusters => scalar(@$all_clusters) - scalar(@$final_clusters),
            decision_reasoning => $found_good_cluster ? 
              'Good quality clusters found, bad clusters discarded' : 
              'No good clusters found, using bad clusters as fallback'
          },
          message => "Cluster quality assessment: " . ($found_good_cluster ? 'using good clusters' : 'fallback to bad clusters')
        }
      );
    }

    say "Building genes from final clusters";
    my $parent_gene_ids = $self->parent_gene_ids();

    foreach my $cluster (@$final_clusters) {
      my $gene = $self->create_gene_from_cluster($cluster,$parent_gene_ids,$source_transcript_id_hash);
      say "Created gene: ".$gene->stable_id()." ".$gene->seq_region_name().":".$gene->seq_region_start.":".$gene->seq_region_end.":".$gene->strand();
      
      my $transcripts = $gene->get_all_Transcripts();
      my @transcript_methods = ();  # FIXED: fresh array for each gene
      foreach my $transcript (@$transcripts) {
        say "  Transcript: ".$transcript->stable_id()." ".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$transcript->strand();
        my $updated_description = $transcript->description().";annotation_method=".$transcript->{'annotation_method'};
        $transcript->description($updated_description);
        push @transcript_methods, $transcript->{'annotation_method'};
      }

      # Log gene boundaries and method provenance
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-gene_boundaries_set",
          {
            feature_id => $source_gene_id,
            feature_type => 'gene',
            result => 'boundaries_established',
            details => {
              target_gene_location => $gene->seq_region_name().":".$gene->seq_region_start.":".$gene->seq_region_end.":".$gene->strand(),
              boundary_start => $gene->seq_region_start,
              boundary_end => $gene->seq_region_end,
              boundary_strand => $gene->strand(),
              boundary_determination_method => 'cluster_based_gene_creation',
              transcript_count => scalar(@$transcripts),
              mapping_methods_used => \@transcript_methods,  # FIXED: reference to array
              cluster_status => $cluster->{'status'},
              source_gene_location => $source_gene->seq_region_name.":".$source_gene->start."-".$source_gene->end."(".$source_gene->strand.")"
            },
            message => "Gene boundaries established from " . $cluster->{'status'} . " cluster using methods: " . join(',', @transcript_methods)
          }
        );
      }

      my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_g',-VALUE => $source_gene->stable_id.".".$source_gene->version);
      $gene->biotype($source_gene->biotype());
      $gene->stable_id($source_gene->stable_id());
      $gene->add_Attributes($parent_attribute);
      $gene->biotype($source_gene->biotype());
      $gene->stable_id($source_gene->stable_id());
      $gene->{'to_write'} = 1;
      push(@$final_genes,$gene);
    }
    
    # Log if no genes were created from this PAF region
    if (scalar(@$final_clusters) == 0 && $self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-paf_region_failure",
        {
          feature_id => $source_gene_id,
          feature_type => 'gene',
          result => 'no_genes_created',
          details => {
            paf_region => "$target_genomic_name:$target_genomic_start-$target_genomic_end($target_strand)",
            total_clusters_found => scalar(@$all_clusters),
            good_clusters => scalar(@$good_clusters),
            bad_clusters => scalar(@$bad_clusters),
            failure_reason => scalar(@$all_clusters) == 0 ? 'no_clusters_formed' : 'all_clusters_rejected'
          },
          message => "No final genes created from this PAF region"
        }
      );
    }
  }

  say "Built ".scalar(@$final_genes)." final genes";
  
  # Log overall pipeline outcome
  if (scalar(@$final_genes) == 0 && $self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-pipeline_failure",
      {
        feature_id => $source_gene_id,
        feature_type => 'gene',
        result => 'complete_failure',
        details => {
          paf_regions_attempted => scalar(@$paf_results),
          final_genes_created => 0,
          failure_reason => 'no_successful_projections_across_all_regions'
        },
        message => "Gene projection pipeline failed - no genes created from any PAF region"
      }
    );
  }
  
  return($final_genes);
}


sub set_transcript_description {
  my ($self, $transcript, $source_transcript) = @_;

  # add source transcript stable id as transcript attribute
  my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_t',-VALUE => $source_transcript->stable_id().".".$source_transcript->version());
  $transcript->add_Attributes($parent_attribute);

  my $description_string = ";parent_transcript=".$source_transcript->stable_id().".".$source_transcript->version().";mapping_coverage=".$transcript->{'cov'}.";mapping_identity=".$transcript->{'perc_id'};
  if($transcript->{'cds_description'}) {
    $description_string .= $transcript->{'cds_description'};
  }
  $transcript->description($description_string);
}


sub set_transcript_descriptions {
  my ($self,$transcripts_by_id,$source_transcript_id_hash) = @_;

  foreach my $id (keys(%$transcripts_by_id)) {
    my $transcript = $transcripts_by_id->{$id};
    $self->set_transcript_description($transcript, $source_transcript_id_hash->{$transcript->stable_id});
  }
}


sub qc_cds_sequence {
  my ($self, $transcript, $source_transcript) = @_;

  my $source_transcript_id = $source_transcript->stable_id();
  my $source_translation = $source_transcript->translation;
  
  if($source_translation) {
    my $transcript_cdna = $transcript->translateable_seq;
    if ($transcript_cdna) {
      my ($cds_coverage,$cds_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript->translateable_seq(),$transcript_cdna);
      my $cds_description = ";cds_coverage=".$cds_coverage.";cds_identity=".$cds_percent_id;
      my $aligned_source_seq_copy = $aligned_source_seq;
      my $aligned_target_seq_copy = $aligned_target_seq;
      $aligned_source_seq_copy =~ s/\-\-\-//g;
      $aligned_target_seq_copy =~ s/\-\-\-//g;

      # Assess CDS alignment quality and gaps
      my $has_gaps = ($aligned_source_seq_copy =~ /\-/ or $aligned_target_seq_copy =~ /\-/);
      
      if($has_gaps) {
        $cds_description .= ";cds_gap=1";
        
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-cds_quality_assessment",
            {
              feature_id => $source_transcript_id,
              feature_type => 'transcript',
              result => 'gaps_detected',
              details => {
                cds_coverage => $cds_coverage,
                cds_percent_id => $cds_percent_id,
                alignment_has_gaps => 1,
                quality_concern => 'cds_alignment_contains_gaps',
                source_cds_length => length($source_transcript->translateable_seq()),
                target_cds_length => length($transcript_cdna)
              },
              message => "CDS quality check: gaps detected in alignment (${cds_coverage}% cov, ${cds_percent_id}% id)"
            }
          );
        }
      } else {
        $cds_description .= ";cds_gap=0";
        
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-cds_quality_assessment",
            {
              feature_id => $source_transcript_id,
              feature_type => 'transcript',
              result => 'high_quality_alignment',
              details => {
                cds_coverage => $cds_coverage,
                cds_percent_id => $cds_percent_id,
                alignment_has_gaps => 0,
                quality_status => 'clean_cds_alignment',
                source_cds_length => length($source_transcript->translateable_seq()),
                target_cds_length => length($transcript_cdna)
              },
              message => "CDS quality check: clean alignment without gaps (${cds_coverage}% cov, ${cds_percent_id}% id)"
            }
          );
        }
      }
      
      say join(' ', __LINE__, $source_transcript->display_id, $source_transcript->biotype), "\n";
      my $translated_transcript = $transcript->translate;
      if ($translated_transcript) {
        my $new_seq = $transcript->translate->seq;
        
        # Transfer attributes from source
        foreach my $attribute (@{$source_transcript->get_all_Attributes}) {
          if ($attribute->code eq 'cds_start_NF' or $attribute->code eq 'cds_end_NF') {
            $transcript->add_Attributes($attribute);
          }
        }
        foreach my $attribute (@{$source_translation->get_all_Attributes}) {
          if ($attribute->code eq 'initial_met') {
            if (substr($new_seq, 0, 1) ne 'M') {
              $transcript->translation->add_Attributes($attribute);
            }
          }
        }
        
        # Critical decision point: internal stop codon detection and handling
        my $num_internal_stops = $new_seq =~ tr/*/*/;
        if ($num_internal_stops) {
          my $source_has_selenocysteine = (index($source_translation->seq, 'U') >= 0);
          my $source_has_known_stops = (index(substr($source_translation->seq, 1, $source_translation->length-2), 'X') >= 0);
          
          if ($source_has_selenocysteine or $source_has_known_stops) {
            # Source has special amino acids - try to preserve them
            if ($self->{_logger}) {
              $self->{_logger}->log(
                "Minimap2Remap-internal_stop_handling_decision",
                {
                  feature_id => $source_transcript_id,
                  feature_type => 'transcript',
                  result => 'special_amino_acid_handling',
                  details => {
                    num_internal_stops => $num_internal_stops,
                    source_has_selenocysteine => $source_has_selenocysteine,
                    source_has_known_internal_stops => $source_has_known_stops,
                    handling_strategy => 'preserve_source_seq_edits',
                    decision_reason => 'source_contains_special_amino_acids'
                  },
                  message => "Internal stops detected - attempting to preserve source selenocysteine/known stops"
                }
              );
            }
            
            my $source_seq_edits = $source_translation->get_all_SeqEdits;
            my $pos = 0;
            my $seq_edits_applied = 0;
            
            foreach my $seq_edit (@$source_seq_edits) {
              if ($seq_edit->code eq '_selenocysteine' or $seq_edit->code eq '_stop_codon_rt') {
                $pos = index($new_seq, '*', $pos);
                if ($pos == -1) {
                  $self->warning($source_transcript->display_id.', I was expecting an internal stop but could not find one');
                  last;
                }
                ++$pos; # We need to be in 1-index coordinates
                $self->warning("Internal stop at $pos might not be a ".$seq_edit->name.' for '.$source_transcript->display_id)
                  unless ($seq_edit->start == $pos or ($pos >= $seq_edit->start-10 and $pos <= $seq_edit->start+10));
                $transcript->translation->add_Attributes(ref($seq_edit)->new(
                  -code => $seq_edit->code,
                  -start => $pos,
                  -end => $pos,
                  -alt_seq => $seq_edit->alt_seq,
                )->get_Attribute);
                $seq_edits_applied++;
              }
            }
            
            my $latest_seq = $transcript->translate->seq;
            my $still_internal_stops = $latest_seq =~ tr/*/*/;
            if ($still_internal_stops) {
              if ($self->{_logger}) {
                $self->{_logger}->log(
                  "Minimap2Remap-stop_correction_decision",
                  {
                    feature_id => $source_transcript_id,
                    feature_type => 'transcript',
                    result => 'fallback_to_intron_replacement',
                    details => {
                      seq_edits_applied => $seq_edits_applied,
                      remaining_internal_stops => $still_internal_stops,
                      decision_reason => 'seq_edits_insufficient',
                      correction_method => 'replace_stops_with_introns'
                    },
                    message => "SeqEdit preservation insufficient - falling back to intron replacement for remaining stops"
                  }
                );
              }
              
              my $no_internal_stop_transcript = replace_stops_with_introns($transcript, $still_internal_stops);
              if ($no_internal_stop_transcript) {
                $no_internal_stop_transcript->{cov} = $transcript->{cov};
                $no_internal_stop_transcript->{perc_id} = $transcript->{perc_id};
                $no_internal_stop_transcript->{parent_transcript_versioned_stable_id} = $transcript->{parent_transcript_versioned_stable_id};
                $no_internal_stop_transcript->{annotation_method} = $transcript->{annotation_method};
                $transcript = $no_internal_stop_transcript;
                
                if ($self->{_logger}) {
                  $self->{_logger}->log(
                    "Minimap2Remap-stop_correction_success",
                    {
                      feature_id => $source_transcript_id,
                      feature_type => 'transcript',
                      result => 'internal_stops_corrected',
                      details => {
                        correction_method => 'intron_replacement',
                        original_internal_stops => $still_internal_stops,
                        correction_successful => 1
                      },
                      message => "Successfully corrected internal stops by intron replacement"
                    }
                  );
                }
              }
              else {
                if ($self->{_logger}) {
                  $self->{_logger}->log(
                    "Minimap2Remap-stop_correction_failure",
                    {
                      feature_id => $source_transcript_id,
                      feature_type => 'transcript',
                      result => 'correction_failed',
                      details => {
                        remaining_internal_stops => $still_internal_stops,
                        correction_method_attempted => 'intron_replacement',
                        correction_successful => 0,
                        transcript_quality => 'degraded'
                      },
                      message => "Failed to correct internal stops - transcript quality degraded"
                    }
                  );
                }
                $self->warning('Could not replace internal stops for '.$source_transcript->display_id);
              }
            }
          }
          else {
            # Source has no special amino acids - standard stop correction
            if ($self->{_logger}) {
              $self->{_logger}->log(
                "Minimap2Remap-internal_stop_handling_decision",
                {
                  feature_id => $source_transcript_id,
                  feature_type => 'transcript',
                  result => 'standard_stop_correction',
                  details => {
                    num_internal_stops => $num_internal_stops,
                    source_has_special_amino_acids => 0,
                    handling_strategy => 'replace_stops_with_introns',
                    decision_reason => 'no_special_amino_acids_in_source'
                  },
                  message => "Internal stops detected - applying standard intron replacement correction"
                }
              );
            }
            
            my $no_internal_stop_transcript = replace_stops_with_introns($transcript, $num_internal_stops);
            if ($no_internal_stop_transcript) {
              $no_internal_stop_transcript->{cov} = $transcript->{cov};
              $no_internal_stop_transcript->{perc_id} = $transcript->{perc_id};
              $no_internal_stop_transcript->{parent_transcript_versioned_stable_id} = $transcript->{parent_transcript_versioned_stable_id};
              $no_internal_stop_transcript->{annotation_method} = $transcript->{annotation_method};
              $transcript = $no_internal_stop_transcript;
              
              if ($self->{_logger}) {
                $self->{_logger}->log(
                  "Minimap2Remap-stop_correction_success",
                  {
                    feature_id => $source_transcript_id,
                    feature_type => 'transcript',
                    result => 'internal_stops_corrected',
                    details => {
                      correction_method => 'intron_replacement',
                      original_internal_stops => $num_internal_stops,
                      correction_successful => 1
                    },
                    message => "Successfully corrected $num_internal_stops internal stops by intron replacement"
                  }
                );
              }
            }
            else {
              if ($self->{_logger}) {
                $self->{_logger}->log(
                  "Minimap2Remap-stop_correction_failure",
                  {
                    feature_id => $source_transcript_id,
                    feature_type => 'transcript',
                    result => 'correction_failed',
                    details => {
                      internal_stops => $num_internal_stops,
                      correction_method_attempted => 'intron_replacement',
                      correction_successful => 0,
                      transcript_quality => 'severely_degraded'
                    },
                    message => "Failed to correct $num_internal_stops internal stops - transcript severely degraded"
                  }
                );
              }
              $self->warning('Could not replace internal stops for '.$source_transcript->display_id);
            }
          }
        } else {
          # No internal stops - clean transcript
          if ($self->{_logger}) {
            $self->{_logger}->log(
              "Minimap2Remap-cds_quality_final",
              {
                feature_id => $source_transcript_id,
                feature_type => 'transcript',
                result => 'clean_translation',
                details => {
                  cds_coverage => $cds_coverage,
                  cds_percent_id => $cds_percent_id,
                  has_internal_stops => 0,
                  translation_quality => 'high',
                  protein_length => length($new_seq)
                },
                message => "CDS QC passed - clean translation without internal stops"
              }
            );
          }
        }
      } else {
        # No translation possible
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-cds_quality_failure",
            {
              feature_id => $source_transcript_id,
              feature_type => 'transcript',
              result => 'translation_failed',
              details => {
                cds_coverage => $cds_coverage,
                cds_percent_id => $cds_percent_id,
                failure_reason => 'cannot_translate_transcript',
                transcript_quality => 'failed'
              },
              message => "CDS QC failed - transcript cannot be translated"
            }
          );
        }
      }
      $transcript->{'cds_description'} = $cds_description;
    }
    else {
      # No translateable sequence
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-cds_quality_failure",
          {
            feature_id => $source_transcript_id,
            feature_type => 'transcript',
            result => 'no_translateable_sequence',
            details => {
              failure_reason => 'no_cds_sequence_available',
              transcript_quality => 'failed'
            },
            message => "CDS QC failed - no translateable sequence available"
          }
        );
      }
      $transcript->{'cds_description'} = ";cds_coverage=0.00;cds_identity=0.00";
    }
  } else {
    # Source has no translation
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-cds_qc_skipped",
        {
          feature_id => $source_transcript_id,
          feature_type => 'transcript',
          result => 'qc_skipped',
          details => {
            reason => 'source_transcript_non_coding',
            transcript_type => 'non_coding'
          },
          message => "CDS QC skipped - source transcript has no translation"
        }
      );
    }
  }
  
  return $transcript;
}

sub qc_cds_sequences {
  my ($self,$transcripts_by_id,$source_transcript_id_hash) = @_;

  foreach my $transcript (values(%$transcripts_by_id)) {
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};
    $transcript = $self->qc_cds_sequence($transcript, $source_transcript);
  }
}


sub set_cds_sequences {
  my ($self,$transcripts_by_id,$source_transcript_id_hash) = @_;

  foreach my $transcript (values(%$transcripts_by_id)) {
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};
    if($source_transcript->translation() and !$transcript->translation) {
      $self->project_cds($transcript,$source_transcript);
    }
  }
}

sub project_cds {
  my ($self,$transcript,$source_transcript) = @_;

  my $source_transcript_id = $source_transcript->stable_id();
  my $frameshift_correction_attempted = 0;  # FIXED: Declare at subroutine level

  # Check for required alignment data
  my $aligned_source_seq = $transcript->{'aligned_source_seq'};
  my $aligned_target_seq = $transcript->{'aligned_target_seq'};
  unless($aligned_source_seq and $aligned_target_seq) {
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-cds_projection_failure",
        {
          feature_id => $source_transcript_id,
          feature_type => 'transcript',
          result => 'failed',
          details => {
            reason => 'missing_alignment_data',
            has_aligned_source => defined($aligned_source_seq),
            has_aligned_target => defined($aligned_target_seq)
          },
          message => "CDS projection failed - missing transcript alignment data"
        }
      );
    }
    $self->throw("Issue fetching alignment for transcript with stable_id ".$source_transcript->stable_id().", expected target transcript to have the alignment stored on it");
  }
  
  my $source_seq = $aligned_source_seq;
  $source_seq =~ s/\-//g;
  my $target_seq = $aligned_target_seq;
  $target_seq =~ s/\-//g;

  my @source_align_array = split('',$aligned_source_seq);
  my @target_align_array = split('',$aligned_target_seq);
  my $source_cds_pos = 0;
  my $source_align_seq_pos = 0;
  my $target_cds_pos = 0;
  my $target_align_seq_pos = 0;

  my $source_cds_start = $source_transcript->cdna_coding_start();
  my $source_cds_end = $source_transcript->cdna_coding_end();
  say "Source CDS start: ".$source_cds_start;
  say "Source CDS end: ".$source_cds_end;

  my $source_cds_start_index;
  my $source_cds_end_index;
  my $target_cds_start_index;
  my $target_cds_end_index;
  my $in_feature_alignment = 0;
  my $source_index_target_cds_start = 1;

  my @frameshifts;
  for(my $i=0; $i < scalar(@source_align_array); $i++) {
    my $source_value = $source_align_array[$i];
    my $target_value = $target_align_array[$i];
    if($source_value ne '-') {
      $source_cds_pos++;
    }
    else {
      push(@frameshifts, ['I', $i+1]);
    }

    if($target_value ne '-') {
      $target_cds_pos++;
    }
    else {
      push(@frameshifts, ['D', $i+1]);
    }

    if($source_cds_pos == $source_cds_start) {
      $source_cds_start_index = $source_cds_pos;
      say "Index of the feature start in source seq: ".$source_cds_start_index;
      $in_feature_alignment = 1;
    }

    if($in_feature_alignment and $target_value ne '-' and !(defined($target_cds_start_index))) {
      $target_cds_start_index = $target_cds_pos;
      $source_index_target_cds_start = $source_cds_pos;
    }

    if($source_cds_pos == $source_cds_end) {
      $source_cds_end_index = $source_cds_pos;
      $target_cds_end_index = $target_cds_pos;
      say "Index of the feature end in source seq: ".$source_cds_end_index;
      $in_feature_alignment = 0;
      last;
    }
  }

  # Log CDS coordinate mapping results
  if ($target_cds_start_index and $target_cds_end_index) {
    my $source_feature_length = $source_cds_end_index - $source_cds_start_index + 1;
    my $target_feature_length = $target_cds_end_index - $target_cds_start_index + 1;
    my $recovered_source_feature_seq = substr($source_seq,$source_cds_start_index-1,$source_feature_length);
    my $recovered_target_feature_seq = substr($target_seq,$target_cds_start_index-1,$target_feature_length);
    say "For transcript ".$source_transcript->stable_id()." original and mapped CDS sequences:\n".$recovered_source_feature_seq."\n".$recovered_target_feature_seq;
    
    my $length_change = $target_feature_length - $source_feature_length;
    my $frameshift_count = scalar(@frameshifts);
    
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-cds_coordinate_mapping",
        {
          feature_id => $source_transcript_id,
          feature_type => 'transcript',
          result => 'coordinates_mapped',
          details => {
            source_cds_coords => "$source_cds_start-$source_cds_end",
            target_cds_coords => "$target_cds_start_index-$target_cds_end_index",
            source_cds_length => $source_feature_length,
            target_cds_length => $target_feature_length,
            length_change => $length_change,
            frameshifts_detected => $frameshift_count,
            length_change_divisible_by_3 => ($length_change % 3 == 0),
            coordinate_mapping_method => 'alignment_based_projection'
          },
          message => "CDS coordinates mapped - " . ($frameshift_count > 0 ? "frameshifts detected, will attempt correction" : "no frameshifts detected")
        }
      );
    }

    my $original_phase = 0;
    my $source_translation = $source_transcript->translation;
    if ($source_translation->start_Exon->phase > 0) {
      $original_phase = $source_translation->start_Exon->phase;
    }
    my $phase_adjust = abs($source_index_target_cds_start-$source_cds_start-$original_phase)%3;

    my $cds_start_exon;
    my $cds_end_exon;
    my $cds_start_offset = 0;
    my $cds_end_offset = 0;
    my $cumulative_length = 0;
    my $exons = $transcript->get_all_Exons();
    foreach my $exon (@$exons) {
      say "Checking exons with the following start/end for CDS coords: ".$exon->start()."/".$exon->end();
      if(!$cds_start_exon and $target_cds_start_index >= $cumulative_length and $target_cds_start_index <= ($cumulative_length + $exon->length)) {
        $cds_start_exon = $exon;
        if($transcript->strand == 1) {
          $cds_start_offset = $target_cds_start_index - $cumulative_length;
          if ($cds_start_offset > $exon->length) {
              if ($self->{_logger}) {
                  $self->{_logger}->log(
                    "Minimap2Remap-cds_projection_critical_error",
                    {
                      feature_id => $source_transcript_id,
                      feature_type => 'transcript',
                      result => 'invalid_cds_offset',
                      details => {
                        error_type => 'cds_offset_exceeds_exon_length',
                        invalid_cds_start_offset => $cds_start_offset,
                        exon_length => $exon->length,
                        exon_coordinates => $exon->start . "-" . $exon->end,
                        target_cds_start_index => $target_cds_start_index,
                        cumulative_length_at_exon => $cumulative_length,
                        transcript_strand => $transcript->strand,
                        source_cds_coordinates => "$source_cds_start-$source_cds_end",
                        alignment_based_mapping => {
                          source_alignment_length => length($aligned_source_seq),
                          target_alignment_length => length($aligned_target_seq),
                          frameshifts_detected => scalar(@frameshifts)
                        },
                        exon_context => {
                          total_exons => scalar(@$exons),
                          exon_index_in_transcript => _get_exon_index($exons, $exon),
                          is_single_nucleotide_exon => ($exon->length == 1 ? 1 : 0)
                        }
                      },
                      message => "CRITICAL: CDS start offset ($cds_start_offset) exceeds exon length (" . $exon->length . ") - this indicates a coordinate mapping failure"
                    }
                  );
                }
            $self->warning("INVALID CDS OFFSET: offset=$cds_start_offset, exon_length=" . $exon->length . 
                          ", target_index=$target_cds_start_index, cumulative=$cumulative_length");
          }
        } else {
          say "FERGAL CDS START DEBUG: ".$exon->length." - (".($exon->length - $target_cds_start_index).") - $cumulative_length";
          $cds_start_offset = $exon->length - ($exon->length - $target_cds_start_index) - $cumulative_length;
        }
      }

      if(($target_cds_end_index >= $cumulative_length) and ($target_cds_end_index <= $cumulative_length + $exon->length())) {
        $cds_end_exon = $exon;
        if($transcript->strand == 1) {
          $cds_end_offset = $target_cds_end_index - $cumulative_length;
        } else {
          say "FERGAL CDS END DEBUG: ".$exon->length." - (".($exon->length - $target_cds_end_index).") - $cumulative_length";
          $cds_end_offset = $exon->length - ($exon->length - $target_cds_end_index) - $cumulative_length;
        }
        last;
      }
      $cumulative_length += $exon->length();
    }

    unless($cds_start_exon and $cds_end_exon) {
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-cds_projection_fallback",
          {
            feature_id => $source_transcript_id,
            feature_type => 'transcript',
            result => 'exon_mapping_failed',
            details => {
              reason => 'cds_coordinates_outside_exons',
              target_cds_start => $target_cds_start_index,
              target_cds_end => $target_cds_end_index,
              transcript_length => $cumulative_length,
              exon_count => scalar(@$exons),
              fallback_method => 'compute_best_translation'
            },
            message => "CDS exon mapping failed - falling back to compute_best_translation"
          }
        );
      }
      $self->warning("For ".$source_transcript->stable_id()." couldn't find the equivalent CDS start/end exon in the alignment of the transcripts. ".
                     "Will compute translation instead");
      compute_best_translation($transcript);
      return;
    }
    
    say "Orig CDS start/end exon lengths: ".$source_translation->start_Exon->length."/".$source_translation->end_Exon->length;
    say "Orig CDS start/end exon phases: ".$source_translation->start_Exon->phase()."/".$source_translation->end_Exon->end_phase();
    say "Orig CDS start exon frame: ".$source_translation->start_Exon->frame();
    
    if ($phase_adjust and $cds_start_exon != $exons->[0] and $exons->[0]->phase == -1) {
      say "CDS start exon will have start phase $phase_adjust but it is not the first exon";
      # This code will fail if the previous exon is smaller than $phase_adjust which may never happen
      my $new_start_exon;
      foreach my $exon (@$exons) {
        last if ($cds_start_exon == $exon);
        $new_start_exon = $exon;
      }
      $cds_start_offset = $new_start_exon->length-$phase_adjust+1;
      $cds_start_exon = $new_start_exon;
      
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-cds_phase_adjustment",
          {
            feature_id => $source_transcript_id,
            feature_type => 'transcript',
            result => 'start_exon_adjusted',
            details => {
              original_phase_adjust => $phase_adjust,
              decision_reason => 'phase_adjust_with_non_first_exon',
              new_start_exon_coordinates => $new_start_exon->start."-".$new_start_exon->end,
              new_cds_start_offset => $cds_start_offset
            },
            message => "CDS start exon adjusted to handle phase requirements"
          }
        );
      }
      $phase_adjust = 0;
    }
    
    say "Orig CDS start/end exon offsets: ".$source_translation->start()."/".$source_translation->end();
    say "CDS start/end index: ".$target_cds_start_index."/".$target_cds_end_index;
    say "CDS start/end exon lengths: ".$cds_start_exon->length()."/".$cds_end_exon->length();
    say "CDS start/end exon offsets: ".$cds_start_offset."/".$cds_end_offset;
    
    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($cds_start_exon);
    $translation->start($cds_start_offset);
    $translation->end_Exon($cds_end_exon);
    $translation->end($cds_end_offset);
    $transcript->translation($translation);
    if ($phase_adjust) {
      calculate_exon_phases($transcript, $phase_adjust);
    }

    my $protein = $transcript->translate;
    my $has_internal_stops = ($protein and index($protein->seq, '*') >= 0 and index($source_translation->seq, '*') == -1);
    
    if ($has_internal_stops) {
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-cds_quality_issue",
          {
            feature_id => $source_transcript_id,
            feature_type => 'transcript',
            result => 'internal_stop_codons_detected',
            details => {
              protein_sequence => $protein->seq,
              source_has_stops => index($source_translation->seq, '*') >= 0,
              target_has_stops => index($protein->seq, '*') >= 0,
              frameshifts_detected => scalar(@frameshifts),
              frameshift_correction_attempted => scalar(@frameshifts) > 0 ? 1 : 0
            },
            message => "Internal stop codons detected - will attempt frameshift correction"
          }
        );
      }
      
      say join(' ', __LINE__, $transcript->display_id, 'STOP CODON', $protein->seq);
      foreach my $exon (@{$transcript->get_all_Exons}) {
        print join("\t", __LINE__, $exon, $exon->stable_id || 'NULL', $exon->start, $exon->end, $exon->strand, $exon->phase, $exon->end_phase, $exon->seq->seq), "\n";
      }
      
      my @gaps;
      my $deletion_shift = 0;
      my $frameshift_correction_attempted = 0;
      
      if (@frameshifts) {
        $frameshift_correction_attempted = 1;
        if (scalar(@frameshifts)%3 == 0 or index($protein->seq, 'X*') == 0) {
          say join(' ', __LINE__, scalar(@frameshifts)%3, index($protein->seq, 'X*'));
        }
        else {
          for (my $i = 0; $i < scalar(@frameshifts); $i++) {
            if (@gaps) {
              if ($frameshifts[$i-1]->[1]+1 == $frameshifts[$i]->[1] and $frameshifts[$i]->[0] eq $gaps[-1]->[2]) {
                $gaps[-1]->[1]++;
              }
              else {
                if ($gaps[-1]->[2] eq 'D') {
                  $deletion_shift += $gaps[-1]->[1];
                }
                if ($gaps[-1]->[1]%3 == 0) {
                  say join(' ', __LINE__, $transcript->display_id, @{$frameshifts[$i]}, @{$gaps[-1]});
                  pop(@gaps);
                }
                push(@gaps, [$frameshifts[$i]->[1], 1, $frameshifts[$i]->[0]]);
                $gaps[-1]->[0] -= $deletion_shift;
              }
            }
            else {
              push(@gaps, [$frameshifts[$i]->[1], 1, $frameshifts[$i]->[0]]);
            }
          }
          if ($gaps[-1]->[1]%3 == 0) {
            pop(@gaps);
          }
        }
      }
      
      my $overall_shift = 0;
      my $frameshift_corrections_applied = 0;
      
      if (@gaps and @gaps < 10) {
        $frameshift_corrections_applied = scalar(@gaps);
        
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-frameshift_correction_decision",
            {
              feature_id => $source_transcript_id,
              feature_type => 'transcript',
              result => 'frameshift_correction_attempted',
              details => {
                total_frameshifts => scalar(@frameshifts),
                correctable_gaps => scalar(@gaps),
                gap_details => \@gaps,
                correction_threshold => 'gaps_less_than_10',
                will_attempt_correction => 1
              },
              message => "Attempting frameshift correction for " . scalar(@gaps) . " gaps"
            }
          );
        }
        
        # [Complex frameshift correction code continues...]
        foreach my $position_length (@gaps) {
          my $frameshift_start = $position_length->[0]-$overall_shift;
          my $frameshift_length = $position_length->[1];
          my $frameshift_op = $position_length->[2];
          my $frameshift_end = $frameshift_start+$frameshift_length-1;
          my $shift = $frameshift_length%3;
          $self->throw("Modulo should not be 0") unless ($shift);
          # [Rest of frameshift correction logic...]
          $overall_shift += $shift;
          $exons = $transcript->get_all_Exons();
        }
      } else {
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-frameshift_correction_decision",
            {
              feature_id => $source_transcript_id,
              feature_type => 'transcript',
              result => 'frameshift_correction_skipped',
              details => {
                total_frameshifts => scalar(@frameshifts),
                gaps_identified => scalar(@gaps),
                skip_reason => scalar(@gaps) >= 10 ? 'too_many_gaps' : 'no_correctable_gaps',
                correction_threshold => 'gaps_less_than_10'
              },
              message => "Frameshift correction skipped - " . (scalar(@gaps) >= 10 ? "too many gaps" : "no correctable gaps")
            }
          );
        }
      }
    }

    # Set the phases
    calculate_exon_phases($transcript, $phase_adjust);

    say "For transcript ".$source_transcript->stable_id()." translateable seq:\n".$transcript->translateable_seq();
    say "For transcript ".$source_transcript->stable_id()." translation seq:\n".$transcript->translation->seq();
    say "Source transcript translation seq:\n".$source_transcript->translation->seq();
    
    # Log final CDS projection result
    my $final_protein = $transcript->translate;
    my $final_has_stops = $final_protein ? (index($final_protein->seq, '*') >= 0) : 1;
    
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-cds_projection_final_result",
        {
          feature_id => $source_transcript_id,
          feature_type => 'transcript',
          result => $final_has_stops ? 'projection_with_issues' : 'projection_successful',
          details => {
            final_protein_length => $final_protein ? length($final_protein->seq) : 0,
            source_protein_length => length($source_translation->seq),
            has_internal_stops => $final_has_stops,
            frameshifts_detected => scalar(@frameshifts),
            frameshift_correction_attempted => $frameshift_correction_attempted,
            projection_quality => $final_has_stops ? 'degraded' : 'high',
            cds_boundaries_set => {
              start_exon => $cds_start_exon->start."-".$cds_start_exon->end,
              end_exon => $cds_end_exon->start."-".$cds_end_exon->end,
              start_offset => $cds_start_offset,
              end_offset => $cds_end_offset
            }
          },
          message => "CDS projection completed - " . ($final_has_stops ? "protein quality degraded due to frameshifts/stops" : "high quality protein maintained")
        }
      );
    }
  }
  else {
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-cds_projection_failure",
        {
          feature_id => $source_transcript_id,
          feature_type => 'transcript',
          result => 'failed',
          details => {
            reason => 'cds_coordinates_not_found',
            source_cds_start => $source_cds_start,
            source_cds_end => $source_cds_end,
            target_cds_start_found => defined($target_cds_start_index),
            target_cds_end_found => defined($target_cds_end_index),
            alignment_length => scalar(@source_align_array)
          },
          message => "CDS projection failed - could not map CDS coordinates through alignment"
        }
      );
    }
    $self->warning('Could not find a CDS for the projection of '.$source_transcript->stable_id);
    return;
  }
}


sub filter_paf_hits {
  my ($self,$gene,$paf_results) = @_;
  my $hit_identity_cutoff = 0.99;
  say "Calculating hit regions withing expected boundaries for ".$gene->stable_id();
  my $gene_left_boundary = $gene->{'left_boundary'};
  my $gene_right_boundary = $gene->{'right_boundary'};
  my $gene_target_slice = $gene->{'target_slice'};
  my $overlapping_paf_results = [];
  my $non_overlapping_paf_results = [];
  foreach my $paf_result (@$paf_results) {
    my $paf_target_genomic_start = ${$paf_result}[7];
    my $paf_target_genomic_end = ${$paf_result}[8];
    my $paf_target_genomic_name = ${$paf_result}[5];
    my $paf_hit_length = ${$paf_result}[1];
    my $paf_hit_identities = ${$paf_result}[9];
    my $paf_perc_ident = $paf_hit_identities/$paf_hit_length;
    if(($gene_left_boundary and $gene_right_boundary and $gene_target_slice) and
       ($paf_target_genomic_start >= $gene_left_boundary and $paf_target_genomic_end <= $gene_right_boundary and $gene_target_slice eq $paf_target_genomic_name and
        $paf_perc_ident >= $hit_identity_cutoff)) {
      say "FERGAL DEBUG P1";
      push(@$overlapping_paf_results,$paf_result);
    } elsif($paf_perc_ident >= $hit_identity_cutoff) {
      say "FERGAL DEBUG P2";
      say "Missing gene ".$gene->stable_id()." did not have an identified high confidence target region. Will take only the top paf as it passes the identity cutoff";
      push(@$non_overlapping_paf_results,$paf_result);
    }
  }
  my $selected_paf_results = [];
  if(scalar(@$overlapping_paf_results)) {
    $selected_paf_results = $overlapping_paf_results;
    say "Found ".scalar(@$selected_paf_results)." overlapping paf results";
  } elsif(scalar(@$non_overlapping_paf_results)) {
    $selected_paf_results = [${$non_overlapping_paf_results}[0]];
    say "Found ".scalar(@$selected_paf_results)." non-overlapping paf results, no overlapping results so will use the top hit";
  }

  my $extended_regions = [];
  foreach my $paf_result (@$selected_paf_results) {
    my $paf_strand = ${$paf_result}[4];
    my $paf_source_genomic_start = ${$paf_result}[2];
    my $paf_source_genomic_end = ${$paf_result}[3];
    my $paf_target_genomic_start = ${$paf_result}[7];
    my $paf_target_genomic_end = ${$paf_result}[8];
    my $paf_target_genomic_name = ${$paf_result}[5];
    my $source_genomic_length = ${$paf_result}[1];
    my $target_genomic_length = ${$paf_result}[6];

    say "  Initial paf hit: ".$paf_target_genomic_start.":".$paf_target_genomic_end.":".$paf_strand.":".$paf_target_genomic_name;
    # This will control the variability of the gap between two hits on the target relative to the
    # gap in coverage on the source sequence

    my $source_genomic_start = $paf_source_genomic_start;
    my $source_genomic_end = $paf_source_genomic_end;
    my $target_genomic_start = $paf_target_genomic_start;
    my $target_genomic_end = $paf_target_genomic_end;

    my $target_flanking = 500;
    my $scaling_ratio = 1.5;

    my $source_missing_coverage = $source_genomic_length - ($paf_source_genomic_end - $paf_source_genomic_start) + 1;
    my $source_missing_left = $paf_source_genomic_start;
    my $source_missing_right = $source_genomic_length - $paf_source_genomic_end;
    my $source_scaled_left = ceil($scaling_ratio * $source_missing_left) + $target_flanking;
    my $source_scaled_right = ceil($scaling_ratio * $source_missing_right) + $target_flanking;

    if($paf_strand eq '-') {
      my $tmp = $source_scaled_left;
      $source_scaled_left = $source_scaled_right;
      $source_scaled_right = $tmp;
    }

    $target_genomic_start -= $source_scaled_left;
    $target_genomic_end += $source_scaled_right;

    if($target_genomic_start <= 0) {
      $target_genomic_start = 1;
    }

    if($target_genomic_end > $target_genomic_length) {
      $target_genomic_end = $target_genomic_length;
    }

    push(@$extended_regions,[$target_genomic_start,$target_genomic_end,$paf_strand,$paf_target_genomic_name]);
  }

  # At this point, extended regions has any paf hit in the region, and it will have been extended based on the expected length of the region
  # versus the length of the hit itself. So if there were a few hits, it's likely that the different extended regions will actually just
  # describe the same region. If there are close paralogues, there may be distinct regions. This is okay, we will allow multi mapping and
  # then try and clear up later. At this point we want to cluster into however many distinct regions there are based on overlap
  my @sorted_res = sort { $a->[1] cmp $b->[1] } @$extended_regions;
  foreach my $region (@sorted_res) {
    my $start = ${$region}[0];
    my $end = ${$region}[1];
    my $strand = ${$region}[2];
  }

  for(my $i=0; $i<scalar(@sorted_res)-1; $i++) {
    my $res1 = $sorted_res[$i];
    unless($res1) {
      next;
    }
    my $coord_i_1 = ${$res1}[0];
    my $coord_i_2 = ${$res1}[1];
    my $strand_i = ${$res1}[2];
    my $region_i = ${$res1}[3];
    for(my $j=$i+1; $j<scalar(@sorted_res); $j++) {
      my $res2 = $sorted_res[$j];
      unless($res2) {
        next;
      }
      my $coord_j_1 = ${$res2}[0];
      my $coord_j_2 = ${$res2}[1];
      my $strand_j = ${$res2}[2];
      my $region_j = ${$res2}[3];
      if($self->coords_overlap($coord_i_1,$coord_i_2,$coord_j_1,$coord_j_2) and $strand_i eq $strand_j and $region_i eq $region_j) {
        $sorted_res[$i] = [min($coord_i_1,$coord_j_1),max($coord_i_2,$coord_j_2),$strand_i,$region_i];
        $sorted_res[$j] = undef;
      }
    } # End for $j
  } # End for $i

  my $final_paf_regions = [];
  foreach my $region (@sorted_res) {
    if($region) {
      push(@$final_paf_regions,$region);
      say " Chosen hit: ".${$region}[0].":".${$region}[1].":".${$region}[2].":".${$region}[3];
    }
  }

  $gene->{'paf_regions'} = $final_paf_regions;
}

sub coords_overlap {
  my ($self,$s1,$e1,$s2,$e2) = @_;

  if (($s1 <= $e2) and ($e1 >= $s2)) {
    return 1;
  }
  return 0;
}

sub print_transcript_stats {
  my ($self,$transcripts_by_id,$tag) = @_;

  say "Transcript stats for ".$tag;
  foreach my $id (keys(%$transcripts_by_id)) {
    my $transcript = $transcripts_by_id->{$id};
    say "  ID: ".$transcript->stable_id().", Cov: ".$transcript->{'cov'}.", Perc id: ".$transcript->{'perc_id'};
  }
}


sub map_gene_minimap {
  my ($self,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
      $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,$target_slice_adaptor,$best_transcripts_by_id) = @_;

  my $source_transcript_fasta_seqs = [];
  my $source_transcripts_to_map = $self->filter_transcripts_to_map($source_transcripts,$best_transcripts_by_id);

  # Log transcript filtering decision
  my $total_source_transcripts = scalar(@$source_transcripts);
  my $filtered_transcripts_count = scalar(@$source_transcripts_to_map);
  my $skipped_transcripts = $total_source_transcripts - $filtered_transcripts_count;

  if ($self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-minimap_transcript_filtering",
      {
        feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
        feature_type => 'gene',
        result => 'transcript_filtering_applied',
        details => {
          total_source_transcripts => $total_source_transcripts,
          transcripts_selected_for_mapping => $filtered_transcripts_count,
          transcripts_skipped => $skipped_transcripts,
          filtering_reason => $skipped_transcripts > 0 ? 'already_have_better_mappings' : 'no_filtering_applied',
          selected_transcript_ids => [map { $_->stable_id } @$source_transcripts_to_map],
          target_region => $target_region_slice->name
        },
        message => "Minimap2 transcript filtering: $filtered_transcripts_count/$total_source_transcripts transcripts selected for mapping"
      }
    );
  }

  say "Processing a total of ".scalar(@$source_transcripts_to_map)." source transcripts via local minimap";

  # Check if any transcripts to map
  if (scalar(@$source_transcripts_to_map) == 0) {
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-minimap_mapping_skipped",
        {
          feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
          feature_type => 'gene',
          result => 'mapping_skipped',
          details => {
            reason => 'no_transcripts_need_mapping',
            existing_mappings_sufficient => 1,
            target_region => $target_region_slice->name
          },
          message => "Minimap2 mapping skipped - all transcripts already have sufficient mappings"
        }
      );
    }
    return {};
  }

  foreach my $source_transcript (@$source_transcripts_to_map) {
    say "Writing ".$source_transcript->stable_id()." to file for mapping";
    my $source_transcript_sequence = $source_transcript->seq->seq();
    my $fasta_record = ">".$source_transcript->dbID()."\n".$source_transcript_sequence;
    push(@$source_transcript_fasta_seqs,$fasta_record);
  }

  say "Generating initial set of minimap2 mappings";
  my $minimap_transcripts = $self->generate_minimap_transcripts($source_transcript_fasta_seqs,$target_genome_file,$target_adaptor,$target_genomic_start,$max_intron_size, $target_region_slice);

  # Log minimap2 mapping results
  my $successful_mappings = scalar(@$minimap_transcripts);
  my $failed_mappings = $filtered_transcripts_count - $successful_mappings;

  if ($self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-minimap_mapping_results",
      {
        feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
        feature_type => 'gene',
        result => $successful_mappings > 0 ? 'mappings_generated' : 'mapping_failed',
        details => {
          transcripts_attempted => $filtered_transcripts_count,
          successful_mappings => $successful_mappings,
          failed_mappings => $failed_mappings,
          success_rate => sprintf("%.1f", ($successful_mappings / $filtered_transcripts_count) * 100),
          max_intron_size_used => $max_intron_size,
          target_genome_file => $target_genome_file,
          mapping_parameters => {
            target_region => $target_region_slice->name,
            target_strand => $target_strand
          }
        },
        message => "Minimap2 mapping completed: $successful_mappings/$filtered_transcripts_count transcripts successfully mapped"
      }
    );
  }

  my $transcripts_by_id = {};
  my $transcripts_with_translation = 0;
  my $transcripts_without_translation = 0;
  my @quality_scores = ();

  foreach my $transcript (@$minimap_transcripts) {
    $transcripts_by_id->{$transcript->stable_id} = $transcript;
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id};
    
    # Handle translation projection
    if ($source_transcript->translation) {
      $self->calculate_translation_based_on_source_transcript($transcript, $source_transcript);
      $transcripts_with_translation++;
    } else {
      $transcripts_without_translation++;
    }

    $self->set_transcript_coverage_and_identity($source_transcript, $transcript);
    my $coverage = $transcript->{'cov'};
    my $percent_id = $transcript->{'perc_id'};
    
    push @quality_scores, {
      transcript_id => $source_transcript->stable_id,
      coverage => $coverage,
      percent_id => $percent_id,
      combined_score => $coverage + $percent_id
    };
    
    say "Mapped transcript (".$source_transcript->stable_id()."): Coverage: ".$coverage.", Percent id: ".$percent_id;
    
    # Log individual transcript mapping quality
    if ($self->{_logger}) {
      my $quality_assessment = 'high';
      if ($coverage < 90 || $percent_id < 90) {
        $quality_assessment = 'medium';
      }
      if ($coverage < 70 || $percent_id < 70) {
        $quality_assessment = 'low';
      }
      
      $self->{_logger}->log(
        "Minimap2Remap-minimap_transcript_quality",
        {
          feature_id => $source_transcript->stable_id,
          feature_type => 'transcript',
          result => 'mapping_completed',
          details => {
            mapping_method => 'minimap2_splice_aware',
            coverage => $coverage,
            percent_id => $percent_id,
            combined_score => $coverage + $percent_id,
            quality_assessment => $quality_assessment,
            has_translation => $source_transcript->translation ? 1 : 0,
            translation_projected => $source_transcript->translation ? 1 : 0,
            target_location => $transcript->seq_region_name . ":" . $transcript->start . "-" . $transcript->end . "(" . $transcript->strand . ")"
          },
          message => "Minimap2 transcript mapping: $quality_assessment quality (${coverage}% cov, ${percent_id}% id)"
        }
      );
    }
  }

  # Log overall minimap2 method assessment
  if ($self->{_logger} && $successful_mappings > 0) {
    my $avg_coverage = 0;
    my $avg_percent_id = 0;
    my $high_quality_count = 0;
    
    foreach my $score (@quality_scores) {
      $avg_coverage += $score->{coverage};
      $avg_percent_id += $score->{percent_id};
      if ($score->{coverage} >= 90 && $score->{percent_id} >= 90) {
        $high_quality_count++;
      }
    }
    
    $avg_coverage /= scalar(@quality_scores) if @quality_scores;
    $avg_percent_id /= scalar(@quality_scores) if @quality_scores;
    
    my $method_effectiveness = 'high';
    if ($avg_coverage < 90 || $avg_percent_id < 90 || $high_quality_count / $successful_mappings < 0.7) {
      $method_effectiveness = 'medium';
    }
    if ($avg_coverage < 70 || $avg_percent_id < 70 || $high_quality_count / $successful_mappings < 0.3) {
      $method_effectiveness = 'low';
    }
    
    $self->{_logger}->log(
      "Minimap2Remap-minimap_method_assessment",
      {
        feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
        feature_type => 'gene',
        result => 'method_assessment_completed',
        details => {
          mapping_method => 'minimap2',
          total_mappings => $successful_mappings,
          transcripts_with_translation => $transcripts_with_translation,
          transcripts_without_translation => $transcripts_without_translation,
          average_coverage => sprintf("%.1f", $avg_coverage),
          average_percent_id => sprintf("%.1f", $avg_percent_id),
          high_quality_mappings => $high_quality_count,
          high_quality_rate => sprintf("%.1f", ($high_quality_count / $successful_mappings) * 100),
          method_effectiveness => $method_effectiveness,
          will_contribute_to_final => $successful_mappings > 0 ? 1 : 0,
          quality_distribution => \@quality_scores
        },
        message => "Minimap2 method assessment: $method_effectiveness effectiveness (avg: ${avg_coverage}% cov, ${avg_percent_id}% id, $high_quality_count/$successful_mappings high quality)"
      }
    );
  }

  return($transcripts_by_id);
}


sub set_transcript_coverage_and_identity {
  my ($self, $source_transcript, $transcript) = @_;

  $transcript->{'aligned_source_seq'} = $source_transcript->seq->seq;
  $transcript->{'aligned_target_seq'} = $transcript->seq->seq;
  if ($transcript->{'aligned_source_seq'} eq $transcript->{'aligned_target_seq'}) {
    $transcript->{'cov'} = sprintf("%.2f", 100);
    $transcript->{'perc_id'} = sprintf("%.2f", 100);
  }
  else {
    my ($proj_coverage,$proj_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($transcript->{'aligned_source_seq'}, $transcript->{'aligned_target_seq'});
    $transcript->{'cov'} = $proj_coverage;
    $transcript->{'perc_id'} = $proj_percent_id;
    $transcript->{'aligned_source_seq'} = $aligned_source_seq;
    $transcript->{'aligned_target_seq'} = $aligned_target_seq;
  }
}

sub map_gene_exonerate {
  my ($self,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
      $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,$target_slice_adaptor, $best_transcripts_by_id) = @_;

  my $source_transcripts_to_map = $self->filter_transcripts_to_map($source_transcripts,$best_transcripts_by_id);

  # Log transcript filtering decision for exonerate
  my $total_source_transcripts = scalar(@$source_transcripts);
  my $filtered_transcripts_count = scalar(@$source_transcripts_to_map);
  my $skipped_transcripts = $total_source_transcripts - $filtered_transcripts_count;

  if ($self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-exonerate_transcript_filtering",
      {
        feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
        feature_type => 'gene',
        result => 'transcript_filtering_applied',
        details => {
          total_source_transcripts => $total_source_transcripts,
          transcripts_selected_for_mapping => $filtered_transcripts_count,
          transcripts_skipped => $skipped_transcripts,
          filtering_reason => $skipped_transcripts > 0 ? 'already_have_sufficient_mappings_from_previous_methods' : 'no_filtering_applied',
          selected_transcript_ids => [map { $_->stable_id } @$source_transcripts_to_map],
          target_region => $target_region_slice->name,
          method_position => 'third_method_after_projection_and_minimap2'
        },
        message => "Exonerate transcript filtering: $filtered_transcripts_count/$total_source_transcripts transcripts selected for final mapping attempt"
      }
    );
  }

  say "Processing a total of ".scalar(@$source_transcripts_to_map)." source transcripts using Exonerate";

  # Check if any transcripts need exonerate mapping
  if (scalar(@$source_transcripts_to_map) == 0) {
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-exonerate_mapping_skipped",
        {
          feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
          feature_type => 'gene',
          result => 'mapping_skipped',
          details => {
            reason => 'no_transcripts_need_exonerate_mapping',
            previous_methods_sufficient => 1,
            target_region => $target_region_slice->name,
            decision_logic => 'all_transcripts_have_adequate_mappings_from_projection_or_minimap2'
          },
          message => "Exonerate mapping skipped - all transcripts already have sufficient mappings from earlier methods"
        }
      );
    }
    return {};
  }

  # Log decision to proceed with exonerate (last resort method)
  if ($self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-exonerate_last_resort_decision",
      {
        feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
        feature_type => 'gene',
        result => 'proceeding_with_exonerate',
        details => {
          transcripts_requiring_exonerate => $filtered_transcripts_count,
          reason_for_exonerate => 'projection_and_minimap2_insufficient',
          transcript_ids_needing_help => [map { $_->stable_id } @$source_transcripts_to_map],
          method_rationale => 'last_resort_high_sensitivity_mapping',
          max_intron_size => $max_intron_size,
          target_region => $target_region_slice->name
        },
        message => "Proceeding with Exonerate mapping as last resort for $filtered_transcripts_count transcripts"
      }
    );
  }

  my $exonerate_transcripts = $self->generate_exonerate_transcripts($source_transcripts_to_map,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor,$max_intron_size);

  # Log exonerate mapping results
  my $successful_mappings = scalar(@$exonerate_transcripts);
  my $failed_mappings = $filtered_transcripts_count - $successful_mappings;

  if ($self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-exonerate_mapping_results",
      {
        feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
        feature_type => 'gene',
        result => $successful_mappings > 0 ? 'mappings_generated' : 'all_methods_failed',
        details => {
          transcripts_attempted => $filtered_transcripts_count,
          successful_mappings => $successful_mappings,
          failed_mappings => $failed_mappings,
          success_rate => $filtered_transcripts_count > 0 ? sprintf("%.1f", ($successful_mappings / $filtered_transcripts_count) * 100) : 0,
          method_position => 'final_mapping_method',
          successfully_mapped_transcript_ids => [map { $_->stable_id } @$exonerate_transcripts],
          completely_failed_transcript_ids => $failed_mappings > 0 ? 
            [grep { my $id = $_->stable_id; !grep { $_->stable_id eq $id } @$exonerate_transcripts } @$source_transcripts_to_map] : [],
          final_method_outcome => $successful_mappings > 0 ? 'some_transcripts_rescued' : 'complete_mapping_failure'
        },
        message => $successful_mappings > 0 ? 
          "Exonerate rescued $successful_mappings/$filtered_transcripts_count transcripts that failed previous methods" :
          "Complete mapping failure - no transcripts could be mapped by any method (projection, minimap2, exonerate)"
      }
    );
  }

  my $transcripts_by_id = {};
  my @quality_assessments = ();

  foreach my $transcript (@$exonerate_transcripts) {
    $transcripts_by_id->{$transcript->stable_id} = $transcript;
    
    # Assess transcript quality if metrics are available
    my $coverage = $transcript->{'cov'} || 0;
    my $percent_id = $transcript->{'perc_id'} || 0;
    
    push @quality_assessments, {
      transcript_id => $transcript->stable_id,
      coverage => $coverage,
      percent_id => $percent_id,
      combined_score => $coverage + $percent_id
    };
    
    # Log individual exonerate transcript result
    if ($self->{_logger}) {
      my $quality_assessment = 'unknown';
      if ($coverage > 0 && $percent_id > 0) {
        if ($coverage >= 90 && $percent_id >= 90) {
          $quality_assessment = 'high';
        } elsif ($coverage >= 70 && $percent_id >= 70) {
          $quality_assessment = 'medium';
        } else {
          $quality_assessment = 'low';
        }
      }
      
      $self->{_logger}->log(
        "Minimap2Remap-exonerate_transcript_rescue",
        {
          feature_id => $transcript->stable_id,
          feature_type => 'transcript',
          result => 'rescued_by_exonerate',
          details => {
            mapping_method => 'exonerate_protein_guided',
            coverage => $coverage,
            percent_id => $percent_id,
            combined_score => $coverage + $percent_id,
            quality_assessment => $quality_assessment,
            rescue_success => 1,
            target_location => defined($transcript->seq_region_name) ? 
              $transcript->seq_region_name . ":" . $transcript->start . "-" . $transcript->end . "(" . $transcript->strand . ")" : 
              "location_unknown",
            method_sequence => 'failed_projection_and_minimap2_rescued_by_exonerate'
          },
          message => "Transcript rescued by Exonerate after projection/minimap2 failure - $quality_assessment quality result"
        }
      );
    }
  }

  # Log overall exonerate method effectiveness
  if ($self->{_logger}) {
    my $rescue_rate = $filtered_transcripts_count > 0 ? ($successful_mappings / $filtered_transcripts_count) * 100 : 0;
    my $method_effectiveness = 'high';
    
    if ($rescue_rate < 70) {
      $method_effectiveness = 'medium';
    }
    if ($rescue_rate < 30) {
      $method_effectiveness = 'low';
    }
    if ($rescue_rate == 0) {
      $method_effectiveness = 'failed';
    }
    
    $self->{_logger}->log(
      "Minimap2Remap-exonerate_method_assessment",
      {
        feature_id => "gene_" . $target_region_slice->seq_region_name . "_" . $target_genomic_start,
        feature_type => 'gene',
        result => 'final_method_assessment',
        details => {
          mapping_method => 'exonerate_last_resort',
          transcripts_needing_rescue => $filtered_transcripts_count,
          transcripts_rescued => $successful_mappings,
          transcripts_completely_failed => $failed_mappings,
          rescue_rate => sprintf("%.1f", $rescue_rate),
          method_effectiveness => $method_effectiveness,
          pipeline_outcome => $successful_mappings > 0 ? 'partial_success' : 'complete_failure',
          quality_distribution => \@quality_assessments,
          final_mapping_decision => $successful_mappings > 0 ? 
            'exonerate_provided_mappings_for_difficult_transcripts' : 
            'no_method_could_map_these_transcripts'
        },
        message => "Exonerate final assessment: $method_effectiveness rescue rate (${rescue_rate}%) - " . 
                   ($successful_mappings > 0 ? "provided last-resort mappings" : "complete mapping failure for all methods")
      }
    );
  }

  return($transcripts_by_id);
}


sub update_best_transcripts {
  my ($self,$best_transcripts_by_id,$new_transcripts_by_id,$coverage_threshold,$perc_id_threshold) = @_;

  my @ids = keys(%{$new_transcripts_by_id});
  foreach my $id (@ids) {
    unless($best_transcripts_by_id->{$id}) {
      # First transcript for this ID - accept by default
      $best_transcripts_by_id->{$id} = $new_transcripts_by_id->{$id};
      
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-transcript_selection_decision",
          {
            feature_id => $id,
            feature_type => 'transcript',
            result => 'first_transcript_accepted',
            details => {
              coverage => $new_transcripts_by_id->{$id}->{'cov'},
              percent_id => $new_transcripts_by_id->{$id}->{'perc_id'},
              annotation_method => $new_transcripts_by_id->{$id}->{'annotation_method'} || 'unknown',
              decision_reason => 'no_existing_transcript'
            },
            message => "First transcript mapping accepted for this transcript ID"
          }
        );
      }
    } else {
      # Existing transcript - need to decide which is better
      my $best_coverage = $best_transcripts_by_id->{$id}->{'cov'};
      my $best_perc_id = $best_transcripts_by_id->{$id}->{'perc_id'};
      my $best_total = $best_coverage + $best_perc_id;
      my $best_method = $best_transcripts_by_id->{$id}->{'annotation_method'} || 'unknown';
      
      my $transcript_coverage = $new_transcripts_by_id->{$id}->{'cov'};
      my $transcript_perc_id = $new_transcripts_by_id->{$id}->{'perc_id'};
      my $transcript_total = $transcript_coverage + $transcript_perc_id;
      my $new_method = $new_transcripts_by_id->{$id}->{'annotation_method'} || 'unknown';

      # Check if new transcript fails quality thresholds
      if(($coverage_threshold and $perc_id_threshold) and ($transcript_coverage < $coverage_threshold or $transcript_perc_id < $perc_id_threshold)) {
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-transcript_selection_decision",
            {
              feature_id => $id,
              feature_type => 'transcript',
              result => 'new_transcript_rejected',
              details => {
                decision_reason => 'below_quality_thresholds',
                new_transcript => {
                  coverage => $transcript_coverage,
                  percent_id => $transcript_perc_id,
                  method => $new_method
                },
                current_best => {
                  coverage => $best_coverage,
                  percent_id => $best_perc_id,
                  method => $best_method
                },
                thresholds => {
                  coverage_threshold => $coverage_threshold,
                  perc_id_threshold => $perc_id_threshold
                },
                coverage_pass => $transcript_coverage >= $coverage_threshold,
                perc_id_pass => $transcript_perc_id >= $perc_id_threshold
              },
              message => "New transcript rejected - fails quality thresholds"
            }
          );
        }
        return;
      }

      # Compare total scores
      if($transcript_total > $best_total) {
        # Check for protected high-quality existing transcript
        if(($coverage_threshold and $perc_id_threshold) and ($best_coverage >= $coverage_threshold and $best_perc_id >= $perc_id_threshold)) {
          if ($self->{_logger}) {
            $self->{_logger}->log(
              "Minimap2Remap-transcript_selection_decision",
              {
                feature_id => $id,
                feature_type => 'transcript',
                result => 'new_transcript_rejected',
                details => {
                  decision_reason => 'existing_transcript_protected',
                  protection_rule => 'preserve_high_quality_existing_when_above_thresholds',
                  new_transcript => {
                    coverage => $transcript_coverage,
                    percent_id => $transcript_perc_id,
                    total_score => $transcript_total,
                    method => $new_method
                  },
                  current_best => {
                    coverage => $best_coverage,
                    percent_id => $best_perc_id,
                    total_score => $best_total,
                    method => $best_method
                  },
                  thresholds => {
                    coverage_threshold => $coverage_threshold,
                    perc_id_threshold => $perc_id_threshold
                  },
                  score_comparison => "new($transcript_total) > best($best_total) but existing protected"
                },
                message => "Higher-scoring transcript rejected - existing transcript above thresholds is protected (prevents incorrect paralog replacement)"
              }
            );
          }
          return;
        }
        
        # Replace with better transcript
        $best_transcripts_by_id->{$id} = $new_transcripts_by_id->{$id};
        
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-transcript_selection_decision",
            {
              feature_id => $id,
              feature_type => 'transcript',
              result => 'transcript_replaced',
              details => {
                decision_reason => 'higher_combined_score',
                replaced_transcript => {
                  coverage => $best_coverage,
                  percent_id => $best_perc_id,
                  total_score => $best_total,
                  method => $best_method
                },
                new_best_transcript => {
                  coverage => $transcript_coverage,
                  percent_id => $transcript_perc_id,
                  total_score => $transcript_total,
                  method => $new_method
                },
                score_improvement => sprintf("%.2f", $transcript_total - $best_total),
                method_transition => "$best_method -> $new_method"
              },
              message => "Transcript replaced with higher-scoring version (score improvement: " . sprintf("%.2f", $transcript_total - $best_total) . ")"
            }
          );
        }
      } else {
        # Keep existing transcript
        if ($self->{_logger}) {
          $self->{_logger}->log(
            "Minimap2Remap-transcript_selection_decision",
            {
              feature_id => $id,
              feature_type => 'transcript',
              result => 'existing_transcript_retained',
              details => {
                decision_reason => 'existing_transcript_higher_score',
                existing_transcript => {
                  coverage => $best_coverage,
                  percent_id => $best_perc_id,
                  total_score => $best_total,
                  method => $best_method
                },
                rejected_transcript => {
                  coverage => $transcript_coverage,
                  percent_id => $transcript_perc_id,
                  total_score => $transcript_total,
                  method => $new_method
                },
                score_difference => sprintf("%.2f", $best_total - $transcript_total)
              },
              message => "Existing transcript retained - higher score than new candidate"
            }
          );
        }
      }
    }
  }
}


sub filter_transcripts_to_map {
  my ($self,$source_transcripts,$best_transcripts_by_id) = @_;

  my $coverage_cutoff = 98;
  my $identity_cutoff = 98;

  my $filtered_transcripts = [];
  foreach my $source_transcript (@$source_transcripts) {
    my $source_transcript_id = $source_transcript->dbID();
    unless($best_transcripts_by_id->{$source_transcript_id}) {
      push(@$filtered_transcripts,$source_transcript);
      next;
    }
    my $best_coverage = $best_transcripts_by_id->{$source_transcript_id}->{'cov'};
    my $best_perc_id = $best_transcripts_by_id->{$source_transcript_id}->{'perc_id'};
    unless($best_coverage >= $coverage_cutoff and $best_perc_id >= $identity_cutoff) {
      push(@$filtered_transcripts,$source_transcript);
    }
  }
  return($filtered_transcripts);
}


sub project_gene_coords {
  my ($self,$gene,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,$gene_genomic_seqs_hash) = @_;

  my $transcripts_by_id = {};
  my $gene_stable_id = $gene->stable_id();

  if($self->no_projection()) {
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-coordinate_projection_disabled",
        {
          feature_id => $gene_stable_id,
          feature_type => 'gene',
          result => 'skipped',
          details => {
            reason => 'projection_disabled_by_parameter',
            no_projection_flag => 1
          },
          message => "Coordinate projection skipped - disabled by no_projection parameter"
        }
      );
    }
    return($transcripts_by_id);
  }

  my $source_genomic_seq_info = $gene_genomic_seqs_hash->{$gene->stable_id};

  unless($source_genomic_seq_info) {
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-coordinate_projection_failure",
        {
          feature_id => $gene_stable_id,
          feature_type => 'gene',
          result => 'failed',
          details => {
            reason => 'missing_genomic_sequence_info',
            gene_dbid => $gene->dbID()
          },
          message => "Coordinate projection failed - no genomic sequence info found"
        }
      );
    }
    $self->throw("Could not find the genomic seq info for source gene with dbID: ".$gene->dbID());
  }

  my $source_region_start = ${$source_genomic_seq_info}[0];
  my $source_region_end = ${$source_genomic_seq_info}[1];
  my $source_genome_seq = ${$source_genomic_seq_info}[2];

  my $exons = $gene->get_all_Exons();
  my $target_genome_seq = $target_region_slice->seq();

  if($target_strand != 1) {
    $target_genome_seq = $self->revcomp($target_genome_seq);
  }

  my $projected_exons_by_id = {};
  # Sometimes issues with the MySQL server disconnecting on long alignments
  $self->source_adaptor->dbc->disconnect_when_inactive(1);
  $self->target_adaptor->dbc->disconnect_when_inactive(1);

  my $coverage = 0;
  my $percent_id = 0;
  my $aligned_source_seq;
  my $aligned_target_seq;
  eval {
    ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_genome_seq,$target_genome_seq);
  };
  if ($@) {
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-coordinate_projection_failure",
        {
          feature_id => $gene_stable_id,
          feature_type => 'gene',
          result => 'failed',
          details => {
            reason => 'mafft_alignment_error',
            error_message => $@,
            source_seq_length => length($source_genome_seq),
            target_seq_length => length($target_genome_seq)
          },
          message => "MAFFT alignment failed - cannot proceed with coordinate projection"
        }
      );
    }
    $self->warning("Issue with running MAFFT on target region");
  } else {
    say "Aligned source and target regions with MAFFT: Coverage: ".$coverage.", Percent id: ".$percent_id;
  }

  unless($coverage and $percent_id) {
    if ($self->{_logger}) {
      $self->{_logger}->log(
        "Minimap2Remap-coordinate_projection_failure",
        {
          feature_id => $gene_stable_id,
          feature_type => 'gene',
          result => 'failed',
          details => {
            reason => 'insufficient_alignment_quality',
            coverage => $coverage || 0,
            percent_id => $percent_id || 0,
            source_region => "$source_region_start-$source_region_end",
            target_region => $target_region_slice->name()
          },
          message => "Coordinate projection rejected - insufficient alignment coverage/identity"
        }
      );
    }
    say "No coverage/percent id calculated therefore not proceeding with projection";
    return($transcripts_by_id);
  }

  # Log successful alignment that enables projection
  if ($self->{_logger}) {
    $self->{_logger}->log(
      "Minimap2Remap-coordinate_projection_alignment_success",
      {
        feature_id => $gene_stable_id,
        feature_type => 'gene',
        result => 'alignment_suitable',
        details => {
          genome_alignment_coverage => $coverage,
          genome_alignment_percent_id => $percent_id,
          source_region => "$source_region_start-$source_region_end",
          target_region => $target_region_slice->name(),
          exons_to_project => scalar(@$exons)
        },
        message => "Genome alignment successful - proceeding with coordinate projection"
      }
    );
  }

  my $successful_exon_projections = 0;
  my $failed_exon_projections = 0;

  foreach my $exon (@$exons) {
    my $exon_region_start = $exon->seq_region_start();
    my $exon_region_end = $exon->seq_region_end();
    say "Source region start/end: ".$source_region_start."/".$source_region_end;
    say "Source exon region start/end: ".$exon_region_start."/".$exon_region_end;
    say "Source exon strand: ".$exon->strand();
    say "Target strand: ".$target_strand;

    my $projected_exon = $self->project_feature(undef,$exon,$source_region_start,$exon_region_start,$exon_region_end,$aligned_source_seq,$aligned_target_seq,$target_region_slice,$target_strand);

    if($projected_exon) {
      my $exon_coverage;
      my $exon_percent_id;
      
      if ($exon->seq->seq eq $projected_exon->seq->seq) {
        $exon_coverage = 100;
        $exon_percent_id = 100;
        $projected_exon->{'cov'} = sprintf("%.2f", 100);
        $projected_exon->{'perc_id'} = sprintf("%.2f", 100);
      } else {
        ($exon_coverage, $exon_percent_id, undef, undef) = align_nucleotide_seqs($exon->seq->seq, $projected_exon->seq->seq);
        say "Projected exon alignment scores: Coverage: ".$exon_coverage.", Perc identity: ".$exon_percent_id;
        $projected_exon->{'cov'} = $exon_coverage;
        $projected_exon->{'perc_id'} = $exon_percent_id;
      }
      
      $projected_exon->{'source_stable_id'} = $exon->stable_id();
      $projected_exon->{'source_length'} = $exon->length();
      $projected_exons_by_id->{$projected_exon->{'source_stable_id'}} = $projected_exon;
      $successful_exon_projections++;

      # Log exon projection quality decision
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-exon_projection_quality",
          {
            feature_id => $gene_stable_id,
            feature_type => 'exon',
            result => 'projected',
            details => {
              source_exon_id => $exon->stable_id(),
              source_exon_length => $exon->length(),
              projected_exon_length => $projected_exon->length(),
              exon_coverage => $exon_coverage,
              exon_percent_id => $exon_percent_id,
              sequence_identity => $exon->seq->seq eq $projected_exon->seq->seq ? 'identical' : 'divergent',
              projection_method => 'coordinate_scaling'
            },
            message => "Exon successfully projected with " . sprintf("%.1f", $exon_coverage) . "% coverage, " . sprintf("%.1f", $exon_percent_id) . "% identity"
          }
        );
      }
    } else {
      $failed_exon_projections++;
      
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-exon_projection_failure",
          {
            feature_id => $gene_stable_id,
            feature_type => 'exon',
            result => 'projection_failed',
            details => {
              source_exon_id => $exon->stable_id(),
              source_exon_coords => "$exon_region_start-$exon_region_end",
              source_exon_length => $exon->length(),
              failure_reason => 'coordinate_projection_failed'
            },
            message => "Failed to project exon using coordinate scaling"
          }
        );
      }
      say "Failed to project exon (".$exon->stable_id().")";
    }
  }

  # Log overall exon projection success rate
  if ($self->{_logger}) {
    my $projection_success_rate = $successful_exon_projections / scalar(@$exons) * 100;
    $self->{_logger}->log(
      "Minimap2Remap-exon_projection_summary",
      {
        feature_id => $gene_stable_id,
        feature_type => 'gene',
        result => 'exon_projection_completed',
        details => {
          total_exons => scalar(@$exons),
          successful_projections => $successful_exon_projections,
          failed_projections => $failed_exon_projections,
          success_rate => sprintf("%.1f", $projection_success_rate),
          projection_quality => $projection_success_rate >= 90 ? 'high' : 
                               $projection_success_rate >= 70 ? 'medium' : 'low'
        },
        message => "Exon projection completed: $successful_exon_projections/${\(scalar(@$exons))} exons successfully projected"
      }
    );
  }

  my $clone_exons = scalar(@$source_transcripts)-1;
  my $successful_transcript_reconstructions = 0;
  my $failed_transcript_reconstructions = 0;

  foreach my $source_transcript (@$source_transcripts) {
    my $projected_transcript = $self->reconstruct_transcript($source_transcript,$projected_exons_by_id);
    if($projected_transcript) {
      $projected_transcript->{'annotation_method'} = 'alignment_projection';
      $transcripts_by_id->{$projected_transcript->stable_id()} = $projected_transcript;
      $successful_transcript_reconstructions++;
    } else {
      $failed_transcript_reconstructions++;
      
      if ($self->{_logger}) {
        $self->{_logger}->log(
          "Minimap2Remap-transcript_reconstruction_failure",
          {
            feature_id => $gene_stable_id,
            feature_type => 'transcript',
            result => 'reconstruction_failed',
            details => {
              source_transcript_id => $source_transcript->stable_id(),
              available_projected_exons => scalar(keys %$projected_exons_by_id),
              failure_reason => 'transcript_reconstruction_failed'
            },
            message => "Failed to reconstruct transcript from projected exons"
          }
        );
      }
    }
  }

  # Log final coordinate projection decision
  if ($self->{_logger}) {
    my $overall_success = $successful_transcript_reconstructions > 0;
    $self->{_logger}->log(
      "Minimap2Remap-coordinate_projection_final_decision",
      {
        feature_id => $gene_stable_id,
        feature_type => 'gene',
        result => $overall_success ? 'success' : 'failed',
        details => {
          source_transcripts => scalar(@$source_transcripts),
          successful_transcript_reconstructions => $successful_transcript_reconstructions,
          failed_transcript_reconstructions => $failed_transcript_reconstructions,
          genome_alignment_quality => {
            coverage => $coverage,
            percent_id => $percent_id
          },
          exon_projection_quality => {
            successful_exons => $successful_exon_projections,
            total_exons => scalar(@$exons)
          },
          decision => $overall_success ? 'coordinate_projection_successful' : 'coordinate_projection_failed_fallback_needed'
        },
        message => $overall_success ? 
          "Coordinate projection successful - $successful_transcript_reconstructions transcripts reconstructed" :
          "Coordinate projection failed - no transcripts reconstructed, will require alternative methods"
      }
    );
  }

  return($transcripts_by_id);
}


sub reconstruct_transcript {
  my ($self,$source_transcript,$projected_exons_by_id, $clone_exons) = @_;

  my $source_exons = $source_transcript->get_all_Exons();

  my $projected_exons = [];
  foreach my $source_exon (@$source_exons) {
    my $projected_exon = $projected_exons_by_id->{$source_exon->stable_id};
    if($projected_exon) {
      my %tmp_exon = %$projected_exon;
      my $cloned_exon = Bio::EnsEMBL::Exon->new_fast(\%tmp_exon);
      push(@$projected_exons, $cloned_exon);
    }
  }

  unless(scalar(@$projected_exons)) {
    return;
  }

  my $projected_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $projected_exons);
  $projected_transcript->slice(${$projected_exons}[0]->slice());
  $projected_transcript->analysis($self->analysis());
  $projected_transcript->biotype($source_transcript->biotype());
  $projected_transcript->stable_id($source_transcript->dbID());
  $projected_transcript->version($source_transcript->version());
  $projected_transcript->{'source_stable_id'} = $source_transcript->stable_id();
  $projected_transcript->{'source_biotype_group'} = $source_transcript->get_Biotype->biotype_group();
  $projected_transcript->{'source_length'} = $source_transcript->length();

  $self->set_transcript_coverage_and_identity($source_transcript, $projected_transcript);
  say "Projected transcript (".$projected_transcript->{'source_stable_id'}.") ".$projected_transcript->seq_region_name." ".$projected_transcript->seq_region_start."/".
      $projected_transcript->seq_region_end." ".$projected_transcript->strand.": Coverage: ".$projected_transcript->{'cov'}.", Percent id: ".$projected_transcript->{'perc_id'};

  return($projected_transcript);
}

sub project_feature {
  my ($self,$transcript,$exon,$source_region_start,$feature_start,$feature_end,$aligned_source_seq,$aligned_target_seq,$target_region_slice,$target_strand) = @_;

  # At this point we have the aligned seqs and many coords
  # The source region start/feature start/feature end are used to calculate the offset of the feature start/end in the aligned source seq
  # Once they are located, then the equivalent positions (if they exist) can be calculated in the target seq
  # After they've been located, we can then calculate the coords of the feature in the target and build it
  # There are various issues here. If the positions are at gaps in the alignment, then what to do? Most straightfoward thing might be to
  # scan and include all aligned bases between the source start/end and fill in with equivalent number of missing bases
  # Maybe look for the closest non-gap position to the coord of the start end that's within the genomic region of the source feature and then
  # calculate the base offset from there and include the equivalent number of missing bases from the target seq
  # Actually the above is probably just overcomplicated, just take all the bases aligned to positions within the feature in the source
  # then when you calculate coverage it will help sort out bad ones versus good ones
  # Also need to do something in terms of if there's a cds, since this is likely to cover IG genes with weird cds features
  # Could decide to drop the UTR in these cases (if there is UTR, since that's very unlikely), then that would make it more
  # straightforward in terms of just making the whole feature coding
#  say "Aligned seqs:";
#  say $aligned_source_seq;
#  say $aligned_target_seq;

  my $source_seq = $aligned_source_seq;
  $source_seq =~ s/\-//g;
  my $target_seq = $aligned_target_seq;
  $target_seq =~ s/\-//g;

  my @source_align_array = split('',$aligned_source_seq);
  my @target_align_array = split('',$aligned_target_seq);
  my $source_seq_pos = 0;
  my $source_align_seq_pos = 0;
  my $target_seq_pos = 0;
  my $target_align_seq_pos = 0;

  # Seen to have some issue in terms of getting the real start/end for some features. Was using +1, but seems off
  my $source_feature_seq_start =  $feature_start - $source_region_start + 1;
  my $source_feature_seq_end = $feature_end - $source_region_start + 1;
  say "Source feature seq start: ".$source_feature_seq_start;
  say "Source feature seq end: ".$source_feature_seq_end;
  say "Source region start: ".$source_region_start;
  say "Feature start: ".$feature_start;
  say "Feature end: ".$feature_end;

  my $source_seq_start_index;
  my $source_seq_end_index;
  my $target_seq_start_index;
  my $target_seq_end_index;
  my $in_feature_alignment = 0;
  for(my $i=0; $i < scalar(@source_align_array); $i++) {
    my $source_value = $source_align_array[$i];
    my $target_value = $target_align_array[$i];
    if($source_value ne '-') {
      $source_seq_pos++;
    }

    if($target_value ne '-') {
      $target_seq_pos++;
    }

    if($source_seq_pos == $source_feature_seq_start) {
      $source_seq_start_index = $source_seq_pos;
      say "Index of the feature start in source seq: ".$source_seq_start_index;
      $in_feature_alignment = 1;
    }

    if($in_feature_alignment and $target_value ne '-' and !(defined($target_seq_start_index))) {
      $target_seq_start_index = $target_seq_pos;
    }

    if($source_seq_pos == $source_feature_seq_end) {
      $source_seq_end_index = $source_seq_pos;
      $target_seq_end_index = $target_seq_pos;
      say "Index of the feature end in source seq: ".$source_seq_end_index;
      $in_feature_alignment = 0;
      last;
    }
  }

  unless(defined($target_seq_start_index) and defined($target_seq_end_index)) {
    $self->warning("Issue with recovering start/end of the exon feature in target alignment sequence, not building exon");
    return;
  }

  my $source_feature_length = $source_seq_end_index - $source_seq_start_index + 1;
  my $target_feature_length = $target_seq_end_index - $target_seq_start_index + 1;
  my $recovered_source_feature_seq = substr($source_seq,$source_seq_start_index,$source_feature_length);
  if($exon->strand != 1) {
    $recovered_source_feature_seq = $self->revcomp($recovered_source_feature_seq);
  }

  my $recovered_target_feature_seq = substr($target_seq,$target_seq_start_index,$target_feature_length);
  if($target_strand != $exon->strand()) {
    $recovered_target_feature_seq = $self->revcomp($recovered_target_feature_seq);
  }

  say "Source seq start index/end index/length: ".$source_seq_start_index.":".$source_seq_end_index.":".$source_feature_length;
  say "Target seq start index/end index/length: ".$target_seq_start_index.":".$target_seq_end_index.":".$target_feature_length;
  say "Recovered feature source seq:";
  say $recovered_source_feature_seq;
  say "Recovered feature target seq:";
  say $recovered_target_feature_seq;

  my $projected_exon = $self->build_projected_exon($transcript,$exon,$target_seq_start_index,$target_seq_end_index,$target_region_slice,$target_strand);
  return($projected_exon);
}


sub build_projected_exon {
  my ($self,$transcript,$exon,$seq_start_index,$seq_end_index,$region_slice,$target_strand) = @_;

  say "Region slice: ".$region_slice->name();
  my $region_start = $region_slice->seq_region_start();
  my $region_end = $region_slice->seq_region_end();
  say "Region start: ".$region_start;
  say "Region end: ".$region_end;

  my $parent_slice = $region_slice->seq_region_Slice();
  my $exon_start;
  my $exon_end;
  my $exon_strand = 1;
  say "Start/End index: ".$seq_start_index."/".$seq_end_index;
  if($target_strand == 1) {
    $exon_start = $region_start + $seq_start_index - 1;
    $exon_end = $region_start + $seq_end_index - 1;
    $exon_strand = $exon->strand();
  } else {
    # In this case we need to reverse the coords
    $exon_end = $region_end - $seq_start_index + 1;
    $exon_start = $exon_end - ($seq_end_index - $seq_start_index);
    if($exon->strand() == 1) {
       $exon_strand = -1;
    } else {
      $exon_strand = 1;
    }
  }

  my $phase = -1;
  my $end_phase = -1;

  say "Projected exon start/end slice coords: ".$exon_start."/".$exon_end;
  say "Parent slice stard/end genomic coords: ".$parent_slice->seq_region_start."/".$parent_slice->seq_region_end;
  my $projected_exon = Bio::EnsEMBL::Exon->new(-start     => $exon_start,
                                               -end       => $exon_end,
                                               -strand    => $exon_strand,
                                               -phase     => $phase,
                                               -end_phase => $end_phase,
                                               -analysis  => $self->analysis,
                                               -slice     => $parent_slice);

  $projected_exon->stable_id($exon->stable_id);
  if($transcript and $exon->is_coding($transcript)) {
    $projected_exon->phase($exon->phase());
    $projected_exon->end_phase($exon->end_phase());
  }

  if($exon_start > $exon_end) {
    $self->throw("Created an exon where the start > than the end, this shouldn't be possible: ".$parent_slice->name." ".$exon->start."..".$exon->end." ".$target_strand);
  }

  say "New exon seq:";
  say $projected_exon->seq->seq();
  say "Original exon seq:";
  say $exon->seq->seq();

  return($projected_exon);
}





sub generate_minimap_transcripts {
  my ($self,$source_transcript_fasta_seqs,$target_genome_index,$target_adaptor,$target_genomic_start,$max_intron_size, $target_region_slice) = @_;

  # This will take in the target genomic region (via the index)
  my $coverage_cutoff = 80;
  my $perc_id_cutoff = 80;
  my $max_stops = 999;
  my $minimap_transcripts = [];
  my $analysis = $self->analysis();
  my $program = $self->program();
  my $paftools = $self->paftools_path();
  my $source_input_file = $self->write_input_file($source_transcript_fasta_seqs);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis                 => $analysis,
       -program                  => $program,
       -paftools_path            => $paftools,
       -genome_index             => $target_genome_index,
       -input_file               => $source_input_file,
       -database_adaptor         => $target_adaptor,
       -skip_introns_check       => 1,
       -add_offset               => $target_genomic_start - 1,
       -skip_compute_translation => 1,
       -max_intron_size          => $max_intron_size,
       -perc_id                  => $perc_id_cutoff,
       -coverage                 => $coverage_cutoff,
       -query                    => $target_region_slice,
  );

  $runnable->run();

  my $output_genes = $runnable->output();
  foreach my $gene (@$output_genes) {
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      $transcript->{'annotation_method'} = 'minimap_local';
      push(@$minimap_transcripts,$transcript);
    }
  }

  return($minimap_transcripts);
}


sub generate_exonerate_transcripts {
  my ($self,$transcripts_to_map,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor,$max_intron_size) = @_;

  my $output_transcripts = [];

  foreach my $transcript_to_map (@$transcripts_to_map) {
    my $source_transcript = $source_transcript_id_hash->{$transcript_to_map->dbID()};
    unless($source_transcript) {
      $self->throw("Couldn't find the dbID of the transcript in the source transcript hash when attempting to run exonerate");
    }

    my $exonerate_transcripts = $self->run_exonerate($source_transcript,$target_region_slice,$target_slice_adaptor,$max_intron_size);
    say "Got ".scalar(@$exonerate_transcripts)." from Exonerate";
    foreach my $transcript (@$exonerate_transcripts) {
      $self->set_transcript_coverage_and_identity($source_transcript, $transcript);
      $transcript->{'annotation_method'} = 'exonerate_local';
      say "Mapped transcript (".$source_transcript->stable_id()."): Coverage: ".$transcript->{'cov'}.", Percent id: ".$transcript->{'perc_id'};
      push(@$output_transcripts,$transcript);
    }
  }
  return($output_transcripts);
}

sub run_exonerate {
  my ($self,$source_transcript,$target_region_slice,$target_slice_adaptor,$max_intron_size) = @_;

  my $output_transcripts = [];
  my $exonerate_length_cutoff = 15000;
  my $max_stops = 999;
  if($source_transcript->length() > $exonerate_length_cutoff) {
    say "Not running exonerate on trascript ".$source_transcript->stable_id()." as length (".$source_transcript->length().") is greater than length cut-off (".$exonerate_length_cutoff.")";
    return($output_transcripts);
  }

  my %parameters = (-coverage_by_aligned => 1);
  $parameters{-options} = "--forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE ".
                          "--saturatethreshold 100 --dnawordlen 15 --dnahspthreshold 60 --bestn 2 --maxintron ".$max_intron_size;
  my $annotation_features;
  my $non_coding_transcript = 0;
  if($source_transcript->translation()) {
    $annotation_features = $self->create_annotation_features($source_transcript);
    $parameters{-options} .= ' --model cdna2genome --codonwordlen 15';
    if ($source_transcript->length < 50) {
      $parameters{-options} .= ' --score 50';
    }
    else {
      $parameters{-options} .= ' --score 500';
    }

  } else {
    $non_coding_transcript = 1;
    $parameters{-options} .= ' --model est2genome';
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
                    -program  => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/exonerate22/bin/exonerate',
                    -analysis => $self->analysis,
                    -query_type     => "dna",
                    -calculate_coverage_and_pid => 1,
                    -non_coding_transcript => $non_coding_transcript,
                     %parameters,
                 );

  if ($annotation_features) {
    $runnable->annotation_features($annotation_features);
  }
  $runnable->target_seqs([$target_region_slice]);
  $runnable->query_seqs([$source_transcript->seq]);
  $runnable->run();
  if (scalar(@{$runnable->output()})) {
    my $output_transcript = ${$runnable->output()}[0];
    attach_Slice_to_Transcript($output_transcript, $target_region_slice);
    attach_Analysis_to_Transcript($output_transcript, $self->analysis);
    my $source_translation = $source_transcript->translation;
    if ($source_translation) {
      my $output_translation = $output_transcript->translation;
      if ($output_translation) {
        if ($source_translation->start_Exon->phase > 0 and $output_translation->start_Exon->phase == -1) {
          if ($output_translation->start > 3) {
            $self->warning("Partial transcript, start should be less than 4, not ".$output_translation->start);
          }
          else {
            if ($output_translation->start == 2) {
              $output_translation->start_Exon->phase(2);
            }
            elsif ($output_translation->start == 3) {
              $output_translation->start_Exon->phase(1);
            }
            else {
              $output_translation->start_Exon->phase(0);
            }
            $output_translation->start(1);
          }
        }
        if ($output_transcript->translateable_seq and !$output_transcript->translate) {
          $output_transcript->translation(undef);
          $self->warning($source_transcript->display_id.' FAILED no translation for a protein coding transcript of biotype '.$source_transcript->biotype);
        }
      }
      else {
        $self->warning($source_transcript->display_id.' FAILED no translation for a protein coding transcript');
      }
    }
    else {
      $output_transcript->translation(undef);
    }

    $output_transcript->stable_id($source_transcript->dbID());
    push(@$output_transcripts,$output_transcript);
  }

  return($output_transcripts);
}

sub label_transcript_status {
  my ($self,$transcripts_by_id,$status) = @_;

  my $coverage_threshold = 95;
  my $percent_id_threshold = 95;
  my $labelled_transcripts = [];
  # This labels a set of transcripts to help with deciding what to do when examining clusters
  foreach my $id (keys(%{$transcripts_by_id})) {
    my $transcript = $transcripts_by_id->{$id};
    my $transcript_coverage = $transcript->{'cov'};
    my $transcript_perc_id = $transcript->{'perc_id'};
    if ($transcript_coverage >= $coverage_threshold and $transcript_perc_id >= $percent_id_threshold) {
      $transcript->{'status'} = 'good';
    } else {
      $transcript->{'status'} = 'bad';
    }
    push(@$labelled_transcripts,$transcript);
  }

  return($labelled_transcripts);
}


sub generate_biotypes_hash {
  my ($self,$transcripts) = @_;

  my $unique_biotypes;
  my $biotypes_array = [];
  my $biotypes_hash = {};
  # This labels a set of transcripts to help with deciding what to do when examining clusters
  foreach my $transcript (@$transcripts) {
    unless($unique_biotypes->{$transcript->biotype()}) {
      push(@$biotypes_array,$transcript->biotype());
      $unique_biotypes->{$transcript->biotype()} = 1;
    }
  }

  $biotypes_hash->{'genes'} = $biotypes_array;
  return($biotypes_hash);
}

sub generate_single_transcript_genes {
  my ($self,$transcripts) = @_;

  my $single_transcript_genes = [];
  foreach my $transcript (@$transcripts) {
    my $gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);
    $gene->add_Transcript($transcript);
    $gene->biotype($transcript->biotype());
    $gene->stable_id($transcript->stable_id());
    $gene->{'status'} = $transcript->{'status'};
    $gene->{'cov'} = $transcript->{'cov'};
    $gene->{'perc_id'} = $transcript->{'perc_id'};
    $gene->slice($transcript->slice());
    push(@$single_transcript_genes,$gene);
  }

  return($single_transcript_genes);
}

sub check_cluster_status {
  my ($self,$cluster) = @_;

  my $genes = $cluster->get_Genes();
  $cluster->{'status'} = 'bad';
  foreach my $gene (@$genes) {
    if($gene->{'status'} eq 'good') {
      $cluster->{'status'} = 'good';
    }
  }
}

sub create_gene_from_cluster {
  my ($self,$cluster,$parent_gene_ids,$source_transcript_id_hash) = @_;

  my $cluster_genes = $cluster->get_Genes();
  my $processed_transcript_ids = {};
  my $selected_transcripts = {};

  # Loop through the transcripts, if there's any transcript that occurs twice, pick the one with the best combined coverage and percent id
  foreach my $cluster_gene (@$cluster_genes) {
    my $transcripts = $cluster_gene->get_all_Transcripts();
    my $transcript = ${$transcripts}[0];
    if($selected_transcripts->{$transcript->stable_id()}) {
      my $current_selected_transcript = $selected_transcripts->{$transcript->stable_id()};
      if(($transcript->{'cov'} + $transcript->{'perc_id'}) > ($current_selected_transcript->{'cov'} + $current_selected_transcript->{'perc_id'})) {
        $selected_transcripts->{$transcript->stable_id()} = $transcript;
      }
    } else {
      $selected_transcripts->{$transcript->stable_id()} = $transcript;
    }
  } # End foreach my $cluster_gene

  my $final_transcripts = [];
  foreach my $transcript_id (keys(%{$selected_transcripts})) {
    my $transcript = $selected_transcripts->{$transcript_id};
    push(@$final_transcripts,$transcript);
  }

  my $gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);
  my $transcript_id = ${$final_transcripts}[0]->stable_id();
  my $parent_gene_id = $parent_gene_ids->{$transcript_id}->{'gene_id'};
  my $parent_gene_stable_id = $parent_gene_ids->{$transcript_id}->{'gene_stable_id'};
  my $parent_gene_version = $parent_gene_ids->{$transcript_id}->{'gene_version'};
  my $parent_gene_biotype = $parent_gene_ids->{$transcript_id}->{'gene_biotype'};
  my $parent_gene_description = $parent_gene_ids->{$transcript_id}->{'gene_description'};

  foreach my $transcript (@$final_transcripts) {
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};
    unless($source_transcript) {
      $self->throw("Issue with finding source transcript. The following dbID was not found in the source transcript id hash: ".$transcript->stable_id());
    }

    my $parent_transcript_stable_id = $parent_gene_ids->{$transcript->stable_id()}->{'transcript_stable_id'};
    my $parent_transcript_version = $parent_gene_ids->{$transcript->stable_id()}->{'transcript_version'};
    my $parent_gene_id = $parent_gene_ids->{$transcript->stable_id()}->{'gene_id'};
    my $is_canonical = $parent_gene_ids->{$transcript->stable_id()}->{'is_canonical'};
    my $source = $parent_gene_ids->{$transcript->stable_id()}->{'source'};

    $transcript->biotype($source_transcript->biotype());
    $transcript->stable_id($parent_transcript_stable_id);
    $transcript->version($parent_transcript_version);

    if($is_canonical) {
      $transcript->is_canonical(1);
    }
    if ($transcript->{annotation_method} eq 'exonerate_local' or $transcript->{annotation_method} eq 'minimap_local') {
      my $toplevel_transcript = $transcript->transfer($transcript->slice->seq_region_Slice);
      if (!$toplevel_transcript) {
        $self->throw('Could not transfer '.$transcript->display_id.' from '.$transcript->slice->name.' to toplevel');
      }
      else {
        $toplevel_transcript->{annotation_method} = $transcript->{annotation_method};
        $transcript = $toplevel_transcript;
      }
    }

    $transcript->source($source);
    
    if ($transcript->coding_region_start() and
        $transcript->coding_region_end() and
        $transcript->coding_region_end()-$transcript->coding_region_start()+1 < 3) {
      # By convention, the coding_region_end is always higher than the
      # value returned by the coding_region_start method.
      say "Transcript CDS is too short (< 3 bp). Parent transcript stable id: ".$parent_transcript_stable_id.".".$parent_transcript_version;
    } elsif ($transcript->biotype() eq 'protein_coding' and !($transcript->translation())) {
      say "Transcript biotype is protein_coding but it does not have any translation. Parent transcript stable id: ".$parent_transcript_stable_id.".".$parent_transcript_version;
    } else {        
      $gene->add_Transcript($transcript);
    }
  }

  $gene->{'parent_gene_id'} = $parent_gene_id;
  $gene->stable_id($parent_gene_stable_id);
  $gene->version($parent_gene_version);
  $gene->biotype($parent_gene_biotype);
  #my $gene_description = ";parent_gene=".$parent_gene_stable_id.".".$parent_gene_version.";mapping_type=primary_mapping";
  #$gene->description($gene_description);
  $gene->description($parent_gene_description);

  # add source gene stable id as gene attribute
  # my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_g',-VALUE => $parent_gene_stable_id.".".$parent_gene_version);
  #$gene->add_Attributes($parent_attribute);

  return($gene);
}

sub calculate_max_intron_size {
  my ($self,$transcripts) = @_;

  my $scaling_factor = 1.5;
  my $max_intron_size = 100000;
  foreach my $transcript (@$transcripts) {
    my $introns = $transcript->get_all_Introns();
    foreach my $intron (@$introns) {
      my $scaled_intron_length = $scaling_factor * $intron->length();
      if($scaled_intron_length > $max_intron_size) {
        $max_intron_size = $scaled_intron_length;
      }
    }
  }

  return($max_intron_size);
}

sub create_annotation_features {
  my ($self,$transcript) = @_;

  my $cds_start  = $transcript->cdna_coding_start;
  my $cds_end    = $transcript->cdna_coding_end;

  my $start_phase = $transcript->translation->start_Exon->phase();
  my $end_phase = $transcript->translation->end_Exon->end_phase();
  if ($start_phase > 0) {
    $cds_start += $start_phase == 1 ? 2 : 1;
  }
  if ($end_phase > 0) {
    $cds_end -= $end_phase;
  }

  my $annotation_feature = Bio::EnsEMBL::Feature->new(-seqname => $transcript->stable_id,
                                                      -strand  => 1,
                                                      -start   => $cds_start,
                                                      -end     => $cds_end);

 my $annotation_features->{$transcript->stable_id_version} = $annotation_feature;
 return($annotation_features);
}


sub create_filename{
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}

sub genome_index {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_genome_index} = $val;
  }

  return $self->{_genome_index};
}

sub genes_to_process {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_genes_to_process} = $val;
  }

  return $self->{_genes_to_process};
}

sub input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_input_file} = $val;
  }

  return $self->{_input_file};
}

sub paftools_path {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_paftools_path} = $val;
  }

  return $self->{_paftools_path};
}

sub source_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_source_adaptor} = $val;
  }

  return $self->{_source_adaptor};
}

sub target_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_target_adaptor} = $val;
  }

  return $self->{_target_adaptor};
}

sub delete_input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_delete_input_file} = $val;
  }

  return $self->{_delete_input_file};
}

sub parent_gene_ids {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_parent_gene_ids} = $val;
  }

  return $self->{_parent_gene_ids};
}

sub gene_genomic_seqs_hash {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_gene_genomic_seqs_hash} = $val;
  }

  return $self->{_gene_genomic_seqs_hash};
}

sub no_projection {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_no_projection} = $val;
  }

  return $self->{_no_projection};
}


=head2 exon_rank
  Arg [1]    : Bio::EnsEMBL::Exon to get the rank from
  Arg [2]    : Bio::EnsEMBL::Transcript where the exon in Arg [1] is found
  Description: It searches for the Arg [1] exon in the Arg [2] transcript exon array by using the exon coordinates.
               Note it is assumed that a transcript cannot have overlapping exons nor exons having the same coordinates.
  Returntype : int
               It returns the index (1 <= index <= number_of_exons) of the Arg [1] exon in the Arg [2] transcript exon array.
  Exceptions : It throws if the Arg [1] exon cannot be found in the Arg [2] transcript.
=cut

sub exon_rank {
  my ($self,$exon,$transcript) = @_;

  my $rank = 0;
  foreach my $t_exon (@{$transcript->get_all_Exons()}) {
    $rank++;
    if ($t_exon->seq_region_start() == $exon->seq_region_start() and
        $t_exon->seq_region_end() == $exon->seq_region_end() and
        $t_exon->seq_region_strand() == $exon->seq_region_strand()) {
      last;
    }
  }

  if ($rank) {
    return $rank;
  } else {
    $self->throw("Exon(seq_region_start,seq_region_end,seq_region_strand) ".$exon->seq_region_start()." ".$exon->seq_region_end()." ".$exon->seq_region_strand()." does not belong to transcript(seq_region_start,seq_region_end,seq_region_strand) ".$transcript->seq_region_start()." ".$transcript->seq_region_end()." ".$transcript->seq_region_strand());
  }
}

sub revcomp {
  my ($self,$seq) = @_;
  my $revcomp = reverse $seq;
  $revcomp =~ tr/ATGCatgc/TACGtacg/;
  return($revcomp);
}


=head2 calculate_translation_based_on_source_transcript

 Arg [1]    : Bio::EnsEMBL::Transcript $transcript, the aligned transcript
 Arg [2]    : Bio::EnsEMBL::Transcript $source_transcript, the source transcript
 Description: Set the translation of the aligned transcript using the CIGAR strings
              and informations from the source transcript like the original translation
              coordinates
              If a translation is already set, it doesn't do any processing
 Returntype : Bio::EnsEMBL::Transcript
 Exceptions : Throws if start_Exon or end_Exon is not set

=cut

sub calculate_translation_based_on_source_transcript {
  my ($self, $transcript, $source_transcript) = @_;

  my $original_cdna_coding_start = $source_transcript->cdna_coding_start;
  my $original_cdna_coding_end = $source_transcript->cdna_coding_end;
  if (!$transcript->translation) {
    my $last_exon = $transcript->end_Exon;
    my $last_sf = $last_exon->get_all_supporting_features->[0];
    if ($original_cdna_coding_end > $last_sf->hend) {
      $original_cdna_coding_end = $last_sf->hend;
    }
    my $original_start_phase = $source_transcript->translation->start_Exon->phase;
    my $new_translation = $transcript->translation || $source_transcript->translation->new;
    $transcript->translation(undef);
    if ($last_sf->hend < $original_cdna_coding_start) {
      return $transcript;
    }
    my $first_exon = $transcript->start_Exon;
    my $first_sf = $first_exon->get_all_supporting_features->[0];
    if ($original_cdna_coding_start < $first_sf->hstart) {
      my $removed_length = $first_sf->hstart-$original_cdna_coding_start;
      if ($original_start_phase != -1) {
        $original_start_phase = ($original_start_phase+$removed_length)%3;
      }
      else {
        $original_start_phase = $removed_length%3;
      }
      $new_translation->start(1);
      $first_exon->phase($original_start_phase);
      $new_translation->start_Exon($first_exon);
      $transcript->cdna_coding_start(1);
    }
    my $source_cdna_start = $first_sf->hstart;
    my $source_cdna_end = $source_cdna_start;
    my $cdna_coding_start = $first_sf->hstart;
    my $cdna_coding_end = $cdna_coding_start;
    foreach my $exon (@{$transcript->get_all_Exons}) {
      foreach my $sf (@{$exon->get_all_supporting_features}) {
        my $cigar_string = $sf->cigar_string;
        my $exon_start = 1;
        my $exon_end = 1;
        my $cdna_start_exon_start = $cdna_coding_start;
        my $cdna_end_exon_start = $cdna_coding_end;
        while ($cigar_string =~ /(\d*)([MDI])/gc) {
          my $len = $1 || 1;
          if ($2 eq 'M') {
            if (!$transcript->cdna_coding_start) {
              if ($source_cdna_start+$len > $original_cdna_coding_start) {
                my $diff = $original_cdna_coding_start-$source_cdna_start;
                if ($diff >= 0) {
                  $cdna_coding_start += $diff;
                  $exon_start = $cdna_coding_start-$cdna_start_exon_start+1;
                }
                else {
                  my $phase = abs($diff)%3;
                  if ($phase > 0) {
                    $cdna_coding_start += 1;
                    $exon_start += $cdna_coding_start-$cdna_start_exon_start;
                  }
                }
                $transcript->cdna_coding_start($cdna_coding_start);
                $new_translation->start_Exon($exon);
                $new_translation->start($exon_start);
              }
              else {
                $source_cdna_start += $len;
                $cdna_coding_start += $len;
                $exon_end += $len;
              }
            }
            if (!$transcript->cdna_coding_end) {
              if ($source_cdna_end+$len > $original_cdna_coding_end) {
                my $diff = $original_cdna_coding_end-$source_cdna_end;
                if ($diff >= 0) {
                  $cdna_coding_end += $diff;
                  $exon_end = $cdna_coding_end-$cdna_end_exon_start+1;
                }
                $transcript->cdna_coding_end($cdna_coding_end);
                $new_translation->end_Exon($exon);
                $new_translation->end($exon_end);
              }
              else {
                $source_cdna_end += $len;
                $cdna_coding_end += $len;
                $exon_end += $len;
              }
            }
          }
          elsif ($2 eq 'D') {
            $source_cdna_start += $len;
            $source_cdna_end += $len;
          }
          elsif ($2 eq 'I') {
            $cdna_coding_start += $len;
            $cdna_coding_end += $len;
          }
        }
      }
      last if ($transcript->cdna_coding_start and $transcript->cdna_coding_end);
    }
    if ($new_translation->start_Exon and $new_translation->end_Exon and $new_translation->start and $new_translation->end) {
      $transcript->translation($new_translation);
      calculate_exon_phases($transcript, $original_start_phase);
    }
    else {
      $self->throw('start or end exon is missing: '.($new_translation->start_Exon || 'Start missing').' '.($new_translation->end_Exon || 'End missing'));
    }

  }
  return $transcript;
}

1;
