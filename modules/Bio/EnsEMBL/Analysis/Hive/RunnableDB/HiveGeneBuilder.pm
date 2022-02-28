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

Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGeneBuilder;

use warnings ;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(Gene_info attach_Analysis_to_Gene_no_ovewrite empty_Gene print_Gene_Transcript_and_Exons);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 
  qw(are_strands_consistent are_phases_consistent calculate_exon_phases
     is_not_folded all_exons_are_valid intron_lengths_all_less_than_maximum exon_overlap features_overlap overlap_length);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(exon_length_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils 
  qw(validate_Translation_coords contains_internal_stops print_Translation print_peptide);
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(cluster_Genes_by_coding_exon_overlap cluster_Genes);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Function  : instatiates a GeneBuilder object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Exceptions: 
  Example   : 

=cut



#sub new {
#  my ($class,@args) = @_;
 # my $self = $class->SUPER::new(@args);
 # $self->read_and_check_config($GENEBUILDER_CONFIG_BY_LOGIC);
 # return $self;
#}

sub param_defaults {
    return {
      post_filter_genes          => '1',
      recovery_overlap_threshold => 0.9,
      skip_readthrough_check => 1,
      load_all_biotypes => 1,
    }
}




sub fetch_input{
  my ($self) = @_;

  my $input_dba = $self->get_database_by_name('source_db');
  my $output_dba = $self->get_database_by_name('target_db');

  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    my $dna_dba = $self->hrdb_get_dba($self->param_required('dna_db'));
    $input_dba->dnadb($dna_dba);
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'output_db');
  $self->create_analysis;

  #fetch sequence
  my $slice = $input_dba->get_SliceAdaptor->fetch_by_name($self->param('iid'));
  $self->query($slice);

  #fetch genes
  if($self->param('load_all_biotypes')) {
    $self->input_genes($slice->get_all_Genes());
  } else {
    $self->get_Genes;
  }

  #print "Have ".@{$self->input_genes}." genes to cluster\n";
  #filter genes

  my $filtered_genes = $self->filter_genes($self->input_genes);
  #print "Have ".@filtered_genes." filtered genes\n";
  #create genebuilder runnable

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder
    ->new(
          -query => $self->query,
          -analysis => $self->analysis,
          -genes => $filtered_genes,
          -output_biotype => $self->OUTPUT_BIOTYPE,
          -max_transcripts_per_cluster => $self->MAX_TRANSCRIPTS_PER_CLUSTER,
          -min_short_intron_len => $self->MIN_SHORT_INTRON_LEN,
          -max_short_intron_len => $self->MAX_SHORT_INTRON_LEN,
          -blessed_biotypes => $self->BLESSED_BIOTYPES,
          -coding_only => $self->CODING_ONLY,
          -skip_readthrough_check => $self->param('skip_readthrough_check'),
         );


  $self->runnable($runnable);

};


sub run {
  my ($self) = @_;
  my $runnable = shift(@{$self->runnable()});
  $runnable->run;
  my $initial_genes = $runnable->output;
  say "Found ".scalar(@$initial_genes)." in initial runnable output";
  unless(scalar(@$initial_genes)) {
    $self->warning("No initial set of output genes created");
    return;
  }

  
  if($self->param('post_filter_genes')) {
    my $output_genes = $self->post_filter_genes($initial_genes);
    unless(scalar(@$output_genes)) {
      $self->throw("The output genes array is empty after running filter genes. This should not happen");
    }
    say "Found ".scalar(@$output_genes)." after post filtering";
    $self->output($output_genes);
  } else {
    $self->output($initial_genes);
  }
}

sub write_output{
  my ($self) = @_;

  my $output_dba = $self->hrdb_get_con('output_db');
  my $ga = $output_dba->get_GeneAdaptor();
  my $sucessful_count = 0;
  logger_info("WRITE OUTPUT have ".@{$self->output}." genes to write");
  foreach my $gene (@{$self->output}){
    attach_Analysis_to_Gene_no_ovewrite($gene, $self->analysis);
    empty_Gene($gene);
    eval{
      $ga->store($gene);
    };
    if($@){
      $self->warning("Failed to write gene ".id($gene)." ".coord_string($gene)." $@");
    }else{
      $sucessful_count++;
      logger_info("STORED GENE ".$gene->dbID);
    }
  }
  if($sucessful_count != @{$self->output}){
    $self->throw("Failed to write some genes");
  }
}

sub post_filter_genes {
  my ($self,$genes) = @_;

  # 1: Find longer models to fix split genes and extend modles from set of recovery dbs
  #    The downsides of this method is that it may join genes by accident. Also it could
  #    add extra transcripts where the difference is relatively minor. However, if the
  #    input recovery set is good then this will work well at fixing common issues with
  #    splitting genes due to layer annotation
#  $genes = $self->recover_long_models($genes);

  # 2: Remove short, suspect transcripts (leaving only the best scoring one if it removes everything)
  #    This will collapse stuff down at each position in terms of 1, 2 or 3 exon models. These models
  #    are most likely to be spurious. The script will remove these models in genes where they have
  #    >= 50 percent unique exons. In general works will with no obvious downside, unless there are
  #    true low exon models with a majority of unique exon structures (unlikely though)
  $genes = $self->remove_bad_transcripts($genes);

  # Haven't implemented the rest yet

  # 3: Remove UTR exons that cover the coding exons of another gene (where the whole exon is non-coding)
  #    Currently disabled as I haven't finished it and haven't come up with a good solution other than
  #    brute force looping through all the exons
  #$genes = $self->trim_utr_overlapping_cds($genes);

  # 4: Look for genes in strange places. For example small genes that are in introns of other genes
  #    Not implemented yet. Sometimes tricky as terminal exons are often misaligned. Might be worth
  #    ignoring any introns bordering a terminal exon when deciding if a gene lies in the intron of another

  return($genes);
}


sub remove_bad_transcripts {
  my ($self,$genes) = @_;

  my $suspect_exon_count = 3;
  my $suspect_exon_size = 30;

  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts();
    if (scalar(@$transcripts) > 1) {
      $gene->flush_Transcripts;
      # First count the frequency of the exons
      my $unique_exon_data;
      my $transcript_index = 0;
      foreach my $transcript (@$transcripts) {
        my $exons = $transcript->get_all_Exons();
        foreach my $exon (@$exons) {
          my $exon_string = ":".$exon->start.":".$exon->end.":";
          if (exists $unique_exon_data->{$exon_string}) {
            $unique_exon_data->{$exon_string}++;
          } else {
            $unique_exon_data->{$exon_string} = 1;
          }
        } # end foreach my $exon
      } # end  foreach my $transcript

      my $highest_score = 0;
      my $highest_scoring_transcript;
      # Now loop back through the transcript to find suspect ones
      foreach my $transcript (@$transcripts) {
        # If we've found a transcript with UTR then assume it's okay
        if($transcript->five_prime_utr || $transcript->three_prime_utr) {
          $gene->add_Transcript($transcript);
          next;
        }

        my $exons = $transcript->get_all_Exons;
        my $exon_count = scalar(@$exons);
        # When there's more than the suspect exon count then assume it's okay
        if($exon_count > $suspect_exon_count) {
          $gene->add_Transcript($transcript);
          next;
        }

        # Now we should a transcript that at least has some dodgy characteristics
        # We need to further investigate to decide if there is actually an issue
        my $unique_exon_count = 0;
        foreach my $exon (@$exons) {
          if ($unique_exon_data->{":".$exon->start.":".$exon->end.":"} == 1) {
            $unique_exon_count++;
          }
        }

        # At this point if we have 2 or more unique exons it implies the transcript
        # is very unusual and should probably be flagged for removal. I haven't thought
        # of a decent system for why we would keep it other than it being the best scoring
        # transcript in the scenario that all transcripts would be removed
        my $transcript_score = $self->transcript_score($transcript);
        if($transcript_score >= $highest_score) {
          $highest_scoring_transcript = $transcript;
          $highest_score = $transcript_score;
        }

        my $unique_fraction = $unique_exon_count / $exon_count;
        if($unique_fraction >= 0.5) {
          say "Skipping transcript because it has ".$suspect_exon_count." or fewer exons and ".$unique_exon_count." are unique";
        }
        else {
          $gene->add_Transcript($transcript);
        }

      } # end 2nd foreach my $transcript

      # If we're going to remove all the transcript then we should at least keep the highest scoring one
      if(scalar(@{$gene->get_all_Transcripts}) == 0) {
        say "Post filtering would remove all transcripts, keeping higest scoring transcript";
        $gene->add_Transcript($highest_scoring_transcript);
      }
    } # end foreach my $gene
  }

  return $genes;
}


=head2 transcript_score

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: Calculate a score for Arg[1] which is the sum of the percentage
              of identity with the hit coverage from each supporting feature.
              It selects the best score if there is more than one supporting
              features.
 Returntype : Int
 Exceptions : None

=cut

sub transcript_score {
  my ($self,$transcript) = @_;

  my $combined_cov_pid = 0;
  foreach my $tsf (@{$transcript->get_all_supporting_features}) {
    $combined_cov_pid = $tsf->hcoverage + $tsf->percent_id if ($tsf->hcoverage + $tsf->percent_id > $combined_cov_pid);
  }
  return $combined_cov_pid;
}


sub recover_long_models {
  my ($self,$genebuilder_genes) = @_;

  my $output_genes = [];

  my $recovery_dbs = $self->param_required('recovery_dbs');
  my $recovery_genes = [];
  say "Preparing to recover gene models";
  foreach my $db_info (@$recovery_dbs) {
    say "Fetching genes for ".$db_info->{'-dbname'};
    my $genes = $self->hrdb_get_dba($db_info)->get_GeneAdaptor->fetch_all_by_Slice($self->query);
    say "Found ".scalar(@$genes)." genes to consider for recovery";
    push(@$recovery_genes,@$genes);
  }

  say "Combining recovery genes under a single biotype";
  foreach my $recovery_gene (@$recovery_genes) {
    $recovery_gene->biotype('recovery');
  }

  my $types_hash;
  $types_hash->{'genebuilder'} = [$self->OUTPUT_BIOTYPE];
  $types_hash->{'recovery'} = ['recovery'];

  say "Clustering genebuilder genes with recovery genes...";
  my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap([@$genebuilder_genes,@$recovery_genes],$types_hash);
  say "...finished clustering genes for recovery";

  say "Found ".scalar(@$clusters)." clusters";
  say "Found ".scalar(@$unclustered)." unclustered genes, will not change";

  foreach my $unclustered (@$unclustered) {
     my $unclustered_genes = $unclustered->get_Genes_by_Set('genebuilder');
     foreach my $unclustered_gene (@$unclustered_genes) {
       push(@$output_genes,$unclustered_gene);
     }
  }

  foreach my $cluster (@$clusters) {
    my $recovery_cluster_genes = $cluster->get_Genes_by_Set('recovery');
    my $genebuilder_cluster_genes = $cluster->get_Genes_by_Set('genebuilder');

    say "In cluster:";
    say "  Recovery genes: ".scalar(@$recovery_cluster_genes);
    say "  Genebuilder genes: ".scalar(@$genebuilder_cluster_genes);
    my $genebuilder_input_genes;
    $genebuilder_input_genes = $self->assess_recovery_overlap($recovery_cluster_genes,$genebuilder_cluster_genes);

    if(scalar(@$genebuilder_input_genes)) {
      say "Inputting ".scalar(@$genebuilder_input_genes)." single transcript genes to new genebuilder run";
      say "Running genebuilder...";
      my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
                      -query => $self->query,
                      -analysis => $self->analysis,
                      -genes => $genebuilder_input_genes,
                      -output_biotype => $self->OUTPUT_BIOTYPE,
                      -max_transcripts_per_cluster => $self->MAX_TRANSCRIPTS_PER_CLUSTER,
                      -min_short_intron_len => $self->MIN_SHORT_INTRON_LEN,
                      -max_short_intron_len => $self->MAX_SHORT_INTRON_LEN,
                      -blessed_biotypes => $self->BLESSED_BIOTYPES,
                      -coding_only => $self->CODING_ONLY,
                      -skip_readthrough_check => 1,
                     );
      $runnable->run;
      say "Created ".scalar(@{$runnable->output})." output genes";
      push(@$output_genes,@{$runnable->output});
    } else {
      say "Not adding models from recovery set so will not run genebuilder again on cluster genes";
      push(@$output_genes,@$genebuilder_cluster_genes);
    }

  }
  return($output_genes);
}


sub assess_recovery_overlap {
  my ($self,$recovery_genes,$genebuilder_genes) = @_;

  my $output_genes = [];

  # Just in case we make multi transcript recovery genes in future, pass it through the single transcript gene method
  my $single_transcript_recovery_genes = $self->make_single_transcript_genes($recovery_genes);
  foreach my $recovery_gene (@$single_transcript_recovery_genes) {
    my $recovery_transcript = shift(@{$recovery_gene->get_all_Transcripts()});
    my $passed_overlap = 0;
    foreach my $genebuilder_gene (@$genebuilder_genes) {
      my $genebuilder_transcripts = $genebuilder_gene->get_all_Transcripts();
      foreach my $genebuilder_transcript (@$genebuilder_transcripts) {
        if($self->pass_overlap($recovery_transcript,$genebuilder_transcript)) {
          $passed_overlap = 1;
	}
      }
    }

    unless($passed_overlap) {
      push(@$output_genes,$recovery_gene);
    }
  }

  push(@$output_genes,@$genebuilder_genes);
  $output_genes = $self->make_single_transcript_genes($output_genes);
  return($output_genes);
}


sub make_single_transcript_genes {
  my ($self,$genes) = @_;

  my $single_transcript_genes =[];

  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts;
    foreach my $transcript (@$transcripts) {
      my $new_gene = new Bio::EnsEMBL::Gene;
      $new_gene->source($gene->source);
      $new_gene->biotype($gene->biotype);
      $new_gene->add_Transcript($transcript);
      push(@$single_transcript_genes,$new_gene);
    } # end foreach my $transcript
  } # end foreach my $gene

  return($single_transcript_genes);
}


sub pass_overlap {
  my ($self,$recovery_transcript,$genebuilder_transcript) = @_;

  my $pass_overlap = 0;
  my $overlap_threshold = $self->param_required('recovery_overlap_threshold');
  my $exon_overlap = exon_overlap($recovery_transcript,$genebuilder_transcript);
  say "Exon overlap: ".$exon_overlap;
  say "Rec Trans: ".$recovery_transcript->seq_region_start.":".$recovery_transcript->seq_region_end.":".
                    $recovery_transcript->strand.":".$recovery_transcript->length;
  say "Genb Trans: ".$genebuilder_transcript->seq_region_start.":".$genebuilder_transcript->seq_region_end.":".
                     $genebuilder_transcript->strand.":".$genebuilder_transcript->length;

  my $exon_non_overlap = abs($recovery_transcript->length()-$genebuilder_transcript->length());
  my $overlap_score = $exon_overlap - $exon_non_overlap;

  my $overlap_fraction = $exon_overlap / $recovery_transcript->length;
  say "Overlap fraction: ".$overlap_fraction;
  say "Overlap score: ".$overlap_score;

  if($overlap_fraction >= $overlap_threshold) {
    $pass_overlap = 1;
  }

  return($pass_overlap);
}


sub trim_utr_overlapping_cds {
  my ($self,$genebuilder_genes) = @_;

  my $output_genes = [];


  my $types_hash;
  $types_hash->{'genes'} = [$self->OUTPUT_BIOTYPE];

  say "Clustering genebuilder genes to check UTRs...";
  my ($clusters, $unclustered) = cluster_Genes($genebuilder_genes,$types_hash);
  say "...finished clustering genes";

  say "Found ".scalar(@$clusters)." clusters";
  say "Found ".scalar(@$unclustered)." unclustered genes, will not change";

  foreach my $unclustered (@$unclustered) {
     my $unclustered_genes = $unclustered->get_Genes();
     foreach my $unclustered_gene (@$unclustered_genes) {
       push(@$output_genes,$unclustered_gene);
     }
  }

  foreach my $cluster (@$clusters) {
    my $cluster_genes = $cluster->get_Genes();

    for(my $i=0; $i<scalar(@$cluster_genes)-1; $i++) {
      my $gene_a = ${$cluster_genes}[$i];
      for(my $j=$i+1; $j<scalar(@$cluster_genes); $j++) {
        my $gene_b = ${$cluster_genes}[$j];
        ($gene_a,$gene_b) = $self->compare_utr_exons($gene_a,$gene_b)
      }
      push(@$output_genes,$gene_a);
    }
    # Push the last gene as for $i runs up to the second last
    push(@$output_genes,${$cluster_genes}[$#$cluster_genes]);
  }

  return($output_genes);
}


sub compare_utr_exons {
#  my ($self,$gene_a,$gene_b) = @_;

#  say "Comparing the following genes:";
#  say "  Gene A: ".$gene_a->seq_region_start.":".$gene_a->seq_region_end.":".$gene_a->strand;
#  say "  Gene B: ".$gene_b->seq_region_start.":".$gene_b->seq_region_end.":".$gene_b->strand;

#  my $transcripts_a = $gene_a->get_all_Transcripts();
 # my $transcripts_b = $gene_b->get_all_Transcripts();
#  foreach my $transcript_a (@$transcripts_a) {
#    my $exons_a = $transcript_a->get_all_Exons();
#    my $coding_exons_a = $transcript_a->get_all_translateable_Exons();
#    my $indices_to_remove
#    foreach my $transcript_b (@$transcripts_b) {
#      unless($transcript_a->five_prime_utr || $transcript_a->three_prime_utr || $transcript_b->five_prime_utr || $transcript_b->three_prime_utr) {
 #       next;
 #     }
 #     # Need to think up a proper algorithm for this
 #   }
 # }

#  return($gene_a,$gene_b);
}

#sub output_db{
#  my ($self, $db) = @_;

#  if($db){
#    $self->param('_output_db',$db);
#  }

#  if(!$self->param('_output_db')){
#    my $db = $self->get_dbadaptor($self->OUTPUT_DB);
#    $self->param('_output_db',$db);
#  }

#  return $self->param('_output_db');
#}

sub get_Genes {
  my ($self) = @_;
  my @genes;

  my $slice = $self->query();

  foreach my $db_name(keys(%{$self->INPUT_GENES})){
    my $biotypes = $self->INPUT_GENES->{$db_name};
    foreach my $biotype(@$biotypes){
      my $genes = $slice->get_all_Genes_by_type($biotype);
      print "Retrieved ".@$genes." of type ".$biotype."\n";
      push(@genes, @$genes);
    }
  }
  $self->input_genes(\@genes);
}


sub input_genes {
  my ($self, $arg) = @_;

  unless($self->param_is_defined('_input_genes')) {
    $self->param('_input_genes',[]);
  }

  if($arg){
    if(!ref($arg) || ref($arg) ne 'ARRAY') {
      $self->throw("Need to pass input genes an arrayref not ".$arg);
    }
    push(@{$self->param('_input_genes')},@$arg);
  }
  return $self->param('_input_genes');
}

sub filter_genes{
  my ($self, $genes) = @_;
  $genes = $self->input_genes if(!$genes);
  print "Have ".@$genes." to filter\n";
  my @filtered;
  GENE:foreach my $gene (@$genes) {
    #throw("Genebuilder only works with one gene one transcript structures")
    #  if(@{$gene->get_all_Transcripts} >= 2);
    my $transcripts = $gene->get_all_Transcripts();
    unless(scalar(@$transcripts)) {
      $self->warning("Likely broken gene as no transcripts were recovered. Skipping");
      next GENE;
    }
    foreach my $transcript(@{$gene->get_all_Transcripts}) {
      my $exons = $transcript->get_all_Exons();
      unless(scalar(@$exons)) {
        $self->warning("Transcript doesn't have exons. Likely a failure during write. Skipping");
        next GENE;
      }

      if($self->validate_Transcript($transcript,$exons)) {
        push(@filtered, $gene);
        next GENE;
      } else {
        print Gene_info($gene)." is invalid skipping\n";
        next GENE;
      }
    }
  }
  return \@filtered;
}

sub validate_Transcript{
  my ($self, $transcript, $exons) = @_;

  my $is_valid = 0;
  #basic transcript validation
  unless(are_strands_consistent($transcript)) {
    print "Transcript has inconsistent strands. ";
    $is_valid++;
  }
  unless(are_phases_consistent($transcript)) {
    print "Transcript has inconsistent exon phases. ";
    $is_valid++;
  }
  unless(is_not_folded($transcript)) {
    print "Transcript seems to be folded (with secondary structure). ";
    $is_valid++;
  }
 EXON:foreach my $exon (@{$exons}){
    if(exon_length_less_than_maximum($exon, $self->MAX_EXON_LENGTH)) {
      next EXON;
    } else {
      print "Exon in transcript exceeds max length. ";
      $is_valid++;
      last EXON;
    }
  }
  if (@$exons > 3 and $transcript->translation and $transcript->translation->length < 100) {
    my $sum = 0;
    my $max = 0;
    foreach my $exon (@$exons) {
      $sum += $exon->length;
      $max = $exon->length if ($max < $exon->length);
    }
    if ($max < 50 and ($sum/@$exons < 20)) {
      print "Transcript has only small exons ";
      ++$is_valid;
    }
  }
  if(contains_internal_stops($transcript)) {
    print "Transcript contains internal stops. ";
    $is_valid++;
  }
  unless(validate_Translation_coords($transcript)) {
    print "Transcript contains invalid translation coords. ";
    $is_valid++;
  }
  return 0 if($is_valid >= 1);
  return 1;
}


#CONFIG METHODS

=head2 INPUT_GENES

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions:
  Example   :

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::GeneBuilder

sub INPUT_GENES {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('INPUT_GENES',$arg);
  }
  return $self->param_required('INPUT_GENES');
}

#sub OUTPUT_DB {
#  my ($self, $arg) = @_;
#  if(defined $arg){
#    $self->param('OUTPUT_DB',$arg);
#  }
#  return $self->param('OUTPUT_DB');
#}

sub OUTPUT_BIOTYPE {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('OUTPUT_BIOTYPE',$arg);
  }
  return $self->param('OUTPUT_BIOTYPE');
}

sub MAX_TRANSCRIPTS_PER_CLUSTER {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('MAX_TRANSCRIPTS_PER_CLUSTER',$arg);
  }
  return $self->param('MAX_TRANSCRIPTS_PER_CLUSTER');
}

sub MIN_SHORT_INTRON_LEN {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('MIN_SHORT_INTRON_LEN',$arg);
  }
  return $self->param('MIN_SHORT_INTRON_LEN');
}

sub MAX_SHORT_INTRON_LEN {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('MAX_SHORT_INTRON_LEN',$arg);
  }
  return $self->param('MAX_SHORT_INTRON_LEN');
}

sub BLESSED_BIOTYPES {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('BLESSED_BIOTYPES',$arg);
  }
  return $self->param('BLESSED_BIOTYPES');
}

sub MAX_EXON_LENGTH {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('MAX_EXON_LENGTH',$arg);
  }
  return $self->param('MAX_EXON_LENGTH');
}

sub CODING_ONLY {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('CODING_ONLY',$arg);
  }
  return $self->param('CODING_ONLY');
}

1;
