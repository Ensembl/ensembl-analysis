=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
     is_not_folded all_exons_are_valid intron_lengths_all_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(exon_length_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils 
  qw(validate_Translation_coords contains_internal_stops print_Translation print_peptide);
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(make_types_hash_with_genes cluster_Genes);

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
      post_filter_genes => '1',
    }
}




sub fetch_input{
  my ($self) = @_;

  my $input_dba = $self->hrdb_get_dba($self->param('layering_output_db'));
  my $output_dba = $self->hrdb_get_dba($self->param('genebuilder_output_db'));

  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));


  if($dna_dba) {
    $input_dba->dnadb($dna_dba);
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($input_dba,'input_db');
  $self->hrdb_set_con($output_dba,'output_db');

  $self->create_analysis;

  #fetch sequence
  my $slice = $input_dba->get_SliceAdaptor->fetch_by_name($self->param('iid'));
  $self->query($slice);


  #fetch genes
  $self->get_Genes;

  #print "Have ".@{$self->input_genes}." genes to cluster\n";
  #filter genes

  my @filtered_genes = @{$self->filter_genes($self->input_genes)};
  #print "Have ".@filtered_genes." filtered genes\n";
  #create genebuilder runnable

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder
    ->new(
          -query => $self->query,
          -analysis => $self->analysis,
          -genes => \@filtered_genes,
          -output_biotype => $self->OUTPUT_BIOTYPE,
          -max_transcripts_per_cluster => $self->MAX_TRANSCRIPTS_PER_CLUSTER,
          -min_short_intron_len => $self->MIN_SHORT_INTRON_LEN,
          -max_short_intron_len => $self->MAX_SHORT_INTRON_LEN,
          -blessed_biotypes => $self->BLESSED_BIOTYPES,
          -coding_only => $self->CODING_ONLY,
         );


  $self->runnable($runnable);

};


sub run {
  my ($self) = @_;
  my $runnable = shift(@{$self->runnable()});
  $runnable->run;
  my $initial_genes = $runnable->output;
  if($self->param('post_filter_genes')) {
    my $output_genes = $self->post_filter_genes($initial_genes);
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

  # This has a few jobs to do
  # 1: find genes where there are multiple transcripts that don't really agree, with very little useful info
  #    This comes out of genblast mainly where it should build a single exon model but tries too hard to find
  #    start/end codons and then builds a 2/3 exon model with spurious terminal exons. We want to reduce these
  #    sites down to a single transcript
  # 2: trim the UTR in cases where the UTR of one gene extends into the coding region of another
  # 3: look for genes in strange places. For example small genes that are in introns of other genes
  # 4: look at overlapping genblast select models and if they are significantly more exon rich than the model
  #    they over lap we want to add in the genblast select model and then rerun genebuilder

  # 1: Remove short, suspect transcripts (leaving only the best scoring one if it removes everything)
  $genes = $self->remove_bad_transcripts($genes);

  # Haven't implemented the rest yet

  # 2: Remove UTR exons that cover the coding exons of another gene (where the whole exon is non-coding)
  #$genes = $self->trim_overlapping_utr($genes);

  # 3: Find genblast_select models
  #$genes = $self->recover_genblast_select($genes);
  return($genes);
}


sub remove_bad_transcripts {
  my ($self,$genes) = @_;

  my $suspect_exon_count = 3;
  my $suspect_exon_size = 30;

  my $filtered_genes = [];
  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts();
    if(scalar(@$transcripts) == 1) {
      push(@$filtered_genes,$gene);
      next;
    }

    # First count the frequency of the exons
    my $unique_exon_data;
    my $transcript_index = 0;
    foreach my $transcript (@$transcripts) {
      my $exons = $transcript->get_all_Exons();
      foreach my $exon (@$exons) {
        my $exon_string = ":".$exon->start.":".$exon->end.":";
        unless($unique_exon_data->{$exon_string}) {
          $unique_exon_data->{$exon_string} = 1;
        } else {
          $unique_exon_data->{$exon_string}++;
	}
      } # end foreach my $exon
    } # end  foreach my $transcript

    my $transcripts_to_keep = [];
    my $highest_score = 0;
    my $highest_scoring_transcript;
    # Now loop back through the transcript to find suspect ones
    foreach my $transcript (@$transcripts) {
      # If we've found a transcript with UTR then assume it's okay
      if($transcript->five_prime_utr || $transcript->three_prime_utr) {
        push(@{$transcripts_to_keep},$transcript);
        next;
      }

      my $exons = $transcript->get_all_Exons;
      my $exon_count = scalar(@$exons);
      # When there's more than the suspect exon count then assume it's okay
      if($exon_count > $suspect_exon_count) {
        push(@{$transcripts_to_keep},$transcript);
        next;
      }

      # Now we should a transcript that at least has some dodgy characteristics
      # We need to further investigate to decide if there is actually an issue
      my $unique_exon_count = 0;
      foreach my $exon (@$exons) {
        my $exon_string = ":".$exon->start.":".$exon->end.":";
        unless($unique_exon_data->{$exon_string} > 1) {
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
        next;
      }

      push(@{$transcripts_to_keep},$transcript);
    } # end 2nd foreach my $transcript

    # If we're going to remove all the transcript then we should at least keep the highest scoring one
    if(scalar(@$transcripts_to_keep) == 0) {
      say "Post filtering would remove all transcripts, keeping higest scoring transcript";
      push(@{$transcripts_to_keep},$highest_scoring_transcript);
    }

    my $new_gene = new Bio::EnsEMBL::Gene;
    $new_gene->source($gene->source);
    $new_gene->biotype($gene->biotype);

    unless(scalar(@$transcripts_to_keep)) {
      $self->throw("In post filtering no transcripts were kept for current gene, this is not right");
    }

    foreach my $transcript (@$transcripts_to_keep) {
      $new_gene->add_Transcript($transcript);
    }

    push(@$filtered_genes,$new_gene);
  } # end foreach my $gene

  unless(scalar(@$filtered_genes) == scalar(@$genes)) {
    $self->throw("The count of the filtered genes (".scalar(@$filtered_genes).") did not match the count of the original set (".scalar(@$genes).")");
  }

  return($filtered_genes);
}


sub transcript_score {
  my ($self,$transcript) = @_;

  my $tsf = pop(@{$transcript->get_all_supporting_features});
  unless($tsf) {
    return(0);
  }

  my $combined_cov_pid = $tsf->hcoverage + $tsf->percent_id;
  return($combined_cov_pid);

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
  GENE:foreach my $gene(@$genes) {
    #throw("Genebuilder only works with one gene one transcript structures")
    #  if(@{$gene->get_all_Transcripts} >= 2);
    foreach my $transcript(@{$gene->get_all_Transcripts}) {
      my $exons = $transcript->get_all_Exons();
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
