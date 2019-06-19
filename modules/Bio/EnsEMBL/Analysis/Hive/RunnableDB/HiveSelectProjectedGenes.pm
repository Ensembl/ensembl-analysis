# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectProjectedGenes

=cut

=head1 DESCRIPTION

HiveSelectProjectedGenes fetches the genes from the WGA and the CESAR2.0 projection output databases and
it selects one of the two overlapping projected transcripts by maximum cov+pid.
It stores the selected single-transcript genes into the output database together with the non-overlapping
projected transcripts from each projection method.   

=head1 OPTIONS

-dna_db        Ensembl database containing the DNA sequences.
-wga_db        Ensembl database containing the WGA projected genes.
-cesar_db      Ensembl database containing the CESAR2.0 projected genes.
-output_db     Ensembl database where the selected projected genes will be stored.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectProjectedGenes;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(empty_Transcript);
use DBI qw(:sql_types);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
# set database connections

  my $self = shift;
  
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  $self->hrdb_set_con($dna_dba,'dna_db');
  
  my $wga_dba = $self->hrdb_get_dba($self->param('wga_db'));
  $self->hrdb_set_con($wga_dba,'wga_db');
  $wga_dba->dnadb($dna_dba);
  
  my $cesar_dba = $self->hrdb_get_dba($self->param('cesar_db'));
  $self->hrdb_set_con($cesar_dba,'cesar_db');
  $cesar_dba->dnadb($dna_dba);
  
  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  $self->hrdb_set_con($output_dba,'output_db');
  $output_dba->dnadb($dna_dba);

  return 1;
}

sub run {
  my $self = shift;
  
  my $wga_ta = $self->hrdb_get_con('wga_db')->get_TranscriptAdaptor();
  my $wga_ga = $self->hrdb_get_con('wga_db')->get_GeneAdaptor();
  my $cesar_ta = $self->hrdb_get_con('cesar_db')->get_TranscriptAdaptor();
  my $cesar_ga = $self->hrdb_get_con('cesar_db')->get_GeneAdaptor();
  my $output_genes;
  my %discarded_wga;

  # loop through all cesar projected transcripts
  foreach my $cesar_t (@{$cesar_ta->fetch_all_by_logic_name('cesar')}) {
    my $orig_cesar_sid = $cesar_t->stable_id();
    my $cesar_sid = $orig_cesar_sid;
    $cesar_sid =~ s/\..*//; # remove stable id version if present
    
    # fetch the equivalent wga projected transcript
    my $wga_t = fetch_transcript_by_stable_id_and_logic_name($wga_ta,$cesar_sid,'project_transcripts');
    
    if ($wga_t) {
      # select the best t between wga and cesar based on cov and pid
      my @wga_sfs = @{$wga_t->get_all_supporting_features()};
      my $wga_sf = shift(@wga_sfs); # fetch any sf as all the sfs have the same coverage and percent id
      my $wga_cov = 0;
      my $wga_pid = 0;
      if ($wga_sf) {
        $wga_cov = $wga_sf->hcoverage();
        $wga_pid = $wga_sf->percent_id();
      }
      my @cesar_sfs = @{$cesar_t->get_all_supporting_features()};
      my $cesar_sf = shift(@cesar_sfs); # fetch any sf as all the sfs have the same coverage and percent id
      my $cesar_cov = 0;
      my $cesar_pid = 0;
      if ($cesar_sf) {
        $cesar_cov = $cesar_sf->hcoverage();
        $cesar_pid = $cesar_sf->percent_id();
      }
      
      if ($wga_cov+$wga_pid > $cesar_cov+$cesar_pid) {
        # select wga if wga cov+pid is greater than cesar cov+pid
        my $selected_gene = fetch_gene_by_transcript_stable_id_and_logic_name($wga_ga,$wga_ta,$cesar_sid,'project_transcripts');
        $selected_gene->flush_Transcripts();
        
        empty_Gene($selected_gene);
        
        $selected_gene->analysis($wga_t->analysis());
        $selected_gene->biotype($wga_t->biotype());
        empty_Transcript($wga_t);
        $selected_gene->add_Transcript($wga_t);
        push(@{$output_genes},$selected_gene);
        
        my $non_selected_gene = fetch_gene_by_transcript_stable_id_and_logic_name($cesar_ga,$cesar_ta,$orig_cesar_sid,'cesar');
        $non_selected_gene->flush_Transcripts();
        
        empty_Gene($non_selected_gene);
        
        $non_selected_gene->analysis($cesar_t->analysis());
        $non_selected_gene->biotype($cesar_t->biotype().'_cesar');
        empty_Transcript($cesar_t);
        $non_selected_gene->add_Transcript($cesar_t);
        push(@{$output_genes},$non_selected_gene);
      } else {
        # select cesar
        my $selected_gene = fetch_gene_by_transcript_stable_id_and_logic_name($cesar_ga,$cesar_ta,$orig_cesar_sid,'cesar');
        $selected_gene->flush_Transcripts();
        
        empty_Gene($selected_gene);
        
        $selected_gene->analysis($cesar_t->analysis());
        $selected_gene->biotype($cesar_t->biotype());
        empty_Transcript($cesar_t);
        $selected_gene->add_Transcript($cesar_t);
        push(@{$output_genes},$selected_gene);
        
        my $non_selected_gene = fetch_gene_by_transcript_stable_id_and_logic_name($wga_ga,$wga_ta,$cesar_sid,'project_transcripts');
        $non_selected_gene->flush_Transcripts();
        
        empty_Gene($non_selected_gene);
        
        $non_selected_gene->analysis($wga_t->analysis());
        $non_selected_gene->biotype($wga_t->biotype().'_wga');
        empty_Transcript($wga_t);
        $non_selected_gene->add_Transcript($wga_t);
        push(@{$output_genes},$non_selected_gene);
        $discarded_wga{$cesar_sid} = 1;
      }
      
    } else {
      # no wga equivalent means cesar is selected
      my $selected_gene = fetch_gene_by_transcript_stable_id_and_logic_name($cesar_ga,$cesar_ta,$orig_cesar_sid,'cesar');
      $selected_gene->flush_Transcripts();
      
      empty_Gene($selected_gene);
      
      $selected_gene->analysis($cesar_t->analysis());
      $selected_gene->biotype($cesar_t->biotype());
      empty_Transcript($cesar_t);
      $selected_gene->add_Transcript($cesar_t);
      push(@{$output_genes},$selected_gene);
    }
  }
  
  # loop through all wga projected transcripts to select the ones which did not have any cesar equivalent
  foreach my $wga_t (@{$wga_ta->fetch_all_by_logic_name('project_transcripts')}) {
    my $wga_t_sid = $wga_t->stable_id();
    if (!(exists $discarded_wga{$wga_t_sid})) {
      my $selected_gene = fetch_gene_by_transcript_stable_id_and_logic_name($wga_ga,$wga_ta,$wga_t_sid,'project_transcripts');
      $selected_gene->flush_Transcripts();
      
      empty_Gene($selected_gene);
      
      $selected_gene->analysis($wga_t->analysis());
      $selected_gene->biotype($wga_t->biotype());
      empty_Transcript($wga_t);
      $selected_gene->add_Transcript($wga_t);
      push(@{$output_genes},$selected_gene);
    }
  }

  $self->output($output_genes);
  return 1;
}

sub write_output {
  my $self = shift;

  my $output_genes = $self->output();
  my $gene_adaptor = $self->hrdb_get_con('output_db')->get_GeneAdaptor();
  foreach my $output_gene (@{$output_genes}) {
    $gene_adaptor->store($output_gene);
  }
  return 1;
}

sub fetch_transcript_by_stable_id_and_logic_name {
  my ($ta,$stable_id,$logic_name) = @_;

  my $constraint = "t.stable_id = ? AND t.analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = ?) ";
  $ta->bind_param_generic_fetch($stable_id,SQL_VARCHAR);
  $ta->bind_param_generic_fetch($logic_name,SQL_VARCHAR);
  my ($transcript) = @{$ta->generic_fetch($constraint)};
  
  return $transcript;
}

sub fetch_gene_by_transcript_stable_id_and_logic_name {
  my ($ga,$ta,$stable_id,$logic_name) = @_;

  my $transcript = fetch_transcript_by_stable_id_and_logic_name($ta,$stable_id,$logic_name);
  my $gene = $ga->fetch_by_transcript_id($transcript->dbID());
  
  return $gene;
}

1;
