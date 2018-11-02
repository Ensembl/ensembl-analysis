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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAEvaluator - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAEvaluator;

use warnings ;
use strict;


use Bio::EnsEMBL::Analysis; 
use Bio::EnsEMBL::Analysis::Runnable::lincRNAEvaluator; 
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(print_Gene_Transcript_and_Exons empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::LincRNA qw(get_genes_of_biotypes_by_db_hash_ref) ;  


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAEvaluator
  Function  : instatiates a lincRNAEvaluator object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAEvaluator 
  Exceptions: 
  Example   : 

=cut


sub fetch_input{
  my ($self) = @_;

  # set up config
  $self->create_analysis;
  
  # Get lincRNA candidate genes and break each of them down into single-transcript genes: 
  my $lincrna_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->LINCRNA_DB);
  my $single_transcript_lincrna_genes = $self->create_single_transcript_genes($lincrna_genes);
  if (@$single_transcript_lincrna_genes) {
    print  'We have '.scalar(@$single_transcript_lincrna_genes).' single_transcript genes (1 gene.. 1 transcript), broken down from '.scalar(@$lincrna_genes)." multi_transcript lincRNA candidates from lincRNAFinder stage.\n";

    # Get ensembl genes as the validation set
    my $validation_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->VALIDATION_DBS);
    print  "\nFetched ".scalar(@$validation_genes)." Ensembl genes as the validation set.\n";

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::lincRNAEvaluator->new(
            -query => $self->query,
            -analysis => $self->analysis,
            -linc_rna_genes => $single_transcript_lincrna_genes,
            -ensembl_genes => $validation_genes,
             );

    $runnable->max_frameshift_intron_len($self->MAX_FRAMESHIFT_INTRON_LEN); #### at6 added 4 Mar
    $runnable->exclude_single_exon_lincRNAs($self->EXCLUDE_SINGLE_EXON_LINCRNAS);  #####at6 added 3 Mar
    $runnable->exclude_artefact_two_exon_lincRNAs($self->EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS);  ####at6 added 3 Mar

    $self->runnable($runnable);
  }
  else {
    $self->input_job->autoflow($self->param('_auto_flow'));
    $self->complete_early('No genes to process');
  }
}

sub run {
  my ($self) = @_;  
  foreach my $runnable (@{$self->runnable}) {
    $runnable->run();

    print  "\nCOLLAPSING SINGLE-TRANSCRIPT lincRNA 'GENES' BY GENEBUILDER TO CREATE FINAL lincRNA GENES.\n";

    #
    # First genebuilder run with unclustered lincRNAs
    #

    print  "\n Running GeneBuilder for ".scalar(@{$runnable->unclustered_ncrnas})." unclustered lincRNAs...\n";
    my $gb = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
            -query => $self->query,
            -analysis => $self->analysis,
            -genes => $runnable->unclustered_ncrnas,     #  \@genes_for_build,
            -output_biotype => $self->FINAL_OUTPUT_BIOTYPE,
            -max_transcripts_per_cluster => $self->MAX_TRANSCRIPTS_PER_CLUSTER,
            -min_short_intron_len => 1,
            -max_short_intron_len => 15,
            -blessed_biotypes => {},
           );
     $gb->run();

     my @output_clustered;
     print  "\n 2GeneBuilder returned ". scalar(@output_clustered) . " lincRNA genes which do not overlap with any other gene ";
     print  "(e.g. not proc_tran, existing lincRNAs or protein_coding genes).\n";
     push @output_clustered, @{$gb->clean_output()};

    #
    # OPTIONAL: second genebuilder run with lincRNAs which cluster with processed transcript
    #

    if ( $self->param('WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS') == 1	 ) {
      print "\n 3Running GeneBuilder for ".scalar(@{$runnable->ncrna_clusters_with_processed_transcript})." lincRNAs which cluster with processed transcript...\n";

      $gb->output_biotype('lincRNA_clusters_with_proc_trans');
      $gb->input_genes($runnable->ncrna_clusters_with_processed_transcript);
      $gb->run();
      print '  GeneBuilder returned '.scalar( @{$gb->output()} )." lincRNA genes which overlapped with processed_transcript genes.\n";
      push @output_clustered, @{$gb->clean_output()};
    }

    #
    # OPTIONAL: third genebuilder run with lincRNAs which cluster with existing lincRNAs in the core/SOURCE_PROTEIN_CODING DB
    #
    if ( $self->param('WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS') == 1 ) {
       print  "\n 5Running GeneBuilder for " . scalar(@{$runnable->ncrna_clusters_with_existing_lincRNAs}). " lincRNAs which cluster with existing lincRNAs...\n";
       $gb->output_biotype('lincRNA_clusters_with_existing_lincRNA');
       $gb->input_genes($runnable->ncrna_clusters_with_existing_lincRNAs);
       $gb->run();
       print  "\n GeneBuilder returned ".scalar( @{$gb->output()} )." lincRNA genes which overlapped with existing lincRNA genes.\n";
       push @output_clustered, @{$gb->clean_output()};
     }

     print  "\n GENEBUILDER RETURNED " .@output_clustered . " lincRNA GENES FOR WRITING (THIS DOES NOT INCLUDE REJECTED lincRNAs). THE TYPES OF lincRNA GENES WRITTEN DEPEND ON THE lincRNAEvaluator CONFIG SETTINGS.\n";
     $self->output( \@output_clustered );  #  This is just to store the lincRNA genes so test_RunnableDB script can find them.

  }
}

sub write_output{
  my ($self) = @_; 

  # $self->create_analysis;
  print  "\nWRITING RESULTS IN OUTPUT DB and/or VALIDATION DB... " . "\n";
  # update genes in the source db which cluster with processed_transcripts or lincRNAs
  # (if requested in the config file)
  if ( $self->param('MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB') == 0  ) { 
    my $db = $self->get_database_by_name('cdna_db');
    my $ga = $db->get_GeneAdaptor;
    my $update_analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($self->OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME);     
    $self->throw('Analysis with logic_name '.$self->OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME.' cannot be found in '.$db->dbname.'@'.$db->host."\n")
      unless (defined $update_analysis && $update_analysis->isa('Bio::EnsEMBL::Analysis'));
    foreach my $runnable (@{$self->runnable}) {
      my $proc_tran_genes_to_update = $runnable->proc_tran_genes_to_update;
      print  "  have " .scalar(@$proc_tran_genes_to_update ) . " processed_transcript genes to update in VALIDATION DB\n";

      for my $ug ( @$proc_tran_genes_to_update ) {
        # before we update the analysis logic_name and gene biotype, we check if the analysis
        # of gene processed_transcript is 'havana' or 'ensembl_havana_gene' to make sure no
        # Ensembl-only proc_trans genes slipped through the net (there shouldn't be any in
        # the first place!)

        my $hav_logic_name_to_match = $self->PROC_TRANS_HAVANA_LOGIC_NAME_STRING;
        if ($ug->analysis->logic_name =~ /$hav_logic_name_to_match/) {
          $ug->analysis($update_analysis);
          $ug->biotype('proc_trans_turned_lincRNA');
          $ga->update($ug);
          print  "  updated gene " . $ug->dbID . "\n";
        }
        else {
          $self->warning('not updating gene '.$ug->biotype.' with dbID '.$ug->dbID.' as it has the wrong analysis logic_name: '
          .$ug->analysis->logic_name.' ( to update, the logic_name should contain the string '.$self->PROC_TRANS_HAVANA_LOGIC_NAME_STRING.')');
        }
      }
    }
  }


  if ( $self->param('MARK_EXISTING_LINCRNA_IN_VALIDATION_DB') == 0  ) {
    my $db = $self->get_dbadaptor($self->UPDATE_SOURCE_DB);
    my $ga = $db->get_GeneAdaptor;
    my $update_analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($self->OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME);

    foreach my $runnable (@{$self->runnable}) {
      my $existing_lincRNA_genes_to_update = $runnable->old_lincRNA_genes_to_update;
      print  "  have " .scalar(@$existing_lincRNA_genes_to_update ) . " existing lincRNA genes to update in VALIDATION DB\n";

      foreach my $old_g (@$existing_lincRNA_genes_to_update ) {
        # No need to check sanity of analysis logic_name here because all "old" (existing) lincRNA genes would have
        # had logic_name "ensembl" or "ensembl_havana_merge". None of them should have logic_name "havana" because
        # Havana does not annotate lincRNAs at the gene level (but only at the transcript level)
        $old_g->biotype('lincRNA_common');
        $old_g->analysis($update_analysis);
        $ga->update($old_g);
        print  '  updated gene '.$old_g->dbID."\n";
      }
    }
  }
 
  my $dba = $self->hrdb_get_dba($self->param('output_db'));
  $self->hrdb_set_con($dba,'output_db');
  my $lincrna_ga  = $self->hrdb_get_con('output_db')->get_GeneAdaptor;
  my $genes_to_write = $self->output;

  if ( $self->WRITE_REJECTED_NCRNAS == 1 ) {
    foreach my $runnable (@{$self->runnable}) {
      print  "Writing rejected genes. They have not been collapsed by GeneBuilder so are still in one-gene-one-transcript format, as generated by lincRNAFinder.\n";
      print  "  mutliple-biotype rejects: " . scalar(@{$runnable->ncrna_clusters_with_multiple_ens_biotypes}) . " genes \n";
      print  "  protein-domain rejects: " . scalar(@{$runnable->ncrna_clusters_with_protein_domain}) . " genes \n";
      print  "  single-biotype rejects " . scalar(@{$runnable->ncrna_clusters_with_single_ens_biotype}) . " genes \n";
      push(@$genes_to_write, @{$runnable->ncrna_clusters_with_single_ens_biotype});
      push(@$genes_to_write, @{$runnable->ncrna_clusters_with_multiple_ens_biotypes});
      push(@$genes_to_write, @{$runnable->ncrna_clusters_with_protein_domain});
    }
  }  

  print  "***HAVE ". scalar(@$genes_to_write) ." GENE(S) TO WRITE IN TOTAL (INCLUDING VALID AND REJECTED lincRNAs).\n";

  my $sucessful_count = 0 ; 
  my $analysis = $self->analysis;
  foreach my $gene(@$genes_to_write){
    my $logic_name_to_be = "Hive_LincRNAEvaluator"; 
    my @t = @{ $gene->get_all_Transcripts}; 
    my $strand_to_use  =  $gene->strand;    
    my $biotype_to_use = "empty_something_is_wrong"; 
    if ( $gene->biotype eq $self->FINAL_OUTPUT_BIOTYPE ) {
      $biotype_to_use = $gene->biotype ; 
      # $logic_name_to_be = "passEval_to_final_check"; # I don't think there is a need to change the logic_name
    } else {
      $biotype_to_use = $t[0]->biotype ;
    }

    ## I will create a new gene without translations and only one transcript to be stored under different analysis ##
    # Make an analysis object (used later for storing stuff in the db)

    my $analysis = Bio::EnsEMBL::Analysis->new(
                                           -logic_name => $logic_name_to_be,
                                           -displayable => 1
                                           );
    my $analysis_adaptor = $dba->get_AnalysisAdaptor();

    # check if the logic name present in my database or if I need to create a new analysis_id, logic_name etc...
    if ($analysis_adaptor->fetch_by_logic_name($logic_name_to_be) ) {
    }else {
        print "# will store the analysis, since it is not exist \n";
        my $description = $logic_name_to_be ;
        my $display_label = $logic_name_to_be;
        $analysis->description($description) if $description;
        $analysis->display_label($display_label) if $display_label;
        $analysis_adaptor->store($analysis);
    }
    
    my $tag      = $gene->display_id; 
    my $gene_linc = Bio::EnsEMBL::Gene->new( 
                                       -analysis => $analysis,
                                       -biotype => $biotype_to_use, 
                                       -strand  => $strand_to_use, 
                                       );
    # print_Gene_Transcript_and_Exons($gene_linc); # if check
    
    foreach my $tr (@t) {
      my $exs      = []; # new array of exons
      $tr->translation(undef);	
  
      foreach my $ex_load (@{ $tr->get_all_Exons } ) {
        my $start_exon     = $ex_load->start();
        my $end_exon       = $ex_load->end();
        my $slice          = $ex_load->slice();
        my $q_strand       = $ex_load->strand; 
  
        my $ex_load        = Bio::EnsEMBL::Exon->new(                     
                                                                -start  => $start_exon,
                                                                -end    => $end_exon,
                                                                -strand => $q_strand, 
                                                                -analysis => $analysis,
                                                                -phase  => -1,
                                                                -end_phase => -1,
                                                                -slice  => $slice
                                               );
        push(@{$exs}, $ex_load); 
      }
      # Now LOAD transcript, gene and finally the db
      my $transcript = Bio::EnsEMBL::Transcript->new(
                                                 -exons => $exs, 
                                                 -analysis => $analysis, 
                                                 -biotype => $biotype_to_use,
                                                 -strand  => $strand_to_use,  
                                                 );
      $gene_linc->add_Transcript($transcript);    
    }

    # $gene_linc->status(undef); 

    # CHECK THIS AGAIN: we are going to have a different analysis for the output of this module. 
    # $gene_linc->analysis($self->analysis);   
    empty_Gene($gene_linc, 1);
    $lincrna_ga->store($gene_linc);
  } 
}



sub output_db{
  my ($self, $db) = @_; 
  if(defined $db && ref($db)=~m/Bio::EnsEMBL::DBSQL/ ){
    $self->{output_db} = $db;
  }
  if(!$self->{output_db}){
    $db = $self->get_dbadaptor($self->FINAL_OUTPUT_DB); 
    $self->{output_db} = $db; 
  }
  return $self->{output_db};
}


sub get_gene_sets {
  my ($self) = @_;
  my @genes; 

   my %sets_to_cluster = %{$self->CLUSTERING_INPUT_GENES};  

    if ( keys %sets_to_cluster != 2 ) { 
        $self->throw("you should only have 2 sets to cluster against - you can't cluster against more sets \n" ); 
    } 
 
  # check if hash-key name is correct : 
  unless ( exists $sets_to_cluster{"SET_1_CDNA"} && 
            exists $sets_to_cluster{"SET_2_PROT"}) {  
   $self->throw( " configuration error - I expect to get 2 sets of genes with names \"SET_1_CDNA\" and \"SET_2_PROT\" - check your config - you can't change these names!!!! \n" ) ; 
  }  
 
  foreach my $set ( keys %sets_to_cluster) { 
    my %this_set = %{$sets_to_cluster{$set}};      
    my @genes_in_set; 
    foreach my $database_db_name ( keys ( %this_set)) {  
       my $set_db = $self->get_dbadaptor($database_db_name);
       my $slice = $self->fetch_sequence($self->input_id, $set_db);  
       my @biotypes = @{$this_set{$database_db_name}}; 
       for my $biotype  ( @biotypes ) { 
          my $genes = $slice->get_all_Genes_by_type($biotype,undef,1);
          if ( @$genes == 0 ) {  
            $self->warning("No genes of biotype $biotype found in $set_db\n"); 
          } 
          print  "Retrieved ".@$genes." of type ".$biotype."\n";
          push @genes_in_set, @$genes; 
      }  
    } 
    my @single_transcript_genes = @{$self->create_single_transcript_genes(\@genes_in_set)};
   $self->all_gene_sets($set,\@genes_in_set); 
  }   
}


sub create_single_transcript_genes{
  my ($self, $genes) = @_; 

  print  "\nCreating single_transcript genes (In other words, one transcript per gene) ....\n" ;
  my @single_transcript_genes; 

 GENE:foreach my $gene(@$genes){ 
 	
 	my @tr = @{$gene->get_all_Transcripts};

    if ( @tr == 1 ) {  
      push @single_transcript_genes, $gene ; 
    } else { 
      foreach my $transcript(@tr ) { 
        my $ng = Bio::EnsEMBL::Gene->new();  
        $ng->add_Transcript($transcript); 
        push @single_transcript_genes, $ng ; 
      }
    }
  }
  return \@single_transcript_genes; 
}






#CONFIG METHODS


# this function checks if everything stored successfully 
sub check_if_all_stored_correctly { 
  my ($self) = @_; 

  my $set_db = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  if($dna_dba) { 
    $set_db->dnadb($dna_dba); 
  } 
  
  my $test_id = $self->param('iid'); 
  my $slice = $self->fetch_sequence($test_id, $set_db);
  print  "check if all genes are fine!! \n" ; 
  my $genes = $slice->get_all_Genes(undef,undef,1) ; 
	return "yes"; 
}




=head2 FINAL_OUTPUT_DB

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::lincRNAEvaluator
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::lincRNAEvaluator

sub FINAL_OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FINAL_OUTPUT_DB', $arg);
  }
  return $self->param('FINAL_OUTPUT_DB');
}  


sub FINAL_OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FINAL_OUTPUT_BIOTYPE', $arg);
  }
  return $self->param('FINAL_OUTPUT_BIOTYPE');
}  


sub VALIDATION_DBS{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('VALIDATION_DBS', $arg);
  }
  return $self->param('VALIDATION_DBS');
} 


sub EXCLUDE_SINGLE_EXON_LINCRNAS {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('EXCLUDE_SINGLE_EXON_LINCRNAS', $arg);
  }
  return $self->param('EXCLUDE_SINGLE_EXON_LINCRNAS');

}

sub MAX_FRAMESHIFT_INTRON_LEN {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MAX_FRAMESHIFT_INTRON_LEN', $arg);
  }
  return $self->param('MAX_FRAMESHIFT_INTRON_LEN');

}

sub EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS', $arg);
  }
  return $self->param('EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS');
}

sub MAX_TRANSCRIPTS_PER_CLUSTER {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('MAX_TRANSCRIPTS_PER_CLUSTER', $arg);
  }
  return $self->param('MAX_TRANSCRIPTS_PER_CLUSTER');
}

sub MARK_EXISTING_LINCRNA_IN_VALIDATION_DB{
  my ($self, $arg) = @_;
  if($arg) {
    $self->param('MARK_EXISTING_LINCRNA_IN_VALIDATION_DB', $arg);
  }
  return $self->param('MARK_EXISTING_LINCRNA_IN_VALIDATION_DB');
}


sub WRITE_REJECTED_NCRNAS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('WRITE_REJECTED_NCRNAS', $arg);
  }
  return $self->param('WRITE_REJECTED_NCRNAS');
}  


sub UPDATE_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('UPDATE_SOURCE_DB', $arg);
  }
  return $self->param('UPDATE_SOURCE_DB');
} 

sub WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS', $arg);
  }
  return $self->param('WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS');
} 

sub WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS', $arg);
  }
  return $self->param('WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS');
}


sub MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB {
  my ($self, $arg) = @_;
  
  if($arg){
    $self->param('MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB', $arg);
  }
  return $self->param('MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB');
} 


sub LINCRNA_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('LINCRNA_DB', $arg);
  }
  return $self->param('LINCRNA_DB');
}  

sub OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME', $arg);
  }
  return $self->param('OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME');
} 


sub PROC_TRANS_HAVANA_LOGIC_NAME_STRING{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->param('PROC_TRANS_HAVANA_LOGIC_NAME_STRING', $arg);
  }
  return $self->param('PROC_TRANS_HAVANA_LOGIC_NAME_STRING');
} 


1;
