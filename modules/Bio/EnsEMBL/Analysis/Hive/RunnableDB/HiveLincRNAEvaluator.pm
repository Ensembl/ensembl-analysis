=head1 LICENSE
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use vars qw(@ISA);
use strict;
use Data::Dumper;


use Bio::EnsEMBL::Hive::Utils ('destringify');
use Bio::EnsEMBL::Analysis; 
# use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
# use Bio::EnsEMBL::Analysis::Config::GeneBuild::lincRNAEvaluator; 
use Bio::EnsEMBL::Analysis::Runnable::lincRNAEvaluator; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange); 
# use Bio::EnsEMBL::Analysis::RunnableDB; 
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(print_Gene_Transcript_and_Exons); 

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
  $self->hive_set_config;
  
  # Get lincRNA candidate genes and break each of them down into single-transcript genes:
  print  "DEBUG:: fetch_input:: dump the object START_FROM" .  "\n"; 
  # get all candidates lincRNA 
  my $lincrna_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->LINCRNA_DB);
  ### ### foreach my $gene (@{$lincrna_genes}) {
  ### ###   print  $gene->display_id, "\t", $gene->adaptor->dbc->dbname, , "\t", $gene->adaptor->dbc->host, "\n";
  ### ### }

  my @single_transcript_lincrna_genes = @{$self->create_single_transcript_genes($lincrna_genes)}; 
  print  "We have ". scalar(@single_transcript_lincrna_genes) . " single_transcript genes (1 gene.. 1 transcript), broken down from " . scalar(@$lincrna_genes) . " multi_transcript lincRNA candidates from lincRNAFinder stage.\n" ;  

  # Get ensembl genes as the validation set  
  my $validation_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->VALIDATION_DBS);  
  print  "\nFetched ". scalar(@$validation_genes) . " Ensembl genes as the validation set.\n";

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::lincRNAEvaluator->new(
          -query => $self->query,
          -analysis => $self->analysis,
          -linc_rna_genes => \@single_transcript_lincrna_genes , 
          -ensembl_genes => $validation_genes,  
           );

  $runnable->max_frameshift_intron_len($self->MAX_FRAMESHIFT_INTRON_LEN); #### at6 added 4 Mar
  $runnable->exclude_single_exon_lincRNAs($self->EXCLUDE_SINGLE_EXON_LINCRNAS);  #####at6 added 3 Mar
  $runnable->exclude_artefact_two_exon_lincRNAs($self->EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS);  ####at6 added 3 Mar

  $self->single_runnable($runnable);     
  # not using RunnableDB SUPER class method "runnable" as the method stores runnable objects in an array.
  # the SUPER "run" method in RunnableDB then loops through the runnables.  As we are implementing the
  # 'run' method specifically for lincRNAEvaluator anyway, we might as well skip the array looping. 
}

sub run {
  my ($self) = @_;  
  $self->single_runnable->run(); 
 
  print  "\nCOLLAPSING SINGLE-TRANSCRIPT lincRNA 'GENES' BY GENEBUILDER TO CREATE FINAL lincRNA GENES.\n";

  #
  # First genebuilder run with unclustered lincRNAs
  #

  print  "\n 1Running GeneBuilder for " . scalar(@{$self->single_runnable->unclustered_ncrnas}). " unclustered lincRNAs...\n";  

  my $gb = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
          -query => $self->query,
          -analysis => $self->analysis,
          -genes => $self->single_runnable->unclustered_ncrnas, #  \@genes_for_build, 
          -output_biotype => $self->FINAL_OUTPUT_BIOTYPE, 
          -max_transcripts_per_cluster => $self->MAX_TRANSCRIPTS_PER_CLUSTER, 
          -min_short_intron_len => 1,
          -max_short_intron_len => 15,
          -blessed_biotypes => {} , 
         );
   $gb->run();

   my @output_clustered; 
   
# print "DEBUG::HIVElincRNaEvaluator::dumper::" . Dumper($gb) . "\n"; 
   push @output_clustered, @{$gb->output()} ;  
   print  "\n 2GeneBuilder returned ". scalar(@output_clustered) . " lincRNA genes which do not overlap with any other gene ";
   print  "(e.g. not proc_tran, existing lincRNAs or protein_coding genes).\n";

  # 
  # OPTIONAL: second genebuilder run with lincRNAs which cluster with processed transcript  
  # 
  
  my $temp_test= "no"; 
  if ( $self->param('WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS') == 1	 ) {
  # if ($temp_test eq "no") { 
     print "\n 3Running GeneBuilder for " . scalar(@{$self->single_runnable->ncrna_clusters_with_processed_transcript}) .  " lincRNAs which cluster with processed transcript...\n" ; 

     $gb = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
            -query => $self->query,
            -analysis => $self->analysis,
            -genes => $self->single_runnable->ncrna_clusters_with_processed_transcript,
            -output_biotype => "lincRNA_clusters_with_proc_trans", 
            -max_transcripts_per_cluster => $self->MAX_TRANSCRIPTS_PER_CLUSTER, 
            -min_short_intron_len => 1,
            -max_short_intron_len => 15,
            -blessed_biotypes => {} , 
           );
     $gb->run();
     # print "DEBUG::HIVElincRNaEvaluator::dumper::" . Dumper($gb) . "\n"; 
 
     print "  GeneBuilder returned ". scalar( @{$gb->output()} ) . " lincRNA genes which overlapped with processed_transcript genes.\n";
     push @output_clustered, @{$gb->output()} ;   
   }

  # 
  # OPTIONAL: third genebuilder run with lincRNAs which cluster with existing lincRNAs in the core/SOURCE_PROTEIN_CODING DB  
  # 
  # if ($temp_test eq "ye") { 
  if ( $self->param('WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS') == 1 ) {
     print  "\n 5Running GeneBuilder for " . scalar(@{$self->single_runnable->ncrna_clusters_with_existing_lincRNAs}). " lincRNAs which cluster with existing lincRNAs...\n" ;
     $gb = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
            -query => $self->query,
            -analysis => $self->analysis,
            -genes => $self->single_runnable->ncrna_clusters_with_existing_lincRNAs,
            -output_biotype => "lincRNA_clusters_with_existing_lincRNA",
            -max_transcripts_per_cluster => $self->MAX_TRANSCRIPTS_PER_CLUSTER,
            -min_short_intron_len => 1,
            -max_short_intron_len => 15,
            -blessed_biotypes => {} ,
           );
     $gb->run();
     print  "\n 6GeneBuilder returned ". scalar( @{$gb->output()} ) . " lincRNA genes which overlapped with existing lincRNA genes.\n";
     push @output_clustered, @{$gb->output()} ;
   }

   print  "\n 7GENEBUILDER RETURNED " .@output_clustered . " lincRNA GENES FOR WRITING (THIS DOES NOT INCLUDE REJECTED lincRNAs). THE TYPES OF lincRNA GENES WRITTEN DEPEND ON THE lincRNAEvaluator CONFIG SETTINGS.\n" ;  
   
   print  "DEBUG::HIVElincRNA::run XXXXXXXXXXXX\n";
   # genes_to_write(); 
   # $self->genes_to_write( \@output_clustered );  # The write_output method takes lincRNA genes from $self->genes_to_write
   # $self->throw("don't let it finish and store!!");    
   $self->output( \@output_clustered );  #  This is just to store the lincRNA genes so test_RunnableDB script can find them.
}



####### START KB15 ADD METHOD
sub get_genes_of_biotypes_by_db_hash_ref { 
  my ($self, $href) = @_;

  my %dbnames_2_biotypes = %$href ; 


  ### ### print  "DEBUG::get_genes_of_biotypes_by_db_hash_ref::Get genes " . scalar(keys %dbnames_2_biotypes) . "\n dumper:" . Dumper(%$href) . "\n"; 


  my @genes_to_fetch;  
  foreach my $db_hash_key ( keys %dbnames_2_biotypes )  {
    # print  "DEBUG::get_genes_of_biotypes_by_db_hash_ref::1 $db_hash_key\n";  # <--- name of the database to use

    my @biotypes_to_fetch = @{$dbnames_2_biotypes{$db_hash_key}};  
    
    # print  "----> " . Dumper(@biotypes_to_fetch) . "<---- \n";
    # foreach my $biotype  ( @biotypes_to_fetch ) { 
    #   print    "----------> $biotype  \n";              # <--- name of the biotype to use
    # } 
    
    
    my $set_db = $self->hrdb_get_dba($self->param($db_hash_key));
   # my $set_db = $self->hrdb_get_dba($self->param('source_protein_coding_db'));


    my $dna_dba = $self->hrdb_get_dba($self->param('reference_db'));
    if($dna_dba) {
      $set_db->dnadb($dna_dba);
    }

    print  "DEBUG::get_genes_of_biotypes_by_db_hash_ref::2 $db_hash_key\n";
    my $test_id = $self->param('iid');
    # my $test_id = "chromosome"; 
    my $slice = $self->fetch_sequence($test_id, $set_db, undef, undef, $db_hash_key)  ;
   # my $slice = $self->fetch_sequence($test_id, $set_db, undef, undef, 'source_protein_coding_db')  ;

    # implementation of fetch_all_biotypes ....  
    my $fetch_all_biotypes_flag ; 
    foreach my $biotype  ( @biotypes_to_fetch ) {   
      if ($biotype=~m/fetch_all_biotypes/ ) {    
        $fetch_all_biotypes_flag = 1 ; 
      }
    }  
    if ( $fetch_all_biotypes_flag ) {  
         print  "fetching ALL biotypes for slice out of db $db_hash_key :\n" ; 
         my $genes = $slice->get_all_Genes(undef,undef,1) ; 
         push @genes_to_fetch, @$genes;  
         my %tmp ; 
         for ( @$genes ) {  
           $tmp{$_->biotype}++; 
         }  
         # foreach ( keys %tmp ) {  
         #   print  "found $_ $tmp{$_}\n" ; 
         # }  
         # print  scalar(@genes_to_fetch) . " genes fetched in total\n" ; 
    } else { 
      foreach my $biotype  ( @biotypes_to_fetch ) { 
         my $genes = $slice->get_all_Genes_by_type($biotype,undef,1);
         if ( @$genes == 0 ) {
           warning("No genes of biotype $biotype found in $set_db\n");
         } 
         # if ( $self->verbose ) { 
           # print  "$db_hash_key [ " . $set_db->dbname  . " ] Retrieved ".@$genes." of type ".$biotype."\n";
           # print  "DEBUG::HiveLincRNA::get_genes_of_biotypes_by_db_hash_ref: " . $db_hash_key . " Retrieved ".@$genes." of biotype ".$biotype."\n";
         # }
         push @genes_to_fetch, @$genes;
      }  
    }  
  } 
  return \@genes_to_fetch;
}
####### END KB15 ADD METHOD




sub write_output{
  my ($self) = @_; 

  print  "\nWRITING RESULTS IN OUTPUT DB and/or VALIDATION DB... " . "\n";



  # update genes in the source db which cluster with processed_transcripts or lincRNAs
  # (if requested in the config file)

  if ( $self->param('MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB') == 0  ) { 

    my @proc_tran_genes_to_update = @{ $self->single_runnable->proc_tran_genes_to_update} ;    
    print  "  have " .scalar(@proc_tran_genes_to_update ) . " processed_transcript genes to update in VALIDATION DB\n" ;  

    # BK add start: 
    my $db = $self->hrdb_get_dba($self->param('source_cdna_db'));
    $self->hrdb_set_con($db,'source_cdna_db');
    my $ga = $self->hrdb_get_con('source_cdna_db')->get_GeneAdaptor;
    # BK add end

    # my $ga = $self->get_dbadaptor($self->UPDATE_SOURCE_DB)->get_GeneAdaptor();
    # my $db = $self->get_dbadaptor($self->UPDATE_SOURCE_DB);    
    my $update_analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($self->OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME);     
    unless ( defined $update_analysis && ref($update_analysis)=~m/Bio::EnsEMBL::Analysis/ ) {  
      throw ( " Analysis with logic_name " . $self->OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME . " can't be found in " . $db->dbname . "\@".$db->host . "\n" );
    }

    for my $ug ( @proc_tran_genes_to_update ) {      

      # before we update the analysis logic_name and gene biotype, we check if the analysis 
      # of gene processed_transcript is 'havana' or 'ensembl_havana_gene' to make sure no
      # Ensembl-only proc_trans genes slipped through the net (there shouldn't be any in
      # the first place!).

      my $hav_logic_name_to_match = $self->PROC_TRANS_HAVANA_LOGIC_NAME_STRING;
      if ( $ug->analysis->logic_name =~m/$hav_logic_name_to_match/ ) {
          $ug->analysis($update_analysis);
          $ug->biotype('proc_trans_turned_lincRNA'); 
          $ga->update($ug);  
          print  "  updated gene " . $ug->dbID . "\n";  
      }else {  
        warning("not updating gene " . $ug->biotype . " with dbID " . $ug->dbID . " as it has the wrong analysis logic_name: " 
        . $ug->analysis->logic_name . " ( to update, the logic_name should contain the string " . $self->PROC_TRANS_HAVANA_LOGIC_NAME_STRING . ")"); 
      } 
    }
  }


  if ( $self->param('MARK_EXISTING_LINCRNA_IN_VALIDATION_DB') == 0  ) {

    my @existing_lincRNA_genes_to_update = @{ $self->single_runnable->old_lincRNA_genes_to_update } ;
    print  "  have " .scalar(@existing_lincRNA_genes_to_update ) . " existing lincRNA genes to update in VALIDATION DB\n" ;
    my $db = $self->get_dbadaptor($self->UPDATE_SOURCE_DB);
    my $ga = $self->get_dbadaptor($self->UPDATE_SOURCE_DB)->get_GeneAdaptor();
    my $update_analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($self->OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME);

    OLD_G: for my $old_g ( @existing_lincRNA_genes_to_update ) {
      # No need to check sanity of analysis logic_name here because all "old" (existing) lincRNA genes would have
      # had logic_name "ensembl" or "ensembl_havana_merge". None of them should have logic_name "havana" because 
      # Havana does not annotate lincRNAs at the gene level (but only at the transcript level)
      $old_g->biotype('lincRNA_common');
      $old_g->analysis($update_analysis);
      $ga->update($old_g);
      print  "  updated gene " . $old_g->dbID . "\n";
    }
  }
 
  # my $lincrna_ga = $self->output_db->get_GeneAdaptor;

## BK add this start
  my $dba = $self->hrdb_get_dba($self->param('lincRNA_output_db'));
  $self->hrdb_set_con($dba,'lincRNA_output_db');

  my $lincrna_ga  = $self->hrdb_get_con('lincRNA_output_db')->get_GeneAdaptor;
## BK add this end





  my @genes_to_write = @{$self->output}; 

  if ( $self->WRITE_REJECTED_NCRNAS == 1 ) {   
    print  "Writing rejected genes. They have not been collapsed by GeneBuilder so are still in one-gene-one-transcript format, as generated by lincRNAFinder.\n"; 
    print  "  mutliple-biotype rejects: " . scalar(@{$self->single_runnable->ncrna_clusters_with_multiple_ens_biotypes}) . " genes \n";
    print  "  protein-domain rejects: " . scalar(@{$self->single_runnable->ncrna_clusters_with_protein_domain}) . " genes \n";
    print  "  single-biotype rejects " . scalar(@{$self->single_runnable->ncrna_clusters_with_single_ens_biotype}) . " genes \n";
   
    push ( @genes_to_write, @{$self->single_runnable->ncrna_clusters_with_single_ens_biotype} );
    push ( @genes_to_write, @{$self->single_runnable->ncrna_clusters_with_multiple_ens_biotypes} );
    push ( @genes_to_write, @{$self->single_runnable->ncrna_clusters_with_protein_domain} );
  }  

  print  "***HAVE ". scalar(@genes_to_write) ." GENE(S) TO WRITE IN TOTAL (INCLUDING VALID AND REJECTED lincRNAs).\n";  

  my $sucessful_count = 0 ; 
  foreach my $gene(@genes_to_write){  
  	print "DEBUG_BK::00gene\n";
    print_Gene_Transcript_and_Exons($gene);

    for ( @{$gene->get_all_Transcripts} ) {    
      print  "DEBUG::Evaluation::write $_ : ".$_->seq_region_start . " " . $_->biotype . " " . $_->translation()->seq() . " " . $_->translation()->start() . " gene: " . $gene->dbID . " " . $gene->biotype . " " . $_->strand . " " . $gene->strand . " " ."\n"; 
    } 
    
    # if (scalar(@t) > 1 ) { 
    #   $self->throw("more than one transcript " )   ; 
    # } 
    my $logic_name_to_be = "lincRNA_set_test_3"; 
    my @t = @{ $gene->get_all_Transcripts}; 
    my $strand_to_use  =  $gene->strand;    
    my $biotype_to_use = "empty_something_is_wrong"; 
    if ( $gene->biotype eq $self->FINAL_OUTPUT_BIOTYPE ) {
      $biotype_to_use = $gene->biotype ; 
      $logic_name_to_be = "final_biotype_set"; 
      print "### NEED TO USE ANOTHER BIOTYPE:: $biotype_to_use AND logic_name  $logic_name_to_be \n ";
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
print "howmany_transcripts:: " . scalar(@t) . "\n" ;

    print_Gene_Transcript_and_Exons($gene);
    
    foreach my $tr (@t) {
    my $exs      = []; # new array of exons
    $tr->translation(undef);	
  
    foreach my $ex_load (@{ $tr->get_all_Exons } ) {
      my $start_exon     = $ex_load->start();
      my $end_exon       = $ex_load->end();
      my $slice          = $ex_load->slice();
      my $q_strand       = $ex_load->strand; 
      # $start_exon        =~ s/-//;
      # $end_exon          =~ s/-//;
      print  'create-exonXXXX: ' . $start_exon . ' ' . $end_exon . ":: " . $q_strand  . "\n";
      
      my $ex_load        = Bio::EnsEMBL::Exon->new(                     
                                                                -start  => $start_exon,
                                                                -end    => $end_exon,
                                                                -strand => $q_strand, 
                                                                -analysis => $analysis,
                                                                -phase  => -1,
                                                                -end_phase => -1,
                                                                -slice  => $slice
                                             );
      print  'create-exon: ' . $start_exon . ' ' . $end_exon . ":: " . $q_strand  . "\n";
      push(@{$exs}, $ex_load); 
    }
      # Now LOAD transcript, gene and finally the db
      my $transcript = Bio::EnsEMBL::Transcript->new(
                                                   -exons => $exs, 
                                                   -analysis => $analysis, 
                                                   -biotype => $biotype_to_use,
                                                   -strand  => $strand_to_use,  
                                                   );
        # my $start_exon = $exs->[0];
        # my $end_exon   = $exs->[-1];
 
      # Set the phases
      # calculate_exon_phases($transcript, 0);
      print "xXXXxx:a bit \n"; 
      $gene_linc->add_Transcript($transcript);    
    }

      print_Gene_Transcript_and_Exons($gene_linc);
  
    
    $gene_linc->status(undef); 
    $gene_linc->analysis($self->analysis);   
    print "exon::01 \n";
    print_Gene_Transcript_and_Exons($gene_linc);
    empty_Gene($gene_linc, 1);
    print "exon::02 \n";
    print_Gene_Transcript_and_Exons($gene_linc);
    print "exon::03 \n";
    
    eval{
      $lincrna_ga->store($gene_linc);
    };
    if($@){
      warning("Failed to write gene ".id($gene_linc)." ".coord_string($gene_linc)." $@");
    }else{
      $sucessful_count++;
      print  "STORED LINCRNA GENE ".$gene_linc->dbID. " "  . $gene_linc->biotype  . " " .  $gene_linc->strand . " " . "  use: $strand_to_use" . "\n";
    }
  } 
  print  $sucessful_count ." genes written to FINAL OUTPUT DB \n"; # . $self->output_db->dbname . " @ ". $self->output_db->host . "\n"  ;   

  if($sucessful_count != @genes_to_write ) { 
    $self->throw("Failed to write some genes");
  }
}

sub output_db{
  my ($self, $db) = @_; 
print  "DEBUG::Evaluator need to change and set output db \n";
  if(defined $db && ref($db)=~m/Bio::EnsEMBL::DBSQL/ ){
    $self->{output_db} = $db;
  }
  if(!$self->{output_db}){
    $db = $self->get_dbadaptor($self->FINAL_OUTPUT_DB); 
    $self->{output_db} = $db; 
  }
  return $self->{output_db};
}


sub strip_phase {
# not used:: Should delete it!  
  my ($transcript_to_strip) = @_;
  my $exon_refs = $$transcript_to_strip->get_all_Exons();
  foreach my $exon (@{$exon_refs}) {
    $exon->phase(-1);
    $exon->end_phase(-1);
  }

}


sub get_gene_sets {
  my ($self) = @_;
  my @genes; 

   my %sets_to_cluster = %{$self->CLUSTERING_INPUT_GENES};  

    if ( keys %sets_to_cluster != 2 ) { 
        throw ("you should only have 2 sets to cluster against - you can't cluster against more sets \n" ); 
    } 
 
  # check if hash-key name is correct : 
  unless ( exists $sets_to_cluster{"SET_1_CDNA"} && 
            exists $sets_to_cluster{"SET_2_PROT"}) {  
   throw( " configuration error - I expect to get 2 sets of genes with names \"SET_1_CDNA\" and \"SET_2_PROT\" - check your config - you can't change these names!!!! \n" ) ; 
  }  
 
  foreach my $set ( keys %sets_to_cluster) { 
    my %this_set = %{$sets_to_cluster{$set}};      
    my @genes_in_set; 
    foreach my $database_db_name ( keys ( %this_set)) {  
       my $set_db = $self->get_dbadaptor($database_db_name);
       #$set_db->disconnect_when_inactive(1);  
       my $slice = $self->fetch_sequence($self->input_id, $set_db);  
       my @biotypes = @{$this_set{$database_db_name}}; 
       for my $biotype  ( @biotypes ) { 
          my $genes = $slice->get_all_Genes_by_type($biotype,undef,1);
          if ( @$genes == 0 ) {  
            warning("No genes of biotype $biotype found in $set_db\n"); 
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


=head2 read_and_check_config

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : hashref from config file
  Function  : call the superclass method to set all the variables and carry
  out some sanity checking
  Returntype: N/A
  Exceptions: throws if certain variables arent set properly
  Example   : 

=cut

sub read_and_check_config{
  my ($self, $hash) = @_;

  $self->SUPER::read_and_check_config($hash);
  #######
  #CHECKS
  ####### 
  foreach my $var(qw(LINCRNA_DB FINAL_OUTPUT_DB FINAL_OUTPUT_BIOTYPE VALIDATION_DBS WRITE_REJECTED_NCRNAS MARK_EXISTING_LINCRNA_IN_VALIDATION_DB)){ 
    throw("RunnableDB::lincRNAEvaluator $var config variable is not defined") if (!defined $self->$var ) ; 
  }
}




# HIVE check
sub hive_set_config {
  my $self = shift;

  # Throw is these aren't present as they should both be defined
  unless($self->param_is_defined('logic_name') && $self->param_is_defined('module')) {
    throw("You must define 'logic_name' and 'module' in the parameters hash of your analysis in the pipeline config file, ".
          "even if they are already defined in the analysis hash itself. This is because the hive will not allow the runnableDB ".
          "to read values of the analysis hash unless they are in the parameters hash. However we need to have a logic name to ".
          "write the genes to and this should also include the module name even if it isn't strictly necessary"
         );
  }

  # Make an analysis object and set it, this will allow the module to write to the output db
  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                           );
  $self->analysis($analysis);

  # Now loop through all the keys in the parameters hash and set anything that can be set
  my $config_hash = $self->param('config_settings');
  foreach my $config_key (keys(%{$config_hash})) {
    if(defined &$config_key) {
      $self->$config_key($config_hash->{$config_key});
    } else {
      throw("You have a key defined in the config_settings hash (in the analysis hash in the pipeline config) that does ".
            "not have a corresponding getter/setter subroutine. Either remove the key or add the getter/setter. Offending ".
            "key:\n".$config_key
           );
    }
  }
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


sub filter_genes_by_translation{ 
  my ($aref) = @_ ;   

  my (@g_with_transl, @no_translation ) ;  

  for my $g ( @$aref ) { 
    my @t = @{ $g->get_all_Transcripts};   
    if (scalar(@t) > 1 ) { 
      throw("more than one transcript " )   ; 
    } 
    for ( @t ) {  
       if(  $_->translation ) {  
         push @g_with_transl , $g; 
       } else {  
         push @no_translation, $g ; 
       } 
     } 
  }  
  return [\@g_with_transl,\@no_translation];
 } 


use vars '$AUTOLOAD';
sub AUTOLOAD {
  my ($self,$val) = @_;
  (my $routine_name=$AUTOLOAD)=~s/.*:://; #trim package name
  $self->{$routine_name}=$val if $val ;
  return $self->{$routine_name} ;
}
sub DESTROY {} # required due to AUTOLOAD





1;
