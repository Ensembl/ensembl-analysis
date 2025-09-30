=head1 LICENSE
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::lincRNAEvaluator - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::lincRNAEvaluator;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::lincRNAEvaluator; 
use Bio::EnsEMBL::Analysis::Runnable::lincRNAEvaluator; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange); 
use Bio::EnsEMBL::Analysis::RunnableDB; 
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild Bio::EnsEMBL::Analysis::RunnableDB);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::lincRNAEvaluator
  Function  : instatiates a lincRNAEvaluator object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::lincRNAEvaluator 
  Exceptions: 
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);  
  $self->read_and_check_config($LINCRNA_EVAL_CONFIG_BY_LOGIC);  
  return $self;
}



sub fetch_input{
  my ($self) = @_;

  # Fetch sequence/slice 

  $self->query($self->fetch_sequence); 

  # Get lincRNA candidate genes and break each of them down into single-transcript genes:

  my $lincrna_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->LINCRNA_DB);
  
  my @single_transcript_lincrna_genes = @{$self->create_single_transcript_genes($lincrna_genes)}; 
  print "Made ". scalar(@single_transcript_lincrna_genes) . " single_transcript genes, broken down from " . scalar(@$lincrna_genes) .
        " multi_transcript lincRNA candidates from lincRNAFinder stage.\n" ;  

  # Get ensembl genes as the validation set  
  
  my $validation_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->VALIDATION_DBS);  
  print "\nFetched ". scalar(@$validation_genes) . " Ensembl genes as the validation set.\n";

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
 
  print "\nCOLLAPSING SINGLE-TRANSCRIPT lincRNA 'GENES' BY GENEBUILDER TO CREATE FINAL lincRNA GENES.\n";

  #
  # First genebuilder run with unclustered lincRNAs
  #

  print "\nRunning GeneBuilder for " . scalar(@{$self->single_runnable->unclustered_ncrnas}). " unclustered lincRNAs...\n" ;  

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

   push @output_clustered, @{$gb->output()} ;  
   print "  GeneBuilder returned ". scalar(@output_clustered) . " lincRNA genes which do not overlap with any other gene ";
   print "(e.g. not proc_tran, existing lincRNAs or protein_coding genes).\n";

  # 
  # OPTIONAL: second genebuilder run with lincRNAs which cluster with processed transcript  
  # 

  if ( $self->WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS == 1 ) { 
     print "\nRunning GeneBuilder for " . scalar(@{$self->single_runnable->ncrna_clusters_with_processed_transcript}). 
      " lincRNAs which cluster with processed transcript...\n" ; 

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
 
     print "  GeneBuilder returned ". scalar( @{$gb->output()} ) . " lincRNA genes which overlapped with processed_transcript genes.\n";
     push @output_clustered, @{$gb->output()} ;   
   }

  # 
  # OPTIONAL: third genebuilder run with lincRNAs which cluster with existing lincRNAs in the core/SOURCE_PROTEIN_CODING DB  
  # 


  if ( $self->WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS == 1 ) {
     print "\nRunning GeneBuilder for " . scalar(@{$self->single_runnable->ncrna_clusters_with_existing_lincRNAs}).
      " lincRNAs which cluster with existing lincRNAs...\n" ;
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
     print "  GeneBuilder returned ". scalar( @{$gb->output()} ) . " lincRNA genes which overlapped with existing lincRNA genes.\n";
     push @output_clustered, @{$gb->output()} ;
   }

   print "\nGENEBUILDER RETURNED " .@output_clustered . " lincRNA GENES FOR WRITING (THIS DOES NOT INCLUDE REJECTED lincRNAs). THE TYPES OF lincRNA GENES WRITTEN DEPEND ON THE lincRNAEvaluator CONFIG SETTINGS.\n" ;  
   $self->genes_to_write( \@output_clustered );  # The write_output method takes lincRNA genes from $self->genes_to_write
   $self->output( \@output_clustered );  #  This is just to store the lincRNA genes so test_RunnableDB script can find them.
}


sub write_output{
  my ($self) = @_; 

  print "\nWRITING RESULTS IN OUTPUT DB and/or VALIDATION DB...\n";

  # update genes in the source db which cluster with processed_transcripts or lincRNAs
  # (if requested in the config file)

  if ( $self->MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB == 1  ) { 

    my @proc_tran_genes_to_update = @{ $self->single_runnable->proc_tran_genes_to_update} ;    
    print "  have " .scalar(@proc_tran_genes_to_update ) . " processed_transcript genes to update in VALIDATION DB\n" ;  
  
    my $ga = $self->get_dbadaptor($self->UPDATE_SOURCE_DB)->get_GeneAdaptor();
    my $db = $self->get_dbadaptor($self->UPDATE_SOURCE_DB);    
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
          print "  updated gene " . $ug->dbID . "\n";  
      }else {  
        warning("not updating gene " . $ug->biotype . " with dbID " . $ug->dbID . " as it has the wrong analysis logic_name: " 
        . $ug->analysis->logic_name . " ( to update, the logic_name should contain the string " . $self->PROC_TRANS_HAVANA_LOGIC_NAME_STRING . ")"); 
      } 
    }
  }


  if ( $self->MARK_EXISTING_LINCRNA_IN_VALIDATION_DB == 1  ) {

    my @existing_lincRNA_genes_to_update = @{ $self->single_runnable->old_lincRNA_genes_to_update } ;
    print "  have " .scalar(@existing_lincRNA_genes_to_update ) . " existing lincRNA genes to update in VALIDATION DB\n" ;
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
      print "  updated gene " . $old_g->dbID . "\n";
    }
  }
 
  my $lincrna_ga = $self->output_db->get_GeneAdaptor;

  my @genes_to_write = @{$self->genes_to_write};

  if ( $self->WRITE_REJECTED_NCRNAS == 1 ) {   
    print "Writing rejected genes. They have not been collapsed by GeneBuilder so are still in one-gene-one-transcript format, as generated by lincRNAFinder.\n"; 
    print "  mutliple-biotype rejects: " . scalar(@{$self->single_runnable->ncrna_clusters_with_multiple_ens_biotypes}) . " genes \n";
    print "  protein-domain rejects: " . scalar(@{$self->single_runnable->ncrna_clusters_with_protein_domain}) . " genes \n";
    print "  single-biotype rejects " . scalar(@{$self->single_runnable->ncrna_clusters_with_single_ens_biotype}) . " genes \n";
   
    push ( @genes_to_write, @{$self->single_runnable->ncrna_clusters_with_single_ens_biotype} );
    push ( @genes_to_write, @{$self->single_runnable->ncrna_clusters_with_multiple_ens_biotypes} );
    push ( @genes_to_write, @{$self->single_runnable->ncrna_clusters_with_protein_domain} );
  }  

  print "***HAVE ". scalar(@genes_to_write) ." GENES TO WRITE IN TOTAL (INCLUDING VALID AND REJECTED lincRNAs).\n";  

  my $sucessful_count = 0 ; 
  foreach my $gene(@genes_to_write){  
    # for ( @{$gene->get_all_Transcripts} ) {  
    #   print "$_ : ".$_->seq_region_start . " " . $_->biotype . "\n"; 
    # }
    $gene->status(undef); 
    $gene->analysis($self->analysis);   
    eval{
      $lincrna_ga->store($gene);
    };
    if($@){
      warning("Failed to write gene ".id($gene)." ".coord_string($gene)." $@");
    }else{
      $sucessful_count++;
      # print "STORED LINCRNA GENE ".$gene->dbID."\n";
    }
  } 

  print $sucessful_count ." genes written to " . $self->output_db->dbname . " @ ".
  $self->output_db->host . "\n"  ;   

  if($sucessful_count != @genes_to_write ) { 
    throw("Failed to write some genes");
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
          print "Retrieved ".@$genes." of type ".$biotype."\n";
          push @genes_in_set, @$genes; 
      }  
    } 
    my @single_transcript_genes = @{$self->create_single_transcript_genes(\@genes_in_set)};
   $self->all_gene_sets($set,\@genes_in_set); 
  }   
}


sub create_single_transcript_genes{
  my ($self, $genes) = @_; 

  print "\nCreating single_transcript genes ....\n" ;
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
  if(defined $arg){
    $self->{'FINAL_OUTPUT_DB'} = $arg;
  }
  return $self->{'FINAL_OUTPUT_DB'};
}  


sub FINAL_OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'FINAL_OUTPUT_BIOTYPE'} = $arg;
  }
  return $self->{'FINAL_OUTPUT_BIOTYPE'};
}  


sub VALIDATION_DBS{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'VALIDATION_DBS'} = $arg;
  }
  return $self->{'VALIDATION_DBS'};
} 


sub EXCLUDE_SINGLE_EXON_LINCRNAS {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'EXCLUDE_SINGLE_EXON_LINCRNAS'} = $arg;
  }
  return $self->{'EXCLUDE_SINGLE_EXON_LINCRNAS'};

}

sub MAX_FRAMESHIFT_INTRON_LEN {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MAX_FRAMESHIFT_INTRON_LEN'} = $arg;
  }
  return $self->{'MAX_FRAMESHIFT_INTRON_LEN'};

}

sub EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS'} = $arg;
  }
  return $self->{'EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS'};
}

sub MAX_TRANSCRIPTS_PER_CLUSTER {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MAX_TRANSCRIPTS_PER_CLUSTER'} = $arg;
  }
  return $self->{'MAX_TRANSCRIPTS_PER_CLUSTER'};
}

sub MARK_EXISTING_LINCRNA_IN_VALIDATION_DB{
  my ($self, $arg) = @_;
  if (defined $arg) {
    $self->{'MARK_EXISTING_LINCRNA_IN_VALIDATION_DB'} = $arg;
  }
  return $self->{'MARK_EXISTING_LINCRNA_IN_VALIDATION_DB'};
}


sub WRITE_REJECTED_NCRNAS{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'WRITE_REJECTED_NCRNAS'} = $arg;
  }
  return $self->{'WRITE_REJECTED_NCRNAS'};
}  


sub UPDATE_SOURCE_DB{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'UPDATE_SOURCE_DB'} = $arg;
  }
  return $self->{'UPDATE_SOURCE_DB'};
} 

sub WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS'} = $arg;
  }
  return $self->{'WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS'};
} 

sub WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS'} = $arg;
  }
  return $self->{'WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS'};
}


sub MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB {
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB'} = $arg;
  }
  return $self->{'MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB'};
} 


sub LINCRNA_DB{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'LINCRNA_DB'} = $arg;
  }
  return $self->{'LINCRNA_DB'};
}  

sub OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME'} = $arg;
  }
  return $self->{'OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME'};
} 


sub PROC_TRANS_HAVANA_LOGIC_NAME_STRING{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'PROC_TRANS_HAVANA_LOGIC_NAME_STRING'} = $arg;
  }
  return $self->{'PROC_TRANS_HAVANA_LOGIC_NAME_STRING'};
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
