#
#
# BioPerl module for GeneBuilder
#
# Cared for by EnsEMBL <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::HavanaAdder

=head1 SYNOPSIS

# This is the main analysis database

    my $genebuilder = new Bio::EnsEMBL::Analysis::Runnable::HavanaAdder
      (
       '-slice'   => $self->query,
       '-input_id' => $self->input_id,
      );



=head1 DESCRIPTION

This module reads your favourite annotations (ensembl, protein_coding,...)
on the one hand, and manually curated plus features on the other hand. 
The product of Bio::EnsEMBL::Analysis::Runnable::HavanaAdder is a combination of
both annotations where redundant transcripts are eliminated.
The resulting transcripts are combined into genes. For more details, follow the list of methods called
by build_Genes() method and the description in each one.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Runnable::HavanaAdder;

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster;

#use Bio::EnsEMBL::Analysis::Config::General     qw (
#							       GB_INPUTID_REGEX
#							      );

use Bio::EnsEMBL::Analysis::Config::HavanaAdder qw (
                                                    GB_ENSEMBL_INPUT_GENETYPE
                                                    GB_HAVANA_INPUT_GENETYPE
                                                    GB_ENSEMBL_PROCESSED_GENETYPE
                                                    GB_HAVANA_PROCESSED_GENETYPE
                                                    GB_ENSEMBL_PSEUDO_GENETYPE
                                                    GB_HAVANA_PSEUDO_GENETYPE
                                                    MERGED_TRANSCRIPT_OUTPUT_TYPE
                                                    HAVANA_LOGIC_NAME
                                                    MERGED_GENE_LOGIC_NAME
                                                    MERGED_TRANSCRIPT_LOGIC_NAME
                                                    HAVANA_GENE_OUTPUT_BIOTYPE
                                                    MERGED_GENE_OUTPUT_BIOTYPE
                                                    ENSEMBL_GENE_OUTPUT_BIOTYPE
                                                    );

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Root);


############################################################

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($slice,$input_id) = $self->_rearrange([qw(SLICE INPUT_ID)],
					      @args);

    $self->throw("Must input a slice to HavanaAdder") unless defined($slice);
    $self->{_final_genes} = [];
    $self->{_gene_types}  = [];

    $self->query($slice);
    $self->gene_types($GB_ENSEMBL_INPUT_GENETYPE);
    $self->gene_types($GB_HAVANA_INPUT_GENETYPE);
  
    $self->input_id($input_id);

    return $self;
}

############################################################

=head2 input_id

 Function: get/set for input id
 Returns : string
 Args    : string (it expects a string of the format chr_name.start_coord-end_coord

=cut
  
sub input_id {
  my ($self,$id) = @_;
  
  if (defined($id)) {
    $self->{_input_id} = $id;
  }
  return $self->{_input_id};
}

############################################################

=head2 build_Genes

 Example    : my @genes = $self->build_Genes
 Description: builds genes. It is like the run method in Runnables. It calls everything that needs to be done.
 Returns    : none
 Args       : none
 Caller     : Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder

=cut

sub build_Genes{
  my ($self) = @_;
  
  print STDERR "Building genes...\n";
  
  # get all genes of type defined in gene_types() on this slice
  $self->get_Genes;

  my @ensembl_coding_transcripts = $self->combined_Transcripts;
  my @ensembl_processed_transcripts = $self->combined_Processed_Transcripts;
  my @ensembl_pseudo_transcripts = $self->combined_PseudoTranscripts;

  my @havana_coding_genes = $self->havana_Coding_Genes;
  my @havana_processed_transcript_genes = $self->havana_Processed_Transcript_Genes;
  my @havana_pseudogene_Genes = $self->havana_Pseudogene_Genes;


  # do a preliminary clustering
  my @preliminary_coding_genes = $self->cluster_into_Genes(\@ensembl_coding_transcripts,\@havana_coding_genes);

  print "Coding gene clusters: ",scalar(@preliminary_coding_genes),"\n";

  my @preliminary_processed_genes = $self->cluster_into_PseudoGenes(\@ensembl_processed_transcripts,\@havana_processed_transcript_genes);

  print "Processed transcript genes clusters: ",scalar(@preliminary_processed_genes),"\n";

  my @preliminary_pseudo_genes = $self->cluster_into_PseudoGenes(\@ensembl_pseudo_transcripts,\@havana_pseudogene_Genes);

  print "Pseudogene clusters: ",scalar(@preliminary_pseudo_genes),"\n";

  push (@preliminary_pseudo_genes,@preliminary_processed_genes);

  my @clustered_gene_set = $self->combine_gene_clusters(\@preliminary_coding_genes,\@preliminary_pseudo_genes);

  print "Total clusters: ",scalar(@clustered_gene_set),"\n";

  # merge redundant ensembl transcripts which match a havana one
  $self->_merge_redundant_transcripts(\@clustered_gene_set);
  #$self->_merge_redundant_transcripts(\@preliminary_coding_genes);

  # make shared exons unique objects
  my @genes =  $self->_make_shared_exons_unique( @clustered_gene_set );
 
  print STDERR scalar(@genes)." genes built\n";
  
  $self->update_biotypes(@genes);

  $self->final_genes( @genes );
}


sub _merge_redundant_transcripts{
  my ($self, $genes) = @_;

  print "Number of genes in clusters: ", scalar(@{$genes}),"\n";

 GENE:
  foreach my $gene(@$genes){
    my @transcripts = @{$gene->get_all_Transcripts};
    my @havana;
    my @ensembl;
    
    #print "number of transcript: ",scalar(@transcripts),"\n";
    # are there any havana transcripts?
    TRANSCRIPT:foreach my $transcript(@transcripts){
      
      if($transcript->biotype =~  /_havana/){
        #print "I'm a havana transcript with biotype: ",$transcript->biotype,"\n";
        push(@havana, $transcript);
        next TRANSCRIPT;
      }
      # }
      push(@ensembl, $transcript);
    }
    if (!scalar(@havana)){
      next GENE;
    }
    if (!scalar(@ensembl)){
      next GENE;
    } 
    #print "havana tran: ",scalar(@havana),"\n";
    #print "ensembl tran: ",scalar(@ensembl),"\n";

    # compare each havana transcript to each ensembl one
    foreach my $ht(@havana){

      # We add an attribute to the havana transcripts that show which supporting features it has
      $self->add_havana_attribute($ht,$ht);

      #print "Deleting havana transcript supporting features\n";
      $ht->flush_supporting_features;

      #print "Number of ensembl transcripts: ",scalar(@ensembl),"\n";
      my $delete_trans = 0;
      my @t_pair;

      foreach my $et(@ensembl){
        
        my $delete_t = $self->are_matched_pair($ht, $et);
        
        # We check all posible match pairs and give preference to the one that shares CDS and UTR
        # This was added so only the best matching havana/ensembl pair is chossen and to avoid a one to many link
        if ($delete_t){
          if ($delete_t == $et){
            #print "I'm here\n";
            $delete_trans = $delete_t;
            @t_pair = ($ht, $et);
          }elsif($delete_trans != $et && $delete_t == $ht){
            $delete_trans = $delete_t;
            @t_pair = ($ht, $et);
          }elsif($delete_trans == 0){
            $delete_trans = $delete_t;
            @t_pair = ($ht, $et);
          }
        }
      }

      if($delete_trans && $delete_trans != 0){
        my $new_bt_0; #biotype
        my $new_bt_1; #biotype
        $self->set_transcript_relation($delete_trans, @t_pair);
        unless ($t_pair[0]->biotype =~/$MERGED_TRANSCRIPT_OUTPUT_TYPE/){
          $new_bt_0 = $t_pair[0]->biotype.$MERGED_TRANSCRIPT_OUTPUT_TYPE;
        }
        unless ($t_pair[1]->biotype =~/$MERGED_TRANSCRIPT_OUTPUT_TYPE/){
          $new_bt_1 = $t_pair[1]->biotype.$MERGED_TRANSCRIPT_OUTPUT_TYPE;
        }
        $t_pair[0]->biotype($new_bt_0);
        $t_pair[1]->biotype($new_bt_1);
        
        # We want to remove the redundant transcript unless both share CDS but have different UTR
        # structure as in that case we annotate both transcripts
        $self->_remove_transcript_from_gene($gene, $delete_trans) unless $delete_trans == 1;         
        
      }else{
        $self->add_ottt_xref($ht);

      }
    }
  }
}

# are_matched_pair check return 4 different possible values:
# return 0 means keep both transcripts as they have different coding region
#          or different exon structure
# return 1 means keep both as they have same coding but different UTR exon structure
# return ($ensembl) means keep havana transcript and remove ensembl 
# return ($havana) means keep ensembl transcript and remove hanana

sub are_matched_pair {
  my($self, $havana, $ensembl) = @_;


  # Fetch all exons in each transcript 
  my @hexons = @{$havana->get_all_Exons};
  my @eexons = @{$ensembl->get_all_Exons};

  my @thexons = @{$havana->get_all_translateable_Exons};
  my @teexons = @{$ensembl->get_all_translateable_Exons};
  
  # Check that the number of exons is the same in both transcripts
  return 0 unless scalar(@hexons) == scalar(@eexons);

  #print "____________________________________\n";
  #print "HAVANA ID: ",$havana->dbID, " ENSEMBL: ",$ensembl->dbID,"\n";
  #print "HAVANA ID: ",$havana->stable_id, " ENSEMBL: ",$ensembl->stable_id,"\n";

  # double check translation coords
  #print "HAVANA TRANS START: ",$havana->translation->genomic_start," END: ",$havana->translation->genomic_end,"\n";
  #print "ENSEMBL TRANS START: ",$ensembl->translation->genomic_start," END: ",$ensembl->translation->genomic_end,"\n";

  my $non_coding_e = 0;
  my $non_coding_h = 0;

  # Check if the transcript is non coding
  if(!@teexons || @teexons == 0){
    $non_coding_e =1 ;
  }
  if(!@thexons || @thexons == 0){
    $non_coding_h =1 ;
  }

  if ($non_coding_h == 1 && $non_coding_e == 1){
    # We check two non coding transcripts. If they have the same structure we keep the one from havana 
    # but if the have same exon strucutre but one is slightly longer we take the longest one of the pair.
    #print "CASE TYPE ONE\n";
    unless ($self->check_internal_exon_structure(\@eexons,\@hexons)){
      #print "BOTH selected\n";
      # CASE 0: the two transcripts have different internal exon structure
      return 0;
    }
    # CASE 1: Havana is longer or both are exactly the same
    if ($self->check_terminal_exon_structure(\@hexons,\@eexons)){
      #print "We keep HAvana\n";
      return $ensembl;
    }else{
      # CASE 2: EnsEMBL is longer than Havana
      #print "We keep Ensembl\n";
      return $havana;
    }


  }elsif($non_coding_h != $non_coding_e){
    # this is a case of a pseudogene overlapping a coding gene so by now we keep both
    print STDERR "Warning Pseudogene and coding overlap for HAVANA ID: ",$havana->dbID, " ENSEMBL: ",$ensembl->dbID,"\n";
    
    # If Havana is coding and ensembl is non coding
    if ($non_coding_h == 0){
      #print "CASE TYPE TWO\n";
      #return 0 unless (scalar(@eexons) == scalar(@thexons));
      unless ((scalar(@eexons) == scalar(@thexons)) && $self->check_internal_exon_structure(\@eexons,\@thexons)){
        # CASE 0: the two transcripts have different internal exon structure
        return 0;
	#Here we may consider removing the Ensembl transcript and keep only the havana coding.
      }
      # CASE 1: If the internal structure is the same we then keep the Havana one.
      return $ensembl;
    }

    # If ensembl is coding and havana is non coding
    if ($non_coding_e == 0){
      # Keep both if they don't have the same coding exon structure 
      return 0 unless scalar(@thexons) == scalar(@teexons);

      # May want to remove the translation and keep ensembl
      # First check if the internal structure of the whole transcripts is conserved
      if(scalar(@hexons) > 1){
        #print "CASE TYPE THREE\n";
        unless ($self->check_internal_exon_structure(\@hexons,\@eexons)){
          # CASE 0: the two transcripts have different internal exon structure we keep both
          return 0;       
        }
      }
      # Now we check if the coding bit is longer than the whole non coding transcript
      # CASE 1: The ensembl transcript is longer or equal so we remove the havana transcript 
      if ($self->check_terminal_exon_structure(\@hexons,\@teexons)){
        #print "CASE CONTROL \n";
        # PRINT SOMETHING USEFULL HERE TO HELP REVIEW
        $ensembl->{translation} = undef;
        $ensembl->biotype($havana->biotype."_e");
        return $havana;
      }else{
      # CASE 2: The havana transcripts is longer in both ends so we remove the ensembl transcript
        #print "CASE CONTROL 2: ",$ensembl,"\n";
        #$ensembl->{_translation_array} = [];
        #$ensembl->biotype($havana->biotype."_e");
        return $ensembl;
      }
    }
    
  }elsif($non_coding_h == 0 && $non_coding_e == 0){
    
    return 0 unless((scalar(@teexons) == scalar(@thexons)) && 
                    $havana->translation->genomic_start == $ensembl->translation->genomic_start &&
                    $havana->translation->genomic_end   == $ensembl->translation->genomic_end);
    
    # special case for single exon genes
    if(scalar(@hexons) == 1){
      #print "SINGLE EXONS!\n";
      
      if ($hexons[0]->start     == $eexons[0]->start &&
          $hexons[0]->end       == $eexons[0]->end &&
          $hexons[0]->strand    == $eexons[0]->strand &&
          $thexons[0]->coding_region_start($havana) == $teexons[0]->coding_region_start($ensembl) &&
          $thexons[0]->coding_region_end($havana) == $teexons[0]->coding_region_end($ensembl)
          ){
        # Both are exactly the same so we delete the Ensembl one unless the Ensembl one is already a merged one <-???? IS not doing this here
        return $ensembl;
        
      }elsif($hexons[0]->start     <= $eexons[0]->start &&
             $hexons[0]->end       >= $eexons[0]->end &&
             $hexons[0]->strand    == $eexons[0]->strand &&
             $eexons[0]->start     == $teexons[0]->coding_region_start($ensembl) &&
             $eexons[0]->end       == $teexons[0]->coding_region_end($ensembl)
             ){
        # Ensembl gene don't have UTR and Havana has then delete Ensembl one
        return $ensembl;
        
      }elsif((($hexons[0]->start    != $eexons[0]->start ||
               $hexons[0]->end       != $eexons[0]->end) &&
              $hexons[0]->strand    == $eexons[0]->strand) &&
             ($eexons[0]->start    != $teexons[0]->coding_region_start($ensembl) ||
              $eexons[0]->end       != $teexons[0]->coding_region_end($ensembl))
             ){
        # Both contain UTR keep ENSEMBL
        return $havana;
        
      }else{
        # We can be here when genes have different UTR start/end and different CDS start/end
        # or when the UTR start/end is the same but the CDS start/end is different
        
        #print "Keep both single exon genes\n";
        return 0;
        
      }
    }
    # if is a multi exons transcript
    else{
      # First we check the internal  coding structure of the transcript where everything has to be exactly equal
      #print "CHECKING INTERNAL EXONS \n";
      for(my $i=1; $i<=($#thexons-1); $i++){
        return 0 unless ($thexons[$i]->start     == $teexons[$i]->start &&
                         $thexons[$i]->end       == $teexons[$i]->end &&
                         $thexons[$i]->strand    == $teexons[$i]->strand 
                         );
      }
      #print "INTERNAL CODING EXONS ARE OK \n";
      
      #  now check the rest of the internal exons that are not coding. This is to find if the UTR exon structure is the same
      for(my $i=1; $i<=($#hexons-1); $i++){
        if($hexons[$i]->start     != $eexons[$i]->start ||
           $hexons[$i]->end       != $eexons[$i]->end ||
           $hexons[$i]->strand    != $eexons[$i]->strand 
                         ){
          print "CASE WITH SAME CDS BUT DIFFERENT UTR STRUCTURE\n";
          print "HAVANA DIFF UTR BOUNDARIES: ".$havana->seq_region_name." - ".$havana->seq_region_start. " - ".$havana->seq_region_end."\n";
          print "ENSEMBL DIFF UTR BOUNDARIES: ".$ensembl->seq_region_name." - ".$ensembl->seq_region_start. " - ".$ensembl->seq_region_end."\n";    
          return 1;
        }
      }
      #print "INTERNAL UTR EXONS ARE OK \n";
      
      # Then check if the first an last exon are the same in both transcripts. If just start and end of UTR are different keep havana one
      # CASE 1: Both coding and UTR are the same, keep Havana and delete Ensembl
      if ($hexons[0]->start     == $eexons[0]->start &&
          $hexons[0]->end       == $eexons[0]->end &&
          $hexons[0]->strand    == $eexons[0]->strand &&
          $hexons[-1]->start    == $eexons[-1]->start &&
          $hexons[-1]->end      == $eexons[-1]->end &&
          $hexons[-1]->strand   == $eexons[-1]->strand 
          ){
        #print "MULTIEXON DELETE ENSEMBL\n";
        return $ensembl;
        
      }elsif (#CASE 2": HAVANA HAS UTR AND ENSEMBL DOESNT, KEEP HAVANA. Forward strand
              $hexons[0]->strand == 1 &&
              $hexons[0]->end       == $eexons[0]->end &&
              $hexons[0]->strand    == $eexons[0]->strand &&
              $hexons[-1]->start    == $eexons[-1]->start &&
              $hexons[-1]->strand   == $eexons[-1]->strand &&
              $eexons[0]->start     == $teexons[0]->coding_region_start($ensembl) &&
              $eexons[-1]->end      == $teexons[-1]->coding_region_end($ensembl) &&
              ($hexons[-1]->end     != $eexons[-1]->end ||
               $hexons[0]->start    != $eexons[0]->start)  
              ){
        #print "MULTIEXON DELETE ENSEMBL\n";      
        return $ensembl;
        
      }elsif (# CASE 3: BOTH ENSEMBL AND HAVANA HAVE UTR BUT WITH DIFFERENT START/END, KEEP Havana. Forward strand
              $hexons[0]->strand == 1 &&
              $hexons[0]->end       == $eexons[0]->end &&
              $hexons[0]->strand    == $eexons[0]->strand &&
              $hexons[-1]->start    == $eexons[-1]->start &&
              $hexons[-1]->strand   == $eexons[-1]->strand &&
              ($eexons[0]->start    != $teexons[0]->coding_region_start($ensembl) ||
               $eexons[-1]->end      != $teexons[-1]->coding_region_end($ensembl)) &&
              ($hexons[-1]->end     != $eexons[-1]->end ||
               $hexons[0]->start     != $eexons[0]->start)
              ){
        #print "MULTIEXON DELETE ENSEMBL\n";
        print "CASE WITH DIFFERENT TERMINAL EXON BOUNDARIES\n";
        print "HAVANA BOUNDARIES: ".$havana->seq_region_name." - ".$havana->seq_region_start. " - ".$havana->seq_region_end."\n";
        print "ENSEMBL BOUNDARIES: ".$ensembl->seq_region_name." - ".$ensembl->seq_region_start. " - ".$ensembl->seq_region_end."\n";   
        return $ensembl;
        
      }elsif (# CASE 4: Same as case 2 but in reverse strand
              $hexons[0]->strand == -1 &&
              $hexons[0]->start     == $eexons[0]->start &&
              $hexons[0]->strand    == $eexons[0]->strand &&
              $hexons[-1]->end      == $eexons[-1]->end &&
              $hexons[-1]->strand   == $eexons[-1]->strand &&
              $eexons[-1]->start    == $teexons[-1]->coding_region_start($ensembl) &&
              $eexons[0]->end       == $teexons[0]->coding_region_end($ensembl) &&
              ($hexons[0]->end      != $eexons[0]->end ||
               $hexons[-1]->start   != $eexons[-1]->start)
              
              ){
        #print "MULTIEXON DELETE ENSEMBL\n";      
        return $ensembl;
        
      }elsif (# CASE 5: Same as case 3 but in reverse strand
              $hexons[0]->strand == -1 &&
              $hexons[0]->start     == $eexons[0]->start &&
              $hexons[0]->strand    == $eexons[0]->strand &&
              $hexons[-1]->end      == $eexons[-1]->end &&
              $hexons[-1]->strand   == $eexons[-1]->strand &&
              ($eexons[-1]->start   != $teexons[-1]->coding_region_start($ensembl) ||
               $eexons[0]->end       != $teexons[0]->coding_region_end($ensembl)) &&
              ($hexons[0]->end      != $eexons[0]->end ||
               $hexons[-1]->start   != $eexons[-1]->start)
              
              ){
        #print "MULTIEXON DELETE HAVANA\n";  
        print "CASE WITH DIFFERENT TERMINAL EXON BOUNDARIES\n";
        print "HAVANA BOUNDARIES: ".$havana->seq_region_name." - ".$havana->seq_region_start. " - ".$havana->seq_region_end."\n";
        print "ENSEMBL BOUNDARIES: ".$ensembl->seq_region_name." - ".$ensembl->seq_region_start. " - ".$ensembl->seq_region_end."\n";    
        return $ensembl;
        
      }else{
        #print "Should I be here?\n";
        #print "Keep MULTIEXON BOTH\n";
        return 1;
        
      }
      
    }
    
    print " WEIRD CASE WE DID NOT THINK ABOUT, CHECK RULES!\n";
    return 0;
  }
}

=head2 check_internal_exon_structure

  Description: check if the start and end of internal exon pairs in two sets of exons is the same
  Return: Returns 0 if they are the same, returns 1 if they are differents

=cut

sub check_internal_exon_structure {

  my ($self, $firstexons, $secondexons) = @_;

  my @exons1 = @{$firstexons};
  my @exons2 = @{$secondexons};

  if (scalar(@exons1) != scalar(@exons2)){
    print "NOOOOOOOOO WHAT HAPPENED HERE?\n";
    print scalar(@exons1), "  vs. ",scalar(@exons2),"\n";
  }

  # We check if the transcript has more than two exon as otherwise we will be checking 
  #only the coordinates in the last/second exon which may produce wrong results 
  if (scalar(@exons1) > 2){
    for(my $i=1; $i<=($#exons1-1); $i++){
      return 0 unless ($exons1[$i]->start   == $exons2[$i]->start &&
                       $exons1[$i]->end     == $exons2[$i]->end &&
                       $exons1[$i]->strand  == $exons2[$i]->strand &&
                       (($exons1[0]->strand == 1 && 
                         $exons1[0]->end     == $exons2[0]->end &&
                         $exons1[-1]->start  == $exons2[-1]->start) ||
                        ($exons1[0]->strand == -1 &&
                         $exons1[0]->start  == $exons2[0]->start &&
                         $exons1[-1]->end   == $exons2[-1]->end))
                       );
    }
    return 1;
  }else{
    return 0 unless (($exons1[0]->strand == 1 && 
                      $exons1[0]->end     == $exons2[0]->end &&
                      $exons1[-1]->start  == $exons2[-1]->start) ||
                     ($exons1[0]->strand == -1 &&
                      $exons1[0]->start  == $exons2[0]->start &&
                      $exons1[-1]->end   == $exons2[-1]->end));
  }
  return 1;
}
  
=head2 check_terminal_exon_structure
  Description: check if begining of transcript and end of transcript coincide by looking at the terminal exon coords.
  Return: return 0 is the first exon set is shorter that the second in both start and end otherwise returns 1.
=cut

sub check_terminal_exon_structure {

  my ($self, $firstexons, $secondexons ) = @_;

  my @exons1 = @{$firstexons};
  my @exons2 = @{$secondexons};

  # I added the following or check "|| $exons1[0] eq $exons1[-1]"
  # to handle single exon genes more efficiently.

  if(($exons1[0]->strand == 1 || $exons1[0] eq $exons1[-1] ) &&
     $exons1[0]->strand == $exons2[0]->strand &&
     $exons1[0]->start  <= $exons2[0]->start &&
     $exons1[-1]->end   >= $exons2[-1]->end){ 
    return 1;
  }elsif($exons1[0]->strand == -1 &&
         $exons1[0]->end    <= $exons2[0]->end &&
         $exons1[-1]->start >= $exons2[-1]->start &&
         $exons1[0]->strand == $exons2[0]->strand){
    return 1;
  }# CASE 2: EnsEMBL is longer than Havana
   elsif(($exons1[0]->strand == 1 || $exons1[0] eq $exons1[-1] )&&
         $exons1[0]->start  >= $exons2[0]->start &&
         $exons1[-1]->end   <= $exons2[-1]->end &&
         $exons1[0]->strand == $exons2[0]->strand){
     return 0;
  }elsif($exons1[0]->strand == -1 &&
         $exons1[0]->end    >= $exons2[0]->end &&
         $exons1[-1]->start <= $exons2[-1]->start &&
         $exons1[0]->strand == $exons2[0]->strand){
    return 0;
  }else{
    return 1;
  }
}

sub add_ottt_xref{
 my($self, $ht) = @_;
 #TEST
 foreach my $entry(@{ $ht->get_all_DBEntries}){
   if ($entry->dbname eq 'Vega_transcript'){
     if($entry->primary_id eq $entry->display_id){
       
       #print "I am adding an OTTT xref to the transcript\n";
       #print "OTTT TO ADD: ",$entry->primary_id,"\n";
       my $xref_ottt = new Bio::EnsEMBL::DBEntry
           (
            -primary_id =>$entry->primary_id,
            -display_id =>$ht->display_id,
            -priority => 1,
            -xref_priority => 0,
            -version => 1,
            -release => 1,
            -dbname => 'OTTT'
            );
       
       $xref_ottt->status("XREF");
       
       $ht->add_DBEntry($xref_ottt);
       #END of TEST
       
     }
   }
 }
}

sub add_ottg_xref{
 my($self, $hg, $ottg) = @_;
 
  print "I am adding an OTTG xref to the transcript\n";
 print "OTTG TO ADD: ",$ottg,"\n";#$entry->primary_id,"\n";
 my $xref_ottg = new Bio::EnsEMBL::DBEntry
     (
      -primary_id =>$ottg,#$entry->primary_id,
      -display_id =>$ottg,#$hg->display_id,
      -priority => 1,
      -xref_priority => 0,
      -version => 1,
      -release => 1,
      -dbname => 'OTTG'
      );
 
 $xref_ottg->status("XREF");
 
 $hg->add_DBEntry($xref_ottg);
 
}


sub set_transcript_relation {
  # $t_pair[0] is the havana transcript and $t_pair[1] is the ensembl transcript
  my($self, $delete_t, @t_pair) = @_;
     
  # If both share CDS and UTR is different in structure and number of exons we still keep both, and we link them via Xref
  if ($delete_t == 1){
    
    # transfer OTTT ID and/or ENST
    foreach my $entry(@{ $t_pair[0]->get_all_DBEntries}){
      if ($entry->dbname eq 'Vega_transcript'){
        if($entry->primary_id eq $entry->display_id){
          
          my $newentry = new Bio::EnsEMBL::DBEntry
              (
               -primary_id => $entry->primary_id,
               -display_id => $entry->display_id,
               -priority => 1,
               -xref_priority => 0,
               -version => 1,
               -release => 1,
               -dbname => 'shares_CDS_with_OTTT'
               );
          
          $newentry->status("XREF");
          
          $t_pair[1]->add_DBEntry($newentry);

          #TEST
          my $xref_ottt = new Bio::EnsEMBL::DBEntry
              (
               -primary_id =>$entry->primary_id,
               -display_id =>$entry->display_id,
               -priority => 1,
               -xref_priority => 0,
               -version => 1,
               -release => 1,
               -dbname => 'OTTT'
               );
          
          #print "OTTT xref to be added here\n";
          
          $xref_ottt->status("XREF");
          
          $t_pair[0]->add_DBEntry($xref_ottt);
          #END of TEST
          
        }
      }
    }

    my $link_attrib = Bio::EnsEMBL::Attribute->new
        (-CODE => 'enst_link',
         -NAME => 'enst link',
         -DESCRIPTION => 'Code to link a OTTT with an ENST when they both share the CDS of ENST',
         -VALUE => $t_pair[1]->dbID);
    
    $t_pair[1]->add_Attributes($link_attrib);

    my $xref_entry = new Bio::EnsEMBL::DBEntry
        (
         -primary_id =>$t_pair[1]->dbID,
         -display_id =>$t_pair[1]->dbID,
         -priority => 1,
         -xref_priority => 0,
         -version => 1,
         -release => 1,
         -dbname => 'shares_CDS_with_ENST'
         );

    $xref_entry->status("XREF");
    
    $t_pair[0]->add_DBEntry($xref_entry);

    #print "OTTT TO ADD: ",$t_pair[0]->stable_id,"\n";

  }
  
  # If transcript to delete is Havana we create an xref for the entry say that the transcript is CDS equal to ENSEMBL
  elsif ($delete_t == $t_pair[0]){
    # transfer OTTT ID and/or ENST
    foreach my $entry(@{ $t_pair[0]->get_all_DBEntries}){
      if ($entry->dbname eq 'Vega_transcript'){
        if($entry->primary_id eq $entry->display_id){
          my $newentry = new Bio::EnsEMBL::DBEntry
              (
               -primary_id => $entry->primary_id,
               -display_id => $entry->display_id,
               -priority => 1,
               -xref_priority => 0,
               -version => 1,
               -release => 1,
               -dbname => 'shares_CDS_with_OTTT'
               );
          
          $newentry->status("XREF");
          
          $t_pair[1]->add_DBEntry($newentry);
        }
      }
    }
    
    # We add a transcript attribute to the ensembl transcript with the start and end coords of the Havana transcript that we will delete
    my $attrib_value = $t_pair[0]->slice->coord_system_name.":".$t_pair[0]->slice->coord_system->version.":".$t_pair[0]->slice->seq_region_name.":".
        $t_pair[0]->start.":".$t_pair[0]->end.":1";
    # print "ATTRIB VALUE:---------- ",$attrib_value,"\n";
    my $attribute = Bio::EnsEMBL::Attribute->new
        (-CODE => 'TranscriptEdge',
         -NAME => 'Transcript Edge',
         -DESCRIPTION => '',
         -VALUE => $attrib_value);
    
    $t_pair[1]->add_Attributes($attribute);
    
    # When we delete a Havana transcript we want to transfer the exon supporting features to the transcript we keep    
    my @delete_e = @{$delete_t->get_all_Exons};
    my @exons    = @{$t_pair[1]->get_all_Exons};
    
    if (scalar(@delete_e) == scalar(@exons)){

        my $e;
        for ($e = 0, $e<scalar(@delete_e), $e++){
          if($delete_e[$e] && $exons[$e]){
            $self->transfer_supporting_evidence($delete_e[$e], $exons[$e]); 
          }
        }
    }
    
    # We add attributes to the havana transcript showing which supporting features where used for the transcript in Havana
    $self->add_havana_attribute($t_pair[0],$t_pair[1]);
    

  }elsif ($delete_t == $t_pair[1]){
    # If the transcript to delete is ENSEMBL we add an xref entry say that both transcripts are exact matches (including UTR)
    foreach my $entry(@{ $t_pair[0]->get_all_DBEntries}){
      if ($entry->dbname eq 'Vega_transcript'){
        if($entry->primary_id eq $entry->display_id){

          my $enstentry = new Bio::EnsEMBL::DBEntry
              ( 
                -primary_id => $entry->primary_id,
                -display_id => $entry->display_id,
                -version => 1,
                -release => 1,
                -priority => 1,
                -xref_priority => 0,
                -dbname => 'shares_CDS_and_UTR_with_OTTT'
                );
          
          $enstentry->status("XREF");
          
          $t_pair[0]->add_DBEntry($enstentry);
        }
      }
    } 
    # Transfer the supporting features both for transcript and exon of the transcript to delete to the transcript we keep
    $self->transfer_supporting_features($delete_t,$t_pair[0]);
   
    # We add attributes to the havana transcript showing which supporting features where used for the transcript in Havana
    #$self->add_havana_attribute($t_pair[0],$t_pair[0]);

  }
}

sub add_havana_attribute{
  my ($self, $transcript, $trans_to_add_attrib) = @_;

  my %evidence;
  my %t_evidence;

  foreach my $tsf (@{$transcript->get_all_supporting_features}){
    $t_evidence{$tsf->hseqname} = 1;
  }

  foreach my $te_key (keys %t_evidence){
    #print "Adding special attrib\n";
    if($te_key->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
      my $attribute = Bio::EnsEMBL::Attribute->new
          (-CODE => 'tp_otter_support',
           -NAME => 'tp otter support',
           -DESCRIPTION => 'Evidence ID that was used as protein transcript supporting feature for building a gene in Vega',
           -VALUE => $te_key);
      
      $trans_to_add_attrib->add_Attributes($attribute);
      
    }
    if($te_key->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){
      my $attribute = Bio::EnsEMBL::Attribute->new
          (-CODE => 'td_otter_support',
           -NAME => 'td otter support',
           -DESCRIPTION => 'Evidence ID that was used as cdna transcript supporting feature for building a gene in Vega',
           -VALUE => $te_key);
      
      $trans_to_add_attrib->add_Attributes($attribute);
      
    }
  }
  
  foreach my $exon (@{$transcript->get_all_Exons}){ 
    foreach my $sf (@{$exon->get_all_supporting_features}){
      $evidence{$sf->hseqname} = 1;
    }
  }

  foreach my $ev_key (keys %evidence){
    #print "Adding special attrib\n";
    if($ev_key->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
      my $attribute = Bio::EnsEMBL::Attribute->new
          (-CODE => 'ep_otter_support',
           -NAME => 'ep otter support',
           -DESCRIPTION => 'Evidence ID that was used as protein exon supporting feature for building a gene in Vega',
           -VALUE => $ev_key);
      
      $trans_to_add_attrib->add_Attributes($attribute);
    }
    if($ev_key->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
      my $attribute = Bio::EnsEMBL::Attribute->new
          (-CODE => 'ed_otter_support',
           -NAME => 'ed otter support',
           -DESCRIPTION => 'Evidence ID that was used as cdna exon supporting feature for building a gene in Vega',
           -VALUE => $ev_key);
      
      $trans_to_add_attrib->add_Attributes($attribute);    
    } 
  } 
}
  
sub transfer_supporting_features{
  my ($self, $delete_t, $transcript) = @_;
  
  #print "TRANSCRIPT IS :  ", $transcript,"\n";
  
  my @exon_features;

  # Delete all the supporting features for the Havana Transcript 
  #$transcript->flush_supporting_features;
  
  my @delete_tsf = @{ $delete_t->get_all_supporting_features };
  #my @transcript_sf = @{ $transcript->get_all_supporting_features };
  
  # print "NUMBER OF TRANSCRIPT SF: ",scalar(@transcript_sf),"\n";
  # print " AND DELETE TSF: ", scalar(@delete_tsf),"\n";
 DTSF: foreach my $dtsf (@delete_tsf){
   next DTSF unless $dtsf->isa("Bio::EnsEMBL::FeaturePair");

   $transcript->add_supporting_features($dtsf);
 }
  
  my @delete_e = @{$delete_t->get_all_Exons};
  my @exons    = @{$transcript->get_all_Exons};
  
  if (scalar(@delete_e) == scalar(@exons)){

      my $e;
      for ($e = 0, $e<scalar(@delete_e), $e++){
        if($delete_e[$e] && $exons[$e]){
          $self->transfer_supporting_evidence($delete_e[$e], $exons[$e]); 
        }
      }
  }
  
#  print "NUMBER AFT ADDITTION: ", scalar(@{ $transcript->get_all_supporting_features }),"\n";
}

sub _remove_transcript_from_gene {
  my ($self, $gene, $trans_to_del)  = @_;
     
#  print "HEY I AAAAAMMMMMMMMMMMMMMMM HEEEEEEEEEEEEEEERRRRRRRRREEEEEEEEEEEEEEE\n";
#  foreach my $entry(@{ $gene->get_all_DBEntries}){
#   print "ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
#  }

  my @newtrans;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans != $trans_to_del) {
      push @newtrans,$trans;
    }
  }

# The naughty bit!
  $gene->{_transcript_array} = [];

  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }

  return scalar(@newtrans);
}


############################################################

sub _make_shared_exons_unique{
  my ( $self, @genes ) = @_;
  my @pruned_genes;
  foreach my $gene ( @genes ){
    #TEST
 #   print "HEY I AAAAAMMMMMMMMMMMMMMMM HEEEEEEEEEEEEEEERRRRRRRRREEEEEEEEEEEEEEE\n";
 #   foreach my $entry(@{ $gene->get_all_DBEntries}){
 #    print "ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
 #   }

    # make different exon objects that are shared between transcripts 
    # ( regarding attributes: start, end, etc )
    # into unique exon objects 
    my $new_gene = $self->prune_Exons($gene);
#    foreach my $entry(@{ $new_gene->get_all_DBEntries}){
#      print "ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
#    }

    push( @pruned_genes, $new_gene );
  }
#  exit;
  return @pruned_genes;
}

############################################################


=head2 get_Genes

 Description: retrieves ensembl and havana gene annotations with supporting evidence. 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub get_Genes {
  my ($self) = @_;
  my @transcripts;
  my @processedtranscripts;
  my @pseudotranscripts;
  my @genes;
  my @processedgenes;
  my @pseudogenes;
  my @hgenes;
  my @hprocessedgenes;
  my @hpseudogenes;
  my $ensemblslice = $self->fetch_sequence($self->input_id, $self->ensembl_db);
  my $havanaslice = $self->fetch_sequence($self->input_id, $self->havana_db);

  # Fetch Ensembl genes
  print STDERR "Fetching ensembl genes\n";  
   
  foreach my $ebiotype (@{$GB_ENSEMBL_INPUT_GENETYPE}){
    EGENE:
    foreach my $egene (@{$ensemblslice->get_all_Genes_by_type($ebiotype)}){
      # Don't add those genes that contain only transcripts imported from HAVANA (this is important during a merge update)
      if ($egene->analysis->logic_name() eq $HAVANA_LOGIC_NAME){
        next EGENE;
      }else{
        push (@genes,$egene);
      } 
    }
  }

 # Fetch Ensembl Processed transcripts
  foreach my $eprocessedbt (@{$GB_ENSEMBL_PROCESSED_GENETYPE}){
 PROCESSED:
    foreach my $eprocessedgene (@{$ensemblslice->get_all_Genes_by_type($eprocessedbt)}){
      # Don't add those genes that contain only transcripts imported from HAVANA (this is important during a merge update)
      if ($eprocessedgene->analysis->logic_name() eq $HAVANA_LOGIC_NAME){
        next PROCESSED;
      }else{
        push (@processedgenes,$eprocessedgene);
      } 
    }
  }

 # Fetch Ensembl pseudogenes
  foreach my $epseudobt (@{$GB_ENSEMBL_PSEUDO_GENETYPE}){
 EPSEUDOGENE:
    foreach my $epseudogene (@{$ensemblslice->get_all_Genes_by_type($epseudobt)}){
      # Don't add those genes that contain only transcripts imported from HAVANA (this is important during a merge update)
      if ($epseudogene->analysis->logic_name() eq $HAVANA_LOGIC_NAME){
        next EPSEUDOGENE;
      }else{
        push (@pseudogenes,$epseudogene);
      } 
    }
  }
 
  print STDERR "Retrieved ".scalar(@genes)." genes of types: ".join(", ",@{$GB_ENSEMBL_INPUT_GENETYPE})."\n";
  print STDERR "Retrieved ".scalar(@processedgenes)." 'processed transcript' genes of types: ".join(", ",@{$GB_ENSEMBL_PROCESSED_GENETYPE})."\n";
  print STDERR "Retrieved ".scalar(@pseudogenes)." pseudogenes of types: ".join(", ",@{$GB_ENSEMBL_PSEUDO_GENETYPE})."\n";

  #Fetch Havana genes

  print STDERR "Fetching havana genes\n";  
  foreach my $hbiotype (@{$GB_HAVANA_INPUT_GENETYPE}){
    foreach my $hgene (@{$havanaslice->get_all_Genes_by_type($hbiotype)}){
  # We change the biotype of the havana genes/transcripts as it could happend to be the same as the ensembl ones
      my $biotype = $hgene->biotype."_havana";
      $hgene->biotype($biotype);
      foreach my $htran (@{$hgene->get_all_Transcripts}) {
        my $tbiotype = $htran->biotype."_havana";
        $htran->biotype($tbiotype);
      }
      push (@hgenes, $hgene);
    }
  }

  print STDERR "Fetching havana 'processed transcript' genes\n";  
  foreach my $hprocessedbiotype (@{$GB_HAVANA_PROCESSED_GENETYPE}){
    foreach my $hprocessedgene (@{$havanaslice->get_all_Genes_by_type($hprocessedbiotype)}){
  # We change the biotype of the havana genes/transcripts as it could happend to be the same as the ensembl ones
      my $processedbiotype = $hprocessedgene->biotype."_havana";
      $hprocessedgene->biotype($processedbiotype);
      foreach my $hprocessedtran (@{$hprocessedgene->get_all_Transcripts}) {
        my $tprocessedbiotype = $hprocessedtran->biotype."_havana";
        $hprocessedtran->biotype($tprocessedbiotype);
      }
      push (@hprocessedgenes, $hprocessedgene);
    }
  }


  #Fetch Havana pseudogenes
  print STDERR "Fetching havana pseudogenes\n";  
  foreach my $hpseudobt (@{$GB_HAVANA_PSEUDO_GENETYPE}){
    foreach my $hpseudogene (@{$havanaslice->get_all_Genes_by_type($hpseudobt)}){
  # We change the biotype of the havana genes/transcripts as it could happend to be the same as the ensembl ones
      my $biotype = $hpseudogene->biotype."_havana";
      $hpseudogene->biotype($biotype);
      foreach my $htran (@{$hpseudogene->get_all_Transcripts}) {
        my $tbiotype = $htran->biotype."_havana";
        $htran->biotype($tbiotype);
      }
      push (@hpseudogenes, $hpseudogene);
    }
  }


  print STDERR "Retrieved ".scalar(@hgenes)." genes of types: ".join(", ",@{$GB_HAVANA_INPUT_GENETYPE})."\n";
  print STDERR "Retrieved ".scalar(@hprocessedgenes)." 'processed transcript' genes of types: ".join(", ",@{$GB_HAVANA_PROCESSED_GENETYPE})."\n";
  print STDERR "Retrieved ".scalar(@hpseudogenes)." pseudogenes of types: ".join(", ",@{$GB_HAVANA_PSEUDO_GENETYPE})."\n";

  # We want to keep the HAVANA genes as they are as they asked us to keep their gene clustering untouched
  $self->havana_Coding_Genes(@hgenes);
  $self->havana_Processed_Transcript_Genes(@hprocessedgenes);
  $self->havana_Pseudogene_Genes(@hpseudogenes);

  #push(@genes, @hgenes);
  #push(@processedgenes, @hprocessedgenes);
  #push(@pseudogenes, @hpseudogenes);

  # We processed the ensembl genes to remove any trace of previous transcript merges in case 
  # you are running a merge update instead of a fresh clean merge
  @transcripts = $self->check_merge_transcript_status(@genes);
  @processedtranscripts = $self->check_merge_transcript_status(@processedgenes);
  @pseudotranscripts = $self->check_merge_transcript_status(@pseudogenes);

  # Join all the gene set together
  #push(@genes, @hgenes);
  
  print STDERR "Finished fetching genes\n";
  $self->combined_Transcripts(@transcripts);
  $self->combined_Processed_Transcripts(@processedtranscripts);
  $self->combined_PseudoTranscripts(@pseudotranscripts);
}

sub check_merge_transcript_status{
  my ($self, @genes) = @_;

  print "Checking premerge gene status\n";

  my @transcripts;
  foreach my $gene(@genes){
  TRANSCRIPT:
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      #First we remove HAVANA only transcripts that are present in merged genes
      if($gene->analysis->logic_name() eq $MERGED_GENE_LOGIC_NAME &&
         $tran->analysis->logic_name() eq $HAVANA_LOGIC_NAME){
        next TRANSCRIPT;
      }elsif($tran->analysis->logic_name() eq $MERGED_TRANSCRIPT_LOGIC_NAME){
        # In case of a merged transcript we want to distinguish the ones that came from HAVANA that have same CDS 
        # but different UTR structure as we want to remove then. This is important for a merge update to avoid 
        # then been wrongly identified as share CDS and UTR 
        my $share_enst = 0;
        my $share_cds_and_utr = 0;
        my @dbentries = @{ $tran->get_all_DBEntries };
        foreach my $dbentry (@dbentries){
          if ($dbentry->dbname eq "shares_CDS_with_ENST"){
            #print "On transcript: ",$tran->dbID," This is a HAVANA shares ENST\n";
            #next TRANSCRIPT;
            $share_enst = 1;
          }
          if ($dbentry->dbname eq "shares_CDS_and_UTR_with_OTTT"){
            $share_cds_and_utr = 1;
            #print "On transcript: ",$tran->dbID," This is a HAVANA shares CDS and UTR\n";
          }
        }
        if ($share_enst == 1 && $share_cds_and_utr == 0){
          next TRANSCRIPT;
        }
      }
      
      $self->flush_xref($tran);
      
      #Check if a transcript is in the discarded genes database before adding it to the merging list.
      if($self->check_transcript_in_discarded_db($tran) != 0){
        #print "Transcript added\n";
        push(@transcripts, $tran);
        
      }
    }
  }
  return @transcripts;
}

sub check_transcript_in_discarded_db{
  my ($self, $tran) = @_;
 
  my @exons = @{$tran->get_all_Exons};

  my $discardedslice = $self->discarded_db->get_SliceAdaptor->fetch_by_region('toplevel',$tran->slice->seq_region_name,$tran->seq_region_start,$tran->seq_region_end);
  #print STDERR "Fetching discarded genes\n"; 
  #print "NUMBER OF DISCARDED GENES: ",scalar(@{$discardedslice->get_all_Genes}),"\n"; 
  DGENE: 
  foreach my $dgene (@{$discardedslice->get_all_Genes}){
    DTRANS:foreach my $dtran (@{$dgene->get_all_Transcripts}){
      my @dexons = @{$dtran->get_all_Exons};
      if(scalar(@exons) == scalar(@dexons)){
        #print "Number of exons: ",scalar(@exons),"\n";
        for (my $i=0; $i < scalar(@exons); $i++){

          if ($exons[$i]->seq_region_start   != $dexons[$i]->seq_region_start ||
              $exons[$i]->strand  != $dexons[$i]->strand ||
              $exons[$i]->seq_region_end     != $dexons[$i]->seq_region_end){
            # if you enter here means that these two transcripts are not the same
            #print "transcript exon coordinates are different\n";
            next DTRANS;
          }
        }
        # If you are here means that both transcripts are the same and $trans must be discarded
        #print "transcript found in discarded db\n";
        return 0;
      }else{
      # if you enter here means that these two transcripts are not the same
        #print "transcript number of exons is different\n";
        next DGENE;
      }
    }
  }
  #If we reach here means that no transcript in the discarded db is the same as our transcript so we keep it
  return 1;
}

sub flush_xref{
  my ($self, $transcript) = @_;

  my @newxrefs;
  #print "THIS IS WHAT NEED EMPTYING: ",$transcript->get_all_DBEntries,"\n";
  foreach my $tran_xref (@{$transcript->get_all_DBEntries}){
    if ($tran_xref->dbname ne "shares_CDS_and_UTR_with_OTTT" &&
        $tran_xref->dbname ne "shares_CDS_with_OTTT" &&
        $tran_xref->dbname ne "shares_CDS_with_ENST" &&
        $tran_xref->dbname ne "OTTT" &&
        $tran_xref->dbname ne "OTTP" &&
        $tran_xref->dbname ne "OTTG"){
      push (@newxrefs, $tran_xref);
    }
  }

# The naughty bit!
  $transcript->{dbentries} = [];
  #$transcript->{display_xref} = [];

  foreach my $newxref (@newxrefs) {
    $transcript->add_DBEntry($newxref);
  }

 # return scalar(@newtrans);
}

###########################################################c

=head2 cluster_Transcripts

 Description : It separates transcripts according to strand and then clusters 
               each set of transcripts by calling _cluster_Transcripts_by_genomic_range()
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster

=cut

sub cluster_Transcripts {
  my ($self,@transcripts) = @_;
 
  my @forward_transcripts;
  my @reverse_transcripts;
 
  foreach my $transcript (@transcripts){
    my @exons = @{ $transcript->get_all_Exons };
    if ( $exons[0]->strand == 1 ){
      push( @forward_transcripts, $transcript );
    }
    else{
      push( @reverse_transcripts, $transcript );
    }
  }
  
  my @forward_clusters;
  my @reverse_clusters;
  
  if ( @forward_transcripts ){
    @forward_clusters = $self->_cluster_Transcripts_by_genomic_range( @forward_transcripts );
  }
  if ( @reverse_transcripts ){
    @reverse_clusters = $self->_cluster_Transcripts_by_genomic_range( @reverse_transcripts );
  }
  my @clusters;
  if ( @forward_clusters ){
    push( @clusters, @forward_clusters);
  }
  if ( @reverse_clusters ){
    push( @clusters, @reverse_clusters);
  }
  return @clusters;
}

############################################################

=head2 _cluster_Transcripts_by_genomic_range

 Description : It clusters transcripts according to genomic overlap
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster

=cut

sub _cluster_Transcripts_by_genomic_range{
  my ($self,@mytranscripts) = @_;
  # first sort the transcripts

  my @transcripts = sort { $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end } @mytranscripts;


  # create a new cluster 
  my $cluster=Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new();
  my $count = 0;
  my @cluster_starts;
  my @cluster_ends;
  my @clusters;
  
  # put the first transcript into these cluster
  $cluster->put_Transcripts( $transcripts[0] );

  $cluster_starts[$count] = $transcripts[0]->start;
  $cluster_ends[$count]   = $transcripts[0]->end;
  
  # store the list of clusters
  push( @clusters, $cluster );
  
  # loop over the rest of the transcripts
 LOOP1:
  for (my $c=1; $c<=$#transcripts; $c++){
    #print STDERR "\nIn cluster ".($count+1)."\n";
    #print STDERR "start: $cluster_starts[$count] end: $cluster_ends[$count]\n";
    #print STDERR "comparing:\n";
    #Bio::EnsEMBL::Analysis::Tools::TranscriptUtils->_print_Transcript( $transcripts[$c] );
    
    if ( !( $transcripts[$c]->end < $cluster_starts[$count] ||
	    $transcripts[$c]->start > $cluster_ends[$count] ) ){
      $cluster->put_Transcripts( $transcripts[$c] );
      
      # re-adjust size of cluster
      if ($transcripts[$c]->start < $cluster_starts[$count]) {
	$cluster_starts[$count] = $transcripts[$c]->start;
      }
      if ( $transcripts[$c]->end > $cluster_ends[$count]) {
	$cluster_ends[$count] =  $transcripts[$c]->end;
      }
    }
    else{
      # else, create a new cluster with this feature
      $count++;
      $cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new();
      $cluster->put_Transcripts( $transcripts[$c] );
      $cluster_starts[$count] = $transcripts[$c]->start;
      $cluster_ends[$count]   = $transcripts[$c]->end;
      
      # store it in the list of clusters
      push(@clusters,$cluster);
    }
  }
  return @clusters;
}

############################################################

=head2 cluster_into_Genes

    Example :   my @genes = $self->cluster_into_Genes(@transcripts);
Description :   it clusters ensembl transcripts and havana genes according to exon overlap.
                It will take care of difficult cases like transcripts within introns.
                It also unify exons that are shared among transcripts.
    Returns :   a beautiful list of geen objects
    Args    :   a list of transcript objects

=cut

sub cluster_into_Genes{
  my ($self, $transcripts_unsorted, $havana_genes) = @_;
  
  my %ottg_xref;

  my $num_trans = scalar(@{$transcripts_unsorted});

  # First clean the coding exon cache in case it has any exons stored from previous called to the cluster_into_Genes function.
  $self->clear_coding_exons_cache;

  my @transcripts_unsorted_translation;

  foreach my $tran(@{$transcripts_unsorted}){
    if ($tran->translation){
      push (@transcripts_unsorted_translation, $tran);
    }
  }

  my @transcripts = sort { $a->coding_region_start <=> $b->coding_region_start ? $a->coding_region_start <=> $b->coding_region_start  : $b->coding_region_end <=> $a->coding_region_end } @transcripts_unsorted_translation;

  # We are going to set the havana genes as initial clusters as we want to keep their structure
  my @clusters;

  foreach my $hav_gene(@{$havana_genes}){
    my $ottg_key;
    foreach my $entry(@{ $hav_gene->get_all_DBEntries}){
      #print "ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
      if ($entry->dbname eq 'Vega_gene'){
        if($entry->primary_id eq $entry->display_id){
          $ottg_key = $entry->primary_id;
        } 
      }
    }

    my @h_clust;
    foreach my $hav_trans(@{$hav_gene->get_all_Transcripts}){
    $ottg_xref{$hav_trans} = $ottg_key;

      push (@h_clust,$hav_trans);
    }
    push(@clusters,\@h_clust);
  }
  #my @clusters = @{$havana_genes};

  # clusters transcripts by whether or not any coding exon overlaps with a coding exon in 
  # another transcript. We will use the set of Havana
  foreach my $tran (@transcripts) {
    
    my @matching_clusters;
  CLUSTER: 
    foreach my $cluster (@clusters) {
      
      #print "Transcript: ",$tran->stable_id," has coding region start: ",$tran->coding_region_start,"\n";
      #print "CLUSTER: ",$cluster,"\n";
      #foreach my $cluster_transcript (@{$cluster->get_all_Transcripts}) {
      foreach my $cluster_transcript (@$cluster) {
        #print $cluster_transcript,"\n";
        if ($cluster_transcript->translation){
          if ($tran->coding_region_end  >= $cluster_transcript->coding_region_start &&
              $tran->coding_region_start <= $cluster_transcript->coding_region_end) {
            
            # foreach my $exon1 (@{$tran->get_all_Exons}) {
            # foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
            my $exons1 = $self->get_coding_exons_for_transcript($tran);
            my $cluster_exons = $self->get_coding_exons_for_transcript($cluster_transcript);
            
            foreach my $exon1 (@{$exons1}) {
              foreach my $cluster_exon (@{$cluster_exons}) {
                
                if ($exon1->overlaps($cluster_exon) && $exon1->strand == $cluster_exon->strand) {
                  push (@matching_clusters, $cluster);
                  next CLUSTER;
                }
              }
            }
          }
        }
      }
    }
    
    if (scalar(@matching_clusters) == 0) {
      my @newcluster;
      push(@newcluster,$tran);
      push(@clusters,\@newcluster);
   
    } 
    elsif (scalar(@matching_clusters) == 1) {
      #print "transcript_id ",$matching_clusters[0]->[0]->dbID,"\n";
      #print "transcript_id db ",$tran->dbID,"\n";
      push @{$matching_clusters[0]}, $tran;
  
    } 
    else {
      print "BEWARE YOU ARE HERE MERGING CLUSTERS\n";
      
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my @merged_cluster;
      my $hav_gene_counter = 0;
      
      # Check how many havana genes are in the cluster
      CLUST:foreach my $clust (@matching_clusters) {
        foreach my $clust_trans (@{$clust}){
          #print "Cluster ", $clust_trans,"\n";
          if ($clust_trans->biotype =~ /havana/){
            my $hav_gene_counter++;
            next CLUST;
          }
        }
      }
      
      # Merge all clustered genes in case there is only one or none havana genes.
      if ($hav_gene_counter <= 1){
        foreach my $clust (@matching_clusters) {
          #print "THIS GENE COUNTER ",$clust,"\n";
          push @merged_cluster, @$clust;
        }
        push @merged_cluster, $tran; 

        # ADD SOME NICE PRINT OUTPUT
        print "MERGING CLUSTER COVERED BY TRANSCRIPT IN: ",$tran->seq_region_name," - ",$tran->seq_region_start," - ",$tran->seq_region_end,"\n";

        #print "content: ",join(' - ',@merged_cluster),"\n";

      }else{
      
        # If there is more than one Havana gene we check which one has a better overlap on the merged transcript
        # and we also check if the havana genes are partial and need to be merged.
        #NEED TO SORT CLUSTER

        @matching_clusters = sort { $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end } @matching_clusters;
        
        # IF MANY CLUSTs equally good then merge all of them. NOPE THIS IS NOT ALLOWED ANY MORE
        my @clust2merge;
        my $cluster_ref = 0;
        my $best_cluster = 0;
        my $best_match = 0;
        my $tran_length = $tran->translate->length;
        foreach my $clust (@matching_clusters) {
          print "More than one Havana gene overlaps. Scoring the best cluster \n";
          #exit;
          my $match_percent = $self->gene_overlap_percent($clust,$tran);
          ## THIS WAS THE BIT TO MERGE TWO HAVANAS 
          #if ($match_percent > 100){
          #  push (@clust2merge, $clust);
          #  next;
          #}
          if ($match_percent > $best_match){
            $best_match = $match_percent;
            $best_cluster = $cluster_ref;
          }
          $cluster_ref++;
        }
        #if (scalar(@clust2merge) >= 1){  
        #  print "Have a perfect overlap with a Havana gene\n";
        #  push (@merged_cluster, @clust2merge);
        #}else{
          print "Only one havana gene gave the best overlap\n";
          push (@merged_cluster, $matching_clusters[$best_cluster]);
        #}
        push (@merged_cluster, $tran);
      }
      
      push @new_clusters,\@merged_cluster;
      # Add back non matching clusters
    MATCHING:
      foreach my $clust (@clusters) {
        
        foreach my $tran_clust (@$clust){
          foreach my $m_clust (@merged_cluster) {
#	foreach my $m_clust (@matching_clusters) {
            if ($tran_clust == $m_clust) {
              next MATCHING;
            }
          }
        }
        push @new_clusters,$clust;
      }
      @clusters =  @new_clusters;
    }
  }
  
  # safety and sanity checks
  $self->check_Clusters(scalar(@transcripts), \@clusters);
  
  # Merge the clustered ensembl transcript with the Havana genes.
  # make and store genes
  #print STDERR scalar(@clusters)." created, turning them into genes...\n";
  my @genes;
  foreach my $cluster(@clusters){
    my $count = 0;
    my $ottg_added = "";
    my $gene = new Bio::EnsEMBL::Gene;
    foreach my $transcript (@$cluster){
      #print "Transcript Stable ID: ",$transcript->dbID,"\n";
      $gene->add_Transcript($transcript);
      if ($ottg_xref{$transcript} && ($ottg_added ne $ottg_xref{$transcript})){
        $self->add_ottg_xref($gene,$ottg_xref{$transcript});
        $ottg_added = $ottg_xref{$transcript};
      }
 
    }
    push( @genes, $gene );
  }

  return @genes;
}

############################################################
=head2 cluster_into_PseudoGenes

    Example :   my @genes = $self->cluster_into_Genes(@transcripts);
Description :   it clusters transcripts into genes according to exon overlap.
                It will take care of difficult cases like transcripts within introns.
                It also unify exons that are shared among transcripts.
    Returns :   a beautiful list of geen objects
    Args    :   a list of transcript objects

=cut

sub cluster_into_PseudoGenes{
  my ($self, $transcripts_unsorted, $havana_non_coding) = @_;

  my %ottg_xref;

  my $num_trans = scalar(@{$transcripts_unsorted});

  my @transcripts = sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end } @{$transcripts_unsorted};

  # I set the havana non coding transcripts as the initial set of cluster
  my @clusters;

  foreach my $hav_gene(@{$havana_non_coding}){

    my $ottg_key;

    foreach my $entry(@{ $hav_gene->get_all_DBEntries}){
      #print "ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
      if ($entry->dbname eq 'Vega_gene'){
        if($entry->primary_id eq $entry->display_id){
          $ottg_key = $entry->primary_id;
        } 
      }
    }

    my @h_clust;
    foreach my $hav_trans(@{$hav_gene->get_all_Transcripts}){
      $ottg_xref{$hav_trans} = $ottg_key;
      push (@h_clust,$hav_trans);
    }
    push(@clusters,\@h_clust);
  }

  # my @clusters = @{$havana_non_coding};

  foreach my $tran (@transcripts) {
    
    my @matching_clusters;
    CLUSTER: 
    foreach my $cluster (@clusters) {
      
     #print "Transcript: ",$tran->stable_id," has start: ",$tran->start,"\n";

     # foreach my $cluster_transcript (@{$cluster->get_all_Transcripts}) {
      foreach my $cluster_transcript (@$cluster) {
        if ($tran->end  >= $cluster_transcript->start &&
            $tran->start <= $cluster_transcript->end) {
          
          # foreach my $exon1 (@{$tran->get_all_Exons}) {
          # foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
          my $exons1 = $tran->get_all_Exons();
          my $cluster_exons = $cluster_transcript->get_all_Exons();

          foreach my $exon1 (@{$exons1}) {
            foreach my $cluster_exon (@{$cluster_exons}) {
              
              if ($exon1->overlaps($cluster_exon) && $exon1->strand == $cluster_exon->strand) {
                push (@matching_clusters, $cluster);
                next CLUSTER;
              }
            }
          }
        }
      }
    }
    
    if (scalar(@matching_clusters) == 0) {
      my @newcluster;
      push(@newcluster,$tran);
      push(@clusters,\@newcluster);
    } 
    elsif (scalar(@matching_clusters) == 1) {
      push @{$matching_clusters[0]}, $tran;
      
    } 
    else {
      print "BEWARE YOU ARE HERE MERGING CLUSTERS\n";      
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my @merged_cluster;
      my $hav_gene_counter = 0;


      # Check how many havana genes are in the cluster
      CLUST:foreach my $clust (@matching_clusters) {
        foreach my $clust_trans (@{$clust}){
          print "Cluster ", $clust_trans,"\n";
          if ($clust_trans->biotype =~ /havana/){
            my $hav_gene_counter++;
            next CLUST;
          }
        }
      }
      
      # Merge all clustered genes in case there is only one or none havana genes.
      if ($hav_gene_counter <= 1){
        foreach my $clust (@matching_clusters) {
          push @merged_cluster, @$clust;
        }
        push @merged_cluster, $tran; 
      }else{
        # If there is more than one Havana gene we check which one has a better overlap on the merger transcript
        # and we also check if the havana genes are partial and need to be merged.
        #NEED TO SORT CLUsTER
        @matching_clusters = sort { $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end } @matching_clusters;
        
        my $cluster_ref = 0;
        my $best_cluster = 0;
        my $best_match = 0;
        my $tran_length = $tran->length;
        foreach my $clust (@matching_clusters) {
          my $match_percent = $self->gene_overlap_percent($clust,$tran);
          if ($match_percent > $best_match){
            $best_match = $match_percent;
            $best_cluster = $cluster_ref;
          }
          $cluster_ref++;
        }
        
        push (@merged_cluster, $matching_clusters[$best_cluster]);
        push (@merged_cluster, $tran);
      }
      
      push @new_clusters,\@merged_cluster;
      # Add back non matching clusters
    MATCHING:
      foreach my $clust (@clusters) {
        
        foreach my $tran_clust (@$clust){
          foreach my $m_clust (@merged_cluster) {
            if ($tran_clust == $m_clust) {
              next MATCHING;
            }
          }
        }
        push @new_clusters,$clust;
      }
      @clusters =  @new_clusters;
    }
  }

  # safety and sanity checks
  $self->check_Clusters(scalar(@transcripts), \@clusters);
  
  # make and store genes
  #print STDERR scalar(@clusters)." created, turning them into genes...\n";
  my @genes;
  foreach my $cluster(@clusters){
    my $count = 0;
    my $ottg_added = "";
    my $gene = new Bio::EnsEMBL::Gene;
    foreach my $transcript (@$cluster){
      #print "Transcript Stable ID: ",$transcript->dbID,"\n";
      $gene->add_Transcript($transcript);
      if ($ottg_xref{$transcript} && ($ottg_added ne $ottg_xref{$transcript})){
        $self->add_ottg_xref($gene,$ottg_xref{$transcript});
        $ottg_added = $ottg_xref{$transcript};
      }
    }
    push( @genes, $gene );
  }
  return @genes;
}

sub combine_gene_clusters {
  my ($self,$preliminary_coding_genes,$preliminary_pseudo_genes) = @_;

 # my @total_genes;
  my @coding_genes = sort { $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end } @{$preliminary_coding_genes};
  my @pseudo_genes = sort { $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end } @{$preliminary_pseudo_genes};

  my @unclestered_pseudos;
     
  CLUSTER:
  foreach my $pseudo_gene(@pseudo_genes){
    if ($pseudo_gene->biotype =~/havana/){
      push(@unclestered_pseudos,$pseudo_gene);
      next CLUSTER;
    }
    my $pseudo_status = 0;
    OVERLAP:
    foreach my $coding_gene (@coding_genes){
      
      if ($coding_gene->end >= $pseudo_gene->start &&
          $coding_gene->start <= $pseudo_gene->end){
        
        my $coding_length = $self->get_coding_length($coding_gene);

        my $cg_exons = $self->get_coding_exons_for_gene($coding_gene);
        
        my $pg_exons = $pseudo_gene->get_all_Exons();
        
        foreach my $cg_exon (@{$cg_exons}) {
          foreach my $pg_exon (@{$pg_exons}) {
            
            if ($cg_exon->overlaps($pg_exon) && $cg_exon->strand == $pg_exon->strand) {
              #check uf the overlap covers at least 10 percent of the coding region of the longest transcript in the gene.
              # This check is a bit experimental.
              if($self->overlap_percent($cg_exon,$pg_exon,$coding_length) > 10){ 
                # Have to add all the transcripts of the pseudo to the gene and remove the pseudogene
                foreach my $p_transcript(@{$pseudo_gene->get_all_Transcripts}){
                  $coding_gene->add_Transcript($p_transcript);
                }
                $pseudo_status = 1;
                next OVERLAP;
              }
            }
          }
        }
      }
    }
    unless ($pseudo_status == 1){ push(@unclestered_pseudos,$pseudo_gene);}
  }

  push(@coding_genes, @unclestered_pseudos);

  return @coding_genes; 
}

sub get_coding_length {
  my ($self, $gene) = @_;

  my $length = 0;

  foreach my $transcript(@{$gene->get_all_Transcripts }){
    if($transcript->translate){
      if($transcript->translate->length > $length){
        $length = $transcript->translate->length;
      }
    }
  }
  return $length;
}

sub gene_overlap_percent {

  my ($self, $gene,$transcript) = @_;
  
  print "GENE CLUST: ",$gene,"\n";

  my $trans_length;
  my @t_exons;
  if ($transcript->translate){
    $trans_length = $transcript->translate->length;
    @t_exons = @{$transcript->get_all_translateable_Exons};
  }else{
    @t_exons = @{$transcript->get_all_Exons};
    $trans_length = $transcript->length;
  }
  
  my $max_overlap = 0;
  my $min_g_t_overlap = 100;
  foreach my $g_trans (@{$gene->get_all_transcripts}){
    my $percent = 0;
    my $g_t_length = $g_trans->length;
    my $g_t_percent = 0;
    foreach my $g_exon (@{g_trans->get_all_Exons}){
      
      foreach my $t_exon (@t_exons){
        if($g_exon->seq_region_start < $t_exon->seq_region_end &&
           $g_exon->seq_region_end   > $t_exon->seq_region_start){
          
          my $low = 0;
          my $high = 0;
          
          if ($g_exon->seq_region_start >= $t_exon->seq_region_start){
            $low = $g_exon->seq_region_start;
          }else{
            $low = $t_exon->seq_region_start;
          }
          if ($g_exon->seq_region_end <= $t_exon->seq_region_end){
            $high = $g_exon->seq_region_end;
          }else{
            $high = $t_exon->seq_region_end;
          }
          
          my $overlap_length = $high-$low;
          $g_t_percent += ($overlap_length/$g_t_length)*100;
          $percent += ($overlap_length/$trans_length)*100;
        }
      }
    }
    # First we check what is the minimun coverage of the havana transcripts in the gene cluster
    # If this value is close to 100% we consider that the ensembl transcript and the gene should
    # be merged no matter what. This is usefull if there are a few partial genes from havana that
    # overlap a transcript in Ensembl
    if ($g_t_percent <$min_g_t_overlap){
      $min_g_t_overlap = $g_t_percent;
    }

    # We also what to know which gene from havana has the best coverage over the ensembl transcript
    # We will use this value in case we have to Havana genes that we don't want to merge and we want
    # to decide which one should merge with the ensembl transcript
    if ($percent >$max_overlap){
      $max_overlap = $percent;
    } 
  }

## This bit of code was to merge havana transcripts that are completely covered by a ensembl transcript. This is not allowed anymore
#  if ($min_g_t_overlap > 98){
#    return 1000;# WE want to merge these two no matter what
#  }else{ 
    return $max_overlap;
#  }
}


sub overlap_percent {

  my ($self, $cg_exon,$pg_exon,$coding_length) = @_;
  my $low = 0;
  my $high = 0;

  if ($cg_exon->seq_region_start >= $pg_exon->seq_region_start){
    $low = $cg_exon->seq_region_start;
  }else{
    $low = $pg_exon->seq_region_start;
  }
  if ($cg_exon->seq_region_end <= $pg_exon->seq_region_end){
    $high = $cg_exon->seq_region_end;
  }else{
    $high = $pg_exon->seq_region_end;
  }

  my $overlap_length = $high-$low;

  my $percent = ($overlap_length/$coding_length)*100;

  return $percent;

}


############################################################

=head2 get_coding_exons_for_transcript

    Example :    my $exons1 = $self->get_coding_exons_for_transcript($tran);
Description :   It returns the coding exons of a transcript and stores 
                them in a hash to safe computer time                
    Returns :   An ArrayRef than contain Exon objects.
    Args    :   a transcript object

=cut

{
  my %coding_exon_cache;

  sub clear_coding_exons_cache {
    %coding_exon_cache = ();
  }


sub get_coding_exons_for_transcript {
    my ($self, $trans) = @_;

    if (exists($coding_exon_cache{$trans})) {
      return $coding_exon_cache{$trans};
    } else {
      my %coding_hash;
      
      next if (!$trans->translation);
      foreach my $exon (@{$trans->get_all_translateable_Exons}) {
        $coding_hash{$exon} = $exon;
      }

      my @coding = sort { $a->start <=> $b->start } values %coding_hash;
      #my @coding = values %coding_hash;

      $coding_exon_cache{$trans} = \@coding;
      return $coding_exon_cache{$trans};
    }
  }
}

=head2 get_coding_exons_for_gene

    Example :    my $exons1 = $self->get_coding_exons_for_gene($gene);
Description :   It returns the coding exons of a transcript and stores 
                them in a hash to safe computer time                
    Returns :   An ArrayRef than contain Exon objects.
    Args    :   a transcript object

=cut


sub get_coding_exons_for_gene {
  my ($self, $gene) = @_;
  
  my @coding;
  
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    next if (!$trans->translation);
    foreach my $exon (@{$trans->get_all_translateable_Exons}) {
      push @coding, $exon;
    }
  }
  
  return \@coding;
  
}



############################################################

sub check_Clusters{
  my ($self, $num_transcripts, $clusters) = @_;
  #Safety checks
  my $ntrans = 0;

  my $cluster_num = 0;

  my %trans_check_hash;
  foreach my $cluster (@$clusters) {
    $ntrans += scalar(@$cluster);

    foreach my $trans (@$cluster) {

      if (defined($trans_check_hash{$trans})) {
        $self->throw("Transcript " . $trans->dbID . " added twice to clusters\n");
      }
      $trans_check_hash{$trans} = 1;
    }
    if (!scalar(@$cluster)) {
      $self->throw("Empty cluster");
    }
  }
  if ($ntrans < $num_transcripts) {
    $self->throw("Not all transcripts have been added into clusters $ntrans and " . $num_transcripts. " \n");
  } 
  #end safety checks
  return;
}


############################################################

sub transcript_high{
  my ($self,$tran) = @_;
  my $high;
  #$tran->sort;
  if ( $tran->start_Exon->strand == 1){
    $high = $tran->end_Exon->end;
  }
  else{
    $high = $tran->start_Exon->end;
  }
  return $high;
}

############################################################

sub transcript_low{
  my ($self,$tran) = @_;
  my $low;
  #$tran->sort;
  if ( $tran->start_Exon->strand == 1){
    $low = $tran->start_Exon->start;
  }
  else{
    $low = $tran->end_Exon->start;
  }
  return $low;
}

############################################################

sub by_transcript_high {
  my $alow;
  my $blow;

  my $ahigh;
  my $bhigh;
  
  # alow and ahigh are the left most and right most coordinates for transcript $a 
  if ($a->start_Exon->strand == 1) {
    $alow  = $a->start_Exon->start;
    $ahigh = $a->end_Exon->end;
  } 
  else {
    $alow  = $a->end_Exon->start;
    $ahigh = $a->start_Exon->end;
  }

  # blow and bhigh are the left most and right most coordinates for transcript $b 
  if ($b->start_Exon->strand == 1) {
    $blow  = $b->start_Exon->start;
    $bhigh = $b->end_Exon->end;
  } 
  else {
    $blow  = $b->end_Exon->start;
    $bhigh = $b->start_Exon->end;
  }

  # return the ascending comparison of the right-most coordinates if they're different
  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } 
  # if they'r equal, return the ascending comparison of the left most coordinate
  else {
    return $alow <=> $blow;
  }
}



############################################################

sub prune_Exons {
  my ($self,$gene) = @_;
  
  my @unique_Exons; 
  
  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_Exon and translation->end_Exon
  
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty
    UNI:foreach my $uni (@unique_Exons) {
	if ($uni->start  == $exon->start  &&
	    $uni->end    == $exon->end    &&
	    $uni->strand == $exon->strand &&
	    $uni->phase  == $exon->phase  &&
	    $uni->end_phase == $exon->end_phase
	   ) {
	  $found = $uni;
	  last UNI;
	}
      }
      if (defined($found)) {
	push(@newexons,$found);
        if($tran->translation){
          if ($exon == $tran->translation->start_Exon){
            $tran->translation->start_Exon($found);
          }
          if ($exon == $tran->translation->end_Exon){
            $tran->translation->end_Exon($found);
          }
        }
      } else {
	push(@newexons,$exon);
	push(@unique_Exons, $exon);
      }
    }          
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
   
    #print "Uniq_tran sid: ",$tran->dbID,"\n";

  }
  return $gene;
}

############################################################

=head2 prune_features

 Description: prunes out duplicated features
 Returntype : array of Bio::EnsEMBL::SeqFeature
 Args       : array of Bio::EnsEMBL::SeqFeature

=cut
    
sub prune_features {
    my ($self,$feature_hash)  = @_;
    my @pruned;
 
  ID:
    foreach my $id (keys %{ $feature_hash }) {
	my @features = @{$feature_hash->{$id}};
	@features = sort {$a->start <=> $b->start} @features;
	
	unless ( @features ){
	    print STDERR "No features here for id: $id\n";
	    next ID;
	}
	while ( @features && !defined $features[0] ){
	    #print STDERR "jumping an undefined feature\n";
	    shift @features;
	}
	
	my $prev = -1;
	
      FEATURE: 
	foreach  my $f (@features) {
	    if ($prev != -1 && $f->hseqname eq $prev->hseqname &&
		$f->start   == $prev->start &&
		$f->end     == $prev->end   &&
		$f->hstart  == $prev->hstart &&
		$f->hend    == $prev->hend   &&
		$f->strand  == $prev->strand &&
		$f->hstrand == $prev->hstrand) 
	    {
		#keep the one with highest score
		if ( $f->score > $prev->score ){
		    $prev->score( $f->score );
		}
		#print STDERR "pruning duplicated feature\n";
		#print STDERR "previous: ".$prev->gffstring."\n";
		#print STDERR "thisone : ".$f->gffstring."\n";
		next FEATURE;
	    } 
	    else {
		push(@pruned,$f);
		$prev = $f;
	    }
	}
    }
    return @pruned;
}

############################################################


############################################################
#
# GETSET METHODS
#
############################################################

# get/set method holding a reference to the db with genewise and combined genes,
# havana genes and discarded genes
# this reference is set in Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder

sub ensembl_db{
 my ($self,$ensembl_db) = @_;
 if ( $ensembl_db ){
   $self->{_ensembl_db} = $ensembl_db;
 }
 
 return $self->{_ensembl_db};
}

sub havana_db{
 my ($self,$havana_db) = @_;
 if ( $havana_db ){
   $self->{_havana_db} = $havana_db;
 }
 
 return $self->{_havana_db};
}

sub discarded_db{
  my ($self, $discarded_db) = @_;

  if ( $discarded_db ){
    $self->{_discarded_db} = $discarded_db;;
  }

  return $self->{_discarded_db};
}

############################################################

sub havana_Coding_Genes {
    my ($self,@genes) = @_;

    if (!defined($self->{_havana_Coding_Genes})) {
        $self->{_havana_Coding_Genes} = [];
    }

    if (scalar @genes > 0) {
	push(@{$self->{_havana_Coding_Genes}},@genes);
    }

    return @{$self->{_havana_Coding_Genes}};
}

sub havana_Processed_Transcript_Genes {
    my ($self,@genes) = @_;

    if (!defined($self->{_havana_Processed_Transcript_Genes})) {
        $self->{_havana_Processed_Transcript_Genes} = [];
    }

    if (scalar @genes > 0) {
	push(@{$self->{_havana_Processed_Transcript_Genes}},@genes);
    }

    return @{$self->{_havana_Processed_Transcript_Genes}};
}

sub havana_Pseudogene_Genes {
    my ($self,@genes) = @_;

    if (!defined($self->{_havana_Pseudogene_Genes})) {
        $self->{_havana_Pseudogene_Genes} = [];
    }

    if (scalar @genes > 0) {
	push(@{$self->{_havana_Pseudogene_Genes}},@genes);
    }

    return @{$self->{_havana_Pseudogene_Genes}};
}

sub combined_Transcripts {
    my ($self,@transcripts) = @_;

    if (!defined($self->{_coding_transcripts})) {
        $self->{_coding_transcripts} = [];
    }

    if (scalar @transcripts > 0) {
	push(@{$self->{_coding_transcripts}},@transcripts);
    }

    return @{$self->{_coding_transcripts}};
}

sub combined_Processed_Transcripts {
    my ($self,@transcripts) = @_;

    if (!defined($self->{_processed_transcripts})) {
        $self->{_processed_transcripts} = [];
    }

    if (scalar @transcripts > 0) {
	push(@{$self->{_processed_transcripts}},@transcripts);
    }

    return @{$self->{_processed_transcripts}};
}


sub combined_PseudoTranscripts {
    my ($self,@pseudotranscripts) = @_;

    if (!defined($self->{_pseudo_transcripts})) {
        $self->{_pseudo_transcripts} = [];
    }

    if (scalar @pseudotranscripts > 0) {
	push(@{$self->{_pseudo_transcripts}},@pseudotranscripts);
    }

    return @{$self->{_pseudo_transcripts}};
}

=head2 update_biotypes

  Description: This check the biotypes of the merged transcript and genes and updates then to reflect the merge

=cut

sub update_biotypes{
  my ($self, @genes) = @_;

  my %pseudobiotypes;
  my %processedbiotypes;
  my %coding_biotypes;

  foreach my $epb (@{$GB_ENSEMBL_PSEUDO_GENETYPE}){
    $pseudobiotypes{$epb}=1;
  }
  foreach my $hpb (@{$GB_HAVANA_PSEUDO_GENETYPE}){
    $pseudobiotypes{$hpb."_havana"}=1;
  }
  foreach my $epb (@{$GB_ENSEMBL_PROCESSED_GENETYPE}){
    $processedbiotypes{$epb}=1;
  }
  foreach my $hpb (@{$GB_HAVANA_PROCESSED_GENETYPE}){
    $processedbiotypes{$hpb."_havana"}=1;
  }

  foreach my $ecb (@{$GB_ENSEMBL_INPUT_GENETYPE}){
    $coding_biotypes{$ecb}=1;
  }
  foreach my $hcb (@{$GB_HAVANA_INPUT_GENETYPE}){
    $coding_biotypes{$hcb."_havana"}=1;
  }

  foreach my $gene (@genes) { 
    my %trans_types;

    my $has_pseudos = 0;
    my $has_processed = 0;
    my $has_coding = 0;
    #$gene->type($GB_GENE_OUTPUT_BIOTYPE);
    # poke the caches
    my %s_pfhash;

    foreach my $tran (@{$gene->get_all_Transcripts}) {
      $trans_types{$tran->biotype} = 1;

      foreach my $pseudobiotype(keys %pseudobiotypes){
        if ($tran->biotype=~ /$pseudobiotype/){
       #   print "HAS PSEUDOS NEW CODE WORKS\n";
          $has_pseudos = 1;
        }
      }
     foreach my $processedbiotype(keys %processedbiotypes){
        if ($tran->biotype=~ /$processedbiotype/){
      #    print "HAS PROCESSED NEW CODE WORKS\n";
          $has_processed = 1;
        }
      }
     foreach my $coding_biotype(keys %coding_biotypes){
        if ($tran->biotype=~ /$coding_biotype/){
         # print "HAS CODING NEW CODE WORKS\n";
          $has_coding = 1;
        }
      }
    }
    
    # MANAGE OUTPUT BIOTYPES BEFORE WRITTING OUTPUT
    my $newbiotype;
    my $biotype_status;
    my $has_havana=0;
    my $has_ensembl=0;
    my $has_merged=0;

    if($has_coding == 1){
      $biotype_status = "protein_coding";
    }elsif($has_processed == 1 && $has_coding == 0){
      $biotype_status = "processed_transcript";
    }elsif($has_pseudos == 1 && $has_coding == 0 && $has_processed == 0){
      $biotype_status = "pseudogene";
    }else{
      print "ERROR: I should not really be here for gene biotype checks\n";
      $biotype_status = "weird_".$gene->biotype;
    }
    
    foreach my $t_biotype (keys %trans_types){
      if($t_biotype =~ /$MERGED_TRANSCRIPT_OUTPUT_TYPE/){
        $has_merged =1;
      }elsif($t_biotype =~ /_havana/){
        $has_havana=1;
      }else{
        $has_ensembl=1;
      }
    }

    if (($has_havana == 1 && $has_ensembl == 1) || $has_merged == 1){
      $newbiotype = $biotype_status.$MERGED_GENE_OUTPUT_BIOTYPE;
      $gene->biotype($newbiotype);
    }elsif($has_havana == 1 && $has_ensembl == 0 && $has_merged == 0){
      $newbiotype = $biotype_status.$HAVANA_GENE_OUTPUT_BIOTYPE;
      $gene->biotype($newbiotype);
    }elsif($has_ensembl == 1 && $has_havana == 0 && $has_merged == 0){
      $newbiotype = $biotype_status.$ENSEMBL_GENE_OUTPUT_BIOTYPE;
      $gene->biotype($newbiotype);
    }else{
      $newbiotype = $biotype_status."weird";
      $gene->biotype($newbiotype);
    }
  } 
}

=head2 final_genes

 Descripton: this holds/returns the final genes produced after clustering transcripts and sharing common exons

=cut

sub final_genes{
  my ($self, @genes) = @_;
  
  if ( @genes ){
    push( @{$self->{_final_genes}}, @genes );
  }
  return @{$self->{_final_genes}};
}

############################################################

=head2 gene_types

 Description: get/set for the type(s) of genes (usually TGE_gw, similarity_genewise and combined_e2g genes) 
              to be used in the genebuilder they get set in new()
              Does not include the ab inition predictions
=cut

sub gene_types {
  my ($self,$type) = @_;

  if (defined($type)) {
     push(@{$self->{_gene_types}},$type);
  }

  return @{$self->{_gene_types}};
}

############################################################

sub features {
  my ($self,@features) = @_;
  
  if (!defined($self->{_feature})) {
    $self->{_feature} = [];
  }
  if ( scalar @features ) {
    push(@{$self->{_feature}},@features);
  }
  return @{$self->{_feature}};
}

############################################################

sub query {
  my ($self,$slice) = @_;
  
  if (defined($slice)) {
    $self->{_query} = $slice;
  }
  return $self->{_query};
}

############################################################

=head2 transfer_supporting_evidence

 Title   : transfer_supporting_evidence
 Usage   : $self->transfer_supporting_evidence($source_exon, $target_exon)
 Function: Transfers supporting evidence from source_exon to target_exon, 
           after checking the coordinates are sane and that the evidence is not already in place.
 Returns : nothing, but $target_exon has additional supporting evidence

=cut

sub transfer_supporting_evidence{
  my ($self, $source_exon, $target_exon) = @_;
  
  my @target_sf = @{$target_exon->get_all_supporting_features};
  #  print "target exon sf: \n";
  #  foreach my $tsf(@target_sf){ print STDERR $tsf; $self->print_FeaturePair($tsf); }
  
  #  print "source exon: \n";
 
  # keep track of features already transferred, so that we do not duplicate
  my %unique_evidence;
  my %hold_evidence;

 SOURCE_FEAT:
  foreach my $feat ( @{$source_exon->get_all_supporting_features}){
    next SOURCE_FEAT unless $feat->isa("Bio::EnsEMBL::FeaturePair");
    
    # skip duplicated evidence objects
    next SOURCE_FEAT if ( $unique_evidence{ $feat } );
    
    # skip duplicated evidence 
    if ( $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } ){
      #print STDERR "Skipping duplicated evidence\n";
      next SOURCE_FEAT;
    }

    #$self->print_FeaturePair($feat);
    
  TARGET_FEAT:
    foreach my $tsf (@target_sf){
      next TARGET_FEAT unless $tsf->isa("Bio::EnsEMBL::FeaturePair");
      
      if($feat->start    == $tsf->start &&
	 $feat->end      == $tsf->end &&
	 $feat->strand   == $tsf->strand &&
	 $feat->hseqname eq $tsf->hseqname &&
	 $feat->hstart   == $tsf->hstart &&
	 $feat->hend     == $tsf->hend){
	
	#print STDERR "feature already in target exon\n";
	next SOURCE_FEAT;
      }
    }
    #print STDERR "from ".$source_exon->dbID." to ".$target_exon->dbID."\n";
    #$self->print_FeaturePair($feat);
    # I may need to add a paranoid check to see that no exons longer than the current one are transferred 
    $target_exon->add_supporting_features($feat);
    $unique_evidence{ $feat } = 1;
    $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } = 1;
  }
}


#fetches sequence from appropriate database

sub fetch_sequence{
  my ($self, $name, $db) = @_;

  my $sa = $db->get_SliceAdaptor; 

  my $slice = $sa->fetch_by_name($name);

  return $slice;
}

1;
