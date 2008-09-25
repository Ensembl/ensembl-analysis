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

use Bio::EnsEMBL::Analysis::Config::GeneBuild::General     qw (
							       GB_INPUTID_REGEX
							      );

use Bio::EnsEMBL::Analysis::Config::HavanaAdder qw (
                                                    GB_ENSEMBL_INPUT_GENETYPE
                                                    GB_HAVANA_INPUT_GENETYPE
                                                    GB_HAVANA_INPUT_TRANSCRIPTTYPES
                                                    MERGED_TRANSCRIPT_OUTPUT_TYPE
                                                    HAVANA_LOGIC_NAME
                                                    MERGED_GENE_LOGIC_NAME
                                                    MERGED_TRANSCRIPT_LOGIC_NAME
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
  my @all_transcripts = $self->combined_Transcripts;
  
  # do a preliminary clustering
  my @preliminary_genes = $self->cluster_into_Genes(@all_transcripts);

  # merge redundant ensembl transcripts which match a havana one
  $self->_merge_redundant_transcripts(\@preliminary_genes);

  # make shared exons unique objects
  my @genes =  $self->_make_shared_exons_unique( @preliminary_genes );
  
  print STDERR scalar(@genes)." genes built\n";
  
  $self->final_genes( @genes );
}


sub _merge_redundant_transcripts{
  my ($self, $genes) = @_;

 GENE:
  foreach my $gene(@$genes){
    my @transcripts = @{$gene->get_all_Transcripts};
    my @havana;
    my @ensembl;
    
    #print "number of transcript: ",scalar(@transcripts),"\n";
    # are there any havana transcripts?
    TRANSCRIPT:foreach my $transcript(@transcripts){
      
      foreach my $htranscript (@{$GB_HAVANA_INPUT_TRANSCRIPTTYPES}){
        if($transcript->biotype eq  $htranscript."_havana"){
          #print "I'm a havana transcript\n";
          push(@havana, $transcript);
          next TRANSCRIPT;
        }
      }
      push(@ensembl, $transcript);
    }
    if (!scalar(@havana)){
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
        
        #if ($delete_t = $self->are_matched_pair($ht, $et)){
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
        
        $self->set_transcript_relation($delete_trans, @t_pair);
        $t_pair[0]->biotype($MERGED_TRANSCRIPT_OUTPUT_TYPE);
        $t_pair[1]->biotype($MERGED_TRANSCRIPT_OUTPUT_TYPE);
        
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
# return 0 means keep both transcript as they have different coding region
# return 1 means keep both as they have same coding but different UTR exon structire
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

  # double check translation coords
  #print "HAVANA TRANS START: ",$havana->translation->genomic_start," END: ",$havana->translation->genomic_end,"\n";
  #print "ENSEMBL TRANS START: ",$ensembl->translation->genomic_start," END: ",$ensembl->translation->genomic_end,"\n";

  return 0 unless($havana->translation->genomic_start == $ensembl->translation->genomic_start &&
                  $havana->translation->genomic_end   == $ensembl->translation->genomic_end);

  # special case for single exon genes
  if(scalar(@hexons == 1)){
    #print "SINGLE EXONS!\n";
 
    if ($hexons[0]->start     == $eexons[0]->start &&
        $hexons[0]->end       == $eexons[0]->end &&
        $hexons[0]->strand    == $eexons[0]->strand &&
        $thexons[0]->coding_region_start($havana) == $teexons[0]->coding_region_start($ensembl) &&
        $thexons[0]->coding_region_end($havana) == $teexons[0]->coding_region_end($ensembl)
        ){
      # Both are exactly the same so we delete the Ensembl one unless the Ensembl one is already a merged one
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
 
   #  now check the rest of the internal exons that are not coding.
    for(my $i=1; $i<=($#hexons-1); $i++){
      return 1 unless ($hexons[$i]->start     == $eexons[$i]->start &&
                       $hexons[$i]->end       == $eexons[$i]->end &&
                       $hexons[$i]->strand    == $eexons[$i]->strand 
                       );
    }
    #print "INTERNAL UTR EXONS ARE OK \n";
    
    # Then check if the first an last exon are the same in both transcripts. If just start and end of UTR are different keep ensembl one
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
      
    }elsif (# CASE 3: BOTH ENSEMBL AND HAVANA HAVE UTR BUT WITH DIFFERENT START/END, KEEP ENSEMBL. Forward strand
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
      #print "MULTIEXON DELETE HAVANA\n";      
      return $havana;
      
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
      return $havana;
      
    }else{
      #print "Should I be here?\n";
      #print "Keep MULTIEXON BOTH\n";
      return 1;
      
    }
    
  }
  
  print " WEIRD CASE WE DID NOT THOUGHT ABOUT, CHECK RULES!\n";
  return 0;
  
}

sub add_ottt_xref{
 my($self, $ht) = @_;
 #TEST
 foreach my $entry(@{ $ht->get_all_DBEntries}){
   if ($entry->dbname eq 'Vega_transcript'){
     if($entry->primary_id eq $entry->display_id){
       
       print "I am adding an OTTT xref to the transcript\n";
       print "OTTT TO ADD: ",$entry->primary_id,"\n";
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
          
          print "OTTT xref to be added here\n";
          
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

    print "OTTT TO ADD: ",$t_pair[0]->stable_id,"\n";

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

  foreach my $exon (@{$transcript->get_all_Exons}){ 
    foreach my $sf (@{$exon->get_all_supporting_features}){
      $evidence{$sf->hseqname} = 1;
    }
  }

  foreach my $ev_key (keys %evidence){
    #print "Adding special attrib\n";
    my $attribute = Bio::EnsEMBL::Attribute->new
        (-CODE => 'otter_support',
         -NAME => 'otter support',
         -DESCRIPTION => 'Evidence ID that was used as supporting feature for building a gene in Vega',
         -VALUE => $ev_key);
    
    $trans_to_add_attrib->add_Attributes($attribute);
    
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
    
    # make different exon objects that are shared between transcripts 
    # ( regarding attributes: start, end, etc )
    # into unique exon objects 
    my $new_gene = $self->prune_Exons($gene);
    push( @pruned_genes, $new_gene );
  }
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
  my @genes;
  my $ensemblslice = $self->fetch_sequence($self->input_id, $self->ensembl_db);
  my $havanaslice = $self->fetch_sequence($self->input_id, $self->havana_db);
  print STDERR "Fetching ensembl genes\n";  
  EGENE: 
  foreach my $egene (@{$ensemblslice->get_all_Genes_by_type($GB_ENSEMBL_INPUT_GENETYPE)}){
    # Don't add those genes that contain only transcripts imported from HAVANA (this is important during a merge update)
    if ($egene->analysis->logic_name() eq $HAVANA_LOGIC_NAME){
      next EGENE;
    }else{
      push (@genes,$egene);
    } 
  }


  print STDERR "Retrieved ".scalar(@genes)." genes of type ".$GB_ENSEMBL_INPUT_GENETYPE."\n";
  print STDERR "Fetching havana genes\n";  
  my @hgenes = @{$havanaslice->get_all_Genes_by_type($GB_HAVANA_INPUT_GENETYPE)};
  print STDERR "Retrieved ".scalar(@hgenes)." genes of type ".$GB_HAVANA_INPUT_GENETYPE."\n";

  # We change the biotype of the havana genes/transcripts as it could happend to be the same as the ensembl ones
  foreach my $hgene(@hgenes){
    my $biotype = $hgene->biotype."_havana";
    $hgene->biotype($biotype);
    foreach my $htran (@{$hgene->get_all_Transcripts}) {
      my $tbiotype = $htran->biotype."_havana";
      $htran->biotype($tbiotype);
    }
  }

  push(@genes, @hgenes);

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
  
  print STDERR "Finished fetching genes\n";
  $self->combined_Transcripts(@transcripts);
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
        print "transcript found in discarded db\n";
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
        $tran_xref->dbname ne  "shares_CDS_with_OTTT" &&
        $tran_xref->dbname ne "shares_CDS_with_ENST" &&
        $tran_xref->dbname ne "OTTT"){
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
Description :   it clusters transcripts into genes according to exon overlap.
                It will take care of difficult cases like transcripts within introns.
                It also unify exons that are shared among transcripts.
    Returns :   a beautiful list of geen objects
    Args    :   a list of transcript objects

=cut

sub cluster_into_Genes{
  my ($self, @transcripts_unsorted) = @_;
  
  my $num_trans = scalar(@transcripts_unsorted);

  # First clean the coding exon cache in case it has any exons stored from previous called to the cluster_into_Genes function.
  $self->clear_coding_exons_cache;

  my @transcripts_unsorted_translation;

  foreach my $tran(@transcripts_unsorted){
    if ($tran->translation){
      push (@transcripts_unsorted_translation, $tran);
    }
  }

  my @transcripts = sort { $a->coding_region_start <=> $b->coding_region_start ? $a->coding_region_start <=> $b->coding_region_start  : $b->coding_region_end <=> $a->coding_region_end } @transcripts_unsorted_translation;
  my @clusters;

  # clusters transcripts by whether or not any coding exon overlaps with a coding exon in 
  # another transcript (came from original prune in GeneBuilder)
  foreach my $tran (@transcripts) {
  # First clean the coding exon cache in case it has any exons stored from previous called to the cluster_into_Genes function.
 # $self->clear_coding_exons_cache;

    my @matching_clusters;
  CLUSTER: 
    foreach my $cluster (@clusters) {
      
     # $self->clear_coding_exons_cache;

     #print "Transcript: ",$tran->stable_id," has coding region start: ",$tran->coding_region_start,"\n";

      foreach my $cluster_transcript (@$cluster) {
        if ($tran->coding_region_end  >= $cluster_transcript->coding_region_start &&
            $tran->coding_region_start <= $cluster_transcript->coding_region_end) {
          
          # foreach my $exon1 (@{$tran->get_all_Exons}) {
          # foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
          my $exons1 = get_coding_exons_for_transcript($tran);
          my $cluster_exons = get_coding_exons_for_transcript($cluster_transcript);

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
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my @merged_cluster;
      foreach my $clust (@matching_clusters) {
        push @merged_cluster, @$clust;
      }
      push @merged_cluster, $tran;
      push @new_clusters,\@merged_cluster;
      # Add back non matching clusters
      foreach my $clust (@clusters) {
        my $found = 0;
      MATCHING: 
	foreach my $m_clust (@matching_clusters) {
          if ($clust == $m_clust) {
            $found = 1;
            last MATCHING;
          }
        }
        if (!$found) {
          push @new_clusters,$clust;
        }
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
    my $gene = new Bio::EnsEMBL::Gene;
    foreach my $transcript (@$cluster){
      #print "Transcript Stable ID: ",$transcript->dbID,"\n";
      $gene->add_Transcript($transcript);
    }
    push( @genes, $gene );
  }
  return @genes;
}

############################################################

=head2 get_coding_exons_for_transcript

    Example :    my $exons1 = $self->get_coding_exons_for_gene($tran);
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
    my ($trans) = @_;

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
  if ($ntrans != $num_transcripts) {
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
	if ($exon == $tran->translation->start_Exon){
	  $tran->translation->start_Exon($found);
	}
	if ($exon == $tran->translation->end_Exon){
	  $tran->translation->end_Exon($found);
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

sub combined_Transcripts {
    my ($self,@transcripts) = @_;

    if (!defined($self->{_genewise_andthelike_transcripts})) {
        $self->{_genewise_andthelike_transcripts} = [];
    }

    if (scalar @transcripts > 0) {
	push(@{$self->{_genewise_andthelike_transcripts}},@transcripts);
    }

    return @{$self->{_genewise_andthelike_transcripts}};
}

############################################################

#=head2 my_genes
#
# Description: this holds and returns the genes that are produced after putting together genewise, combined and
#              processed_supporte_ab_initio predictions and removing the redundant set, giving priority
#              to long CDSs + UTR
#
#=cut
#
#
#sub my_genes {
#  my ($self,@genes) = @_;
#  
#  unless($self->{_my_genes}){
#    $self->{_my_genes} = [];
#  }
#
#  if (@genes){
#    push(@{$self->{_my_genes}},@genes);
#  }
#  return @{$self->{_my_genes}};
#}

############################################################

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

#=head2 predictions
#
# Description: get/set for the PredictionTranscripts. It is  set in new()
#
#=cut
#
#sub predictions {
#  my ($self,@predictions) = @_;
#
#  if(!$self->{_predictions}){
#    $self->{_predictions} = [];
#  }
#  if ( @predictions ) {
#     push(@{$self->{_predictions}},@predictions);
#  }
#  return @{$self->{_predictions}};
#}

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
