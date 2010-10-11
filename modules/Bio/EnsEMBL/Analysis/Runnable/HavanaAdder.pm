
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
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Analysis::Config::HavanaAdder qw (
  ENSEMBL_INPUT_CODING_TYPE
  HAVANA_INPUT_CODING_TYPE
  ENSEMBL_INPUT_PROCESSED_TYPE
  HAVANA_INPUT_PROCESSED_TYPE
  ENSEMBL_INPUT_PSEUDO_TYPE
  HAVANA_INPUT_PSEUDO_TYPE
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
use warnings;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Root);


############################################################

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ( $slice, $input_id ) = $self->_rearrange( [qw(SLICE INPUT_ID)], @args );

  $self->throw("Must input a slice to HavanaAdder") unless defined($slice);
  $self->{_final_genes} = [];
  $self->{_gene_types}  = [];

  $self->query($slice);
  $self->gene_types($ENSEMBL_INPUT_CODING_TYPE);
  $self->gene_types($HAVANA_INPUT_CODING_TYPE);

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
  my ( $self, $id ) = @_;

  if ( defined($id) ) {
    $self->{_input_id} = $id;
  }
  return $self->{_input_id};
}

############################################################

=head2 build_Genes

 Example    : my @genes = $self->build_Genes
 Description: Builds genes. It is like the run method in Runnables. It
              calls everything that needs to be done.
 Returns    : none
 Args       : none
 Caller     : Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder

=cut

sub build_Genes {
  my ($self) = @_;

  print STDERR "Building genes...\n";

  # get all genes of type defined in gene_types() on this slice
  $self->get_Genes;

  my @ensembl_coding_transcripts    = @{ $self->combined_Transcripts };
  my @ensembl_processed_transcripts = @{ $self->combined_Processed_Transcripts };
  my @ensembl_pseudo_transcripts    = @{ $self->combined_PseudoTranscripts };

  my @havana_coding_genes               = @{ $self->havana_Coding_Genes };
  my @havana_processed_transcript_genes = @{ $self->havana_Processed_Transcript_Genes };
  my @havana_pseudogene_Genes           = @{ $self->havana_Pseudogene_Genes };

  # do a preliminary clustering
  my @preliminary_coding_genes = @{
    $self->cluster_into_Genes( \@ensembl_coding_transcripts,
                               \@havana_coding_genes, 1 )};

  print "\nCoding gene clusters: ", scalar(@preliminary_coding_genes), "\n";

  my @preliminary_processed_genes = @{
    $self->cluster_into_Genes( \@ensembl_processed_transcripts,
                               \@havana_processed_transcript_genes ) };

  print "\nProcessed transcript genes clusters: ",
    scalar(@preliminary_processed_genes), "\n";

  my @preliminary_pseudo_genes = @{
    $self->cluster_into_Genes( \@ensembl_pseudo_transcripts,
                               \@havana_pseudogene_Genes ) };

  print "\nPseudogene clusters: ", scalar(@preliminary_pseudo_genes), "\n";

  push( @preliminary_pseudo_genes, @preliminary_processed_genes );

  my @clustered_gene_set = @{
    $self->combine_gene_clusters( \@preliminary_coding_genes,
                                  \@preliminary_pseudo_genes ) };

  print "\nTotal clusters: ", scalar(@clustered_gene_set), "\n";

  # merge redundant ensembl transcripts which match a havana one
  $self->_merge_redundant_transcripts( \@clustered_gene_set );

  # make shared exons unique objects
  my @genes = @{ $self->_make_shared_exons_unique( \@clustered_gene_set ) };

  print STDERR scalar(@genes) . " genes built\n";

  $self->update_biotypes(@genes);

  $self->final_genes( \@genes );
} ## end sub build_Genes


sub _merge_redundant_transcripts {
  my ( $self, $genes ) = @_;

  print "\nNumber of genes in clusters: ", scalar( @{$genes} ), "\n";

GENE:
  foreach my $gene ( @{$genes} ) {
    my @transcripts = @{ $gene->get_all_Transcripts };
    my @havana;
    my @ensembl;

    #print "Number of transcript: ", scalar(@transcripts),"\n";

    # Classify transcripts according to their origin (Havana or Ensembl)
  TRANSCRIPT: foreach my $transcript (@transcripts) {

      if ( $transcript->biotype =~ /_hav/ ) {
        #print "I'm a havana transcript with biotype: ",$transcript->biotype .  "\t" . $transcript->dbID . "\n";
        push( @havana, $transcript );
        next TRANSCRIPT;
      }

      print "\t\tin _merge_redundant_transcripts: "
        . "adding this transcript to ensembl array: "
        . $transcript->dbID . " ("
        . $transcript->biotype . ")\n";

      push( @ensembl, $transcript );
    }
    if ( !scalar(@havana) ) {
      next GENE;
    }

    if ( !scalar(@ensembl) ) {
      # As this is a havana only gene I still want to add an OTTT
      # entry to the transcripts and an OTTG to the gene.
      $self->add_ottg_xref( $gene, $gene->stable_id );

      foreach my $hav_trans (@havana) {
        $self->add_ottt_xref($hav_trans);
      }

      next GENE;
    }
    print "Havana tran: ",  scalar(@havana),  "\n";
    print "Ensembl tran: ", scalar(@ensembl), "\n";

    # Compare each havana transcript to each ensembl one
    foreach my $ht (@havana) {
      print "looking at havana trans: " . $ht->dbID . "\n";

      # We add an attribute to the havana transcripts that shows which
      # supporting features have been used to build the Exons. This
      # will allow us to distinguish them when it comes to the web
      # display.
      print "\nAdding attributes to havana transcript " . $ht->dbID . "\n";
      $self->add_havana_attribute( $ht, $ht );

      #print "Deleting havana transcript supporting features\n";

      # Original supporting features in havana transcripts are removed
      # as Havana build transcripts based on exon supporting features
      # so the transcripts supporting features are artificial features
      # that are not accurate.
      $ht->flush_supporting_features;

      print "\nNumber of ensembl transcripts: ", scalar(@ensembl),"\n";
      my $delete_trans = 0;
      my @t_pair;

      foreach my $et (@ensembl) {
        print "looking at ensembl trans: " . $et->dbID . "\n";

        my $delete_t = $self->are_matched_pair( $ht, $et );

        # We check all possible match pairs and give preference to the
        # one that shares CDS and UTR. This was added so only the best
        # matching havana/ensembl pair is chosen and to avoid a one
        # to many link

        # TODO: have a closer look at the code here.
        if ($delete_t) {
          if ( $delete_t == $et ) {
            print "--->>> I'm deleting matching transcript: " .  $delete_t->dbID . "\n";
            $delete_trans = $delete_t;
            @t_pair = ( $ht, $et );
          } elsif ( $delete_trans != $et && $delete_t == $ht ) {
            print "deleted_trans != et but == ht " . $delete_t->dbID . "\n";
            $delete_trans = $delete_t;
            @t_pair = ( $ht, $et );
          } elsif ( $delete_trans == 0 ) {
            $delete_trans = $delete_t;
            print "the other case, keeping both? " . $et->dbID . " and " .  $ht->dbID . "\n";
            @t_pair = ( $ht, $et );
          }
        }
      }

      if ( $delete_trans && $delete_trans != 0 ) {
        my $new_bt_0;    #biotype
        my $new_bt_1;    #biotype
        $self->set_transcript_relation( $delete_trans, @t_pair );

        # We flag the transcript biotype to show that one of the
        # transcripts is a merged one.
        unless ( $t_pair[0]->biotype =~ /$MERGED_TRANSCRIPT_OUTPUT_TYPE/ ) {
          $new_bt_0 = $t_pair[0]->biotype . $MERGED_TRANSCRIPT_OUTPUT_TYPE;
        }
        unless ( $t_pair[1]->biotype =~ /$MERGED_TRANSCRIPT_OUTPUT_TYPE/ ) {
          $new_bt_1 = $t_pair[1]->biotype . $MERGED_TRANSCRIPT_OUTPUT_TYPE;
        }
        $t_pair[0]->biotype($new_bt_0);
        $t_pair[1]->biotype($new_bt_1);

        # We want to remove the redundant transcript unless both share
        # CDS but have different UTR structure as in that case we
        # annotate both transcripts.
        $self->_remove_transcript_from_gene( $gene, $delete_trans )
          unless $delete_trans == 1;

      } else {
        $self->add_ottt_xref($ht);
      }
    } ## end foreach my $ht (@havana)

    $gene->recalculate_coordinates;

  } ## end foreach my $gene (@$genes)
} ## end sub _merge_redundant_transcripts

# are_matched_pair check return 4 different possible values:
# return 0 means keep both transcripts as they have different coding region
#          or different exon structure
# return 1 means keep both as they have same coding but different UTR exon structure
# return ($ensembl) means keep havana transcript and remove ensembl 
# return ($havana) means keep ensembl transcript and remove hanana

sub are_matched_pair {
  my ( $self, $havana, $ensembl ) = @_;

  # Fetch all exons in each transcript
  my @hexons = @{ $havana->get_all_Exons };
  my @eexons = @{ $ensembl->get_all_Exons };

  my $non_coding_h = 0;

  # Keep track if the havana transcript is coding or not via its exons
  my @thexons = @{ $havana->get_all_translateable_Exons };
  if ( !@thexons || @thexons == 0 ) {
    $non_coding_h = 1;
  }

  my $non_coding_e   = 0;
  my $coding_non_met = 0;
  my @teexons        = [];

  print "____________________________________\n";
  print "HAVANA ID: ",$havana->dbID, " ENSEMBL: ",$ensembl->dbID,"\n";
  print "HAVANA ID: ",$havana->stable_id, " ENSEMBL: ",$ensembl->stable_id,"\n";

  print "\nEnsembl trans dbID: " . $ensembl->dbID . "\tbiotype: " . $ensembl->biotype. "\n";
  print "no. of exons: " . scalar(@{ $ensembl->get_all_Exons }) . "\n";
  print "no. of coding exons: " . scalar(@{ $ensembl->get_all_translateable_Exons }) . "\n";

  print "\nHavana trans dbID: " . $havana->dbID . "\tbiotype: " .  $havana->biotype . "\n";
  print "no. of exons: " . scalar(@{ $havana->get_all_Exons }) . "\n";
  print "no. of coding exons: " . scalar(@{ $havana->get_all_translateable_Exons }) . "\n";

  my $complete_5        = 0;
  my $incomplete        = 0;
  my $full_length_cdnas = 0;

  if ( @{ $ensembl->get_all_translateable_Exons } ) {

    my $cdna         = $ensembl->spliced_seq;
    my $coding_start = $ensembl->cdna_coding_start;
    my $coding_end   = $ensembl->cdna_coding_end;
    my $cdna_start   = substr( $cdna, $coding_start - 1, 3 );
    my $cdna_end     = substr( $cdna, $coding_end - 3, 3 );

    if ( $cdna_start eq 'ATG' ) {
      print "The ensembl transcript start with ATG\n";

      if ( $cdna_end eq 'TAG' || $cdna_end eq 'TGA' || $cdna_end eq 'TAA' ) {
        print "Ensembl transcript has start and stop codon: " .  $ensembl->dbID . "\n";
        $full_length_cdnas++;
        @teexons = @{ $ensembl->get_all_translateable_Exons };
      } else {
        # Flag the transcripts with start but no end for later deletion.

        $non_coding_e = 2;
        $complete_5++;
        #print "Have a MET but not a stop:  " .  $ensembl->dbID . "\n";
      }
    } else {
      $non_coding_e = 2;
      $incomplete++;
      #print "Not full-length transcript to delete: " . $ensembl->dbID . "\n";
    }
  } else {
    $non_coding_e = 1;
    print "\nNon-coding ensembl trans id " . $ensembl->dbID . "\n";
  }

  # Check that the number of exons is the same in both transcripts
  if (scalar(@hexons) != scalar(@eexons)) {
    print "et and ht don't have same number of exons, skipping the comparison\n";
    return 0;
  }
  #return 0 unless scalar(@hexons) == scalar(@eexons);

  if ( $non_coding_h == 1 && $non_coding_e == 1 ) {
    print "\n===>>> ensembl and havana are both non-coding <<<===\n";
    # We check two non coding transcripts. If they have the same
    # structure we keep the one from havana but if the have same exon
    # strucutre but one is slightly longer we take the longest one of
    # the pair.

    if ( !$self->check_internal_exon_structure( \@eexons, \@hexons ) ) {
      print "value is: " . $self->check_internal_exon_structure( \@eexons, \@hexons ) . "\n";
      print "NON-CODING CASE 0 - BOTH selected\n";

      # CASE 0: the two transcripts have different internal exon structure
      return 0;
    }
    # CASE 1: Havana is longer or both are exactly the same
    print "NON-CODING CASE 1 - Havana longer or both are the same length\n";
    if ( $self->check_terminal_exon_structure( \@hexons, \@eexons ) ) {
      print "value is: " . $self->check_internal_exon_structure( \@eexons, \@hexons ) . "\n";
      print "We keep Havana - havana longer\n";
      return $ensembl;
    } else {
      # CASE 2: EnsEMBL is longer than Havana
      print "NON-CODING CASE 2: Ensembl is longer but we keep Havana\n";
      return $ensembl;
    }
  } elsif ( $non_coding_h != $non_coding_e ) {
    if ( $non_coding_e == 2 ) {
      print "===>>> incomplete ensembl models, will not be merged\n";
      # Any transcript coming here should be deleted as we don't want to have
      # to have them in the merge set.

      print "Not full_length ensembl transcript, will be deleted: " . $ensembl->dbID . "\n";
      print "\tFull-length: $full_length_cdnas\n"
        . "\t5' complete: $complete_5\n"
        . "\tIncomplete: $incomplete\n";

      return $ensembl;
    }

    print "\n===>>> coding and non-coding overlap <<<===\n";
    # This is a case of a pseudogene overlapping a coding gene so by
    # now we keep both.

    print "Warning Pseudogene and coding overlap for"
        . " HAVANA ID: " . $havana->dbID
        . " ENSEMBL: "   . $ensembl->dbID
        . "\n";

    # If Havana is coding and ensembl is non coding
    if ( $non_coding_h == 0 ) {
      print "\n===>>> ensembl non-coding and havana coding <<<===\n";

      unless ( ( scalar(@eexons) == scalar(@thexons) )
              && $self->check_internal_exon_structure( \@eexons, \@thexons ) )
      {
        # CASE 0: the two transcripts have different internal exon structure
        print "HAV CODING vs ENS NON-CODING CASE 0 - "
            . "different internal exon structure, both selected\n";

        print "value is: " . $self->check_internal_exon_structure( \@eexons, \@hexons ) . "\n";
        return 0;
      }
      # Here we may consider removing the Ensembl transcript and keep
      # only the havana coding???

      # CASE 1: If the internal structure is the same we then keep the
      #         Havana one.
      print "HAV CODING vs ENS NON-CODING CASE 1 - internal structure the same, keep havana\n";

      print "deleting ensembl " . $ensembl->biotype . "\n";
      return $ensembl;
    }

    # If ensembl is coding and havana is non coding, the ensembl
    # transcript is turned to non-coding.
    if ( $non_coding_e == 0 ) {

      print "===>>> ensembl is coding, havana is non-coding <<<===\n";

      # Keep both if they don't have the same coding exon structure
      return 0 unless scalar(@hexons) == scalar(@eexons);

      # Check that the ensembl transcript is not a CCDS model, if so, don't
      # make it to non-coding and keep both transcripts.

      if ( $self->check_transcript_in_external_db('ccds', $ensembl ) == 0 ) {
        print "ENS CODING vs HAV NON-CODING CASE 2 - keeping ensembl regardless\n";
        print "The Ensembl model " . $ensembl->dbID . " is a CCDS transcript, "
            . "will not convert it to non-coding.\n";
        return 0;
      }

      # May want to remove the translation and keep ensembl?
      #
      # First check if the internal structure of the whole
      # transcripts is conserved.
      if ( scalar(@hexons) > 1 ) {
        unless ( $self->check_internal_exon_structure( \@hexons, \@eexons ) )
        {
          # CASE 0: The two transcripts have different internal exon
          #         structure we keep both
          print "ENS CODING vs HAV NON-COSDING CASE 3 - "
              . "different internal exon structure, both selected\n";
          return 0;
        }
      }
      # Now we check if the coding bit is longer than the whole non
      # coding transcript.

      # CASE 1: The ensembl transcript is longer or equal so we remove
      #         the havana transcript. This is not true any more, we
      #         keep the havana transcript regardless of length. Since the
      #         ensembl model is deleted, there's no need to change its
      #         biotype or remove its translation.
      if ( $self->check_terminal_exon_structure( \@hexons, \@teexons ) ) {
        print "ENS CODING vs HAV NON-COSDING CASE 4 - "
            . "removing ensembl model as it is equal or longer "
            . "than the Havana model but Havana models are prioritised.\n";

        #print "making ensembl model " . $ensembl->dbID . " to non-coding\n";

        # Do we need to do this? NO!
        #$ensembl->{translation} = undef;
        #$ensembl->biotype($havana->biotype."_e");
        return $ensembl;

      } else {
        # CASE 2: The havana transcripts is longer in both ends so we
        #         remove the ensembl transcript. Again, there's no need to
        #         change the biotype of the model or remove its translation.

        print "ENS CODING vs HAV NON-COSDING CASE 5 - "
            . "havana transcript is longer in both ends, "
            . "removing ensembl: " . $ensembl->dbID . "\n";
        #print "changing the biotype from: " $ensembl->biotype . " to " .  $havana->biotype . "_e\n";

        #$ensembl->{_translation_array} = [];
        #$ensembl->biotype($havana->biotype."_e");
        return $ensembl;
      }
    } ## end if ( $non_coding_e == ...

  } elsif ( $non_coding_h == 0 && $non_coding_e == 0 ) {

    return 0
      unless ( ( scalar(@teexons) == scalar(@thexons) )
               && $havana->translation->genomic_start == $ensembl->translation->genomic_start
               && $havana->translation->genomic_end == $ensembl->translation->genomic_end );

    # Special case for single exon genes
    if ( scalar(@hexons) == 1 ) {
      #print "SINGLE EXONS!\n";

      if (    $hexons[0]->start == $eexons[0]->start
           && $hexons[0]->end == $eexons[0]->end
           && $hexons[0]->strand == $eexons[0]->strand
           && $thexons[0]->coding_region_start($havana) == $teexons[0]->coding_region_start($ensembl)
           && $thexons[0]->coding_region_end($havana) == $teexons[0]->coding_region_end($ensembl) )
      {
        # Both are exactly the same so we delete the Ensembl one
        # unless the Ensembl one is already a merged one <-???? IS not doing
        # this here.
        return $ensembl;

      } elsif ($hexons[0]->start <= $eexons[0]->start
            && $hexons[0]->end >= $eexons[0]->end
            && $hexons[0]->strand == $eexons[0]->strand
            && $eexons[0]->start == $teexons[0]->coding_region_start($ensembl)
            && $eexons[0]->end == $teexons[0]->coding_region_end($ensembl) )
      {
        # Ensembl gene don't have UTR and Havana has then delete Ensembl one
        return $ensembl;

      }
      elsif ( ( (
                 $hexons[0]->start != $eexons[0]->start
              || $hexons[0]->end != $eexons[0]->end )
            && $hexons[0]->strand == $eexons[0]->strand )
          && ( $eexons[0]->start != $teexons[0]->coding_region_start($ensembl)
              || $eexons[0]->end != $teexons[0]->coding_region_end($ensembl) )
        )
      {
        # Both contain UTR keep ENSEMBL
        return $havana;

      } else {
        # We can be here when genes have different UTR start/end and
        # different CDS start/end or when the UTR start/end is the
        # same but the CDS start/end is different.

        #print "Keep both single exon genes\n";
        return 0;

      }
    } ## end if ( scalar(@hexons) ==...

    # if is a multi exons transcript
    else {
      # First we check the internal coding structure of the transcript
      # where everything has to be exactly equal

      #print "CHECKING INTERNAL EXONS \n";
      for ( my $i = 1 ; $i <= ( $#thexons - 1 ) ; $i++ ) {
        return 0
          unless (    $thexons[$i]->start == $teexons[$i]->start
                   && $thexons[$i]->end == $teexons[$i]->end
                   && $thexons[$i]->strand == $teexons[$i]->strand );
      }
      #print "INTERNAL CODING EXONS ARE OK \n";

      # Now check the rest of the internal exons that are not coding.
      # This is to find if the UTR exon structure is the same
      for ( my $i = 1 ; $i <= ( $#hexons - 1 ) ; $i++ ) {
        if (    $hexons[$i]->start != $eexons[$i]->start
             || $hexons[$i]->end != $eexons[$i]->end
             || $hexons[$i]->strand != $eexons[$i]->strand )
        {
          #print "CASE WITH SAME CDS BUT DIFFERENT UTR STRUCTURE\n";
          #print "HAVANA DIFF UTR BOUNDARIES: ".$havana->seq_region_name." - ".$havana->seq_region_start. " - ".$havana->seq_region_end."\n";
          #print "ENSEMBL DIFF UTR BOUNDARIES: ".$ensembl->seq_region_name." - ".$ensembl->seq_region_start. " - ".$ensembl->seq_region_end."\n";
          return 1;
        }
      }
      #print "INTERNAL UTR EXONS ARE OK \n";

      # Then check if the first an last exon are the same in both
      # transcripts. If just start and end of UTR are different keep havana
      # one.
      #
      # CASE 1: Both coding and UTR are the same, keep Havana and delete Ensembl
      if (    $hexons[0]->start == $eexons[0]->start
           && $hexons[0]->end == $eexons[0]->end
           && $hexons[0]->strand == $eexons[0]->strand
           && $hexons[-1]->start == $eexons[-1]->start
           && $hexons[-1]->end == $eexons[-1]->end
           && $hexons[-1]->strand == $eexons[-1]->strand )
      {
        #print "MULTIEXON DELETE ENSEMBL\n";
        return $ensembl;

      } elsif
        ( #CASE 2": HAVANA HAS UTR AND ENSEMBL DOESNT, KEEP HAVANA. Forward strand
           $hexons[0]->strand == 1
        && $hexons[0]->end == $eexons[0]->end
        && $hexons[0]->strand == $eexons[0]->strand
        && $hexons[-1]->start == $eexons[-1]->start
        && $hexons[-1]->strand == $eexons[-1]->strand
        && $eexons[0]->start == $teexons[0]->coding_region_start($ensembl)
        && $eexons[-1]->end == $teexons[-1]->coding_region_end($ensembl)
        && (    $hexons[-1]->end != $eexons[-1]->end
             || $hexons[0]->start != $eexons[0]->start ) )
      {
        #print "MULTIEXON DELETE ENSEMBL\n";
        return $ensembl;

      } elsif
        ( # CASE 3: BOTH ENSEMBL AND HAVANA HAVE UTR BUT WITH DIFFERENT START/END, KEEP Havana. Forward strand
           $hexons[0]->strand == 1
        && $hexons[0]->end == $eexons[0]->end
        && $hexons[0]->strand == $eexons[0]->strand
        && $hexons[-1]->start == $eexons[-1]->start
        && $hexons[-1]->strand == $eexons[-1]->strand
        && (   $eexons[0]->start != $teexons[0]->coding_region_start($ensembl)
            || $eexons[-1]->end != $teexons[-1]->coding_region_end($ensembl) )
        && (    $hexons[-1]->end != $eexons[-1]->end
             || $hexons[0]->start != $eexons[0]->start ) )
      {
        #print "MULTIEXON DELETE ENSEMBL\n";
        #print "CASE WITH DIFFERENT TERMINAL EXON BOUNDARIES\n";
        #print "HAVANA BOUNDARIES: ".$havana->seq_region_name." - ".$havana->seq_region_start. " - ".$havana->seq_region_end."\n";
        #print "ENSEMBL BOUNDARIES: ".$ensembl->seq_region_name." - ".$ensembl->seq_region_start. " - ".$ensembl->seq_region_end."\n";
        return $ensembl;

      } elsif (    # CASE 4: Same as case 2 but in reverse strand
           $hexons[0]->strand == -1
        && $hexons[0]->start == $eexons[0]->start
        && $hexons[0]->strand == $eexons[0]->strand
        && $hexons[-1]->end == $eexons[-1]->end
        && $hexons[-1]->strand == $eexons[-1]->strand
        && $eexons[-1]->start == $teexons[-1]->coding_region_start($ensembl)
        && $eexons[0]->end == $teexons[0]->coding_region_end($ensembl)
        && (    $hexons[0]->end != $eexons[0]->end
             || $hexons[-1]->start != $eexons[-1]->start )

        )
      {
        #print "MULTIEXON DELETE ENSEMBL\n";
        return $ensembl;

      } elsif (    # CASE 5: Same as case 3 but in reverse strand
           $hexons[0]->strand == -1
        && $hexons[0]->start == $eexons[0]->start
        && $hexons[0]->strand == $eexons[0]->strand
        && $hexons[-1]->end == $eexons[-1]->end
        && $hexons[-1]->strand == $eexons[-1]->strand
        && ( $eexons[-1]->start != $teexons[-1]->coding_region_start($ensembl)
             || $eexons[0]->end != $teexons[0]->coding_region_end($ensembl) )
        && (    $hexons[0]->end != $eexons[0]->end
             || $hexons[-1]->start != $eexons[-1]->start )

        )
      {
        #print "MULTIEXON DELETE HAVANA\n";
        #print "CASE WITH DIFFERENT TERMINAL EXON BOUNDARIES\n";
        #print "HAVANA BOUNDARIES: ".$havana->seq_region_name." - ".$havana->seq_region_start. " - ".$havana->seq_region_end."\n";
        #print "ENSEMBL BOUNDARIES: ".$ensembl->seq_region_name." - ".$ensembl->seq_region_start. " - ".$ensembl->seq_region_end."\n";
        return $ensembl;

      } else {
        print "\n\tShould I be here?\n";
        print "\tKeep MULTIEXON BOTH: havana "
          . $havana->dbID
          . " and ensembl: "
          . $ensembl->dbID . "\n";
        print "\tHAVANA BOUNDARIES: "
          . $havana->seq_region_name . " - "
          . $havana->seq_region_start . " - "
          . $havana->seq_region_end . "\n";
        print "\tENSEMBL BOUNDARIES: "
          . $ensembl->seq_region_name . " - "
          . $ensembl->seq_region_start . " - "
          . $ensembl->seq_region_end . "\n";
        print "\n";
        return 1;
      }

    } ## end else [ if ( scalar(@hexons) ==...

    print "WEIRD CASE WE DID NOT THINK ABOUT, CHECK RULES!\n";
    return 0;
  } ## end elsif ( $non_coding_h == ...
} ## end sub are_matched_pair

=head2 check_internal_exon_structure

  Description: Check if the start and end of internal exon pairs in two sets of exons is the same
  Return     : Returns 0 if they are the same, returns 1 if they are differents

=cut

sub check_internal_exon_structure {

  my ( $self, $firstexons, $secondexons ) = @_;

  my @exons1 = @{$firstexons};
  my @exons2 = @{$secondexons};

  if ( scalar(@exons1) != scalar(@exons2) ) {
    print "NOOOOOOOOO WHAT HAPPENED HERE?\n";
    print scalar(@exons1), "  vs. ", scalar(@exons2), "\n";
  }

  # We check if the transcript has more than two exon as otherwise
  # we will be checking only the coordinates in the last/second exon
  # which may produce wrong results.
  if ( scalar(@exons1) > 2 ) {
    for ( my $i = 1 ; $i <= ( $#exons1 - 1 ) ; $i++ ) {
      return 0
        unless (    $exons1[$i]->start == $exons2[$i]->start
                 && $exons1[$i]->end == $exons2[$i]->end
                 && $exons1[$i]->strand == $exons2[$i]->strand
                 && ( (    $exons1[0]->strand == 1
                        && $exons1[0]->end == $exons2[0]->end
                        && $exons1[-1]->start == $exons2[-1]->start )
                      || (    $exons1[0]->strand == -1
                           && $exons1[0]->start == $exons2[0]->start
                           && $exons1[-1]->end == $exons2[-1]->end ) ) );
    }
    return 1;
  } else {
    return 0
      unless ( (    $exons1[0]->strand == 1
                 && $exons1[0]->end == $exons2[0]->end
                 && $exons1[-1]->start == $exons2[-1]->start )
               || (    $exons1[0]->strand == -1
                    && $exons1[0]->start == $exons2[0]->start
                    && $exons1[-1]->end == $exons2[-1]->end ) );
  }
  return 1;
} ## end sub check_internal_exon_structure

=head2 check_terminal_exon_structure

  Description: Checks if beginning and end of transcripts
               coincide by looking at the terminal exon coords.
  Return     : Returns 0 is the first exon set is shorter that the second in
               both start and end otherwise returns 1.
=cut

sub check_terminal_exon_structure {

  my ( $self, $firstexons, $secondexons ) = @_;

  my @exons1 = @{$firstexons};
  my @exons2 = @{$secondexons};

  # I added the following or check "|| $exons1[0] eq $exons1[-1]"
  # to handle single exon genes more efficiently.

  if (    ( $exons1[0]->strand == 1 || $exons1[0] eq $exons1[-1] )
       && $exons1[0]->strand == $exons2[0]->strand
       && $exons1[0]->start <= $exons2[0]->start
       && $exons1[-1]->end >= $exons2[-1]->end )
  {
    return 1;
  } elsif (    $exons1[0]->strand == -1
            && $exons1[0]->end <= $exons2[0]->end
            && $exons1[-1]->start >= $exons2[-1]->start
            && $exons1[0]->strand == $exons2[0]->strand )
  {
    return 1;
  }
  # CASE 2: EnsEMBL is longer than Havana
  elsif (    ( $exons1[0]->strand == 1 || $exons1[0] eq $exons1[-1] )
          && $exons1[0]->start >= $exons2[0]->start
          && $exons1[-1]->end <= $exons2[-1]->end
          && $exons1[0]->strand == $exons2[0]->strand )
  {
    return 0;
  } elsif (    $exons1[0]->strand == -1
            && $exons1[0]->end >= $exons2[0]->end
            && $exons1[-1]->start <= $exons2[-1]->start
            && $exons1[0]->strand == $exons2[0]->strand )
  {
    return 0;
  } else {
    return 1;
  }
} ## end sub check_terminal_exon_structure

sub add_ottt_xref {
  my ( $self, $ht ) = @_;

  foreach my $entry ( @{ $ht->get_all_DBEntries } ) {
    if ( $entry->dbname eq 'Vega_transcript' ) {
      if ( $entry->primary_id eq $entry->display_id ) {

        #print "I am adding an OTTT xref to the transcript\n";
        #print "OTTT TO ADD: ",$entry->primary_id,"\n";
        my $xref_ottt =
          new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                     -display_id    => $ht->display_id,
                                     -priority      => 1,
                                     -xref_priority => 0,
                                     -version       => 1,
                                     -release       => 1,
                                     -dbname        => 'OTTT' );

        $xref_ottt->status("XREF");
        $ht->add_DBEntry($xref_ottt);
      }
    }
  }
} ## end sub add_ottt_xref

sub add_ottg_xref {
  my ( $self, $hg, $ottg ) = @_;

  #print "I am adding an OTTG xref to the transcript or Gene\n";

  # Don't want to save the OTTG xref twice so check if it's already stored.
  foreach my $entry ( @{ $hg->get_all_DBEntries } ) {
    if ( $entry->dbname eq 'OTTG' && ( $entry->primary_id eq $ottg ) ) {
      return 0;
    }
  }
  #print "OTTG TO ADD: ", $ottg, "\n";
  my $xref_ottg =
    new Bio::EnsEMBL::DBEntry( -primary_id    => $ottg,
                               -display_id    => $ottg,
                               -priority      => 1,
                               -xref_priority => 0,
                               -version       => 1,
                               -release       => 1,
                               -dbname        => 'OTTG' );

  $xref_ottg->status("XREF");

  $hg->add_DBEntry($xref_ottg);
} ## end sub add_ottg_xref


sub set_transcript_relation {
  # $t_pair[0] is the havana transcript and $t_pair[1] is the ensembl transcript
  my ( $self, $delete_t, @t_pair ) = @_;

  # If both share CDS and UTR is different in structure and number of
  # exons we still keep both, and we link them via Xref.
  if ( $delete_t == 1 ) {

    # transfer OTTT ID and/or ENST
    foreach my $entry ( @{ $t_pair[0]->get_all_DBEntries } ) {
      if ( $entry->dbname eq 'Vega_transcript' ) {
        if ( $entry->primary_id eq $entry->display_id ) {

          my $newentry =
            new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                       -display_id    => $entry->display_id,
                                       -priority      => 1,
                                       -xref_priority => 0,
                                       -version       => 1,
                                       -release       => 1,
                                       -dbname => 'shares_CDS_with_OTTT' );

          $newentry->status("XREF");

          $t_pair[1]->add_DBEntry($newentry);

          #TEST
          my $xref_ottt =
            new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                       -display_id    => $entry->display_id,
                                       -priority      => 1,
                                       -xref_priority => 0,
                                       -version       => 1,
                                       -release       => 1,
                                       -dbname        => 'OTTT' );

          #print "OTTT xref to be added here\n";

          $xref_ottt->status("XREF");

          $t_pair[0]->add_DBEntry($xref_ottt);
          #END of TEST

        } ## end if ( $entry->primary_id...
      } ## end if ( $entry->dbname eq...
    } ## end foreach my $entry ( @{ $t_pair...

    my $link_attrib =
      Bio::EnsEMBL::Attribute->new(
      -CODE => 'enst_link',
      -NAME => 'enst link',
      -DESCRIPTION => 'Code to link a OTTT with an ENST when they both share the CDS of ENST',
      -VALUE => $t_pair[1]->dbID );

    $t_pair[1]->add_Attributes($link_attrib);

    my $xref_entry =
      new Bio::EnsEMBL::DBEntry( -primary_id    => $t_pair[1]->dbID,
                                 -display_id    => $t_pair[1]->dbID,
                                 -priority      => 1,
                                 -xref_priority => 0,
                                 -version       => 1,
                                 -release       => 1,
                                 -dbname        => 'shares_CDS_with_ENST' );

    $xref_entry->status("XREF");

    $t_pair[0]->add_DBEntry($xref_entry);

    #print "OTTT TO ADD: ",$t_pair[0]->stable_id,"\n";

  } ## end if ( $delete_t == 1 )

  # If transcript to delete is havana we create an xref for the entry
  # saying that the transcript is CDS equal to ensembl.
  elsif ( $delete_t == $t_pair[0] ) {
    # transfer OTTT ID and/or ENST
    foreach my $entry ( @{ $t_pair[0]->get_all_DBEntries } ) {
      if ( $entry->dbname eq 'Vega_transcript' ) {
        if ( $entry->primary_id eq $entry->display_id ) {
          my $newentry =
            new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                       -display_id    => $entry->display_id,
                                       -priority      => 1,
                                       -xref_priority => 0,
                                       -version       => 1,
                                       -release       => 1,
                                       -dbname => 'shares_CDS_with_OTTT' );

          $newentry->status("XREF");

          $t_pair[1]->add_DBEntry($newentry);
        }
      }
    }

    # We add a transcript attribute to the ensembl transcript with the
    # start and end coords of the Havana transcript that we will delete
    my $attrib_value =
        $t_pair[0]->slice->coord_system_name . ":"
      . $t_pair[0]->slice->coord_system->version . ":"
      . $t_pair[0]->slice->seq_region_name . ":"
      . $t_pair[0]->start . ":"
      . $t_pair[0]->end . ":1";
    print "\tATTRIB VALUE:---------- ", $attrib_value, "\n";
    my $attribute =
      Bio::EnsEMBL::Attribute->new( -CODE        => 'TranscriptEdge',
                                    -NAME        => 'Transcript Edge',
                                    -DESCRIPTION => '',
                                    -VALUE       => $attrib_value );

    $t_pair[1]->add_Attributes($attribute);

    # When we delete a Havana transcript we want to transfer the exon
    # supporting features to the transcript we keep
    my @delete_e = @{ $delete_t->get_all_Exons };
    my @exons    = @{ $t_pair[1]->get_all_Exons };

    if ( scalar(@delete_e) == scalar(@exons) ) {

      my $e;
      for ( $e = 0, $e < scalar(@delete_e), $e++ ) {
        if ( $delete_e[$e] && $exons[$e] ) {
          $self->transfer_supporting_evidence( $delete_e[$e], $exons[$e] );
        }
      }
    }

    # We add attributes to the havana transcript showing which
    # supporting features where used for the transcript in Havana
    $self->add_havana_attribute( $t_pair[0], $t_pair[1] );

  } elsif ( $delete_t == $t_pair[1] ) {
    # If the transcript to delete is ENSEMBL we add an xref entry say
    # that both transcripts are exact matches (including UTR)
    foreach my $entry ( @{ $t_pair[0]->get_all_DBEntries } ) {
      if ( $entry->dbname eq 'Vega_transcript' ) {
        if ( $entry->primary_id eq $entry->display_id ) {

          my $enstentry =
            new Bio::EnsEMBL::DBEntry(
                                     -primary_id    => $entry->primary_id,
                                     -display_id    => $entry->display_id,
                                     -version       => 1,
                                     -release       => 1,
                                     -priority      => 1,
                                     -xref_priority => 0,
                                     -dbname => 'shares_CDS_and_UTR_with_OTTT'
            );

          $enstentry->status("XREF");

          $t_pair[0]->add_DBEntry($enstentry);
        }
      }
    }
    # Transfer the supporting features both for transcript and exon of
    # the transcript to delete to the transcript we keep
    $self->transfer_supporting_features( $delete_t, $t_pair[0] );

    # We add attributes to the havana transcript showing which
    # supporting features where used for the transcript in Havana
    #
    #$self->add_havana_attribute($t_pair[0],$t_pair[0]);

  } ## end elsif ( $delete_t == $t_pair...
} ## end sub set_transcript_relation

sub add_havana_attribute {
  my ( $self, $transcript, $trans_to_add_attrib ) = @_;

  my %evidence;
  my %t_evidence;

  foreach my $tsf ( @{ $transcript->get_all_supporting_features } ) {
    $t_evidence{ $tsf->hseqname } = 1;
  }

  foreach my $te_key ( keys %t_evidence ) {
    #print "Adding special attrib\n";

    if ( $te_key->isa("Bio::EnsEMBL::DnaPepAlignFeature") ) {
      my $attribute =
        Bio::EnsEMBL::Attribute->new(
        -CODE => 'tp_ott_support',
        -NAME => 'otter protein transcript support',
        -DESCRIPTION => 'Evidence ID that was used as protein transcript supporting feature for building a gene in Vega',
        -VALUE => $te_key );

      $trans_to_add_attrib->add_Attributes($attribute);

    }

    if ( $te_key->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
      my $attribute =
        Bio::EnsEMBL::Attribute->new(
        -CODE => 'td_ott_support',
        -NAME => 'otter dna transcript support',
        -DESCRIPTION => 'Evidence ID that was used as cdna transcript supporting feature for building a gene in Vega',
        -VALUE => $te_key );

      $trans_to_add_attrib->add_Attributes($attribute);

    }
  } ## end foreach my $te_key ( keys %t_evidence)

  foreach my $exon ( @{ $transcript->get_all_Exons } ) {
    foreach my $sf ( @{ $exon->get_all_supporting_features } ) {
      $evidence{ $sf->hseqname } = 1;
    }
  }

  foreach my $ev_key ( keys %evidence ) {
    #print "Adding special attrib\n";
    if ( $ev_key->isa("Bio::EnsEMBL::DnaPepAlignFeature") ) {
      my $attribute =
        Bio::EnsEMBL::Attribute->new(
        -CODE => 'ep_ott_support',
        -NAME => 'otter protein exon support',
        -DESCRIPTION => 'Evidence ID that was used as protein exon supporting feature for building a gene in Vega',
        -VALUE => $ev_key );

      $trans_to_add_attrib->add_Attributes($attribute);
    }

    if ( $ev_key->isa("Bio::EnsEMBL::DnaPepAlignFeature") ) {
      my $attribute =
        Bio::EnsEMBL::Attribute->new(
        -CODE => 'ed_ott_support',
        -NAME => 'otter dna exon support',
        -DESCRIPTION => 'Evidence ID that was used as cdna exon supporting feature for building a gene in Vega',
        -VALUE => $ev_key );

      $trans_to_add_attrib->add_Attributes($attribute);
    }
  } ## end foreach my $ev_key ( keys %evidence)
} ## end sub add_havana_attribute

sub transfer_supporting_features {
  my ( $self, $delete_t, $transcript ) = @_;

  #print "TRANSCRIPT IS :  ", $transcript,"\n";

  my @exon_features;

  # Delete all the supporting features for the Havana Transcript
  #$transcript->flush_supporting_features;

  my @delete_tsf = @{ $delete_t->get_all_supporting_features };
  #my @transcript_sf = @{ $transcript->get_all_supporting_features };

  # print "NUMBER OF TRANSCRIPT SF: ",scalar(@transcript_sf),"\n";
  # print " AND DELETE TSF: ", scalar(@delete_tsf),"\n";
DTSF: foreach my $dtsf (@delete_tsf) {
    next DTSF unless $dtsf->isa("Bio::EnsEMBL::FeaturePair");

    $transcript->add_supporting_features($dtsf);
  }

  my @delete_e = @{ $delete_t->get_all_Exons };
  my @exons    = @{ $transcript->get_all_Exons };

  if ( scalar(@delete_e) == scalar(@exons) ) {

    my $e;
    for ( $e = 0, $e < scalar(@delete_e), $e++ ) {
      if ( $delete_e[$e] && $exons[$e] ) {
        $self->transfer_supporting_evidence( $delete_e[$e], $exons[$e] );
      }
    }
  }

  # print "NUMBER AFT ADDITTION: ", scalar(@{ $transcript->get_all_supporting_features }),"\n";
} ## end sub transfer_supporting_features

sub _remove_transcript_from_gene {
  my ( $self, $gene, $trans_to_del ) = @_;

  #  foreach my $entry(@{ $gene->get_all_DBEntries}){
  #   print "ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
  #  }

  my @newtrans;
  foreach my $trans ( @{ $gene->get_all_Transcripts } ) {
    if ( $trans != $trans_to_del ) {
      push @newtrans, $trans;
    }
  }

  # The naughty bit!
  $gene->{_transcript_array} = [];

  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }

  return scalar(@newtrans);
} ## end sub _remove_transcript_from_gene


############################################################

sub _make_shared_exons_unique {
  my ( $self, $genes ) = @_;
  my @pruned_genes;
  foreach my $gene (@$genes) {

    # Make different exon objects that are shared between transcripts
    # (regarding attributes: start, end, etc) into unique exon
    # objects.
    my $new_gene = $self->prune_Exons($gene);
    push( @pruned_genes, $new_gene );
  }
  return \@pruned_genes;
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
  my $ensemblslice =
    $self->fetch_sequence( $self->input_id, $self->ensembl_db );
  my $havanaslice =
    $self->fetch_sequence( $self->input_id, $self->havana_db );

  # Fetch Ensembl genes
  print STDERR "Fetching ensembl genes\n";

  foreach my $ebiotype ( @{$ENSEMBL_INPUT_CODING_TYPE} ) {
  EGENE:
    foreach my $egene ( @{ $ensemblslice->get_all_Genes_by_type($ebiotype,undef,1) } ) {
      $egene->load();
      # Don't add those genes that contain only transcripts imported
      # from HAVANA (this is important during a merge update)

      if ( $egene->analysis->logic_name() eq $HAVANA_LOGIC_NAME ) {
        next EGENE;
      } else {
        #print "in get_Genes: coding genes " . $egene->stable_id . "\n";
        push( @genes, $egene );
      }
    }
  }

  # Fetch Ensembl Processed transcripts
  foreach my $eprocessedbt ( @{$ENSEMBL_INPUT_PROCESSED_TYPE} ) {
  PROCESSED:
    foreach my $eprocessedgene ( @{ $ensemblslice->get_all_Genes_by_type($eprocessedbt,undef,1) } ) {
      $eprocessedgene->load();
      # Don't add those genes that contain only transcripts imported
      # from HAVANA (this is important during a merge update)
      if ( $eprocessedgene->analysis->logic_name() eq $HAVANA_LOGIC_NAME ) {
        next PROCESSED;
      } else {
        #print "in get_Genes: Processed transcripts  " . $eprocessedgene->stable_id . "\n";
        push( @processedgenes, $eprocessedgene );
      }
    }
  }

  # Fetch Ensembl pseudogenes
  foreach my $epseudobt ( @{$ENSEMBL_INPUT_PSEUDO_TYPE} ) {
  EPSEUDOGENE:
    foreach my $epseudogene ( @{ $ensemblslice->get_all_Genes_by_type($epseudobt,undef,1) } ) {
      $epseudogene->load();
      # Don't add those genes that contain only transcripts imported
      # from HAVANA (this is important during a merge update)
      if ( $epseudogene->analysis->logic_name() eq $HAVANA_LOGIC_NAME ) {
        next EPSEUDOGENE;
      } else {
        #print "in get_Genes: pseudogenes " . $epseudogene->stable_id . "\n";
        push( @pseudogenes, $epseudogene );
      }
    }
  }

  print STDERR "Retrieved "
    . scalar(@genes)
    . " genes of types: "
    . join( ", ", @{$ENSEMBL_INPUT_CODING_TYPE} ) . "\n";

  print STDERR "Retrieved "
    . scalar(@processedgenes)
    . " 'processed transcript' genes of types: "
    . join( ", ", @{$ENSEMBL_INPUT_PROCESSED_TYPE} ) . "\n";

  print STDERR "Retrieved "
    . scalar(@pseudogenes)
    . " pseudogenes of types: "
    . join( ", ", @{$ENSEMBL_INPUT_PSEUDO_TYPE} ) . "\n";

  # Fetch Havana genes
  print STDERR "Fetching havana genes\n";
  foreach my $hbiotype ( @{$HAVANA_INPUT_CODING_TYPE} ) {
    foreach my $hgene ( @{ $havanaslice->get_all_Genes_by_type($hbiotype,undef,1) } ) {
      $hgene->load();
      # We change the biotype of the havana genes/transcripts as it
      # could happend to be the same as the ensembl ones
      my $biotype = $hgene->biotype . "_hav";
      $hgene->biotype($biotype);
      foreach my $htran ( @{ $hgene->get_all_Transcripts } ) {
        my $tbiotype = $htran->biotype . "_hav";
        $htran->biotype($tbiotype);
      }
      push( @hgenes, $hgene );
    }
  }

  print STDERR "Fetching havana 'processed transcript' genes\n";
  foreach my $hprocessedbiotype ( @{$HAVANA_INPUT_PROCESSED_TYPE} ) {
    foreach my $hprocessedgene ( @{ $havanaslice->get_all_Genes_by_type($hprocessedbiotype,undef,1) } ) {
      $hprocessedgene->load();
      # We change the biotype of the havana genes/transcripts as it
      # could happend to be the same as the ensembl ones
      my $processedbiotype = $hprocessedgene->biotype . "_hav";
      $hprocessedgene->biotype($processedbiotype);

      foreach my $hprocessedtran ( @{ $hprocessedgene->get_all_Transcripts } )
      {
        my $tprocessedbiotype = $hprocessedtran->biotype . "_hav";
        $hprocessedtran->biotype($tprocessedbiotype);

      }
      push( @hprocessedgenes, $hprocessedgene );
    }
  }

  #Fetch Havana pseudogenes
  print STDERR "Fetching havana pseudogenes\n";
  foreach my $hpseudobt ( @{$HAVANA_INPUT_PSEUDO_TYPE} ) {
    foreach my $hpseudogene ( @{ $havanaslice->get_all_Genes_by_type($hpseudobt,undef,1) } ) {
      $hpseudogene->load();
      # We change the biotype of the havana genes/transcripts as it
      # could happend to be the same as the ensembl ones
      my $biotype = $hpseudogene->biotype . "_hav";
      $hpseudogene->biotype($biotype);
      foreach my $htran ( @{ $hpseudogene->get_all_Transcripts } ) {
        my $tbiotype = $htran->biotype . "_hav";
        $htran->biotype($tbiotype);
      }
      push( @hpseudogenes, $hpseudogene );
    }
  }

  print STDERR "Retrieved "
    . scalar(@hgenes)
    . " genes of types: "
    . join( ", ", @{$HAVANA_INPUT_CODING_TYPE} ) . "\n";

  print STDERR "Retrieved "
    . scalar(@hprocessedgenes)
    . " 'processed transcript' genes of types: "
    . join( ", ", @{$HAVANA_INPUT_PROCESSED_TYPE} ) . "\n";

  print STDERR "Retrieved "
    . scalar(@hpseudogenes)
    . " pseudogenes of types: "
    . join( ", ", @{$HAVANA_INPUT_PSEUDO_TYPE} ) . "\n";

  # We want to keep the HAVANA genes as they are as they asked us to
  # keep their gene clustering untouched
  $self->havana_Coding_Genes( \@hgenes );
  $self->havana_Processed_Transcript_Genes( \@hprocessedgenes );
  $self->havana_Pseudogene_Genes( \@hpseudogenes );

  # We processed the ensembl genes to remove any trace of previous
  # transcript merges in case you are running a merge update instead
  # of a fresh clean merge.
  @transcripts = @{ $self->check_merge_transcript_status( \@genes ) };

  @processedtranscripts =
    @{ $self->check_merge_transcript_status( \@processedgenes ) };

  @pseudotranscripts =
    @{ $self->check_merge_transcript_status( \@pseudogenes ) };

  # Join all the gene set together
  print STDERR "Finished fetching genes\n";
  $self->combined_Transcripts( \@transcripts );
  $self->combined_Processed_Transcripts( \@processedtranscripts );
  $self->combined_PseudoTranscripts( \@pseudotranscripts );
} ## end sub get_Genes

sub check_merge_transcript_status {
  my ( $self, $genes ) = @_;

  #print "Checking premerge gene status\n";

  my @transcripts;
  foreach my $gene (@$genes) {
  TRANSCRIPT:
    foreach my $tran ( @{ $gene->get_all_Transcripts } ) {

     #First we remove HAVANA only transcripts that are present in merged genes
      if (    $gene->analysis->logic_name() eq $MERGED_GENE_LOGIC_NAME
           && $tran->analysis->logic_name() eq $HAVANA_LOGIC_NAME )
      {
        next TRANSCRIPT;
      } elsif (
              $tran->analysis->logic_name() eq $MERGED_TRANSCRIPT_LOGIC_NAME )
      {
        # In case of a merged transcript we want to distinguish the
        # ones that came from HAVANA that have same CDS but different
        # UTR structure as we want to remove them. This is important
        # for a merge update to avoid them been wrongly identified as
        # share CDS and UTR
        my $share_enst        = 0;
        my $share_cds_and_utr = 0;
        my @dbentries         = @{ $tran->get_all_DBEntries };
        foreach my $dbentry (@dbentries) {
          if ( $dbentry->dbname eq "shares_CDS_with_ENST" ) {

            #print "On transcript: ",$tran->dbID," This is a HAVANA shares ENST\n";
            $share_enst = 1;
          }
          if ( $dbentry->dbname eq "shares_CDS_and_UTR_with_OTTT" ) {
            $share_cds_and_utr = 1;

            #print "On transcript: ",$tran->dbID," This is a HAVANA shares CDS and UTR\n";
          }
        }
        if ( $share_enst == 1 && $share_cds_and_utr == 0 ) {
          next TRANSCRIPT;
        }
      }

      $self->flush_xref($tran);

      # Check if a transcript is in the discarded genes database
      # before adding it to the merging list.
      if ( $self->check_transcript_in_external_db('discarded', $tran) != 0 ) {
        #print "Transcript added\n";
        push( @transcripts, $tran );

      }
    } ## end foreach my $tran ( @{ $gene...
  } ## end foreach my $gene (@$genes)
  return \@transcripts;
} ## end sub check_merge_transcript_status

sub check_transcript_in_external_db {
  my ( $self, $dbname, $trans ) = @_;

  my @exons   = @{ $trans->get_all_Exons };
  my @t_exons = @{ $trans->get_all_translateable_Exons };

  #print "--->>> exons: " . scalar(@exons) . "\n";
  #print "--->>> t_exons: " . scalar(@t_exons) . "\n";

  my $ext_slice;
  if ( $dbname =~ /discarded/ ) {
    $ext_slice =
      $self->discarded_db->get_SliceAdaptor->fetch_by_region( 'toplevel',
                           $trans->slice->seq_region_name,
                           $trans->seq_region_start, $trans->seq_region_end );
  } else {    # external_db = ccds_db
    $ext_slice =
      $self->ccds_db->get_SliceAdaptor->fetch_by_region( 'toplevel',
                     $trans->slice->seq_region_name,
                     $trans->seq_region_start, $trans->seq_region_end );

  }

EXT_GENE:
  foreach my $ext_gene ( @{ $ext_slice->get_all_Genes } ) {
  EXT_TRANS:
    foreach my $ext_trans ( @{ $ext_gene->get_all_Transcripts } ) {
      my @ext_exons = @{ $ext_trans->get_all_Exons };
      if ( $dbname =~ /discarded/ ) {
        if ( scalar(@exons) == scalar(@ext_exons) ) {
          #print "Number of exons: ",scalar(@exons),"\n";
          for ( my $i = 0 ; $i < scalar(@exons) ; $i++ ) {
            if (
              $exons[$i]->seq_region_start != $ext_exons[$i]->seq_region_start
              || $exons[$i]->strand != $ext_exons[$i]->strand
              || $exons[$i]->seq_region_end != $ext_exons[$i]->seq_region_end )
            {
              # If you enter here means that these two transcripts are
              # not the same.

              #print "transcript exon coordinates are different\n";
              next EXT_TRANS;
            }
          }
          #print "transcript found in discarded db\n";
          return 0;
        } else {
          #print "discarded db: number of exons is different\n";
          next EXT_GENE;
        }
      } else {
        #print "comparing ccds: " . $ext_trans->stable_id() . " vs trans: " .  $trans->dbID . " (" . $trans->stable_id() . ")\n";
        if (@t_exons) {
          #print "--->>> ccds exons: " . scalar(@ext_exons) . "\n";
          #print "--->>> trans t_exons: " . scalar(@t_exons) . "\n";

          if ( scalar(@t_exons) == scalar(@ext_exons) ) {
            #print $trans->dbID . " :: " . $trans->stable_id() . "\n";
            for ( my $i = 0 ; $i < scalar(@t_exons) ; $i++ ) {
              #print "te start: " . $t_exons[$i]->coding_region_start($trans) . " te_end: " . $t_exons[$i]->coding_region_end($trans) . "\n";
              #print "ext exon start: " .  $ext_exons[$i]->seq_region_start . " ext exon end: " .  $ext_exons[$i]->seq_region_end . "\n";

              if ( $t_exons[$i]->coding_region_start($trans) != $ext_exons[$i]->seq_region_start
                   || $t_exons[$i]->strand != $ext_exons[$i]->strand
                   || $t_exons[$i]->coding_region_end($trans) != $ext_exons[$i]->seq_region_end )
              {
                next EXT_TRANS;
              }
            }
            #print "\t--->>> transcript found in ccds db\n";
            return 0;
          } else {
            #print "ccds db: number of (translatable) exons is different\n";
            next EXT_GENE;
          }
        } else {
          if ( scalar(@exons) == scalar(@ext_exons) ) {
            for ( my $i = 0 ; $i < scalar(@exons) ; $i++ ) {
              #print "exon start: " . $exons[$i]->seq_region_start . " vs ext_exon start: " . $ext_exons[$i]->seq_region_start . "\n";
              if ( $exons[$i]->seq_region_start != $ext_exons[$i]->seq_region_start
                   || $exons[$i]->strand != $ext_exons[$i]->strand
                   || $exons[$i]->seq_region_end != $ext_exons[$i]->seq_region_end )
              {
                next EXT_TRANS;
              }
            }
            # If you are here means that both transcripts are the same.

            #print "transcript found in ccds db\n";
            return 0;
          } else {

            # If you enter here means that these two transcripts are not
            # the same.

            #print "ccds db: number of (non-coding) exons is different\n";
            next EXT_GENE;
          }
        } ## end else [ if (@t_exons)
      } ## end else [ if ( $dbname =~ /discarded/)
    } ## end foreach my $ext_trans ( @{ ...
  } ## end foreach my $ext_gene ( @{ $ext_slice...

  # If we reach here it means that no transcript in the external db
  # is the same as our transcript so we keep it.

  #print "transcript not found in external db\n";
  return 1;
} ## end sub check_transcript_in_external_db

sub flush_xref {
  my ( $self, $transcript ) = @_;

  my @newxrefs;
  #print "THIS IS WHAT NEEDS EMPTYING: ",$transcript->get_all_DBEntries,"\n";
  foreach my $tran_xref ( @{ $transcript->get_all_DBEntries } ) {
    if (    $tran_xref->dbname ne "shares_CDS_and_UTR_with_OTTT"
         && $tran_xref->dbname ne "shares_CDS_with_OTTT"
         && $tran_xref->dbname ne "shares_CDS_with_ENST"
         && $tran_xref->dbname ne "OTTT"
         && $tran_xref->dbname ne "OTTP"
         && $tran_xref->dbname ne "OTTG"
         && $tran_xref->dbname ne "Vega_Transcript"
         && $tran_xref->dbname ne "Vega_Gene" )
    {
      push( @newxrefs, $tran_xref );
    }
  }

  # The naughty bit!
  $transcript->{dbentries} = [];
  #$transcript->{display_xref} = [];

  foreach my $newxref (@newxrefs) {
    $transcript->add_DBEntry($newxref);
  }

} ## end sub flush_xref

###########################################################c

=head2 cluster_Transcripts

 Description : It separates transcripts according to strand and then clusters 
               each set of transcripts by calling _cluster_Transcripts_by_genomic_range()
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster

=cut

sub cluster_Transcripts {
  my ( $self, $transcripts ) = @_;

  my @forward_transcripts;
  my @reverse_transcripts;

  foreach my $transcript (@$transcripts) {
    my @exons = @{ $transcript->get_all_Exons };
    if ( $exons[0]->strand == 1 ) {
      push( @forward_transcripts, $transcript );
    } else {
      push( @reverse_transcripts, $transcript );
    }
  }

  my @forward_clusters;
  my @reverse_clusters;

  if (@forward_transcripts) {
    @forward_clusters =
      @{ $self->_cluster_Transcripts_by_genomic_range( \@forward_transcripts )
      };
  }
  if (@reverse_transcripts) {
    @reverse_clusters =
      @{ $self->_cluster_Transcripts_by_genomic_range( \@reverse_transcripts )
      };
  }
  my @clusters;
  if (@forward_clusters) {
    push( @clusters, @forward_clusters );
  }
  if (@reverse_clusters) {
    push( @clusters, @reverse_clusters );
  }
  return \@clusters;
} ## end sub cluster_Transcripts

############################################################

=head2 _cluster_Transcripts_by_genomic_range

 Description : It clusters transcripts according to genomic overlap
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster

=cut

sub _cluster_Transcripts_by_genomic_range {
  my ( $self, $mytranscripts ) = @_;

  # First sort the transcripts
  my @transcripts = sort {
    $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end
  } @$mytranscripts;

  # Create a new cluster
  my $cluster =
    Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new();
  my $count = 0;
  my @cluster_starts;
  my @cluster_ends;
  my @clusters;

  # Put the first transcript into these cluster
  $cluster->put_Transcripts( [ $transcripts[0] ] );

  $cluster_starts[$count] = $transcripts[0]->start;
  $cluster_ends[$count]   = $transcripts[0]->end;

  # Store the list of clusters
  push( @clusters, $cluster );

  # Loop over the rest of the transcripts
LOOP1:
  for ( my $c = 1 ; $c <= $#transcripts ; $c++ ) {
    #print STDERR "\nIn cluster ".($count+1)."\n";
    #print STDERR "start: $cluster_starts[$count] end: $cluster_ends[$count]\n";
    #print STDERR "comparing:\n";

    if (
         !(    $transcripts[$c]->end < $cluster_starts[$count]
            || $transcripts[$c]->start > $cluster_ends[$count] ) )
    {
      $cluster->put_Transcripts( [ $transcripts[$c] ] );

      # Re-adjust size of cluster
      if ( $transcripts[$c]->start < $cluster_starts[$count] ) {
        $cluster_starts[$count] = $transcripts[$c]->start;
      }
      if ( $transcripts[$c]->end > $cluster_ends[$count] ) {
        $cluster_ends[$count] = $transcripts[$c]->end;
      }
    } else {
      # Else, create a new cluster with this feature
      $count++;
      $cluster =
        Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new();
      $cluster->put_Transcripts( [ $transcripts[$c] ] );
      $cluster_starts[$count] = $transcripts[$c]->start;
      $cluster_ends[$count]   = $transcripts[$c]->end;

      # Store it in the list of clusters
      push( @clusters, $cluster );
    }
  } ## end for ( my $c = 1 ; $c <=...
  return \@clusters;
} ## end sub _cluster_Transcripts_by_genomic_range

sub combine_gene_clusters {
  my ( $self, $preliminary_coding_genes, $preliminary_pseudo_genes ) = @_;

  my @coding_genes = sort {
    $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end
  } @{ $preliminary_coding_genes };

  my @pseudo_genes = sort {
    $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end
  } @{ $preliminary_pseudo_genes };

  my @unclestered_pseudos;
  my $is_pseudo_havana = 0;

CLUSTER:
  foreach my $pseudo_gene (@pseudo_genes) {
    print "Pseudogene biotype: ", $pseudo_gene->biotype,"\tstart: " .  $pseudo_gene->seq_region_start . "\n";
    if ( $pseudo_gene->biotype =~ /hav/ ) {
      $is_pseudo_havana = 1;
    }
    my $pseudo_status = 0;
  OVERLAP:
    foreach my $coding_gene (@coding_genes) {

      #print "--->>> \$is_pseudo_havana = $is_pseudo_havana\n";
      #print "--->>> coding gene biotype: " . $coding_gene->biotype . "\n";

      # I need to add this line as I will be making some coding genes
      # UNDEF as I add them to pseudogenes.
      next OVERLAP unless $coding_gene;

      foreach my $trans (@{ $coding_gene->get_all_Transcripts }) {

        if ($trans->biotype =~ /hav/) {
          next OVERLAP;
        }

        #print "checking if trans " . $trans->dbID . " is CCDS\n";
        if ( $self->check_transcript_in_external_db('ccds', $trans) == 0 ) {
          print "Keeping ensembl transcript as it is CCDS " . $trans->dbID . "\n";
          next OVERLAP;
        }
      }

      if ( $is_pseudo_havana == 1 && $coding_gene->biotype =~ /hav/ ) {
        print "Jumping over havana coding gene! " . $coding_gene->dbID .
              " biotype: " . $coding_gene->biotype . "\n";
        next OVERLAP;
      }

      if (    $coding_gene->end >= $pseudo_gene->start
           && $coding_gene->start <= $pseudo_gene->end )
      {
        #print "--->>> coding gene smaller: " . $coding_gene->biotype . "\n";
        #print "--->>> \$is_pseudo_havana = $is_pseudo_havana\n";

        my $coding_length = $self->get_coding_length($coding_gene);

        my $cg_exons = $self->get_coding_exons_for_gene($coding_gene);
        my $pg_exons = $pseudo_gene->get_all_Exons();

        foreach my $cg_exon ( @{$cg_exons} ) {
          foreach my $pg_exon ( @{$pg_exons} ) {

            if (    $cg_exon->overlaps($pg_exon)
                 && $cg_exon->strand == $pg_exon->strand )
            {
              # Check if the overlap covers at least 10 percent of the
              # coding region of the longest transcript in the gene.

              # NOTE!!! This check is a bit experimental.
              if ( $self->overlap_percent( $cg_exon, $pg_exon, $coding_length ) > 10 ) {
                if ( $is_pseudo_havana == 1 ) {
                  print "I'm looking into merging a pseudogene\n";

                  # As the pseudogene is Havana I will add the coding
                  # transcripts to the pseudo and remove the translation.
                  foreach my $c_transcript ( @{ $coding_gene->get_all_Transcripts } ) {
                    print "\ncombine_gene_clusters --->>> changing coding trans to non-coding: "
                      . $c_transcript->dbID . " (" . $c_transcript->biotype . ")\n";

                    #print "NOTE!!! removing the translation, adding trans to pseudogene, "
                    #    . "and changing the biotype to: " . $pseudo_gene->biotype . "_ens_merged\n";

                    $c_transcript->{translation} = undef;
                    $c_transcript->biotype($pseudo_gene->biotype .  "_ens_merged");
                    $pseudo_gene->add_Transcript($c_transcript);
                  }
                  $coding_gene = undef;
                  next OVERLAP;
                } else {

                  # Have to add all the transcripts of the pseudo to
                  # the gene and remove the pseudogene.
                  foreach my $p_transcript (
                                      @{ $pseudo_gene->get_all_Transcripts } )
                  {
                    #print "\ncombine_gene_clusters --->>> adding pseudo transcript to the coding gene: " . $p_transcript->dbID . "\n";
                    $coding_gene->add_Transcript($p_transcript);
                  }
                  $pseudo_status = 1;
                  #SMJS Need to look closer at this (used to be next OVERLAP;)
                  print " !!!!!!!!! Done noncoding-coding merge - jumping to next noncoding gene\n";
                  next CLUSTER;
                }
              } ## end if ( $self->overlap_percent...
            } ## end if ( $cg_exon->overlaps...
          } ## end foreach my $pg_exon ( @{$pg_exons...
        } ## end foreach my $cg_exon ( @{$cg_exons...
      } ## end if ( $coding_gene->end...
    } ## end foreach my $coding_gene (@coding_genes)
    unless ( $pseudo_status == 1 ) {
      push( @unclestered_pseudos, $pseudo_gene );
    }
  } ## end foreach my $pseudo_gene (@pseudo_genes)

  my @final_clustered_genes;

  foreach my $coding (@coding_genes) {
    if ($coding) {
      #print "YOUR CODING:", $coding,"\n";
      push( @final_clustered_genes, $coding );
    }
  }

  push( @final_clustered_genes, @unclestered_pseudos );

  return \@final_clustered_genes;
} ## end sub combine_gene_clusters

sub get_coding_length {
  my ( $self, $gene ) = @_;

  my $length = 0;

  foreach my $transcript ( @{ $gene->get_all_Transcripts } ) {
    if ( $transcript->translate ) {
      if ( $transcript->translate->length > $length ) {
        $length = $transcript->translate->length;
      }
    }
  }
  return $length;
}

sub clust_overlap_percent {
  my ( $self, $clust, $transcript ) = @_;

  my $trans_length;
  my @t_exons;
  if ( $transcript->translate ) {
    $trans_length = $transcript->translate->length;
    @t_exons      = @{ $transcript->get_all_translateable_Exons };
  } else {
    @t_exons      = @{ $transcript->get_all_Exons };
    $trans_length = $transcript->length;
  }

  my $max_overlap     = 0;
  my $min_g_t_overlap = 100;
  foreach my $g_trans ( @{$clust} ) {
    my $percent     = 0;
    my $g_t_length  = $g_trans->length;
    my $g_t_percent = 0;
    foreach my $g_exon ( @{ $g_trans->get_all_Exons } ) {

      foreach my $t_exon (@t_exons) {
        if (    $g_exon->seq_region_start < $t_exon->seq_region_end
             && $g_exon->seq_region_end > $t_exon->seq_region_start )
        {

          my $low  = 0;
          my $high = 0;

          if ( $g_exon->seq_region_start >= $t_exon->seq_region_start ) {
            $low = $g_exon->seq_region_start;
          } else {
            $low = $t_exon->seq_region_start;
          }
          if ( $g_exon->seq_region_end <= $t_exon->seq_region_end ) {
            $high = $g_exon->seq_region_end;
          } else {
            $high = $t_exon->seq_region_end;
          }

          my $overlap_length = $high - $low;
          $g_t_percent += ( $overlap_length/$g_t_length )*100;
          $percent     += ( $overlap_length/$trans_length )*100;
        }
      }
    } ## end foreach my $g_exon ( @{ $g_trans...

    # First we check what is the minimun coverage of the havana
    # transcripts in the gene cluster If this value is close to 100%
    # we consider that the ensembl transcript and the gene should
    # be merged no matter what. This is usefull if there are a few
    # partial genes from havana that overlap a transcript in Ensembl
    if ( $g_t_percent < $min_g_t_overlap ) {
      $min_g_t_overlap = $g_t_percent;
    }

    # We also want to know which gene from havana has the best
    # coverage over the ensembl transcript. We will use this value
    # in case we have to Havana genes that we don't want to merge
    # and we want to decide which one should merge with the ensembl
    # transcript
    if ( $percent > $max_overlap ) {
      $max_overlap = $percent;
    }
  } ## end foreach my $g_trans ( @{$clust...
  return $max_overlap;
} ## end sub clust_overlap_percent


sub overlap_percent {

  my ( $self, $cg_exon, $pg_exon, $coding_length ) = @_;
  my $low  = 0;
  my $high = 0;

  if ( $cg_exon->seq_region_start >= $pg_exon->seq_region_start ) {
    $low = $cg_exon->seq_region_start;
  } else {
    $low = $pg_exon->seq_region_start;
  }
  if ( $cg_exon->seq_region_end <= $pg_exon->seq_region_end ) {
    $high = $cg_exon->seq_region_end;
  } else {
    $high = $pg_exon->seq_region_end;
  }

  my $overlap_length = $high - $low;

  my $percent = ( $overlap_length/$coding_length )*100;

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
    my ( $self, $trans ) = @_;

    if ( exists( $coding_exon_cache{$trans} ) ) {
      return $coding_exon_cache{$trans};
    } else {
      my %coding_hash;

      next if ( !$trans->translation );
      foreach my $exon ( @{ $trans->get_all_translateable_Exons } ) {
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
  my ( $self, $gene ) = @_;

  my @coding;

  foreach my $trans ( @{ $gene->get_all_Transcripts } ) {
    next if ( !$trans->translation );
    foreach my $exon ( @{ $trans->get_all_translateable_Exons } ) {
      push @coding, $exon;
    }
  }

  return \@coding;

}



############################################################

sub check_Clusters {
  my ( $self, $num_transcripts, $clusters ) = @_;

  # Safety checks
  my $ntrans = 0;

  my $cluster_num = 0;

  my %trans_check_hash;
  foreach my $cluster (@$clusters) {
    $ntrans += scalar(@$cluster);

    foreach my $trans (@$cluster) {

      if ( defined( $trans_check_hash{$trans} ) ) {
        $self->throw(
              "Transcript " . $trans->dbID . " already exist in clusters.\n" );
      }
      $trans_check_hash{$trans} = 1;
    }
    if ( !scalar(@$cluster) ) {
      $self->throw("Empty cluster");
    }
  }
  if ( $ntrans < $num_transcripts ) {
    $self->throw( "Not all transcripts have been added into clusters: "
                . "$ntrans and " . $num_transcripts . " \n" );
  }
  #end safety checks
  return;
} ## end sub check_Clusters


############################################################

sub prune_Exons {
  my ( $self, $gene ) = @_;

  my @unique_Exons;

  # Keeping track of all unique exons found so far to avoid making
  # duplicates. Need to be very careful about translation->start_Exon
  # and translation->end_Exon

  foreach my $tran ( @{ $gene->get_all_Transcripts } ) {
    my @newexons;
    foreach my $exon ( @{ $tran->get_all_Exons } ) {
      my $found;
      #always empty
    UNI: foreach my $uni (@unique_Exons) {
        if (    $uni->start == $exon->start
             && $uni->end == $exon->end
             && $uni->strand == $exon->strand
             && $uni->phase == $exon->phase
             && $uni->end_phase == $exon->end_phase )
        {
          $found = $uni;
          last UNI;
        }
      }
      if ( defined($found) ) {
        push( @newexons, $found );
        if ( $tran->translation ) {
          if ( $exon == $tran->translation->start_Exon ) {
            $tran->translation->start_Exon($found);
          }
          if ( $exon == $tran->translation->end_Exon ) {
            $tran->translation->end_Exon($found);
          }
        }
      } else {
        push( @newexons,     $exon );
        push( @unique_Exons, $exon );
      }
    } ## end foreach my $exon ( @{ $tran...
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }

    #print "Uniq_tran sid: ",$tran->dbID,"\n";

  } ## end foreach my $tran ( @{ $gene...
  return $gene;
} ## end sub prune_Exons

############################################################

=head2 prune_features

 Description: prunes out duplicated features
 Returntype : array of Bio::EnsEMBL::SeqFeature
 Args       : array of Bio::EnsEMBL::SeqFeature

=cut

sub prune_features {
  my ( $self, $feature_hash ) = @_;
  my @pruned;

ID:
  foreach my $id ( keys %{$feature_hash} ) {
    my @features = @{ $feature_hash->{$id} };
    @features = sort { $a->start <=> $b->start } @features;

    unless (@features) {
      print STDERR "No features here for id: $id\n";
      next ID;
    }
    while ( @features && !defined $features[0] ) {
      #print STDERR "jumping an undefined feature\n";
      shift @features;
    }

    my $prev = -1;

  FEATURE:
    foreach my $f (@features) {
      if (    $prev != -1
           && $f->hseqname eq $prev->hseqname
           && $f->start == $prev->start
           && $f->end == $prev->end
           && $f->hstart == $prev->hstart
           && $f->hend == $prev->hend
           && $f->strand == $prev->strand
           && $f->hstrand == $prev->hstrand )
      {
        #keep the one with highest score
        if ( $f->score > $prev->score ) {
          $prev->score( $f->score );
        }
        #print STDERR "pruning duplicated feature\n";
        #print STDERR "previous: ".$prev->gffstring."\n";
        #print STDERR "thisone : ".$f->gffstring."\n";
        next FEATURE;
      } else {
        push( @pruned, $f );
        $prev = $f;
      }
    }
  } ## end foreach my $id ( keys %{$feature_hash...
  return \@pruned;
} ## end sub prune_features
############################################################

=head2 transfer_supporting_evidence

 Title   : transfer_supporting_evidence
 Usage   : $self->transfer_supporting_evidence($source_exon, $target_exon)
 Function: Transfers supporting evidence from source_exon to target_exon, 
           after checking the coordinates are sane and that the evidence is not already in place.
 Returns : nothing, but $target_exon has additional supporting evidence

=cut

sub transfer_supporting_evidence {
  my ( $self, $source_exon, $target_exon ) = @_;

  my @target_sf = @{ $target_exon->get_all_supporting_features };

  # keep track of features already transferred, so that we do not duplicate
  my %unique_evidence;
  my %hold_evidence;

SOURCE_FEAT:
  foreach my $feat ( @{ $source_exon->get_all_supporting_features } ) {
    next SOURCE_FEAT unless $feat->isa("Bio::EnsEMBL::FeaturePair");

    # skip duplicated evidence objects
    next SOURCE_FEAT if ( $unique_evidence{$feat} );

    # skip duplicated evidence
    if ( $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }
         { $feat->hstart }{ $feat->hend } )
    {
      #print STDERR "Skipping duplicated evidence\n";
      next SOURCE_FEAT;
    }

    #$self->print_FeaturePair($feat);

  TARGET_FEAT:
    foreach my $tsf (@target_sf) {
      next TARGET_FEAT unless $tsf->isa("Bio::EnsEMBL::FeaturePair");

      if (    $feat->start == $tsf->start
           && $feat->end == $tsf->end
           && $feat->strand == $tsf->strand
           && $feat->hseqname eq $tsf->hseqname
           && $feat->hstart == $tsf->hstart
           && $feat->hend == $tsf->hend )
      {

        #print STDERR "feature already in target exon\n";
        next SOURCE_FEAT;
      }
    }
    #print STDERR "from ".$source_exon->dbID." to ".$target_exon->dbID."\n";
    #$self->print_FeaturePair($feat);
    # I may need to add a paranoid check to see that no exons longer than the current one are transferred
    $target_exon->add_supporting_features($feat);
    $unique_evidence{$feat} = 1;
    $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }
      { $feat->hstart }{ $feat->hend } = 1;
  } ## end foreach my $feat ( @{ $source_exon...
} ## end sub transfer_supporting_evidence

=head2 update_biotypes

  Description: This check the biotypes of the merged transcript and genes and updates then to reflect the merge

=cut

sub update_biotypes {
  my ( $self, @genes ) = @_;

  my %pseudobiotypes;
  my %processedbiotypes;
  my %coding_biotypes;

  foreach my $epb ( @{$ENSEMBL_INPUT_PSEUDO_TYPE} ) {
    $pseudobiotypes{$epb} = 1;
  }
  foreach my $hpb ( @{$HAVANA_INPUT_PSEUDO_TYPE} ) {
    $pseudobiotypes{ $hpb . "_hav" } = 1;
  }
  foreach my $epb ( @{$ENSEMBL_INPUT_PROCESSED_TYPE} ) {
    $processedbiotypes{$epb} = 1;
  }
  foreach my $hpb ( @{$HAVANA_INPUT_PROCESSED_TYPE} ) {
    $processedbiotypes{ $hpb . "_hav" } = 1;
  }

  foreach my $ecb ( @{$ENSEMBL_INPUT_CODING_TYPE} ) {
    $coding_biotypes{$ecb} = 1;
  }
  foreach my $hcb ( @{$HAVANA_INPUT_CODING_TYPE} ) {
    $coding_biotypes{ $hcb . "_hav" } = 1;
  }

  foreach my $gene (@genes) {
    my %trans_types;

    my $has_pseudos   = 0;
    my $has_processed = 0;
    my $has_coding    = 0;
    #$gene->type($GB_GENE_OUTPUT_BIOTYPE);
    # poke the caches
    my %s_pfhash;

    foreach my $tran ( @{ $gene->get_all_Transcripts } ) {
      $trans_types{ $tran->biotype } = 1;

      foreach my $pseudobiotype ( keys %pseudobiotypes ) {
        if ( $tran->biotype =~ /$pseudobiotype/ ) {
          #print "HAS PSEUDOS NEW CODE WORKS\n";
          $has_pseudos = 1;
        }
      }
      foreach my $processedbiotype ( keys %processedbiotypes ) {
        if ( $tran->biotype =~ /$processedbiotype/ ) {
          #print "HAS PROCESSED NEW CODE WORKS\n";
          $has_processed = 1;
        }
      }
      foreach my $coding_biotype ( keys %coding_biotypes ) {
        if ( $tran->biotype =~ /$coding_biotype/ ) {
          #print "HAS CODING NEW CODE WORKS\n";
          $has_coding = 1;
        }
      }
    }

    # MANAGE OUTPUT BIOTYPES BEFORE WRITTING OUTPUT
    my $newbiotype;
    my $biotype_status;
    my $has_havana  = 0;
    my $has_ensembl = 0;
    my $has_merged  = 0;

    if ( $has_coding == 1 ) {
      $biotype_status = "protein_coding";
    } elsif ( $has_processed == 1 && $has_coding == 0 ) {
      $biotype_status = "processed_transcript";
    } elsif ( $has_pseudos == 1 && $has_coding == 0 && $has_processed == 0 ) {
      $biotype_status = "pseudogene";
    } else {
      print "ERROR: I should not really be here for gene biotype checks\n";
      $biotype_status = "weird_" . $gene->biotype;
    }

    foreach my $t_biotype ( keys %trans_types ) {
      if ( $t_biotype =~ /$MERGED_TRANSCRIPT_OUTPUT_TYPE/ ) {
        $has_merged = 1;
      } elsif ( $t_biotype =~ /_hav/ ) {
        $has_havana = 1;
      } else {
        $has_ensembl = 1;
      }
    }

    if ( ( $has_havana == 1 && $has_ensembl == 1 ) || $has_merged == 1 ) {
      $newbiotype = $biotype_status . $MERGED_GENE_OUTPUT_BIOTYPE;
      $gene->biotype($newbiotype);
    } elsif ( $has_havana == 1 && $has_ensembl == 0 && $has_merged == 0 ) {
      $newbiotype = $biotype_status . $HAVANA_GENE_OUTPUT_BIOTYPE;
      $gene->biotype($newbiotype);
    } elsif ( $has_ensembl == 1 && $has_havana == 0 && $has_merged == 0 ) {
      $newbiotype = $biotype_status . $ENSEMBL_GENE_OUTPUT_BIOTYPE;
      $gene->biotype($newbiotype);
    } else {
      $newbiotype = $biotype_status . "weird";
      $gene->biotype($newbiotype);
    }
  } ## end foreach my $gene (@genes)
} ## end sub update_biotypes

############################################################


############################################################
#
# GETSET METHODS
#
############################################################

# get/set method holding a reference to the db with genewise and combined genes,
# havana genes and discarded genes
# this reference is set in Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder

sub ensembl_db {
  my ( $self, $ensembl_db ) = @_;
  if ($ensembl_db) {
    $self->{_ensembl_db} = $ensembl_db;
  }

  return $self->{_ensembl_db};
}

sub havana_db {
  my ( $self, $havana_db ) = @_;
  if ($havana_db) {
    $self->{_havana_db} = $havana_db;
  }

  return $self->{_havana_db};
}

sub discarded_db {
  my ( $self, $discarded_db ) = @_;

  if ($discarded_db) {
    $self->{_discarded_db} = $discarded_db;
  }

  return $self->{_discarded_db};
}


sub ccds_db {
  my ($self, $ccds_db) = @_;
  if ($ccds_db) {
    $self->{_ccds_db} = $ccds_db;
  }

  return $self->{_ccds_db};
}

############################################################

sub havana_Coding_Genes {
  my ( $self, $g ) = @_;
  my @genes;
  if ( defined($g) ) {
    @genes = @$g;
  }

  if ( !defined( $self->{_havana_Coding_Genes} ) ) {
    $self->{_havana_Coding_Genes} = [];
  }

  if ( scalar @genes > 0 ) {
    push( @{ $self->{_havana_Coding_Genes} }, @genes );
  }

  return $self->{_havana_Coding_Genes};
}

sub havana_Processed_Transcript_Genes {
  my ( $self, $g ) = @_;
  my @genes;
  if ( defined($g) ) {
    @genes = @$g;
  }

  if ( !defined( $self->{_havana_Processed_Transcript_Genes} ) ) {
    $self->{_havana_Processed_Transcript_Genes} = [];
  }

  if ( scalar @genes > 0 ) {
    push( @{ $self->{_havana_Processed_Transcript_Genes} }, @genes );
  }

  return $self->{_havana_Processed_Transcript_Genes};
}

sub havana_Pseudogene_Genes {
  my ( $self, $g ) = @_;
  my @genes;
  if ( defined($g) ) {
    @genes = @$g;
  }

  if ( !defined( $self->{_havana_Pseudogene_Genes} ) ) {
    $self->{_havana_Pseudogene_Genes} = [];
  }

  if ( scalar @genes > 0 ) {
    push( @{ $self->{_havana_Pseudogene_Genes} }, @genes );
  }

  return $self->{_havana_Pseudogene_Genes};
}

sub combined_Transcripts {
  my ( $self, $t ) = @_;
  my @transcripts;

  if ( defined($t) ) {
    @transcripts = @$t;
  }

  if ( !defined( $self->{_coding_transcripts} ) ) {
    $self->{_coding_transcripts} = [];
  }

  if ( scalar @transcripts > 0 ) {
    push( @{ $self->{_coding_transcripts} }, @transcripts );
  }

  return $self->{_coding_transcripts};
}

sub combined_Processed_Transcripts {
  my ( $self, $t ) = @_;
  my @transcripts;

  if ( defined($t) ) {
    @transcripts = @$t;
  }

  if ( !defined( $self->{_processed_transcripts} ) ) {
    $self->{_processed_transcripts} = [];
  }

  if ( scalar @transcripts > 0 ) {
    push( @{ $self->{_processed_transcripts} }, @transcripts );
  }

  return $self->{_processed_transcripts};
}


sub combined_PseudoTranscripts {
  my ( $self, $pt ) = @_;
  my @pseudotranscripts;

  if ( defined($pt) ) {
    @pseudotranscripts = @$pt;
  }

  if ( !defined( $self->{_pseudo_transcripts} ) ) {
    $self->{_pseudo_transcripts} = [];
  }

  if ( scalar @pseudotranscripts > 0 ) {
    push( @{ $self->{_pseudo_transcripts} }, @pseudotranscripts );
  }
  return $self->{_pseudo_transcripts};
}


=head2 final_genes

 Descripton: this holds/returns the final genes produced after clustering transcripts and sharing common exons

=cut

sub final_genes {
  my ( $self, $g ) = @_;
  my @genes;
  if ( defined($g) ) {
    @genes = @$g;
  }

  if (@genes) {
    push( @{ $self->{_final_genes} }, @genes );
  }
  return $self->{_final_genes};
}

############################################################

=head2 gene_types

 Description: get/set for the type(s) of genes (usually TGE_gw, similarity_genewise and combined_e2g genes) 
              to be used in the genebuilder they get set in new()
              Does not include the ab inition predictions
=cut

sub gene_types {
  my ( $self, $type ) = @_;

  if ( defined($type) ) {
    push( @{ $self->{_gene_types} }, $type );
  }

  return $self->{_gene_types};
}

############################################################

sub features {
  my ( $self, $f ) = @_;
  my @features = @$f;

  if ( !defined( $self->{_feature} ) ) {
    $self->{_feature} = [];
  }
  if ( scalar @features ) {
    push( @{ $self->{_feature} }, @features );
  }
  return $self->{_feature};
}

############################################################

sub query {
  my ( $self, $slice ) = @_;

  if ( defined($slice) ) {
    $self->{_query} = $slice;
  }
  return $self->{_query};
}

#fetches sequence from appropriate database

sub fetch_sequence {
  my ( $self, $name, $db ) = @_;

  my $sa = $db->get_SliceAdaptor;

  my $slice = $sa->fetch_by_name($name);

  return $slice;
}


sub sort_clusters {
  my ($self, @matching_clusters) = @_;

  my @start_end_mc;
  for ( my $i = 0 ; $i < scalar(@matching_clusters) ; $i++ ) {
    my $mc = $matching_clusters[$i];
    my $start;
    my $end;
    foreach my $mt (@$mc) {
      if ( !defined($start) || $mt->start < $start ) {
        $start = $mt->start;
      }
      if ( !defined($end) || $mt->end > $end ) {
        $end = $mt->end;
      }
    }
    my @arr = ( $start, $end, $mc );
    $start_end_mc[$i] = \@arr;
    #print "start_end_mc[$i] before  = " . $start_end_mc[$i]->[0] . " " . $start_end_mc[$i]->[1] . " " . $start_end_mc[$i]->[2] . "\n";
  }

  @start_end_mc =
    sort { $a->[0] <=> $b->[0] ? $a->[0] <=> $b->[0] : $b->[1] <=> $a->[1] }
    @start_end_mc;

  my @sorted_clusters;
  foreach my $sem (@start_end_mc) {
    #print "sem = " . $sem->[0] . " " . $sem->[1] . " " . $sem->[2] . "\n";
    push @sorted_clusters, $sem->[2];
  }
  @matching_clusters = @sorted_clusters;
  return @matching_clusters;
} ## end sub sort_clusters


############################################################

=head2 cluster_into_Genes

  Arg [1]    : Listref of Bio::EnsEMBL::Transcript objects
  Arg [2]    : Listref of Bio::EnsEMBL::Transcript objects
  Arg [3]    : (optional) Int - flag for indicating transcripts are coding
  Example    : my @genes = @{ cluster_into_Genes(\@ensembl_transcripts, \@havana_transcripts ) };
  Description: Clusters ensembl transcripts and havana genes according
               It will take care of difficult cases like transcripts
               within introns. It also unify exons that are shared
               among transcripts.

              By default, the transcripts are assumed to be
              non-coding, either pseudogenes or processed transcripts.
              If clustering coding transcripts, set the $coding flag
              to one.
    Returns : Listref of Bio::EnsEMBL::Gene objects

=cut

sub cluster_into_Genes {
  my ( $self, $transcripts_unsorted, $havana_genes, $coding ) = @_;

  my $num_trans = scalar( @{$transcripts_unsorted} );
  my %ottg_xref;

  my @transcripts;
  if ($coding) {
    @transcripts = @{ $self->sort_transcripts( $transcripts_unsorted, 1 ) };
  } else {
    @transcripts = @{ $self->sort_transcripts($transcripts_unsorted) };
  }

  # We are going to set the havana genes as initial clusters as we
  # want to keep their structure.
  my @clusters;

  my @read_thru_trans;
  my @pure_read_thru_genes;
  print "Havana coding gene cluster size: "
    . scalar( @{$havana_genes} ) . "\n";
  foreach my $hav_gene ( @{$havana_genes} ) {

    #print "\$hav_gene->dbID " . $hav_gene->dbID . "\tbiotype: " . $hav_gene->biotype . "\n";

    my $ottg_key;
    foreach my $entry ( @{ $hav_gene->get_all_DBEntries } ) {
      #print "ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
      if ( $entry->dbname eq 'Vega_gene' ) {
        if ( $entry->primary_id eq $entry->display_id ) {
          $ottg_key = $entry->primary_id;
        }
      }
    }

    my @h_clust;
    my $read_thru;
    my $read_through_trans = 0;
    foreach my $hav_trans ( @{ $hav_gene->get_all_Transcripts } ) {
      #print "Havana trans id: ". $hav_trans->dbID . "\n";
      $ottg_xref{$hav_trans} = $ottg_key;

      # We want to keep track of the read_through transcripts. If a
      # gene has only read_through transcripts, then it's put aside
      # and not included in the merge process. The gene is then added
      # to the final cluster at the end.
      $read_thru = 0;
      my @hav_attrib = @{ $hav_trans->get_all_Attributes() };
      foreach my $attrib (@hav_attrib) {
        if ( $attrib->name eq 'readthrough transcript' ) {
          print "\nHave a trans with readthrough attrib: "
            . $hav_trans->dbID . "\n";
          $read_thru = 1;
          $read_through_trans++;
        }
      }
    }
    if ( $read_thru == 1 ) {

      # If the number of transcripts is the same as the read_through
      # counter then that gene is a pure read_through gene and should
      # be kept separate from rest of the Havana genes.
      if (
        $read_through_trans == scalar( @{ $hav_gene->get_all_Transcripts } ) )
      {
        #print "Pushing " . $hav_gene->dbID . " in pure set\n";
        my @prt_clust = @{ $hav_gene->get_all_Transcripts };
        push( @pure_read_thru_genes, \@prt_clust );
      } else {
        #print "Pushing " . $hav_gene->dbID . " in general set\n";
        push( @h_clust,  @{ $hav_gene->get_all_Transcripts } );
        push( @clusters, \@h_clust );
      }
    } else {
      push( @h_clust,  @{ $hav_gene->get_all_Transcripts } );
      push( @clusters, \@h_clust );
    }

  } ## end foreach my $hav_gene ( @{$havana_genes...

  print "There are " . scalar(@clusters) . " havana pseudogene clusters\n";
  print "There are " . scalar(@transcripts) . " ensembl transcripts\n";
  print "There are "
    . scalar(@pure_read_thru_genes)
    . " read_through havana gene clusters\n";

  # Clusters transcripts by whether or not any coding exon overlaps
  # with a coding exon in another transcript. We will use the set of
  # Havana.
  foreach my $tran (@transcripts) {
    print "Transcript: ", $tran->stable_id, " biotype ", $tran->biotype, "\n";

    my @matching_clusters;
  CLUSTER:
    foreach my $cluster (@clusters) {

      #print "havana transcript: ", $cluster->[0]->stable_id, "\t"
      #     . "dbID: ", $cluster->[0]->dbID, "\n";
      if ($coding) {
      CLUSTER_TRANS:
        foreach my $cluster_transcript (@$cluster) {
          #print "\$cluster_transcript->dbID "
          #  . $cluster_transcript->dbID
          #  . "\t\$cluster_transcript->biotype "
          #  . $cluster_transcript->biotype . "\n";

          if ( $cluster_transcript->translation ) {
            if ( $tran->coding_region_end >=
                    $cluster_transcript->coding_region_start
                 && $tran->coding_region_start <=
                 $cluster_transcript->coding_region_end )
            {

              my $exons1 = $self->get_coding_exons_for_transcript($tran);
              my $cluster_exons =
                $self->get_coding_exons_for_transcript($cluster_transcript);

              foreach my $exon1 ( @{$exons1} ) {
                foreach my $cluster_exon ( @{$cluster_exons} ) {

                  if (    $exon1->overlaps($cluster_exon)
                       && $exon1->strand == $cluster_exon->strand )
                  {
                    push( @matching_clusters, $cluster );
                    next CLUSTER;
                  }
                }
              }
            }
          }
        } ## end foreach my $cluster_transcript...
      } else {
        # If clustering pseudogenes or processed transcripts

        foreach my $cluster_transcript (@$cluster) {
          #print "\$cluster_transcript->dbID " . $cluster_transcript->dbID
          #    . "\t\$cluster_transcript->biotype " . $cluster_transcript->biotype . "\n";

          if (    $tran->end >= $cluster_transcript->start
               && $tran->start <= $cluster_transcript->end )
          {

            my $exons1        = $tran->get_all_Exons();
            my $cluster_exons = $cluster_transcript->get_all_Exons();

            foreach my $exon1 ( @{$exons1} ) {
              foreach my $cluster_exon ( @{$cluster_exons} ) {

                if (    $exon1->overlaps($cluster_exon)
                     && $exon1->strand == $cluster_exon->strand )
                {
                  push( @matching_clusters, $cluster );
                  next CLUSTER;
                }
              }
            }
          }
        }
      } ## end else [ if ($coding)
    } ## end foreach my $cluster (@clusters)

    if ( scalar(@matching_clusters) == 0 ) {
      print "\nNo matching clusters, new cluster being created with trans: "
        . $tran->dbID . "\n";

      my @newcluster;

      push( @newcluster, $tran );
      push( @clusters,   \@newcluster );
      print "\nSize of all clusters: " . scalar(@clusters) . "\n";
    } elsif ( scalar(@matching_clusters) == 1 ) {
      print "\nHave one trans in the matching clusters, "
        . $tran->dbID
        . " adding it on to the first cluster in matching_clusters\n";

      print "transcript_id ", $matching_clusters[0]->[0]->dbID, "\n";
      print "adding this one to the above: ", $tran->dbID, "\n";

      push @{ $matching_clusters[0] }, $tran;
    } else {
      print "\nBEWARE YOU ARE HERE MERGING CLUSTERS\n";
      print "Have >1 trans in the matching clusters:\n";

      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my @merged_cluster;
      my $hav_gene_counter = 0;

      # Check how many havana genes are in the cluster
    CLUST: foreach my $matching_cluster (@matching_clusters) {

        foreach my $clust_trans ( @{$matching_cluster} ) {
          print "\tPseudogene cluster -  transcript dbID: "
            . $clust_trans->dbID . " ("
            . $clust_trans->biotype . ")\n";
          if ( $clust_trans->biotype =~ /hav/ ) {
            $hav_gene_counter++;
            next CLUST;
          }
        }
      }

      # Merge all clustered genes in case there is only one or none havana genes.

      #print "\$hav_gene_counter: $hav_gene_counter\n";
      if ( $hav_gene_counter <= 1 ) {
        foreach my $cluster_to_be_merged (@matching_clusters) {
          push @merged_cluster, @{$cluster_to_be_merged};
        }
        push @merged_cluster, $tran;

        # ADD SOME NICE PRINT OUTPUT
        print "MERGING CLUSTER COVERED BY TRANSCRIPT IN: ",
          $tran->seq_region_name, " - ", $tran->seq_region_start, " - ",
          $tran->seq_region_end, "\n";

        print "content: ", join( ' - ', @merged_cluster ), "\n";
      } else {

        # If there is more than one Havana gene we check which one has
        # a better overlap on the merger transcript and we also check
        # if the havana genes are partial and need to be merged.

        # Sort the clusters

        #print "Before call matching_clusters size = " .  scalar(@matching_clusters) . "\n";

        @matching_clusters = $self->sort_clusters(@matching_clusters);

        my $cluster_ref  = 0;
        my $best_cluster = 0;
        my $best_match   = 0;

        my $tran_length;
        if ($coding) {
          $tran->translate->length;
        } else {
          $tran_length = $tran->length;
        }

        foreach my $overlapping_cluster (@matching_clusters) {
          print "\nMore than one Havana gene overlaps. Scoring the best cluster \n";
          print "Trans id: "
            . $tran->dbID
            . "\tbiotype: "
            . $tran->biotype . "\n";

          my $match_percent =
            $self->clust_overlap_percent( $overlapping_cluster, $tran );

          if ( $match_percent > $best_match ) {
            $best_match   = $match_percent;
            $best_cluster = $cluster_ref;
          }
          $cluster_ref++;
        }
        print "Only one havana gene gave the best overlap\n\n";
        print "Havana trans id: "
          . $matching_clusters[$best_cluster]->[0]->dbID
          . "\tbiotype: "
          . $matching_clusters[$best_cluster]->[0]->biotype . "\n";

        push( @merged_cluster, @{ $matching_clusters[$best_cluster] } );
        push( @merged_cluster, $tran );
      } ## end else [ if ( $hav_gene_counter...

      push( @new_clusters, \@merged_cluster );

      # Add back non matching clusters
    MATCHING:
      foreach my $clust (@clusters) {

        foreach my $tran_clust (@$clust) {
          foreach my $m_clust (@merged_cluster) {
            if ( $tran_clust == $m_clust ) {
              next MATCHING;
            }
          }
        }
        push( @new_clusters, $clust );
      }
      @clusters = @new_clusters;
    } ## end else [ if ( scalar(@matching_clusters...
  } ## end foreach my $tran (@transcripts)

  #print "size of pure_read_thru_genes: " . scalar(@pure_read_thru_genes) .  "\n";
  #print "size of clusters: " . scalar(@clusters) .  "\n";

  # Adding the read_through genes that we don't want to merge in
  # to the final cluster.
  push( @clusters, @pure_read_thru_genes );

  # safety and sanity checks
  $self->check_Clusters( scalar(@transcripts), \@clusters );

  # Merge the clustered ensembl transcript with the Havana genes.

  # make and store genes
  #print STDERR scalar(@clusters)." created, turning them into genes...\n";
  my @genes;
  foreach my $cluster (@clusters) {
    my $count         = 0;
    my $ottg_added    = "";
    my $biotype_added = 0;
    my $gene_biotype  = "";
    my $gene          = new Bio::EnsEMBL::Gene;

    foreach my $transcript (@$cluster) {
      if ( $biotype_added == 0 && $transcript->biotype =~ /hav/ ) {
        $gene_biotype = $transcript->biotype;
      } elsif ( $biotype_added == 0 ) {
        $gene_biotype = $transcript->biotype;
      }
      #print "Transcript Stable ID: ",$transcript->dbID,"\n";
      #print "This is the transcript_biotype ",$transcript->biotype,"\n";
      #print "Transcript dbID: "
      #  . $transcript->dbID . "\t("
      #  . $transcript->biotype . ")\n";

      $gene->add_Transcript($transcript);
      if ( $ottg_xref{$transcript}
           && ( $ottg_added ne $ottg_xref{$transcript} ) )
      {
        # Need to add the OTTG as the gene stable_id otherwise it will
        # be undef.
        $gene->stable_id( $ottg_xref{$transcript} );
        $self->add_ottg_xref( $gene, $ottg_xref{$transcript} );
        $ottg_added = $ottg_xref{$transcript};
      }
    }
    $gene->biotype($gene_biotype);
    push( @genes, $gene );
  } ## end foreach my $cluster (@clusters)

  return \@genes;

} ## end sub cluster_into_Genes


sub sort_transcripts {
  my ( $self, $transcripts_unsorted, $coding ) = @_;

  my @transcripts;

  if ($coding) {
    # First clean the coding exon cache in case it has any exons stored
    # from previous called to the cluster_into_Genes function.
    $self->clear_coding_exons_cache;

    my @transcripts_unsorted_translation;
    foreach my $tran ( @{$transcripts_unsorted} ) {
      if ( $tran->translation ) {
        push( @transcripts_unsorted_translation, $tran );
      }
    }

    @transcripts = sort {
          $a->coding_region_start <=> $b->coding_region_start
        ? $a->coding_region_start <=> $b->coding_region_start
        : $b->coding_region_end <=> $a->coding_region_end
    } @transcripts_unsorted_translation;
  } else {
    @transcripts = sort {
          $a->start <=> $b->start
        ? $a->start <=> $b->start
        : $b->end <=> $a->end
    } @{$transcripts_unsorted};
  }
  return \@transcripts;
} ## end sub sort_transcripts


1;
