
=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

  Questions may also be sent to the Ensembl help desk
  at <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::HavanaAdder

=head1 SYNOPSIS

  This is the main analysis database

  my $genebuilder =
    new Bio::EnsEMBL::Analysis::Runnable::HavanaAdder(
                                                 '-slice'    => $self->query,
                                                 '-input_id' => $self->input_id,
    );


=head1 DESCRIPTION

  This module reads your favourite annotations (ensembl,
  protein_coding,...) on the one hand, and manually
  curated plus features on the other hand. The product of
  Bio::EnsEMBL::Analysis::Runnable::HavanaAdder is a combination of both
  annotations where redundant transcripts are eliminated. The resulting
  transcripts are combined into genes. For more details, follow the list
  of methods called by build_Genes() method and the description in each
  one.

  The rest of the documentation details each of the object
  methods. Internal methods are usually preceded with a _

=head1 METHODS

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
  HAVANA_INPUT_NONCODING_TYPE
  ENSEMBL_INPUT_NONCODING_TYPE
);
  #ENSEMBL_INPUT_PROCESSED_TYPE no need to add this

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

  Arg        : String (it expects a string of the format chr_name.start_coord-end_coord)
  Description: Getter/setter for input id
  Returns    : String

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

  Arg        : none
  Example    : my @genes = $self->build_Genes
  Description: Builds genes. It is like the run method in Runnables. It
               calls everything that needs to be done.
  Returns    : none
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

#  foreach my $g (@genes) {
#    foreach my $t ( @{ $g->get_all_Transcripts } ) {
#      print "DEBUG: trans: " . $t->dbID . "\t"
#        . "nr of tsf: " . scalar( @{ $t->get_all_supporting_features } ) . "\n";
#      foreach my $e ( @{ $t->get_all_Exons } ) {
#        print "DEBUG: nr esf for exon " . $e->dbID . " :: "
#          . scalar( @{ $e->get_all_supporting_features } ) . "\n";
#      }
#    }
#  }
  print STDERR scalar(@genes) . " genes built\n";

  $self->update_gene_biotypes(@genes);

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

      if ( $transcript->biotype =~ /_hav$/ ) {
        #print "DEBUG: I'm a havana transcript with biotype: ",
        #  $transcript->biotype . "\t" . $transcript->dbID . "\n";
        push( @havana, $transcript );
        next TRANSCRIPT;
      }

      print "\n\n\t\tDEBUG: in _merge_redundant_transcripts: "
        . "adding this transcript to ensembl array: "
        . $transcript->dbID . " ("
        . $transcript->biotype . ")\n";

      push( @ensembl, $transcript );
    }
    #print "DEBUG: in _merge_redundant_transcripts: "
    #    . "bfr checking the array size: " . scalar(@ensembl) . "\n";
    if ( !scalar(@havana) ) {
      next GENE;
    }

    if ( !scalar(@ensembl) ) {
      # As this is a havana only gene I still want to add an OTTT
      # entry to the transcripts and an OTTG to the gene.
      $self->add_ottg_xref( $gene, $gene->stable_id, $gene->version );

      foreach my $hav_trans (@havana) {
        $self->add_ottt_xref($hav_trans);
      }

      next GENE;
    }
    print "Havana tran: ",  scalar(@havana),  "\n";
    print "Ensembl tran: ", scalar(@ensembl), "\n";

    # Compare each havana transcript to each ensembl one
    foreach my $ht (@havana) {
      print "\n===> Looking at havana trans: " . $ht->dbID . "\n";

      foreach my $et (@ensembl) {

        print "\n===> Looking at ensembl trans: " . $et->dbID . "\n";

        my $delete_t = $self->are_matched_pair( $gene, $ht, $et );

        foreach my $ht_attrib (@{$ht->get_all_Attributes()}) {
          if ($ht_attrib->code eq 'gene_cluster') {
            print "hav transcript is in a gene_cluster, "
                . "will remove the ensembl transcript. delete_t = 2.\n";
            $delete_t = 2;
          }
        }


        my @t_pair;
        my $delete_trans = 0;

        # We check all possible match pairs and give preference to the
        # one that shares CDS and UTR. This was added so only the best
        # matching havana/ensembl pair is chosen and to avoid a one
        # to many link

        # delete_t can be 0 (don't merge, trans completely different),
        #                 1 (don't merge, UTR structures different),
        #               $et (delete Ens trans) or
        #               $ht (delete Hav trans)

        #print "DEBUG: \tdelete_t was $delete_t "
        #    . "when just returned from are_matched_pair.\n";

        # if delete_t is 1 (don't merge), et, ht or 2 (Ens is incomplete, remove)
        if ($delete_t) {
          #print "DEBUG: delete_t is: " . $delete_t . "\n";
          if ( $delete_t == $et ) {
            print "DEBUG: --->>> are_matched_pair returned Ens trans " . $et->dbID . 
                  " which will be merged to the Hav trans and then deleted.\n";
            $delete_trans = $et;
            @t_pair = ( $ht, $et );
          } elsif ( $delete_t == $ht ) {
            print "DEBUG: --->>> are_matched_pair returned Hav trans " . $ht->dbID .
                  " which will be merged to the Ens trans and then deleted.\n";
            $delete_trans = $ht;
            @t_pair = ( $ht, $et );
          } elsif ( $delete_t == 2) {
            print "DEBUG: --->>> are_matched_pair returned value 2, Ens trans ". $et->dbID .
                  " is incomplete and will be deleted without merging with Hav trans.\n";
            $delete_trans = $et;
          } else {
            # We get here when $delete_t = 1 (et and ht differ in UTR structure) 
            print "DEBUG: --->>> are_matched_pair returned value 1, "
                . "Ens trans and Hav trans differ in UTR structure, not merging.\n";
            $delete_trans = $delete_t;
            @t_pair = ( $ht, $et );
          }
        } else {
          print "DEBUG: --->>> are_matched_pair returned value 0, "
              . " Ens trans and Hav trans are completely different, not merging.\n";
          $delete_trans = $delete_t;
        }

        if ( $delete_trans && $delete_trans !=1 && $delete_t !=2) {

          print "DEBUG: --->>> _merge_redundant_transcripts - "
              . " in the tricky if: delete_t: " . $delete_t
              . "\tdelete_trans: " . $delete_trans . "\n";

          # When delete_t is NOT 0, 1 or 2, i.e. it is either et or
          # ht, we want to merge sth)

          # For each pair of Ens + Hav transcripts we merge, we need to:
          #
          # (1) transfer supporting features from the to-be-delete trans onto the kept trans
          #     and setting correct xrefs using the "set_transcript_relation" method;
          # (2) add a biotype suffix to merged transcripts to mark them as merged;
          # (3) finally delete the Ens trans OR Hav trans from the pair as appropriate;

          # 1. Transfer supporting features:
          $self->set_transcript_relation( $delete_trans, @t_pair );

          # 2. Flag merged transcripts with biotype suffix if they
          #    haven't got any yet:
          my $new_bt_0;    #new biotype string for merged Hav trans
          my $new_bt_1;    #new biotype string for merged Ens trans

          # biotype suffix for Hav
          if ( $t_pair[0]->biotype !~ /$MERGED_TRANSCRIPT_OUTPUT_TYPE$/ ) {
            $new_bt_0 = $t_pair[0]->biotype . $MERGED_TRANSCRIPT_OUTPUT_TYPE;
            $t_pair[0]->biotype($new_bt_0);
          }

          # biotype suffix for Ens
          if ( $t_pair[1]->biotype !~ /$MERGED_TRANSCRIPT_OUTPUT_TYPE$/ ) {
            $new_bt_1 = $t_pair[1]->biotype . $MERGED_TRANSCRIPT_OUTPUT_TYPE;
            $t_pair[1]->biotype($new_bt_1);
          }

          # 3. Finally remove the Ens or Hav transcript which has been
          #    merged into its counterpart:
          $self->_remove_transcript_from_gene( $gene, $delete_trans )
            unless $delete_trans == 1;
        }
        else {
          #print "DEBUG: --->>> _merge_redundant_transcripts - "
          #    . "in the tricky else: delete_t: " . $delete_t
          #    . "\tdelete_trans: " . $delete_trans . "\n";

          # delete_t = 0 means Ens and Hav trans are completely different, hence no merge.
          # delete_t = 1 means Ens and Hav trans share the same CDS but differ in UTRs
          # delete_t = 2 means we need to remove Ens trans because its CDS is incomplete

          # When delete_t is 0, 1, or 2, we do not merge transcripts. 

          # For all three cases, we add OTTT xref to the Havanna
          # transcript. And we do nothing else for delete_t = 0.

          # If delete_t = 1, we want to set xref connections between the Ens and
          # Hav transcripts by calling set_transcript_relation to show they share 
          # the same coding region.
 
          # If delete_t = 2, we do not merge but we ned to delete Ens
          # transcript.

          $self->add_ottt_xref($ht);

          if ($delete_t == 1) {
            $self->set_transcript_relation( $delete_t, @t_pair );
          } 
          if ($delete_t == 2) {
            $self->_remove_transcript_from_gene( $gene, $et );
          }
        }
      }
    } ## end foreach my $ht (@havana)
    $gene->recalculate_coordinates;
  } ## end foreach my $gene (@$genes)
} ## end sub _merge_redundant_transcripts

sub gene_has_assembly_error_attribute {

  my $gene = shift;

  my @attribs = @{$gene->get_all_Attributes('NoTransRefError')};
  return (scalar(@attribs) > 0);
}

sub are_matched_pair {

  # This method checks pairs of Ens + Hav transcripts and returns 4
  # different possible values:
  # return 0 means keep both transcripts as they have different coding
  #          region or different exon structure
  # return 1 means keep both as they have same coding but different
  #          UTR exon structure
  # return ($ensembl) means keep havana transcript and remove ensembl 
  # return ($havana) means keep ensembl transcript and remove hanana

  my ( $self, $gene, $havana, $ensembl ) = @_;

  # Fetch all exons in each transcript
  my @hexons = @{ $havana->get_all_Exons };
  my @eexons = @{ $ensembl->get_all_Exons };

######## DEBUGGING TESTS ############
#  my $ens_sf = 0;
#  foreach my $ens_e (@eexons) {
#    print "DEBUG (are_matched_pair): " . $ens_e->dbID . " :: "
#        . scalar( @{ $ens_e->get_all_supporting_features } ) . "\n";
#    $ens_sf += scalar(@{ $ens_e->get_all_supporting_features } );
#  }
#
#  my $hav_sf = 0;
#  foreach my $hav_e (@hexons) {
#    $hav_sf += scalar(@{ $hav_e->get_all_supporting_features } );
#  }
######## DEBUGGING TESTS ############

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

  #  print "no. of ens_esf: " . $ens_sf . "\n";
  #  print "no. of hav_esf: " . $hav_sf . "\n";

  ###
  ### Before we start comparing each pair of Ensembl and Havana transcripts,
  ### we tag Ensembl models which are (1) non-coding ($non_coding_e = 1) and
  ### (2) non-CCDS and incomplete ($non_coding_e = 2). 
  ### All other Ensembl transcripts will have $non_coding_e set at "0".
  ###

  my $complete_5        = 0;
  my $incomplete        = 0;
  my $full_length_cdnas = 0;

  my $ccds_check_result =  $self->check_transcript_in_external_db('ccds', $ensembl );

  if ( @{ $ensembl->get_all_translateable_Exons } &&  $ccds_check_result == 1 ) {
    my $cdna         = $ensembl->spliced_seq;
    my $coding_start = $ensembl->cdna_coding_start;
    my $coding_end   = $ensembl->cdna_coding_end;
    my $cdna_start   = substr( $cdna, $coding_start - 1, 3 );
    my $cdna_end     = substr( $cdna, $coding_end - 3, 3 );

    if ( $cdna_start eq 'ATG' ) {
      print "The ensembl transcript starts with ATG\n";

      if ( $cdna_end eq 'TAG' || $cdna_end eq 'TGA' || $cdna_end eq 'TAA' ) {
        print "Ensembl transcript has start and stop codon: " .  $ensembl->dbID . "\n";
        $full_length_cdnas++;
        @teexons = @{ $ensembl->get_all_translateable_Exons };
      } else {
        # Flag the transcripts with start but no end for later deletion.
        $non_coding_e = 2;
        $complete_5++;
        #print "DEBUG: Have a MET but not a stop:  " .  $ensembl->dbID . "\n";
      }
    } else {
      $non_coding_e = 2;
      $incomplete++;
      #print "DEBUG: Not full-length transcript to delete: " . $ensembl->dbID . "\n";
    }
  } elsif ( @{ $ensembl->get_all_translateable_Exons } && $ccds_check_result == 0 ) {
    print "\nEnsembl transcript is coding and is part of CCDS. Not checking cds completeness.\n";
    @teexons = @{ $ensembl->get_all_translateable_Exons };
  } else {
    $non_coding_e = 1;
    print "\nNon-coding ensembl trans id " . $ensembl->dbID . "\n";
  }

  ###
  ### Now check the pairs!
  ###

  # Check that the number of exons is the same in both transcripts
  if (scalar(@hexons) != scalar(@eexons)) {
    print "Different number of exons in ensembl and havana transcripts, "
        . "skipping the comparison\n";
    return 0;
  }

  my %e_input_noncoding_type;
  foreach my $noncoding_type ( @{$ENSEMBL_INPUT_NONCODING_TYPE} ) {
    # note only pure non_coding (not including pseudo)
    $e_input_noncoding_type{$noncoding_type} = 1;
  }

  if ( $non_coding_h == 1 && $non_coding_e == 1 ) {
    print "\n===>>> ensembl and havana are both non-coding <<<===\n";
    # CASE 4:
    # if both of them are single-exon transcripts
    # and the Ensembl transcript is in ENSEMBL_INPUT_NONCODING_TYPE (ncRNAs)
    # and their length is the same,
    # we keep the Ensembl transcript (because we want to keep Ensembl ncRNA biotypes)
    if ( (scalar(@hexons) == 1) # already checked that the number of exons is the same
      && ($e_input_noncoding_type{$ensembl->biotype})
      && ($ensembl->length == $havana->length) ) {
      print "NON-CODING CASE 4 - single-exon transcripts, same length, Ensembl ncRNA (keep Ensembl)";
      return $havana;

    # We check two non coding transcripts. If they have the same
    # structure we keep the one from havana but if the have same exon
    # structure but one is slightly longer we take the longest one of
    # the pair.
    } elsif ( !$self->check_internal_exon_structure( \@eexons, \@hexons ) ) {
        #print "DEBUG: value is: " . $self->check_internal_exon_structure( \@eexons, \@hexons ) . "\n";
        print "NON-CODING CASE 0 - BOTH selected\n";

        # CASE 0: the two transcripts have different internal exon structure
        return 0;
    }
    # CASE 1: Havana is longer or both are exactly the same
    print "NON-CODING CASE 1 - Havana longer or both are the same length\n";
    if ( $self->check_terminal_exon_structure( \@hexons, \@eexons ) ) {
      #print "DEBUG: value is: " . $self->check_terminal_exon_structure( \@eexons, \@hexons ) . "\n";
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
      # Any transcript coming here should be deleted as we don't want
      # to have them in the merge set.

      print "Not full_length ensembl transcript, will be deleted: " . $ensembl->dbID . "\n";
      print "\tFull-length: $full_length_cdnas\n"
          . "\t5' complete: $complete_5\n"
          . "\tIncomplete: $incomplete\n";

      return 2;
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

        #print "DEBUG: value is: " . $self->check_internal_exon_structure( \@eexons, \@hexons ) . "\n";
        return 0;
      }

      # CASE 1: If the internal structure is the same we then keep the
      #         Havana one.
      print "HAV CODING vs ENS NON-CODING CASE 1 - "
          . "internal structure the same, keeing havana\n";

      return $ensembl;
    }

    if ( $non_coding_e == 0 ) {

      print "===>>> ensembl is coding, havana is non-coding <<<===\n";

      # Keep both if they don't have the same coding exon structure
      return 0 unless scalar(@hexons) == scalar(@eexons);

      # CASE: GENOME ASSEMBLY ERROR
      # Check if there is an havana gene attrib for the havana transcript that
      # means we should keep the ensembl protein_coding model because havana annotation
      # is based on a sequence with a genome assembly error.
      # Don't make Ensembl to non-coding and keep both transcripts.

      if (gene_has_assembly_error_attribute($gene)) {
        print "ENS CODING vs HAV NON-CODING CASE 'GENOME ASSEMBLY ERROR' - "
            . "keeping ensembl (and havana) regardless because\n";
        print "the Havana model " . $ensembl->dbID
            . " is based on a region where there is a genome assembly error.\n";
        return 0;
      }
      # Check that the ensembl transcript is not a CCDS model.

      # CASE 0: Ensembl transcript belongs to CCDS set, don't
      # make it to non-coding and keep both transcripts.

      # TODO: In the future, this CCDS check should refer to an
      # internal flagging system in this module, instead of referring
      # to the CCDS DB.

      if ( $self->check_transcript_in_external_db('ccds', $ensembl ) == 0 ) {
        print "ENS CODING vs HAV NON-CODING CASE 0 - "
            . "keeping ensembl regardless because\n";
        print "the Ensembl model " . $ensembl->dbID
            . " is a CCDS transcript will not convert it to non-coding.\n";
        return 0;
      }

      # If Ensembl model is not CCDS...

      # First check if the internal structure of the whole
      # transcripts is conserved.
      if ( scalar(@hexons) > 1 ) {
        unless ( $self->check_internal_exon_structure( \@hexons, \@eexons ) )
        {
          # CASE 1: The two transcripts have different internal exon
          #         structure we keep both
          print "ENS CODING vs HAV NON-CODING CASE 1 - "
              . "different internal exon structure, both selected\n";
          return 0;
        }
      }

      # Now we check if the coding bit is longer than the whole non
      # coding transcript.

      # CASE 2: The ensembl transcript is longer or equal but we still
      #         discard the Ensembl model. We keep the havana transcript
      #         regardless of length. Since the ensembl model is deleted, 
      #         there's *no* need to turn it into a non-coding model by
      #         changing its biotype or remove its translation.

      if ( $self->check_terminal_exon_structure( \@hexons, \@teexons ) ) {
        print "ENS CODING vs HAV NON-CODING CASE 2 - "
            . "removing ensembl model as it is equal or longer "
            . "than the Havana model but Havana models are prioritised.\n";

        #print "DEBUG: making ensembl model " . $ensembl->dbID . " to non-coding\n";

        return $ensembl;

      } else {

      # CASE 3: The havana transcript is longer at both ends so we
      #         remove the ensembl transcript. Again, there's no need to
      #         change the biotype of the model or remove its translation.

        print "ENS CODING vs HAV NON-CODING CASE 3 - "
            . "havana transcript is longer in both ends, "
            . "removing ensembl: " . $ensembl->dbID . "\n";

        #print "DEBUG: changing the biotype from: "
        #  . $ensembl->biotype . " to "
        #  . $havana->biotype . "_e\n";

        return $ensembl;
      }
    } ## end if ( $non_coding_e == ...

  } elsif ( $non_coding_h == 0 && $non_coding_e == 0 ) {

    return 0
      unless ( ( scalar(@teexons) == scalar(@thexons) )
               && $havana->translation->genomic_start == $ensembl->translation->genomic_start
               && $havana->translation->genomic_end   == $ensembl->translation->genomic_end );

    # Special case for single exon genes
    if ( scalar(@hexons) == 1 ) {
      #print "DEBUG: SINGLE EXONS!\n";

      if (    $hexons[0]->start  == $eexons[0]->start
           && $hexons[0]->end    == $eexons[0]->end
           && $hexons[0]->strand == $eexons[0]->strand
           && $thexons[0]->coding_region_start($havana) == $teexons[0]->coding_region_start($ensembl)
           && $thexons[0]->coding_region_end($havana)   == $teexons[0]->coding_region_end($ensembl) )
      {
        # Both are exactly the same so we delete the Ensembl one
        # unless the Ensembl one is already a merged one <-???? IS not doing
        # this here.
        return $ensembl;

      } elsif ($hexons[0]->start  <= $eexons[0]->start
            && $hexons[0]->end    >= $eexons[0]->end
            && $hexons[0]->strand == $eexons[0]->strand
            && $eexons[0]->start  == $teexons[0]->coding_region_start($ensembl)
            && $eexons[0]->end    == $teexons[0]->coding_region_end($ensembl) )
      {
        # Ensembl gene don't have UTR and Havana has then delete Ensembl one
        return $ensembl;

      }
      elsif ( ( (
                 $hexons[0]->start != $eexons[0]->start
              || $hexons[0]->end   != $eexons[0]->end )
            && $hexons[0]->strand  == $eexons[0]->strand )
          && ( $eexons[0]->start   != $teexons[0]->coding_region_start($ensembl)
              || $eexons[0]->end   != $teexons[0]->coding_region_end($ensembl) )
        )
      {
        # Both contain UTR keep ENSEMBL
        return $havana;

      } else {
        # We can be here when genes have different UTR start/end and
        # different CDS start/end or when the UTR start/end is the
        # same but the CDS start/end is different.

        #print "DEBUG: Keep both single exon genes\n";
        return 0;

      }
    } ## end if ( scalar(@hexons) ==...

    # if is a multi exons transcript
    else {
      # First we check the internal coding structure of the transcript
      # where everything has to be exactly equal

      #print "DEBUG: CHECKING INTERNAL EXONS \n";
      for ( my $i = 1 ; $i <= ( $#thexons - 1 ) ; $i++ ) {
        return 0
          unless (    $thexons[$i]->start == $teexons[$i]->start
                   && $thexons[$i]->end == $teexons[$i]->end
                   && $thexons[$i]->strand == $teexons[$i]->strand );
      }
      #print "DEBUG: INTERNAL CODING EXONS ARE OK \n";

      # Now check the rest of the internal exons that are not coding.
      # This is to find if the UTR exon structure is the same
      for ( my $i = 1 ; $i <= ( $#hexons - 1 ) ; $i++ ) {
        if (    $hexons[$i]->start != $eexons[$i]->start
             || $hexons[$i]->end != $eexons[$i]->end
             || $hexons[$i]->strand != $eexons[$i]->strand )
        {
          print "DEBUG: CASE WITH SAME CDS BUT DIFFERENT UTR STRUCTURE\n";
          print "DEBUG: HAVANA DIFF UTR BOUNDARIES: "
            . $havana->seq_region_name . " - "
            . $havana->seq_region_start . " - "
            . $havana->seq_region_end . "\n";
          print "DEBUG: ENSEMBL DIFF UTR BOUNDARIES: "
            . $ensembl->seq_region_name . " - "
            . $ensembl->seq_region_start . " - "
            . $ensembl->seq_region_end . "\n";
          return 1;
        }
      }
      #print "DEBUG: Internal utr exons are ok.\n";

      # Then check if the first an last exon are the same in both
      # transcripts. If just start and end of UTR are different keep
      # havana one.
      #
      # CASE 1: Both coding and UTR are the same, keep Havana and
      #         delete Ensembl
      if (    $hexons[0]->start == $eexons[0]->start
           && $hexons[0]->end == $eexons[0]->end
           && $hexons[0]->strand == $eexons[0]->strand
           && $hexons[-1]->start == $eexons[-1]->start
           && $hexons[-1]->end == $eexons[-1]->end
           && $hexons[-1]->strand == $eexons[-1]->strand )
      {
        #print "DEBUG: Multiexon delete ensembl\n";
        return $ensembl;

      } elsif
        ( # CASE 2: Havana has utr and ensembl doesnt, keep havana.
          #         Forward strand
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
        #print "DEBUG: Multiexon delete ensembl\n";
        return $ensembl;

      } elsif
        ( # CASE 3: Both ensembl and havana have utr but with
          #         different start/end, keep havana. Forward strand
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
        print "DEBUG: Case 3: MULTIEXON DELETE ENSEMBL\n";
        #print "DEBUG: CASE WITH DIFFERENT TERMINAL EXON BOUNDARIES\n";
        #print "DEBUG: HAVANA BOUNDARIES: ".$havana->seq_region_name." - ".$havana->seq_region_start. " - ".$havana->seq_region_end."\n";
        #print "DEBUG: ENSEMBL BOUNDARIES: ".$ensembl->seq_region_name." - ".$ensembl->seq_region_start. " - ".$ensembl->seq_region_end."\n";
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
        print "DEBUG: Case 4: MULTIEXON DELETE ENSEMBL\n";
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
        print "DEBUG: Case 5: different terminal exon boundaries, delete ens\n";
        #print "DEBUG: Havana boundaries: "
        #  . $havana->seq_region_name . " - "
        #  . $havana->seq_region_start . " - "
        #  . $havana->seq_region_end . "\n";
        #print "DEBUG: Ensembl boundaries: "
        #  . $ensembl->seq_region_name . " - "
        #  . $ensembl->seq_region_start . " - "
        #  . $ensembl->seq_region_end . "\n";
        return $ensembl;

      } else {
        print "DEBUG: \tShould I be here?\n";
        print "DEBUG: \tKeep multiexon both: havana "
          . $havana->dbID
          . " and ensembl: "
          . $ensembl->dbID . "\n";
        print "DEBUG: \tHavana boundaries: "
          . $havana->seq_region_name . " - "
          . $havana->seq_region_start . " - "
          . $havana->seq_region_end . "\n";
        print "DEBUG: \tEnsembl boundaries: "
          . $ensembl->seq_region_name . " - "
          . $ensembl->seq_region_start . " - "
          . $ensembl->seq_region_end . "\n";
        print "\n";
        return 1;
      }

    } ## end else [ if ( scalar(@hexons) ==...

    print "Weird case we did not think about, check rules!\n";
    return 0;
  } ## end elsif ( $non_coding_h == ...
} ## end sub are_matched_pair

=head2 check_internal_exon_structure

  Arg        : None
  Description: Check if the start and end of internal exon pairs in
               two sets of exons is the same
  Return     : Returns 0 if they are the same, returns 1 if they are
               differents

=cut

sub check_internal_exon_structure {

  my ( $self, $firstexons, $secondexons ) = @_;

  my @exons1 = @{$firstexons};
  my @exons2 = @{$secondexons};

  if ( scalar(@exons1) != scalar(@exons2) ) {
    print "NOOOOOOOOO WHAT HAPPENED HERE?\n";
    print "DEBUG: ", scalar(@exons1), "  vs. ", scalar(@exons2), "\n";
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

  Arg        : None
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

  my $found = 0;
  foreach my $entry ( @{ $ht->get_all_DBEntries } ) {
    if ( $entry->dbname eq 'Vega_transcript' ) {
      if ( $entry->primary_id eq $entry->display_id ) {

        #print "DEBUG: I am adding an OTTT xref to the transcript\n";
        #print "DEBUG: OTTT TO ADD: ",$entry->primary_id,"\n";
        my $xref_ottt =
          new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                     -display_id    => $ht->display_id,
                                     -priority      => 1,
                                     -xref_priority => 0,
                                     -version       => $entry->version,
                                     -release       => 1,
                                     -dbname        => 'OTTT' );

        $xref_ottt->status("XREF");
        $ht->add_DBEntry($xref_ottt);
        $found = 1;
      } else {
        if (!$found) {
          my $xref_ottt =
            new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                      -display_id    => $ht->display_id,
                                      -priority      => 1,
                                      -xref_priority => 0,
                                      -version       => $entry->version,
                                      -release       => 1,
                                      -dbname        => 'OTTT' );

          $xref_ottt->status("XREF");
          $ht->add_DBEntry($xref_ottt);
          warning("This entry doesn't have a display_id, it's from a patch. "
            . "Will set entry->primary_id and ht->display_id with "
            . "corresponding entry version\n"
            . $entry->primary_id . " version: " . $entry->version . "\n");
          $found = 1;
        }
      }
    }
  }
} ## end sub add_ottt_xref

sub add_ottg_xref {
  my ( $self, $hg, $ottg, $ottg_version ) = @_;

  #print "Creating OTTG: " . $ottg. "  " . $ottg_version. "\n";

  # Don't want to save the OTTG xref twice so check if it's already stored.
  foreach my $entry ( @{ $hg->get_all_DBEntries } ) {
    if ( $entry->dbname eq 'OTTG' && ( $entry->primary_id eq $ottg ) ) {
      return 0;
    }
  }
  #print "DEBUG: OTTG TO ADD: ", $ottg, "\n";
  my $xref_ottg =
    new Bio::EnsEMBL::DBEntry( -primary_id    => $ottg,
                               -display_id    => $ottg,
                               -priority      => 1,
                               -xref_priority => 0,
                               -version       => $ottg_version,
                               -release       => 1,
                               -dbname        => 'OTTG' );

  $xref_ottg->status("XREF");

  $hg->add_DBEntry($xref_ottg);
} ## end sub add_ottg_xref


sub set_transcript_relation {
  # $t_pair[0] is the havana transcript and $t_pair[1] is the ensembl transcript
  # $delete_t here is passed from $delete_t or $delete_trans from _merge_redundant_transcripts.

  # $delete_t here can be $et (delete Ens trans), 
  #                       $ht (delete Hav trans),
  #                         1 (keep both as CDS is the same but UTR structures differ)

  # print "DEBUG (set_transcript_relation)\n";
  my ( $self, $delete_t, @t_pair ) = @_;

  # If both share CDS and UTR is different in structure and number of
  # exons we still keep both, and we link them via Xref.
  if ( $delete_t == 1 ) {
    #print "DEBUG (set_transcript_relation) we're keeping both transcripts. "
    #    . "transferring only the xrefs inbtw\n";

    # transfer OTTT ID and/or ENST
    my $found = 0;
    foreach my $entry ( @{ $t_pair[0]->get_all_DBEntries } ) {
      if ( $entry->dbname eq 'Vega_transcript' ) {
        if ( $entry->primary_id eq $entry->display_id ) {

          my $newentry =
            new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                       -display_id    => $entry->display_id,
                                       -priority      => 1,
                                       -xref_priority => 0,
                                       -version       => $entry->version,
                                       -release       => 1,
                                       -dbname => 'shares_CDS_with_OTTT' );

          $newentry->status("XREF");
          $t_pair[1]->add_DBEntry($newentry);

          my $xref_ottt =
            new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                       -display_id    => $entry->display_id,
                                       -priority      => 1,
                                       -xref_priority => 0,
                                       -version       => $entry->version,
                                       -release       => 1,
                                       -dbname        => 'OTTT' );

          # print "DEBUG: OTTT xref to be added here\n";

          $xref_ottt->status("XREF");
          $t_pair[0]->add_DBEntry($xref_ottt);
          $found = 1;
        } else {
          if (!$found) {
            my $newentry =
              new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                        -display_id    => $entry->primary_id,,
                                        -priority      => 1,
                                        -xref_priority => 0,
                                        -version       => $entry->version,
                                        -release       => 1,
                                        -dbname => 'shares_CDS_with_OTTT' );

            $newentry->status("XREF");
            $t_pair[1]->add_DBEntry($newentry);

            my $xref_ottt =
              new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                        -display_id    => $entry->primary_id,,
                                        -priority      => 1,
                                        -xref_priority => 0,
                                        -version       => $entry->version,
                                        -release       => 1,
                                        -dbname        => 'OTTT' );

            $xref_ottt->status("XREF");
            $t_pair[0]->add_DBEntry($xref_ottt);
            $found = 1;
          }
        }
      } ## end if ( $entry->dbname eq...
    } ## end foreach my $entry ( @{ $t_pair...


    my $link_attrib =
      Bio::EnsEMBL::Attribute->new(
      -CODE => 'enst_link',
      -NAME => 'enst link',
      -DESCRIPTION => 'Code to link a OTTT with an ENST when they both share the CDS of ENST',
      -VALUE => $t_pair[1]->dbID );

    $t_pair[1]->add_Attributes($link_attrib);

    # print "DEBUG: OTTT TO ADD: ",$t_pair[0]->stable_id,"\n";

  } ## end if ( $delete_t == 1 ), where Ens and Hav trans share CDS but differ in UTR structure

  # If transcript to delete is havana we create an xref for the entry
  # saying that the transcript is CDS equal to ensembl.

  elsif ( $delete_t == $t_pair[0] ) {
    #print "DEBUG (set_transcript_relation) del trans is havana: "
    #    . $delete_t->dbID . "\n";

    # transfer OTTT ID and/or ENST
    my $found = 0;
    foreach my $entry ( @{ $t_pair[0]->get_all_DBEntries } ) {
      if ( $entry->dbname eq 'Vega_transcript' ) {
        if ( $entry->primary_id eq $entry->display_id ) {
          my $newentry =
            new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                       -display_id    => $entry->display_id,
                                       -priority      => 1,
                                       -xref_priority => 0,
                                       -version       => $entry->version,
                                       -release       => 1,
                                       -dbname => 'shares_CDS_with_OTTT' );

          $newentry->status("XREF");
          $t_pair[1]->add_DBEntry($newentry);
          $found = 1;
        } else {
          if (!$found) {
            my $newentry =
              new Bio::EnsEMBL::DBEntry( -primary_id    => $entry->primary_id,
                                        -display_id    => $entry->primary_id,,
                                        -priority      => 1,
                                        -xref_priority => 0,
                                        -version       => $entry->version,
                                        -release       => 1,
                                        -dbname => 'shares_CDS_with_OTTT' );

            $newentry->status("XREF");
            $t_pair[1]->add_DBEntry($newentry);
            $found = 1;
          }
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
    $self->transfer_exon_support($delete_t, $t_pair[1]);

  } # end else ( $delete_t == $tpair[0]), i.e. deleting Havana transcript
 
  elsif ( $delete_t == $t_pair[1] ) {
    #print "DEBUG (set_transcript_relation) del trans is ensembl: "
    #    . $delete_t->dbID . "\n";

    # If the transcript to delete is ENSEMBL we add an xref entry say
    # that both transcripts are exact matches (including UTR)
    my $found = 0;
    foreach my $entry ( @{ $t_pair[0]->get_all_DBEntries } ) {
      if ( $entry->dbname eq 'Vega_transcript' ) {
        if ( $entry->primary_id eq $entry->display_id ) {

          my $enstentry =
            new Bio::EnsEMBL::DBEntry(
                                     -primary_id    => $entry->primary_id,
                                     -display_id    => $entry->display_id,
                                     -version       => $entry->version,
                                     -release       => 1,
                                     -priority      => 1,
                                     -xref_priority => 0,
                                     -dbname => 'shares_CDS_and_UTR_with_OTTT'
            );
          $enstentry->status("XREF");
          $t_pair[0]->add_DBEntry($enstentry);
          $found = 1;
        } else {
          if (!$found) {
            my $enstentry =
              new Bio::EnsEMBL::DBEntry(
                                      -primary_id    => $entry->primary_id,
                                      -display_id    => $entry->display_id,
                                      -version       => $entry->version,
                                      -release       => 1,
                                      -priority      => 1,
                                      -xref_priority => 0,
                                      -dbname => 'shares_CDS_and_UTR_with_OTTT'
              );
            $enstentry->status("XREF");
            $t_pair[0]->add_DBEntry($enstentry);
            $found = 1;
          }
        }
      }
    }

    # Transfer the supporting features both for transcript and exon of
    # the transcript to delete to the transcript we keep
    $self->transfer_transcript_and_exon_supporting_features( $delete_t, $t_pair[0] );

  } ## end elsif ( $delete_t == $t_pair[1] ), i.e. deleting Ensembl transcript
} ## end sub set_transcript_relation


sub transfer_transcript_and_exon_supporting_features {
  my ( $self, $delete_t, $transcript ) = @_;

  my @exon_features;

  my @delete_tsf = @{ $delete_t->get_all_supporting_features };
  my @transcript_sf = @{ $transcript->get_all_supporting_features };

  #print "DEBUG (transfer_transcript_and_exon_supporting_features) "
  #    . "\# of tsf bfr addition: ",scalar(@transcript_sf),"\n";
  #print "DEBUG (transfer_transcript_and_exon_supporting_features) "
  #    . "and delete tsf: ", scalar(@delete_tsf),"\n";

  DTSF: foreach my $dtsf (@delete_tsf) {
    if ( !$dtsf->isa("Bio::EnsEMBL::FeaturePair") ) {
      next DTSF;
    }
    $transcript->add_supporting_features($dtsf);
  }

  #print "DEBUG (transfer_transcript_and_exon_supporting_features) "
  #  . "after adding tsf to trans: "
  #  . scalar( @{ $transcript->get_all_supporting_features } ) . "\n";

  #print "DEBUG (transfer_transcript_and_exon_supporting_features) about to transfer exon sf\n";
  $self->transfer_exon_support($delete_t, $transcript);

} ## end sub transfer_transcript_and_exon_supporting_features

sub _remove_transcript_from_gene {
  my ( $self, $gene, $trans_to_del ) = @_;

  #  foreach my $entry(@{ $gene->get_all_DBEntries}){
  #   print "DEBUG: ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
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

  Arg        : None
  Description: Retrieves ensembl and havana gene annotations with supporting evidence. 
  ReturnType : None, but $self->combined_Transcripts is filled

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
        #print "DEBUG: in get_Genes: coding genes " . $egene->stable_id . "\n";
        push( @genes, $egene );
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
        #print "DEBUG: in get_Genes: pseudogenes " . $epseudogene->stable_id . "\n";
        push( @pseudogenes, $epseudogene );
      }
    }
  }

  print STDERR "Fetching ensembl 'non-coding' genes\n";
  foreach my $e_noncoding_biotype ( @{$ENSEMBL_INPUT_NONCODING_TYPE} ) {
    foreach my $e_noncoding_gene ( @{ $ensemblslice->get_all_Genes_by_type($e_noncoding_biotype,undef,1) } ) {
      $e_noncoding_gene->load();
      push( @processedgenes, $e_noncoding_gene);
    }
  }

  # no need to add this yet since we don't have ensembl processed_transcript types
  #print STDERR "Fetching ensembl 'processed transcript' genes\n";
  #foreach my $eprocessedbiotype ( @{$ENSEMBL_INPUT_PROCESSED_TYPE} ) {
  #  foreach my $eprocessedgene ( @{ $ensemblslice->get_all_Genes_by_type($eprocessedbiotype,undef,1) } ) {
  #    $eprocessedgene->load();
  #    push( @processedgenes, $eprocessedgene );
  #  }
  #}

  print STDERR "Retrieved "
    . scalar(@genes)
    . " genes of types: "
    . join( ", ", @{$ENSEMBL_INPUT_CODING_TYPE} ) . "\n";

  print STDERR "Retrieved "
    . scalar(@pseudogenes)
    . " pseudogenes of types: "
    . join( ", ", @{$ENSEMBL_INPUT_PSEUDO_TYPE} ) . "\n";

  print STDERR "Retrieved "
    . scalar(@processedgenes)
    . " non-coding genes of types: "
    . join( ", ", @{$ENSEMBL_INPUT_NONCODING_TYPE} ) . "\n";
    #. " processed-transcript and non-coding genes of types: "
    #. join( ", ", @{$ENSEMBL_INPUT_PROCESSED_TYPE}, @{$ENSEMBL_INPUT_NONCODING_TYPE} )
    #. "\n";

  # Fetch Havana genes
  print STDERR "Fetching havana genes\n";
  foreach my $hbiotype ( @{$HAVANA_INPUT_CODING_TYPE} ) {
    foreach my $hgene ( @{ $havanaslice->get_all_Genes_by_type($hbiotype,undef,1) } ) {
      $hgene->load();

      # We change the biotype of the havana genes/transcripts as it
      # could happend to be the same as the ensembl ones
      #my $biotype = $hgene->biotype . "_hav";
      #$hgene->biotype($biotype);
      foreach my $htran ( @{ $hgene->get_all_Transcripts } ) {
        my $tbiotype = $htran->biotype . "_hav";
        $htran->biotype($tbiotype);
      }
      push( @hgenes, $hgene );
    }
  }

  print STDERR "Fetching havana 'non-coding' genes\n";
  foreach my $h_noncoding_biotype ( @{$HAVANA_INPUT_NONCODING_TYPE} ) {
    foreach my $h_noncoding_gene ( @{ $havanaslice->get_all_Genes_by_type($h_noncoding_biotype,undef,1) } ) {
      $h_noncoding_gene->load();
      # We change the biotype of the havana genes/transcripts as it
      # could happend to be the same as the ensembl ones
      #my $noncoding_biotype = $h_noncoding_gene->biotype . "_hav";
      #$h_noncoding_gene->biotype($noncoding_biotype);

      foreach my $h_noncoding_trans ( @{ $h_noncoding_gene->get_all_Transcripts } )
      {
        my $trans_noncoding_biotype = $h_noncoding_trans->biotype . "_hav";
        $h_noncoding_trans->biotype($trans_noncoding_biotype);
      }
      push( @hprocessedgenes, $h_noncoding_gene);
    }
  }

  print STDERR "Fetching havana 'processed transcript' genes\n";
  foreach my $hprocessedbiotype ( @{$HAVANA_INPUT_PROCESSED_TYPE} ) {
    foreach my $hprocessedgene ( @{ $havanaslice->get_all_Genes_by_type($hprocessedbiotype,undef,1) } ) {
      $hprocessedgene->load();
      # We change the biotype of the havana genes/transcripts as it
      # could happend to be the same as the ensembl ones
      #my $processedbiotype = $hprocessedgene->biotype . "_hav";
      #$hprocessedgene->biotype($processedbiotype);

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
      #my $biotype = $hpseudogene->biotype . "_hav";
      #$hpseudogene->biotype($biotype);
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
    . " processed-transcript and non-coding genes of types: "
    . join( ", ", @{$HAVANA_INPUT_PROCESSED_TYPE}, @{$HAVANA_INPUT_NONCODING_TYPE} )
    . "\n";

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
  @transcripts = @{ $self->flush_transcript_xrefs( \@genes ) };

  @processedtranscripts =
    @{ $self->flush_transcript_xrefs( \@processedgenes ) };

  @pseudotranscripts =
    @{ $self->flush_transcript_xrefs( \@pseudogenes ) };

  # Join all the gene set together
  print STDERR "Finished fetching genes\n";
  $self->combined_Transcripts( \@transcripts );
  $self->combined_Processed_Transcripts( \@processedtranscripts );
  $self->combined_PseudoTranscripts( \@pseudotranscripts );
} ## end sub get_Genes

sub flush_transcript_xrefs {
  my ( $self, $genes ) = @_;

  #print "DEBUG: Checking premerge gene status\n";

  my @transcripts;
  foreach my $gene (@$genes) {
  TRANSCRIPT:
    foreach my $tran ( @{ $gene->get_all_Transcripts } ) {

      $self->flush_xref($tran);

      # add transcript to the merging list
      push( @transcripts, $tran );

    } ## end foreach my $tran ( @{ $gene...
  } ## end foreach my $gene (@$genes)
  return \@transcripts;
} ## end sub check_merge_transcript_status

sub check_transcript_in_external_db {
  my ( $self, $dbname, $trans ) = @_;

  my @exons   = @{ $trans->get_all_Exons };
  my @t_exons = @{ $trans->get_all_translateable_Exons };

  #print "DEBUG: --->>> exons: " . scalar(@exons) . "\n";
  #print "DEBUG: --->>> t_exons: " . scalar(@t_exons) . "\n";

  my $ext_slice;

  # external_db = ccds_db
  my $ccds_db = $self->ccds_db();
  if ( defined($ccds_db) ) {
    $ext_slice =
      $ccds_db->get_SliceAdaptor()->fetch_by_region( 'toplevel',
                                     $trans->slice()->seq_region_name(),
                                     $trans->seq_region_start(),
                                     $trans->seq_region_end() );
  }
  else {
    return 1;
  }


EXT_GENE:
  foreach my $ext_gene ( @{ $ext_slice->get_all_Genes } ) {
  EXT_TRANS:
    foreach my $ext_trans ( @{ $ext_gene->get_all_Transcripts } ) {
      my @ext_exons = @{ $ext_trans->get_all_Exons };
      my @ext_t_exons = @{ $ext_trans->get_all_translateable_Exons };

      #print "DEBUG: comparing ccds: " . $ext_trans->stable_id()
      #    . " vs trans: " . $trans->dbID . " ("
      #    . $trans->stable_id() . ")\n";
      if (@t_exons) {
        #print "DEBUG: --->>> ccds exons: " . scalar(@ext_exons) . "\n";
        #print "DEBUG: --->>> trans t_exons: " . scalar(@t_exons) . "\n";

        if ( scalar(@t_exons) == scalar(@ext_exons) ) {
          #print "DEBUG: " . $trans->dbID . " :: " . $trans->stable_id() . "\n";
          for ( my $i = 0 ; $i < scalar(@t_exons) ; $i++ ) {

            # Work out Ens coding exon start and end positions and convert
            # them to genomic coords to compare to external (CCDS) transcript

            my $exon_start_in_cds = $t_exons[$i]->coding_region_start($trans);
            my $exon_end_in_cds = $t_exons[$i]->coding_region_end($trans);
            my @genomic_coords = $trans->cdna2genomic( $exon_start_in_cds, $exon_end_in_cds );

            #print "DEBUG: te_start in slice: " . $t_exons[$i]->coding_region_start($trans)
            #  . " te_end in slice: " . $t_exons[$i]->coding_region_end($trans) . "\n";
            #print "DEBUG te_start in seq_region: " . $genomic_coords[0]->start
            #  . " te_end in seq_region: " . $genomic_coords[0]->end ."\n";
            #print "DEBUG: ext exon start: " . $ext_exons[$i]->seq_region_start
            #  . " ext exon end: " . $ext_exons[$i]->seq_region_end . "\n";

            ### CAUTION!!! THE FOLLOWING COMPARISON ONLY WORKS IF HavanaAdder
            ### IS RUN ON A WHOLE CHROMOSOME, OR ELSE "coding_region_start/end"
            ### WILL NEVER BE THE SAME AS seq_region_start/end!

            if ( $genomic_coords[0]->start != $ext_exons[$i]->seq_region_start
                 || $t_exons[$i]->strand != $ext_exons[$i]->strand
                 || $genomic_coords[0]->end != $ext_exons[$i]->seq_region_end )

            {
              #print "DEBUG: number of translateable exons matched "
              #    . "but not all exon boundaries matched.  Check next CCDS model...\n";

              # CCDS is one-gene-one-transcript. "next EXT_GENE" is equivalent to
              # next EXT_GENE;
              next EXT_GENE;
            }
          }
          # print "DEBUG: \t--->>> transcript " . $trans->display_id
          #     . " found in ccds db\n";
          return 0;
        } else {
          # print "DEBUG: ccds db: number of (translatable) exons is "
          #    . "different between " .  $ext_trans->stable_id()
          #    . " and ". $trans->display_id . "\n";
          next EXT_GENE;
        }
      } else {
        if ( scalar(@exons) == scalar(@ext_exons) ) {
          for ( my $i = 0 ; $i < scalar(@exons) ; $i++ ) {
            #print "DEBUG: exon start: " . $exons[$i]->seq_region_start
            #  . " vs ext_exon start: " . $ext_exons[$i]->seq_region_start . "\n";
            if ( $exons[$i]->seq_region_start != $ext_exons[$i]->seq_region_start
                 || $exons[$i]->strand != $ext_exons[$i]->strand
                 || $exons[$i]->seq_region_end != $ext_exons[$i]->seq_region_end )
            {
              next EXT_TRANS;
            }
          }
          # If you are here means that both transcripts are the same.

          #print "DEBUG: transcript found in ccds db\n";
          return 0;
        } else {

          # If you enter here means that these two transcripts are not
          # the same.

          #print "DEBUG: ccds db: number of (non-coding) exons is different\n";
          next EXT_GENE;
        }
      } ## end else [ if (@t_exons)
    } ## end foreach my $ext_trans ( @{ ...
  } ## end foreach my $ext_gene ( @{ $ext_slice...

  # If we reach here it means that no transcript in the external db
  # is the same as our transcript so we keep it.

  #print "DEBUG: transcript not found in external db\n";
  return 1;
} ## end sub check_transcript_in_external_db

sub flush_xref {
  my ( $self, $transcript ) = @_;

  my @newxrefs;
  #print "DEBUG: THIS IS WHAT NEEDS EMPTYING: ",$transcript->get_all_DBEntries,"\n";
  foreach my $tran_xref ( @{ $transcript->get_all_DBEntries } ) {
    if (    $tran_xref->dbname ne "shares_CDS_and_UTR_with_OTTT"
         && $tran_xref->dbname ne "shares_CDS_with_OTTT"
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

  Arg        : Array of Bio::EnsEMBL::Transcript
  Description: It separates transcripts according to strand and then clusters 
               each set of transcripts by calling _cluster_Transcripts_by_genomic_range()
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

  Arg        : Array of Bio::EnsEMBL::Transcript
  Description: It clusters transcripts according to genomic overlap
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
    #print "\nDEBUG: In cluster ".($count+1)."\n"
    #    . "start: $cluster_starts[$count] end: $cluster_ends[$count]\n"
    #    . "comparing:\n";

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
    print "Pseudogene biotype: ", $pseudo_gene->biotype,
      "\tstart: " . $pseudo_gene->seq_region_start . "\n";

    foreach my $gene_attrib (@{$pseudo_gene->get_all_Attributes()}) {
      if ($gene_attrib->code eq 'hav_gene_type') {
        #print "Have a blessed gene, use it to merge if possible\n";
        $is_pseudo_havana = 1;
        last;
      }
    }

    if ( $pseudo_gene->biotype =~ /hav/ ) {
      $is_pseudo_havana = 1;
    }
    my $pseudo_status = 0;
  OVERLAP:
    foreach my $coding_gene (@coding_genes) {

      #print "DEBUG (combine_gene_clusters) "
      #    . "\$is_pseudo_havana = $is_pseudo_havana\n";

      # I need to add this line as I will be making some coding genes
      # UNDEF as I add them to pseudogenes.
      next OVERLAP unless $coding_gene;

      foreach my $trans (@{ $coding_gene->get_all_Transcripts }) {

        if ($trans->biotype =~ /hav/) {
          next OVERLAP;
        }

        #print "DEBUG: checking if trans " . $trans->dbID . " is CCDS\n";
        if ( $self->check_transcript_in_external_db('ccds', $trans) == 0 ) {
          print "Keeping ensembl transcript as it is CCDS " . $trans->dbID . "\n";

          # TODO:
          # Can add a flag here to keep track of the Ensembl
          # transcript that needs to be saved.
          #
          # The flagging system can be used later in the
          # are_matched_pair method when Ensembl CCDS transcripts
          # needs "protection" again. Currently the are_matched_pair
          # method checks against the CCDS DB to confirm the
          # protection status of models. Using the internal system
          # should speed things up.
          next OVERLAP;
        }
      }

      if ( $is_pseudo_havana == 1 && $coding_gene->biotype =~ /hav/ ) {
        #print "Jumping over havana coding gene! " . $coding_gene->dbID .
        #      " biotype: " . $coding_gene->biotype . "\n";
        next OVERLAP;
      }

      if (    $coding_gene->end >= $pseudo_gene->start
           && $coding_gene->start <= $pseudo_gene->end )
      {
        #print "DEBUG (combine_gene_clusters) "
        #    . "coding gene smaller: " . $coding_gene->biotype . "\n";
        #print "DEBUG (combine_gene_clusters) "
        #    . "\$is_pseudo_havana = $is_pseudo_havana\n";
        print "combine_gene_clusters\n";
        print "coding_gene dbID: ".$coding_gene->dbID."\n" if ($coding_gene->dbID);
        print "coding_gene start: ".$coding_gene->start."\n";
        print "coding_gene end: ".$coding_gene->end."\n";
        print "pseudo_gene dbID: ".$pseudo_gene->dbID."\n" if ($coding_gene->dbID);
        print "pseudo_gene start: ".$pseudo_gene->start."\n";
        print "pseudo_gene end: ".$pseudo_gene->end."\n";

        my $coding_length = $self->get_coding_length($coding_gene);
        print "coding length: $coding_length\n";
        my $cg_exons = $self->get_coding_exons_for_gene($coding_gene);
        print "coding gene exons (translateable): $cg_exons\n";
        my $pg_exons = $pseudo_gene->get_all_Exons();
        print "pseudo gene exons: $pg_exons\n";

        foreach my $cg_exon ( @{$cg_exons} ) {
          foreach my $pg_exon ( @{$pg_exons} ) {

            print "pseudo_gene exon start: ".$cg_exon->start."\n";
            print "pseudo_gene exon end: ".$cg_exon->end."\n";
            print "pseudo_gene exon start: ".$pg_exon->start."\n";
            print "pseudo_gene exon end: ".$pg_exon->end."\n";

            if (    $cg_exon->overlaps($pg_exon)
                 && $cg_exon->strand == $pg_exon->strand )
            {
              # Check if the overlap covers at least 10 percent of the
              # coding region of the longest transcript in the gene.
              print "cg_exon overlaps pg_exon\n";

              # NOTE!!! This check is a bit experimental.
              if ( $self->overlap_percent( $cg_exon, $pg_exon, $coding_length ) > 10 ) {
                if ( $is_pseudo_havana == 1 ) {
                  print "I'm looking into merging a pseudogene\n";

                  # As the pseudogene is Havana I will add the coding
                  # transcripts to the pseudo and remove the translation.
                  foreach my $c_transcript ( @{ $coding_gene->get_all_Transcripts } ) {
                    #print "\nDEBUG (combine_gene_clusters) "
                    #  . "changing coding trans to non-coding: "
                    #  . $c_transcript->dbID . " (" . $c_transcript->biotype . ")\n";

                    #print "\nDEBUG (combine_gene_clusters) "
                    #  . "removing the translation, adding trans to pseudogene, "
                    #  . "and changing the biotype to: "
                    #  . $pseudo_gene->biotype . "_ens\n";

                    # if the hav gene has NoTransRefError attrib, keep ens protein_coding biotype
                    if (gene_has_assembly_error_attribute($pseudo_gene)) {
                      print "But I found NoTransRefError gene attrib in Hav gene ".$pseudo_gene->stable_id." . Ens transcript " . $c_transcript->dbID . " will be kept as ".$c_transcript->biotype .  "_ens\n";
                      $c_transcript->biotype($c_transcript->biotype .  "_ens");
                    } else {
                      #print "And I did not found any assembly error attribute\n";
                      $c_transcript->{translation} = undef;
                      $c_transcript->biotype($pseudo_gene->biotype .  "_ens");
                    }
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
                    #print "DEBUG: \ncombine_gene_clusters --->>> "
                    #  . "adding pseudo transcript to the coding gene: "
                    #  . $p_transcript->dbID . "\n";
                    $coding_gene->add_Transcript($p_transcript);
                  }
                  $pseudo_status = 1;
                  #SMJS Need to look closer at this (used to be next OVERLAP;)
                  print "DEBUG: !!!!!!!!! Done noncoding-coding merge - jumping to next noncoding gene\n";
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


{
  my %coding_exon_cache;

  sub clear_coding_exons_cache {
    %coding_exon_cache = ();
  }

=head2 get_coding_exons_for_transcript

  Arg        : Bio::EnsEMBL::Transcript object
  Description: It returns the coding exons of a transcript and
               stores them in a hash to safe computer time.
  Returns    : An array ref than contain exon objects.
  Example    : my $exons = $self->get_coding_exons_for_transcript($transcript);

=cut
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

      $coding_exon_cache{$trans} = \@coding;
      return $coding_exon_cache{$trans};
    }
  }
}


##################################################################

=head2 get_coding_exons_for_gene

  Example    : my $exons1 = $self->get_coding_exons_for_gene($gene);
  Description: It returns the coding exons of a transcript and stores
               them in a hash to safe computer time
  Returns    : An ArrayRef than contain Exon objects.
  Args       : Bio::EnsEMBL::Transcript object

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
          #print "DEBUG: prune_Exons: found two copies of the same exon, "
          #    . "tranferring esf from exon to ini now\n";
          $self->transfer_supporting_evidence_between_exons($exon, $uni);
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

    #print "DEBUG: Uniq_tran sid: ",$tran->dbID,"\n";

  } ## end foreach my $tran ( @{ $gene...
  return $gene;
} ## end sub prune_Exons

############################################################

=head2 transfer_supporting_evidence_between_exons

 Example    : $self->transfer_supporting_evidence_between_exons($source_exon, $target_exon)
 Description: Transfers supporting evidence from source_exon to target_exon, 
              after checking the coordinates are sane and that the evidence is not already in place.
 Returns    : Nothing, but $target_exon has additional supporting evidence

=cut

sub transfer_supporting_evidence_between_exons {
  my ( $self, $source_exon, $target_exon ) = @_;

  my @target_sf = @{ $target_exon->get_all_supporting_features };

  # Keep track of features already transferred, so that we do not duplicate
  my %unique_evidence;
  my %hold_evidence;

  #print "DEBUG (transfer_supporting_evidence_between_exons): "
  #  . "bfr source exon sf: " . scalar( @{ $source_exon->get_all_supporting_features } )
  #  . " vs target exon sf: " . scalar( @{ $target_exon->get_all_supporting_features } ) . "\t\t";

SOURCE_FEAT:
  foreach my $feat ( @{ $source_exon->get_all_supporting_features } ) {
    if ( !$feat->isa("Bio::EnsEMBL::FeaturePair") ) {
      next SOURCE_FEAT;
    }

    # Skip duplicated evidence objects
    next SOURCE_FEAT if ( $unique_evidence{$feat} );

    # Skip duplicated evidence
    if ( $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }
         { $feat->hstart }{ $feat->hend } )
    {
      # Skipping duplicated evidence
      next SOURCE_FEAT;
    }

  TARGET_FEAT:
    foreach my $tsf (@target_sf) {

      # Skip comparison between source and target features if the
      # target one isn't a daf or paf feature pair
      next TARGET_FEAT unless $tsf->isa("Bio::EnsEMBL::FeaturePair");

      # If supp feat in source exon is identical to existing supp feat
      # in the target exon, we don't want to transfer source feat onto
      # target exon and there's no need to check other existing supp
      # feat in the target exon.
      if (    $feat->start    == $tsf->start
           && $feat->end      == $tsf->end
           && $feat->strand   == $tsf->strand
           && $feat->hseqname eq $tsf->hseqname
           && $feat->hstart   == $tsf->hstart
           && $feat->hend     == $tsf->hend )
      {
        next SOURCE_FEAT;
      }
    }
    #print "DEBUG (transfer_supporting_evidence_between_exons): from "
    #  . $source_exon->dbID . " to "
    #  . $target_exon->dbID . "\n";

    # I may need to add a paranoid check to see that no exons longer
    # than the current one are transferred
    $target_exon->add_supporting_features($feat);
    $unique_evidence{$feat} = 1;
    $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }
      { $feat->hstart }{ $feat->hend } = 1;
  } ## end foreach my $feat ( @{ $source_exon...

  #print "DEBUG: no. of target exon sf afr: "
  #    . scalar( @{ $target_exon->get_all_supporting_features } ) . "\n";

} ## end sub transfer_supporting_evidence_between_exons

############################################################################

=head2 update_gene_biotypes

  Example    : $self->update_gene_biotypes(@genes)
  Description: For each gene, it checks the biotypes of its associated merged 
               transcript(s) and then set the gene biotype to reflect the
               combination of transcript biotypes.
  Returns    : Nothing, but each gene object in the @genes array will have its
               biotype attribute set.

=cut

sub update_gene_biotypes {
  my ( $self, @genes ) = @_;

  my %pseudobiotypes;
  my %processedbiotypes;
  #my %noncoding_biotypes;
  my %coding_biotypes;

  foreach my $epb ( @{$ENSEMBL_INPUT_PSEUDO_TYPE} ) {
    $pseudobiotypes{$epb} = 1;
  }
  foreach my $hpb ( @{$HAVANA_INPUT_PSEUDO_TYPE} ) {
    $pseudobiotypes{ $hpb . "_hav" } = 1;
  }
  foreach my $hpb ( @{$HAVANA_INPUT_PROCESSED_TYPE} ) {
    $processedbiotypes{ $hpb . "_hav" } = 1;
  }
  #foreach my $hpb ( @{$HAVANA_INPUT_NONCODING_TYPE} ) {
  #  $noncoding_biotypes{ $hpb . "_hav" } = 1;
  #}
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
    #my $has_noncoding = 0;
    my $has_coding    = 0;

    # Keeping track of the Havana non-coding gene biotypes, we don't
    # want to change these so they are 'blessed'.
    my $blessed_biotype;

    foreach my $tran ( @{ $gene->get_all_Transcripts } ) {
      $trans_types{ $tran->biotype } = 1;

      foreach my $pseudobiotype ( keys %pseudobiotypes ) {
        if ( $tran->biotype =~ /$pseudobiotype/ ) {
          $has_pseudos = 1;
        }
      }
      foreach my $processedbiotype ( keys %processedbiotypes ) {
        if ( $tran->biotype =~ /$processedbiotype/ ) {
          $has_processed = 1;
        }
      }
      #foreach my $noncoding_biotype ( keys %noncoding_biotypes ) {
      #  if ( $tran->biotype =~ /$noncoding_biotype/ ) {
      #    $has_noncoding = 1;
      #  }
      #}
      foreach my $coding_biotype ( keys %coding_biotypes ) {
        if ( $tran->biotype =~ /$coding_biotype/ ) {
          $has_coding = 1;
        }
      }
    }

    foreach my $gene_attrib (@{$gene->get_all_Attributes()}) {
      if ($gene_attrib->code eq 'hav_gene_type') {
        foreach my $non_coding_type (@{$HAVANA_INPUT_NONCODING_TYPE}) {
          if ($non_coding_type eq $gene_attrib->value) {
            print "Have a blessed non-coding gene with biotype: "
              . $gene_attrib->value . " vs "
              . $non_coding_type
              . " will not change the biotype\n";
            $blessed_biotype = $gene_attrib->value;
          }
        }
      }
    }

    # MANAGE OUTPUT BIOTYPES BEFORE WRITING OUTPUT
    my $newbiotype;
    my $biotype_status;
    my $has_havana  = 0;
    my $has_ensembl = 0;
    my $has_merged  = 0;

    if ( $has_coding == 1 ) {
      $biotype_status = "protein_coding";
    } elsif ( $has_pseudos == 1 && $has_coding == 0 ) {
      $biotype_status = "pseudogene";
    } elsif ( $has_processed == 1 && $has_coding == 0 && $has_pseudos == 0 ) {
      $biotype_status = "processed_transcript";
    } elsif ( $has_processed == 1 && $blessed_biotype) {
      $biotype_status = "processed_transcript";
    } elsif ($blessed_biotype) {
      #print "DEBUG: update_gene_biotype: type to keep: $blessed_biotype\n";
      $biotype_status = $blessed_biotype;
    } else {
      print "ERROR: I should not really be here for gene biotype checks\n";
      $biotype_status = "weird_" . $gene->biotype;
    }

    foreach my $t_biotype ( keys %trans_types ) {
      # Be careful with merged transcript as they can have either
      # $MERGED_TRANSCRIPT_OUTPUT_TYPE suffix or hard-coded
      # "_ens" suffix. The latter comes from special cases
      # where an Ensembl coding gene overlaps with a Havana non-coding
      # gene, and at least 10% of the longest transcript of the
      # Ensembl gene overlaps with the Havana one. In these cases,
      # we turn the Ensembl coding gene into pseudogene, keep the
      # longest Ensembl transcript, strip off the Ensembl transcript's
      # translation and then add "_ens" to transcript.
      if (    $t_biotype =~ /$MERGED_TRANSCRIPT_OUTPUT_TYPE$/
           || $t_biotype =~ /_ens$/ )
      {
        $has_merged = 1;
      } elsif ( $t_biotype =~ /_hav$/ ) {
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
} ## end sub update_gene_biotypes

############################################################


############################################################
#
# GETSET METHODS
#
############################################################

# get/set method holding a reference to the db with genewise and combined genes,
# havana genes
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

 Descripton: This holds/returns the final genes produced after
             clustering transcripts and sharing common exons

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

  Description: Get/set for the type(s) of genes (usually TGE_gw,
               similarity_genewise and combined_e2g genes) to be used in the
               genebuilder they get set in new() Does not include the ab inition
               predictions

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
    #print "DEBUG: start_end_mc[$i] before  = " . $start_end_mc[$i]->[0]
    #    . " " . $start_end_mc[$i]->[1] . " " . $start_end_mc[$i]->[2] . "\n";
  }

  @start_end_mc =
    sort { $a->[0] <=> $b->[0] ? $a->[0] <=> $b->[0] : $b->[1] <=> $a->[1] }
    @start_end_mc;

  my @sorted_clusters;
  foreach my $sem (@start_end_mc) {
    #print "DEBUG: sem = " . $sem->[0] . " " . $sem->[1] . " " . $sem->[2] . "\n";
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
  Description: Clusters ensembl transcripts and havana genes according to exon-overlap.
               It will take care of difficult cases like transcripts
               within introns. It also unify exons that are shared
               among transcripts.

               By default, the transcripts are assumed to be
               non-coding, either pseudogenes or processed transcripts.
               If clustering coding transcripts, set the $coding flag
               to one.
  Returns    : Listref of Bio::EnsEMBL::Gene objects

=cut

sub cluster_into_Genes {
  my ( $self, $transcripts_unsorted, $havana_genes, $coding ) = @_;

  my $num_trans = scalar( @{$transcripts_unsorted} );
  my %ottg_xref;
  my %ottg_type;
  my %gene_attribs_to_keep;

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

  print "Havana";
  print " coding" if ($coding); 
  print " non-coding" if (!$coding);
  print " gene cluster size: "
    . scalar( @{$havana_genes} ) . "\n";
  foreach my $hav_gene ( @{$havana_genes} ) {

    #print "DEBUG: \$hav_gene->dbID " . $hav_gene->dbID
    #    . "\tbiotype: " . $hav_gene->biotype . " " . $hav_gene->stable_id."\n";

    my $hav_stable_id = $hav_gene->stable_id();
    my $hav_type = $hav_gene->biotype();
    $ottg_type{$hav_stable_id} = $hav_type;

    my @hav_gene_attrib = @{ $hav_gene->get_all_Attributes() };

    # Want to transfer the ncrna_host and NoTransRefError gene attributes from vega gene
    # to merged gene.
    foreach my $gene_attrib (@hav_gene_attrib) {
      if ( ($gene_attrib->name() eq 'ncrna_host') or
           ($gene_attrib->code() eq 'NoTransRefError') ) {
        print "Found gene_attrib to keep: ". $gene_attrib->code() ."\n";
        push @{$gene_attribs_to_keep{$hav_stable_id}},$gene_attrib;
      }
    }
    my ($ottg_key, $ottg_version);
      my $found = 0;
    foreach my $entry ( @{ $hav_gene->get_all_DBEntries } ) {
      #print "DEBUG: ENTRY: ",$entry->dbname," : ",$entry->primary_id,"\n";
      if ( $entry->dbname eq 'Vega_gene' ) {
        #print "primary id: " . $entry->primary_id . "\tdisplay_id: "
        #    . $entry->display_id . "\n";
        if ( $entry->primary_id eq $entry->display_id ) {
          $ottg_key = $entry->primary_id;
          $ottg_version = $entry->version;
          $found = 1;
        } else {
          if (!$found) {
            $ottg_key = $entry->primary_id;
            $ottg_version = $entry->version;
            warning("This entry doesn't have a display_id, it's from a patch. "
              . "Will set ottg_key = primary_id with corresponding version\n"
              . $ottg_key . " version: " . $ottg_version . "\n");
            $found = 1;
          }
        }
      }
    }

    my $read_thru;
    my $read_through_trans = 0;

    foreach my $hav_trans ( @{ $hav_gene->get_all_Transcripts } ) {
      #print "DEBUG: Havana trans id: ". $hav_trans->dbID . "\n";
      $ottg_xref{$hav_trans} = $ottg_key . "_" . $ottg_version;
      #print "ottg_hav_trans: " . $ottg_xref{$hav_trans} . "\n";

      $read_thru = 0;
      my @hav_attrib = @{ $hav_trans->get_all_Attributes() };
      foreach my $attrib (@hav_attrib) {
        if ( $attrib->name eq 'readthrough transcript' ) {
          print "Have a Hav readthrough transcript:  "
            . $hav_trans->stable_id . "\n";
          $read_thru = 1;
          $read_through_trans++;
        }
      }
    }

    my @h_clust;
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

  print "There are " . scalar(@clusters) . " havana clusters\n";
  print "There are " . scalar(@transcripts) . " ensembl transcripts\n";
  print "There are "
    . scalar(@pure_read_thru_genes)
    . " read_through havana gene clusters\n";

  # Clusters transcripts by whether or not any coding exon overlaps
  # with a coding exon in another transcript. We will use the set of
  # Havana.

  foreach my $tran (@transcripts) {

    print "\n==========\nLooking at Ensembl transcript: ", $tran->stable_id, " biotype ", $tran->biotype, "\n";
    my @matching_clusters;

    CLUSTER: foreach my $cluster (@clusters) {

      #print "DEBUG: cluster transcript: ", $cluster->[0]->stable_id, "\t"
      #    . "dbID: ", $cluster->[0]->dbID, "\n";

      if ($coding) {
        foreach my $cluster_transcript (@$cluster) {
           # print "DEBUG: \$cluster_transcript->stable_id "
           # . $cluster_transcript->stable_id
           # . "\t\$cluster_transcript->biotype "
           # . $cluster_transcript->biotype . "\n";

          if ( $cluster_transcript->translation ) {
            if ( $tran->coding_region_end >=
                    $cluster_transcript->coding_region_start
                 && $tran->coding_region_start <=
                 $cluster_transcript->coding_region_end )
            {

              #print "DEBUG: we have overlap of coding trans!\n";
              my $exons1 = $self->get_coding_exons_for_transcript($tran);
              my $cluster_exons =
                $self->get_coding_exons_for_transcript($cluster_transcript);

              foreach my $exon1 ( @{$exons1} ) {
                foreach my $cluster_exon ( @{$cluster_exons} ) {

                  if (    $exon1->overlaps($cluster_exon)
                       && $exon1->strand == $cluster_exon->strand )
                  {
                    #print "DEBUG: we have overlap of coding exons, adding "
                    #    . "the whole cluster to matching clusters\n";
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
          #print "DEBUG: \$cluster_transcript->dbID " . $cluster_transcript->dbID
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
      print "DEBUG: Have one trans in the matching clusters, "
        . $tran->stable_id
        . " adding it on to the first cluster in matching_clusters\n";

      print "DEBUG: adding "
        . $tran->stable_id . " to "
        . $matching_clusters[0]->[0]->stable_id
        . " in matching_cluster\n";

      push @{ $matching_clusters[0] }, $tran;
    } else {
      print "DEBUG: BEWARE YOU ARE HERE MERGING CLUSTERS\n";
      #print "DEBUG: Have >1 trans in the matching clusters:\n";

      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my @merged_cluster;
      my $hav_gene_counter = 0;

      # Check how many havana genes are in the cluster
    CLUST: foreach my $matching_cluster (@matching_clusters) {

        foreach my $clust_trans ( @{$matching_cluster} ) {
          #print "\tCluster -  transcript dbID: "
          #  . $clust_trans->dbID . " ("
          #  . $clust_trans->biotype . ")\n";
          if ( $clust_trans->biotype =~ /hav/ ) {
            $hav_gene_counter++;
            next CLUST;
          }
        }
      }

      # Merge all clustered genes in case there is only one or none havana genes.

      #print "DEBUG: \$hav_gene_counter: $hav_gene_counter\n";
      if ( $hav_gene_counter <= 1 ) {
        foreach my $cluster_to_be_merged (@matching_clusters) {
          push @merged_cluster, @{$cluster_to_be_merged};
        }
        push @merged_cluster, $tran;

        # Add some nice print output
        print "MERGING CLUSTER COVERED BY TRANSCRIPT IN: ",
          $tran->seq_region_name, " - ", $tran->seq_region_start, " - ",
          $tran->seq_region_end, "\n";

        #print "DEBUG: content: ", join( ' - ', @merged_cluster ), "\n";
      } else {

        # If there is more than one Havana gene we check which one has
        # a better overlap on the merger transcript and we also check
        # if the havana genes are partial and need to be merged.

        # Sort the clusters

        #print "DEBUG: Before call matching_clusters size = "
        #    .  scalar(@matching_clusters) . "\n";

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
          print "DEBUG: More than one Havana gene overlaps. Scoring the best cluster \n";
          print "DEBUG: Overlapped Hav transcript is " . ${$overlapping_cluster}[0]->stable_id
              . " biotype "                            . ${$overlapping_cluster}[0]->biotype . "\n";
          print "DEBUG: Trans id: " . $tran->dbID      . "\tbiotype: " . $tran->biotype . "\n";

          my $match_percent =
            $self->clust_overlap_percent( $overlapping_cluster, $tran );

          if ( $match_percent > $best_match ) {
            $best_match   = $match_percent;
            $best_cluster = $cluster_ref;
          }
          $cluster_ref++;
        }
        print "DEBUG: Only one havana gene gave the best overlap. Hav trans id: "
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

  #print "DEBUG: size of clusters: " . scalar(@clusters) .  "\n";

  # Adding the read_through genes that we don't want to merge in
  # to the final cluster.
  print "size of pure_read_thru_genes: " . scalar(@pure_read_thru_genes) .  "\n";
  push( @clusters, @pure_read_thru_genes );

  print "size of clusters: " . scalar(@clusters) .  "\n";

  # safety and sanity checks
  $self->check_Clusters( scalar(@transcripts), \@clusters );

  # Merge the clustered ensembl transcript with the Havana genes.

  # make and store genes
  #print "DEBUG: scalar(@clusters)." created, turning them into genes...\n";
  my @genes;
  foreach my $cluster (@clusters) {
    my $count         = 0;
    my $ottg_added    = "";
    my $gene_biotype  = "";
    my $gene          = new Bio::EnsEMBL::Gene;

    foreach my $transcript (@$cluster) {
      if ( $transcript->biotype =~ /hav/ ) {
        $gene_biotype = $transcript->biotype;
      } else {
        $gene_biotype = $transcript->biotype;
      }
      #print "DEBUG: Transcript Stable ID: " . $transcript->stable_id
      #  . "\ttranscript dbID: " . $transcript->dbID
      #  . "\t(" . $transcript->biotype . ")\n";

      my (@ottg_stable_id_array, $ottg_stable_id, $ottg_version);
      if ( defined($ottg_xref{$transcript})) {
        @ottg_stable_id_array = split('_', $ottg_xref{$transcript});
        $ottg_stable_id = $ottg_stable_id_array[0];
        $ottg_version   = $ottg_stable_id_array[1];
      }

      $gene->add_Transcript($transcript);
      if ( $ottg_stable_id
           && ( $ottg_added ne $ottg_stable_id ) )
      {
        # Need to add the OTTG as the gene stable_id otherwise it will
        # be undef.
        $gene->stable_id( $ottg_stable_id );
        $gene->version( $ottg_version );
        $self->add_ottg_xref( $gene, $ottg_stable_id, $ottg_version );
        $ottg_added = $ottg_stable_id;
      }
    }
    if ( exists($ottg_type{$ottg_added}) ) {
      my $blessed_type = $ottg_type{$ottg_added};
      #print "DEBUG: found ID: $blessed_type, keep biotype and "
      #  . "add attrib to the gene\n";

      my $attribute =
        Bio::EnsEMBL::Attribute->new( -CODE        => 'hav_gene_type',
                                      -NAME        => 'Havana gene biotype',
                                      -DESCRIPTION => 'Gene biotype assigned by Havana',
                                      -VALUE       => $blessed_type);

      $gene->add_Attributes($attribute);

      $gene->biotype($blessed_type);
    } else {
      $gene->biotype($gene_biotype);
    }

    # Keeping the ncrna_host and NoTransRefError attributes
    if ( exists($gene_attribs_to_keep{$ottg_added} ) ) {
      foreach my $attrib (@{$gene_attribs_to_keep{$ottg_added}}) {
        $gene->add_Attributes($attrib);
      }
    }

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


#########################################################################

=head2 transfer_exon_support

  Arg[1]     : Bio::EnsEMBL::Transcript object (source)
  Arg[2]     : Bio::EnsEMBL::Transcript object (target)
  Example    : $self->transfer_exon_support($source_transcript, $recipient_transcript)
  Description: Transfers supporting evidence from the source transcript's exons to
               the recipient transcript's exons.
  Returns    : Nothing, but the recipient transcript will have new exon-level supporting
               evidence added to its exons from the source transcript's exons.

=cut


sub transfer_exon_support {
  my ( $self, $trans_to_delete, $trans_to_keep ) = @_;

  my @delete_e = @{ $trans_to_delete->get_all_Exons };
  my @exons    = @{ $trans_to_keep-> get_all_Exons };

############# DEBUGGING TESTS #########
#
#  my ( $del_sf, $trans_sf ) = (0, 0);
#  foreach my $del_e ( @{ $trans_to_delete->get_all_Exons } ) {
#    print "DEBUG (transfer_exon_support): trans_to_del: " . $trans_to_delete->dbID
#      . " exon_id: " . $del_e->dbID . " :: "
#      . scalar( @{ $del_e->get_all_supporting_features } ) . "\n";
#    $del_sf += scalar( @{ $del_e->get_all_supporting_features } );
#  }
#
#  foreach my $trans_e ( @{ $trans_to_keep->get_all_Exons } ) {
#    print "DEBUG (transfer_exon_support): trans_to_keep: " . $trans_to_keep->dbID
#      . " exon_id: " . $trans_e->dbID . " :: "
#      . scalar( @{ $trans_e->get_all_supporting_features } ) . "\n";
#     $trans_sf += scalar( @{ $trans_e->get_all_supporting_features } );
#  }
#  print "DEBUG (transfer_exon_support): \$del_sf: $del_sf\n";
#  print "DEBUG (transfer_exon_support): \$trans_sf: $trans_sf\n";
#
############# DEBUGGING TESTS #########

  if ( scalar(@delete_e) == scalar(@exons) ) {
    for ( my $e = 0; $e < scalar(@delete_e); $e++ ) {
      $self->transfer_supporting_evidence_between_exons( $delete_e[$e], $exons[$e] );
    }
  } else {
    #print "DEBUG (transfer_exon_support): number of exons don't match, not doing anything.\n";
  }

} ## end sub transfer_exon_support


1;
