
=head1 NAME

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils - utilities for gene objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(clone_Gene);

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils 
  
  to get all methods

=head1 DESCRIPTION

All methods in this class should take a Bio::EnsEMBL::Gene
object as their first argument.

The methods provided should carry out some standard 
functionality for said objects such as printing info, and 
cloning

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;

use strict;
use warnings;
use Exporter;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             print_Gene 
             print_Gene_Transcript_and_Exons
             clone_Gene 
             Gene_info
             prune_Exons
             get_one2one_orth_for_gene_in_other_species
             get_one2one_homology_for_gene_in_other_species 
             get_transcript_with_longest_CDS
             attach_Slice_to_Gene
             attach_Analysis_to_Gene
             attach_Analysis_to_Gene_no_support
             fully_load_Gene
             empty_Gene
             compute_6frame_translations
             convert_to_single_transcript_gene
            );

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(print_Transcript clone_Transcript get_evidence_ids attach_Slice_to_Transcript fully_load_Transcript empty_Transcript attach_Analysis_to_Transcript attach_Analysis_to_Transcript_no_support print_Transcript_and_Exons);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(run_translate add_ORF_to_transcript ); 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw (id coord_string empty_Object);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils;
use Bio::EnsEMBL::Gene;



=head2 print_Gene

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : prints out information about the gene object 
  passed in, it will also iterate down the gene structure 
  printing information about transcripts, exons etc
  Returntype: n/a
  Exceptions: n/a
  Example   : print_Gene($gene);

=cut




sub print_Gene{
  my $gene = shift;
  print Gene_info($gene)."\n";
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    my $indent = "\t";
    print_Transcript($transcript, $indent);
  }
}


sub print_Gene_Transcript_and_Exons{
  my $gene = shift;
  print Gene_info($gene)."\n";
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    my $indent = "\t";
    print_Transcript_and_Exons($transcript, $indent);
  }
}


=head2 clone_Gene

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : to produce a copy of the given gene object
  which can be altered without changing the original object or 
  its children to the original object
  Returntype: Bio::EnsEMBL::Gene
  Exceptions: none
  Example   : my $newgene = clone_Gene($gene);

=cut



sub clone_Gene{
  my ($gene, $clone_xrefs) = @_; 

  $clone_xrefs = 1 if(!defined($clone_xrefs));
  my $newgene = Bio::EnsEMBL::Gene->new();
  $newgene->dbID($gene->dbID);
  $newgene->biotype($gene->biotype);
  $newgene->analysis($gene->analysis);
  $newgene->stable_id($gene->stable_id);
  $newgene->version($gene->version);
  if ($clone_xrefs){
    foreach my $DBEntry (@{$gene->get_all_DBEntries}){
      $newgene->add_DBEntry($DBEntry);
    } 
    $newgene->display_xref($gene->display_xref) ; 
  }  
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    my $newtranscript = clone_Transcript($transcript, $clone_xrefs);
    $newgene->add_Transcript($newtranscript);
  }
  $newgene->slice($gene->slice);
  return $newgene;
}



=head2 Gene_info

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : returns a string of information about the gene
  Returntype: n/a
  Exceptions: none
  Example   : print Gene_info($gene)."\n";

=cut




sub Gene_info{
  my ($gene) = @_;
  my $coord_string = coord_string($gene);
  my $id = id($gene);
  return "GENE: id ".$id." ".$coord_string." biotype ".$gene->biotype;
}




=head2 prune_Exons

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : remove duplicate exons between Transcripts but
  ensure translation and other exons are maintained
  Returntype: 
  Exceptions: 
  Example   : 

=cut



sub prune_Exons{
  my ($gene) = @_;
  my @unique_exons;
  # keep track of all unique exons found so far to avoid making 
  # duplicates need to be very careful about 
  # translation->start_Exon and translation->end_Exon 
  #  
  my $cloned_gene = clone_Gene($gene); 
  foreach my $tran (@{$cloned_gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty
    UNI:foreach my $uni (@unique_exons) {
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
        push(@unique_exons, $exon);
      }
    }          
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }

  return $cloned_gene;
}


=head2 get_transcript_with_longest_CDS

  Arg [0]   : Bio::EnsEMBL::Gene
  Function  : returns Array containing the longest transcript-obj, the length of all translateable exons and 
              the number of exons for this transcript 
  Returntype: Array longest transcript, max. length, max nr of exons for this transcript 
  Exceptions: 
  Example   : 

=cut


sub get_transcript_with_longest_CDS {
  my $gene = shift;

  my $maxlen = 0;
  my $nexonmax = 0;
  my $longest; 

  foreach my $trans (@{$gene->get_all_Transcripts}) {
    my $len = 0;
    if ($trans->translation) {
      my @trans_exons = @{$trans->get_all_translateable_Exons()};
      foreach my $exon (@trans_exons) {
        $len+= $exon->length;
      }
      if ($len > $maxlen) {
        $maxlen = $len;
        $longest = $trans;
        $nexonmax = scalar(@{$trans->get_all_Exons});
      }
    }
  }
  if (!defined($longest)) {
    print "Didn't find a longest transcript for gene ". $gene->stable_id . " type " . $gene->type . "\n";
  }
  return ($longest,$maxlen,$nexonmax);
}


=head2 get_one2one_orth_for_gene_in_other_species 

  Arg [0]   : Bio::EnsEMBL::Gene
  Arg [1]   : String species_name  ("Homo sapiens") 
  Description: Returns the one2one orthologue as a Bio::EnsEMBL::Gene object in the
               specified species or undef if there is non or more than one orthologue 
  Returntype: Bio::EnsEMBL::Gene object 
  Exceptions: none
  Example   : $one2one_orth_mus = get_one2one_orth_for_gene_in_other_species($gene, "Mus musculus"); 

=cut


sub get_one2one_orth_for_gene_in_other_species {
    my ($gene, $other_species) = @_ ;


    # this is a list of homolgy-desription terms which are used in ensembl_compara databases 
    # versions v18 - v 47. They are stored in homology.description table. 


    my @orthology_description_dictionary = qw (  
                                                PIP
                                                SEED
                                                BRH
                                                DWGA
                                                MBRH
                                                RHS
                                                UBRH
                                                YoungParalogues
                                                ortholog_many2many
                                                ortholog_one2many
                                                ortholog_one2one
                                                apparent_ortholog_one2one
                                                between_species_paralog
                                                within_species_paralog
                                               ) ; 

   my %orthology_dict; 
   @orthology_dict{@orthology_description_dictionary} = 1;

   my $one2one_homologue_in_other_species = undef; 

   my @homologies_found ; 

    foreach my $homolog_to_check  ( @{ $gene->get_all_homologous_Genes()} ) {
       my ($check_homg, $check_homology, $check_species ) = @$homolog_to_check ;
 
        unless ( exists $orthology_dict{$check_homology->description}){ 
          throw("The homology-description : ".$check_homology->description. " is an unknow, probably new description. This".
                " is probably due to changes in compara. Please add the new description to the \@orthology_description_dictionary ".
                " in the GeneUtils.pm module and make sure that the  retrieval of one2one orthologues is not affected\n") ; 
        }

       if ($check_species eq $other_species) {  
         push @homologies_found , [$check_homology,$check_homg]; 
         #print "Homology found : " .  $check_homology->description. "\t" . $check_homg->stable_id . "\n" ;    

         # now this is a tricky one ... depending on the string we need to make a decision if it's a one2one or not. 
         # ortholgos can be retrieved by gene-trees as well as the homology pipeline. 
         # compara changes this string ( stored in the compara homology-table column description

         if ( $check_homology->description=~m/ortholog_one2one/){ 
  
            # a ortholog_one2one-relation has already been found - this is wrong, or compara changed ...          
            
            if ($one2one_homologue_in_other_species ) {   
              my $err_str = "We have more than one ortholog_one2one-relationship found for : " . $gene->stable_id ."\t" . $check_species. "\n" ;   

              for my $hg ( @homologies_found ) {  
                my ( $homology, $homolog ) = @$hg  ; 
                $err_str .= $homology->description ."\t" . $homolog->stable_id . "\n" ; 
              } 
              throw("$err_str"); 
            } 
            $one2one_homologue_in_other_species = $check_homg ;
         } elsif ( $check_homology->description=~m/(PIP|SEED|BRH|DWGA|MBRH|RHS|UBRH)/){ 
           if ( $one2one_homologue_in_other_species ) {   
             print "more than one homology found, so i assume that it's not a one2one homolog\n" ; 
             return undef ; 
           } 
           $one2one_homologue_in_other_species = $check_homg ;        
         }
       } 
    }  
#    if ( scalar (@homologies_found ) > 0 ) { 
#      print "\nRelationships identified for : " . $gene->stable_id . " :\n" ; 
#      for my $hg ( @homologies_found ) {   
#         my ( $homology, $homolog ) = @$hg  ; 
#         print "\t" .  $homology->description ."\t" . $homolog->stable_id . "\n" ; 
#      } 
#      print "\n" ;   
#    }  
    return $one2one_homologue_in_other_species;
}



sub get_one2one_homology_for_gene_in_other_species {
    my ($gene, $other_species) = @_ ;

   print "using method get_one2one_homology_for_gene_in_other_species analysis_2008_01_21   \n" ; 

    my $one2one_homologue_gene_in_other_species = undef;
    my $homology = undef ;
    my $more_than_one = undef ;
    my @ortholog_one2one;

    # problem with the old routine is that it only works with old compara 
    # versions ( like with versions where there are homology-classes like 
    # UBRH, MBRH etc. ) - if we use the new, tree-based compara pipeline we 
    # get other homology descriptions. 


    foreach my $homolog_to_check  ( @{ $gene->get_all_homologous_Genes()} ) {
       my ($check_homologue_gene, $check_homology, $check_species ) = @$homolog_to_check ;

       if ($check_species eq $other_species) {

        print "homology_key\t$other_species\t".$check_homology->description."\t" .
        $gene->stable_id  . "\t" .  $check_homologue_gene->stable_id  ."\n" ;



        if ( $check_homology->description =~m/ortholog_one2one/ ||
             $check_homology->description =~m/apparent_ortholog_one2one/)  {
          push @ortholog_one2one, $check_homology ;
        }

        if ( $one2one_homologue_gene_in_other_species ) {
          $one2one_homologue_gene_in_other_species = undef;
          $more_than_one = 1 ;
        }
        $one2one_homologue_gene_in_other_species = $check_homologue_gene ;
        $homology = $check_homology ;
      }
    }

    if ( $more_than_one ) {
      if ( scalar(@ortholog_one2one) eq 1 && $ortholog_one2one[0]){
        return $ortholog_one2one[0] ;
      } elsif ( scalar(@ortholog_one2one) > 1) {
        throw("weird more than one ortholog_one2one found\n");
      } else {
        #print "no orthologue_one2one found\n" ;         
        $homology = undef ;
      }
    }
    return $homology;
}





=head2 attach_Slice_to_Gene

  Arg [1]   : Bio::EnsEMBL::Gene
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : attach given slice to gene and its object hierachy
  Returntype: n/a
  Exceptions: 
  Example   : 

=cut


sub attach_Slice_to_Gene{
  my ($gene, $slice) = @_;
  $gene->slice($slice);
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    attach_Slice_to_Transcript($transcript, $slice);
  }
}


sub attach_Analysis_to_Gene{
  my ($gene, $analysis) = @_;
  $gene->analysis($analysis);
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    attach_Analysis_to_Transcript($transcript, $analysis);
  }
}

sub attach_Analysis_to_Gene_no_support{
  my ($gene, $analysis) = @_;
  $gene->analysis($analysis);
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    attach_Analysis_to_Transcript_no_support($transcript, $analysis);
  }
}
=head2 fully_load_Gene

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : descend object hierachy to ensure is fully loaded from the database
  Returntype: Bio::EnsEMBL::Gene
  Exceptions: 
  Example   : 

=cut


sub fully_load_Gene{
  my ($gene, $keep_xrefs) = @_;
  $keep_xrefs = 1 if(!defined($keep_xrefs));
  foreach my $t(@{$gene->get_all_Transcripts}){
    fully_load_Transcript($t, $keep_xrefs);
  }
  $gene->analysis;
  $gene->get_all_DBEntries if($keep_xrefs);
  $gene->get_all_Attributes;
  $gene->stable_id;
  return $gene;
}


=head2 empty_Gene

  Arg [1]   : Bio::EnsEMBL::Gene
  Arg [2]   : Boolean, whether to remove the stable id
  Arg [3]   : Boolean, whether to remove the xrefs
  Function  : detached gene and its object hierachy from database
  Returntype: Bio::EnsEMBL::Gene
  Exceptions: 
  Example   : 

=cut



sub empty_Gene{
  my ($gene, $remove_stable_id, $remove_xrefs) = @_;
  fully_load_Gene($gene);
  foreach my $t(@{$gene->get_all_Transcripts}){
    empty_Transcript($t, $remove_stable_id, $remove_xrefs);
  }
  $gene->display_xref(undef) if($remove_xrefs);
  empty_Object($gene, $remove_stable_id);
  return $gene;
}


=head2 compute_6frame_translations

  Arg [1]   : Bio::EnsEMBL::Gene 
  Function  : computes all possible 6-frame-translations for all transcripts of a gene 
              and returns a new Bio::EnsEMBL::Gene object with one Transcript added for each 
              translation found; used to check if ncRNA's can be translated + contain protein_domains ...  
  Returntype: Bio::EnsEMBL::Gene
  Exceptions: Warns in unable to create translation
  Example   : $new_gene = compute_6frame_translations($old_gene) ; 

=cut



sub compute_6frame_translations{
  my ($gene ) = @_;

  my @tr = @{ $gene->get_all_Transcripts};

  my @new_transcripts ;

  TRANSCRIPTS: for my $transcript ( @tr ) {

    my @met_predictions = @{run_translate ($transcript, 1)};
    my @nomet_predictions = @{run_translate ($transcript)};

    ORF: for my $orf ( @met_predictions, @nomet_predictions ) {
      # create new transcript for every ORF 
      my $nt = new Bio::EnsEMBL::Transcript( -EXONS => $transcript->get_all_Exons) ; 
      $nt = add_ORF_to_transcript($orf, $nt) ;  
      push @new_transcripts, $nt ;  
    } 
 }
    
 my $new_gene = Bio::EnsEMBL::Gene->new();
 $new_gene->biotype("");
 
 for my $nt ( @new_transcripts ) { 
   $new_gene->add_Transcript($nt) ;
 }  
  return $new_gene ;
}


1;


sub convert_to_single_transcript_gene {
  my ($multitrans_gene) = @_;

  my @singletrans_gene ;

  foreach my $trans ( @{ $multitrans_gene->get_all_Transcripts } ) {
  my $newgene = Bio::EnsEMBL::Gene->new();

  $newgene->add_Transcript($trans) ;
  $newgene->biotype($multitrans_gene->biotype) ;

  push @singletrans_gene, $newgene;
  }

return \@singletrans_gene ;
}
