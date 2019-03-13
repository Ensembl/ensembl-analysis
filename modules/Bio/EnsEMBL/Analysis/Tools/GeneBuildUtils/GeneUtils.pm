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

please send any questions to http://lists.ensembl.org/mailman/listinfo/dev

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
             attach_Slice_to_Gene
             attach_Analysis_to_Gene
             attach_Analysis_to_Gene_no_support
             attach_Analysis_to_Gene_no_ovewrite
             clone_Gene 
             compute_6frame_translations
             convert_to_single_transcript_gene
             empty_Gene 
             filter_Genes_by_Exon_count
             fully_load_Gene
             get_one2one_orth_for_gene_in_other_species
             get_one2one_homology_for_gene_in_other_species  
             get_readthroughs_count 
             get_single_Exon_Genes 
             get_multi_Exon_Genes
             get_transcript_with_longest_CDS
             Gene_info
             print_Gene 
             print_Gene_Transcript_and_Exons
             prune_Exons
             remove_Transcript_from_Gene
             validate_store
            );

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(print_Transcript clone_Transcript get_evidence_ids attach_Slice_to_Transcript fully_load_Transcript empty_Transcript attach_Analysis_to_Transcript attach_Analysis_to_Transcript_no_support attach_Analysis_to_Transcript_no_overwrite print_Transcript_and_Exons);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(run_translate add_ORF_to_transcript compute_6frame_translations_for_transcript); 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw (id seq_region_coord_string empty_Object);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils;
use Bio::EnsEMBL::Gene;



=head2 filter_Genes_by_Exon_count 

  Arg [1]   : Arrayref of Bio::EnsEMBL::Gene objects with 1 or more transcripts 
  Function  : filters Genes into 2 subsets - genes which only 'have' single-exon-transcripts and 
              Genes with transcripts with more than one exon. 
  Returntype: n/a
  Exceptions: n/a
  Example   : my ( $single_exon_genes, $multi_exon_genes ) = filter_Genes_by_Exon_count(\@genes);

=cut


sub filter_Genes_by_Exon_count {
  my ($genes ) = @_;

  my (@single_exon_genes,@multi_exon_genes) ;

  for my $g ( @$genes ) {
    my $max_nr_exons = 0 ;
    for my $t ( @{$g->get_all_Transcripts} ) {
      my $exons = scalar( @{$t->get_all_Exons} );
      if ( $max_nr_exons < $exons ) {
        $max_nr_exons = $exons;
      }
    }
    if ( $max_nr_exons == 1 ) {
      push @single_exon_genes, $g ;
    } elsif ( $max_nr_exons > 1 ) {
      push @multi_exon_genes, $g;
    } else {
      throw (" Gene ".$g->dbID. " does not have any exons !");
    }
  }
  return ( \@single_exon_genes , \@multi_exon_genes ) ;
}


sub get_single_Exon_Genes {
  my ($genes ) = @_;
  my ($single_exon_genes , $multi) = filter_Genes_by_Exon_count($genes);  
  return $single_exon_genes ; 
} 

sub get_multi_Exon_Genes {
  my ($genes ) = @_;
  my ($single, $multi_exon_genes ) = filter_Genes_by_Exon_count($genes); 
  return $multi_exon_genes ;  
}


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
  $newgene->source($gene->source);
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
  my $coord_string = seq_region_coord_string($gene);
  my $id = $gene->display_id;
  my $logic_name = $gene->analysis ? $gene->analysis->logic_name : 'NO_ANALYSIS';
  return "GENE: id ".$id." ".$coord_string." biotype ".$gene->biotype.' '.$logic_name;
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


=head2 attach_Analysis_to_Gene_no_ovewrite

 Arg [1]    : Bio::EnsEMBL::Gene
 Arg [2]    : Bio::EnsEMBL::Analysis
 Description: Attach a Arg[2] to Arg[1] and its sub objects like transcripts
              unless the analysis is already set
 Returntype : None
 Exceptions : Throws if Arg[2] is not a Bio::EnsEMBL::Analysis

=cut

sub attach_Analysis_to_Gene_no_ovewrite {
  my ($gene, $analysis) = @_;

  throw('You need a Bio::EnsEMBL::Analysis object not a "'.ref($analysis).'"')
    unless ($analysis and ref($analysis) eq 'Bio::EnsEMBL::Analysis');
  $gene->analysis($analysis) unless ($gene->analysis);
  foreach my $transcript (@{$gene->get_all_Transcripts}) {
    attach_Analysis_to_Transcript_no_overwrite($transcript, $analysis);
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

  # and more!
  $gene->canonical_transcript();
  #$gene->canonical_annotation();
  #$gene->get_all_alt_alleles();
  # etc

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
  if ($remove_xrefs) {
    $gene->display_xref(undef);
    # It is naughty to go into the gene object
    # but we have no API method to do this:
    $gene->{'dbentries'} = [];
  }
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
    push @new_transcripts, @{compute_6frame_translations_for_transcript($transcript)}; 
  } 
    
  my $new_gene = Bio::EnsEMBL::Gene->new();
  $new_gene->biotype($gene->biotype);
 
 for my $nt ( @new_transcripts ) { 
   $new_gene->add_Transcript($nt) ;
 }  
  return $new_gene ;
}



=head2 convert_to_single_transcript_gene 

  Arg [1]   : Bio::EnsEMBL::Gene 
  Function  : converts all Transcripts of a Gene into new Bio::EnsEMBL::Gene objects and 
              returns them as array reference 
  Returntype: arraryref of Bio::EnsEMBL::Gene
  Exceptions: 
  Example   : @single_transcript_genes = @{ convert_to_single_transcript_gene($gene)} 

=cut

sub convert_to_single_transcript_gene { 
  my ($gene) = @_; 

  my @single_transcript_genes;

  my @tr = @{$gene->get_all_Transcripts};
  if ( scalar(@tr) == 1 ) {
     push @single_transcript_genes, $gene ;
  } else {
     foreach my $transcript(@tr) {
       my $ng = Bio::EnsEMBL::Gene->new();
       $ng->add_Transcript($transcript);
       $ng->biotype($gene->biotype) ; 
       push @single_transcript_genes, $ng ;     
     }
  }
  return \@single_transcript_genes;
}

=head2 get_readthroughs_count

  Arg [1]   : Bio::EnsEMBL::Gene 
  Function  : returns the number of transcripts having a readthrough transcript attribute
  Returntype: int
  Exceptions: none
  Example   : $num_readthrough_transcripts = get_readthroughs_count($gene)

=cut

sub get_readthroughs_count { 
  my ($gene) = @_; 

  my $readthrough_transcripts = 0;
  foreach my $t (@{$gene->get_all_Transcripts()}) {
    if (scalar(@{$t->get_all_Attributes('readthrough_tra')}) > 0) {
      $readthrough_transcripts++;
    }
  }
  return $readthrough_transcripts;
}


=head2 remove_Transcript_from_Gene

 Arg [1]    : Bio::EnsEMBL::Gene, gene to remove a transcript from
 Arg [2]    : Bio::EnsEMBL::Transcript, transcript to remove from the gene
 Arg [3]    : Hashref of String (optional), Hash of biotypes which are not to be removed
 Description: Remove a transcript from a gene.
              You should use this method if the objects are not in a database or if you
              do not want the gene coordinates to be updated in the source database.
              Otherwise you might want to use: $gene->remove_Transcript($transcript)
 Returntype : Boolean, 1 if the transcript has been removed, 0 if the transcript is blessed
 Exceptions : None

=cut

sub remove_Transcript_from_Gene {
  my ($gene, $transcript_to_delete, $blessed_biotypes) = @_;

  if ($blessed_biotypes and exists $blessed_biotypes->{$transcript_to_delete->biotype}) {
    return 0;
  }
  my $transcripts = $gene->get_all_Transcripts;
  $gene->flush_Transcripts;
  foreach my $transcript (@$transcripts) {
    if ($transcript != $transcript_to_delete) {
      $gene->add_Transcript($transcript);
    }
  }
  return 1;
}


=head2 validate_store

 Arg [1]    : Bio::EnsEMBL::Gene, the gene you tried to store
 Arg [2]    : Bio::EnsEMBL::Gene, the stored copy of the gene read from the output db via the dbID

 Description: Compare a gene in memory to the stored copy of the gene in the output db. We have seen issues with
              the API where sometimes the gene is incompletely stored, but this is not caught as an issue by the
              core API.
 Returntype : Boolean, 1 if the gene has been stored correctly (down to the exon level), 0 if it hasn't
 Exceptions : None

=cut

sub validate_store {
  my ($gene1, $gene2) = @_;

  unless($gene1 && $gene2) {
    warning("One or both of the input genes were not passed in");
    return(0);
  }

  my $transcripts1 = $gene1->get_all_Transcripts();
  my $transcripts2 = $gene2->get_all_Transcripts();
  my $transcripts_count1 = scalar(@$transcripts1);
  my $transcripts_count2 = scalar(@$transcripts2);

  unless($transcripts_count1 && $transcripts_count2 && ($transcripts_count1 == $transcripts_count2)) {
    warning("Transcript counts do not match between the gene in memory and the stored gene");
    return(0);
  }

  for(my $i=0; $i<$transcripts_count1; $i++) {
    my $transcript1 = $$transcripts1[$i];
    my $transcript2 = $$transcripts2[$i];
    my $exons1 = $transcript1->get_all_Exons();
    my $exons2 = $transcript2->get_all_Exons();

    my $exons_count1 = scalar(@$exons1);
    my $exons_count2 = scalar(@$exons2);

    unless($exons_count1 && $exons_count2 && ($exons_count1 == $exons_count2)) {
      warning("Exon counts do not match between the gene in memory and the stored gene.\nStored transcript id: ".$transcript2->dbID."\n".
              "Exon count in memory: ".$exons_count1."\nExon count in db: ".$exons_count2);
      return(0);
    }

    if($transcript1->translation) {
      unless($transcript2->translation) {
        warning("The transcript in memory has a translation but the stored transcript does not");
        return(0);
      }

      my $translation1 = $transcript1->translate->seq;
      my $translation2 = $transcript2->translate->seq;
      unless($translation1 eq $translation2) {
        warning("The translation from the transcript in memory does not match the stored transcript.\nTranslation in memory:\n".$translation1.
                "\nTranslation in db:\n".$translation2);
        return(0);
      }
    }
  }
  return 1;
}

1;
