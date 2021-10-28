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
             clean_utrs
            );

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(print_Transcript clone_Transcript get_evidence_ids attach_Slice_to_Transcript fully_load_Transcript empty_Transcript attach_Analysis_to_Transcript attach_Analysis_to_Transcript_no_support attach_Analysis_to_Transcript_no_overwrite print_Transcript_and_Exons);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(run_translate add_ORF_to_transcript compute_6frame_translations_for_transcript calculate_sequence_content);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw (id seq_region_coord_string empty_Object);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(make_types_hash cluster_Genes);
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


sub clean_utrs {
  my ($genes, $min_size_utr_exon, $ratio_5prime_utr, $ratio_3prime_utr, $ratio_same_transcript, $ratio_max_allowed_difference, $ratio_expansion,
      $minimum_expanding_number_for_single_transcript, $ratio_transcript_fragment, $ratio_exon_expansion, $ratio_utrs, $store_rejected) = @_;

  $min_size_utr_exon //= 30;
  $ratio_5prime_utr //= .3;
  $ratio_3prime_utr //= .6;
  $ratio_same_transcript //= .02;
  $ratio_max_allowed_difference //= .05;
  $ratio_expansion //= 3;
  $minimum_expanding_number_for_single_transcript //= 2;
  $ratio_transcript_fragment //= 3;
  $ratio_exon_expansion //= 2;
  $ratio_utrs //= 2;
  $store_rejected //= 0;
  my @rejected;
  my ($clusters, $unclustered) = cluster_Genes($genes, make_types_hash($genes, undef, 'set1'));
  foreach my $cluster (@$clusters) {
    my @overlapping_genes = sort {$a->start <=> $b->start || $a->end <=> $b->end} @{$cluster->get_Genes_by_Set('set1')};
    for (my $gene_index = 0; $gene_index <= $#overlapping_genes; $gene_index++) {
      my $gene = $overlapping_genes[$gene_index];
      my $transcripts = $gene->get_all_Transcripts;
      for (my $next_gene_index = 0; $next_gene_index <= $#overlapping_genes; $next_gene_index++) {
        my $next_gene = $overlapping_genes[$next_gene_index];
        if ($gene_index != $next_gene_index) {
          if ($gene->overlaps_local($next_gene)) {
            my $change_happened = 0;
            foreach my $transcript (@$transcripts) {
              if ($transcript->overlaps_local($next_gene)) {
                my $cds_start_genomic = $transcript->coding_region_start;
                my $cds_end_genomic = $transcript->coding_region_end;
                my $cds_start_index = 0;
                my $cds_end_index = 0;
                my %overlapping_exons;
                my %exons_to_delete;
                my @exons = sort {$a->start <=> $b->start} @{$transcript->get_all_Exons};
                my $count = 0;
                foreach my $exon (@exons) {
                  if ($cds_start_genomic >= $exon->start and $cds_start_genomic <= $exon->end) {
                    $cds_start_index = $count;
                  }
                  if ($cds_end_genomic >= $exon->start and $cds_end_genomic <= $exon->end) {
                    $cds_end_index = $count;
                  }
                  if ($exon->overlaps_local($next_gene)) {
                    $overlapping_exons{$exon->start.':'.$exon->end} = $exon;
                  }
                  ++$count;
                }
                if (scalar(keys %overlapping_exons)) {
                  my $new_utr_exon_start = 0;
                  my $new_utr_exon_end = 0;
                  foreach my $next_transcript (@{$next_gene->get_all_Transcripts}) {
                    foreach my $cds_exon (@{$next_transcript->get_all_CDS}) {
                      foreach my $utr_exon (values %overlapping_exons) {
                        if ($cds_exon->overlaps_local($utr_exon)) {
                          $exons_to_delete{$utr_exon->start.':'.$utr_exon->end} = $utr_exon;
                          if (!exists $overlapping_exons{$cds_exon->start.':'.$cds_exon->end}) {
                            if ($utr_exon->start <= $cds_start_genomic and $utr_exon->end >= $cds_start_genomic
                                and $utr_exon->start <= $next_transcript->coding_region_end and $utr_exon->end >= $next_transcript->coding_region_end) {
                              $new_utr_exon_start = $cds_exon->end;
                            }
                            if ($utr_exon->start <= $cds_end_genomic and $utr_exon->end >= $cds_end_genomic
                                and $utr_exon->start <= $next_transcript->coding_region_start and $utr_exon->end >= $next_transcript->coding_region_start) {
                              $new_utr_exon_end = $cds_exon->start;
                            }
                          }
                        }
                      }
                    }
                  }
                  if (scalar(keys %exons_to_delete)) {
                    $change_happened = 1;
                    my $translation;
                    if ($new_utr_exon_end and $transcript->strand == -1) {
                      $translation = $transcript->translation;
                    }
                    elsif ($new_utr_exon_start and $transcript->strand == 1) {
                      $translation = $transcript->translation;
                    }
                    $transcript->flush_Exons;
                    $transcript->flush_IntronSupportingEvidence;
                    my $start_index = 0;
                    my $end_index = $#exons;
                    for (my $index = $cds_start_index; $index >= 0; $index--) {
                      if (exists $exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}) {
                        if ($exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}->start <= $cds_start_genomic and $exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}->end >= $cds_start_genomic) {
                          $start_index = $index;
                        }
                        last;
                      }
                      else {
                        $start_index = $index;
                      }
                    }
                    for (my $index = $cds_end_index; $index <= $#exons; $index++) {
                      if (exists $exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}) {
                        if ($exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}->start <= $cds_end_genomic and $exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}->end >= $cds_end_genomic) {
                          $end_index = $index;
                        }
                        last;
                      }
                      else {
                        $end_index = $index;
                      }
                    }
                    foreach my $exon (@exons[$start_index..$end_index]) {
                      if ($new_utr_exon_start >= $exon->start and $new_utr_exon_start <= $exon->end) {
                        if ($exon->strand == -1) {
                          $new_utr_exon_start = $cds_start_genomic-int(($cds_start_genomic-$exon->start)*$ratio_3prime_utr);
                        }
                        else {
                          $new_utr_exon_start = $cds_start_genomic-int(($cds_start_genomic-$exon->start)*$ratio_5prime_utr);
                        }
                        if ($cds_start_genomic-$new_utr_exon_start < $min_size_utr_exon) {
                          $new_utr_exon_start = $cds_start_genomic-$min_size_utr_exon;
                        }
                        if ($new_utr_exon_start < $exon->start) {
                          $new_utr_exon_start = $exon->start;
                        }
                        $translation->start($cds_start_genomic-$new_utr_exon_start+1) if ($translation);
                        $exon->start($new_utr_exon_start);
                      }
                      if ($new_utr_exon_end >= $exon->start and $new_utr_exon_end <= $exon->end) {
                        if ($exon->strand == -1) {
                          $new_utr_exon_end = $cds_end_genomic+int(($exon->end-$cds_end_genomic)*$ratio_5prime_utr);
                        }
                        else {
                          $new_utr_exon_end = $cds_end_genomic+int(($exon->end-$cds_end_genomic)*$ratio_3prime_utr);
                        }
                        if ($new_utr_exon_end-$cds_end_genomic < $min_size_utr_exon) {
                          $new_utr_exon_end = $cds_end_genomic+$min_size_utr_exon;
                        }
                        if ($new_utr_exon_end > $exon->end) {
                          $new_utr_exon_end = $exon->end;
                        }
                        $translation->start($new_utr_exon_end-$cds_end_genomic+1) if ($translation);
                        $exon->end($new_utr_exon_end);
                      }
                      $transcript->add_Exon($exon);
                    }
                    if ($translation and $translation->start_Exon == $translation->end_Exon) {
                      if ($transcript->strand == -1) {
                        $translation->start($translation->start_Exon->end-$cds_end_genomic+1);
                        $translation->end($translation->start_Exon->end-$cds_start_genomic+1);
                      }
                      else {
                        $translation->start($cds_start_genomic-$translation->start_Exon->start+1);
                        $translation->end($cds_end_genomic-$translation->start_Exon->start+1);
                      }
                    }
                  }
                }
              }
            }
            if ($change_happened) {
              my %hashes;
              my $transcripts = $gene->get_all_Transcripts;
              foreach my $transcript (@$transcripts) {
                my $id = '';
                foreach my $exon (@{$transcript->get_all_Exons}) {
                  $id .= join(':', $exon->start, $exon->end, $exon->phase, $exon->end_phase);
                }
                if ($transcript->translation) {
                  $id .= $transcript->translation->start.':'.$transcript->translation->end;
                }
                push(@{$hashes{$id}}, $transcript);
              }
              if (scalar(keys %hashes) != @$transcripts) {
                $gene->flush_Transcripts;
                foreach my $item (values %hashes) {
                  $gene->add_Transcript($item->[0]);
                }
              }
              else {
                $gene->recalculate_coordinates;
              }
              throw($gene->display_id.' has no transcript') unless (@{$gene->get_all_Transcripts});
            }
          }
        }
      }
    }
  }
  my @extra_genes;
  foreach my $gene (@$genes) {
    my @transcripts = sort {$a->end-$a->start <=> $b->end-$b->start } @{$gene->get_all_Transcripts};
    my $tcount = 0;
    foreach my $transcript (@transcripts) {
      ++$tcount if ($transcripts[0]->coding_region_start <= $transcript->coding_region_end and $transcripts[0]->coding_region_end >= $transcript->coding_region_start);
    }
    if ($tcount != scalar(@transcripts)) {
      my %data;
      foreach my $transcript (@transcripts) {
        my $stable_id = $transcript->display_id;
        $data{$stable_id}->{cds_length} = $transcript->translation->length;
        $data{$stable_id}->{cds_content} = calculate_sequence_content($transcript->translation);
        if (($transcript->coding_region_end-$transcript->coding_region_end)*$ratio_utrs < ($transcript->end-$transcript->start)) {
          my $exons = $transcript->get_all_Exons;
          my $coding_start = $transcript->coding_region_start;
          my $coding_end = $transcript->coding_region_end;
          $transcript->flush_Exons;
          $transcript->flush_IntronSupportingEvidence;
          foreach my $exon (@$exons) {
            if ($exon->start <= $coding_end and $exon->end >= $coding_start) {
              $transcript->add_Exon($exon);
            }
          }
        }
      }
      my @genes;
      my %bridging_transcripts;
      foreach my $transcript (@transcripts) {
        my $current_gene;
        foreach my $cluster_gene (reverse @genes) {
          if ($transcript->coding_region_start <= $cluster_gene->{_gb_coding_end} and $transcript->coding_region_end >= $cluster_gene->{_gb_coding_start}) {
            if ($current_gene) {
              $bridging_transcripts{$transcript->display_id} = $transcript;
              last;
            }
            else {
              $current_gene = $cluster_gene;
              if ($transcript->coding_region_start < $cluster_gene->{_gb_coding_start}) {
                $cluster_gene->{_gb_coding_start} = $transcript->coding_region_start;
              }
              if ($transcript->coding_region_end > $cluster_gene->{_gb_coding_end}) {
                $cluster_gene->{_gb_coding_end} = $transcript->coding_region_end;
              }
            }
          }
        }
        if (!exists $bridging_transcripts{$transcript->display_id}) {
          if ($current_gene) {
            $current_gene->add_Transcript($transcript);
          }
          else {
            $current_gene = Bio::EnsEMBL::Gene->new();
            $current_gene->add_Transcript($transcript);
            $current_gene->analysis($transcript->analysis);
            $current_gene->biotype('protein_coding');
            $current_gene->{_gb_coding_start} = $transcript->coding_region_start;
            $current_gene->{_gb_coding_end} = $transcript->coding_region_end;
            push(@genes, $current_gene);
          }
        }
      }
      if (scalar(keys %bridging_transcripts) == 1) {
        my ($bridging_transcript) = values %bridging_transcripts;
        my $max_allowed_difference = int($bridging_transcript->translation->length*$ratio_max_allowed_difference);
        my $bridging_stable_id = $bridging_transcript->display_id;
        my $remove_transcript = 0;
        my %bridging_exon_seen;
        my $bridging_exons = $bridging_transcript->get_all_Exons;
        foreach my $new_gene (@genes) {
          foreach my $new_transcript (@{$new_gene->get_all_Transcripts}) {
            my $new_stable_id = $new_transcript->display_id;
            if ($data{$bridging_stable_id}->{cds_length}/$data{$new_stable_id}->{cds_length} > 1-$ratio_same_transcript
                  and $data{$bridging_stable_id}->{cds_length}/$data{$new_stable_id}->{cds_length} < 1+$ratio_same_transcript) {
              my $bridge_value = 0;
              my $new_value = 0;
              my $diff = 0;
              $remove_transcript = 1;
              foreach my $key (keys %{$data{$new_stable_id}->{cds_content}}) {
                $bridge_value += $data{$bridging_stable_id}->{cds_content}->{$key} || 0;
                $new_value += $data{$new_stable_id}->{cds_content}->{$key};
                $diff = abs($bridge_value-$new_value) if (abs($bridge_value-$new_value) > $diff);
                if ($diff > $max_allowed_difference) {
                  $remove_transcript = 0;
                }
              }
            }
            else {
              my $index = -1;
              foreach my $bridging_exon (@$bridging_exons) {
                ++$index;
                next if (exists $bridging_exon_seen{$bridging_exon});
                if ($bridging_exon->start <= $new_transcript->coding_region_end and $bridging_exon->end >= $new_transcript->coding_region_start) {
                  $bridging_exon_seen{$bridging_exon} = [$bridging_exon, $new_gene, $index];
                }
              }
            }
          }
        }
        if (keys %bridging_exon_seen) {
          my $previous_gene;
          my $previous_rank;
          foreach my $item (sort {$a->[2] <=> $b->[2]} values %bridging_exon_seen) {
            if (defined $previous_gene and $previous_gene != $item->[1]) {
              if ($previous_rank+1 == $item->[2]) {
                $remove_transcript = 1;
                last;
              }
            }
            $previous_gene = $item->[1];
            $previous_rank = $item->[2];
          }
        }
        if ($remove_transcript) {
          $gene->flush_Transcripts;
          my $first_gene = shift(@genes);
          foreach my $t (@{$first_gene->get_all_Transcripts}) {
            $gene->add_Transcript($t);
          }
          foreach my $new_gene (@genes) {
            push(@extra_genes, $new_gene);
          }
          if ($store_rejected) {
            $bridging_transcript->biotype('readthrough');
            my $readthrough = Bio::EnsEMBL::Gene->new();
            $readthrough->add_Transcript($bridging_transcript);
            $readthrough->analysis($bridging_transcript->analysis);
            $readthrough->biotype($bridging_transcript->biotype);
            push(@rejected, $readthrough);
          }
        }
      }
      elsif (scalar(keys %bridging_transcripts) > 1) {
        my @no_fragments;
        $gene->flush_Transcripts;
        foreach my $small_cluster (@genes) {
          my $nof_gene = Bio::EnsEMBL::Gene->new();
          foreach my $small_transcript (@{$small_cluster->get_all_Transcripts}) {
            $nof_gene->add_Transcript($small_transcript) if (uc(substr($small_transcript->translateable_seq, 0, 3)) eq 'ATG');
          }
          if ($nof_gene->get_all_Transcripts and @{$nof_gene->get_all_Transcripts}) {
            $nof_gene->analysis($small_cluster->analysis);
            $nof_gene->biotype($small_cluster->biotype);
            push(@no_fragments, $nof_gene);
          }
          else {
          }
        }
        if (@no_fragments) {
          if (@no_fragments == 1) {
            foreach my $nof_transcript (@{$no_fragments[0]->get_all_Transcripts}, values %bridging_transcripts) {
              $gene->add_Transcript($nof_transcript);
            }
          }
          else {
            foreach my $bridging_transcript (sort {$a->end-$a->start <=> $b->end-$b->start} values %bridging_transcripts) {
              my $current_nof_gene;
              my $bridging = 0;
              foreach my $nof_gene (@no_fragments) {
                if ($bridging_transcript->coding_region_start <= $nof_gene->end and $bridging_transcript->coding_region_end >= $nof_gene->start) {
                  if ($current_nof_gene) {
                    $bridging = 1;
                  }
                  else {
                    $current_nof_gene = $nof_gene;
                  }
                }
              }
              if ($bridging) {
                if ($store_rejected) {
                  $bridging_transcript->biotype('readthrough');
                  my $readthrough = Bio::EnsEMBL::Gene->new();
                  $readthrough->add_Transcript($bridging_transcript);
                  $readthrough->analysis($bridging_transcript->analysis);
                  $readthrough->biotype($bridging_transcript->biotype);
                  push(@rejected, $readthrough);
                }
              }
              else {
                if (!$current_nof_gene) {
                  $current_nof_gene = Bio::EnsEMBL::Gene->new();
                  $current_nof_gene->analysis($bridging_transcript->analysis);
                  $current_nof_gene->biotype($bridging_transcript->biotype);
                  push(@no_fragments, $current_nof_gene);
                }
                $current_nof_gene->add_Transcript($bridging_transcript);
              }
            }
            my $first_gene = shift(@no_fragments);
            foreach my $t (@{$first_gene->get_all_Transcripts}) {
              $gene->add_Transcript($t);
            }
            foreach my $new_gene (@no_fragments) {
              push(@extra_genes, $new_gene);
            }
          }
        }
        else {
          foreach my $bridging_transcript (values %bridging_transcripts) {
            $gene->add_Transcript($bridging_transcript);
          }
        }
      }
    }
  }
  foreach my $gene (@$genes, @extra_genes) {
    my @transcripts = sort {$a->end-$a->start <=> $b->end-$b->start } @{$gene->get_all_Transcripts};
    my @genes;
    my %expanding_transcripts;
    my $current_gene_size = $transcripts[0]->end-$transcripts[0]->start;
    $gene->flush_Transcripts;
    foreach my $transcript (@transcripts) {
      if ($current_gene_size*$ratio_expansion > $transcript->end-$transcript->start) {
        $current_gene_size = $transcript->end-$transcript->start;
        $gene->add_Transcript($transcript);
      }
      else {
        $expanding_transcripts{$transcript->display_id} = $transcript
      }
    }
    if (@{$gene->get_all_Transcripts} == 1 and scalar(keys %expanding_transcripts) > $minimum_expanding_number_for_single_transcript) {
      foreach my $expanding_transcript (values %expanding_transcripts) {
        $gene->add_Transcript($expanding_transcript);
      }
    }
    elsif (scalar(keys %expanding_transcripts) >= @{$gene->get_all_Transcripts}) {
      my $expansion_added = 0;
      my $exon_count = 0;
      my %possible_fragments;
      foreach my $t (@{$gene->get_all_Transcripts}) {
        $exon_count = scalar(@{$t->get_all_Exons}) if ($exon_count < scalar(@{$t->get_all_Exons}));
        $possible_fragments{$t} = 1;
      }
      foreach my $expanding_transcript (values %expanding_transcripts) {
        if ($gene->length*$ratio_expansion >= ($expanding_transcript->coding_region_end-$expanding_transcript->coding_region_start+1)) {
          my $expanding_exons = $expanding_transcript->get_all_Exons;
          $expanding_transcript->flush_Exons;
          $expanding_transcript->flush_IntronSupportingEvidence;
          my %utr_5p = map {$_->start.':'.$_->end => $_ } @{$expanding_transcript->get_all_five_prime_UTRs};
          foreach my $expanding_exon (@$expanding_exons) {
            $expanding_transcript->add_Exon($expanding_exon) unless (exists $utr_5p{$expanding_exon->start.':'.$expanding_exon->end});
          }
          $expanding_exons = $expanding_transcript->get_all_Exons;
          $expanding_transcript->flush_Exons;
          $expanding_transcript->flush_IntronSupportingEvidence;
          my %utr_3p = map {$_->start.':'.$_->end => $_ } @{$expanding_transcript->get_all_three_prime_UTRs};
          foreach my $expanding_exon (@$expanding_exons) {
            $expanding_transcript->add_Exon($expanding_exon) unless (exists $utr_3p{$expanding_exon->start.':'.$expanding_exon->end});
          }
          $gene->add_Transcript($expanding_transcript);
        }
        elsif ($exon_count < scalar(@{$expanding_transcript->get_all_CDS})) {
          $gene->add_Transcript($expanding_transcript);
          ++$expansion_added;
        }
        else {
          if ($store_rejected) {
            my $expanding_gene = Bio::EnsEMBL::Gene->new();
            $expanding_gene->add_Transcript($expanding_transcript);
            $expanding_gene->analysis($expanding_transcript->analysis);
            $expanding_gene->biotype('expanding');
            push(@rejected, $expanding_gene);
          }
        }
      }
      if ($expansion_added >= scalar(keys %expanding_transcripts) and scalar(keys %possible_fragments)*$ratio_transcript_fragment < $expansion_added) {
        my $all_transcripts = $gene->get_all_Transcripts;
        $gene->flush_Transcripts;
        foreach my $all_transcript (@$all_transcripts) {
          if (exists $possible_fragments{$all_transcript}) {
            if ($store_rejected) {
              my $fragment_gene = Bio::EnsEMBL::Gene->new();
              $fragment_gene->add_Transcript($all_transcript);
              $fragment_gene->analysis($all_transcript->analysis);
              $fragment_gene->biotype('fragment');
              push(@rejected, $fragment_gene);
            }
          }
          else {
            $gene->add_Transcript($all_transcript);
          }
        }
      }
    }
    else {
      my $expansion_added = 0;
      my $exon_count = 0;
      my %possible_fragments;
      foreach my $t (@{$gene->get_all_Transcripts}) {
        $exon_count = scalar(@{$t->get_all_Exons}) if ($exon_count < scalar(@{$t->get_all_Exons}));
        $possible_fragments{$t} = 1;
      }
      $exon_count *= $ratio_exon_expansion;
      foreach my $expanding_transcript (values %expanding_transcripts) {
        my $process_transcript = 1;
        my $is_new_gene;
        my $genomic_start = $expanding_transcript->coding_region_start;
        my $genomic_end = $expanding_transcript->coding_region_end;
        my $expanding_exons = $expanding_transcript->get_all_Exons;
        $expanding_transcript->flush_Exons;
        $expanding_transcript->flush_IntronSupportingEvidence;
        foreach my $expanding_exon (@$expanding_exons) {
          if ($expanding_exon->start <= $genomic_end and $expanding_exon->end >= $genomic_start
             or $expanding_exon->overlaps_local($gene)) {
            $expanding_transcript->add_Exon($expanding_exon);
          }
        }
        if ($expanding_transcript->overlaps_local($gene)) {
          if ($gene->length*$ratio_expansion >= ($expanding_transcript->end-$expanding_transcript->start+1)) {
            $gene->add_Transcript($expanding_transcript);
            $process_transcript = 0;
          }
        }
        else {
          $is_new_gene = $expanding_transcript->get_all_Exons;
        }
        if ($process_transcript) {
          if ($exon_count < scalar(@{$expanding_transcript->get_all_CDS})) {
            $gene->add_Transcript($expanding_transcript);
            ++$expansion_added;
            $process_transcript = 0;
          }
        }
        if ($process_transcript) {
          if ($store_rejected) {
            my $expanding_gene = Bio::EnsEMBL::Gene->new();
            $expanding_gene->add_Transcript($expanding_transcript);
            $expanding_gene->analysis($expanding_transcript->analysis);
            $expanding_gene->biotype('expanding');
            push(@rejected, $expanding_gene);
          }
        }
      }
      if ($expansion_added >= scalar(keys %expanding_transcripts) and scalar(keys %possible_fragments)*$ratio_transcript_fragment < $expansion_added) {
        my $all_transcripts = $gene->get_all_Transcripts;
        $gene->flush_Transcripts;
        foreach my $all_transcript (@$all_transcripts) {
          if (exists $possible_fragments{$all_transcript}) {
            if ($store_rejected) {
              my $fragment_gene = Bio::EnsEMBL::Gene->new();
              $fragment_gene->add_Transcript($all_transcript);
              $fragment_gene->analysis($all_transcript->analysis);
              $fragment_gene->biotype('fragment');
              push(@rejected, $fragment_gene);
            }
          }
          else {
            $gene->add_Transcript($all_transcript);
          }
        }
      }
    }
  }
}


1;
