
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
             clone_Gene 
             Gene_info
             prune_Exons
             get_one2one_orth_for_gene_in_other_species
             get_transcript_with_longest_CDS
            );

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(print_Transcript clone_Transcript get_evidence_ids);

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw (id coord_string);
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
  my ($gene) = @_;
  my $newgene = Bio::EnsEMBL::Gene->new();
  $newgene->dbID($gene->dbID);
  $newgene->biotype($gene->biotype);
  $newgene->analysis($gene->analysis);
  $newgene->stable_id($gene->stable_id);
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    my $newtranscript = clone_Transcript($transcript);
    $newgene->add_Transcript($newtranscript);
  }
  $newgene->slice($newgene->slice);
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
  return "GENE: id ".$id." ".$coord_string." logic_name ".$gene->analysis->logic_name;
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


sub get_transcript_with_longest_CDS
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
        $nexonmax = scalar($trans->get_all_Exons)
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
  Description: Returns the one2one orthologue as a bIO::EnsEMBL::Gene object in the
               specified species or undef if there is non or more than one orthologue 
  Returntype: Bio::EnsEMBL::Gene object 
  Exceptions: none
  Example   : $one2one_orth_mus = get_one2one_orth_for_gene_in_other_species($gene, "Mus musculus"); 

=cut


sub get_one2one_orth_for_gene_in_other_species {
    my ($gene, $other_species) = @_ ;

    my $one2one_homologue_in_other_species = undef;

    foreach my $homolog_to_check  ( @{ $gene->get_all_homologous_Genes()} ) {
       my ($check_homg, $check_homology, $check_species ) = @$homolog_to_check ;

       if ($check_species eq $other_species) {
        if ($one2one_homologue_in_other_species ) {
          $one2one_homologue_in_other_species = undef;
          last;
        }
        $one2one_homologue_in_other_species = $check_homg ;
      }
    }
    return $one2one_homologue_in_other_species;
}












1;
