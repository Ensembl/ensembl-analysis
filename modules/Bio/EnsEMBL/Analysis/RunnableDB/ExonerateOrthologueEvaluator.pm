=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateOrthologueEvaluator - 

=head1 SYNOPSIS

my $exonerate2genes = Bio::EnsEMBL::Analysis::RunnableDB::ExonerateOrthologueEvaluator->new(
                              -db         => $refdb,
			      -analysis   => $analysis_obj,
			      -database   => $EST_GENOMIC,
			      -query_seqs => \@sequences,
			      -query_type => 'dna',
			     );

$exonerate2genes->fetch_input();
$exonerate2genes->run();
$exonerate2genes->output();
$exonerate2genes->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::Exonerate2Genes and is used 
in the Orthologue-evaluation process. Before writing Gene structures in the ORTHOLOGUE_DB
it checks if these structures already exists in either core or the ORTHOLOGUE_DB
itself. 


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/ExonerateOrthologueEvaluator.pm,v $
# $Version: $
package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateOrthologueEvaluator;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes ; 
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene; 
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator; 

use vars qw(@ISA);
@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes);

sub run{
  my ($self) = @_; 
  print "starting run method \n" ; 
  $self->SUPER::run(); 
  print "run method\n" ;    
   
  # load connection to core where the 'old' genes are stored in 
 
  Bio::EnsEMBL::Registry->load_all($$MAIN_CONFIG{LOCATION_OF_COMPARA_REGISTRY_FILE},1,1);
  my $storedb = Bio::EnsEMBL::Registry->get_DBAdaptor($$MAIN_CONFIG{QUERY_SPECIES},'core') ;   
  throw ("Could not get dba\n") unless $storedb ;   

  # compare against CORE and OUTPUT  

  my $outdb = $self->get_output_db;
  my $outdb_ga = $outdb->get_GeneAdaptor;  
  my $store_ga = $storedb->get_GeneAdaptor;   
  my $sa = $outdb->get_SliceAdaptor;  
   
#  my $test_slice = 'chromosome:Btau_3.1:22:30177919:30506937:1'  ; 
#  my $slice = $dba->get_SliceAdaptor->fetch_by_name($test_slice) ; 
#  my $ga = $dba->get_GeneAdaptor();  
#  my @new_genes = @{$ga->fetch_all_by_Slice($slice)};  


  my @dbas = ( $outdb_ga, $store_ga ) ; 

  my @new_genes = @{ $self->output } ;  
  my %bt;
  $bt{$outdb_ga}="_out";
  $bt{$store_ga}="_core";

  for my $ga (@dbas ) { 
    for my $gene ( @new_genes ) { 
      my $slice = $sa->fetch_by_name($gene->slice->name) ;  
      my @stored_genes = @{$ga->fetch_all_by_Slice($slice)};  
      for my $sg( @stored_genes) {     
        if ( compare_Genes($gene, $sg ) ) {  
          $gene->biotype("redundant".$bt{$ga});  
          # print "gene is the same\n" ; 
        } 
      }
    } 
  }
}


sub compare_Genes {
  my ($gene1,$gene2) = @_; 

  # quit if genes do not have genomic overlap 
  
  # start-------gene1------end   start--------gene2----------end
  #  

  if ($gene1->end < $gene2->start || $gene1->start > $gene2->end) {
    # Failed extents check
    return 0;
  }
    # overlap check based on all (noncoding + coding) Exons 

    my %exon; 
    my $all_exon_overlap;
    foreach my $exon1 (@{$gene1->get_all_Exons}){ 
      print $exon1->hashkey."\n" ; 
      foreach my $exon2 (@{$gene2->get_all_Exons}){
        if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
          #print "Passed exon overlap check (noncod. + cod. exons checked)  - returning 1\n";
          $exon{$exon1->hashkey} = 1;  
          #print 'hashkeys are the same' if ( $exon1->hashkey eq $exon2->hashkey ) ; 
          #return 1;
        }
      }
    } 
    if ( scalar (keys %exon)  == scalar (@{$gene1->get_all_Exons}) ){ 
      return 1 ; 
    } 

  return 0;
}




sub write_output{
  my ($self,@output) = @_;

  my $outdb = $self->get_output_db;
  my $gene_adaptor = $outdb->get_GeneAdaptor;  

  unless (@output){
    @output = @{$self->output};
  }
  
  my $fails = 0;
  my $total = 0;
  foreach my $gene (@output){
    eval {
      $gene_adaptor->store($gene);
    };    
    if ($@){
      warning("Unable to store gene!!\n$@");
      $fails++;
    }
    $total++;
  }
  if ($fails > 0) {
    throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }
}


1;
