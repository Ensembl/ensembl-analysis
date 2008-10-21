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

Bio::EnsEMBL::Analysis::Runnable::HaplotypeProjection

=head1 SYNOPSIS

# This is the main analysis database

    my $genebuilder = new Bio::EnsEMBL::Analysis::Runnable::HaplotypeProjection
      (       
       '-hap_slice' => $self->query,
       '-slice'   => $self->target,
       '-input_id' => $self->input_id,
      );



=head1 DESCRIPTION

This module aligned the genomic sequence of a haplotype region and the corresponding
reference chromosome region and projects the gene annotations

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Runnable::HaplotypeProjection;

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

use Bio::EnsEMBL::Analysis::Config::HaplotypeProjection qw (
                                                 
                                                    );

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Root);


############################################################

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($hap_slice,$slice,$input_id) = $self->_rearrange([qw(HAP_SLICE SLICE INPUT_ID)],
					      @args);

    $self->throw("Must input a reference slice and a haplotype slice to HaplotypeProjection") unless (defined($slice) && defined($hap_slice));
    $self->{_final_genes} = [];
    $self->{_gene_types}  = [];

    $self->query($hap_slice);
    $self->target($slice);

    #$self->gene_types($GB_ENSEMBL_INPUT_GENETYPE);
  
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


=head2 create_alignment

 Description: retrieves ensembl and havana gene annotations with supporting evidence. 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub create_alignment {
my $target_slice = $self->target;
my $query_slice =  $self->query;

}

############################################################


=head2 filter_alignment

 Description: retrieves ensembl and havana gene annotations with supporting evidence. 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub filter_alignment {

}

############################################################


=head2 make_map_regions

 Description: retrieves ensembl and havana gene annotations with supporting evidence. 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub make_map_regions {

} 

############################################################


=head2 project_genes

 Description: retrieves ensembl and havana gene annotations with supporting evidence. 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub project_genes {

# Call get genes on target slice to get the genes to project.

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

###########################################################c

=head2 cluster_Transcripts

 Description : It separates transcripts according to strand and then clusters 
               each set of transcripts by calling _cluster_Transcripts_by_genomic_range()
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster

=cut


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

sub discarded_db{
  my ($self, $discarded_db) = @_;

  if ( $discarded_db ){
    $self->{_discarded_db} = $discarded_db;;
  }

  return $self->{_discarded_db};
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

sub target {
  my ($self,$slice) = @_;
  
  if (defined($slice)) {
    $self->{_target} = $slice;
  }
  return $self->{_target};
}

#fetches sequence from appropriate database

sub fetch_sequence{
  my ($self, $name, $db) = @_;

  my $sa = $db->get_SliceAdaptor; 

  my $slice = $sa->fetch_by_name($name);

  return $slice;
}

1;
