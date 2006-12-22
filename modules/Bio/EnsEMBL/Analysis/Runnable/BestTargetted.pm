# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
=head1 NAME 

Bio::EnsEMBL::Analysis::Runnable::BestTargetted

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::BestTargetted->new(
      -query => $slice,
     );
  $runnable->run;
  my @genes = @{$runnable->output};

=head1 DESCRIPTION

BestTargetted combines gene-structures from two different targetted runs 
generating a new set containing the best genes from both sets. 


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut



package Bio::EnsEMBL::Analysis::Runnable::BestTargetted;

use strict;
use warnings; 
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(clone_Transcript);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);




sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

    
  my( $primary_logic_name, $secondary_logic_name, $seqfetcher, $verbose, $genes ) = 
      rearrange([qw(
                     PRIMARY_LOGIC_NAME
                     SECONDARY_LOGIC_NAME
                     SEQFETCHER
                     VERBOSE
                     GENES
                    )], @args);
  
  $self->seqfetcher($seqfetcher) if defined $seqfetcher;
  $self->primary($primary_logic_name) if defined $primary_logic_name;
  $self->secondary($secondary_logic_name) if defined $secondary_logic_name;
  $self->verbose($verbose) if defined $verbose;
  $self->all_genes($genes) if defined $genes;

  return $self ; 
}


=head2 run

   Arg : none
   Function  : Runs BestTargetted
               - compares overlapping transcripts from both sets, constructed from the same evidence 
               and stores the best model
   Returnval : none 

=cut  


sub run { 
  my ($self) = @_ ; 
  
  my %outgenes;
  
   
  # get all genes using both logic names on slice and cluster them
 
  my $set1_name = 'primary';
  my $set2_name = 'secondary';
  
  my $set1_logic = $self->$set1_name;
  my $set2_logic = $self->$set2_name;
  
  
  my %logic_hash;
  $logic_hash{$set1_name} = [$set1_logic];
  $logic_hash{$set2_name} = [$set2_logic];
  
     
  my @allgenes = @{$self->all_genes};
 
  my ($clusters, $non_clusters) = cluster_Genes(\@allgenes, \%logic_hash ) ; 

  if ($self->verbose){
    print @$non_clusters ." non_clusters and ".@$clusters ." clusters\n";
  }
  
  foreach my $non_cluster(@$non_clusters){
    
    my @genes = $non_cluster->get_Genes();
    
    foreach my $gene (@genes){
      
      if ($self->verbose){
        print "storing non_cluster_gene ".$gene->dbID."\n";
      }
      
      $outgenes{$gene} = make_outgenes($self->analysis, $gene);
    }
  }
  
  my @twoways;  

  foreach my $cluster (@$clusters) {
    my @inc_sets = @{$cluster->get_sets_included};
    
   
    if (scalar @inc_sets == 2) {
      push @twoways, $cluster;
    
    } elsif (scalar @inc_sets == 1) { 
      foreach my $gene ($cluster->get_Genes_by_Set($inc_sets[0])) {
        
        if ($self->verbose){
          print "storing single_set_cluster_gene ".$gene->dbID." ".$inc_sets[0]."\n";
        }
        $outgenes{$gene} = make_outgenes($self->analysis, $gene);
      }
    } 
  }
  
  #print "Got " . scalar(@twoways) . " twoways\n";


  CLUSTER: foreach my $cluster (@twoways) {
    my @genes = $cluster->get_Genes_by_Set($set1_name);
    my @comp_genes = $cluster->get_Genes_by_Set($set2_name);

    my %hadmatch;
    my %unmatched_pairs = ();
      

    my $cluster_match = 0;
    foreach my $gene (@genes) {
      
      foreach my $comp_gene (@comp_genes) {
        
        TRANS: foreach my $trans (@{$gene->get_all_Transcripts}) {
		      
          #find the protein this transcript was made from:
          my $protein = undef;
          my @tsfs = @{$trans->get_all_supporting_features};
    	    if (scalar(@tsfs) == 1) {
              $protein = $tsfs[0]->hseqname;
    	    }else{
            print STDERR "warning transcript".$trans->dbID  ." has ".scalar(@tsfs)." tsfs\n";
            exit;
          }
          

          if (!defined($trans->translation)) {
            print STDERR $trans->stable_id." does not translate\n";
				    next TRANS;
          }
          my @exons = sort {$a->start <=> $b->start} @{$trans->get_all_translateable_Exons};


          COMPTRANS: foreach my $comp_trans (@{$comp_gene->get_all_Transcripts}) {
            #only want to compare transcripts made from the same protein:
            
            my $comp_protein = undef;
            my @comp_tsfs = @{$comp_trans->get_all_supporting_features};
    	      if (scalar(@comp_tsfs) == 1) {
                $comp_protein = $comp_tsfs[0]->hseqname;
    	      }else{
              print STDERR "warning transcript".$comp_trans->dbID  ." has ".scalar(@comp_tsfs)." tsfs\n";
              exit;
            }
            
            next if ($protein ne $comp_protein); #only compare transcripts made from same protein
            
            
            #only compare those which overlap
            unless ($comp_trans->start <= $trans->end && $comp_trans->end >= $trans->start){
              next COMPTRANS;
             
            }
            
            

            if (!defined($comp_trans->translation)) {
              print STDERR $comp_trans->stable_id ."does not have a translation\n";
				      next COMPTRANS;
            }
            my @comp_exons = sort {$a->start <=> $b->start} @{$comp_trans->get_all_translateable_Exons};

            if (scalar(@comp_exons) != scalar(@exons)){
              $unmatched_pairs{$trans}{'comp_gene'} = $comp_gene;
              $unmatched_pairs{$trans}{'protein'} = $protein;
              $unmatched_pairs{$comp_trans}{'gene'} = $gene;
              next COMPTRANS;
            }


            foreach my $exon (@exons) {
              my $comp_exon = shift @comp_exons;


              if ($exon->start != $comp_exon->start ||
                $exon->end != $comp_exon->end ||
                ($exon->strand != $comp_exon->strand) ) {
                
                #have 2 transcripts made from the same protein but they do not match - store these for further investigation
                
                $unmatched_pairs{$trans}{'comp_gene'} = $comp_gene;
                $unmatched_pairs{$trans}{'protein'} = $protein;
                $unmatched_pairs{$comp_trans}{'gene'} = $gene;
                next COMPTRANS;
              }
            }
            
            if ($self->verbose){
              print "storing trans ".$trans->dbID." exact_match\n";
              print "exact_match ".$trans->dbID." trans and ".$comp_trans->dbID." comp_trans\n";
            }
            
            if (! exists $outgenes{$gene}){
              $outgenes{$gene} = make_outgenes($self->analysis, $gene, $trans);
            }
				    
            
            $hadmatch{$trans->dbID} = 1;
            $hadmatch{$comp_trans->dbID} = 1;
          }
        }
      }
      
      
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        
        if (!$hadmatch{$trans->dbID}) {
        
          if (exists $unmatched_pairs{$trans}){
          
            
            my $protein_id = $unmatched_pairs{$trans}{'protein'};
            my $comp_gene = $unmatched_pairs{$trans}{'comp_gene'};
            my $comp_trans = ${$comp_gene->get_all_Transcripts}[0];

            #compare the translations of the 2 transcripts to the protein sequence:

            my $seq = $trans->translate->seq;
            my $comp_seq = $comp_trans->translate->seq;
            my $pep = $self->seqfetcher->get_Seq_by_acc($protein_id)->seq;


            if ($seq eq $pep){
              
              if ($self->verbose){
                print "storing trans ".$trans->dbID." trans_matches_protein_seq\n";
                print "trans_matches_protein_seq ".$trans->dbID." trans and ".$comp_trans->dbID." comp_trans\n";
              }
              if (! exists $outgenes{$gene}){
                $outgenes{$gene} = make_outgenes($self->analysis, $gene, $trans);
              }
	

            }elsif ($comp_seq eq $pep){
              if ($self->verbose){
                print "storing comp_trans ".$comp_trans->dbID." comp_trans_matches_protein_seq\n";
                print "comp_trans_matches_protein_seq ".$trans->dbID." trans and ".$comp_trans->dbID." comp_trans\n";
              }
              if (! exists $outgenes{$gene}){
                $outgenes{$gene} = make_outgenes($self->analysis, $comp_gene, $comp_trans);
              }

            }elsif ($seq eq $comp_seq){
              if ($self->verbose){
                print "storing trans ".$trans->dbID." trans_matches_comp_trans_seq\n";
                print "trans_matches_comp_trans_seq ".$trans->dbID." trans and ".$comp_trans->dbID." comp_trans\n";
              }
              if (! exists $outgenes{$gene}){
                $outgenes{$gene} = make_outgenes($self->analysis, $gene, $trans);
              }
	
            }else{

              #the seqs are different and don't match the protein - use exonerate to determine which is more similar
              
              my $target =  Bio::Seq->new( -display_id => $protein_id, -seq => $pep);
              my $query1 =  Bio::Seq->new( -display_id => $trans->dbID, -seq => $seq);
              my $query2 =  Bio::Seq->new( -display_id => $comp_trans->dbID, -seq => $comp_seq);
              
              
              my $gene_score = shift @{$self->run_exonerate($query1, $target)};
              my $comp_gene_score = shift @{$self->run_exonerate($query2, $target)};
              

              if ($gene_score >= $comp_gene_score){
                if ($self->verbose){
                  print "storing trans ".$trans->dbID." exonerate_score\n";
                  print "exonerate ".$trans->dbID." ".$gene_score." and ".$comp_trans->dbID." ".$comp_gene_score."\n";
                }
                if (! exists $outgenes{$gene}){
                  $outgenes{$gene} = make_outgenes($self->analysis, $gene);
                }
              }else{
                if ($self->verbose){
                  print "storing comp_trans ".$comp_trans->dbID." exonerate_score\n";
                  print "exonerate ".$trans->dbID." ".$gene_score." and ".$comp_trans->dbID." ".$comp_gene_score."\n";
                }
                if (! exists $outgenes{$comp_gene}){
                  $outgenes{$comp_gene} = make_outgenes($self->analysis, $comp_gene);
                }
              }
            }
          }else{
            if ($self->verbose){
              print "storing trans ".$trans->dbID." unpaired_transcript\n";
            }
            if (! exists $outgenes{$gene}){
              $outgenes{$gene} = make_outgenes($self->analysis, $gene, $trans);
            }
          }
        }
      }
    }
      #repeat to identify comp_trans which are in a cluster but there is no trans made from same thing
     
    foreach my $comp_gene (@comp_genes) {
      foreach my $comp_trans (@{$comp_gene->get_all_Transcripts}) {
        if (!$hadmatch{$comp_trans->dbID} && !exists $unmatched_pairs{$comp_trans}) {
          
          if ($self->verbose){
            print "storing comp_trans ".$comp_trans->dbID." unpaired_transcript\n";
          }
          if (! exists $outgenes{$comp_gene}){
            $outgenes{$comp_gene} = make_outgenes($self->analysis, $comp_gene, $comp_trans);
          }
        }
      }  
    }    
  }
  
  my @outgenes;
  foreach my $gene ( keys %outgenes ) {
    push @outgenes, $outgenes{$gene};
  }
  $self->output(\@outgenes);
  
  
  return ; 
}


=head2 make_outgenes

  Arg [0]   : ref to a Bio::EnsEMBL::Analysis 
  Arg [1]   : ref to a Bio::EnsEMBL::Gene 
  Arg [2]   : ref to a Bio::EnsEMBL::Transcript (optional)
  Function  : Clones the gene object and assigns the new analysis to it
  Returntype: new Bio::EnsEMBL::Gene 
  Example   : 

=cut




sub make_outgenes {
  my ($analysis, $gene, $trans) = @_;
  
  my $outgene = copy_empty_gene($analysis, $gene);

  if (not defined $trans){
    my @trans = @{$gene->get_all_Transcripts};
  
    if (scalar @trans > 1){
      print STDERR "gene ".$$gene->dbID." has >1 transcript\n";
    }else{
      $trans = $trans[0];
    }
  }

  my $outtrans = clone_Transcript($trans); 
  $outtrans->analysis($analysis);
  $outtrans->dbID(undef);
  
 
  #update the transcript_supporting_feature analysis_ids to match the transcript
  foreach my $sf(@{$outtrans->get_all_supporting_features}){
    $sf->dbID(undef);
    $sf->analysis($analysis);
  }
  
  #update the supporting_feature analysis_ids to match the transcript
  foreach my $exon (@{$outtrans->get_all_Exons}){
    $exon->dbID(undef);
    foreach my $sf (@{$exon->get_all_supporting_features}){
      $sf->dbID(undef);
      $sf->analysis($analysis);
    }
  }
  
  $outgene->add_Transcript($outtrans);
 
  return $outgene;
}

=head2 copy_empty_gene

  Arg [0]   : ref to a Bio::EnsEMBL::Analysis
  Arg [1]   : ref to a Bio::EnsEMBL::Gene 
  Function  : makes a gene and sets biotype and analysis
  Returntype: new Bio::EnsEMBL::Gene

=cut



sub copy_empty_gene {
  my ($analysis, $gene) = @_;
  
  my $newgene = new Bio::EnsEMBL::Gene;
  $newgene->biotype( $gene->biotype);
  $newgene->analysis($analysis);
  
  return $newgene;
}

=head2 run_exonerate

  Arg [0]   : query  Bio::Seq object
  Arg [1]   : target Bio::Seq object
  Function  : calls exonerate->run
  Returntype: score - integer

=cut
      


sub run_exonerate {
  my ($self, $query, $target) = @_;
  
  my $exonerate = Bio::EnsEMBL::Analysis::Runnable::BestTargettedExonerate->new
      (
       -program     => $self->program,
       -analysis    => $self->analysis,
       -query_seqs  => [$query],
       -target_seqs => [$target],
       -options     => "--model affine:local",
       );
       
  
  $exonerate->run;
  my ($score) = $exonerate->output;
  
  return $score;
}

########################
sub seqfetcher {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_seqfetcher} = $val;
  }

  return $self->{_bt_seqfetcher};
}

sub primary {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_primary} = $val;
  }

  return $self->{_bt_primary};
}

sub secondary {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_secondary} = $val;
  }

  return $self->{_bt_secondary};
}

sub verbose {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_verbose} = $val;
  }

  return $self->{_bt_verbose};
}

sub all_genes {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_all_genes} = $val;
  }

  return $self->{_all_genes};
}




##########################################
### local class
##########################################

package Bio::EnsEMBL::Analysis::Runnable::BestTargettedExonerate;

use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseExonerate);


sub target_type {
  my ($self) = @_;

  return 'protein';
}


sub query_type {
  my ($self) = @_;
  
  return 'protein';
}



sub parse_results {
  my ($self, $fh) = @_;

  my $score = 0;
  
  while (<$fh>){
    if ($_=~/^RESULT:/){
      my @tmp = split/\s+/, $_;
      if ($score < $tmp[9]){
        $score = $tmp[9]; #return the best score
      }
    }
  }  

  $self->output([$score]);
}


1; 

