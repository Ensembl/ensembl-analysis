=head1 LICENSE

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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::UTR_Builder - 

=head1 SYNOPSIS

my $utrbuilder_runnable = new Bio::EnsEMBL::Analysis::RunnableDB::UTR_Builder(
  -db        => $db,
  -input_id  => $input_id
 );
 $utrbuilder_runnable->fetch_input();
 $utrbuilder_runnable->run();
 $utrbuilder_runnable->output();
 $utrbuilder_runnable->write_output(); #writes to DB

=head1 DESCRIPTION

This is the new version of the UTR-addition procedure.
It combines predictions made from proteins with cDNA alignments to add UTR regions to the
gene models. It can also inlcude ESTs and ditags. It uses code from Coalescer/Consensus to produce
score for the alternative models and chose the best option.
It also includes ("look-for-both") code to correct the phases of the transcripts unless they are "blessed"
and inculdes the option to check for predefined protein/cDNA pairing as a first step,
looking for NM/NPentries in a GeneBank file.

Config files to set-up are
   Bio::EnsEMBL::Analysis::Config::GeneBuild::UTR_Builder
   Bio::EnsEMBL::Analysis::Config::Databases
   Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus (just a copy of the example file)
   Bio::EnsEMBL::Analysis::Config::GeneBuild::KillListFilter

=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are
usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::UTR_Builder;

use warnings ;
use vars qw(@ISA);
use strict;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw( throw warning verbose);
use Bio::EnsEMBL::Analysis::RunnableDB;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(
								clone_Exon
								transfer_supporting_evidence
								validate_Exon_coords
							       );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(
								      is_Transcript_sane
								      all_exons_are_valid
								      intron_lengths_all_less_than_maximum
								      set_start_codon
								      set_stop_codon
								      clone_Transcript
								      has_no_unwanted_evidence
								     );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(
								       validate_Translation_coords
								       compute_translation
								       contains_internal_stops
								       print_Translation
								      );
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::KillList::KillList;

#CONFIG FILES
use Bio::EnsEMBL::Analysis::Config::GeneBuild::UTR_Builder;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::KillListFilter;

@ISA = qw ( Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild );


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::UTR_Builder
  Function  : instatiates object & check config file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::UTR_Builder

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($UTR_BUILDER_CONFIG_BY_LOGIC) ;
  $self->read_and_check_config_value($TRANSCRIPT_CONSENSUS_CONFIG_BY_LOGIC, undef,['MIN_CONSENSUS', 'END_EXON_PENALTY', 'EST_OVERLAP_PENALTY', 'SHORT_INTRON_PENALTY', 'SHORT_EXON_PENALTY', 'UTR_PENALTY', 'VERBOSE']) ;

  return $self;
}

my $totalgenes = 0;

=head2 fetch_input

  Arg [1]    : none
  Description: Get all raw data needed from the databases
  Returntype : none

=cut

sub fetch_input{
  my ($self) = @_;

  my $slice = $self->fetch_sequence(undef, $self->db);
  if(!$slice){ throw "can't fetch input slice!\n" }
  $self->query($slice);

  my @databases = @{ $self->INPUT_DBS } ;
  # fetch all general input genes

  foreach my $db ( @databases) {
    my $dba = $self->get_dbadaptor($db) ;
    $dba->dnadb($self->db) ;
    my $slice = $self->fetch_sequence($self->input_id, $dba ) ;
    foreach my $input_genetype (@{ $self->INPUT_GENETYPES }) {
      my $input_genes = $slice->get_all_Genes_by_type($input_genetype);
      print STDERR "got " . scalar(@{$input_genes}) . " $input_genetype genes [ ".
      $dba->dbc->dbname() . "@" . $dba->dbc->host . " ]\n";
      $self->gw_genes( $input_genes );
    }
  }

  # get blessed genes
  my $blessed_slice;
  my @blessed_genes;
  my $blessed_type;
  if($self->BLESSED_DB and scalar(@{$self->BLESSED_GENETYPES})){
	my $blessed_db = get_db_adaptor_by_string($self->BLESSED_DB,1);
    if ($self->BLESSED_DB){
      #fetch blessed genes from blessed db
      $blessed_slice = $blessed_db->get_SliceAdaptor->fetch_by_name($self->input_id);
    }
    else{
      #fetch blessed genes from gw db
      $blessed_slice = $self->query;
    }
    foreach my $bgt( @{$self->BLESSED_GENETYPES} ){
      my $blessed_genes = $blessed_slice->get_all_Genes_by_type($bgt);
      print STDERR "got " . scalar(@{$blessed_genes}) . " $bgt genes [ ".
	$blessed_slice->db->dbname() . "@" . $blessed_slice->db->host ." ]\n";

      $self->blessed_genes( $blessed_slice, $blessed_genes );
      $blessed_type .= $bgt."";
    }
    # store all blessed type names for VIP treatment
    $self->{'blessed_type'} = ($blessed_type.$self->BLESSED_UTR_GENETYPE);
  $self->BLESSED_DB->dbc->disconnect_when_inactive(1);
  }

  # are there any genes here at all?
  if(!scalar @{$self->gw_genes} and !scalar @{$self->blessed_genes}){
    warn ("No genewise or blessed genes fetched from database(s). Are you sure there should no genes?\n");
    return 0;
  }

  # get cdnas
  if(!(defined($self->cDNA_GENETYPE) and $self->CDNA_DB)){
    warn("cDNA_GENETYPE and CDNA_DB are both undefined in your config. Are you sure you are not using any cDNAs as UTR evidence!?\n");
  }
  else{
    foreach my $dbname (@{$self->CDNA_DB}) {
       my $dbh = $self->get_dbadaptor($dbname) ;
       my $cdna_vc = $dbh->get_SliceAdaptor->fetch_by_name($self->input_id);
       $self->_cdna_slice($cdna_vc);
       foreach my $cdna_type (@{$self->cDNA_GENETYPE}){
          my @cdna_genes = @{$self->_cdna_slice->get_all_Genes_by_type($cdna_type)};
          print STDERR "got ".scalar(@cdna_genes)." ".$cdna_type." cDNAs.\n" if $self->VERBOSE;
          $dbh->dbc->disconnect_when_inactive(1);

      # filter cdnas
          my $filtered_cdna = $self->_filter_cdnas(\@cdna_genes, 0);
      
          $self->cdna_genes($filtered_cdna);
          print STDERR "got " . scalar(@{$filtered_cdna}) . " cDNAs after filtering.\n" if $self->VERBOSE;
       }
    }
    print STDERR "got " . scalar(@{$self->cdna_genes}) . " filtered cDNAs accross all databases.\n";
  }

  # get ESTs
  if(defined($self->EST_GENETYPE) && $self->EST_DB){
    my $est_vc    = $self->EST_DB->get_SliceAdaptor->fetch_by_name($self->input_id);
    my @est_genes = @{$est_vc->get_all_Genes_by_type($self->EST_GENETYPE)};
    print STDERR "got " . scalar(@est_genes) . " ".$self->EST_GENETYPE." ESTs.\n";
    if($self->FILTER_ESTS){
      my $filtered_ests = $self->_filter_cdnas(\@est_genes, 1);
      $self->ests($filtered_ests);
    }
    else{
      $self->ests(\@est_genes);
    }
    print STDERR "got " . scalar(@{$self->ests}) . " filtered ESTs.\n";
  }

  # get ditags
  my ($dfa, $ditag_slice);
  my @ditags;
  if((scalar @{$self->DITAG_TYPE_NAMES}) && $self->DITAG_DB){
	my $ditag_db = get_db_adaptor_by_string($self->DITAG_DB,1);
    $dfa         = $ditag_db->get_DitagFeatureAdaptor;
    $ditag_slice = $ditag_db->get_SliceAdaptor->fetch_by_name($self->input_id);

    foreach my $ditag_type (@{$self->DITAG_TYPE_NAMES}) {
      my @type_ditags = @{$dfa->fetch_pairs_by_Slice($ditag_slice, $ditag_type)};
      print STDERR "got " . scalar(@type_ditags) . " ".$ditag_type." ditags.\n" if $self->VERBOSE;
      push(@ditags, @type_ditags);
    }
    $self->DITAG_DB->dbc->disconnect_when_inactive(1);
  }
  if(scalar @ditags){
    @ditags = sort {($a->{'start'} <=> $b->{'start'}) or ($a->{'end'} <=> $b->{'end'})} @ditags;
    $self->ditags(\@ditags);
    print STDERR "got " . scalar(@ditags) ." sorted ditags of all specified types.\n";
  }
  else{
    print STDERR "not using Ditags.\n";
  }
  # db disconnections
  foreach my $db ( @databases ) {
     my $dba = $self->get_dbadaptor($db) ;
     $dba->dnadb($self->db) ;
     $dba->dbc->disconnect_when_inactive(1) ;
  }

  if ($self->CDNA_DB) {
    foreach my $dbname (@{$self->CDNA_DB}) {
       my $dbh = $self->get_dbadaptor($dbname) ;
       $dbh->dbc->disconnect_when_inactive(1) ;
    }
  }
  $self->EST_DB->dbc->disconnect_when_inactive(1)
    if($self->EST_DB);

  # set evidence sets for Coalescer code
  my (@est_biotypes, @cdna_logicnames);
  push(@est_biotypes, $self->EST_GENETYPE);
  my @simgw_biotypes  = @{ $self->INPUT_GENETYPES };
  push(@cdna_logicnames, @{$self->cDNA_GENETYPE});
  $self->{evidence_sets} = {
			    'est'   => \@est_biotypes,
			    'simgw' => \@simgw_biotypes,
			    'cdna'  => \@cdna_logicnames,
			   };

  # prune genes during filtering?
  # (usually not used, as too many good things are thrown out.)
  $self->prune($self->PRUNE_GENES);

  # get rid of identical models?
  # (usually not used, as we want the same number of output as input genes)
  $self->{'remove_redundant'} = 0;

  #get killlist entries for cDNAs to ignore
  my $kill_list_type = "cDNA";
  my $kill_list = $self->populate_kill_list($kill_list_type);
  $self->kill_list($kill_list);

  if($self->LOOK_FOR_KNOWN){
    #read KnowUTR pairing into hash
    $self->create_predefined_pairing($self->KNOWNUTR_FILE);
  }

}


=head2 run

  Description: general run method
  Returntype : none

=cut

sub run {
  my ($self) = @_;

  # filter the genewise predictions to remove fragments.
  # Make sure we keep the fragments to let the
  # genebuilder do the final sorting out, but don't add UTRs to them.
  print STDERR "\nTOTAL GENES IN: ".(scalar @{$self->gw_genes} + scalar @{$self->blessed_genes})."\n";
  my $ininumber = scalar @{$self->gw_genes};
  my $filtered_genes = $self->filter_genes($self->gw_genes);
  print STDERR "filtered unblessed (presumably genewise) from $ininumber to ".(scalar @$filtered_genes). " genes.\n" if $filtered_genes;

  #for blessed genes, only check region and existing UTR
  $ininumber = scalar @{$self->blessed_genes};
  my $filtered_blessed_genes = $self->filter_genes($self->blessed_genes, 1);
  print STDERR "filtered blessed from $ininumber to ".(scalar @$filtered_blessed_genes). " genes.\n" if $filtered_blessed_genes;

  # match gene-model to cDNAs/ESTs for genewise genes by calling run_matching.
  # run_matching is called for blessed genes only if there are blessed genes in the first place.

  $self->run_matching($filtered_genes, $self->UTR_GENETYPE, 0);

  if (scalar @$filtered_blessed_genes) {
    $self->run_matching($filtered_blessed_genes, $self->BLESSED_UTR_GENETYPE,  1);
  }

  print STDERR "\n\n### Have altogether ".(scalar @{$self->combined_genes})." +UTR genes & ".
               (scalar @{$self->unmatched_genes})." noUTR genes for final output.\n";

  # check translation, start and stop
  my @remapped;
  push( @remapped, @{$self->remap_genes( $self->combined_genes,  undef )} );
  print STDERR "Remapped genes which have UTRs added.\n";
  push( @remapped, @{$self->remap_genes( $self->unmatched_genes, $self->EXTEND_BIOTYPE_OF_UNCHANGED_GENES )} );
  print STDERR "Remapped genes which can't have UTRs added.\n";

  #results are stored in "output"
  $self->output(@remapped);
  print "TOTAL GENES OUT: ".(scalar @remapped)."\n";
  if($self->VERBOSE){
    foreach ( @remapped ) {
      print "final biotype : " . $_->biotype . "\n";
    }
  }

}


=head2 filter_genes

  Arg [1]    : ref to array of Bio::EnsEMBL::Gene
  Description: We do not want to add UTRs to fragments of genes
  otherwise they get boosted in importance during the final gene build
  and produce rubbish alternative transcripts and in the worst cases
  can misjoin clusters especially if there is some
  misassembled/duplicated sequence. One solution is to cluster the
  genewises at this stage in much the same way as we do during the
  GeneBuilder run, and attempt to add UTRs ony to the best ones eg
  fullest genomic extent.
  Also checks for existing UTR, these genes are just passed on
  Returntype : ref to array of Bio::EnsEMBL::Gene

=cut

sub filter_genes {
  my ($self, $genesref, $blessed) = @_;

  my @tested_genes;
  my $pruned_CDS;

  #run tests: get rid of long-intron genes, etc.
 GENES:
  foreach my $test_gene (@$genesref){
    if( $test_gene ){
      my $test_transcript = $test_gene->get_all_Transcripts->[0];

      #would like to avoid slow sorting call - but the "get_all_Transcripts" does not do this job properly...
      #$test_transcript->sort;

      if($test_transcript->three_prime_utr or $test_transcript->five_prime_utr){
	print STDERR "Gene has UTR already. Skipping.\n" if $self->VERBOSE;
	$self->unmatched_genes($test_gene);
	next GENES;
      }

      # We test if the transcript is fully contained in the slice and also check transcript's sanity.
      # The treatment of blessed and unblessed models differ in the sanity checks.
      # Blessed models get a less stringent sanity checks (i.e. we don't test exon/intron lengths)
      # they're allowed to have sometimes weird structures.

      # Models which fall off a slice will be picked up by another slice, so they'll be accounted
      # for somehow. There's no need to push them into the unmatched_genes hash.

      # Models which didn't pass the sanity check test will be captured in the unmatched_genes hash.  No UTRs
      # will be added to them and they'll be written to the output DB as they were.
      

      if($blessed){
	if($test_transcript->start < 1 && $test_transcript->end > 1){
	  print STDERR "ignoring blessed gene at seq_region/chr ".$test_transcript->seq_region_name . " (" . 
                       $test_transcript->seq_region_start . "-" . $test_transcript->seq_region_end . 
                       ") as it falls off the slice by its lower end\n" if $self->VERBOSE;
	} 
	elsif ( is_Transcript_sane($test_transcript)
            && ($test_transcript->start < $self->query->length)
            && ($test_transcript->end > 1)
            && has_no_unwanted_evidence($test_transcript) ){

 	  push(@tested_genes, $test_gene);  # test_gene not falling off slice + is sane.
	} else {
          $self->unmatched_genes($test_gene); # test_gene doesn't fall off slice but isn't sane
        }
      } # end if $blessed

      elsif (!$blessed) {
        if($test_transcript->start < 1 && $test_transcript->end > 1){
          print STDERR "ignoring unblessed gene at seq_region/chr ".$test_transcript->seq_region_name . " (" . 
                        $test_transcript->seq_region_start . "-" . $test_transcript->seq_region_end . 
                        ") as it falls off the slice by its lower end\n" if $self->VERBOSE;
	}        
	elsif( is_Transcript_sane($test_transcript)
	    && all_exons_are_valid($test_transcript, $self->MAX_EXON_LENGTH)
	    && intron_lengths_all_less_than_maximum($test_transcript, $self->MAX_INTRON_LENGTH)
	    && ($test_transcript->start < $self->query->length)
	    && ($test_transcript->end > 1)
	    && has_no_unwanted_evidence($test_transcript) ) {

	  push(@tested_genes, $test_gene);  # test_gene not falling off slice + is sane.
	} else{
	  $self->unmatched_genes($test_gene);  # test_gene doesn't fall off slice but isn't sane
	}
      } # end if !$blessed
    }  # end test_gene/GENES block

    else{
      throw "\nERROR: No genes passed to the filter_genes method from the run method via genesref(?). Got no genes to filter!!\n\n";
    }
  }
  $genesref = \@tested_genes;

  #cluster genes (stores internally)
  if(!$blessed){
    # all the genes are single transcript at this stage, cluster them
    my $clustered_genes = $self->cluster_CDS($genesref);

    #puning if desired
    if($self->prune){
      $genesref = $self->prune_CDS($clustered_genes);
    }

  }

  return $genesref;
}


=head2 cluster_CDS

  Arg [1]    : ref to array of Bio::EnsEMBL::Gene
  Description: Separates CDSs according to strand and then clusters
               each set by calling Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils:cluster_Genes
  Returntype : ref to array of Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster

=cut

sub cluster_CDS {
  my ($self, $CDSref) = @_;

  my @forward_CDS;
  my @reverse_CDS;

  foreach my $gene (@$CDSref){
    my @exons = @{ $gene->get_all_Exons };
    if ( $exons[0]->strand == 1 ){
      push( @forward_CDS, $gene );
    }
    else{
      push( @reverse_CDS, $gene );
    }
  }


  my (@clusters, @forward_clusters, @reverse_clusters);
  my ($forward_clusters, $forward_clusters_rest);
  my ($reverse_clusters, $reverse_clusters_rest);

  if ( @forward_CDS ){
    ($forward_clusters, $forward_clusters_rest) = cluster_Genes( \@forward_CDS, $self->{evidence_sets} );
    push(@forward_clusters, @$forward_clusters);
    push(@forward_clusters, @$forward_clusters_rest);
  }
  if ( @reverse_CDS ){
    ($reverse_clusters, $reverse_clusters_rest) = cluster_Genes( \@reverse_CDS, $self->{evidence_sets} );
    push(@reverse_clusters, @$reverse_clusters);
    push(@reverse_clusters, @$reverse_clusters_rest);
  }


  #store the two cluster sets and return combined set
  if ( scalar @forward_clusters ){
    #store for later
    $self->forward_genewise_clusters(\@forward_clusters);
    push( @clusters, @forward_clusters);
  }
  if ( scalar @reverse_clusters ){
    #store for later
    $self->reverse_genewise_clusters(\@reverse_clusters);
    push( @clusters, @reverse_clusters);
  }

  return \@clusters;
}


=head2 prune_CDS

  Arg [1]    : ref to array of GeneClusters
  Description: rejects duplicate CDS, makes sure they are kept on
               unmatched_genes or they (and their supporting evidence)
               will be lost for the rest of the build
               Note: pruning will not be used unless specified otherwise
  Returntype : ref to array of Bio::EnsEMBL::Gene

=cut

sub prune_CDS {
  my ($self, $gene_cluster_ref) = @_;
  my $cluster_count = 0;
  my @pruned_transcripts;
  my @pruned_genes;

 CLUSTER:
  foreach my $gene_cluster ( @$gene_cluster_ref ){
    $cluster_count++;
    my @unpruned_genes;

    #check for multi-transcript genes
    foreach my $unpruned_gene (@{ $gene_cluster->get_Genes }){
      if(scalar(@{$unpruned_gene->get_all_Transcripts}) > 1){
        throw($unpruned_gene->dbID . " has >1 transcript - can't handle it yet \n");
      }
      push(@unpruned_genes, $unpruned_gene);
    }

    # sort the unpruned_genes by total exon length of their single transcripts - there are no UTRs to worry about yet
    @unpruned_genes = sort{$b->get_all_Transcripts->[0]->length <=> $a->get_all_Transcripts->[0]->length} @unpruned_genes;

    # do we really just want to take the first transcript only?
    # If there's a very long single exon gene we will lose any underlying multi-exon transcripts
    # this may increase problems with the loss of valid single exon genes as mentioned below. 
    # it's a balance between keeping multi exon transcripts and losing single exon ones
    my $maxexon_number = 0;

    foreach my $gene (@unpruned_genes){
      my $exon_number = scalar(@{$gene->get_all_Transcripts->[0]->get_all_Exons});
      if ( $exon_number > $maxexon_number ){
	$maxexon_number = $exon_number;
      }
    }

    if ($maxexon_number == 1){ # ie the longest transcript is a single exon one
      # take the longest:
      push (@pruned_genes, $unpruned_genes[0]);
      print STDERR "VAC: found single_exon_transcript: " .$unpruned_genes[0]->dbID. "\n" if $self->VERBOSE;
      shift @unpruned_genes;
      $self->unmatched_genes(@unpruned_genes);
      next CLUSTER;
    }

    # otherwise we need to deal with multi exon transcripts and reject duplicates.
    # links each exon in the transcripts of this cluster with a hash of other exons it is paired with
    my %pairhash;

    # allows retrieval of exon objects by exon->id - convenience
    my %exonhash;

    # prune redundant transcripts

  GENE:
    foreach my $gene (@unpruned_genes) {

      my @exons = @{$gene->get_all_Transcripts->[0]->get_all_Exons};

      my $i     = 0;
      my $found = 1;

      # 10.1.2002 VAC we know there's a potential problem here - single exon transcripts which are in a
      # cluster where the longest transcript has > 1 exon are not going to be considered in
      # this loop, so they'll always be marked "transcript already seen"
      # How to sort them out? If the single exon overlaps an exon in a multi exon transcript then
      # by our rules it probably ought to be rejected the same way transcripts with shared exon-pairs are.
      # Tough one.

    EXONS:
      for ($i = 0; $i < $#exons; $i++) {
	my $foundpair = 0;
	my $exon1 = $exons[$i];
	my $exon2 = $exons[$i+1];

        # go through the exon pairs already stored in %pairhash. 
        # If there is a pair whose exon1 overlaps this exon1, and 
        # whose exon2 overlaps this exon2, then these two transcripts are paired

      EXONHASH:
        foreach my $first_exon_id (keys %pairhash) {
          my $first_exon = $exonhash{$first_exon_id};

          foreach my $second_exon_id (keys %{$pairhash{$first_exon}}) {
            my $second_exon = $exonhash{$second_exon_id};

            if ( $exon1->overlaps($first_exon) && $exon2->overlaps($second_exon) ) {
              $foundpair = 1;
	      last EXONHASH;
              # this method allows a transcript to be covered by exon pairs
              # from different transcripts, rejecting possible
              # splicing variants. Needs rethinking
		
            }
          }
        }

	
	if ($foundpair == 0) {	# ie this exon pair does not overlap with a pair yet found in another transcript
	  $found = 0;		# ie currently this transcript is not paired with another

	  # store the exons so they can be retrieved by id
	  $exonhash{$exon1} = $exon1;
	  $exonhash{$exon2} = $exon2;

	  # store the pairing between these 2 exons
	  $pairhash{$exon1}{$exon2} = 1;
	}
      }				# end of EXONS

      # decide whether this is a new transcript or whether it has already been seen
      if ($found == 0) {
	push(@pruned_genes, $gene);
      }
      elsif ($found == 1 && $#exons == 0){
        $self->unmatched_genes($gene);
      }
      else {
      # dont know what this is about...
        $self->unmatched_genes($gene);
	if ( $gene == $unpruned_genes[0] ){
	}
      }
    } # end of this gene
  } #end CLUSTER

  return \@pruned_genes;
}


=head2 run_matching

  Arg [1]    : ref to array of Bio::EnsEMBL::Gene
  Arg [2]    : string representing genetype to be associated with UTR-modified genes
  Arg [3]    : bool indicating that we are NOT trying to map to predefined cDNA (blessed genes)
  Description: main function which matches CDS (genewise) predictions to cDNA alignments
               considering only terminal CDS exons, looking for "knowns" first, optionally
               using ESTs and ditags.
  Returntype : none

=cut

sub run_matching{
  my ($self, $genesref, $combined_genetype, $blessed) = @_;

  # merge exons with frameshifts into a big exon
  my @merged_genes = ();

  #maybe we should not do merging with blessed genes at all?!
  @merged_genes = $self->_merge_genes($genesref, $blessed);
    
  print STDERR "\n### Got ". scalar(@merged_genes) . " blessed genes which have been filtered and frameshift-exon-merged. " if ($blessed);
  print STDERR "\n### Got ". scalar(@merged_genes) . " unblessed genes which have been filtered and frameshift-exon-merged. " if (!$blessed);
  print STDERR "Working on them one by one ...\n";
  # sort genewises by exonic length and genomic length
  @merged_genes = sort {
                         my $result = ( ($b->get_all_Transcripts->[0]->end - $b->get_all_Transcripts->[0]->start + 1) <=>
					($a->get_all_Transcripts->[0]->end - $b->get_all_Transcripts->[0]->start + 1) );

			 unless ($result){
			   return ( ($b->get_all_Transcripts->[0]->length) <=>
				    ($a->get_all_Transcripts->[0]->length) );
			 }
			 return $result;
		       } @merged_genes;

  # find one-2-one matching between proteins and cdnas
 CDS:
  foreach my $cds (@merged_genes){

    print STDERR "\n-------------- next cds: ".$cds->seq_region_start."-".$cds->seq_region_end." ---------------\n";

    # should be only 1 transcript
    my @cds_tran  = @{$cds->get_all_Transcripts};			
    my @cds_exons = @{$cds_tran[0]->get_all_Exons}; # ordered array of exons
    my $strand    = $cds_exons[0]->strand;
    my $cdna_match;
    my $usingKnown = 0;

    if($cds_exons[$#cds_exons]->strand != $strand){
      warn("first and last cds exons have different strands - can't make a sensible combined gene\n");
      # get and store unmerged version of cds
      my $unmerged_cds = $self->retrieve_unmerged_gene($cds);
      $self->unmatched_genes($unmerged_cds);
      next CDS;
    }

    #Look for a pre-defined pairing between a protein and a cDNA the gene was build from.
    #for a genewise model, look for pre-defined pairing in the KNOWN_UTR gpff file.
    #for a blessed model, look for pre-defined pairing via the NM_* xrefs in the xref table of BLESSED_DB

    my $predef_match = undef;
    my $combined_transcript = undef;

    if($self->LOOK_FOR_KNOWN){

      $predef_match = $self->check_for_predefined_pairing($cds, $blessed);

      if(defined $predef_match){

	$combined_transcript = $self->combine_genes($cds, $predef_match, $blessed);
        if (defined $combined_transcript) {
          #$combined_transcript->sort;

          # just check combined transcript works before throwing away the original transcript
          if ($blessed) {
            if (  $combined_transcript
	      && is_Transcript_sane($combined_transcript)
	      && validate_Translation_coords($combined_transcript, 1)
	      && !contains_internal_stops($combined_transcript)
	      && $combined_transcript->translate
	      && has_no_unwanted_evidence($combined_transcript) ){

	    # make sure combined transcript doesn't misjoin any genewise clusters
	      if($self->find_cluster_joiners($combined_transcript)){
                print STDERR "Found a cluster_joiner! Can't use this blessed combined transcript although it's made from predefined cDNA: " 
                              if $self->VERBOSE;
                print STDERR " start " . $combined_transcript->seq_region_start. " end ".$combined_transcript->seq_region_end."\n\n" if $self->VERBOSE;
	        $combined_transcript = undef;   # combined_transcript not used
 	      } else {
                print STDERR "Blessed combined transcript with predefined cDNA is sane and isn't cluster_joiner.\n" if $self->VERBOSE; 
              }
            } else {
              print STDERR "Blessed combined transcript at seq_region/chr ". $combined_transcript->seq_region_name . " " . $combined_transcript->start .
                           "-" . $combined_transcript->end . " did not pass sanity check.\n\n" if $self->VERBOSE;
              $combined_transcript = undef;
            }
          } # end if ($blessed), checking blessed combined transcript

          elsif (!$blessed) {
            if (  $combined_transcript
	      && is_Transcript_sane($combined_transcript)
	      && all_exons_are_valid($combined_transcript, $self->MAX_EXON_LENGTH, 1)
	      && intron_lengths_all_less_than_maximum($combined_transcript, $self->MAX_INTRON_LENGTH)
	      && validate_Translation_coords($combined_transcript, 1)
	      && !contains_internal_stops($combined_transcript)
	      && $combined_transcript->translate
	      && has_no_unwanted_evidence($combined_transcript) ){

	    # make sure combined transcript doesn't misjoin any genewise clusters
	      if($self->find_cluster_joiners($combined_transcript)){
                print STDERR "Found a cluster_joiner! Can't use this unblessed combined transcript although it's made from predefined cDNA: " 
                             if $self->VERBOSE;
                print STDERR " start " . $combined_transcript->seq_region_start. " end ".$combined_transcript->seq_region_end."\n\n" if $self->VERBOSE;
                $combined_transcript = undef;
	      } else {
                print STDERR "Unblessed combined transcript with predefined cDNA is sane and isn't cluster_joiner.\n" if $self->VERBOSE;
              }
            }
	    else {
	      print STDERR "Unblessed combined transcript at seq_region/chr". $combined_transcript->seq_region_name . " " . $combined_transcript->start .
                           "-" . $combined_transcript->end . " did not pass sanity check.\n\n" if $self->VERBOSE;
	      $combined_transcript = undef;
            }
          } 
        }
        else {
            print STDERR "Predefined cDNA was found but can't be used as there were problems combining it with protein-coding region. ".
                         "Will fall back to more generic evidence.\n";
        }
      }
    }
   
    # If at this point a pre-defined cDNA gene model is found for the coding model
    # and a combined transcript is successfully made, we move on to the next coding
    # model. Otherwise, we need to look at my generic evidence (other cDNAs or ESTs
    # in the vicinity of the coding model) and score the evidence to choose the best UTR
    # evidence.
 
    if($predef_match && $combined_transcript){
      print STDERR "Used predefined cDNA successfully!\n" if $self->VERBOSE;
      $cdna_match = $predef_match;
      $usingKnown = 1;
    }
    else{

      # find matching cdnas using scoring of all evidence available
      my ($matching_cdnas, $utr_length_hash, $UTR_side_indicator_hash) = $self->match_protein_to_cdna($cds, 0);
      # Note: $utr_length_hash is never used once returned.

      if(!(scalar @$matching_cdnas)){
	warn("Could not identify any matching cDNA!\n");        
        my $unmerged_cds;
        # try to use ESTs instead only if EST_DB is defined
	if ($self->EST_DB) {
	  #we could use coalescer code to produce joined ESTs (find longest path, etc.)?  Not implemented yet.
	  ($matching_cdnas, $utr_length_hash, $UTR_side_indicator_hash) = $self->match_protein_to_cdna($cds, 1);
	  if(!(scalar @$matching_cdnas)){
            $unmerged_cds = $self->retrieve_unmerged_gene($cds);
	    warn("Could not identify any matching ESTs either!\n");
            $self->unmatched_genes($unmerged_cds);
	    next CDS;
          }
        } else {
          $unmerged_cds = $self->retrieve_unmerged_gene($cds);
          $self->unmatched_genes($unmerged_cds);
          next CDS;
        }
      }

      ## scoring code...

      print STDERR "Now clustering UTR evidence and scoring them using TranscriptConsensus.\n";

      #convert genes to extended objects
      my $matching_extended_cdnas = $self->convert_to_extended_genes($matching_cdnas);

      #cluster matching cDNA or EST genes
      #call cluster method from ClusterUtils
      # print "EXTENDED : ",join('   ',@{$matching_extended_cdnas}),"\n";
      #foreach my $keys (keys %{$self->{evidence_sets}}){
      #  print "Evidence: ",@{$self->{evidence_sets}{$keys}},"\n";
      #}
      my ($clusters, $non_clusters) = cluster_Genes( $matching_extended_cdnas, $self->{evidence_sets} ) ;

      #store genes seperated by strand #USAGE??
      my $genes_by_strand;
      my @possible_transcripts;
      foreach my $gene (@$matching_extended_cdnas){
	push @{$genes_by_strand->{$gene->strand}}, $gene;
      }

      #apply collapsing method from TranscriptConsensus (creates the scores)
      my @cluster_to_use = ();
      if(scalar @$clusters){
	push(@cluster_to_use, @$clusters);
      }
      if(scalar @$non_clusters){
	push(@cluster_to_use, @$non_clusters);
      }

      foreach my $cluster (@cluster_to_use){
	print STDERR "EVIDENCE SET CLUSTER: ".$cluster->start." ".$cluster->end." ".$cluster->strand."\n" if $self->VERBOSE;
	my $tc = Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus->new (-analysis => $self->analysis) ;
        $tc->{min_consensus} = $self->MIN_CONSENSUS ;
        $tc->{end_exon_penalty} = $self->END_EXON_PENALTY ;
        $tc->{est_overlap_penalty} = $self->EST_OVERLAP_PENALTY ;
        $tc->{short_intron_penalty} = $self->SHORT_INTRON_PENALTY ;
        $tc->{short_exon_penalty} = $self->SHORT_EXON_PENALTY ;
        $tc->{utr_penalty} = $self->UTR_PENALTY ;
        $tc->{verbose} = $self->VERBOSE ;

        my $collapsed_cluster = $tc->collapse_cluster($cluster, $genes_by_strand);
	my $potential_genes = $cluster->get_Genes;
	foreach my $potential_gene ( @$potential_genes ){
	  foreach my $potential_trans (@{$potential_gene->get_all_Transcripts}) {

	    #create an extended transcript (with scores)
	    my @exon_array = ( @{$potential_trans->get_all_Exons} );
	    my $new_transcript = $tc->transcript_from_exons(\@exon_array, undef, $collapsed_cluster) ;

            #add supporting features
            $new_transcript->add_supporting_features(@{$potential_trans->get_all_supporting_features});

	    #add a ditag score
	    print STDERR "new collapsed transcript score :".$new_transcript->score()."\n" if $self->VERBOSE;
	    my ($ditag_score, $index) = $self->score_ditags($self->ditags, $new_transcript, 0);
	    $new_transcript->score($new_transcript->score() + $ditag_score);

	    #add a UTR length score
	    my $UTR_score = $self->calculate_UTR_score($cds_tran[0], $new_transcript);
	    $new_transcript->score($new_transcript->score() + $UTR_score);

	    push(@possible_transcripts, $new_transcript);
	  }
	}
      } #clusters

      #get the highest-scoring transcripts that survived collapsing and scoring

      @possible_transcripts = sort { $a->score <=> $b->score } @possible_transcripts if @possible_transcripts;
      print STDERR "\nGot ". scalar@possible_transcripts . " possible transcripts as UTR evidence after clustering and TranscriptConsensus collapsing.\n";

      my ($cdnaname, $genename);
      my $round = 0;
    POS:
      while(my $chosen_transcript = pop @possible_transcripts){

        print STDERR "Checking individual possible cDNA/EST transcript as UTR evidence: ".
                     ${$chosen_transcript->get_all_supporting_features}[0]-> hseqname ."..... \n" if $self->VERBOSE;
        my $chosen_feats = $chosen_transcript->get_all_supporting_features;

        foreach my $feature (@$chosen_feats) {
          print $feature->hseqname."\n";
        }

	#make a gene from the candidate cDNA (survivied TranscriptConsensus)
	$cdna_match = Bio::EnsEMBL::Gene->new;
	$cdna_match->slice($chosen_transcript->slice);
	$cdna_match->add_Transcript($chosen_transcript);
	$cdna_match->analysis($chosen_transcript->analysis);
        $cdna_match->get_all_Transcripts->[0]->add_supporting_features(@$chosen_feats);
	#combine (not predefined) candidate cDNA gene object with the cds gene object

	$combined_transcript = $self->combine_genes($cds, $cdna_match, $blessed);

	# Check combined transcript for sanity and non-cluster-joiners before throwing away the original transcript
	
	if ( defined($combined_transcript)){
	  #$combined_transcript->sort;
          if ($blessed) {
            if( is_Transcript_sane($combined_transcript) && has_no_unwanted_evidence($combined_transcript) ){
              # make sure combined transcript doesn't misjoin any genewise clusters
	      if($self->find_cluster_joiners($combined_transcript)){
	        print STDERR "Found a cluster_joiner! Can't use this blessed combined transcript: ";
	        print STDERR " start " . $combined_transcript->seq_region_start. " end ".$combined_transcript->seq_region_end."\n\n";
  	        #dont use this one
	        $combined_transcript = undef;
	      } else {
                print STDERR "Blessed combined transcript is sane and isn't cluster_joiner.\n";
              }
            } else{
              print STDERR "Blessed combined transcript at seq_region/chr". $combined_transcript->seq_region_name . " " . 
                           $combined_transcript->seq_region_start . "-" . $combined_transcript->seq_region_end . " did not pass sanity check.\n\n";
              $combined_transcript = undef;
	    }
          }
          elsif (!$blessed) {
            if(is_Transcript_sane($combined_transcript)
              && all_exons_are_valid($combined_transcript, $self->MAX_EXON_LENGTH, 1)
              && intron_lengths_all_less_than_maximum($combined_transcript, $self->MAX_INTRON_LENGTH)
              && has_no_unwanted_evidence($combined_transcript)
             ){
              if($self->find_cluster_joiners($combined_transcript)){
                print STDERR "Found a cluster_joiner! Can't use this unblessed combined transcript: ";
                print STDERR " start ". $combined_transcript->seq_region_start. " end ".$combined_transcript->seq_region_end."\n\n";
                #dont use this one
                $combined_transcript = undef;
              } else {
                print STDERR "Unblessed combined transcript is sane and isn't a cluster_joiner.\n";
              }
            } else{
              print STDERR "Unblessed combined transcript at seq_region/chr". $combined_transcript->seq_region_name . " " . 
                            $combined_transcript->seq_region_start . "-" . $combined_transcript->seq_region_end . " did not pass sanity check.\n\n";
              $combined_transcript = undef;
            }
          }
	} # close if (defined $combined_transcript) loop that throws away insane transcripts or cluster joiners


	if(defined $combined_transcript){
	  #leave the loop if $combined_transcript survived sanity and cluster-join check
	  last POS;
	} 
      }
    } # end  POS find match
    
    if (defined $combined_transcript){
      #transfer evidence
      print "Combined transcript checked to be OK. Now transferring evidence onto it.......\n";
      $combined_transcript = $self->_transfer_evidence($combined_transcript, $cdna_match);

      #set biotype
      my $genetype;
      if ( $self->EXTEND_ORIGINAL_BIOTYPE && (length($self->EXTEND_ORIGINAL_BIOTYPE)>0)) {
	my $sep = "";
	if($usingKnown){
	  $sep = "_"  unless ( $self->KNOWN_UTR_GENETYPE =~ m/^_/ );
	  $genetype = $combined_transcript->biotype. $sep. $self->KNOWN_UTR_GENETYPE;
	}
	else{
	  $sep = "_"  unless ( $self->UTR_GENETYPE =~ m/^_/ );
	  $genetype = $combined_transcript->biotype. $sep. $self->UTR_GENETYPE;
	}
      }
      else {
	if($usingKnown){
	  $genetype = $self->KNOWN_UTR_GENETYPE;
	}
	else{
	  $genetype = $combined_genetype;
	}
      }


      #print STDERR "MAKING_GENE FROM ".." AND ".$cdna_match->hit_name."\n";

      print "Combined transcript has supporting features added to it.  Now making combined gene from combined transcript....\n";

      # my $combined_genes = $self->make_gene($genetype, $combined_transcript, $blessed);  # this line of code messed up passing of $blessed!
      my $combined_genes = $self->make_gene($genetype, [$combined_transcript], $blessed);  # make_gene method is expecting an array of transcript refs!

      my $combined_gene;

      #check phases, etc.- if it is not a blessed gene
      if(!$blessed){
	$combined_gene = $self->look_for_both($combined_genes->[0]);
	my $combined_transcript2 = $combined_gene->get_all_Transcripts->[0];
	#$combined_transcript2->sort;

	 if(! (is_Transcript_sane($combined_transcript2)
	       && all_exons_are_valid($combined_transcript2, $self->MAX_EXON_LENGTH, 1)
	       && intron_lengths_all_less_than_maximum($combined_transcript2, $self->MAX_INTRON_LENGTH)
	       && validate_Translation_coords($combined_transcript2, 1)
	       && !contains_internal_stops($combined_transcript2)
	       && $combined_transcript2->translate
	       && has_no_unwanted_evidence($combined_transcript2) )
	   ){

	   #revert to previous version if tests failed
	   $combined_gene = $combined_genes->[0];
	 }
      }
      else{
	#use blessed genes without modifications
	$combined_gene = $combined_genes->[0];
  
      }

      #store as combined
      $self->combined_genes($combined_gene);
      print STDERR "RESULT-CombinedGene: ".$combined_gene->seq_region_name." ".$combined_gene->seq_region_start.
	    "-".$combined_gene->seq_region_end."\n";
    }
    else{
      #retrieving unmerged if no UTR could be defined
      my $unmerged_cds = $self->retrieve_unmerged_gene($cds);
      $self->unmatched_genes($unmerged_cds);
      print STDERR "RESULT-UnmergedGene: ".$unmerged_cds->seq_region_name." ".$unmerged_cds->seq_region_start.
	    "-".$unmerged_cds->seq_region_end."\n";
      next CDS;
    }

  }

  print STDERR "\n### At the end of the matching unblessed (presumably genewise) genes with cDNAs/ESTs, I can add UTR to ".
        (scalar @{$self->combined_genes})." gene(s). ". (scalar @{$self->unmatched_genes})." gene(s) are without UTRs.\n" if !$blessed;
  print STDERR "\n### At the end of the matching blessed genes with cDNAs/ESTs, including the unblessed genes I can add UTR to ".
        (scalar @{$self->combined_genes})." gene(s). ".  (scalar @{$self->unmatched_genes})." gene(s) are without UTRs.\n" if $blessed;


  if(!$blessed && ($self->{'remove_redundant'})){
    #remove redundatant models from the umatched group
    my $unique_unmatched_genes = $self->remove_redundant_models($self->combined_genes, $self->unmatched_genes);

    #replace all unmatched genes
    $self->{'_unmatched_genes'} = [];
    $self->unmatched_genes(@{$unique_unmatched_genes});
    print STDERR "without redundant models, we have ".(scalar @{$self->combined_genes})." combined_genes genes".
      " and ".(scalar @{$self->unmatched_genes})." unmatched_genes\n" if $self->VERBOSE;
  }

}


=head2 calculate_UTR_score

  Arg [1]    : genewise transcript
  Arg [2]    : transcript adding UTR region
  Description: calculate a score favouring transcripts that add UTR on both sides
               Could also favour longer UTRs if desired
  Returntype : int score

=cut

sub calculate_UTR_score {
  my ($self, $gw_gene, $matching_gene) = @_;

  my $BONUS_FOR_BOTH_SIDES = 1;
  my $UTR_score = 0;
  if(($gw_gene->start - $matching_gene->start) && ($matching_gene->end - $gw_gene->end)){
    $UTR_score += $BONUS_FOR_BOTH_SIDES;
  }
  #...to be extended if needed

  return $UTR_score;
}


=head2 score_ditags

  Arg [1]    : array ref with mapped ditags to analyse
  Arg [2]    : transcript to analyse
  Arg [3]    : (optional) int index, where to start looking in the ditag-array
  Description: score ditags in the region of a exon, transcripts or est
               that might support it
               Give a higher score to those features that match well
               to ditags and have a high tag_count.
  Returntype : int score: high value indicates many/well matching ditags
               int index: updated array index
  Exceptions : none

=cut

sub score_ditags{
  my ($self, $ditags, $feature, $index) = @_;

  my $ditag_score    = 0;
  if(!$index){ $index = 0 }
  #check start & end seperately, use best matching
  for(my $i = $index; $i < (scalar @$ditags); $i++){
    my $ditag = $ditags->[$i];


    if(($ditag->{'start'} < $feature->end) && ($ditag->{'end'} > $feature->start)){
      my $start_distance = abs($ditag->{'start'} - $feature->start);
      my $end_distance   = abs($ditag->{'end'}   - $feature->end);

      if(($start_distance < $self->DITAG_WINDOW) && ($end_distance < $self->DITAG_WINDOW)){
        #matching ditag; produce a score favoring those transcripts,
        #that have a high ditag count &/| perfectly positioned ditags.

        #magic equation to generate a score
        my $score = 0;
        $score += ($ditag->{'tag_count'} / (($start_distance +1) * 0.9));
        $score += ($ditag->{'tag_count'} / (($end_distance +1) * 0.9));

        if($score > 0){ $ditag_score += $score; }
      }
    }

    if($ditag->{'end'} < $feature->start){
      $index++;
    }
    elsif($ditag->{'start'} > $feature->start){
      last;
    }
  }
  print STDERR " returning ditag score $ditag_score.\n" if $self->VERBOSE;

  return($ditag_score, $index);
}


=head2 check_for_predefined_pairing

  Arg [1]    : Bio::EnsEMBL::Gene CDS gene
  Arg [2]    : Boolean to indicate blessed gene
  Description: check whether there is a specific cDNA assigned to the given gene model
  Returntype : Bio::EnsEMBL::Gene cDNA gene

=cut

sub check_for_predefined_pairing {
  my ($self, $gene, $blessed) = @_;

  my $cdna = undef;
  my $protein_id;
  my $cdna_id;

  # get hash with cDNAs: $cdna_evidence{$evidence->hseqname()} = $cdna-gene;
  my $cdna_evidence = $self->_cdna_evidence();


  # get the protein id, the gene was built from;
  # do this for blessed genes also, as this should be most reliable
 EXON:
  foreach my $exon(@{$gene->get_all_Exons}){
    my @feat = @{$exon->get_all_supporting_features};
    foreach my $feat(@feat){
      $protein_id = $feat->hseqname;
      last EXON if (defined $protein_id);
    }
  }

  if ((!defined $protein_id || $protein_id eq '') && $blessed){
    #try to use RefSeq xrefs for the blessed genes
    #The CCDS genes are stored with NM-entries in the xref table,
    #so these can be used directly as cDNA-ids.
    #Might have to be adjusted if this changes.
    $protein_id = 'blessed';
    my @xrefs = @{$gene->get_all_DBLinks()};

    print STDERR "Blessed gene has " . (scalar(@xrefs)) . " xrefs of all sorts. Now checking for useful NM_ xrefs...\n" if $self->VERBOSE;

    if(scalar @xrefs){
      foreach my $xref (@xrefs){
	if($xref->display_id =~ "^NM_"){
	  $cdna_id = $xref->display_id;
	  $cdna_id =~ s/\.\S+$//;
	  print STDERR "have xref $cdna_id for blessed gene\n";
	  last;
	}
      }
    }
    else{ print STDERR "no NM_ xrefs for this blessed gene\n"; }
  }
  else{
    if (!defined $protein_id || $protein_id eq ''){
      print STDERR "Couldn't find protein accession for unblessed gene.\n";
      return undef;
    }

    # Using the protein accession of the targetted gene, determine the
    # corresponding cDNA id from the pre-loaded hash.
    # The matching won't work if (i) the protein accession is not from RefSeq NP*
    # but Uniprot or (ii) the NP protein accession doesn't exist in the
    # pre-loaded hash (e.g. the protein has been suppressed by NCBI).
    
    # Protein version info is removed before looking up the pre-loaded hash

    $protein_id =~ s/\.\S+$//;
    $cdna_id = $self->get_cdna_id_from_protein_id($protein_id);
  }

  if (!defined $cdna_id || $cdna_id eq ''){
    print STDERR "No predefined cDNA found for $protein_id.\n";
    return undef;
  }
  else{
    print STDERR "Found predefined cDNA $cdna_id for $protein_id\n";
  }

  #make sure it's not on the kill list
  if(defined ($self->kill_list()->{$cdna_id})){
    print STDERR "Skipping cDNA " . $cdna_id . " as it's present in kill list\n";
    return undef;
  }

  #get the gene for this cDNA
  my $cdna_gene = $cdna_evidence->{$cdna_id};

  #check if they really are overlapping to avoid disappointment when combining
  if($cdna_gene){
    if($self->VERBOSE){ print STDERR "Found cDNA gene object for $cdna_id among the filtered cDNAs.\n" }
    if( !(($gene->seq_region_name  eq $cdna_gene->seq_region_name)
	  && ($gene->strand  == $cdna_gene->strand)
	  && ($gene->seq_region_start  < $cdna_gene->seq_region_end)
	  && ($cdna_gene->seq_region_start < $gene->seq_region_end)
	  && ($gene->get_all_Transcripts->[0]->get_all_Exons->[0]->strand 
	          == $cdna_gene->get_all_Transcripts->[0]->get_all_Exons->[0]->strand)) ){

      print STDERR "cDNA not overlapping unblessed coding model properly, can't use cDNA.";
      $cdna_gene = undef;
    }
  }
  else {
    print STDERR "Unable to fetch cDNA gene object for $cdna_id. Maybe the cDNA alignment is not present in one of the input databases.\n";
    return undef;
  }

  return $cdna_gene;
}


=head2 get_cdna_id_from_protein_id

  Arg [1]    : String protein-id
  Description: Find the cDNA-id a given protein-id is linked to as a "known" evidence
  Returntype : String cDNA-id

=cut

sub get_cdna_id_from_protein_id {
  my ($self, $protein_id) = @_;

  my $cdna_id = undef;
  if(defined($self->_known_pairs()->{$protein_id})){
    $cdna_id = $self->_known_pairs()->{$protein_id};
  }

  return $cdna_id;
}


=head2 create_predefined_pairing

  Arg [1]    : none
  Description: Read GenBank sequence file linking NP proteins to a specific NM cDNA
               to use as "known" UTR evidence 
               Source : ftp://ftp.ncbi.nlm.nih.gov/refseq/<SPECIES>/mRNA_Prot/<species>.protein.gpff.gz
  Returns    : nothing

=cut

sub create_predefined_pairing {
  my ($self, $known_utr_file) = @_;

  print STDERR "\nParsing GenBank file $known_utr_file for KnowUTR pairing.\n" if $self->VERBOSE;

  open(REFSEQ, "<$known_utr_file") or die "Can't open ".$known_utr_file.": $@\n";

  my $cdna_id;
  my $protein_id;
  my %predefined_pairing = ();

  while(<REFSEQ>){
    #use only these 2 fields of the file:
    next unless /^VERSION|DBSOURCE/;

    if(/VERSION/){

      next unless /(NP\S+)/;

      if(defined $protein_id){
	die("previous protein_id [$protein_id] has not been cleared out\n");
      }
      if(defined $cdna_id){
	die("previous cdna_id [$cdna_id] has not been cleared out ...\n");
      }
      $protein_id = $1;
    }

    if(/DBSOURCE/){
      # don't want NCs or NGs, etc.
      if(/(NC\_\S+)|(NG\_\S+)|(XM\S+)|(AC\_\S+)/){
	$cdna_id = undef;
	$protein_id = undef;
	next;
      }
      if(!defined $protein_id){
	die("something very wrong - no protein_id for $_\n");
      }
      if (defined $cdna_id){
	die("previous cdna_id [$cdna_id] has not been cleared out ...\n");
      }

      next unless /(NM\S+)/;

      $cdna_id = $1;
      #remove version info
      $cdna_id =~ s/\.\S+$//;
      $protein_id =~ s/\.\S+$//;
      if(exists($predefined_pairing{$protein_id})){ print STDERR $protein_id." allready exists!\n"; }
      $predefined_pairing{$protein_id} = $cdna_id;

      $cdna_id    = undef;
      $protein_id = undef;
    }
  }

  close REFSEQ;
  $self->_known_pairs(\%predefined_pairing);
}


=head2 find_cluster_joiners

  Arg [1]    : ref to array of GeneClusters
  Description: screens UTR modified transcripts for those that
               potentially misjoin genewise clusters, and remove them
  Returntype : ref to array of Bio::EnsEMBL::Gene 

=cut

sub find_cluster_joiners{
  my ($self, $transcript) = @_;

  my @clusters;
  my $transcript_start;
  my $transcript_end;
  my $overlaps_previous_cluster = 0;
  my $matching_clusters = 0;

  if ($transcript->get_all_Exons->[0]->strand == 1){
    @clusters = @{$self->forward_genewise_clusters};
    $transcript_start = $transcript->start_Exon->start;
    $transcript_end   = $transcript->end_Exon->end;
  }
  else{
    @clusters = @{$self->reverse_genewise_clusters};
    $transcript_start = $transcript->end_Exon->start;
    $transcript_end = $transcript->start_Exon->end;
  }

  print STDERR "Checking transcript for cluster joining: seq_region/chr ". $transcript->seq_region_name . ", ".
        $transcript->seq_region_start. "-". $transcript->seq_region_end. ".\n" if $self->VERBOSE;
 
  if(!scalar(@clusters)){
    print STDERR "Odd, not even a single cluster of combined_transcript!? Should be at least one even if there are no cluster-joiners.\n" if $self->VERBOSE;
    return 0;
  }

  #NEW CODE: CLUSTERING WITH EXON COORDINATES
  foreach my $cluster(@clusters){
    if( _overlapping_genes($transcript, $cluster) ){
      $matching_clusters++;
    }
    if($matching_clusters>1){
      print STDERR "transcript joins clusters - you should discard it\n" if $self->VERBOSE;
      return 1;
    }
  }

  return 0;
}


=head2 _overlapping_genes

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster
  Description: Check whether any of the exons of the transcript
               overlap any of the exons of the cluster.
  Returntype : boolean, true if overlap exists

=cut

sub _overlapping_genes {
  my ($gene1, $cluster) = @_;

  # quit if genes do not have genomic overlap
  if ($gene1->end < $cluster->start || $gene1->start > $cluster->end) {
    return 0;
  }
  # overlap check based on all (noncoding + coding) Exons
  foreach my $exon1 (@{$gene1->get_all_Exons}){
    foreach my $exon2 (@{$cluster->get_all_Exons}){
      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	return 1;
      }
    }
  }

  return 0;
}


=head2 _filter_cdnas

  Arg [1]    : ref to array of Bio::EnsEMBL::Gene
  Arg [2]    : flag whether these genes are EST genes rather than cDNA genes
  Description: This method checks the cDNA and EST genes
               translation is not checked as these genescome without
               translation
  Returntype : ref to array of filtered Bio::EnsEMBL::Gene

=cut

sub _filter_cdnas{
  my ($self, $cdna_arrayref, $ests) = @_;
  my @newcdna;
  my %cdna_evidence;

 cDNA_GENE:
  foreach my $cdna (@{$cdna_arrayref}) {

  cDNA_TRANSCRIPT:
    foreach my $tran (@{$cdna->get_all_Transcripts}) {
      #$tran->sort;

      if(!$ests){
	#store cDNA gene evidence for later
	foreach my $evidence (@{ $tran->get_all_supporting_features }){
	  #remove version number?!
	  my $evidence_name = $evidence->hseqname();
	  $evidence_name =~ s/\.\S+$//;
	  #print STDERR "evidence: ".$evidence->hseqname()." => ".$evidence_name."\n";
	  $cdna_evidence{$evidence_name} = $cdna;
	}
      }

      # rejecting on basis of intron length may not be valid here
      # - it may not be that simple in the same way as it isn;t that simple in Targetted & Similarity builds
      next cDNA_TRANSCRIPT unless ( is_Transcript_sane($tran)
				    && all_exons_are_valid($tran, $self->MAX_EXON_LENGTH, 1)
				    && intron_lengths_all_less_than_maximum($tran, $self->MAX_INTRON_LENGTH)
				    && has_no_unwanted_evidence($tran) );
      push(@newcdna, $cdna);
    }
  }

  if(!$ests){
    #store cDNAs for later use
    $self->_cdna_evidence(\%cdna_evidence);
  }

  return \@newcdna;
}


=head2  write_output

  Arg [1]    : none
  Description: Do some cleaning up uf the result and store the genes to the database
  Returntype : none

=cut

sub write_output {
  my($self) = @_;

  # write genes in the output database
  my $output_db = get_db_adaptor_by_string($self->OUTPUT_DB,1);
  my $gene_adaptor = $output_db->get_GeneAdaptor;

  print STDERR "Have ".scalar (@{$self->output})."(".$totalgenes.") genes to write\n";

 GENE:
  foreach my $gene (@{$self->output}) {

    if(!$gene->analysis ||
       $gene->analysis->logic_name ne $self->analysis->logic_name){
      $gene->analysis($self->analysis);
    }

    # double check gene coordinates
    $gene->recalculate_coordinates;

    #As all genes/exons are stored as new in the target db
    #it's save to remove the old adaptor & old dbID here,
    #to avoid warnings from the store function.
    foreach my $exon (@{$gene->get_all_Exons}) {
      $exon->adaptor(undef);
      $exon->dbID(undef);
    }
    $gene->dbID(undef) ;
    $gene->adaptor(undef);
    for ( @{$gene->get_all_Transcripts} ) {
      for ( @{$_->get_all_supporting_features}){
        $_->dbID(undef);
        $_->adaptor(undef);
      }
      $_->dbID(undef);
      $_->adaptor(undef);
    }

    eval {
      $gene_adaptor->store($gene); 
      print STDERR "wrote gene dbID " . $gene->dbID  . "\t".$gene->biotype . "\n"  ; 
    };
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE " . $gene->dbID. "of type " . $gene->type  . "\n\n$@\n\nSkipping this gene\n";
    }
  }

}


=head2  convert_to_extended_genes

  Arg [1]    : ref to array of Bio::EnsEMBL::Gene
  Description: converts transcripts to Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended
  Returntype : ref to array of TranscriptExtended
  Exceptions : throws if gene has multiple transcripts

=cut

sub convert_to_extended_genes {
  my ($self, $genes) =@_;

  my @new_genes;
  for my $g (@$genes){ 
    #my $st  = _get_evidence_set( $g->biotype );
    #my $set = $g->biotype;
    my ( $set ) = $self->get_evidence_set ( $g->biotype ) ;
    for my $t (@{$g->get_all_Transcripts}){
      throw ("gene has more than one trancript - only processing 1-gene-1-transcript-genes")
	if (@{$g->get_all_Transcripts}>1);

      # conversion 
      my $gene_from_pt = Bio::EnsEMBL::Gene->new(
                         -start    => $t->start ,
                         -end      => $t->end ,
                         -strand   => $t->strand ,
                         -slice    => $t->slice ,
                         -biotype  => $g->biotype,
                         -analysis => $t->analysis,
                        );

      my $new_tr = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended->new(
                    -BIOTYPE  => $g->biotype ,
                    -ANALYSIS => $t->analysis ,
                   );
      $new_tr->ev_set($set);

      $new_tr->add_supporting_features(@{$t->get_all_supporting_features});

      my @pt_exons  = @{$t->get_all_Exons} ;
      for (my $i=0 ; $i<scalar(@pt_exons) ; $i++) { 
	# converting Bio::EnsEMBL::PredictionExon into ExonExtened (ISA Bio::EnsEMBL::Exon)
	my $pte = $pt_exons[$i];
	bless $pte,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended";
	$pte->biotype($g->biotype);
	$pte->ev_set($set);
	$pte->end_phase(0);
	$pte->phase(0);
	$pte->next_exon($pt_exons[$i+1]) ;
        if ($i == 0 ) {
	  $pte->prev_exon(0);
	} else {
	  $pte->prev_exon($pt_exons[$i-1]);
	}
	$pte->transcript($new_tr) ; 
	$pte->analysis($t->analysis) ; 
      };

      # Extending the Bio::EnsEMBL::Transcript object by ev_set methods
      for (@pt_exons) {
	$new_tr->add_Exon($_);
      }
      #$new_tr->sort;
      $gene_from_pt->add_Transcript($new_tr);
      push @new_genes , $gene_from_pt;
    }
  }

  return \@new_genes;
}


=head2 make_gene

  Arg [1]    : string representing genetype to be associated with genes
  Arg [2]    : an array of Bio::EnsEMBL::Transcript
  Description: Constructs Bio::EnsEMBL::Gene objects from UTR-modified transcripts
  Returntype : none; new genes are stored in self->combined_genes
  Exceptions : throws when missing genetype or analysis

=cut

sub make_gene{
  my ($self, $genetype, $transcripts, $blessed) = @_;
  
  unless ( $genetype ){
    throw("You must define UTR_GENETYPE in Bio::EnsEMBL::Analysis::Conf::GeneBuild::UTR_Builder");
  }

  # an analysis should be passed in via the RunnableDB.m parent class:
  my $analysis = $self->analysis;
  unless ($analysis){
    throw("You have to pass an analysis to this RunnableDB through new()");
  }

  my @genes;
  my $count=0;

  foreach my $trans(@$transcripts){
    #$trans->sort;
    if (!$blessed) {
      unless ( is_Transcript_sane($trans)
  	     && intron_lengths_all_less_than_maximum($trans, $self->MAX_INTRON_LENGTH)
	     && all_exons_are_valid($trans, $self->MAX_EXON_LENGTH, 1)
	     && has_no_unwanted_evidence($trans) ){
        print STDERR "\nUnblessed transcript rejected when trying to make combined_gene from combined_transcript.\n";
        return;
      }
    } elsif ($blessed) {
      unless ( is_Transcript_sane($trans)
             && has_no_unwanted_evidence($trans) ){
        print STDERR "\nBlessed transcript rejected when trying to make combined_gene from combined_transcript.\n";
        return;
      }
    }
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->biotype($genetype);
    $trans->biotype($genetype);
    $trans->analysis($analysis);
    $gene->add_Transcript($trans);
    $gene->analysis($analysis);

    # do not modify the analysis of the supporting features
    # they should be the original ones: cdna, targetted_genewise or similarity_genewise

    if($self->validate_gene($gene)){
      push (@genes,$gene);
      $count++;
    }
    else{
      print STDERR "\nCould not validate gene!\n";
    }
  }

  return(\@genes);
}


=head2 match_protein_to_cdna

  Arg [1]    : Bio::EnsEMBL::Gene
  Arg [2]    : int is_est
  Description: this method tries to find the cdnas that can be merged with the genewise genes.
             Basically it checks for exact internal splice matching in the 5' and 3' exons of the genewise gene.
             In order to match a starting genewise exon with an cDNA exon, we need to have 
             a. exactly coinciding exon boundaries
             b. either cdna exon has start <= genewise exon start, 
             OR cdna exon start lies within $exon_slop bp of genewise exon start AND 
             the cdna transcript will add exra UTR exons. 
             Substitute "end" for "start" for 3prime ends of transcripts  
             BUT do not allow through any cDNA that will result just in a 
             shortened peptide without additional UTR exons.
  Returntype : ref to array of Bio::EnsEMBL::Gene, ref to hash relating Bio::EnsEMBL::Gene to UTR length
  Exceptions : 

=cut

sub match_protein_to_cdna {
  my ( $self, $gw, $is_est ) = @_;

  print STDERR "Looking at ESTs.\n" if $is_est;
  my ( %UTR_hash, %UTR_side_indicator_hash );
  my @other_genes;

  my @matching_e2g;
  my @gw_tran = @{ $gw->get_all_Transcripts };

  my @gw_exons = @{ $gw_tran[0]->get_all_Exons };
  my $strand   = $gw_exons[0]->strand;
  if ( $gw_exons[$#gw_exons]->strand != $strand ) {
    warn( "first and last gw exons have different strands " .
          "- eventually can't make a sensible combined gene (gw dbID " .
          $gw_tran[0]->dbId . ")" );
    return undef;
  }
  if (@gw_exons) {
    if ( $strand == 1 ) {
      @gw_exons = sort { $a->start <=> $b->start } @gw_exons;
    }
    else {
      @gw_exons = sort { $b->start <=> $a->start } @gw_exons;
    }
  }
  else {
    warn( "gw gene without exons: " . $gw->dbID . ", skipping it" );
    return undef;
  }

  my $exon_slop = 20;

  my $cds_length = length( $gw_tran[0]->translateable_seq() );

  my @genes;
  if ($is_est) {
    @genes = @{ $self->ests };
  }
  else {
    @genes = @{ $self->cdna_genes };
  }

cDNA:
  foreach my $e2g (@genes) {

    my @egtran = @{ $e2g->get_all_Transcripts };

    # if the cDNA or EST ($egtran[0]) is on the kill list, we skip it

    my $egtran_hitname =
      ${ $egtran[0]->get_all_supporting_features }[0]->hseqname;

    $egtran_hitname =~ s/\.\d//;

    if ( defined( $self->kill_list()->{$egtran_hitname} ) ) {
      print STDERR "Skipping cDNA " . $egtran_hitname .
        " as it's present in kill list\n";
      next cDNA;
    }

    my @eg_exons = @{ $egtran[0]->get_all_Exons };

    $strand = $eg_exons[0]->strand;
    if ( $eg_exons[$#eg_exons]->strand != $strand ) {
      warn( "first and last e2g exons have different strands" .
            " - skipping e2g transcript (dbID " . $egtran[0]->dbID .
            ")" );
      next cDNA;
    }
    if ( $strand == 1 ) {
      @eg_exons = sort { $a->start <=> $b->start } @eg_exons;
    }
    else {
      @eg_exons = sort { $b->start <=> $a->start } @eg_exons;
    }

    my $fiveprime_match  = 0;
    my $threeprime_match = 0;
    my $left_exon;
    my $right_exon;
    my $left_diff  = 0;
    my $right_diff = 0;

    # Lets deal with single exon genes first
    if ( $#gw_exons == 0 ) {
      foreach my $current_exon (@eg_exons) {
        if ( $current_exon->strand != $gw_exons[0]->strand ) {
          next cDNA;
        }

        # don't yet deal with genewise leakage for single exon genes
        if ( $gw_exons[0]->end <= $current_exon->end &&
             $gw_exons[0]->start >= $current_exon->start )
        {
          $fiveprime_match  = 1;
          $threeprime_match = 1;

          $left_exon  = $current_exon;
          $right_exon = $current_exon;
          $left_diff  = $gw_exons[0]->start - $current_exon->start;
          $right_diff = $current_exon->end - $gw_exons[0]->end;
        }
      }

    }
    else {    # Now the multi exon genewises

    cDNA_EXONS:
      foreach my $current_exon (@eg_exons) {
        if ( $current_exon->strand != $gw_exons[0]->strand ) {
          next cDNA;
        }

        if ( $gw_exons[0]->strand == 1 ) {

          #FORWARD:

          # 5prime
          if (
            $gw_exons[0]->end == $current_exon->end &&
            # either e2g exon starts before genewise exon
            ( $current_exon->start <= $gw_exons[0]->start ||
              # or e2g exon is a bit shorter but there
              # are spliced UTR exons as well
              ( abs( $current_exon->start - $gw_exons[0]->start ) <=
                $exon_slop && $current_exon != $eg_exons[0] ) ) )
          {
            $fiveprime_match = 1;
            $left_exon       = $current_exon;
            $left_diff       = $gw_exons[0]->start - $current_exon->start;
          }

          # 3prime
          elsif (
            $gw_exons[$#gw_exons]->start == $current_exon->start &&
            # either e2g exon ends after genewise exon
            ( $current_exon->end >= $gw_exons[$#gw_exons]->end ||
              # or there are UTR exons to be added
              ( abs( $current_exon->end - $gw_exons[$#gw_exons]->end )
                <= $exon_slop && $current_exon != $eg_exons[$#eg_exons]
              ) ) )
          {
            $threeprime_match = 1;
            $right_exon       = $current_exon;
            $right_diff       = $current_exon->end - $gw_exons[0]->end;
          }
        } ## end if ( $gw_exons[0]->strand...)

        elsif ( $gw_exons[0]->strand == -1 ) {

          #REVERSE:

          # 5prime
          if (
            $gw_exons[0]->start == $current_exon->start &&
            # either e2g exon ends after gw exon
            ( $current_exon->end >= $gw_exons[0]->end ||
              # or there are UTR exons to be added
              ( abs( $current_exon->end - $gw_exons[0]->end ) <=
                $exon_slop && $current_exon != $eg_exons[0] ) ) )
          {
            $fiveprime_match = 1;
            $right_exon      = $current_exon;
            $right_diff      = $current_exon->end - $gw_exons[0]->end;
          }

          #3prime
          elsif (
            $gw_exons[$#gw_exons]->end == $current_exon->end &&
            # either e2g exon starts before gw exon
            ( $current_exon->start <= $gw_exons[$#gw_exons]->start ||
              # or there are UTR exons to be added
              ( abs( $current_exon->start - $gw_exons[$#gw_exons]->start
                ) <= $exon_slop &&
                $current_exon != $eg_exons[$#eg_exons] ) ) )
          {
            $threeprime_match = 1;
            $left_exon        = $current_exon;
            $left_diff = $gw_exons[0]->start - $current_exon->start;
          }
        } ## end elsif ( $gw_exons[0]->strand...)
      } ## end cDNA_EXONS: foreach my $current_exon (@eg_exons)
    } ## end else [ if ( $#gw_exons == 0 )]

    # can match either end, or both
    if ( $fiveprime_match || $threeprime_match ) {
      my ( $UTR_length, $left_UTR_length, $right_UTR_length ) =
        $self->_compute_UTRlength( $egtran[0], $left_exon,
                                   $left_diff, $right_exon,
                                   $right_diff );

      my $UTR_diff = $egtran[0]->length;

      #make sure CDS is not much smaller than UTR
      if ( ( $cds_length*10 ) > $UTR_diff ) {
        $UTR_hash{$e2g}                = $UTR_length;
        $UTR_side_indicator_hash{$e2g} = 1;

        if ( $self->VERBOSE ) {
          print STDERR "considering cDNA " . $e2g->seq_region_start .
            "-" . $e2g->seq_region_end . "[cds_len*10 " .
            ( $cds_length*10 ) . " vs UTR_diff " . $UTR_diff . "]\n";
        }

        push( @matching_e2g, $e2g );
      }
      else {
        if ( $self->VERBOSE ) {
          print STDERR "didnt pass UTR length check [cds_len*10" .
            ( $cds_length*10 ) . " vs UTR_diff " . $UTR_diff . "]: " .
            $e2g->seq_region_start . "-" . $e2g->seq_region_end . "\n";
        }
      }
    } ## end if ( $fiveprime_match ...)

  } ## end cDNA: foreach my $e2g (@genes)

  print STDERR "\nmatch_protein_to_cdna is returning " .
    ( scalar @matching_e2g ) . " UTR evidence candidates.\n";

  return ( \@matching_e2g, \%UTR_hash, \%UTR_side_indicator_hash );
} ## end sub match_protein_to_cdna


=head2 _compute_UTRlength

  Arg [1]    : eg-transcript
  Arg [2]    : first matching exon
  Arg [3]    : (int) UTR bases of first matching exon
  Arg [4]    : last matching exon
  Arg [5]    : (int) UTR bases of first matching exon
  Description: add up genomic extend of UTR regions
  Returntype : int (basepairs) total UTR-length, left UTR-length, right UTR-lenght
  Exceptions : none

=cut

sub _compute_UTRlength{
 my ($self, $transcript, $left_exon, $left_diff, $right_exon, $right_diff) = @_;

 my $strand = $transcript->start_Exon->strand;
 my @exons  = sort { $a->start <=> $b->start } @{ $transcript->get_all_Exons };

 my $UTRlength  = 0;
 my $in_UTR     = 1;
 my $start_flag = 0;
 my $left_UTR   = $left_diff;
 my $right_UTR  = $right_diff;

 foreach my $exon ( @exons ){
   if ( defined $left_exon && $exon == $left_exon ){
     $UTRlength += $left_diff;
     $in_UTR     = 0;
   }
   elsif( defined $right_exon && $exon == $right_exon ){
     $UTRlength += $right_diff;
     $in_UTR     = 1;
   }
   elsif( $in_UTR == 1 ){
     $UTRlength += $exon->length;
     if(!$start_flag){
       $left_UTR  += $UTRlength;
       $start_flag = 1;
     }
     else{
       $right_UTR  += $UTRlength;
     }
   }
 }

 return($UTRlength, $left_UTR, $right_UTR);
}


=head2 _merge_genes

  Arg [1]    : ref to array of Bio::EnsEMBL::Gene
  Description: merges adjacent exons if they are frameshifted; stores
               component exons as subSeqFeatures of the merged exon
  Returntype : ref to array of Bio::EnsEMBL::Gene

=cut

sub _merge_genes {
  my ($self, $genesref, $blessed) = @_;
  my @merged;
  my $count = 1;

 UNMERGED_GENE:
  foreach my $unmerged (@{$genesref}){
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->dbID($unmerged->dbID);
    foreach my $dblink (@{$unmerged->get_all_DBLinks()}){
    # print STDERR "Adding xref ". $dblink->display_d . " from unmerged to merged gene......\n";
    $gene->add_DBEntry($dblink);
    }

    my @pred_exons;
    my $ecount = 0;

    # order is crucial
    
    my @trans = @{$unmerged->get_all_Transcripts};
    if(scalar(@trans) != 1) { 
      throw("Gene with dbID " . $unmerged->dbID .
		   " has NO or more than one related transcript, where a 1-gene-to-1-transcript-relation is assumed.".
		   " Check preceding analysis \n");
    }
    #$trans[0]->sort;
    # check the sanity of the transcript only if it's not blessed. (We don't care about long introns etc in blessed models.)
    if (!$blessed) {
      if( !is_Transcript_sane($trans[0])
             || !all_exons_are_valid($trans[0], $self->MAX_EXON_LENGTH, 1)
             || !has_no_unwanted_evidence($trans[0]) 
             || !intron_lengths_all_less_than_maximum($trans[0], $self->MAX_INTRON_LENGTH) ){
        print STDERR "Transcript of gene (dbID ". $unmerged->dbID .") at seq_region/chr ". $unmerged->seq_region_name.", ".
             $unmerged->seq_region_start." to ".$unmerged->seq_region_end ." did NOT pass sanity check! (See warning above)\n";
        next UNMERGED_GENE;
      }
    }
    ### we follow here 5' -> 3' orientation ###

    my $cloned_translation = new Bio::EnsEMBL::Translation;

    my @unmerged_exons = @{$trans[0]->get_all_Exons};
    my $strand   = $unmerged_exons[0]->strand; 
    my $previous_exon;

  EXON:
    foreach my $exon(@unmerged_exons){
	
	## genewise frameshift? we merge here two exons separated by max 10 bases into a single exon
	#if ($ecount && $pred_exons[$ecount-1]){
	#  $previous_exon = $pred_exons[$ecount-1];
	#}
	
	$ecount++;
	
	my $separation = 0;
	my $merge_it   = 0;
	
	## we put every exon following a frameshift into the first exon before the frameshift
	## following the ordering 5' -> 3'
	if (defined($previous_exon)){
	
	  #print STDERR "previous exon: ".$previous_exon->start."-".$previous_exon->end."\n";
	  #print STDERR "current exon : ".$exon->start."-".$exon->end."\n";
	  if ($strand == 1){
	    $separation = $exon->start - $previous_exon->end - 1;
	  }
	  elsif( $strand == -1 ){
	    $separation = $previous_exon->end - $exon->start - 1;
	  }
	  if ($separation <=10){
	    $merge_it = 1;
	  }	
	}
	
    # NEVER MERGE EXONS WITH FRAMESHIFT AS IT GENERATES SUPPORTING FEATURES BEYOND THE END OF THE EXON
    $merge_it = 0;

	if ( defined($previous_exon) && $merge_it == 1){
	  # combine the two
	
	  # the first exon (5'->3' orientation always) is the containing exon,
	  # which gets expanded and the other exons are added into it
      
          # print STDERR "Merging exon of length " . $exon->length . " with exon of length " . $previous_exon->length ."\n";

          if ($strand == 1) {
            $previous_exon->end($exon->end);
          } else {
            $previous_exon->start($exon->start);
          }

          # print STDERR "Addding exon " . $exon->seq_region_start . " " . $exon->seq_region_end . " to previous exon ".
          #      $previous_exon->seq_region_start . " ". $previous_exon->seq_region_end. " as sub_SeqFeature\n"; 

	  $previous_exon->add_sub_SeqFeature($exon,'');
	
	  # if this is end of translation, keep that inf:
	  if ( $exon == $trans[0]->translation->end_Exon ){
	    $cloned_translation->end_Exon( $previous_exon );
	    $cloned_translation->end($trans[0]->translation->end);
	  }
	
	  my %evidence_hash;
	  foreach my $sf( @{$exon->get_all_supporting_features}){
	    if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	      next;
	    }
	    #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
	    $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
	    $previous_exon->add_supporting_features($sf);
	  }
	  next EXON;
	}
	else{
          # This is the majority of cases (no frameshift exons to merge), avoid the print statement below

	  # make a new Exon - clone using Bio::EnsEMBL::Analysis::Tools::GeneB uildUtils::ExonUtils
	  my $cloned_exon = clone_Exon($exon);
   
          # print "No fs exons to merge. Adding exon " . $exon->seq_region_start . " " . $exon->seq_region_end . " as sub_SeqFeature to cloned exon ".
          #      $cloned_exon->seq_region_start . " ". $cloned_exon->seq_region_end. "\n"; 

	  $cloned_exon->add_sub_SeqFeature($exon,'');
	
	  # if this is start/end of translation, keep that info:
	  if ( $exon == $trans[0]->translation->start_Exon ){
	    $cloned_translation->start_Exon( $cloned_exon );
	    $cloned_translation->start($trans[0]->translation->start);
	  }
	  if ( $exon == $trans[0]->translation->end_Exon ){
	    $cloned_translation->end_Exon( $cloned_exon );
	    $cloned_translation->end($trans[0]->translation->end);
	  }
	
	  $previous_exon = $cloned_exon;
	  push(@pred_exons, $cloned_exon);
	}
      }

    # transcript
    my $merged_transcript   = new Bio::EnsEMBL::Transcript;
    $merged_transcript->dbID($trans[0]->dbID);

    foreach my $pe(@pred_exons){
	$merged_transcript->add_Exon($pe);
    }

    #$merged_transcript->sort;
    $merged_transcript->translation($cloned_translation);

    my @seqeds = @{$trans[0]->translation->get_all_SeqEdits};
    if (scalar(@seqeds)) {
      print "Copying sequence edits\n" if $self->VERBOSE; 
      foreach my $se (@seqeds) {
	my $new_se =
	  Bio::EnsEMBL::SeqEdit->new(
				     -CODE    => $se->code,
				     -NAME    => $se->name,
				     -DESC    => $se->description,
				     -START   => $se->start,
				     -END     => $se->end,
				     -ALT_SEQ => $se->alt_seq
				    );
	my $attribute = $new_se->get_Attribute();
	$cloned_translation->add_Attributes($attribute);
      }
    }
    my @support = @{$trans[0]->get_all_supporting_features};
    if (scalar(@support)) {
      $merged_transcript->add_supporting_features(@support);
    }
    $merged_transcript->add_Attributes(@{$trans[0]->get_all_Attributes});

    # and gene
    $gene->add_Transcript($merged_transcript);
    $gene->biotype($unmerged->biotype);
    # my @merged_gene_xrefs = @{$gene->get_all_DBLinks()};     # Check if fs-merged genes (esp blessed ones) still retain xrefs or not
    # print "Frameshift-merged gene has ".scalar(@merged_gene_xrefs). " xrefs.\n";
    push(@merged, $gene);
    $count++;

    # store match between merged and original gene so we can easily retrieve the latter if we need to
    $self->merged_unmerged_pairs($gene,$unmerged);

  } # end UNMERGED_GENE

  return @merged;
}


=head2 combine_genes

  Arg [1]    : genewise gene
  Arg [2]    : cDNA gene
  Description: combine gene with matching cDNA
  Returntype : Bio::EnsEMBL::Transcript

=cut

sub combine_genes{
  my ($self, $gw, $e2g, $blessed) = @_; 
  my $modified_peptide = 0;
  my @combined_transcripts = ();

  # should be only 1 transcript
  my @gw_tran   = @{$gw->get_all_Transcripts};
  #$gw_tran[0]->sort;

  my @gw_exons  = @{$gw_tran[0]->get_all_Exons}; # ordered array of exons

  my @egtran    = @{$e2g->get_all_Transcripts};
  #$egtran[0]->sort;
  my @e2g_exons = @{$egtran[0]->get_all_Exons}; # ordered array of exons

  # clone transcript
  my $newtranscript = new Bio::EnsEMBL::Transcript;

  if ( $self->EXTEND_ORIGINAL_BIOTYPE && (length($self->EXTEND_ORIGINAL_BIOTYPE)>0) ) {
    # original biotype of gene / transcript will only be extended, not overwritten
    $newtranscript->biotype($gw->biotype) ;
  }
  foreach my $exon(@gw_exons){
    my $new_exon = clone_Exon($exon);
    if($exon->sub_SeqFeature && scalar($exon->sub_SeqFeature) > 1 ){
      my @sf    = sort {$a->start <=> $b->start} $exon->sub_SeqFeature;
      foreach my $sf (@sf) {
        my $cloned_sf = clone_Exon($sf);
        # print STDERR "Adding subsf ". $cloned_sf->seq_region_start . "-" . $cloned_sf->seq_region_end . " to cloned Exon of length " . 
        #       $exon->length . " (" . $exon->seq_region_start . "-" . $exon->seq_region_end. ") before cds and UTR evid are combined.\n";
        $new_exon->add_sub_SeqFeature($cloned_sf,'');
      }
    }
    $newtranscript->add_Exon($new_exon);
  }

  my @support = @{$gw_tran[0]->get_all_supporting_features};
  if (scalar(@support)) {
    $newtranscript->add_supporting_features(@support);
  }
  $newtranscript->add_Attributes(@{$gw_tran[0]->get_all_Attributes});

  my $translation = new Bio::EnsEMBL::Translation;

  $translation->start($gw_tran[0]->translation->start);
  $translation->end($gw_tran[0]->translation->end);
  $translation->start_Exon($gw_tran[0]->translation->start_Exon);
  $translation->end_Exon($gw_tran[0]->translation->end_Exon);

  my @seqeds = @{$gw_tran[0]->translation->get_all_SeqEdits};
  if (scalar(@seqeds)) {
    print "Copying sequence edits\n" if $self->VERBOSE;
    foreach my $se (@seqeds) {
      my $new_se =
              Bio::EnsEMBL::SeqEdit->new(
                -CODE    => $se->code,
                -NAME    => $se->name,
                -DESC    => $se->description,
                -START   => $se->start,
                -END     => $se->end,
                -ALT_SEQ => $se->alt_seq
              );
      my $attribute = $new_se->get_Attribute();
      $translation->add_Attributes($attribute);
    }
  }

  $newtranscript->translation($translation);
  $newtranscript->translation->start_Exon($newtranscript->start_Exon);
  $newtranscript->translation->end_Exon($newtranscript->end_Exon);

  my $eecount = 0;
  my $modified_peptide_flag = 0;

 EACH_E2G_EXON:
  foreach my $ee (@e2g_exons){

      # check strands are consistent
      if ($ee->strand != $gw_exons[0]->strand){
	  warn("coding (gw or blessed) and e2g exons have different strands - can't combine genes\n") ;
	  return undef;
      }

      # single exon genewise prediction?
      if(scalar(@gw_exons) == 1) {
	  ($newtranscript, $modified_peptide_flag) = $self->transcript_from_single_exon_genewise( $ee,
												  $gw_exons[0],
												  $newtranscript,
												  $eecount,
												  @e2g_exons);
      }

      else {
	  ($newtranscript, $modified_peptide_flag) = $self->transcript_from_multi_exon_genewise($ee,
												$newtranscript,
												$translation,
												$eecount,
												$gw,
												$e2g)
	  }


      if ( $modified_peptide_flag ){
	$modified_peptide = 1;
      }

      # increment the exon
      $eecount++;

  } # end of EACH_E2G_EXON

  print STDERR "Combined transcript just created in combine_genes method has ".  scalar@{$newtranscript->get_all_Exons()}  . " exon(s), translation length "
               . $newtranscript->translation->length()."\n";
               

  #don't modify translation of blessed genes
  my $biotype = $gw->biotype;
  if(defined($newtranscript) && $modified_peptide && (($self->{'blessed_type'}) =~ m/$biotype/)){
    print STDERR "translation of blessed gene would need to be modified - not using combined gene.\n";
    return undef;
  }

  ##############################
  # expand merged exons
  ##############################

  # the new transcript is made from a merged genewise gene
  # check the transcript and expand frameshifts in all but original 3' gw_exon
  # (the sub_SeqFeatures have been flushed for this exon)
  if (defined($newtranscript)){
    #$newtranscript->sort;

    foreach my $ex (@{$newtranscript->get_all_Exons}){

      # print STDERR "Looking at exon of length " . $ex->length . "\n";
      if($ex->sub_SeqFeature && scalar($ex->sub_SeqFeature) > 1 ){
	my @sf    = sort {$a->start <=> $b->start} $ex->sub_SeqFeature;
        print STDERR "Unpacking a merged exon of length ". $ex->length . " at seq_region/chr ". $ex->seq_region_name . 
                     " (" . $ex->seq_region_start . "-" . $ex->seq_region_end .")\n" if $self->VERBOSE;
 
	my $first = shift(@sf);

	$ex->end($first->end);

	# add back the remaining component exons
	foreach my $s(@sf){
	  $newtranscript->add_Exon($s);
	  #$newtranscript->sort;
	}
	# flush the sub_SeqFeatures
	$ex->flush_sub_SeqFeature;
      }
    }

    # check that the resulting transcript:

    if (!$blessed) {
      unless( is_Transcript_sane($newtranscript)
  	    && all_exons_are_valid($newtranscript, $self->MAX_EXON_LENGTH, 1)
	    && intron_lengths_all_less_than_maximum($newtranscript, $self->MAX_INTRON_LENGTH)
	    && has_no_unwanted_evidence($newtranscript) ){
        print STDERR "problems with this unblessed combined_transcript inside combine_genes method (no predefined cDNA), return undef\n";
        return undef;
      }
    } elsif ($blessed) {
      unless( is_Transcript_sane($newtranscript)
            && has_no_unwanted_evidence($newtranscript) ){
        print STDERR "problems with this blessed combined_transcript inside combine_genes method (no predefined cDNA), return undef\n";
        return undef;
      }   
    }

    #check the translation
    unless( validate_Translation_coords($newtranscript, 1)
	    && !contains_internal_stops($newtranscript)
	    && $newtranscript->translate){
      print STDERR "problems with this translation from the combined_transcript inside combine_genes method, return undef\n";
      return undef;
    }

    # check translation is the same as for the genewise gene we built from
    my $foundtrans = 0;

    # the genewise translation can be modified due to a disagreement in a
    # splice site with cdnas. This can happen as neither blast nor genewise can
    # always find very tiny exons.
    # we then recalculate the translation using genomewise:

    my $newtrans;
    if ( $modified_peptide ){
      my $strand = $newtranscript->start_Exon->strand;

      $newtrans = $self->_recalculate_translation($newtranscript, $strand);

      unless( validate_Translation_coords($newtrans, 1) && is_Transcript_sane($newtrans)
	      && all_exons_are_valid($newtrans, $self->MAX_EXON_LENGTH, 1)
	      && intron_lengths_all_less_than_maximum($newtrans, $self->MAX_INTRON_LENGTH) ){
	print STDERR "problems with this genomewise alternative model, returning original transript.\n";
	$newtrans = $newtranscript;
      }
    }
    else{
      $newtrans = $newtranscript;
    }
    return $newtrans;
  }  # if defined $newtranscript
  else{
    warn("No combination could be built from coding model and designated cDNA/EST evidence.\n");
    return undef;
  }
}

=head2 transcript_from_single_exon_genewise

  Arg [1]    :
  Description: This method will actually do the combination of both
               cdna and genewise gene.  Note that if there is a match
               on one end but not on the other, the code will extend
               one end, but will leave the other as it is in the
               genewise genes. This will explit cdna matches that look
               fine on one end and we disregard the mismatching part.
  Returntype :

=cut

sub transcript_from_single_exon_genewise {
  my ($self, $eg_exon, $gw_exon, $orig_transcript, $exoncount, @e2g_exons) = @_;

  # clone the orig_transcript which just got passed in. It clones the translation too!
  # Will work on this cloned transcript (and its underlying cloned translation) throughout this method.
  # While cloning transcript, we manually cloned the sub_SeqFeatures, because whie clone_Transcript
  # calls clone_Exon, clone_Exon method doesn't handle sub_SeqFeatures.

  my $transcript = clone_Transcript($orig_transcript);
  foreach my $orig_exon(@{$orig_transcript->get_all_Exons()}){
    foreach my $new_exon(@{$transcript->get_all_Exons()}){
      if ( ($orig_exon->start == $new_exon->start) && ($orig_exon->end == $new_exon->end) ) {
        if($orig_exon->sub_SeqFeature && scalar($orig_exon->sub_SeqFeature) > 1 ){
          my @orig_sf    = sort {$a->start <=> $b->start} $orig_exon->sub_SeqFeature;
          foreach my $orig_sf (@orig_sf) {
          my $cloned_sf=clone_Exon($orig_sf);
          # print STDERR "Adding subsf ". $orig_sf->seq_region_start . "-" . $orig_sf->seq_region_end . " to cloned Exon of length " . 
          #    $new_exon->length . " (" . $new_exon->seq_region_start . "-" . $new_exon->seq_region_end. ") in transcript_from_single_exon_genewise.\n";
              $new_exon->add_sub_SeqFeature($cloned_sf,'');
          }
        }
      }  # close if (exon start and end matches)
    } # close foreach new exon
  }  # close foreach original exon

  # Note that $translation is from the cloned transcript, not the transcript we originally passed into this method
  my $translation = $transcript->translation;

  # save out current translation end - we will need this if we have to unmerge frameshifted exons later
  my $orig_tend = $translation->end;

  # stay with being strict about gw vs e2g coords - may change this later ...
  # the overlapping e2g exon must at least cover the entire gw_exon
  if ($gw_exon->start >= $eg_exon->start && $gw_exon->end <= $eg_exon->end){
	
    my $egstart = $eg_exon->start;
    my $egend   = $eg_exon->end;
    my $gwstart = $gw_exon->start;
    my $gwend   = $gw_exon->end;

    # modify the coordinates of the first exon in $newtranscript
    my $ex = $transcript->start_Exon;
	
    $ex->start($eg_exon->start);
    $ex->end($eg_exon->end);
    
    # need to explicitly set the translation start & end exons here.
    $translation->start_Exon($ex);
	
    # end_exon may be adjusted by 3' coding exon frameshift expansion. Ouch.
    $translation->end_Exon($ex);

    # need to deal with translation start and end this time - varies depending on strand
	
    #FORWARD:
    if($gw_exon->strand == 1){
      my $diff = $gwstart - $egstart;
      my $tstart = $translation->start;
      my $tend = $translation->end;

      #print STDERR "diff: ".$diff." translation start : ".$tstart." end: ".$tend."\n";
      #print STDERR "setting new translation to start: ".($tstart+$diff)." end: ".($tend+$diff)."\n";
      $translation->start($tstart + $diff);
      $translation->end($tend + $diff);

      if($translation->start < 0){
	warn("Forward strand: setting very dodgy translation start: " . $translation->start.  "\n");
	return(undef, 1);
      }
      if($translation->end > $translation->end_Exon->length){
	warn("Forward strand: setting dodgy translation end: " . $translation->end .
		   " exon_length: " . $translation->end_Exon->length . "\n");
	return(undef, 1);
      }
    }

    #REVERSE:
    elsif($gw_exon->strand == -1){
      my $diff   = $egend - $gwend;
      my $tstart = $translation->start;
      my $tend   = $translation->end;
      $translation->start($tstart+$diff);
      $translation->end($tend + $diff);

      if($translation->start < 0){
	warn("Forward strand: setting very dodgy translation start: " . $translation->start.  "\n");
	return(undef, 1);
      }
      if($translation->end > $translation->end_Exon->length){
	warn("Forward strand: setting dodgy translation end: " . $translation->end .
		   " exon_length: " . $translation->end_Exon->length . "\n");
	return(undef, 1);
      }
    }
	
    # expand frameshifted single exon genewises back from one exon to multiple exons
    if(defined($ex->sub_SeqFeature) && (scalar($ex->sub_SeqFeature) > 1)){
      print STDERR "frameshift in a single exon genewise\n";
      my @sf = $ex->sub_SeqFeature;

      # save the "current" start and end of frameshift(fs)-merged exon before we bring back (unpack) the individual fs exons
      # NOTE: $ex is the first exon (well, only exon, the fs-merged exon) of the cloned transcript we made at the beginning of this method

      my $cstart = $ex->start;
      my $cend   = $ex->end;
      my $exlength = $ex->length;
      # print STDERR "Current start and end of fs-merged exon are: " . $ex->seq_region_start . " " . $ex->seq_region_end. "\n";

      # Bring back the first exon by modifying the end position of the merged exon
      # (which is already part of the cloned transcript):

      my $first_sf = shift(@sf);
      $ex->end($first_sf->end);
      $transcript->translation->start_Exon($ex);
      # print STDERR "End of first unpacked exon is ". $first_sf->seq_region_end . " (genomic).\n ";

      # Bring back the last exon by engineering a new exon from the last sub_SeqFeature
      # and then adding this new exon to the cloned transcript:

      my $last_sf = pop(@sf);
      $last_sf->end($cend);
      $transcript->add_Exon($last_sf);
      $transcript->translation->end_Exon($last_sf); 
      $transcript->translation->end($orig_tend);  # Putting translation end_Exon isn't enough. we need to say where the translation should end. 

      # print STDERR "End of last unpacked exon is ". $last_sf->seq_region_end ." (genomic).\n";
      # print STDERR "tln end exon end is ".$translation->end_Exon->end."\n";
	
      # If there are more than two frameshift exons, what we had done above is to
      # deal with the start and end exons only.  Need to get any remaining exons:
      foreach my $s(@sf){
        print STDERR "Adding unpacked exon to a coding+UTR combined transcript: ". $s->seq_region_start . " " . $s->seq_region_end . "\n" if $self->VERBOSE;
        $transcript->add_Exon($s);
      }
      
      #$transcript->sort;
      # flush the sub_SeqFeatures
      $ex->flush_sub_SeqFeature;

      if(defined($ex->sub_SeqFeature)){
        my @sf_after = $ex->sub_SeqFeature;
        foreach my $sf_after (@sf_after) {
          print STDERR "BAD!! subseqfeature remained with combined transcript after flush: ". $sf_after->seq_region_start . " " . $sf_after->seq_region_end."\n" if $self->VERBOSE;
        }
      }
    }

    # After unpacking all frameshift-merged exons, the transcript only contains its
    # original exons, where terminal exons could have been extended by bits of UTRs.
    # Still, if the coding model contained 3 exons to start with, it still contains
    # 3 exons at this point.
    # But the UTR can be much longer, i.e. beyond the terminal exons of the coding
    # model. So we need to add back E2G exons, both at 5' and 3':

    $self->add_5prime_exons($transcript, $exoncount, @e2g_exons);
    $self->add_3prime_exons($transcript, $exoncount, @e2g_exons);
  }
  return ($transcript,0); 
}


=head2 transcript_from_multi_exon_genewise

  Arg [1]    : 
  Description: This method will actually do the combination of both
               cdna and genewise gene.  Note that if there is a match
               on one end but not on the other, the code will extend
               one end, but will leave the other as it is in the
               genewise genes. This will explit cdna matches that look
               fine on one end and we disregard the mismatching part.
  Returntype : 

=cut

sub transcript_from_multi_exon_genewise {
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;

  # $current_exon is the exon one the e2g_transcript we are in at the moment
  # $exoncount is the position of the e2g exon in the array

  # save out current translation->end - we'll need it if we have to expand 3prime exon later
  # my $orig_tend = $translation->end;

  my @gwtran  = @{$gw_gene->get_all_Transcripts};
  #$gwtran[0]->sort;
  my @gwexons = @{$gwtran[0]->get_all_Exons};

  my @egtran  = @{$eg_gene->get_all_Transcripts};
  #$egtran[0]->sort;
  my @egexons = @{$egtran[0]->get_all_Exons};

  # in order to match a starting genewise exon with an e2g exon, we need to have
  # a. exactly coinciding exon ends
  # b. exon starts lying within $exon_slop bp of each other.
  # previously we had required e2g start to be strictly <= gw start, but this will lose us some valid UTRs
  # substitute "end" for "start" for 3' ends of transcripts

  # compare to the first genewise exon
  if($gwexons[0]->strand == 1){
    return $self->transcript_from_multi_exon_genewise_forward($current_exon, $transcript, $translation, 
							      $exoncount, $gw_gene, $eg_gene);
  }
  elsif( $gwexons[0]->strand == -1 ){
    return $self->transcript_from_multi_exon_genewise_reverse($current_exon, $transcript, $translation, 
							      $exoncount, $gw_gene, $eg_gene);
  }
}


=head2 transcript_from_multi_exon_genewise_forward

  Arg [1]    : 
  Description: 
  Returntype : 

=cut

sub transcript_from_multi_exon_genewise_forward {
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;

  my $modified_peptide = 0;

  my @gwtran  = @{$gw_gene->get_all_Transcripts};
  #$gwtran[0]->sort;
  my @gwexons = @{$gwtran[0]->get_all_Exons};

  my @egtran  = @{$eg_gene->get_all_Transcripts};
  #$egtran[0]->sort;
  my @egexons = @{$egtran[0]->get_all_Exons};

  # save out current translation->end - we'll need it if we have to expand 3prime exon later
  my $orig_tend = $translation->end;

  my $exon_slop = 20;

  #5_PRIME:
  if (#they have a coincident end
      $gwexons[0]->end == $current_exon->end &&

      # either e2g exon starts before genewise exon
      ($current_exon->start <= $gwexons[0]->start ||

       # or e2g exon is a bit shorter but there are spliced UTR exons as well
       (abs($current_exon->start - $gwexons[0]->start) <= $exon_slop &&
	$current_exon != $egexons[0]))){

    my $current_start = $current_exon->start;
    my $gwstart       = $gwexons[0]->start;

    # this exon will be the start of translation, convention: phase = -1
    my $ex = $transcript->start_Exon;
    $ex->phase(-1);

    # modify the coordinates of the first exon in $newtranscript if
    # e2g is larger on this end than gw.
    if ( $current_exon->start < $gwexons[0]->start ){
      $ex->start($current_exon->start);
    }
    elsif( $current_exon->start == $gwexons[0]->start ){
      $ex->start($gwstart);
      $ex->phase($gwexons[0]->phase);
    }
    # if the e2g exon starts after the gw exon,
    # modify the start only if this e2g exon is not the first of the transcript
    elsif(  $current_start > $gwstart && $exoncount != 0 ) {
      $ex->start($current_exon->start);
    }

    # add all the exons from the est2genome transcript, previous to this one
    transfer_supporting_evidence($current_exon, $ex);
    $self->add_5prime_exons($transcript, $exoncount, @egexons);

    # fix translation start
    if($gwstart >= $current_start){
      # take what it was for the gw gene, and add on the extra
      my $tstart = $translation->start;
      print STDERR "Forward 5': original translation start: $tstart "; ##
      $tstart += ($gwstart - $current_start);
      $translation->start($tstart);
      print STDERR "re-setting translation start to: $tstart\n"; ##
    }

    # only trust a smaller cdna exon if it is not the first of the transcript
    # (it could be a truncated cdna)
    elsif($gwstart < $current_start && $exoncount != 0){

      $modified_peptide = 1;
      #print STDERR "SHORTENING GENEWISE TRANSLATION\n";
      # genewise has leaked over the start. Tougher call - we need to take into account the 
      # frame here as well
      #print STDERR "gw exon starts: $gwstart < new start: $current_start\n";
      #print STDERR "modifying exon, as cdna exon is not the first of transcript-> exoncount = $exoncount\n";

      # $diff is the number of bases we chop from the genewise exon
      my $diff   = $current_start - $gwstart;
      my $tstart = $translation->start;
      warn("this is a case where gw translation starts at $tstart > 1") if ($tstart>1);
      print STDERR "gw translation start: ".$tstart."\n";
      #print STDERR "start_exon: ".$translation->start_Exon->start.
      #"-".$translation->start_Exon->end.
      #" length: ".($translation->start_Exon->end - $translation->start_Exon->start + 1).
      #" phase: ".$translation->start_Exon->phase.
      #" end_phase: ".$translation->start_Exon->end_phase."\n";

      if($diff % 3 == 0) {
	# we chop exactily N codons from the beginning of translation
	$translation->start(1);
      }
      elsif ($diff % 3 == 1) {
	# we chop N codons plus one base
	$translation->start(3);
      }
      elsif ($diff % 3 == 2) {
	# we chop N codons plus 2 bases
	$translation->start(2);
      }
      else {
	$translation->start(1);
	warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");
      }
    }

    else{
      #print STDERR "gw exon starts: $gwstart > new start: $current_start";
      #print STDERR "but cdna exon is the first of transcript-> exoncount = $exoncount, so we don't modify it\n";
    }
    throw("setting very dodgy translation start: " . $translation->start.  "\n")
      unless $translation->start > 0;

  } # end 5' exon

  # 3_PRIME:
  elsif (# they have coincident start
	 $gwexons[$#gwexons]->start == $current_exon->start &&
	
	 # either e2g exon ends after genewise exon
	 ($current_exon->end >= $gwexons[$#gwexons]->end ||
	
	  # or we allow to end before if there are UTR exons to be added
	  (abs($current_exon->end - $gwexons[$#gwexons]->end) <= $exon_slop &&
	   $current_exon != $egexons[$#egexons]))){

    my $end_translation_shift = 0;

    # modify the coordinates of the last exon in $newtranscript
    # e2g is larger on this end than gw.
    my $ex = $transcript->end_Exon;

    # this exon is the end of translation, convention: end_phase = -1
    $ex->end_phase(-1);

    if ( $current_exon->end > $gwexons[$#gwexons]->end ){
      $ex->end($current_exon->end);
    }
    elsif( $current_exon->end == $gwexons[$#gwexons]->end ){
      $ex->end($gwexons[$#gwexons]->end);
      $ex->end_phase($gwexons[$#gwexons]->end_phase);
    }
    # if the e2g exon ends before the gw exon,
    # modify the end only if this e2g exon is not the last of the transcript
    elsif ( $current_exon->end < $gwexons[$#gwexons]->end && $exoncount != $#egexons ){
	
      $modified_peptide = 1;
      #print STDERR "SHORTENING GENEWISE TRANSLATION\n";
      ## fix translation end iff genewise has leaked over - will need truncating
      my $diff   = $gwexons[$#gwexons]->end - $current_exon->end;
      #print STDERR "diff: $diff\n";
      my $tend   = $translation->end;

      my $gw_exon_length   = $gwexons[$#gwexons]->end - $gwexons[$#gwexons]->start + 1;
      my $cdna_exon_length = $current_exon->end - $current_exon->start + 1;
      #print STDERR "gw exon length  : $gw_exon_length\n";
      #print STDERR "cdna exon length: $cdna_exon_length\n";

      my $length_diff = $gw_exon_length - $cdna_exon_length;
      #print STDERR "length diff: ".$length_diff."\n"; # should be == diff

      $ex->end($current_exon->end);

      if($diff % 3 == 0) {
	# we chop exactily N codons from the end of the translation
	# so it can end where the cdna exon ends
	$translation->end($cdna_exon_length);
	$end_translation_shift = $length_diff;
      }
      elsif ($diff % 3 == 1) {
	# we chop N codons plus one base
	# it should end on a full codon, so we need to end translation 2 bases earlier:
	$translation->end($cdna_exon_length - 2);
	$end_translation_shift = $length_diff + 2;
      }
      elsif ($diff % 3 == 2) {
	# we chop N codons plus 2 bases
	# it should end on a full codon, so we need to end translation 1 bases earlier:
	$translation->end($cdna_exon_length - 1);
	$end_translation_shift = $length_diff + 1;
      }
      else {
	# absolute genebuild paranoia 8-)
	$translation->end($cdna_exon_length);
	warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");
      }
      #print STDERR "Forward: translation end set to : ".$translation->end."\n";

    }
    # need to explicitly set the translation end exon for translation to work out
    my $end_ex = $transcript->end_Exon;
    $translation->end_Exon($end_ex);

    # strand = 1
    my $expanded = $self->expand_3prime_exon($ex, $transcript, 1);

    if($expanded){
      # set translation end to what it originally was in the unmerged genewise gene
      # taking into account the diff
      #print STDERR "Forward: expanded 3' exon, re-setting end of translation from ".$translation->end." to orig_end ($orig_tend)- ( length_diff + shift_due_to_phases ) ($end_translation_shift)".($orig_tend - $end_translation_shift)."\n";
      $translation->end($orig_tend - $end_translation_shift);
    }

    # finally add any 3 prime e2g exons
    transfer_supporting_evidence($current_exon, $ex);
    $self->add_3prime_exons($transcript, $exoncount, @egexons);

  } # end 3' exon

  return ($transcript,$modified_peptide);
}

=head2 transcript_from_multi_exon_genewise_reverse

  Arg [1]    : 
  Description: 
  Returntype : 
  Exceptions : 
  Example    : 

=cut

sub transcript_from_multi_exon_genewise_reverse {
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;

  my $modified_peptide = 0;
  my @gwtran  = @{$gw_gene->get_all_Transcripts};
  #$gwtran[0]->sort;
  my @gwexons = @{$gwtran[0]->get_all_Exons};

  my @egtran  = @{$eg_gene->get_all_Transcripts};
  #$egtran[0]->sort;
  my @egexons = @{$egtran[0]->get_all_Exons};

  # save out current translation->end - we'll need it if we have to expand 3prime exon later
  my $orig_tend = $translation->end;

  my $exon_slop = 20;

  # 5_PRIME:
  if ($gwexons[0]->start == $current_exon->start &&
      # either e2g exon ends after gw exon
      ($current_exon->end >= $gwexons[0]->end ||
       # or there are UTR exons to be added
       (abs($current_exon->end - $gwexons[0]->end) <= $exon_slop &&
	$current_exon != $egexons[0]))){

    # sort out translation start
    my $tstart = $translation->start;
    if($current_exon->end >= $gwexons[0]->end){
      # take what it was for the gw gene, and add on the extra
      $tstart += $current_exon->end - $gwexons[0]->end;
      $translation->start($tstart);
    }
    elsif( $current_exon->end < $gwexons[0]->end && $current_exon != $egexons[0] ){
      # genewise has leaked over the start. Tougher call - we need to take into account the
      # frame here as well
      $modified_peptide = 1;
      #print STDERR "SHORTENING GENEWISE TRANSLATION\n";
      #print STDERR "In Reverse strand. gw exon ends: ".$gwexons[0]->end." > cdna exon end: ".$current_exon->end."\n";
      #print STDERR "modifying exon, as cdna exon is not the first of transcript-> exoncount = $exoncount\n";

      my $diff    = $gwexons[0]->end - $current_exon->end;
      my $gwstart = $gwexons[0]->end;
      my $current_start = $current_exon->end;
      my $tstart  = $translation->start;

      #print STDERR "start_exon: ".$translation->start_Exon->start.
      #  "-".$translation->start_Exon->end.
      #    " length: ".($translation->start_Exon->end - $translation->start_Exon->start + 1).
      #      " phase: ".$translation->start_Exon->phase.
      #	" end_phase: ".$translation->start_Exon->end_phase."\n";

      if    ($diff % 3 == 0) { $translation->start(1); }
      elsif ($diff % 3 == 1) { $translation->start(3); }
      elsif ($diff % 3 == 2) { $translation->start(2); }
      else {
	$translation->start(1);
	warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");}
    }

    throw("setting very dodgy translation start: " . $translation->start.  "\n")
      unless $translation->start > 0;

    # this exon is the start of translation, convention: phase = -1
    my $ex = $transcript->start_Exon;
    $ex->phase(-1);

    # modify the coordinates of the first exon in $newtranscript
    if ( $current_exon->end > $gwexons[0]->end){

      ## HERE WAS THE PROBLEM WHICH CHANGED THE SOURCE COORDINATES! ##
      $ex->end($current_exon->end);
      $ex->phase(-1);

    }
    elsif (  $current_exon->end == $gwexons[0]->end){
      $ex->end($gwexons[0]->end);
      $ex->phase($gwexons[0]->phase);
    }
    elsif (  ($current_exon->end < $gwexons[0]->end) && ($current_exon != $egexons[0]) ){
      $ex->end($current_exon->end);
    }

    # need to explicitly set the translation start exon for translation to work out
    $translation->start_Exon($ex);

    transfer_supporting_evidence($current_exon, $ex);
    $self->add_5prime_exons($transcript, $exoncount, @egexons);

  }
  # end 5' exon

  # 3_PRIME:
  elsif ($gwexons[$#gwexons]->end == $current_exon->end &&
	 # either e2g exon starts before gw exon
	 ($current_exon->start <= $gwexons[$#gwexons]->start ||
	  # or there are UTR exons to be added
	  (abs($current_exon->start - $gwexons[$#gwexons]->start) <= $exon_slop &&
	   $current_exon != $egexons[$#egexons]))){
    my $end_translation_shift = 0;

    # this exon is the end of translation, convention: end_phase = -1
    my $ex = $transcript->end_Exon;
    $ex->end_phase(-1);

    # modify the coordinates of the last exon in $newtranscript
    if ( $current_exon->start < $gwexons[$#gwexons]->start ){
      # no need to modify translation->end as the 'end' of this exon has not changed
      $ex->start($current_exon->start);
      $ex->end_phase(-1);
    }
    elsif( $current_exon->start == $gwexons[$#gwexons]->start){
      $ex->start($gwexons[$#gwexons]->start);
      $ex->end_phase($gwexons[$#gwexons]->end_phase);
    }

    # if the e2g exon starts after the gw exon,
    # modify the end only if this e2g exon is not the last of the transcript
    elsif ( $current_exon->start > $gwexons[$#gwexons]->start && $exoncount != $#egexons ){

      $modified_peptide = 1;
      #print STDERR "SHORTENING GENEWISE TRANSLATION\n";
      #print STDERR "In Reverse strand: gw exon start: ".$gwexons[$#gwexons]->start." < cdna exon start: ".$current_exon->start."\n";
      #print STDERR "modifying exon, as cdna exon is not the last of transcript-> exoncount = $exoncount, and #egexons = $#egexons\n";

      ## adjust translation
      my $diff   = $current_exon->start - $gwexons[$#gwexons]->start;
      #print STDERR "diff: $diff\n";
      my $tend   = $translation->end;
	
      my $gw_exon_length   = $gwexons[$#gwexons]->end - $gwexons[$#gwexons]->start + 1;
      my $cdna_exon_length = $current_exon->end - $current_exon->start + 1;
      #print STDERR "gw exon length  : $gw_exon_length\n";
      #print STDERR "cdna exon length: $cdna_exon_length\n";
	
      my $length_diff = $gw_exon_length - $cdna_exon_length;

      # modify the combined exon coordinate to be that of the cdna
      $ex->start($current_exon->start);

      if($diff % 3 == 0) {
	# we chop exactily N codons from the end of the translation
	# so it can end where the cdna exon ends
	$translation->end($cdna_exon_length);
	$end_translation_shift = $length_diff;
      }
      elsif ($diff % 3 == 1) {
	# we chop N codons plus one base
	# it should end on a full codon, so we need to end translation 2 bases earlier:
	$translation->end($cdna_exon_length - 2);
	$end_translation_shift = $length_diff + 2;
      }
      elsif ($diff % 3 == 2) {
	# we chop N codons plus 2 bases
	# it should end on a full codon, so we need to end translation 1 bases earlier:
	$translation->end($cdna_exon_length - 1);
	$end_translation_shift = $length_diff + 1;
      }
      else {
	# absolute genebuild paranoia 8-)
	$translation->end($cdna_exon_length);
	warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");
      }

    }	
    # strand = -1
    my $expanded = $self->expand_3prime_exon($ex, $transcript, -1);

    # need to explicitly set the translation end exon for translation to work out
    my $end_ex = $transcript->end_Exon;
    $translation->end_Exon($end_ex);

    if($expanded){
      # set translation end to what it originally was in the unmerged genewise gene
      #print STDERR "Reverse: expanded 3' exon, re-setting translation exon ".$translation->end.
      # " to original end( $orig_tend ) - shifts_due_to_phases_etc ( $end_translation_shift ) :".
      # ($orig_tend - $end_translation_shift)."\n";
      $translation->end($orig_tend - $end_translation_shift);
    }
    transfer_supporting_evidence($current_exon, $ex);
    $self->add_3prime_exons($transcript, $exoncount, @egexons);

  } # end 3' exon

  return ($transcript, $modified_peptide);
}

=head2 add_5prime_exons

  Arg [1]    : 
  Description: 
  Returntype : 

=cut

sub add_5prime_exons {
  my ($self, $transcript, $exoncount, @e2g_exons) = @_;

  # add all the exons from the est2genome transcript, previous to this one
  # db handle will be screwed up, need to mak new exons from these
  my $c = 0;
  my $modified = 0;
  while($c < $exoncount){
    my $newexon = new Bio::EnsEMBL::Exon;
    my $oldexon = $e2g_exons[$c];
    $newexon->start($oldexon->start);
    $newexon->end($oldexon->end);
    $newexon->strand($oldexon->strand);

    # these are all 5prime UTR exons
    $newexon->phase(-1);
    $newexon->end_phase(-1);
    $newexon->slice($oldexon->slice);
    my %evidence_hash;
    #print STDERR "adding evidence at 5':\n";
    foreach my $sf( @{$oldexon->get_all_supporting_features} ){
      if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	next;
      }
      $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
      #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
      $newexon->add_supporting_features($sf);
    }
    print STDERR "Adding 5prime UTR exon " . $newexon->start . " " . $newexon->end . "\n" if $self->VERBOSE;
    $transcript->add_Exon($newexon);
    $modified = 1;
    $c++;
  }

  if ($modified == 1){
    $transcript->translation->start_Exon->phase(-1);
  }
}

# $exon is the terminal exon in the genewise transcript, $transcript. We need
# to expand any frameshifts we merged in the terminal genewise exon. 
# The expansion is made by putting $exon to be the last (3' end) component, so we modify its
# start but not its end. The rest of the components are added. The translation end will have to be modified,
# this happens in the method _transcript_from_multi_exon....

=head2 expand_3prime_exon

  Arg [1]    : 
  Description: $exon is the terminal exon in the genewise transcript,
               $transcript. We need to expand any frameshifts we
               merged in the terminal genewise exon.  The expansion is
               made by putting $exon to be the last (3 end)
               component, so we modify its start but not its end. The
               rest of the components are added. The translation end
               will have to be modified, this happens in the method
               _transcript_from_multi_exon....
  Returntype : 

=cut

sub expand_3prime_exon{
  my ($self, $exon, $transcript, $strand) = @_;

  if(defined($exon->sub_SeqFeature) && (scalar($exon->sub_SeqFeature) > 1)){
    #print STDERR "expanding 3'prime frameshifted exon $exon in strand $strand: ".
    #$exon->start."-".$exon->end." phase: ".$exon->phase." end_phase: ".$exon->end_phase."\n";
    my @sf = $exon->sub_SeqFeature;

    my $last = pop(@sf);
    #print STDERR "last component: ".$last->start."-".$last->end." phase ".$last->phase." end_phase ".$last->end_phase."\n";

    #print STDERR "setting exon $exon start: ".$last->start." phase: ".$last->phase."\n";  
    $exon->start($last->start); # but don't you dare touch the end!
    $exon->dbID($last->dbID);
    $exon->phase($last->phase);

    # add back the remaining component exons
    foreach my $s(@sf){
      #print STDERR "adding exon: ".$s->start."-".$s->end."\n";
      $transcript->add_Exon($s);
      #$transcript->sort;
    }
    # flush the sub_SeqFeatures so we don't try to re-expand later
    $exon->flush_sub_SeqFeature;
    return 1;
  }

  # else, no expansion
  return 0;
}

=head2 add_3prime_exons

  Arg [1]    : 
  Description: $exoncount tells us which position in the array of e2g
               exons corresponds to the end of the genewise transcript
               so we add back exons 3 to that position.  $exon and
               $transcript are references to Exon and Transcript
               objects.
  Returntype : 

=cut

sub add_3prime_exons {
  my ($self, $transcript, $exoncount, @e2g_exons) = @_;
  # need to deal with frameshifts - 3' exon is a special case as its end might have changed

  # add all the exons from the est2genome transcript, subsequent to this one
  my $c = $#e2g_exons;
  my $modified = 0;
  while($c > $exoncount){
    my $newexon = new Bio::EnsEMBL::Exon;
    my $oldexon = $e2g_exons[$c];
    $newexon->start($oldexon->start);
    $newexon->end($oldexon->end);
    $newexon->strand($oldexon->strand);
	
    # these are all exons with UTR:
    $newexon->phase(-1);
    $newexon->end_phase(-1);
    $newexon->slice($oldexon->slice);
    print STDERR "adding evidence in 3':\n" if $self->VERBOSE;
    my %evidence_hash;
    foreach my $sf( @{$oldexon->get_all_supporting_features }){
      if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	next;
      }
      $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
      #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
      $newexon->add_supporting_features($sf);
    }
    print STDERR "Adding 3prime UTR exon " . $newexon->start . " " . $newexon->end . "\n" if $self->VERBOSE;
    $transcript->add_Exon($newexon);
    $modified = 1;
    $c--;
  }

  if ($modified == 1){
    $transcript->translation->end_Exon->end_phase(-1);
  }
}


=head2 remap_genes

  Description: strictly speaking this is no longer remapping anything...
               checks translation, set start/stop, can set biotype
  Returntype : ref to array of Bio::EnsEMBL::Gene

=cut

sub remap_genes {
  my ($self, $genes, $biotype_suffix) = @_;

  my @remapped_genes = ();
  my $blessed_type   = $self->{'blessed_type'};
  print STDERR "remapping ".scalar @$genes." genes.\n" if $self->VERBOSE;

 GENE:
  foreach my $gene (@$genes) {

    #leave the blessed genes alone
    my $biotype = $gene->biotype;
    #force a centain biotype?
    if($biotype_suffix && (length($biotype_suffix) > 0)){
      $gene->biotype($gene->biotype().$biotype_suffix);
    }
    if(defined $blessed_type && $blessed_type =~ m/$biotype/){
      print STDERR "not remapping ".$biotype."\n" if $self->VERBOSE;
      push(@remapped_genes, $gene);
      next GENE;
    }
    print STDERR "remapping ".$gene->biotype."\n" if $self->VERBOSE;

    my @t = @{$gene->get_all_Transcripts};
    my $tran = $t[0];

    # check that it translates
    unless(validate_Translation_coords($tran, 1)
	   && !contains_internal_stops($tran)
	   && $tran->translate){
      print STDERR "\nERROR: Gene at ".$gene->seq_region_name." ".$gene->seq_region_start."-".$gene->seq_region_end." ".
	    $gene->seq_region_strand." doesn't translate!\n\n";
      push(@remapped_genes, $gene);  ##added
      next GENE;
    }

    foreach my $transcript ( @{$gene->get_all_Transcripts} ){
      if($biotype){
	$transcript->biotype($gene->biotype);
	$transcript->analysis($gene->analysis);
      }
      # set start and stop codons
      set_start_codon($transcript);
      set_stop_codon($transcript);
    }

    push(@remapped_genes, $gene);
  }
  return \@remapped_genes;
}


=head2 validate_gene

  Arg [1]    : Bio::EnsEMBL::Gene
  Description: checks start and end coordinates of each exon of each transcript are sane
  Returntype : 1 if gene is valid, otherwise zero

=cut

sub validate_gene {
  my ($self, $gene) = @_;

  # should be only a single transcript
  my @transcripts = @{$gene->get_all_Transcripts};
  if(scalar(@transcripts) != 1) {
    print STDERR "Rejecting gene - should have one transcript, not " . scalar(@transcripts) . "\n";
    return 0;
  }

  foreach my $transcript(@transcripts){
    foreach my $exon(@{$transcript->get_all_Exons}){
      unless( validate_Exon_coords($exon, 1) ){
	print STDERR "Rejecting gene because of invalid exon\n";
	return 0;
      }
    }
  }

  return 1;
}


=head2 _recalculate_translation

 Arg[1]      : Bio::EnsEMBL::Transcript
 Arg[2]      : strand of the transcript
 Description : a transcript is used as evidence for genomewise
              to recalculate the ORF. The idea is to use this when
              the peptide has been shortened, due to a genewise model
              being incompatible with the cdna splicing. This can happen when 
              genewise cannot find very short exons
              and attaches them to one of the flanking exons.
              We tell genomewise to keep the splice boundaries pretty much
              static, so that we preserve the original splicing structure.
 Returntype  : Bio::EnsEMBL::Transcript

=cut

sub _recalculate_translation {
  my ($self, $mytranscript, $strand) = @_;

  my $this_is_my_transcript = clone_Transcript($mytranscript);

  compute_translation($mytranscript);

  # check that everything is sane:
  unless(validate_Translation_coords($mytranscript, 1)
	 && contains_internal_stops($mytranscript)
	 && $mytranscript->translate){
    print STDERR "Problem with the translation as it didn't pass sanity check. Returning the original transcript\n";
    return $this_is_my_transcript;
  }
  return $mytranscript;
}


=head2 _transfer_evidence

  Arg [1]    : reference to Bio::EnsEMBL::Tanscript $combined_transcript
  Arg [2]    : reference to Bio::EnsEMBL::Transcript $cdna_transcript
  Description: transfers cdna evidence to combined transcript
  Returntype: Bio::EnsEMBL::Transcript

=cut

sub _transfer_evidence {
  my ($self, $combined_transcript, $cdna_transcript) = @_;

  my $first_support_id;
  foreach my $combined_exon(@{$combined_transcript->get_all_Exons}){
    foreach my $cdna_exon(@{$cdna_transcript->get_all_Exons}){
      # exact match or overlap?

# exact match
#      if($combined_exon->start  == $cdna_exon->start &&
#         $combined_exon->end    == $cdna_exon->end &&
#         $combined_exon->strand == $cdna_exon->strand){
#         print STDERR "exact match " . $combined_exon->dbID . " with " . $cdna_exon->dbID . "; transferring evidence\n";
#         transfer_supporting_evidence($cdna_exon, $combined_exon);
#      }

      # overlap - feature boundaries may well be wonky
      if($combined_exon->overlaps($cdna_exon)){
	if($combined_exon->strand != $cdna_exon->strand){
	  print STDERR "OVERLAPPING-BUT-DIFFERENT_STRANDS!\n";
	}
	else{
	  transfer_supporting_evidence($cdna_exon, $combined_exon);

	}
      }
    }
  }

  my $cdna_trans = $cdna_transcript->get_all_Transcripts()->[0];
  foreach my $tsf (@{$cdna_trans->get_all_supporting_features}) {
    print STDERR "adding supporting feature: ".$tsf->hseqname."\n" if $self->VERBOSE;
    $combined_transcript->add_supporting_features($tsf);
  }
  return $combined_transcript;
}


=head2 look_for_both

  Arg [1]    : Bio::EnsEMBL::Gene
  Description: a copy of Steve's look_for_both-script,
               checks phases, etc.
  Returntype : Bio::EnsEMBL::Gene

=cut

sub look_for_both {
  my ($self, $gene) = @_;

  my $time = time;
  my $nupdated_start = 0;
  my $nupdated_end = 0;
  my $metcnt = 1;
  my $maxterdist = $self->MAX_CDS_EXTEND;

  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans->translation) {
      my $tln = $trans->translation;
      my $coding_start = $trans->cdna_coding_start;
      my $orig_coding_start = $coding_start;
      #$trans->sort;
      my $cdna_seq = uc($trans->spliced_seq);
      my @pepgencoords = $trans->pep2genomic(1,1);
      if(scalar(@pepgencoords) > 2) {
	print STDERR "pep start does not map cleanly\n";
	goto TLNEND; # I swore I'd never use this - this code desperately needs a rewrite
      }
      my $pepgenstart = $pepgencoords[0]->start;
      my $pepgenend   = $pepgencoords[$#pepgencoords]->end;

      unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
	print STDERR "pep start maps to gap\n";
	goto TLNEND; # I swore I'd never use this - this code desperately needs a rewrite
      }
      unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
	print STDERR "pep start (end of) maps to gap\n";
	goto TLNEND; # I swore I'd never use this - this code desperately needs a rewrite
      }
  
      print STDERR "Pep genomic location = " . $pepgenstart . " " . $pepgenend . "\n" if $self->VERBOSE;
      
      my $startseq= substr($cdna_seq,$coding_start-1,3);
      print STDERR "cdna seq for pep start = " . $startseq . "\n" if $self->VERBOSE;
      if ($startseq ne "ATG") {
	if ($coding_start > 3) {
	  my $had_stop = 0;
	  while ($coding_start > 3 && !$had_stop) {
	    my $testseq = substr($cdna_seq,$coding_start-4,3);
	    if ($testseq eq "ATG") {
	      print_Translation($trans) if $self->VERBOSE;

	      my @coords = $trans->cdna2genomic($coding_start-3,$coding_start-1,$gene->strand);
	      my $new_start;
	      my $new_end;
	      if(scalar(@coords) > 2) {
		throw "Shouldn't happen - new coding start maps to >2 locations in genome - I'm out of here\n";
	      } elsif (scalar(@coords) == 2) {
		print STDERR "WOW ISN'T NATURE HORRIBLE: new coding start crosses intron\n";
		print STDERR "coord[0] = " . $coords[0]->start . " " . $coords[0]->end ."\n";
		print STDERR "coord[1] = " . $coords[1]->start . " " . $coords[1]->end ."\n";
		if ($gene->strand == 1) {
		  $new_start = $coords[0]->start;
		  $new_end   = $coords[$#coords]->end;
		} else {
		  $new_start = $coords[0]->end;
		  $new_end   = $coords[$#coords]->start;
		}
	      } else {
		$new_start = $coords[0]->start;
		$new_end   = $coords[0]->end;
	      }

	      unless($coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
		print STDERR "Shouldn't happen - new start maps to gap - I'm out of here\n";
		next;
	      }
	      unless($coords[$#coords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
		print STDERR "Shouldn't happen - new start (end of) maps to gap - I'm out of here\n";
		next;
	      }
          
	      print STDERR "genomic pos for new start = " . $new_start . " " . $new_end . "\n" if $self->VERBOSE;
	      
	      if ($new_end - $new_start == 2) {
  
		$nupdated_start++;
  
		my $newstartexon;
		foreach my $exon (@{$trans->get_all_Exons}) {
		  if ($exon->end >= $new_start && $exon->start <= $new_start) {
		    $newstartexon = $exon;
		    last;
		  }
		}
  
		if ($newstartexon == $tln->start_Exon) {
    
		  if ($tln->start_Exon->strand == 1) {
		    $tln->start($new_start - $tln->start_Exon->start + 1);
		  } else {
		    $tln->start($tln->start_Exon->end - $new_end + 1);
		  }

		  # NAUGHTY, but hey I should have to do this - I've changed the translation after all
		  $trans->{'transcript_mapper'} = undef;
		  $trans->{'coding_region_start'} = undef;
		  $trans->{'coding_region_end'} = undef;
		  $trans->{'cdna_coding_start'} = undef;
		  $trans->{'cdna_coding_end'} = undef;
		  
		} else {
		  # find exon
		  if (!defined($newstartexon)) {
		    print STDERR "Failed finding new start exon - how can this be?\n";
		    next;
		  }
		  # create a copy of if and of current start exon (because of phase change)
		  my $copyexon = new Bio::EnsEMBL::Exon(
							-start  => $tln->start_Exon->start,
							-end    => $tln->start_Exon->end,
							-strand => $gene->strand,
						       );
		  my $copynewstartexon = new Bio::EnsEMBL::Exon(
								-start  => $newstartexon->start,
								-end    => $newstartexon->end,
								-strand => $gene->strand,
							       );
		  
		  # $copyexon->phase(0);
		  $copyexon->end_phase($tln->start_Exon->end_phase);
		  $copyexon->contig($tln->start_Exon->contig);
		  if ($tln->start_Exon->stable_id) {
		    $copyexon->stable_id($tln->start_Exon->stable_id . "MET" . $metcnt++);
		    $copyexon->created($time);
		    $copyexon->modified($time);
		    $copyexon->version(1);
		  }
    
		  $copynewstartexon->phase($newstartexon->phase);
		  # $copynewstartexon->end_phase(0);
		  $copynewstartexon->contig($newstartexon->contig);
		  if ($newstartexon->stable_id) {
		    $copynewstartexon->stable_id($newstartexon->stable_id . "MET" . $metcnt++);
		    $copynewstartexon->created($time);
		    $copynewstartexon->modified($time);
		    $copynewstartexon->version(1);
		  }
    
		  # TODO evidence
        
		  if ($copynewstartexon->strand == 1) {
		    $tln->start($new_start - $copynewstartexon->start + 1);
		  } else {
		    $tln->start($copynewstartexon->end - $new_end + 1);
		  }
  
		  # Replace exons in transcript, and fix phases
  
		  my @newexons;
		  my $inrange = 0;
		  foreach my $exon (@{$trans->get_all_Exons}) {
		    if ($inrange) {
		      $exon->phase( $newexons[$#newexons]->end_phase );
		      $exon->end_phase(($exon->length + $exon->phase) % 3);
		    }
		    if ($exon == $tln->start_Exon) {
		      $copyexon->phase( $newexons[$#newexons]->end_phase );
  
		      push @newexons,$copyexon;
		      $inrange = 0;
		    } elsif ($exon == $newstartexon) {
		      push @newexons,$copynewstartexon;
		      $copynewstartexon->end_phase(($exon->length - $tln->start + 1)%3);
		      print STDERR "Setting end_phase on new start exon to " . $copynewstartexon->end_phase . 
			    " l = " . $exon->length . " ts = " . $tln->start . "\n" if $self->VERBOSE;
		      $inrange = 1;
		    } else {
		      push @newexons,$exon;
		    }
		  }
  
  
		  $trans->flush_Exons;
		  foreach my $exon (@newexons) {
		    $trans->add_Exon($exon);
		  }
  
		  # Reset translation start exon
		  if ($tln->end_Exon == $tln->start_Exon) {
		    $tln->end_Exon($copyexon);
		  }
		  $tln->start_Exon($copynewstartexon);
  
		}
		print_Translation($trans);
	      } else {
		print STDERR "Across exons - not handling this\n";
	      }

	      last;
	    } else {
	      if ($testseq =~ /TAA/ or $testseq =~ /TGA/ or $testseq =~ /TAG/) {
		$had_stop = 1;
	      } else {
		$coding_start -= 3;
	      }
	    }
	  }
	} else {
	  print STDERR "Coding region starts between the 1st and 3rd base of the transcript.  Coding start codon isn't ATG ".
                       "but a max of 3 bases upstream is not enough to search for the next nearest ATG. NOT looking into genomic\n"if $self->VERBOSE;
	}
      }
      
    TLNEND:
      {
	my $coding_end = $trans->cdna_coding_end;
	my $orig_coding_end = $coding_end;
	
	#$trans->sort;
    
	my $peplen = $trans->translate->length;
	
	my @pepgencoords = $trans->pep2genomic($peplen,$peplen);
	
	if(scalar(@pepgencoords) > 2) {
	  print STDERR "pep end does not map cleanly\n";
	  next;
	}
    
	my $pepgenstart = $pepgencoords[0]->start;
	my $pepgenend   = $pepgencoords[$#pepgencoords]->end;
	
	unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
	  print STDERR "pep end maps to gap\n";
	  next;
	}
	unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
	  print STDERR "pep start (end of) maps to gap\n";
	  next;
	}
    
	#print "End Pep genomic location = " . $pepgenstart . " " . $pepgenend . "\n";
	
	my $endseq= substr($cdna_seq,$coding_end-3,3);
	my $cdnalen = length($cdna_seq);
  
	#print "cdna seq for pep end = " . $endseq . "\n";
	my $longendseq= substr($cdna_seq,$coding_end-6,12);
	#print "long end seq (-3 to len 12) = $longendseq\n";
	
	#          if (!($endseq ne "TGA" and $endseq ne "TAA" and $endseq ne "TAG")) {
	#            print "Has end " . $trans->translateable_seq . "\n";
	#          } 
	if ($endseq ne "TGA" and $endseq ne "TAA" and $endseq ne "TAG") {
	  if (($cdnalen-$coding_end) > 3) {
	    while (($cdnalen-$coding_end) > 0 && ($coding_end-$orig_coding_end) <= $maxterdist) {
	      my $testseq = substr($cdna_seq,$coding_end,3);
	      #print "Test seq = $testseq\n" if $self->VERBOSE ; 
	     
	      if ($testseq eq "TGA" or $testseq eq "TAA" or $testseq eq "TAG") {

		my @coords = $trans->cdna2genomic($coding_end+1,$coding_end+3,$gene->strand);
		my $new_start;
		my $new_end;
		if(scalar(@coords) > 2) {
		  throw("new end does not map cleanly\n");
		} elsif (scalar(@coords) == 2) {
		  print STDERR "WOW ISN'T NATURE HORRIBLE: new end crosses intron\n";
		  print STDERR "coord[0] = " . $coords[0]->start . " " . $coords[0]->end ."\n";
		  print STDERR "coord[1] = " . $coords[1]->start . " " . $coords[1]->end ."\n";
		  if ($gene->strand == 1) {
		    $new_start = $coords[0]->start;
		    $new_end   = $coords[$#coords]->end;
		  } else {
		    $new_start = $coords[0]->end;
		    $new_end   = $coords[$#coords]->start;
		  }
		} else {
		  $new_start = $coords[0]->start;
		  $new_end   = $coords[0]->end;
		}
		
		unless($coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
		  print STDERR "new start maps to gap\n";
		  next;
		}
		unless($coords[$#coords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
		  print STDERR "new start (end of) maps to gap\n";
		  next;
		}
		
		if ($new_end - $new_start == 2) {
    
		  #print "Sequence of genomic pos of new end = " . $slice->subseq($new_start,$new_end,$gene->strand) . "\n";
		  $nupdated_end++;
    
		  my $newendexon;
		  foreach my $exon (@{$trans->get_all_Exons}) {
		    if ($exon->end >= $new_start && $exon->start <= $new_start) {
		      $newendexon = $exon;
		      last;
		    }
		  }
    
		  if ($newendexon == $tln->end_Exon) {
		    if ($tln->end_Exon->strand == 1) {
		      $tln->end($new_end - $tln->end_Exon->start + 1);
		    } else {
		      $tln->end($tln->end_Exon->end - $new_start + 1);
		    }
		    
		    # NAUGHTY, but hey I should have to do this - I've changed the translation after all
		    $trans->{'transcript_mapper'} = undef;
		    $trans->{'coding_region_start'} = undef;
		    $trans->{'coding_region_end'} = undef;
		    $trans->{'cdna_coding_start'} = undef;
		    $trans->{'cdna_coding_end'} = undef;
		    
		  } else {
		    # find exon
		    if (!defined($newendexon)) {
		      print STDERR  "Failed finding new end exon - how can this be?\n";
		      next;
		    }
		    # create a copy of if and of current end exon (because of phase change)
		    my $copyexon = new Bio::EnsEMBL::Exon(
							  -start  => $tln->end_Exon->start,
							  -end    => $tln->end_Exon->end,
							  -strand => $gene->strand,
							 );
		    my $copynewendexon = new Bio::EnsEMBL::Exon(
								-start  => $newendexon->start,
								-end    => $newendexon->end,
								-strand => $gene->strand,
							       );
      
		    $copyexon->phase($tln->end_Exon->phase);
		    $copyexon->end_phase($tln->end_Exon->end_phase);
		    $copyexon->contig($tln->end_Exon->contig);
		    if ($tln->end_Exon->stable_id) {
		      $copyexon->stable_id($tln->end_Exon->stable_id . "TER" . $metcnt++);
		      $copyexon->created($time);
		      $copyexon->modified($time);
		      $copyexon->version(1);
		    }
      
		    $copynewendexon->phase($newendexon->phase);
		    # $copynewendexon->end_phase(0);
		    $copynewendexon->contig($newendexon->contig);
		    if ($newendexon->stable_id) {
		      $copynewendexon->stable_id($newendexon->stable_id . "TER" . $metcnt++);
		      $copynewendexon->created($time);
		      $copynewendexon->modified($time);
		      $copynewendexon->version(1);
		    }
      
		    # TODO evidence
          
		    if ($copynewendexon->strand == 1) {
		      $tln->end($new_end - $copynewendexon->start + 1);
		    } else {
		      $tln->end($copynewendexon->end - $new_start + 1 );
  
		      my $tercodon = $copynewendexon->seq->subseq($copynewendexon->end - $new_start-1, $copynewendexon->end - $new_start +1);
		      #reverse($tercodon);
		      #$tercodon =~ tr /ACGT/TGCA/;
  		      
		    }
    
		    # Replace exons in transcript, and fix phases
		    my @newexons;
		    my $inrange = 0;
		    foreach my $exon (@{$trans->get_all_Exons}) {
		      if ($inrange) {
			print STDERR "in range exon before phase = " . $exon->phase . " endphase " . $exon->end_phase . "\n" if $self->VERBOSE;
			$exon->phase( $newexons[$#newexons]->end_phase );
			$exon->end_phase(($exon->length + $exon->phase) % 3);
			print STDERR "in range exon after phase = " . $exon->phase . " endphase " . $exon->end_phase . "\n" if $self->VERBOSE;
		      }
		      if ($exon == $tln->end_Exon) {
			my $phase = $exon->phase;
			if ($phase == -1) {
			  $phase = 0;
			}
			if ($exon == $tln->start_Exon) {
			  $copyexon->end_phase(($exon->length - $tln->start + 1)%3);
			} else {
			  $copyexon->end_phase(($exon->length + $exon->phase)%3);
			}
			print STDERR "Setting end_phase on old end exon to " . $copyexon->end_phase . " l = " . $exon->length . "\n" if $self->VERBOSE;
    
			push @newexons,$copyexon;
			$inrange = 1;
		      } elsif ($exon == $newendexon) {
			$copynewendexon->phase( $newexons[$#newexons]->end_phase );
			$copynewendexon->end_phase( -1);
  
			push @newexons,$copynewendexon;
			$inrange = 0;
		      } else {
			push @newexons,$exon;
		      }
		    }
    
		    $trans->flush_Exons;
		    foreach my $exon (@newexons) {
		      $trans->add_Exon($exon);
		    }
		    
		    # Reset translation start exon
		    if ($tln->end_Exon == $tln->start_Exon) {
		      $tln->start_Exon($copyexon);
		    }
		    $tln->end_Exon($copynewendexon);
		    
		  }
		  print_Translation($trans);
		  # print "translateable seq = \n";
		  # print $trans->translateable_seq . "\n"; 
		} else {
		  print STDERR "Across exons - not handling this\n" if $self->VERBOSE;
		}
		last;
	      }
	      $coding_end += 3;
	    }
	  } else {
            print STDERR "Coding region ends between the 3rd last to the last base of the transcript.  Stop codon isn't TGG, TGA or TAG ".
                       "but a max of 3 bases downstream is not enough to search for the next nearest stop codon. NOT looking into genomic\n"if $self->VERBOSE;

	    print STDERR "Not enough bases downstream - NOT looking into genomic\n" if $self->VERBOSE;
	  }
	}
      }
    }

    # These lines force loads from the database to stop attempted lazy
    # loading during the write (which fail because they are to the wrong
    # db)
    $trans->get_all_supporting_features();

    my @exons= @{$trans->get_all_Exons};
    my $get = $trans->translation;
    $trans->translation(undef);

    foreach my $exon (@exons) {
      $exon->stable_id;
      $exon->contig($gene->slice);
      $exon->get_all_supporting_features;
    }
  }

  return($gene);
}


=head2 forward_genewise_clusters

  Arg [1]    : 
  Description: get/set for genewise clusters
  Returntype : 

=cut

sub forward_genewise_clusters{
  my ($self, $cluster_ref) = @_;

  if (!defined($self->{'_forward_genewise_clusters'})) {
    $self->{'_forward_genewise_clusters'} = [];
  }

  if (defined $cluster_ref && scalar(@{$cluster_ref})) {
    push(@{$self->{'_forward_genewise_clusters'}},@{$cluster_ref});
    # store them sorted
    @{$self->{'_forward_genewise_clusters'}} = sort { $a->start <=> $b->start } @{$self->{'_forward_genewise_clusters'}};
  }

  return $self->{'_forward_genewise_clusters'};
}

=head2 reverse_genewise_clusters

  Arg [1]    : 
  Description: get/set for genewise clusters
  Returntype : 
  Exceptions : 
  Example    : 

=cut

sub reverse_genewise_clusters{
  my ($self, $cluster_ref) = @_;

  if (!defined($self->{'_reverse_genewise_clusters'})) {
    $self->{'_reverse_genewise_clusters'} = [];
  }

  if (defined $cluster_ref && scalar(@{$cluster_ref})) {
    push(@{$self->{'_reverse_genewise_clusters'}},@{$cluster_ref});
    # store them sorted
    @{$self->{'_reverse_genewise_clusters'}} = sort { $a->start <=> $b->start } @{$self->{'_reverse_genewise_clusters'}};
  }

  return $self->{'_reverse_genewise_clusters'};
}


=head2 cdna_genes

  Arg [1]    : 
  Description: get/set for e2g gene array
  Returntype : 

=cut

sub cdna_genes {
  my ($self, $genes) = @_;

  if (!defined($self->{'_cdna_genes'})) {
    $self->{'_cdna_genes'} = [];
  }

  if (defined $genes && scalar(@{$genes})) {
    push(@{$self->{'_cdna_genes'}},@{$genes});
  }

  return ($self->{'_cdna_genes'});
}


=head2 prune

  Arg [1]    : (optional) bool
  Description: get/set for option to prune genes

=cut

sub prune {
  my ($self, $bool) = @_;

  if (!defined($self->{'_prune_genes'})) {
    $self->{'_prune_genes'} = undef;
  }

  if (defined $bool) {
    $self->{'_prune_genes'} = $bool;
  }

  return ($self->{'_prune_genes'});
}


=head2 ests

  Arg [1]    : (optional) ref to array with ests
  Description: get/set for ests
  Returntype : array ref with ST objects

=cut

sub ests {
  my ($self, $ests) = @_;

  if (!defined($self->{'_ests'})) {
    $self->{'_ests'} = [];
  }

  if (defined $ests && scalar(@{$ests})) {
    push(@{$self->{'_ests'}},@{$ests});
  }

  return($self->{'_ests'});
}


=head2 ditags

  Arg [1]    : (optional) ref to array with ditags
  Description: get/set ditags
  Returntype : array ref with ditag objects
  Exceptions : none

=cut

sub ditags {
  my ($self, $ditags) = @_;

  if (!defined($self->{'_ditags'})) {
    $self->{'_ditags'} = [];
  }

  if (defined $ditags && scalar(@{$ditags})) {
    push(@{$self->{'_ditags'}}, @{$ditags});
  }

  return($self->{'_ditags'});
}


=head2 gw_genes

  Arg [1]    : 
  Description: get/set for genewise gene array
  Returntype : 
  Exceptions : 
  Example    : 

=cut

sub gw_genes {
  my ($self, $genes) = @_;
  if (!defined($self->{'_gw_genes'})) {
    $self->{'_gw_genes'} = [];
  }

  if (defined $genes && scalar(@{$genes})) {
    push(@{$self->{'_gw_genes'}},@{$genes});
  }

  return $self->{'_gw_genes'};
}

=head2 blessed_genes

  Arg [1]    : 
  Description: get/set for blessed gene array
               clones the genes
  Returntype : 

=cut

sub blessed_genes {
  my ($self, $slice, $genes) = @_;

  if (!defined($self->{'_blessed_genes'})) {
    $self->{'_blessed_genes'} = [];
  }

  if (defined $slice && defined $genes && scalar(@{$genes})) {

    # split input genes into one transcript per gene; keep type the same
  OLDGENE:
    foreach my $gene(@{$genes}){
      foreach my $transcript(@{$gene->get_all_Transcripts}){

	# make a new gene
	my $newgene = new Bio::EnsEMBL::Gene;
	$newgene->biotype($gene->biotype);
	#preserve the xref and dbID of the blessed genes!
	foreach my $dblink (@{$gene->get_all_DBLinks()}){
	  # print STDERR "Adding xref: ".$dblink->display_id." from blessed to cloned-blessed gene.....\n";
	  $newgene->add_DBEntry($dblink);
	}
        $newgene->dbID($gene->dbID);

	# clone transcript
	my $newtranscript = new Bio::EnsEMBL::Transcript;
	$newtranscript->slice($slice);

	# clone translation
	my $newtranslation = new Bio::EnsEMBL::Translation;
        my @seqeds = @{$transcript->translation->get_all_SeqEdits};
        if (scalar(@seqeds)) {
          print "Copying sequence edits\n" if $self->VERBOSE;
          foreach my $se (@seqeds) {
            my $new_se =
              Bio::EnsEMBL::SeqEdit->new(
                -CODE    => $se->code,
                -NAME    => $se->name,
                -DESC    => $se->description,
                -START   => $se->start,
                -END     => $se->end,
                -ALT_SEQ => $se->alt_seq
              );
            my $attribute = $new_se->get_Attribute();
            $newtranslation->add_Attributes($attribute);
          }
        }
        my @support = @{$transcript->get_all_supporting_features};
        if (scalar(@support)) {
          $newtranscript->add_supporting_features(@support);
        }
        $newtranscript->add_Attributes(@{$transcript->get_all_Attributes});

	$newtranscript->translation($newtranslation);

	foreach my $exon(@{$transcript->get_all_Exons}){
	  # clone the exon
	  my $newexon = clone_Exon($exon);
	  # get rid of stable ids
	  $newexon->stable_id('');

	  # if this is start/end of translation, keep that info:
	  if ( $exon == $transcript->translation->start_Exon ){
		$newtranslation->start_Exon( $newexon );
		$newtranslation->start($transcript->translation->start);
	    }
	    if ( $exon == $transcript->translation->end_Exon ){
		$newtranslation->end_Exon( $newexon );
		$newtranslation->end($transcript->translation->end);
	    }
	  $newtranscript->add_Exon($newexon);

	  #add sf
	  $newexon->add_supporting_features(@{$exon->get_all_supporting_features});
	}

	$newgene->add_Transcript($newtranscript);
	push(@{$self->{'_blessed_genes'}}, $newgene);	

      }
    }
  }

  return $self->{'_blessed_genes'};
}

=head2 combined_genes

  Arg [1]    : ref to Bio::EnsEMBL::Gene
  Description: get/set for combined gene array
  Returntype : 

=cut

sub combined_genes {
  my ($self, $genesref) = @_;

  if (!defined($self->{'_combined_genes'})) {
    $self->{'_combined_genes'} = [];
  }
  if (defined $genesref) {
    push(@{$self->{'_combined_genes'}}, $genesref);
  }
#  if (defined $genesref && scalar(@{$genesref})) {
#    push(@{$self->{'_combined_genes'}},@{$genesref});
#  }

  return $self->{'_combined_genes'};
}

=head2 unmatched_genes

  Arg [1]    : array of genes
  Description: get/set for unmatched gene array
  Returntype : none

=cut

sub unmatched_genes {
  my ($self, @genes) = @_;

  if (!defined($self->{'_unmatched_genes'})) {
    $self->{'_unmatched_genes'} = [];
  }

  if(@genes){
    push(@{$self->{'_unmatched_genes'}},@genes);
  }


  return $self->{'_unmatched_genes'};
}


=head2 populate_kill_list

  Arg [1]    : String mol-type
  Description: read ids to ignore from KillList db
  Returntype : ref to hash

=cut

sub populate_kill_list {
  my ($self, $type) = @_;

  my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => $type);
  my %kill_list = %{$kill_list_object->get_kill_list()};

  return \%kill_list;
}

=head2 kill_list

  Arg [1]    : optional hash ref with kill list
  Description: get/set kill list object
  Returntype : ref to hash

=cut

sub kill_list {
  my ($self, $kill_list) = @_;

  if($kill_list){
    $self->{'_kill_list'} = $kill_list;
  }

  return $self->{'_kill_list'};

}


=head2 retrieve_unmerged_gene

  Arg [1]    : Bio::EnsEMBL::Gene
  Description: Returns unmeregd, frameshifted version of input gene
  Returntype : Bio::EnsEMBL::Gene

=cut

sub retrieve_unmerged_gene{
  my ($self, $merged_gene) = @_;

  my %pairs = %{$self->merged_unmerged_pairs()};
  if(!exists $pairs{$merged_gene}){
    print STDERR "Can't retrieve unmerged \n";
    print STDERR "for gene ".$merged_gene->dbID." \n";
  }

  return $pairs{$merged_gene};
}

=head2 merged_unmerged_pairs

  Arg [1]    : 
  Description: get/set for pairs of frameshift merged and unmerged
               genes. Key is merged gene, value is unmerged

  Returntype : 

=cut

sub merged_unmerged_pairs {
  my ($self, $merged_gene, $unmerged_gene) = @_;

  if (!defined($self->{'_merged_unmerged_pairs'})) {
    $self->{'_merged_unmerged_pairs'} = {};
  }

  if ($unmerged_gene && $merged_gene) {
    $self->{'_merged_unmerged_pairs'}{$merged_gene}= $unmerged_gene;
  }

  # hash ref
  return $self->{'_merged_unmerged_pairs'};
}

=head2 modified_unmodified_pairs

  Arg [1]    : 
  Description: get/set for pairs of UTR modified and unmodified genes.
               Key is modified gene, value is unmodified
  Returntype : 

=cut

sub modified_unmodified_pairs {
  my ($self, $modified_gene, $unmodified_gene) = @_;

  if (!defined($self->{'_modified_unmodified_pairs'})) {
    $self->{'_modified_unmodified_pairs'} = {};
  }

  if ($unmodified_gene && $modified_gene) {
    $self->{'_modified_unmodified_pairs'}{$modified_gene}= $unmodified_gene;
  }

  # hash ref
  return $self->{'_modified_unmodified_pairs'};
}


=head2 output

  Arg [1]    : optional RunnableDB output
  Description: get/set for result of Analysis
  Returntype : array ref with genes

=cut

sub output{
  my ($self,@genes) = @_;

  if (!defined($self->{'_output'})) {
    $self->{'_output'} = [];
  }

  if(@genes){
    push(@{$self->{'_output'}},@genes);
  }

  return $self->{'_output'};
}


=head2 _cdna_slice

  Arg [1]    : 
  Description: 
  Returntype : 

=cut

sub _cdna_slice {
  my ($self, $slice) = @_;

  if($slice){
    $self->{'_cdna_slice'} = $slice;
  }

  return $self->{'_cdna_slice'};
}

=head2 _cdna_evidence

  Arg [1]    : 
  Description: 
  Returntype : 

=cut

sub _cdna_evidence {
  my ($self, $cdna_evidence) = @_;

  if($cdna_evidence){
    $self->{'_cdna_evidence'} = $cdna_evidence;
  }

  return $self->{'_cdna_evidence'};
}

=head2 _known_pairs

  Arg [1]    : optional hash ref
  Description: store the links between NM and NP entries
  Returntype : hashref

=cut

sub _known_pairs {
  my ($self, $hashref) = @_;

  if (defined $hashref) {
    $self->{_known_pairs} = $hashref;
  }

  return $self->{_known_pairs};
}

######## database access functions ###########

=head2 CDNA_DB

  Arg [1]    : optional Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set for db storing exonerate alignments of cDNAs 
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : Throws if input isn't DBAdaptor or db parameters haven't been defined

=cut

sub CDNA_DB {
  my( $self, $val ) = @_;

  if ($val){
     $self->{_cdna_dbs} = $val;
  }

  return $self->{_cdna_dbs};
}

=head2 EST_DB

  Arg [1]    : optional Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set for db storing exonerate alignments of ESTs 
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : Throws if input isn't DBAdaptor or db parameters haven't been defined

=cut

sub EST_DB {
  my( $self, $est_db ) = @_;

  if ($est_db){
    $self->{_est_db} = get_db_adaptor_by_string($est_db,1);
  }

  return $self->{_est_db};
}

=head2 DITAG_DB

  Arg [1]    : optional Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set for db storing exonerate alignments of ditags
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : Throws if input isn't DBAdaptor or db parameters haven't been defined

=cut

sub DITAG_DB {
    my( $self, $ditag_db ) = @_;

    if ($ditag_db){
      $self->{_ditag_db} = $ditag_db;
    }

    if((defined $self->{_ditag_db}) && (!$self->{_ditag_db}) && (scalar @{$self->DITAG_TYPE_NAMES})){
      throw("You have defined DITAG_TYPE_NAMES, but no DITAG_DB parameters.\n");
    }
    return $self->{_ditag_db};
}

=head2 BLESSED_DB

  Arg [1]    : optional Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set for db storing blessed gene structures
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : Throws if input isn't DBAdaptor or db parameters haven't been defined

=cut

sub BLESSED_DB {
  my( $self, $blessed_db ) = @_;

  if ($blessed_db){
    $self->{_blessed_db} = $blessed_db;
  }

  return $self->{_blessed_db};
}

=head2 OUTPUT_DB

  Arg [1]    : optional Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set for db in which UTR modified genes should be stored
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : Throws if input isn't DBAdaptor or db parameters haven't been defined

=cut

sub OUTPUT_DB {
  my( $self, $output_db ) = @_;

  if ($output_db){
    $self->{_output_db} = $output_db;
  }
  if(!$self->{_output_db}){
    throw("Please define database parameters for output db.\n");
  }

  return $self->{_output_db};
}

####  other config variable get/set methods  ######

=head2 INPUT_GENETYPES

  Arg [1]    : optional parameter
  Description: setter / getter for config vars
  Returntype : config value

=cut

sub INPUT_GENETYPES {
  my ($self, $input_genetypes) = @_;

  if($input_genetypes){
    $self->{'_input_genetypes'} = $input_genetypes;
  }

  return $self->{'_input_genetypes'};
}

sub BLESSED_GENETYPES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_blessed_genetypes} = $value;
  }
  return $self->{_blessed_genetypes};
}

sub VERBOSE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_verbose} = $value;
  }
  return $self->{_verbose};
}

sub MAX_EXON_LENGTH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_max_exon_length} = $value;
  }
  return $self->{_max_exon_length};
}

sub MAX_CDS_EXTEND {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_max_cds_extend} = $value;
  }
  return $self->{_max_cds_extend};
}

sub DITAG_TYPE_NAMES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_ditag_type_names} = $value;
  }
  return $self->{_ditag_type_names};
}

sub cDNA_GENETYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_cdna_genetype} = $value;
  }
  return $self->{_cdna_genetype};
}

sub EXTEND_ORIGINAL_BIOTYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_extend_original_biotype} = $value;
  }
  return $self->{_extend_original_biotype};
}

sub EXTEND_BIOTYPE_OF_UNCHANGED_GENES {
  my ($self,$value) = @_;

  if (defined $value) {
	  $value = "_".$value  unless ( $value =~ m/^_/ );
    $self->{_extend_biotype_of_unchanged_genes} = $value;
  }
  return $self->{_extend_biotype_of_unchanged_genes};
}   

sub KNOWN_UTR_GENETYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_known_utr_genetype} = $value;
  }
  return $self->{_known_utr_genetype};
}

sub UTR_GENETYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_utr_genetype} = $value;
  }
  return $self->{_utr_genetype};
}

sub PRUNE_GENES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_prune_genes} = $value;
  }
  return $self->{_prune_genes};
}

sub FILTER_ESTS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_filter_ests} = $value;
  }
  return $self->{_filter_ests};
}

sub MAX_INTRON_LENGTH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_max_intron} = $value;
  }
  return $self->{_max_intron};
}

sub EST_GENETYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_est_genetype} = $value;
  }
  return $self->{_est_genetype};
}

sub BLESSED_UTR_GENETYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_blessed_utr_genetype} = $value;
  }
  return $self->{_blessed_utr_genetype};
}

sub LOOK_FOR_KNOWN {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_look_for_known} = $value;
  }
  return $self->{_look_for_known};
}

sub DITAG_WINDOW {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_ditag_window} = $value;
  }
  return $self->{_ditag_window};
}

sub KNOWNUTR_FILE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_knownutr_file} = $value;
  }
  return $self->{_knownutr_file};
}

sub INPUT_GENES {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'INPUT_GENES'} = $arg ;
  }
  return $self->{'INPUT_GENES'} ;
}

sub ABINITIO_SETS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'ABINITIO_SETS'} = $arg ;
  }
  return $self->{'ABINITIO_SETS'} ;
}

sub SIMGW_SETS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'SIMGW_SETS'} = $arg ;
  }
  return $self->{'SIMGW_SETS'} ;
}

sub EST_SETS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'EST_SETS'} = $arg ;
  }
  return $self->{'EST_SETS'} ;
}

sub AB_INITIO_LOGICNAMES {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'AB_INITIO_LOGICNAMES'} = $arg ;
  }
  return $self->{'AB_INITIO_LOGICNAMES'} ;
}

sub OUTPUT_DATABASE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'OUTPUT_DATABASE'} = $arg ;
  }
  return $self->{'OUTPUT_DATABASE'} ;
}

sub FILTER_SINGLETONS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'FILTER_SINGLETONS'} = $arg ;
  }
  return $self->{'FILTER_SINGLETONS'} ;
}

sub FILTER_NON_CONSENSUS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'FILTER_NON_CONSENSUS'} = $arg ;
  }
  return $self->{'FILTER_NON_CONSENSUS'} ;
}

sub ADD_UTR {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'ADD_UTR'} = $arg ;
  }
  return $self->{'ADD_UTR'} ;
}

sub MIN_CONSENSUS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'MIN_CONSENSUS'} = $arg ;
  }
  return $self->{'MIN_CONSENSUS'} ;
}

sub UTR_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'UTR_PENALTY'} = $arg ;
  }
  return $self->{'UTR_PENALTY'} ;
}

sub END_EXON_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'END_EXON_PENALTY'} = $arg ;
  }
  return $self->{'END_EXON_PENALTY'} ;
}

sub EST_OVERLAP_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'EST_OVERLAP_PENALTY'} = $arg ;
  }
  return $self->{'EST_OVERLAP_PENALTY'} ;
}

sub SHORT_INTRON_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'SHORT_INTRON_PENALTY'} = $arg ;
  }
  return $self->{'SHORT_INTRON_PENALTY'} ;
}

sub SHORT_EXON_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'SHORT_EXON_PENALTY'} = $arg ;
  }
  return $self->{'SHORT_EXON_PENALTY'} ;
}

sub GOOD_PERCENT {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'GOOD_PERCENT'} = $arg ;
  }
  return $self->{'GOOD_PERCENT'} ;
}

sub GOOD_BIOTYPE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'GOOD_BIOTYPE'} = $arg ;
  }
  return $self->{'GOOD_BIOTYPE'} ;
}

sub BAD_BIOTYPE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'BAD_BIOTYPE'} = $arg ;
  }
  return $self->{'BAD_BIOTYPE'} ;
}

sub SMALL_BIOTYPE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'SMALL_BIOTYPE'} = $arg ;
  }
  return $self->{'SMALL_BIOTYPE'} ;
}

sub RNASEQ_INTRON_NAME {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'RNASEQ_INTRON_NAME'} = $arg ;
  }
  return $self->{'RNASEQ_INTRON_NAME'} ;
}

sub RNASEQ_INTRON_DB {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'RNASEQ_INTRON_DB'} = $arg ;
  }
  return $self->{'RNASEQ_INTRON_DB'} ;
}

sub INPUT_DBS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'INPUT_DBS'} = $arg ;
  }
  return $self->{'INPUT_DBS'} ;
}


1;
