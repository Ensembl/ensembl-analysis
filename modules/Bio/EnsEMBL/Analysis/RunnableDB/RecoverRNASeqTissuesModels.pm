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

Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels - 

=head1 SYNOPSIS


=head1 DESCRIPTION

This module helps to remove RNASeq models coming from the merged set
but that are very different from the models build from the tissues
alone. Most of the time the models from the tissues would be better.
It only change the biotype for the models

=head1 METHODS

=cut


# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/RecoverRNASeqTissuesModels.pm,v $
# $Version: $
package Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils; 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::RecoverRNASeqTissuesModels;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::DB::Sam;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels
  Function  : instatiates a RecoverRNASeqTissuesModels object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels
  Exceptions: 
  Example   : 

=cut


sub new {
  my ($class,@args) = @_ ;
  my $self = $class->SUPER::new(@args) ;
  $self->read_and_check_config($RECOVER_RNASEQ_MODELS_CONFIG_BY_LOGIC) ;

  return $self ;
}


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels
  Function  : fetch input (gene-objects) out of different databases specified
              in the RecoverRNASeqTissuesModels config.
  Returntype: void
  Exceptions: none
  Example   : 

=cut

sub fetch_input {
  my ($self) = @_ ;

  $self->throw("No input id") unless defined($self->input_id) ;
  my ($merged_genes, $merged_biotypes) = $self->fetch_all_genes($self->MERGED_SET_DATABASE, $self->MERGED_BIOTYPES);
  $self->merged_genes($merged_genes);
  $self->MERGED_BIOTYPES($merged_biotypes);
  my ($tissues_genes, $tissues_biotypes) = $self->fetch_all_genes($self->TISSUES_SET_DATABASE, $self->TISSUES_BIOTYPES);
  $self->tissues_genes($tissues_genes);
  $self->TISSUES_BIOTYPES($tissues_biotypes);
  $self->_disconnect_connections;
}


=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels
  Function  : Run the analysis
  Returntype: void
  Exceptions: none
  Example   : 

=cut

sub run {
  my ($self) = @_;

  my $set1_name = 'merged_set';
  my $set2_name = 'tissues_set';
  my $set1_only_cluster_nr = 0;
  my $set2_only_cluster_nr = 0;
  my $set1_orphan_gene_nr = 0;
  my $set2_orphan_gene_nr = 0;
  my $twoway_cluster_nr = 0;
  my $set1_genes_in_twoway_clust_nr = 0;
  my $set2_genes_in_twoway_clust_nr = 0;
  my ($curr_sr_start) = $self->input_id =~ /[^:]+:[^:]+:[^:]+:(\d+)/;

  my %types_hash;
  $types_hash{$set1_name} = $self->MERGED_BIOTYPES;
  $types_hash{$set2_name} = $self->TISSUES_BIOTYPES;

  my ($clustered,$unclustered) = cluster_Genes( [@{$self->merged_genes}, @{$self->tissues_genes}] , \%types_hash , $self->CODING_ONLY ,  $self->IGNORE_STRAND);

  my $clustered_count = scalar @$clustered;
  my $unclustered_count = scalar @$unclustered;
  print STDOUT "clustered_count is $clustered_count and unclustered_count is $unclustered_count\n";
  ##### First we look at single-set clusters which contain multiple set1 models exclusively
  ##### or multiple set2 models exclusively.


  if ($self->CHECK_SINGLE) {
    my $single_set_clusters = get_single_clusters($clustered);

    my $single_set_count = scalar @$single_set_clusters;
    if (scalar @$single_set_clusters) {
      logger_info("\n\tLooking at $set1_name-only or $set2_name-only $single_set_count clusters.\n");
 
      foreach my $single_set_clust ( @$single_set_clusters) {

        my $oneway_clust_sr_start = $curr_sr_start + $single_set_clust->start;
        my $oneway_clust_sr_end = $curr_sr_start + $single_set_clust->end;
 
        my @sets = @{ $single_set_clust->get_sets_included() } ; 
        for my $set ( @sets ) {   
          if ( $set eq $set1_name ) {
            print("Skipping single-cluster from $set1_name...\n");
            #$set1_only_cluster_nr ++ ;
            #$self->report_single_set_clusters($single_set_clust, $set1_name, $self->input_id, $oneway_clust_sr_start, $oneway_clust_sr_end);
          } elsif ( $set eq $set2_name ) {
            print("Looking at single-cluster from $set2_name...\n");
            $set2_only_cluster_nr ++;
            $self->report_single_set_clusters($single_set_clust, $set2_name, $self->input_id, $oneway_clust_sr_start, $oneway_clust_sr_end);

            # get the gene with the longest translation and write it to the output db
            print STDOUT "\n";
          }
        }
        $self->output($single_set_clust->get_Genes);
      }
    }
    else {
      print STDOUT "\tNo single-set clusters.\n";
    }
  } # End if ($self->CHECK_SINGLE)
  

  ##### Second, we look at unclustered genes, i.e. one set1 gene or one set2 gene
  ##### per locus with no corresponding genes from the same locus in another set.

  if ($self->CHECK_ORPHAN) {
    my (@set1_orphan_genes, @set2_orphan_genes);
    foreach my $orphan_cluster(@$unclustered) {
      push ( @set1_orphan_genes, @{$orphan_cluster->get_Genes_by_Set($set1_name)} );
      push ( @set2_orphan_genes, @{$orphan_cluster->get_Genes_by_Set($set2_name)} );
    }

    logger_info("\n\tLooking at $set2_name lone genes....\n");

    my $set1_count_orphan_genes = scalar(@set1_orphan_genes);
    if ( $set1_count_orphan_genes ) {
      print("Looking at $set1_count_orphan_genes orphan genes in $set1_name.\n");
      $set1_orphan_gene_nr += scalar(@set1_orphan_genes);
      $self->report_orphan_genes(\@set1_orphan_genes, $set1_name);
      $self->output(\@set1_orphan_genes);
    } else {
      print STDOUT "\tNo $set1_name orphan genes.\n";
    }
    my $set2_count_orphan_genes = scalar(@set2_orphan_genes);
    if ( $set2_count_orphan_genes ) {
      print("Looking at $set2_count_orphan_genes orphan genes in $set2_name.\n");
      $set2_orphan_gene_nr += scalar(@set2_orphan_genes);
      $self->report_orphan_genes(\@set2_orphan_genes, $set2_name);
      $self->output(\@set2_orphan_genes);
    } else {
      print STDOUT "\tNo $set2_name orphan genes.\n";
    }
  } # End if ($self->CHECK_ORPHAN)


  ##### Third, we look at cases where set1 and set2 genes overlap
  #####
  if ($self->CHECK_TWOWAY) {
    logger_info("\n\tLooking at clusters containing both $set1_name and $set2_name genes.\n");
    my $twoway_clusters = get_twoway_clusters($clustered);
    print STDOUT "\t" . scalar@$twoway_clusters . " two_way clusters found.\n\n";

    if ( scalar(@$twoway_clusters) )  {
      my $bam = Bio::DB::Sam->new(
          -bam => $self->INTRON_BAM_FILE,
          -autoindex => 1,
      );
      my $cluster_counter = 0;
      TWO_WAY_CLUSTER: foreach my $twoway_clust (@$twoway_clusters) {
        $twoway_cluster_nr ++;
        $cluster_counter ++;
        my $twoway_clust_sr_start = $curr_sr_start  + $twoway_clust->start;
        my $twoway_clust_sr_end = $curr_sr_start + $twoway_clust->end;

        logger_info("\t\tLooking at Cluster number $cluster_counter: ".$self->input_id." $twoway_clust_sr_start - $twoway_clust_sr_end.\n");

        my @set1_genes_overlapped_with_set2_genes = @{$twoway_clust->get_Genes_by_Set($set1_name)};
        my @set2_genes_overlapped_with_set1_genes = @{$twoway_clust->get_Genes_by_Set($set2_name)};

        my $set1_genes_count = scalar(@set1_genes_overlapped_with_set2_genes);
        my $set2_genes_count = scalar(@set2_genes_overlapped_with_set1_genes);

        $set1_genes_in_twoway_clust_nr += $set1_genes_count;
        $set2_genes_in_twoway_clust_nr += $set2_genes_count;
        my $set1_start = 10000000000;
        my $set1_end = 0;
        my $set2_start = 10000000000;
        my $set2_end = 0;
        my %h_set1_boundaries;
        my %h_set2_boundaries;
        foreach my $gene (@set1_genes_overlapped_with_set2_genes) {
            $set1_start = $gene->start if ($gene->start < $set1_start);
            $set1_end = $gene->end if ($gene->end > $set1_end);
            $h_set1_boundaries{$gene->start}++;
            $h_set1_boundaries{$gene->end}++;
        }
        foreach my $gene (@set2_genes_overlapped_with_set1_genes) {
            $set2_start = $gene->start if ($gene->start < $set2_start);
            $set2_end = $gene->end if ($gene->end > $set2_end);
            $h_set2_boundaries{$gene->start}++;
            $h_set2_boundaries{$gene->end}++;
        }
        my $is_set2_better = 0;
        if ($set1_start != $set2_start or $set1_end != $set2_end) {
            if (scalar(keys %h_set2_boundaries) == 2) {
                $is_set2_better += check_intron_support($bam, $set1_genes_overlapped_with_set2_genes[0], $set2_genes_overlapped_with_set1_genes[0]);
            }
            else {
                foreach my $model (@set2_genes_overlapped_with_set1_genes) {
                    $is_set2_better += check_intron_support($bam, $set1_genes_overlapped_with_set2_genes[0], $model);
                }
            }
            if ($is_set2_better) {
                print STDERR 'You might want to delete this model: ', $set1_genes_overlapped_with_set2_genes[0]->display_id, 'Start: ',$set1_genes_overlapped_with_set2_genes[0]->start ,' End: ', $set1_genes_overlapped_with_set2_genes[0]->end,' Strand: ', $set1_genes_overlapped_with_set2_genes[0]->strand, "\n";
                $self->_change_biotype(\@set1_genes_overlapped_with_set2_genes, $self->DELETE_BIOTYPE);
                $self->output(\@set2_genes_overlapped_with_set1_genes);
            }
        }
        $self->output(\@set1_genes_overlapped_with_set2_genes);
      }
    }
  }  # end if ($self->CHECK_TWOWAY)

  print STDOUT "\n\n";
  print STDOUT "**************************** SUMMARY **********************************\n";
  print STDOUT "***********************************************************************\n";

  if ($self->CHECK_SINGLE) {
    print STDOUT "Found $set1_only_cluster_nr clusters with $set1_name genes only.\n";
    print STDOUT "Found $set2_only_cluster_nr clusters with $set2_name genes only.\n";
  }
  if ($self->CHECK_ORPHAN) {
    print STDOUT "Found $set1_orphan_gene_nr $set1_name orphan (unclustered) genes.\n";
    print STDOUT "Found $set2_orphan_gene_nr $set2_name orphan (unclustered) genes.\n";
  }
  if ($self->CHECK_TWOWAY) {
    print STDOUT "Found $twoway_cluster_nr two_way clusters containing both $set1_name and $set2_name genes, ".
              "involving $set1_genes_in_twoway_clust_nr $set1_name genes and $set2_genes_in_twoway_clust_nr $set2_name genes.\n";
  }
}


=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels
  Function  : Write the genes in the OUTPUT_DATABASE database
  Returntype: void
  Exceptions: none
  Example   : 

=cut

sub write_output {
    my ($self) = @_;

    my $dbadaptor = $self->get_dbadaptor($self->OUTPUT_DATABASE);
    my $gene_adaptor = $dbadaptor->get_GeneAdaptor;
    foreach my $gene (@{$self->output}) {
        $gene->load;
        $gene_adaptor->store($gene);
    }
}


=head2 _change_biotype

  Arg [1]   : Listref of Bio::EnsEMBL::Gene
  Arg [2]   : String, the new biotype
  Function  : Change the biotype for a list of Bio::EnsEMBL::Gene
  Returntype: void
  Exceptions: none
  Example   : 

=cut

sub _change_biotype {
    my ($self, $genes, $biotype) = @_;

    foreach my $gene (@$genes) {
        $gene->biotype($biotype.'_'.$gene->biotype);
    }
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub fetch_all_genes {
    my ($self, $dbs, $biotypes) = @_;

    my @db_adaptors;
    my @genes;
    my @biotypes;
    foreach my $dbname (@{$dbs}) {
        my $dbadaptor;
        if ($self->_get_db($dbname)) {
            $dbadaptor = $_;
        }
        else {
            $dbadaptor = $self->get_dbadaptor($dbname);
            $self->_cache_db($dbname, $dbadaptor);
        }
        my $slice = $dbadaptor->get_SliceAdaptor->fetch_by_name($self->input_id);
        if (scalar(@{$biotypes->{$dbname}})) {
          push(@genes, @{_fetch_genes_by_specified_biotypes($biotypes->{$dbname}, $slice)});
          push(@biotypes, @{$biotypes->{$dbname}});
        }  else {
          my ($fetched_biotypes, $fetched_genes) = _fetch_genes_of_all_biotypes ($slice);
          push(@biotypes, @$fetched_biotypes);
          push(@genes, @$fetched_genes);
        }
    }
    return \@genes, \@biotypes;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub check_intron_support {
    my ($bam, $set1_gene, $set2_gene) = @_;
    my @set1_exons = sort {$a->start <=> $b->start} @{$set1_gene->get_all_Exons()};
    my @set2_exons = sort {$a->start <=> $b->start} @{$set2_gene->get_all_Exons()};
    my $set1_count_support = 0;
    my $set2_count_support = 0;
    for (my $i = 0; $i < (scalar(@set1_exons)-1); $i++) {
        $set1_count_support += count_suppport($bam, $set1_gene->seq_region_name, $set1_exons[$i]->end, $set1_exons[$i+1]->start);
    }
    for (my $i = 0; $i < (scalar(@set2_exons)-1); $i++) {
        $set2_count_support += count_suppport($bam, $set2_gene->seq_region_name, $set2_exons[$i]->end, $set2_exons[$i+1]->start);
    }
    warn( 'Set1: ', $set1_count_support, ' Set2: ', $set2_count_support, "\n") if ($set2_count_support > $set1_count_support);
    return 1 if ($set2_count_support > $set1_count_support);
    return 0;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub count_suppport {
    my ($bam, $seq_region_name, $exon1_end, $exon2_start) = @_;
    my $segment = $bam->segment($seq_region_name, $exon1_end, $exon2_start);
    throw("Bam file segment not found for slice " .  $seq_region_name . "\n") unless $segment;
    my $count = 0;
    my $iterator = $segment->features(-iterator=>1);
    while (my $read = $iterator->next_seq) {
        my $ugfs = ungapped_features($read);
        foreach my $ugf (@{$ugfs}) {
            next unless ($ugf->[1] == $exon1_end and $ugf->[2] == $exon2_start);
            ++$count;
        }
    }
    return $count;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub ungapped_features {
  my ($read) = @_;
  my @ugfs;
  
  my $string = $read->cigar_str;
  my $start = $read->start;
  my $end = $read->end;
  my @pieces = ( $string =~ /(\d*[MDN])/g );
  foreach my $piece (@pieces) {
    my ($length) = ( $piece =~ /^(\d*)/ );
    if( $length eq "" ) { $length = 1 }
    if( $piece =~ /M$/ ) {
      $start += $length;
    } else {
      my ( $qstart, $qend);
      $qstart = $start-1;
      $start += $length;
      $qend = $start;
      
      my $ugf;
      $ugf->[0] = $read->query->name;
      $ugf->[1] = $qstart;
      $ugf->[2] = $qend;
      push @ugfs, $ugf;
    }
  }
  return \@ugfs;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub report_single_set_clusters {
  my ($self, $single_set_clust, $set_name, $sr_name, $clust_start, $clust_end) = @_;
  print STDOUT $set_name."-only-cluster: $sr_name\t$clust_start\t$clust_end.\n";

  my @this_set_only_cluster_genes = @ {$single_set_clust->get_Genes};
  foreach my $this_set_only_clust_gene(@this_set_only_cluster_genes) {
    print STDOUT $set_name."-only-cluster gene: ";
    $self->print_gene_info($this_set_only_clust_gene);
  }
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub report_orphan_genes {
  my ($self, $orphan_genes, $set_name) = @_;
  print STDOUT "\t".scalar(@$orphan_genes) . " $set_name orphan genes found.\n" ;
  for my $orphan_gene( @$orphan_genes ) {  
    print STDOUT "$set_name orphan: ";
    $self->print_gene_info($orphan_gene);

    # write gene to the output db
  }
  $self->output($orphan_genes);
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub print_gene_info {
  my ($self, $gene) = @_;
  my $transcripts = $gene->get_all_Transcripts;

  my $hit_name = '';
  if (defined ${${$transcripts}[0]->get_all_supporting_features}[0]) {
    $hit_name =  ${${$transcripts}[0]->get_all_supporting_features}[0]->hseqname;
  }

    print STDOUT $gene->display_id. " ";

  print STDOUT $gene->slice->seq_region_name . " " . $gene->seq_region_start . " " . $gene->seq_region_end . " " . $gene->biotype .
        " max_cod_ex_cnt: ". max_coding_exon_count($gene);

  if ($self->CHECK_STOP) {
    my $n_trans_with_stop = 0;
    my $total_stops_in_gene = 0;

    ($n_trans_with_stop, $total_stops_in_gene) = check_internal_stops($transcripts);
    print STDOUT " $n_trans_with_stop trans_have_stops, total_nr_of_stops_in_gene: $total_stops_in_gene";
  }
  if ($self->CHECK_FRAMESHIFT) {
    my $fshift_in_gene_cnt = check_frameshift_intron($transcripts);
    print STDOUT " $fshift_in_gene_cnt trans_have_fshift";
  }
  print STDOUT " ".$hit_name.".";
  print STDOUT "\n";
}

sub max_coding_exon_count {
  my ($gene) = @_;
  my $max_coding_ex_cnt = 0;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    my $ex_cnt = scalar(@{$trans->get_all_translateable_Exons});
    if ($ex_cnt > $max_coding_ex_cnt) {
      $max_coding_ex_cnt = $ex_cnt;
    }
  }
  return $max_coding_ex_cnt;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub check_frameshift_intron {
  my ($transcripts) = @_;
  my $prev_exon=undef;
  my $n_trans_with_fshift ++;

  foreach my $transcript(@$transcripts) {
    FEX: foreach my $exon (@{$transcript->get_all_translateable_Exons}) {
      if (defined($prev_exon)) {
        my $intron_len;
        if ($transcript->strand == 1) {
           $intron_len = $exon->start - $prev_exon->end -1;
        } else {
           $intron_len = $prev_exon->start - $exon->end -1;
        }

        if ($intron_len <= 17) {
          $n_trans_with_fshift++;
          last FEX;
        }
      }
      $prev_exon = $exon;
    }
  }
  return $n_trans_with_fshift;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub check_internal_stops {
  my ($transcripts) = @_;
  my $stoppy_trans_cnt = 0;
  my $total_stops_cnt_in_gene = 0;

  foreach my $trans (@$transcripts) {
    if ($trans->translation) {
      my $translation_len = $trans->translate->length();
      my $tlnseq = $trans->translate->seq;
      if ($tlnseq =~ /\*/) {
        $stoppy_trans_cnt++;
        $total_stops_cnt_in_gene += ($tlnseq =~ tr/*/*/);
      }
    }
  }
  return ($stoppy_trans_cnt, $total_stops_cnt_in_gene);
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub _cache_db {
    my ($self, $dbname, $dbadaptor) = @_;

    $self->{'_db_cache'}{$dbname} = $dbadaptor;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub _get_db {
    my ($self, $dbname) = @_;

    return values %{$self->{'_db_cache'}} unless $dbname;
    return $self->{'_db_cache'}{$dbname} if (exists $self->{'_db_cache'}{$dbname});
    return undef;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub _disconnect_connections {
    my ($self) = @_;

    foreach my $db ($self->_get_db) {
        $db->dbc->disconnect_when_inactive(1);
    }
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub _fetch_genes_by_specified_biotypes {
  my ($biotypes, $slice) = @_;
  my @genes;
  foreach my $biotype (@$biotypes) {
    push ( @genes , @{$slice->get_all_Genes_by_type($biotype)} );
  }
  return \@genes;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub _fetch_genes_of_all_biotypes {
  my ($slice) = @_;
  my $genes = $slice->get_all_Genes();
  my %all_biotypes;
  foreach my $gene(@$genes) {
    $all_biotypes{$gene->biotype} = 1;
  }
  my @biotypes = keys(%all_biotypes);
  return (\@biotypes, $genes);
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub merged_genes {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_merged_genes'} = $arg ;
  }
  return $self->{'_merged_genes'} ;
}


=head2 fetch_all_genes

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub tissues_genes {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_tissues_genes'} = $arg ;
  }
  return $self->{'_tissues_genes'} ;
}

=head2 INPUT_GENES 

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype:
  Exceptions: 
  Example   : 

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::RecoverRNASeqTissuesModels


sub VERBOSE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_VERBOSE'} = $arg ;
  }
  return $self->{'_VERBOSE'} ;
}

sub IGNORE_STRAND {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_IGNORE_STRAND'} = $arg ;
  }
  return $self->{'_IGNORE_STRAND'} ;
}

sub CODING_ONLY {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_CODING_ONLY'} = $arg ;
  }
  return $self->{'_CODING_ONLY'} ;
}

sub CHECK_TWOWAY {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_CHECK_TWOWAY'} = $arg ;
  }
  return $self->{'_CHECK_TWOWAY'} ;
}

sub CHECK_STOP {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_CHECK_STOP'} = $arg ;
  }
  return $self->{'_CHECK_STOP'} ;
}

sub CHECK_FRAMESHIFT {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_CHECK_FRAMESHIFT'} = $arg ;
  }
  return $self->{'_CHECK_FRAMESHIFT'} ;
}

sub CHECK_SINGLE {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_CHECK_SINGLE'} = $arg ;
  }
  return $self->{'_CHECK_SINGLE'} ;
}

sub CHECK_ORPHAN {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_CHECK_ORPHAN'} = $arg ;
  }
  return $self->{'_CHECK_ORPHAN'} ;
}

sub DELETE_BIOTYPE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_DELETE_BIOTYPE'} = $arg ;
  }
  return $self->{'_DELETE_BIOTYPE'} ;
}

sub MERGED_SET_DATABASE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_MERGED_SET_DATABASE'} = $arg ;
  }
  return $self->{'_MERGED_SET_DATABASE'} ;
}

sub TISSUES_UNIPROT_DATABASE {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_TISSUES_UNIPROT_DATABASE'} = $arg ;
  }
  return $self->{'_TISSUES_UNIPROT_DATABASE'} ;
}

sub MERGED_BIOTYPES {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_MERGED_BIOTYPES'} = $arg ;
  }
  return $self->{'_MERGED_BIOTYPES'} ;
}

sub TISSUES_SET_DATABASE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_TISSUES_SET_DATABASE'} = $arg ;
  }
  return $self->{'_TISSUES_SET_DATABASE'} ;
}

sub MERGED_UNIPROT_DATABASE {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'_MERGED_UNIPROT_DATABASE'} = $arg ;
  }
  return $self->{'_MERGED_UNIPROT_DATABASE'} ;
}

sub TISSUES_BIOTYPES {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_TISSUES_BIOTYPES'} = $arg ;
  }
  return $self->{'_TISSUES_BIOTYPES'} ;
}

sub OUTPUT_DATABASE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_OUTPUT_DATABASE'} = $arg ;
  }
  return $self->{'_OUTPUT_DATABASE'} ;
}

sub INTRON_BAM_FILE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_INTRON_BAM_FILE'} = $arg ;
  }
  return $self->{'_INTRON_BAM_FILE'} ;
}


1;
