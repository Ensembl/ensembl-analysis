=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels - 

=head1 SYNOPSIS


=head1 DESCRIPTION

This module helps to remove RNASeq models coming from the merged set
but that are very different from the models build from the tissues
alone. Most of the time the models from the tissues would be better.
It only change the biotype for the models

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::RecoverRNASeqTissuesModels;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(clone_Gene empty_Gene); 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::RecoverRNASeqTissuesModels;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::DB::HTS;

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
            my $bestgene = $self->get_longest_translation_gene($single_set_clust->get_Genes);
            $self->output([$bestgene]);
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
#    my @set2_orphan_genes;
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
      my $bam = Bio::DB::HTS->new(
          -bam => $self->INTRON_BAM_FILE,
          -autoindex => 1,
      );
      my $cluster_counter = 0;
      foreach my $twoway_clust (@$twoway_clusters) {
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
        if ($set1_start != $set2_start or $set1_end != $set2_end) {
            my $is_set2_better = 0;
            my $bestgene = $self->get_longest_translation_gene(\@set2_genes_overlapped_with_set1_genes);
            my $use_recovered_gene = 0;
            foreach my $set1_gene (@set1_genes_overlapped_with_set2_genes) {
                $is_set2_better = check_intron_support($bam, $set1_gene, $bestgene, $self->INTRON_ALLOWANCE);
                if ($is_set2_better > 0) {
                    print STDERR 'You might want to delete this model: ', $set1_gene->display_id, 'Start: ',$set1_gene->seq_region_start ,' End: ', $set1_gene->seq_region_end,' Strand: ', $set1_gene->strand, "\n";
                    if ( !($bestgene->seq_region_end < $set1_gene->seq_region_start or $set1_gene->seq_region_end  < $bestgene->seq_region_start )) {
##                        print STDERR 'DEBUG DELETE ', $bestgene->seq_region_end, ' < ', $set1_gene->seq_region_start, ' OR ', $set1_gene->seq_region_end , ' < ', $bestgene->seq_region_start, "\n";
                        $set1_gene->biotype($self->DELETE_BIOTYPE);
                    }
#                    print STDERR 'DEBUG DELETED ', $bestgene->seq_region_end, ' < ', $set1_gene->seq_region_start, ' OR ', $set1_gene->seq_region_end , ' < ', $bestgene->seq_region_start, "\n";
                    $use_recovered_gene = 1;
                }
                elsif ($is_set2_better == 0) {
                    if ($bestgene->get_all_Transcripts->[0]->translation->length > $set1_gene->get_all_Transcripts->[0]->translation->length) {
                        print STDERR 'DEBUG: WARNING set2 translation is longer than set1', "\n";
                        $set1_gene->biotype($self->DELETE_BIOTYPE);
                        $use_recovered_gene = 1;
                    }
                }
            }
            $self->output([$bestgene]) if ($use_recovered_gene);
        }
        $self->output(\@set1_genes_overlapped_with_set2_genes);
        $self->output(\@set2_genes_overlapped_with_set1_genes);
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
        empty_Gene($gene);
        $gene_adaptor->store($gene);
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

    my @genes;
    my @biotypes;
    foreach my $dbname (@{$dbs}) {
        my $dbadaptor = $self->_get_db($dbname);
        if (!$dbadaptor) {
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


=head2 check_intron_support

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
    my ($bam, $set1_gene, $set2_gene, $intron_allowance) = @_;
    my @set1_exons = sort {$a->start <=> $b->start} @{$set1_gene->get_all_Exons()};
    my @set2_exons = sort {$a->start <=> $b->start} @{$set2_gene->get_all_Exons()};
    my $set1_count_support = 0;
    my $set2_count_support = 0;
    my @count_support_set1;
    my @count_support_set2;
    for (my $i = 0; $i < (scalar(@set1_exons)-1); $i++) {
        my $count_support = count_support($bam, $set1_gene->seq_region_name, $set1_exons[$i]->seq_region_end, $set1_exons[$i+1]->seq_region_start);
        push(@count_support_set1, $count_support);
        $set1_count_support += $count_support;
    }
    for (my $i = 0; $i < (scalar(@set2_exons)-1); $i++) {
        my $count_support = count_support($bam, $set2_gene->seq_region_name, $set2_exons[$i]->seq_region_end, $set2_exons[$i+1]->seq_region_start);
        push(@count_support_set2, $count_support);
        $set2_count_support += $count_support;
    }
    my $is_modified = 0;
    if (@count_support_set1 > 1) {
        if ($count_support_set1[0] < $count_support_set1[1]/100) {
            $set1_count_support -= $count_support_set1[0];
            shift(@set1_exons);
            $is_modified = 1;
        }
        if ($count_support_set1[$#count_support_set1] < $count_support_set1[$#count_support_set1-1]/100) {
            $set1_count_support -= $count_support_set1[$#count_support_set1];
            pop(@set1_exons);
            $is_modified = 1;
        }
    }
    if (@count_support_set2 > 1) {
        if ($count_support_set2[0] < $count_support_set2[1]/100) {
            $set2_count_support -= $count_support_set2[0];
            shift(@set2_exons);
        }
        if ($count_support_set2[$#count_support_set2] < $count_support_set2[$#count_support_set2-1]/100) {
            $set2_count_support -= $count_support_set2[$#count_support_set2];
            pop(@set2_exons);
        }
    }
    warning( 'Set1: ', $set1_count_support, ' Set2: ', $set2_count_support, "\n") if ($set2_count_support > $set1_count_support);
    return 1 if ($set2_count_support > $set1_count_support);
    if (scalar(@set1_exons) == scalar(@set2_exons)+1) {
#        print STDERR 'DEBUG check_intron_support: set1 ',scalar(@set1_exons), ' == set2+1 ', scalar(@set2_exons), "\n";
        if ($set2_count_support+($set2_count_support/$intron_allowance) > $set1_count_support) {
#            print STDERR 'DEBUG check_intron_support: better model set1 ', $set1_count_support, ' == set2 ', $set2_count_support, "\n";
            return 1;
        }
    }
    if ($set1_count_support == $set2_count_support) {
        if ($is_modified) {
            print STDERR 'DEBUG check_intron_support: support count is equal but set1 has been modified, so set2 should be better', "\n";
            return 1;
        }
        return 0;
    }
    return -1;
}


=head2 count_support

  Arg [1]   : Listref of String
  Arg [2]   : Hashref, the key is the name of the database
               the values are the biotypes to retrieve
  Function  : Fetch all the genes for all the database specified
              if no biotype is given, it will fetch all the genes
  Returntype: Listref of Bio::EnsEMBL::Gene and String of biotypes
  Exceptions: none
  Example   : 

=cut

sub count_support {
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
    print STDERR $count;
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

sub get_longest_translation_gene {
  my ($self, $genes) = @_;

  print STDOUT "\nLooking for gene with longest translation...\n";

  my $num_groups = $self->NUM_GROUPS;
  my $biotype = $self->GOOD_BIOTYPE;
  my $max_exon_gene;
  my $max_exon_length = 0;
  my $max_exon_cnt = 0;
  my $longest_translation_gene;
  my $longest_exons_gene;
  my $max_translation_length = 0;

  my $num_genes = scalar(@$genes);
  my %most_abundant;
  my %introns_score;
  foreach my $this_set_only_clust_gene (@$genes) {
#  return ($id, $ex_cnt, $translation_length, $exons_length, $max_translation);
    my ($id, $ex_cnt, $translation_length, $exons_length, $max_translation) = getting_max($this_set_only_clust_gene, \%introns_score);
    if ($ex_cnt > $max_exon_cnt) {
      $max_exon_cnt = $ex_cnt;
      $max_exon_gene = $this_set_only_clust_gene;
    }

    if ($translation_length > $max_translation_length) {
      $max_translation_length = $translation_length;
      $longest_translation_gene = $this_set_only_clust_gene;
    }

    if ($exons_length > $max_exon_length) {
      $max_exon_length = $exons_length;
      $longest_exons_gene = $this_set_only_clust_gene;
    }
    $most_abundant{$id}{count}++;
    push(@{$most_abundant{$id}{genes}}, $this_set_only_clust_gene);
  }
  my @ids = sort {$most_abundant{$b}{count} <=> $most_abundant{$a}{count}} keys %most_abundant;
  my $longestID = $longest_translation_gene->dbID;
  my $maxexonID = $max_exon_gene->dbID;
  my $longexonID = $longest_exons_gene->dbID;
  my $selected_gene;
#  print STDERR 'DEBUG: ', $longest_translation_gene->display_id, '/', $longestID, '@', $max_translation_length , ' -- ', $max_exon_gene->display_id, '/', $maxexonID, '@', $max_exon_cnt , ' -- ', $longest_exons_gene->display_id, '/', $longexonID, '@', $max_exon_length, "\n";
  if (@ids > 1) {
      if (@ids > 3) {
          $biotype = 'fragmented';
          my $num_ids = @ids;
#          print STDERR 'DEBUG: FRAG ', $num_ids, "\n";
          if ($most_abundant{$ids[0]}{count} >= ($num_ids*0.66)) {
              $biotype = 'hqfragmented';
#              print STDERR 'DEBUG: hqfragmented: ', $most_abundant{$ids[0]}{genes}[0]->dbID, "\t", $most_abundant{$ids[0]}{count}, "\n";
          }
#          foreach my $id (@ids) {
#              print STDERR '   ID: ', $id, ': ', $most_abundant{$id}{count}, "\n";
#          }
      }
      else {
          if ($num_genes < $num_groups/10) {
              $biotype = 'verylow';
#              print STDERR 'DEBUG: Very low ', $num_genes, ' but you have ', $num_groups, "\n";
          }
          elsif ($num_genes < $num_groups/2) {
              $biotype = 'low';
#              print STDERR 'DEBUG: Low ', $num_genes, ' but you have ', $num_groups, "\n";
          }
          elsif ($most_abundant{$ids[0]}{count} < ($num_genes/2)) {
              $biotype = 'weak';
#            print STDERR 'DEBUG: Weak ', $most_abundant{$ids[0]}{genes}[0]->dbID, ': ', $most_abundant{$ids[0]}{count}, "\n";
          }
          elsif ($most_abundant{$ids[0]}{count} == $most_abundant{$ids[1]}{count}) {
              print STDERR 'INFO: ', $most_abundant{$ids[0]}{count}, ' => ', $most_abundant{$ids[0]}{genes}[0]->dbID, ' & ', $most_abundant{$ids[1]}{genes}[0]->dbID, "\n";
          }
      }
  }
  if ($longestID == $maxexonID and $longestID == $longexonID and (grep {$_->dbID == $longexonID} @{$most_abundant{$ids[0]}{genes}} or ($most_abundant{$ids[1]}{count} > $most_abundant{$ids[0]}{count}/2 and grep {$_->dbID == $longexonID} @{$most_abundant{$ids[1]}{genes}})) ) {
#    print STDERR 'DEBUG: the longest model is significant', "\n";
    $biotype = 'longest';
  }
  elsif (grep {$_->dbID == $longestID} @{$most_abundant{$ids[0]}{genes}}) {
      if ($longestID != $longexonID) {
          if (grep {$_->dbID == $longexonID} @{$most_abundant{$ids[0]}{genes}}) {
              if (max_translation($longest_translation_gene) eq max_translation($longest_exons_gene)) {
                  $selected_gene = $longest_exons_gene;
              }
          }
          else {
#        print STDOUT "WARNING: gene with longest translation $max_translation_length (gene id $longestID) is different from gene with maximum exon count $max_exon_cnt (gene id $maxexonID). Gene with longest translation is (always) chosen.";
            $biotype .= "_warning";
              foreach my $gene (@{$most_abundant{$ids[0]}{genes}}) {
                  $selected_gene = $gene if ($gene->start < $longest_translation_gene->start and $gene->end >= $longest_translation_gene->end);
              }
#              print STDERR 'DEBUG: WARNING longest exon is not in the most abundant', "\n";
          }
      }
      else {
#        print STDOUT "GREAT! gene with longest translation $max_translation_length (gene id $longestID) is the same as gene with maximum exon count $max_exon_cnt (gene id $maxexonID).";
        $biotype .= "_great";
      }
  }
  else {
#              print STDERR 'DEBUG: ', $longestID, ' <=> ', $longexonID, "\n";
      if (grep {$_->dbID == $longexonID} @{$most_abundant{$ids[0]}{genes}}) {
          $selected_gene = $longest_exons_gene;
#          print STDERR 'DEBUG: using gene with maximum number of exons ', $selected_gene->dbID, ' => ', $longest_exons_gene->dbID, "\n";
      }
      else {
#          print STDERR 'DEBUG: COUNT: ',$most_abundant{$ids[0]}{count} , ' <==> ',$most_abundant{$ids[1]}{count} , "\n";
          if ($most_abundant{$ids[1]}{count} > $most_abundant{$ids[0]}{count}*$self->ABUNDANCE_THRESHOLD) {
              $selected_gene = $longest_translation_gene;
#              print STDERR 'DEBUG: using gene with longest translation ', $selected_gene->dbID, "\n";
          }
          elsif ($most_abundant{$ids[0]}{count} > $num_groups/10) {
              if (grep {$_->dbID == $maxexonID} @{$most_abundant{$ids[0]}{genes}}) {
                  $selected_gene = $max_exon_gene;
#                  print STDERR 'DEBUG: using gene with max number of exons ', $selected_gene->dbID, "\n";
              }
              else {
                  if ($longestID == $longexonID and $longestID == $maxexonID) {
                      $selected_gene = $longest_translation_gene;
#                      print STDERR 'DEBUG: using gene with longest translation as selected is significant', $selected_gene->dbID, "\n";
                  }
                  else {
                  $selected_gene = $most_abundant{$ids[0]}{genes}[0];
#                  print STDERR 'DEBUG: may have been using gene with longest translation ', $selected_gene->dbID, "\n";
                  }
              }
          }
          else {
              $selected_gene = $most_abundant{$ids[0]}{genes}[0];
#              print STDERR 'DEBUG: using gene with maximum number of copies ', $selected_gene->dbID, ' => ', $most_abundant{$ids[0]}{genes}[0]->dbID, "\n";
          }
      }
  }
  if ($num_genes < $num_groups/2) {
    $biotype .= '_low';
  }
  $selected_gene = $longest_translation_gene unless ($selected_gene);
  my $min_start = 100000000000;
  my $max_start = 0;
  my $min_score = 10000;
  my $min_score_start = 0;
  my $full_score = 0;
  my $num_introns = 0;
  foreach my $transcript (@{$selected_gene->get_all_Transcripts}) {
      foreach my $ise (@{$transcript->get_all_IntronSupportingEvidence}) {
          if (exists $introns_score{$ise->start.':'.$ise->end}) {
              ++$num_introns;
              my $start = $ise->start;
#              print STDERR 'DEBUG SCORE ', $start, ' ^ ', $introns_score{$start.':'.$ise->end}, "\n";
              $min_start = $start if ($min_start > $start);
              $max_start = $start if ($max_start < $start);
              if ($introns_score{$start.':'.$ise->end} < $min_score) {
                  $min_score_start = $start;
                  $min_score = $introns_score{$start.':'.$ise->end};
              }
              $full_score += $introns_score{$start.':'.$ise->end};
          }
      }
  }
  if ($min_score < $full_score/($num_introns*100) and $min_score_start > $min_start and $min_score_start < $max_start) {
#      print STDERR 'DEBUG SPLITTING ', $min_score, ' < ', ($full_score/($num_introns*100)), "\n";
      if ($selected_gene->dbID != $longestID and grep {$_->dbID == $longestID} @{$most_abundant{$ids[0]}{genes}}) {
          $selected_gene = $longest_translation_gene;
#          print STDERR 'DEBUG: SPLITTING using gene with longest translation ', $selected_gene->dbID, "\n";
      }
      elsif ($selected_gene->dbID != $longexonID and grep {$_->dbID == $longexonID} @{$most_abundant{$ids[0]}{genes}}){
          $selected_gene = $longest_exons_gene;
#          print STDERR 'DEBUG: SPLITTING using gene with longest exons ', $selected_gene->dbID, "\n";
      }
      elsif ($selected_gene->dbID != $maxexonID and grep {$_->dbID == $maxexonID} @{$most_abundant{$ids[0]}{genes}}) {
          $selected_gene = $max_exon_gene;
#          print STDERR 'DEBUG: SPLITTING using gene with max exons ', $selected_gene->dbID, "\n";
      }
      else {
          $selected_gene = $most_abundant{$ids[0]}{genes}[0];
#          print STDERR 'DEBUG: SPLITTING using most abundant gene ', $selected_gene->dbID, "\n";
      }
  }

  my $bestgene = clone_Gene($selected_gene);
  $bestgene->biotype($biotype);
  $bestgene->analysis($self->analysis);
  $bestgene->source($self->analysis->logic_name);
# Updating the depth score to our new model
  foreach my $transcript (@{$bestgene->get_all_Transcripts}) {
      foreach my $ise (@{$transcript->get_all_IntronSupportingEvidence}) {
          if (exists $introns_score{$ise->start.':'.$ise->end}) {
              $ise->score($introns_score{$ise->start.':'.$ise->end});
          }
      }
  }
  print STDERR 'CHOSEN ', $bestgene->display_id, ' / ', $bestgene->dbID, ' @ ', $bestgene->seq_region_start, ' - ', $bestgene->seq_region_end, "\n";

  return $bestgene;
}

sub max_translation {
  my ($gene) = @_;
  my $max_transl_len = 0;
  my $max_translation;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    my $trans_len = $trans->translation->length;
    if ($trans_len > $max_transl_len) {
      $max_transl_len = $trans_len;
      $max_translation = $trans->translation;
    }
  }
  return $max_translation;
}

sub getting_max {
  my ($gene, $href_introns_score) = @_;
  my $id;
  my $max_num_coding_exons = 0;
  my $max_length_exons = 0;
  my $max_length_translation = 0;
  my $max_translation;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    my $num_coding_exons = 0;
    foreach my $coding_exon (@{$trans->get_all_translateable_Exons}) {
        $num_coding_exons++;
        $id .= $coding_exon->start.':'.$coding_exon->end.':';
    }
    if ($num_coding_exons > $max_num_coding_exons) {
      $max_num_coding_exons = $num_coding_exons;
    }
    my $length_exons = 0;
    foreach my $exon (@{$trans->get_all_Exons}) {
      $length_exons += $exon->length;
    }
    if ($length_exons > $max_length_exons) {
        $max_length_exons = $length_exons;
    }
    my $trans_len = $trans->translation->length;
    if ($trans_len > $max_length_translation) {
      $max_length_translation = $trans_len;
      $max_translation = $trans->translation;
    }
    foreach my $ise (@{$trans->get_all_IntronSupportingEvidence}) {
        my $id = $ise->seq_region_start.':'.$ise->seq_region_end;
        if (exists $href_introns_score->{$id}) {
            $href_introns_score->{$id} += $ise->score;
        }
        else {
            $href_introns_score->{$id} = $ise->score;
        }
    }
  }
  return ($id, $max_num_coding_exons, $max_length_translation, $max_length_exons, $max_translation);
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

sub NUM_GROUPS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_NUM_GROUPS'} = $arg ;
  }
  return $self->{'_NUM_GROUPS'} ;
}

sub GOOD_BIOTYPE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_GOOD_BIOTYPE'} = $arg ;
  }
  return $self->{'_GOOD_BIOTYPE'} ;
}

sub INTRON_ALLOWANCE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_INTRON_ALLOWANCE'} = $arg ;
  }
  return $self->{'_INTRON_ALLOWANCE'} ;
}

sub ABUNDANCE_THRESHOLD {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'_ABUNDANCE_THRESHOLD'} = $arg ;
  }
  return $self->{'_ABUNDANCE_THRESHOLD'} ;
}


1;
