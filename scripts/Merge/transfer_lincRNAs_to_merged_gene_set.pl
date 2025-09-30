#!/usr/bin/env perl
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

=pod

=head1 NAME 

transfer_lincRNAs_to_merged_gene_set.pl

=head1 DESCRIPTION

This script transfers lincRNAs from an existing genebuild (e.g. release version
X) to an updated/patched genebuild (release version X+1) given the underlying
genome assembly has not changed.  The updated/patched genebuild won't be 
complete for handover unless lincRNAs have been imported.

The "patched" genebuild referred to here is often an updated ensembl-havana
merged gene set.

The script starts by fetching lincRNAs from the current core DB and genes from
the patched core DB.  It then clusters the lincRNAs and patched genes with one
another.  Depending on the outcome of the clustering analysis, lincRNAs will
be treated differently:

(1) If a lincRNA does not overlap with any patched gene and the lincRNA has Ensembl supporting evidence
    (if it hasn't, ignore because if HAVANA have annotated it, it will have been there after the merge),
    then copy the lincRNA to the patched core DB, delete HAVANA supporting evidence (if any) and
    label as ensembl_lincrna.

(2) If a lincRNA overlaps with a protein_coding gene, we discard the lincRNA
    (not copied to the patched core DB).

(3) If a lincRNA overlaps with a processed_transcript gene and the lincRNA has Ensembl supporting evidence
    (if it hasn't, ignore because if HAVANA have annotated it, it will have been there after the merge),
    and the processed_transcript gene corresponds to 'havana' (we ignore ensembl_havana_gene because
    they were merged as a result of combining hav processed_transcript with ens protein_coding, so 
    the lincrna should have been discarded against the ens protein_coding)
    We also do *not* copy the lincRNA model across to the patched DB, but we do the following
    changes:
      (i) gene logic name to "ensembl_havana_lincrna", to indicate agreement between
    ensembl and havana.
     (ii) gene biotype to "lincRNA".
    (iii) transcript logic_name to "ensembl_havana_lincrna".


=head1 OPTIONS

=head2 DB connection:

=head3 DB containing lincRNAs

        -dbuser         Username to connect to DB containing lincRNAs

        -dbpass         Password to connect to DB containing lincRNAs 
                        (optional, in most cases, all you need is read-only)

        -dbhost         Where the lincRNA DB is

        -dbname         lincRNA DB Name

        -dbport         port to connect to lincRNA DB


=head3 DB containing patched gene set

        -newdbuser      Username to connect to DB containing patched gene set

        -newdbpass      Password to connect to DB containing patched gene set.
                        Optional. However, you must use this option if you want
                        to write lincRNA genes directly into this DB while the
                        script is run.

        -newdbhost      Where the DB containing patched gene set is

        -newdbname      Pathced gene set DB Name

        -newdbport      Port to connect to patched gene set DB

        -----------------------------------------------------------------------

=head2 Fetching genes for clustering:

=head3 Compulsory options:

        -set1_name      Name of your lincRNA gene set. Default set as "lincRNA".

        -set2_name      Name of your patched gene set which does not contain
                        lincRNAs.  Default set as "merged".

        -coordsystem    The coordinate system from which your seq_region slices
                        are to be fetched from (e.g. "chromosome", "toplevel",
                        "supercontig")

        -path           The assembly version name (e.g. "GRCh37", "NCBIM37")


=head3 Optional flags:
 
        -biotypes       Optional.  Use this if you want to fetch specific
                        biotypes of genes from the patched set for clustering.
                        Biotypes should be separated by commas but no white
                        spaces, e.g. protein_coding,pseudogene. 
                        Otherwise, the script fetches all genes from patched
                        DB by default.

        -non_ref        A boolean flag to indicate if genes from non-reference
                        genome sequences (e.g. haplotypic regions in human) 
                        should be retrieved too.  This option is set to false
                        by default (i.e. only fetching from reference genome 
                        sequence).
 
        -duplicate_sr   A boolean flag to indicate whether duplicated reference
                        seq_regions (e.g. PAR on Y-chr)  should be retrieved. 
                        This option is set to false by default (i.e. not 
                        fetching duplicated seq_regions)

        -----------------------------------------------------------------------

=head2 Output options:

        -verbose        use this option to get more print statements to follow
                        the script.  Set to 0 (not verbose) by default.

=head3 For a dry run

        -gsi_to_copy    The name of a text file containing gene stable IDs of
                        lincRNAs worth copying.  This allows checking of the 
                        output of the script before any changes to DBs are being
                        made. Either use this option together with -sql_update, 
                        or use "-write" on its own.

        -sql_update     The name of a text file containing SQL update statements
                        which change logic_name of processed_transcript genes 
                        (in patched gene set) to "ensembl_havana_lincrna" if they 
                        overlap with lincRNAs. Use this option together with
                        -gsi_to_copy.  Don't use this option if you're using
                        "-write".

=head3 For a "wet" run (writing directly to patched DB when the script is run)

        -write          A boolean flag to indicate whether lincRNAs should be 
                        copied directly into the patched DB and whether patched 
                        gene logic_name should be changed directly when it 
                        overlaps with a lincRNA. This option is set to false by
                        default.  Only turn on this flag if you're absolutely 
                        sure that the script has identified the correct set of 
                        lincRNAs for transfer.
                        Don't use this flag with "gsi_to_copy" or "sql_update".

                        Also, make sure the output logic_name "ensembl_havana_
                        lincrna" already exists in your output DB before you
                        run the script.

        -----------------------------------------------------------------------



=head1 EXAMPLES

=head2 Example 1

To compare lincRNAs with all genes from the patched set, print out stable IDs
of lincRNAs worth copying + SQL statements for changing logic_name of some
patched genes:

perl transfer_lincRNAs_to_merged_gene_set.pl -dbuser read_only -dbhost live_host\
-dbname lincRNA_DB -newdbuser read_only -newdbhost patched_DB_host\
-newdbname patched_genes_DB -coordsystem chromosome -path NCBIM37 \
-gsi_to_copy lincRNAs_to_copy.GSI -sql_update proc_trans_logic_name_update.sql \
-verbose >& compare_lincRNA_vs_patched_genes.log

=head2 Example 2

To compare lincRNAs with only protein_coding genes in patched set, and copy all
"fit" lincRNAs into the patched DB:

perl transfer_lincRNAs_to_merged_gene_set.pl -dbuser read_only -dbhost live_host\
-dbname lincRNA_DB -newdbuser read_only -newdbhost patched_DB_host\
-newdbname patched_genes_DB -coordsystem chromosome -path NCBIM37 \
-biotypes protein_coding -write -verbose >& compare_lincRNA_vs_patched_genes.log

=cut


use strict;
use warnings; 
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Getopt::Long qw(:config no_ignore_case);

$|=1;

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname);
my ($newdbhost, $newdbuser, $newdbpass, $newdbport, $newdbname);
my ($vegadbhost, $vegadbuser, $vegadbpass, $vegadbport, $vegadbname);

my $set1_name = 'lincRNA';
my $set2_name = 'merged';

my @biotypes;
my $coordsystem;
my $path;
my $non_ref = 0;
my $duplicate_sr = 0;
my $gsi_to_copy = undef;
my $sql_update = undef;
my $write = 0;
my $verbose = 0;


GetOptions( 'dbhost|host|h=s'      => \$dbhost,
            'dbname|db|D=s'      => \$dbname,
            'dbuser|user|u=s'      => \$dbuser,
            'dbpass|pass|p:s'      => \$dbpass,
            'dbport|port|P=s'      => \$dbport,
            'newdbhost=s'   => \$newdbhost,
            'newdbname=s'   => \$newdbname,
            'newdbuser=s'   => \$newdbuser,
            'newdbpass=s'   => \$newdbpass,
            'newdbport=s'   => \$newdbport,
            'vegadbname=s'  => \$vegadbname,
            'vegadbhost=s'  => \$vegadbhost,
            'vegadbport=s'  => \$vegadbport,
            'vegadbuser=s'  => \$vegadbuser,
            'set1_name:s'   => \$set1_name,
            'set2_name:s'   => \$set2_name,
            'biotypes:s'    => \@biotypes,
            'coordsystem|cs_name=s' => \$coordsystem,
            'path|cs_version=s'        => \$path,
            'non_ref!'      => \$non_ref,
            'duplicate_sr!' => \$duplicate_sr,
            'write!'        => \$write,
            'gsi_to_copy:s' => \$gsi_to_copy,
            'sql_update:s'  => \$sql_update,
            'verbose!'      => \$verbose,
) or throw("Failed to get opts");

my $db =
  Bio::EnsEMBL::DBSQL::DBAdaptor->new( -dbname => $dbname,
                                       -user   => $dbuser,
                                       -pass   => $dbpass,
                                       -host   => $dbhost,
                                       -port   => $dbport, );

my $new_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                              -dbname => $newdbname,
                                              -user => $newdbuser,
                                              -pass => $newdbpass,
                                              -host => $newdbhost,
                                              -port => $newdbport,
                                             );

my $vega_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                              -dbname => $vegadbname,
                                              -user => $vegadbuser,
                                              -pass => $vegadbpass,
                                              -host => $vegadbhost,
                                              -port => $vegadbport,
                                             );

if ($gsi_to_copy && $write) {
  die print "You should only use either the -gsi_to_copy flag for a dry run, or -write flag if ".
            "you're sure about writing data into your patched gene set DB, but not both flags.\n";
}

if ($sql_update && $write) {
  die print "You should only use either -sql_update for a dry run or -write flag if you ".
  "are sure about writing data into your patched gene set DB.\n";
}

if ($gsi_to_copy) {
  open (OUT, ">$gsi_to_copy") || die print "Cannot open output file to store gene stable IDs.\n";
}

if ($sql_update) {
  open (OUT2, ">$sql_update") || die print "Cannot open output file to store sql update statements.\n";
}

my $gene_adaptor = $new_db->get_GeneAdaptor;
my $trans_adaptor = $new_db->get_TranscriptAdaptor;

if ($write) {
  print "\nWATCH OUT! lincRNAs transfer and analysis logic_name "
      . "update for some genes will be performed directly on "
      . "the 2nd set (patched set) core DB.\n";
} else {
  throw("This script only works with -write option.");
}

# Fetch all the slices required. Can fetch non_reference (e.g. halpotypic
# regions) as well as duplicated regions (e.g. PAR of Y-chromosome) 
# depending on commandline flags

print "\nFetching slices from current core DB.\n" if ($verbose);

if ($non_ref) {
  print "Slices fetched will include non-reference sequences.\n" if ($verbose);
}

if ($duplicate_sr) {
  print "Slices fetched will include duplicated regions (e.g. PAR on Y-chromosome).\n" 
  if ($verbose);
}

my $slices =
  $db->get_SliceAdaptor->fetch_all( $coordsystem, $path, $non_ref,
                                    $duplicate_sr );

my $total_lincRNAs                     = 0;
my $lincRNA_cluster_nr                 = 0;
my $lincRNAs_to_copy_counter           = 0;
my $lincRNA_overlap_coding_counter     = 0;
my $lincRNA_overlap_proc_trans_counter = 0;
my $removed_genes_nr                   = 0;

if ( scalar(@biotypes) ) {
  if ($verbose) {
    print "\nWill fetch a subset of set2 (patched) genes by biotype.\n"
  }
  @biotypes = split( /,/, join( ',', @biotypes ) );
  print "\tBiotypes: " if ($verbose);
  foreach my $biotype (@biotypes) {
    print $biotype. "  " if ($verbose);
  }
  print "\n" if ($verbose);
}

# For each slice, get all lincRNAs from existing core DB and required merged
# genes from the new core DB.  If the -biotype flag was used on the command
# line, the script assumes that the merged genes need to be filtered by
# biotype. If -biotype flag was not used, the script assumes that we want
# to compare *all* merged genes with lincRNAs.

SLICE: foreach my $slice (@$slices) {

  print "\n***** Working on slice: " . $slice->name . "\n";

  my @existing_lincRNAs = @{ $slice->get_all_Genes_by_type('lincRNA') };

  my @merged_genes_for_comparison;

  if ( scalar(@biotypes) ) {
    foreach my $biotype (@biotypes) {
      push( @merged_genes_for_comparison,
            @{ $slice->get_all_Genes_by_type($biotype) } );
    }
  } else {
    @merged_genes_for_comparison =
      @{ $new_db->get_SliceAdaptor->fetch_by_name( $slice->name )->get_all_Genes() };
  }

  print "\t" . scalar(@existing_lincRNAs) . " $set1_name genes found\n";
  $total_lincRNAs += scalar(@existing_lincRNAs);
  print "\t"
    . scalar(@merged_genes_for_comparison)
    . " $set2_name genes found\n";

  if ( scalar(@existing_lincRNAs) == 0 ) {
    print "\tThis slice contains no lincRNA genes at all, skip slice.\n"
      if ($verbose);
    next SLICE;
  }

  my ( $clustered, $unclustered ) =
    @{
    simple_cluster_Genes( \@existing_lincRNAs,           $set1_name,
                          \@merged_genes_for_comparison, $set2_name
    ) };

  # First we look at single-set clusters which contain multiple lincRNAs
  # exclusively, without any merged genes. We want to copy each of the
  # clustered lincRNAs to the new DB if they have Ensembl supporting evidence.

  print "\n\tLooking at lincRNA-only clusters and unclustered lincRNA genes.\n" if ($verbose);

  my $single_set_clusters         = get_single_clusters($clustered);
  my $lincRNA_cluster_on_slice_nr = 0;

  foreach my $single_set_clust (@$single_set_clusters) {
    my @sets = @{ $single_set_clust->get_sets_included() };
    for my $set (@sets) {
      if ( $set eq $set1_name ) {
        $lincRNA_cluster_on_slice_nr++;
        $lincRNA_cluster_nr++;
      }
    }
  }

  # all unclustered genes are stored in GeneCluster objects too (1
  # cluster 1 gene).
  my @orphan_clusters_from_both_sets = @$unclustered;
  my @orphan_lincRNA_genes;
  foreach my $orphan_cluster (@orphan_clusters_from_both_sets) {
    push( @orphan_lincRNA_genes,
          @{ $orphan_cluster->get_Genes_by_Set($set1_name) } );
  }
  print "\t\t"
    . scalar(@orphan_lincRNA_genes)
    . " unclustered orphan lincRNA genes found.\n";

  my @all_single_set_clust_lincRNA_genes;
  if ( $lincRNA_cluster_on_slice_nr >= 1 ) {
    @all_single_set_clust_lincRNA_genes = @{ get_oneway_clustering_genes_of_set( $clustered, $set1_name ) };
    print "\t\t"
      . scalar(@all_single_set_clust_lincRNA_genes)
      . " single-set clustered lincRNA genes found in $lincRNA_cluster_nr clusters.\n";
  } elsif ( $lincRNA_cluster_on_slice_nr == 0 ) {
    print "\t\tThere are no single-set clustered lincRNA genes "
        . "found on this slice.\n";
  }

  my @lincRNA_genes;
  if (@all_single_set_clust_lincRNA_genes) {
    push(@lincRNA_genes,@all_single_set_clust_lincRNA_genes);
  }
  if (@orphan_lincRNA_genes) {
    push(@lincRNA_genes,@orphan_lincRNA_genes);
  }

  if (@lincRNA_genes) {
    foreach my $linc_RNA_gene (@lincRNA_genes) {

      # Need to check if the gene is still present in the new database.
      if (!check_in_new_db($linc_RNA_gene)) {
        $removed_genes_nr++;
        #next;
      }

      $linc_RNA_gene->load;
      if (has_ensembl_evidence($linc_RNA_gene)) {
        my $ensembl_only_gene = $linc_RNA_gene;
        if (has_havana_evidence($linc_RNA_gene)) {
          $ensembl_only_gene = delete_havana_evidence($linc_RNA_gene); # includes: analysis IDs for gene and transcript have been updated to ensembl_lincrna
        }
        $lincRNAs_to_copy_counter++;
        if ($write) {
          my $new_dbID = $gene_adaptor->store($ensembl_only_gene);
          print "Stored ensembl lincRNA gene with new dbID $new_dbID.\n";
        } else {
          print OUT $ensembl_only_gene->stable_id . "\n";
        }
      }
    }
  }

  # Second, we look at cases where lincRNAs overlap with merged genes.
  # If lincRNA overlaps with a protein_coding gene, we won't copy
  # the lincRNA.
  # Third, if lincRNA overlaps with a processed_transcript gene and
  # the lincRNA has Ensembl supporting evidence, we change
  # the biotype of the processed_transcript gene to "lincRNA" and
  # the logic_name to "ensembl_havana_lincrna". We also change the
  # logic_name of all underlying transcripts to "ensembl_havana_lincrna".

  print "\n\tLooking at clusters containing both lincRNA and merged genes.\n"
    if ($verbose);

  my $lincRNA_merged_clusters = get_twoway_clusters($clustered);
  print "\t\t"
    . scalar @$lincRNA_merged_clusters
    . " two_way clusters found.\n";
  if ( scalar(@$lincRNA_merged_clusters) ) {
    my $cluster_counter = 0;
  TWO_WAY_CLUSTER: foreach my $twoway_clust (@$lincRNA_merged_clusters) {
      $cluster_counter++;
      print "\t\tLooking at Cluster number $cluster_counter: "
        . $twoway_clust->start . " "
        . $twoway_clust->end . "\n";
      my @lincRNA_overlapped_with_merged =
        @{ $twoway_clust->get_Genes_by_Set($set1_name) };
      my @overlapped_merged_genes =
        @{ $twoway_clust->get_Genes_by_Set($set2_name) };

      foreach my $merged_gene (@overlapped_merged_genes) {
        if ( $merged_gene->biotype =~ /protein_coding/ ) {
          $lincRNA_overlap_coding_counter++;
          my $problem_lincRNA_gsi =
            $lincRNA_overlapped_with_merged[0]->stable_id;
          my $coding_merged_gene_gsi = $merged_gene->stable_id;
          print "\t\tlincRNA ($problem_lincRNA_gsi) overlaps with "
            . "a protein_coding merged gene ($coding_merged_gene_gsi). "
            . "All lincRNA(s) in this cluster will not be transferred.\n";
          next TWO_WAY_CLUSTER;
        } elsif ( $merged_gene->biotype =~ /processed_transcript/ and
                  has_ensembl_evidence($lincRNA_overlapped_with_merged[0]) and
                  $merged_gene->analysis->logic_name eq 'havana' # ignoring 'ensembl_havana_gene'
                ) {
          $lincRNA_overlap_proc_trans_counter++;
          print "\t\tOne or more lincRNAs with Ensembl supporting evidence in this cluster overlap "
              . "with a processed_transcript merged gene.\n";
          print "\t\tProc_trans gene stable ID: "
            . $merged_gene->stable_id
            . " , analysis logic_name is "
            . $merged_gene->analysis->logic_name . ".\n";

          my $analysis_adap = $new_db->get_AnalysisAdaptor();
          my $new_analysis_obj =
            $analysis_adap->fetch_by_logic_name('ensembl_havana_lincrna');
          my $new_analysis_id = $new_analysis_obj->dbID;

          my @merged_trans = @{ $merged_gene->get_all_Transcripts };

          if ($write) {
            $merged_gene->analysis($new_analysis_obj);
            $merged_gene->biotype("lincRNA");
            $gene_adaptor->update($merged_gene);

            foreach my $merged_trans (@merged_trans) {
              $merged_trans->analysis($new_analysis_obj);
              $trans_adaptor->update($merged_trans);
            }

            print "\t\tChanged logic_name of processed_transcript gene "
              . $merged_gene->stable_id
              . " to 'ensembl_havana_lincrna' and updated its biotype to lincRNA. ";
            print "Also updated logic_name of underlying transcripts.\n";

            #my $new_dbID = $gene_adaptor->store($merged_gene);
            #print "Stored ens-hav lincRNA gene with new dbID $new_dbID.\n";

          } else {
            print "\t\tNeed to change logic_name of proc_trans gene to "
              . $merged_gene->stable_id
              . " to 'ensembl_havana_lincrna'.\n";
            print OUT2 "update gene set analysis_id = $new_analysis_id "
                . "where gene_id = "
              . $merged_gene->dbID . ";\n";
            print OUT2 "update gene set biotype = 'lincRNA' where gene_id = "
              . $merged_gene->dbID . ";\n";
            foreach my $merged_tr (@merged_trans) {
              print OUT2 "update transcript "
                . "set analysis_id = $new_analysis_id "
                . "where transcript_id = " . $merged_tr->dbID . ";\n";
            }
          }
        } ## end elsif ( $merged_gene->biotype...
      } ## end foreach my $merged_gene (@overlapped_merged_genes)
    } ## end foreach my $twoway_clust (@$lincRNA_merged_clusters)
  } ## end if ( scalar(@$lincRNA_merged_clusters...
} ## end foreach my $slice (@$slices)

sub has_ensembl_evidence {
  my $gene = shift;
  my @exons = @{ $gene->get_all_Exons };
  my $ensembl_evidence = 0;

  EXON: foreach my $exon (@exons) {
    my @evidence = @{ $exon->get_all_supporting_features };
    foreach my $supporting_feature (@evidence) {
      my $analysis = $supporting_feature->analysis;
      if ($analysis->logic_name !~ /havana/) {
        $ensembl_evidence = 1;
        last EXON;
      }
    }
  }
  return $ensembl_evidence;
}

sub has_havana_evidence {
  my $gene = shift;
  my @exons = @{ $gene->get_all_Exons };
  my $havana_evidence = 0;

  EXON: foreach my $exon (@exons) {
    my @evidence = @{ $exon->get_all_supporting_features };
    foreach my $supporting_feature (@evidence) {
      my $analysis = $supporting_feature->analysis;
      if ($analysis->logic_name =~ /havana/) {
        $havana_evidence = 1;
        last EXON;
      }
    }
  }
  return $havana_evidence;
}

sub delete_havana_evidence {
# note that this sub only works for lincrna
  my $gene = shift;

  my $analysis_adaptor = $new_db->get_AnalysisAdaptor();
  my $new_analysis = $analysis_adaptor->fetch_by_logic_name('ensembl_lincrna');

  my @transcripts = @{ $gene->get_all_Transcripts };
  my @ensembl_transcripts;

  TRANSCRIPT: foreach my $transcript (@transcripts) {
    my @ensembl_exons;

    EXON: foreach my $exon (@{ $transcript->get_all_Exons }) {
      my @evidence = @{ $exon->get_all_supporting_features };
      my @ensembl_evidence;

      FEATURE: foreach my $supporting_feature (@evidence) {
        my $analysis = $supporting_feature->analysis;
        if ($analysis->logic_name !~ /havana/) {
          push(@ensembl_evidence,$supporting_feature);
        }
      } # end of feature

      $exon->flush_supporting_features;
      if (@ensembl_evidence) {
        $exon->add_supporting_features(@ensembl_evidence);
        push(@ensembl_exons,$exon);
      }
    } # end of exon

    $transcript->flush_Exons;
    if (@ensembl_exons) {
      foreach my $ensembl_exon (@ensembl_exons) {
        $transcript->add_Exon($ensembl_exon); # includes recalculate_coordinates
      }
      $transcript->analysis($new_analysis);
      push(@ensembl_transcripts,$transcript);
    }
  } # end of transcript

  $gene->{'_transcript_array'} = []; # flush transcripts
  if (@ensembl_transcripts) {
    foreach my $ensembl_transcript (@ensembl_transcripts) {
      $gene->add_Transcript($ensembl_transcript); # includes recalculate_coordinates
    }
  }
  $gene->analysis($new_analysis);
  return $gene;
}

sub check_in_new_db {
  my $gene = shift;

  my $gene_id = $gene->stable_id();
  my $new_gene = $gene_adaptor->fetch_by_stable_id($gene_id);
  print "checking if the gene exists in new db\n";
  if (!defined($new_gene)) {
    print "Gene $gene_id has been removed in the new merge database\n";
    return 0;
  }
  return 1;
}


sub check_in_vega {
  my $gene = shift;

  my $vega_ga = $vega_db->get_GeneAdaptor();

  if (    $gene->analysis->logic_name eq 'ensembl_havana_lincrna'
       || $gene->analysis->logic_name eq 'havana' )
  {
    print "The gene " . $gene->stable_id . " is a merged or havana one, "
      . "let's see if it's still in the Vega db\n";

    foreach my $vega_xref ( @{ $gene->get_all_DBEntries() } ) {
      if ( $vega_xref->primary_id =~ /^OTT/ ) {
        print "Found a vega xref\n";
        my $vega_gene =
          $vega_ga->fetch_by_stable_id( $vega_xref->primary_id );
        if ( !defined($vega_gene) ) {
          print "Couldn't find the gene, it's been removed.\n";
          return 0;
        }
      }
    }
  }
  return 1;
}

print "\n\n";
print "**************************** SUMMARY **********************************\n";
print "***********************************************************************\n";
print "Started with $total_lincRNAs lincRNAs.\n";
print "$removed_genes_nr lincRNA genes have been removed.\n";
print "$lincRNAs_to_copy_counter lincRNA genes are worth copying.\n";
print "$lincRNA_overlap_coding_counter lincRNA genes overlapped with protein_coding genes, discarded.\n\n";
print "$lincRNA_overlap_proc_trans_counter overlapped with processed_transcript (PT) genes.\n";
print "Remember to change PT genes' logic_name to 'ensembl_havana_lincrna' and biotype to 'lincRNA'.\n";
print "Transcript logic_name should be 'ensembl_havana_lincrna'. Do not change transcript biotype.\n";

close OUT;
