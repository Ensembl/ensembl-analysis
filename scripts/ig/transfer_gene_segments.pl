# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

transfer_ig_genes.pl

=head1 DESCRIPTION

This script transfers IG gene segments from a source ("query") DB to an output 
("target") DB. In the process, if the "query" gene segment overlaps at the exon
level with a gene (of a specified biotype) in the target DB, the "query" gene will
*overwrite* the "target" one. 

Note: the script expects the query DB to contain *no other gene types* except 
Ig gene segments.


=head1 OPTIONS

=head2 DB connection:

=head3 DB containing query (source) Ig gene segments

        -qydbuser        Read-only username to connect to query DB
                        
        -qydbhost        Where the query DB is

        -qydbname        query DB name

        -qydbport        Port to connect to query DB


=head3 DB containing target (output) genes (some of them to be overwritten by "query" genes)

        -tgdbuser        Username to connect to target DB with write permissions

        -tgdbpass        Password to connect to target DB

        -tgdbhost        Where the target DB is

        -tgdbname        target DB name
 
        -tgdbport        Port to connect to target DB


=head3 DB containing DNA sequences. You must provide such details if your query
       and/or target DB contains no DNA sequence.

        -dnadbuser       Read-only username to connect to DNA DB
 
        -dnadbhost       Where the DNA DB is

        -dnadbname       DNA DB name

        -dnadbport       Port to connect to DNA DB

        -path            The version of the assembly you want to fetch DNA sequences
                         from (e.g "GRCh37", "NCBIM37")

        -----------------------------------------------------------------------

=head2 Analysis and output options:

        -tgbiotypes     A list of biotypes of target genes which should be compared
                        against query genes *and* can be replaced by query genes if 
                        substantial exon overlap is found between the query and target. 

                        The list should be comma-separated with no whitespace, 
                        e.g."protein_coding,pseudogene"
 
        -sub_biotypes   Optional. A boolean flag to indicate whether a biotype 
                        substitution check needs to be done. Switched off by default. 

                        The check mainly concerns copying  Ig gene segments from 
                        Vega DB (query DB) to a target DB which already contains some
                        Ig gene segments from Ensembl. Vega Ig gene segments always 
                        have priority over Ensembl ones, so Vega Ig models will 
                        overwrite Ensembl ones. However, Ensembl Ig models 
                        sometimes have more detailed biotypes (e.g. "IG_V_gene"
                        instead of simply "IG_gene"), and we want to keep the more
                        elaborate biotypes.  Therefore, the script can check if the 
                        target gene has a biotype starting with "IG*", and if yes, 
                        the script will assign the original target's biotype onto
                        the overwriting Vega Ig model.
                        
        -verbose        A boolean option to get more print statements to follow
                        the script.  Set to 0 (not verbose) by default.


=cut


#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Digest::MD5 qw(md5_hex);

my ( $dbname,     $dbhost,      $dbuser,      $dbport, 
     $tgdbname,   $tgdbhost,    $tgdbuser,    $tgdbport,    $tgdbpass,
     $path,       @tg_biotypes, %tg_biotypes, $tg_biotypes,
     $dna_dbuser, $dna_dbhost, $dna_dbname, $dna_dbport );

$dbuser = 'ensro' ;
$dbport = 3306 ;
$tgdbuser = 'ensro' ;
$tgdbport = 3306 ;
$dna_dbuser = 'ensro' ;
$dna_dbport = 3306 ;
my $sub_biotypes = 0;
my $verbose = 0;

GetOptions( 'qydbname|dbname|db|D=s'   => \$dbname,
            'qydbuser|dbuser|user|u=s'   => \$dbuser,
            'qydbhost|dbhost|host|h=s'   => \$dbhost,
            'qydbport|dbport|port|P=s'   => \$dbport,
            'tgdbname=s'   => \$tgdbname,
            'tgdbuser=s'   => \$tgdbuser,
            'tgdbhost=s'   => \$tgdbhost,
            'tgdbport=s'   => \$tgdbport,
            'tgdbpass=s'   => \$tgdbpass,
            'tgbiotypes=s' => \$tg_biotypes,
            'dnadbuser:s'  => \$dna_dbuser,
            'dnadbhost:s'  => \$dna_dbhost,
            'dnadbname:s'  => \$dna_dbname,
            'dnadbport:s'  => \$dna_dbport,
            'sub_biotypes!'=> \$sub_biotypes,
            'verbose!'      => \$verbose,
            'path=s'       => \$path );

my $qy_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
        '-dbname' => $dbname,
        '-host' => $dbhost,
        '-user' => $dbuser,
        '-port' => $dbport,
        );

my $tg_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
        '-dbname' => $tgdbname,
        '-host' => $tgdbhost,
        '-user' => $tgdbuser,
        '-port' => $tgdbport,
        '-pass' => $tgdbpass
        );

my $qyDB_has_DNA = check_if_DB_contains_DNA($qy_db);
my $tgDB_has_DNA = check_if_DB_contains_DNA($tg_db);

if ($qyDB_has_DNA == 0 || $tgDB_has_DNA == 0) {
  if (!defined $dna_dbname || !defined $dna_dbhost) {
  throw ("You must provide both -dnadbname and -dnadbhost on the commandline to connect to ".
         "DNA_DB or else the code will die when trying to fully load genes.");
  }

  my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new( -dbname => $dna_dbname,
         -host   => $dna_dbhost,
         -user   => $dna_dbuser,
         -port   => $dna_dbport,
         -path   => $path
        );

  if ($qyDB_has_DNA == 0) {
    print "Attaching DNA_DB to " . $qy_db->dbc->dbname . "\n";
    $qy_db->dnadb($dnadb);
  }
  if ($tgDB_has_DNA ==0) {
    print "Attaching DNA_DB to " . $tg_db->dbc->dbname . "\n";
    $tg_db->dnadb($dnadb);
  }
}

print "\n";

if ($tg_biotypes) {
  @tg_biotypes = split (",",$tg_biotypes);
	print "This is your array of transcript biotypes which can be overwritten by Ig models should there be exon overlap : ",join(" - ",@tg_biotypes),"\n";
  map { $tg_biotypes{$_} = 1 } @tg_biotypes;
} else {
  map { $tg_biotypes{$_} = 1 } ('protein_coding', 'pseudogene');
}

$verbose and print STDERR "\nCurrent target db gene summary:\n" . &gene_stats_string . "\n";

my (@genes, %genes_by_slice);

$verbose and print STDERR "Fetching new genes...\n";

if (@ARGV) {
  foreach my $slid (@ARGV) {
    my $sl = $qy_db->get_SliceAdaptor->fetch_by_name($slid);
    my $tl_sl = $qy_db->get_SliceAdaptor->fetch_by_region('toplevel',
                                                       $sl->seq_region_name);
    my @genes = @{$sl->get_all_Genes};
    @genes = map { $_->transfer($tl_sl) } @genes;
    push @{$genes_by_slice{$sl->seq_region_name}}, @genes;
  }

} else {
  @genes = @{$qy_db->get_GeneAdaptor->fetch_all};
}

#########################################################
# Fully load query DB (source) Ig genes to keep their 
# protein annotation features  
#########################################################
my ($new_sid_hash, $pep_feat_hash);

$verbose and print STDERR "Fully loading genes...\n";
$pep_feat_hash = &fully_load_genes(\@genes);

########################################################
# Find relationships between new genes and old geneset
########################################################
$verbose and print STDERR "Comparing new genes to current...\n";
my ($genes_to_delete_hash, $stable_id_event_hash);

if ($sub_biotypes){
  ($genes_to_delete_hash, $stable_id_event_hash) = 
    &compare_new_genes_with_current_genes(\@genes, 1);
}else{
  ($genes_to_delete_hash, $stable_id_event_hash) =
    &compare_new_genes_with_current_genes(\@genes, 0);
}
########################################################
# remove old genes and store new ones
########################################################
$verbose and print STDERR "Storing new genes...\n";
@genes = @{&store_genes($pep_feat_hash, \@genes)};

$verbose and print "Removing interfering old genes...\n";
foreach my $g (values %$genes_to_delete_hash) {
  $tg_db->get_GeneAdaptor->remove($g);
}


$verbose and printf(STDERR "Done (%s removed, %d added)\n%s\n",
                    scalar(keys %$genes_to_delete_hash),
                    scalar(@genes),
                    &gene_stats_string);


##############################################################

sub fully_load_genes {
  my $glist = shift;
  
  my $prot_feat_hash = {};

  foreach my $g (@$glist) {
    foreach my $t (@{$g->get_all_Transcripts}) {
      foreach my $e (@{$t->get_all_Exons}) {
        $e->get_all_supporting_features;
      }
      my $tr = $t->translation;
      if (defined $tr) {
        # adaptor cascade for gene does not handle protein
        # features. In order to do that here, we have to keep
        # the the translation object so that we can obtain
        # its new dbID after storage
        $prot_feat_hash->{$tr->dbID} = [$tr];
        foreach my $pf (@{$tr->get_all_ProteinFeatures}) {
          push @{$prot_feat_hash->{$tr->dbID}}, $pf;
        }
      }
      $t->get_all_supporting_features;
    }
  }

  return $prot_feat_hash;
}

######################################################################

######################################################################

sub store_genes {
  my ($prot_feat_hash, $glist) = @_;

  my $g_adap = $tg_db->get_GeneAdaptor;
  my $p_adap = $tg_db->get_ProteinFeatureAdaptor;
  my $dbea   = $tg_db->get_DBEntryAdaptor;

  foreach my $g (@$glist) {
    $g_adap->store($g);
  }

  # At this point, all translations should have new dbIDs.
  # We can now store the protein features
  foreach my $old_dbid (keys %$prot_feat_hash) {
    my ($trn, @feats) = @{$prot_feat_hash->{$old_dbid}};

    my $new_dbid = $trn->dbID;
    foreach my $f (@feats) {
      $p_adap->store($f, $new_dbid);
    }
  }

  # finally, refetch the stored genes from the target database
  # so that they we are working exclusively with target databases
  # adaptors from here on in
  my @tg_genes;
  foreach my $g (@$glist) {
    my $ng = $g_adap->fetch_by_dbID($g->dbID);
    
    foreach my $t (@{$ng->get_all_Transcripts}) {
      $t->get_all_supporting_features;
      $t->translation;
      foreach my $e (@{$t->get_all_Exons}) {
        $e->get_all_supporting_features;
      }
    }
    push @tg_genes, $ng;
  }

  return \@tg_genes;
}

######################################################################

sub compare_new_genes_with_current_genes {
  my ($glist, $sub_biotypes) = @_;

  my (%by_slice, %genes_to_remove, %stable_id_event);
  
  map { push @{$by_slice{$_->slice->seq_region_name}}, $_ } @$glist;

  foreach my $sr_id (keys %by_slice) {
    my @g = sort { $a->start <=> $b->start } @{$by_slice{$sr_id}};
    
    my $tg_tl_slice = $tg_db->get_SliceAdaptor->fetch_by_region('toplevel',
                                                                 $sr_id);
     
    foreach my $g (@g) {
      # Need to fully load the genes
      $g->load();
      my $oslice = $tg_db->get_SliceAdaptor->fetch_by_region('toplevel',
                                                             $sr_id,
                                                             $g->start,
                                                             $g->end);

      my @ogenes = map { $_->transfer($tg_tl_slice) } @{$oslice->get_all_Genes};
      @ogenes = grep { exists($tg_biotypes{$_->biotype}) } @ogenes;
      # print "YOUR HAVE THIS ORIGINAL GENES: ",scalar(@ogenes),"\n"; 
      if (@ogenes) {
        foreach my $t (@{$g->get_all_Transcripts}) {
          my @exons = @{$g->get_all_Exons};
          
          foreach my $og (@ogenes) {
            # Fully load the genes
            $og->load();

            foreach my $ot (@{$og->get_all_Transcripts}) {
              my $has_exon_overlap = 0;
              
              PAIR: foreach my $oe (@{$ot->get_all_Exons}) {
                foreach my $e (@exons) {
                  if ($e->strand == $oe->strand and
                      $e->overlaps($oe)) {
                    $has_exon_overlap = 1;
                    last PAIR;
                  }
                }
              }
              
              if ($has_exon_overlap) {
                # transcripts $t and $ot have exon overlap
                # therefore:
                #   delete gene $og (and all transcripts)
                #   map $og to $g
                #   map $ot to $t
                #   map $ot->translation to $t->translation

                # Do a biotype check here if "sub_biotypes" flag
                # is set to 1.
                # If the original/target gene (to be replaced) has 
                # an "IG*" biotype, we keep the target biotype as 
                # it's already quite specific.
                # i.e. we do replace the target gene *model* with 
                # the query *model* but the replaced gene will
                # retain its original biotype.

                if ($sub_biotypes == 1){
                  if ($og->biotype =~/IG/){
                    $g->biotype($og->biotype);
                    $t->biotype($ot->biotype);
                  }
                }
                $genes_to_remove{$og->dbID} = $og;
              }
            }
          }
        }
      }
    }
  }

  return (\%genes_to_remove, \%stable_id_event);
}

######################################################################

sub gene_stats_string {
  my $summary_string = "";

  my $sql = "select biotype, count(*) from gene group by biotype order by biotype";
  my $st = $tg_db->dbc->prepare($sql);
  $st->execute;
  while(my ($tp, $c) = $st->fetchrow_array) {
    $summary_string .= sprintf("%-20s %5d\n", $tp, $c);
  }
  $st->finish;

  return $summary_string;
}

#######################################################################

sub check_if_DB_contains_DNA {
  my ($db) = @_;
  my $sql_command = "select count(*) from dna";
  my $sth = $db->dbc->prepare($sql_command);
  $sth->execute();
  my @dna_array = $sth->fetchrow_array;
  if ($dna_array[0] > 0) {
    print "Your DB ". $db->dbc->dbname ." contains DNA sequences. No need to attach a ". 
          "DNA_DB to it.\n" if ($verbose);
    return 1;
  } else { 
    print "Your DB ". $db->dbc->dbname ." does not contain DNA sequences.\n"
          if ($verbose);
    return 0;
  }
}

