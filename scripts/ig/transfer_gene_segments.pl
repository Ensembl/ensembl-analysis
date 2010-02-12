#!/usr/local/ensembl/bin/perl

use strict;
use warnings;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Digest::MD5 qw(md5_hex);

my ( $dbname,      $dbhost,                 $dbuser,
     $dbport,      $dbpass,                 $tgdbname,
     $tgdbhost,    $tgdbuser,               $tgdbport,
     $tgdbpass,    @tg_biotypes,            %tg_biotypes,
     $tg_biotypes, $patch,                  $old_db_name,
     $new_db_name, $new_mapping_session_id, $verbose, );

$dbuser = 'ensro';
$tgdbuser = 'ensro';

$dbport = 3306;
$tgdbport = 3306;

my $sub_biotypes;

GetOptions( 'qydbname=s'       => \$dbname,
            'qydbuser=s'       => \$dbuser,
            'qydbhost=s'       => \$dbhost,
            'qydbport=s'       => \$dbport,
            'qydbpass=s'       => \$dbpass,
            'tgdbname=s'       => \$tgdbname,
            'tgdbuser=s'       => \$tgdbuser,
            'tgdbhost=s'       => \$tgdbhost,
            'tgdbport=s'       => \$tgdbport,
            'tgdbpass=s'       => \$tgdbpass,
            'tgbiotype=s'      => \$tg_biotypes,
            'patch'            => \$patch,
            'old_db_name=s'    => \$old_db_name,
            'new_db_name=s'    => \$new_db_name,
            'new_session_id=s' => \$new_mapping_session_id,
            'sub_biotypes'     => \$sub_biotypes,
            'verbose'          => \$verbose, );

if ($patch) {
  die "You must provide an old_db_name and a new_db_name"
      if not defined $old_db_name or not defined $new_db_name;
}

my $qy_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
        '-dbname' => $dbname,
        '-host' => $dbhost,
        '-user' => $dbuser,
        '-port' => $dbport,
        '-pass' => $dbpass
        );

my $tg_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
        '-dbname' => $tgdbname,
        '-host' => $tgdbhost,
        '-user' => $tgdbuser,
        '-port' => $tgdbport,
        '-pass' => $tgdbpass
        );

if ($tg_biotypes) {
  @tg_biotypes = split (",",$tg_biotypes);
	print "This is your array: ",join(" - ",@tg_biotypes),"\n";
  map { $tg_biotypes{$_} = 1 } @tg_biotypes;
} else {
  map { $tg_biotypes{$_} = 1 } ('protein_coding', 'pseudogene');
}




$verbose and print STDERR "Current target db gene summary:\n" . gene_stats_string . "\n";

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
# Assign stable ids to new genes and copy them over
#########################################################
my ($new_sid_hash, $pep_feat_hash);

$verbose and print STDERR "Fully loading genes...\n";
$pep_feat_hash = fully_load_genes(\@genes);

if ($patch) {
  $verbose and print STDERR "Setting up new stable ids...\n";
  set_stable_ids(\@genes);
}

########################################################
# Find relationships between new genes and old geneset
########################################################
$verbose and print STDERR "Comparing new genes to current...\n";
my ($genes_to_delete_hash, $stable_id_event_hash);

if ($sub_biotypes){
  ($genes_to_delete_hash, $stable_id_event_hash) = 
    compare_new_genes_with_current_genes(\@genes, 1);
}else{
  ($genes_to_delete_hash, $stable_id_event_hash) =
    compare_new_genes_with_current_genes(\@genes, 0);
}
######################################################################
# remove old genes and store new ones
######################################################################
$verbose and print STDERR "Storing new genes...\n";
@genes = @{store_genes($pep_feat_hash, \@genes)};

$verbose and print "Removing interfering old genes...\n";
foreach my $g (values %$genes_to_delete_hash) {
  $tg_db->get_GeneAdaptor->remove($g);
}


if ($patch) {
  ######################################################################
  # update relevant stable id tables with relationships
  ######################################################################
  $verbose and print STDERR "Populating the stable_id_event table...\n";
  my $map_session_id = populate_stable_id_event([values %$genes_to_delete_hash],
                                                 \@genes,
                                                 $stable_id_event_hash);
  
  ######################################################################
  # for translations no longer present in the database, add to
  # peptide/gene archive
  ######################################################################
  $verbose and print STDERR "Populating the peptide/gene archive tables...\n";
  populate_gene_archive([values %$genes_to_delete_hash], 
                         \@genes,
                         $map_session_id);
}


$verbose and printf(STDERR "Done (%s removed, %d added)\n%s\n",
                    scalar(keys %$genes_to_delete_hash),
                    scalar(@genes),
                    gene_stats_string);


######################################################################

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

sub record_stable_ids {
  my $glist = shift;

  my %stable_ids;

  foreach my $g (@$glist) {
    $stable_ids{$g->stable_id} = 'gene';
    foreach my $t (@{$g->get_all_Transcripts}) {
      $stable_ids{$t->stable_id} = 'transcript';
      if ($t->translation) {
        $stable_ids{$t->translation->stable_id} = 'translation';
      }
    }
  }

  return \%stable_ids;
}

######################################################################

sub set_stable_ids {
  my $glist = shift;

  sub increment_stable_id {
    my $sid = shift;
    
    my ($prefix,$suffix) = $sid =~ /([a-zA-Z]+)([0-9]+)/;
    return sprintf "%s%011d", $prefix, $suffix+1;
  }
  
  my $tm = time;

  my (%sids);

  my $archive_sql = "select type, max(old_stable_id) from stable_id_event group by type";

  my $st = $tg_db->dbc->prepare($archive_sql);
  $st->execute;
  while(my ($type, $max) = $st->fetchrow_array) {
    $sids{$type} = $max;
  }
  $st->finish;

  foreach my $tb ('gene', 'transcript', 'translation', 'exon') {
    my $sql = "select max(stable_id) from " . $tb . "_stable_id";
    $st = $tg_db->dbc->prepare($sql);
    $st->execute;
    while(my ($max) = $st->fetchrow_array) {
      if (not exists $sids{$tb} or ($max gt $sids{$tb})) {
        $sids{$tb} = $max;
      }
    }
    $st->finish;

    $sids{$tb} = increment_stable_id($sids{$tb});
  }

  foreach my $g (@$glist) {
    my $gid = $sids{gene}; $sids{gene} = increment_stable_id($gid);
    $g->stable_id($gid); 
    $g->created_date($tm);
    $g->modified_date($tm);
    $g->version(1);

    my %exon_ids; 
    foreach my $e (@{$g->get_all_Exons}) {
      my $eid = $sids{exon}; $sids{exon} = increment_stable_id($eid);
      $exon_ids{$e->dbID} = $eid;      
    }
    foreach my $t (@{$g->get_all_Transcripts}) {
      my $tid = $sids{transcript}; $sids{transcript} = increment_stable_id($tid);
      $t->stable_id($tid);
      $t->created_date($tm);
      $t->modified_date($tm);
      $t->version(1);
      foreach my $e (@{$t->get_all_Exons}) {
        $e->stable_id($exon_ids{$e->dbID});
        $e->created_date($tm);
        $e->modified_date($tm);
        $e->version(1);
      }
      if ($t->translation) {
        my $tr = $t->translation;
        my $trid = $sids{translation}; $sids{translation} = increment_stable_id($trid);
        $tr->stable_id($trid);
        $tr->created_date($tm);
        $tr->modified_date($tm);
        $tr->version(1);
      }
    }
  }
}

######################################################################

sub store_genes {
  my ($prot_feat_hash, $glist) = @_;

  #foreach my $g (@g) {
  #  printf "GENE %d = %s\n", $g->dbID, $g->stable_id;
  #  foreach my $t (@{$g->get_all_Transcripts}) {
  #    printf " TRANSCRIPT = %d = %s\n", $t->dbID, $t->stable_id;
  #    foreach my $e (@{$t->get_all_Exons}) {
  #      printf "  EXON %d = %s\n", $e->dbID, $e->stable_id;
  #    }
  #    if ($t->translation) {
  #      printf " TRANSLATION %d = %s\n", $t->translation->dbID, $t->translation->stable_id;
  #    }
  #  }
  #}

  my $g_adap = $tg_db->get_GeneAdaptor;
  my $p_adap = $tg_db->get_ProteinFeatureAdaptor;

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
      my $oslice = $tg_db->get_SliceAdaptor->fetch_by_region('toplevel',
                                                             $sr_id,
                                                             $g->start,
                                                             $g->end);
      
      my @ogenes = map { $_->transfer($tg_tl_slice) } @{$oslice->get_all_Genes};
      @ogenes = grep { exists($tg_biotypes{$_->biotype}) } @ogenes;
      #print "YOUR HAVE THIS ORIGINAL GTENES: ",scalar(@ogenes),"\n";
      if (@ogenes) {
        foreach my $t (@{$g->get_all_Transcripts}) {
          my @exons = @{$g->get_all_Exons};
          
          foreach my $og (@ogenes) {
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
                #   map $ot->translation to $t->translationi
								if ($sub_biotypes == 1){
								  if ($og->biotype =~/IG/){
                    $g->biotype($og->biotype);
								  }
								}
                if ($patch) {
                  $genes_to_remove{$og->stable_id} = $og;
                  $stable_id_event{gene}->{$og->stable_id}->{$g->stable_id} = 1;
                  $stable_id_event{transcript}->{$ot->stable_id}->{$t->stable_id} = 1;
                  if ($t->translation and $ot->translation) {
                    $stable_id_event{translation}->{$ot->translation->stable_id}->{$t->translation->stable_id} = 1;
                  }
                } else {
                  $genes_to_remove{$og->dbID} = $og;
                }
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

sub populate_stable_id_event {
  my ($del_genes, $new_genes, $stable_id_event) = @_;

  my ($prev_session_id, $this_session_id) = create_new_mapping_session;

  # add the deleted genes/transcripts/translations to the peptide/gene archive
  # add (NULL, new id) entries to stable_id_event for each new gene
  # add (id, id) entries to stable_event for each unaffected gene
  # add (deleted_id, new_id) entries to stable_id_event for old2new relationships

  my $deleted_id_hash = record_stable_ids($del_genes);
  my $new_id_hash = record_stable_ids($new_genes);


  my $sql = "SELECT new_stable_id, new_version, type ";
  $sql .= "FROM stable_id_event ";
  $sql .= "WHERE mapping_session_id = $prev_session_id ";
  $sql .= "AND new_stable_id IS NOT NULL";

  my %current;

  my $st = $tg_db->dbc->prepare($sql);
  $st->execute;
  while (my ($cur_id, $cur_ver, $type) = $st->fetchrow_array) {
    if (not exists $current{$cur_id}) {
      $current{$cur_id} = {
        type => $type,
        version => $cur_ver,
      };
    } else {
      if ($cur_ver > $current{$cur_id}->{version}) {
        $current{$cur_id}->{version} = $cur_ver;
      }
    }
  }
  $st->finish;

  $sql = "INSERT into stable_id_event VALUES(?,?,?,?,?,?,?)";
  $st = $tg_db->dbc->prepare($sql);

  #
  # first, catalog which new ids have been mapped to
  #
  my %mapped_new;
  foreach my $tp (keys %$stable_id_event) {
    foreach my $oid (keys %{$stable_id_event->{$tp}}) {
      foreach my $nid (keys %{$stable_id_event->{$tp}->{$oid}}) {
        $mapped_new{$nid} = 1;
      }
    }
  }

  #
  # now add (null, id) entries for new ids that do not map to old
  #
  foreach my $id (keys %$new_id_hash) {
    if (not exists $mapped_new{$id}) {
      my $tp = $new_id_hash->{$id};
      $st->execute(undef, 0, $id, 1, $this_session_id, $tp, 0);
    }
  }

  #
  # now add records for the mapped ids, existing and new
  #
  foreach my $id (keys %current) {
    my ($tp, $ver) = ($current{$id}->{type}, $current{$id}->{version});

    # Two cases
    # 1. the id has not been deleted: add entry ($id, $id)
    # 2. the gene has been deleted:
    #   - if mapped, add entry ($id, null) to represent deletion
    #   - otherwise, add entry ($id, $other)

    if (not exists $deleted_id_hash->{$id}) {
      $st->execute($id, $ver, $id, $ver, $this_session_id, $tp, 1);
    } else {
      if (exists $stable_id_event->{$tp}->{$id}) {
        foreach my $nid (keys %{$stable_id_event->{$tp}->{$id}}) {
          $st->execute($id, $ver, $nid, 1, $this_session_id, $tp, 1);
        }
      } else {
        $st->execute($id, $ver, undef, 0, $this_session_id, $tp, 0);
      }
    }
  }
  $st->finish;

  return $this_session_id;
}

######################################################################

sub populate_gene_archive {
  my ($deleted_genes, $new_genes, $mapping_session_id) = @_;

  my %new_digests;

  foreach my $g (@$new_genes) {
    foreach my $t (@{$g->get_all_Transcripts}) {
      if ($t->translation) {
        $new_digests{md5_hex($t->translate->seq)} = 1;
      }
    }
  }
  
  my %archive;

  foreach my $g (@$deleted_genes) {
    foreach my $t (@{$g->get_all_Transcripts}) {
      if ($t->translation) {
        my $md5 = md5_hex($t->translate->seq);
        if (not exists $new_digests{$md5}) {
          if (not exists $archive{$md5}) {
            $archive{$md5}->{peptide} = $t->translate->seq;
          }
          push @{$archive{$md5}->{gene_archive}}, {
            gid => $g->stable_id,
            gid_ver => $g->version,
            tid => $t->stable_id,
            tid_ver => $t->version,
            trid => $t->translation->stable_id,
            trid_ver => $t->translation->version,
          };
        }
      }
    }
  }

  my $sql1 = "SELECT peptide_archive_id FROM peptide_archive WHERE md5_checksum = ?";
  my $sql2 = "INSERT into peptide_archive(md5_checksum, peptide_seq) VALUES(?,?)";

  my $st1 = $tg_db->dbc->prepare($sql1);
  my $st2 = $tg_db->dbc->prepare($sql2);
  foreach my $md5 (keys %archive) {
    $st1->execute($md5);
    while(my ($pid) = $st1->fetchrow_array) {
      $archive{$md5}->{peptide_archive_id} = $pid;
    }
    if (not exists $archive{$md5}->{peptide_archive_id}) {
      $st2->execute($md5, $archive{$md5}->{peptide});
      $archive{$md5}->{peptide_archive_id} = $st2->{mysql_insertid};
    }
  }
  $st1->finish;
  $st2->finish;

  $sql1 = "INSERT into gene_archive VALUES(?,?,?,?,?,?,?,?)";
  $st1 = $tg_db->dbc->prepare($sql1);
  foreach my $md5 (keys %archive) {
    foreach my $en (@{$archive{$md5}->{gene_archive}}) {
      $st1->execute($en->{gid},
                   $en->{gid_ver},
                   $en->{tid},
                   $en->{tid_ver},
                   $en->{trid},
                   $en->{trid_ver},
                   $archive{$md5}->{peptide_archive_id},
                   $mapping_session_id)
    }
  }
  $st1->finish;
}

######################################################################

sub create_new_mapping_session {

  my ($last_session_id, 
      $last_session_rel, 
      $last_session_ass,
      $new_session_id);

  my $sql = "select mapping_session_id, new_release, new_assembly from mapping_session";
  my $st = $tg_db->dbc->prepare($sql);
  $st->execute;
  while (my ($id, $rel_num, $ass) = $st->fetchrow_array) {
    if ($rel_num =~ /^\s*(\d+)\s*$/) {
      my $rel = $1;
      if (not defined $last_session_rel or
          $rel > $last_session_rel) {
        $last_session_id = $id;
        $last_session_rel = $rel;
        $last_session_ass = $ass;
      }
    }
  }
  $st->finish;

  if (defined $new_mapping_session_id) {
    $sql = "insert into mapping_session ";
    $sql .= " values(?,?,?,?,?,?,?, now())";
    $st = $tg_db->dbc->prepare($sql);
    $st->execute($new_mapping_session_id,
                 $old_db_name, 
                 $new_db_name, 
                 $last_session_rel,
                 $last_session_rel + 1,
                 $last_session_ass,
                 $last_session_ass);
    $new_session_id = $new_mapping_session_id;
  } else {
    $sql = "insert into mapping_session ";
    $sql .= "(old_db_name, new_db_name, old_release, new_release, old_assembly, new_assembly, created)";
    $sql .= " values(?,?,?,?,?,? now())";
    $st = $tg_db->dbc->prepare($sql);
    $st->execute($old_db_name, 
                 $new_db_name, 
                 $last_session_rel,
                 $last_session_rel + 1,
                 $last_session_ass,
                 $last_session_ass);
    $new_session_id = $st->{mysql_insertid}; 
  }    

  return ($last_session_id, $new_session_id);
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
