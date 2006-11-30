#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::IgSegBuilder

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::RunnableDB::IgSegBuilder;

use vars qw(@ISA);
use strict;

# Object preamble

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::IgSegBuilder;

use Bio::EnsEMBL::Gene;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

############################################################

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);    
  
  $self->read_and_check_config($IGSEG_CONFIG_BY_LOGIC);
           
  # other stuff here
  
  return $self;
}


############################################################

sub fetch_input {
  my($self) = @_;
  
  my $trandb = $self->get_dbadaptor($self->TRANDB_DATABASES_NAME);
  eval {
    my $dnadb = $self->get_dbadaptor($self->DNADB_DATABASES_NAME);
    $trandb->dnadb($dnadb);
  };

  my $slice = $trandb->get_SliceAdaptor->fetch_by_name($self->input_id);
  my $tlslice = $trandb->get_SliceAdaptor->fetch_by_region('toplevel',
                                                           $slice->seq_region_name);
  $self->query($tlslice);

  my (@lv, @d, @j, @c);

  if ($self->LV_LOGICS) {
    foreach my $logic (@{$self->LV_LOGICS}) {
      foreach my $t (@{$slice->get_all_Transcripts(1, $logic)}) {
        push @lv, $t->transfer($tlslice);
      }
    }
  }
  if ($self->D_LOGICS) {
    foreach my $logic (@{$self->D_LOGICS}) {
      foreach my $t (@{$slice->get_all_Transcripts(1, $logic)}) {
        push @d, $t->transfer($tlslice);
      }
    }
  }
  if ($self->J_LOGICS) {
    foreach my $logic (@{$self->J_LOGICS}) {
      foreach my $t (@{$slice->get_all_Transcripts(1, $logic)}) {
        push @j, $t->transfer($tlslice); 
      }
    }
  }
  if ($self->C_LOGICS) {
    foreach my $logic (@{$self->C_LOGICS}) {
      foreach my $t (@{$slice->get_all_Transcripts(1, $logic)}) {
        push @c, $t->transfer($tlslice);
      }
    }
  }

  $self->LV_transcripts(\@lv);
  $self->D_transcripts(\@d);
  $self->J_transcripts(\@j);
  $self->C_transcripts(\@c);
}

#############################################################

sub write_output {
  my ($self) = @_;

  my $out_db = $self->get_dbadaptor($self->OUTPUT_DB_NAME);
  my $ga = $out_db->get_GeneAdaptor;

  foreach my $g (@{$self->output}) {
    $ga->store($g);
  }
}

#############################################################

sub run {
  my ($self) = @_;

  my @all_genes;

  # LV segments...

  my @lv_genes = @{$self->cluster_by_genomic_overlap($self->LV_transcripts)};
  @lv_genes = @{$self->transfer_to_local_slices(\@lv_genes)};
  foreach my $g (@lv_genes) {
    $g = $self->adjust_gene_3($g, "CAC");
    $g = $self->prune_LV_transcripts($g);

    # extend/trim translation to exclude stops
    push @all_genes, $g;
  }


  # J segments...

  my @j_genes = @{$self->cluster_by_genomic_overlap($self->J_transcripts)};
  @j_genes = @{$self->transfer_to_local_slices(\@j_genes)};
  foreach my $g (@j_genes) {
    $g = $self->adjust_gene_3($g, "GT");
    $g = $self->adjust_gene_5($g, "GTG");
    $g = $self->prune_D_J_transcripts($g);
    # extend/trim translation to exclude stops    

    push @all_genes, $g;
  }

  # C segments...

  my @c_genes = @{$self->cluster_by_genomic_overlap($self->C_transcripts)};
  @c_genes = @{$self->transfer_to_local_slices(\@c_genes)};
  foreach my $g (@c_genes) {
    $g = $self->adjust_gene_5($g, "AG");
    $g = $self->prune_C_transcripts($g);

    push @all_genes, $g;
  }

  # D segments. Don't do anything for now...

  my @d_genes = @{$self->cluster_by_genomic_overlap($self->D_transcripts)};
  @d_genes = @{$self->transfer_to_local_slices(\@d_genes)};
  foreach my $g (@d_genes) {
    $g = $self->adjust_gene_5($g, "GTG");
    $g = $self->adjust_gene_3($g, "CAC");;
    $g = $self->prune_D_J_transcripts($g);

    # extend/trim translation to exclude stops    
    push @all_genes, $g;
  }

  my @final_genes;
  foreach my $g (@all_genes) {
    $g = $g->transfer($self->query);
    $g->analysis($self->analysis);
    $g->biotype($self->OUTPUT_BIOTYPE);
    $self->prune_Exons($g);

    push @final_genes, $g;
  }

  $self->output(\@final_genes);
}

############################################################
# helper methods
############################################################

sub cluster_by_genomic_overlap {
  my ($self, $tran) = @_;

  my %tran_groups;
  foreach my $t (@$tran) {
    my $sr = $t->slice->seq_region_name;
    my $str = $t->strand;

    push @{$tran_groups{$sr . ":" . $str}}, $t;
  }

  my @genes;

  foreach my $grp (values %tran_groups) {
    my @c;
    foreach my $t (sort {$a->start <=> $b->start} @$grp) {
      if (not @c or $t->start > $c[-1]->end + 1) {
        push @c, Bio::EnsEMBL::Gene->new();
      }
      $c[-1]->add_Transcript($t);
    }
    push @genes, @c;
  }

  return \@genes;
}


########################################################

# the following method transfers the genes to a set
# of smaller slices; some slices are reverse strand
# to ensure that all gene and sequence is forward
# strand (makes motif checking easier)
#
sub transfer_to_local_slices {
  my ($self, $genes) = @_;
  
  my @groups;
  foreach my $g (sort { $a->start <=> $b->start } @$genes) {
    if (not @groups or 
        $groups[-1]->{strand} != $g->strand or
        $g->start > $groups[-1]->{end} + 100000) {
      push @groups, {
        start => $g->start,
        end   => $g->end, 
        strand => $g->strand,
        genes => [],
      };
    } elsif ($g->end > $groups[-1]->{end}) {
      $groups[-1]->{end} = $g->end;
    }
    push @{$groups[-1]->{genes}}, $g;
  }

  my @retgenes;
  foreach my $grp (@groups) {
    my $gstart = $grp->{start} - 10; 
    $gstart = 1 if $gstart < 1;
    my $gend   = $grp->{end} + 10; 
    $gstart = $self->query->end if $gstart > $self->query->end;
    
    my $local_slice = $self->get_dbadaptor($self->DNADB_DATABASES_NAME)->
        get_SliceAdaptor->fetch_by_region($self->query->coord_system->name,
                                          $self->query->seq_region_name,
                                          $gstart,
                                          $gend,
                                          $grp->{strand});
    print "TRANSFERRINF TO LOCAL SLICE: ", $local_slice->start, " ", $local_slice->end, "\n";
    push @retgenes, map { $_->transfer($local_slice) } @{$grp->{genes}}; 
  }

  return \@retgenes;
}

#############################################################
sub prune_Exons {
  my ($self, $gene) = @_;

  my %gene_exons;

  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @t_exons;
    foreach my $e (@{$tran->get_all_Exons}) {
      my $k = join(":",
                   $e->start,
                   $e->end,
                   $e->strand,
                   $e->phase,
                   $e->end_phase);
      if (not exists $gene_exons{$k}) {
        $gene_exons{$k} = $e;
      }
      push @t_exons, $gene_exons{$k};

      if ($e == $tran->translation->start_Exon) {
        $tran->translation->start_Exon($gene_exons{$k});
      }
      if ($e == $tran->translation->end_Exon) {
        $tran->translation->end_Exon($gene_exons{$k});
      }
    }
    
    $tran->flush_Exons;
    map { $tran->add_Exon($_) } @t_exons;
  }
  return $gene;
}


#############################################################

sub prune_LV_transcripts {
  my ($self, $gene) = @_;

  my @ftran;  
  my @tran = @{$gene->get_all_Transcripts};

  #
  # if there are transcripts with a single intron, remove others
  #
  foreach my $t (@tran) {
    if (scalar(@{$t->get_all_Introns}) == 1) {
      push @ftran, $t;
    }
  }
  if (@ftran) {
    @tran = @ftran; 
    @ftran = ();
  }
  
  #
  # if there are transcripts with all-consensus introns, remove others
  #
  my $seq = uc($gene->slice->seq);

  foreach my $t (@tran) {
    my $all_consensus = 1;
    foreach my $i (@{$t->get_all_Introns}) {
      my $left = substr($seq, $i->start - 1, 2);
      my $right = substr($seq, $i->end - 2, 2);
 
      if ($left ne 'GT' or $right ne 'AG') {
        $all_consensus = 0;
        last;
      }
    }
    push @ftran, $t if $all_consensus;
  }
  if (@ftran) {
    @tran = @ftran;
    @ftran = ();
  }

  #
  # if there are transcripts that begin with ATG, remove others
  #
  foreach my $t (@tran) {
    my $start = substr($seq, $t->start - 1, 3);
    push @ftran, $t if $start eq 'ATG';
  }
  if (@ftran) {
    @tran = @ftran;
    @ftran = ();
  }

  #
  # finally, remove redundant transcripts; only those with
  # unique introns are kept
  #
  @tran = sort _by_total_exon_length @tran;

  foreach my $t (@tran) {
    my @intr = @{$t->get_all_Introns};
    my $redundant = 0;
    foreach my $kt (@ftran) {
      my @kintr = @{$kt->get_all_Introns};
      my $matches = 0;
      foreach my $i (@intr) {
        foreach my $ki (@kintr) {
          if ($i->start == $ki->start and
              $i->end   == $ki->end) {
            $matches++;
            last;
          }
        }
      }
      if ($matches == scalar(@intr)) {
        $redundant = 1;
        last;
      }
    }
    if (not $redundant) {
      push @ftran, $t;
    }
  }

  return Bio::EnsEMBL::Gene->new(-transcripts => \@ftran);
}

########################################################

sub prune_C_transcripts {
  my ($self, $gene) = @_;

  my @tran = @{$gene->get_all_Transcripts};
  @tran = sort _by_total_exon_length @tran;

  my (@newtran, @intron_lists, %all_exons);
  foreach my $tran (@tran) {
    my @exons = sort {$a->start <=> $b->start} @{$tran->get_all_Exons};
    my @introns = sort {$a->start <=> $b->start} @{$tran->get_all_Introns};
    my $is_redundant = 0;
    foreach my $il (@intron_lists) {
      my $all_found = 1;
      foreach my $ti (@introns) {
        my $found_in_other = 0;
        foreach my $oi (@$il) {
          if ($ti->start == $oi->start and
              $ti->end   == $oi->end) {
            $found_in_other = 1;
            last;
          }
        }
        if (not $found_in_other) {
          $all_found = 0;
        }        
      }
      if ($all_found) {
        $is_redundant = 1;
        last;
      }
    }
    if (not $is_redundant) {
      push @newtran, $tran;
      push @intron_lists, \@introns;
      map { $all_exons{$_->start . ":" . $_->end} = $_} @exons;
    } else {
      # all introns have been seen in another transcript,
      # but perhaps the terminal exons contribute something
      # unique?
      my $all_covered = 1;
      foreach my $tex ($exons[0], $exons[-1]) {
        my $covered = 0;
        foreach my $eid (keys %all_exons) {
          if ($all_exons{$eid}->start <= $tex->start and
              $all_exons{$eid}->end   >= $tex->end) {
            # the exon is completely covered by another
            $covered = 1;
            last;
          }
        }
        if (not $covered) {
          $all_covered = 0;
          last;
        }
      }
      if (not $all_covered) {
        push @newtran, $tran;
        push @intron_lists, \@introns;
        map { $all_exons{$_->start . ":" . $_->end} = $_} @exons;
      }
    }
  }

  return Bio::EnsEMBL::Gene->new(-transcripts => \@newtran);
}

########################################################

sub prune_D_J_transcripts {
  my ($self, $gene) = @_;

  my @tran = @{$gene->get_all_Transcripts};
  @tran = grep { scalar(@{$_->get_all_Exons}) == 1 } @tran;

  my ($best) = sort { $a->length <=> $b->length } @tran;

  return Bio::EnsEMBL::Gene->new(-transcripts => [$best]);
}

#########################################################

sub adjust_gene_5 {
  my ($self, $gene, $motif) = @_;

  my $seq = uc($gene->slice->seq);

  # we will consider adding/removing a maximum of
  # 2 bps to the transcript start.

  my (@adjusted_trans);
  
  foreach my $t (@{$gene->get_all_Transcripts}) {
    my $reg_end = $t->start + 1;
    my $reg_start = $t->start - 2 - length($motif);
    my $subseq = substr($seq, $reg_start - 1, $reg_end - $reg_start + 1);
    print "SUBSEQ = $subseq\n";

    # favour extension so take the first match
    if ((my $subseqoff = index($subseq, $motif)) >= 0) {
      my $newstart = $subseqoff + length($motif) + $reg_start;
      $self->adjust_transcript($t, $t->start - $newstart, 0);
    }
    push @adjusted_trans, $t;
  }

  return Bio::EnsEMBL::Gene->new(-transcripts => \@adjusted_trans);
}


sub adjust_gene_3 {
  my ($self, $gene, $motif) = @_;

  my $seq = uc($gene->slice->seq);

  # we will consider adding/removing a maximum of
  # 2 bps to the transcript end. 

  my (@adjusted_trans);
  
  foreach my $t (@{$gene->get_all_Transcripts}) {
    my $tname = $t->get_all_supporting_features->[0]->hseqname;

    my $reg_start = $t->end - 1;
    my $reg_end = $t->end + 2 + length($motif);
    my $subseq = substr($seq, $reg_start - 1, $reg_end - $reg_start + 1);
    
    # favour extension so take the last match
    if ((my $subseqoff = index($subseq, $motif)) >= 0) {
      while((my $anotheroff = index($subseq, $motif, $subseqoff + 1)) >= 0) {
        $subseqoff = $anotheroff;
      }
      my $newend = $subseqoff + $reg_start - 1;
      $self->adjust_transcript($t, 0, $newend - $t->end);
    }
    push @adjusted_trans, $t;
  }

  return Bio::EnsEMBL::Gene->new(-transcripts => \@adjusted_trans);

}

# offsets in the following are how much to EXTEND transcript
# by at each end (-> negative offsets are truncations)

sub adjust_transcript {
  my ($self, $tran, $offset5, $offset3) = @_;

  my @exons = @{$tran->get_all_Exons};

  if ($offset5) {
    $exons[0]->start($exons[0]->start - $offset5);
    $exons[0]->phase( ($exons[0]->phase - $offset5) % 3 );
    if ($tran->translation) {
      $tran->translation->end($tran->translation->end + $offset5);
    }
    $tran->recalculate_coordinates;
  }

  if ($offset3) {
    $exons[-1]->end($exons[-1]->end + $offset3);    
    $exons[-1]->end_phase( ($exons[-1]->end_phase + $offset3) % 3 );
    if ($tran->translation) {
      $tran->translation->end($tran->translation->end + $offset3);
    }
    $tran->recalculate_coordinates;
  }
}


#######

sub _by_total_exon_length {
  my $alen = 0;
  my $blen = 0;
  foreach my $e (@{$a->get_all_Exons}) {
    $alen += $e->length;
  }
  foreach my $e (@{$b->get_all_Exons}) {
    $blen += $e->length;
  }
  
  return $blen <=> $alen;
}



############################################################
# containers
############################################################

sub LV_transcripts{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_lv_transcripts} = $value;
  }
  return $self->{_lv_transcripts};
}


sub D_transcripts{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_d_transcripts} = $value;
  }
  return $self->{_d_transcripts};
}


sub J_transcripts{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_j_transcripts} = $value;
  }
  return $self->{_j_transcripts};
}


sub C_transcripts{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_c_transcripts} = $value;
  }
  return $self->{_c_transcripts};
}



############################################################
# config holders
############################################################

sub LV_LOGICS {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_lv_logic} = $value;
  }
  return $self->{_lv_logic};
}


sub J_LOGICS {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_j_logic} = $value;
  }
  return $self->{_j_logic};

}

sub D_LOGICS {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_d_logic} = $value;
  }
  return $self->{_d_logic};

}

sub C_LOGICS {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_c_logic} = $value;
  }
  return $self->{_c_logic};
}


sub DNADB_DATABASES_NAME {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_dnadb_name} = $val;
  }

  return $self->{_dnadb_name};
}


sub TRANDB_DATABASES_NAME {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_transcriptdb_name} = $val;
  }

  return $self->{_transcriptdb_name};
}


sub OUTPUTDB_DATABASES_NAME {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_outputdb_name} = $val;
  }

  return $self->{_outputdb_name};
}


sub OUTPUT_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_output_biotype} = $val;
  }

  return $self->{_output_biotype};
}



1;
