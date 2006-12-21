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

use Bio::SeqIO;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

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

  my $slice = $trandb->get_SliceAdaptor->fetch_by_name($self->input_id);
  my $tlslice = $trandb->get_SliceAdaptor->fetch_by_region('toplevel',
                                                           $slice->seq_region_name);
  $self->query($tlslice);

  my (@lv, @d, @j, @c);

  my @groups = 
      (
       { logics => $self->LV_LOGICS, result => \@lv, },
       { logics => $self->D_LOGICS,  result => \@d,  },
       { logics => $self->J_LOGICS,  result => \@j,  },
       { logics => $self->C_LOGICS,  result => \@c,  },
       );

  foreach my $grp (@groups) {
    my ($logics, $result) = ($grp->{logics}, $grp->{result});

    if ($logics) {
      foreach my $logic (@$logics) {
        foreach my $t (@{$slice->get_all_Transcripts(1, $logic)}) {
          map { $_->get_all_supporting_features } ($t, @{$t->get_all_Exons});
          $t = $t->transfer($tlslice);


          if ($t->coding_region_start > $t->start or $t->coding_region_end < $t->end) {
            my @e = @{$t->get_all_translateable_Exons};
            my $tr = Bio::EnsEMBL::Translation->new(-start_Exon => $e[0],
                                                    -end_Exon => $e[-1],
                                                    -start => 1,
                                                    -end => $e[-1]->length);
            my @sfs = @{$t->get_all_supporting_features};
            $t = Bio::EnsEMBL::Transcript
                ->new(-analysis => $t->analysis,
                      -exons    => $t->get_all_translateable_Exons);
            $t->add_supporting_features(@sfs);
            $t->translation($tr);
            
            if ($t->get_all_Exons->[0]->phase < 0) {
              $t->get_all_Exons->[0]->phase(0);
            }
            if ($t->get_all_Exons->[-1]->end_phase < 0) {
              $t->get_all_Exons->[-1]->end_phase(0);
            }
          }
          push @{$result}, $t;
        }
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

  my $out_db = $self->get_dbadaptor($self->OUTPUTDB_DATABASES_NAME);
  my $ga = $out_db->get_GeneAdaptor;

  foreach my $g (@{$self->output}) {
    $ga->store($g);
  }
}

#############################################################

sub run {
  my ($self) = @_;

  my @genes;

  # LV segments...

  my @lv_genes = @{$self->cluster_by_genomic_overlap($self->LV_transcripts)};
  @lv_genes = @{$self->transfer_to_local_slices(\@lv_genes)};
  foreach my $g (@lv_genes) {
    $self->adjust_gene_3($g, "CAC");
    $self->trim_gene_translation3($g);
    $g = $self->prune_LV_transcripts($g);
    if (@{$g->get_all_Transcripts}) {
      $g = $g->transfer($self->query);
      $g->biotype($self->LV_OUTPUT_BIOTYPE);
      map { $_->biotype($self->LV_OUTPUT_BIOTYPE) } @{$g->get_all_Transcripts};
      push @genes, $g;
    }
  }
  @lv_genes = @genes; @genes = ();


  # C segments...

  my @c_genes = @{$self->cluster_by_genomic_overlap($self->C_transcripts)};
  @c_genes = @{$self->transfer_to_local_slices(\@c_genes)};
  foreach my $g (@c_genes) {
    $self->adjust_gene_5($g, "AG");
    $self->set_gene_stop_codon($g);
    $g = $self->prune_C_transcripts($g);
    if (@{$g->get_all_Transcripts}) {
      $g = $g->transfer($self->query);
      $g->biotype($self->C_OUTPUT_BIOTYPE);
      map { $_->biotype($self->C_OUTPUT_BIOTYPE) } @{$g->get_all_Transcripts};
      push @genes, $g;
    }
  }
  @c_genes = @genes; @genes = ();


  # D segments...

  my @d_genes = @{$self->cluster_by_genomic_overlap($self->D_transcripts)};
  @d_genes = grep { $self->gene_close_to_others($_, 2000000, \@lv_genes, \@c_genes) } @d_genes;
  @d_genes = @{$self->transfer_to_local_slices(\@d_genes)};

  foreach my $g (@d_genes) {
    $self->adjust_gene_5($g, "GTG");
    $self->adjust_gene_3($g, "CAC");;
    #$self->trim_gene_translation5($g);
    #$self->trim_gene_translation3($g);
    $g = $self->prune_D_J_transcripts($g);
    if (@{$g->get_all_Transcripts}) {
      $g = $g->transfer($self->query);
      $g->biotype($self->D_OUTPUT_BIOTYPE);
      map { $_->biotype($self->D_OUTPUT_BIOTYPE) } @{$g->get_all_Transcripts};
      push @genes, $g;
    }
  }
  @d_genes = @genes; @genes = ();


  # J segments...

  my @j_genes = @{$self->cluster_by_genomic_overlap($self->J_transcripts)};
  @j_genes = grep { $self->gene_close_to_others($_, 2000000, \@lv_genes, \@c_genes) } @j_genes;
  @j_genes = @{$self->transfer_to_local_slices(\@j_genes)};
  foreach my $g (@j_genes) {
    $self->adjust_gene_5($g, "GTG");
    $self->adjust_gene_3($g, "GT");
    #$self->trim_gene_translation5($g);
    #$self->trim_gene_translation3($g);
    $g = $self->prune_D_J_transcripts($g);
    if (@{$g->get_all_Transcripts}) {
      $g = $g->transfer($self->query);
      $g->biotype($self->J_OUTPUT_BIOTYPE);
      map { $_->biotype($self->J_OUTPUT_BIOTYPE) } @{$g->get_all_Transcripts};
      push @genes, $g;
    }
  }
  @j_genes = @genes; @genes = ();

  foreach my $g (@lv_genes, @c_genes, @j_genes, @d_genes) {
    $self->prune_Exons($g);
    $g->analysis($self->analysis);

    push @genes, $g;
  }
  @genes = sort { $a->start <=> $b->start } @genes;

  $self->output(\@genes);
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
    
    my $local_slice = $self->get_dbadaptor($self->TRANDB_DATABASES_NAME)->
        get_SliceAdaptor->fetch_by_region($self->query->coord_system->name,
                                          $self->query->seq_region_name,
                                          $gstart,
                                          $gend,
                                          $grp->{strand});
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
      } else {
        $gene_exons{$k}->
            add_supporting_features(@{$e->get_all_supporting_features});
      }
      push @t_exons, $gene_exons{$k};

      if ($tran->translation) {
        if ($e == $tran->translation->start_Exon) {
          $tran->translation->start_Exon($gene_exons{$k});
        }
        if ($e == $tran->translation->end_Exon) {
          $tran->translation->end_Exon($gene_exons{$k});
        }
      }
    }
    
    $tran->flush_Exons;
    map { $tran->add_Exon($_) } @t_exons;
  }
}


#############################################################
sub prune_LV_transcripts {
  my ($self, $gene) = @_;

  my @ftran;  
  my @tran = @{$gene->get_all_Transcripts};
  @tran = grep { $_->translate->seq !~ /\*/ } @tran;

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
  @tran = grep { $_->translate->seq !~ /\*/ } @tran;
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

  my ($best) = sort _by_total_exon_length @{$gene->get_all_Transcripts};

  if (scalar(@{$best->get_all_Exons}) > 1) {
    my ($exon, $min, $max);
    foreach my $e (@{$best->get_all_Exons}) {
      if (not defined $exon) {
        $exon = $e;
        $min = $e->start;
        $max = $e->end;
      } else {
        $min = $e->start if $e->start < $min;
        $max = $e->end   if $e->end > $max;
      }
    }
    $exon->start($min);
    $exon->end($max);

    $best->flush_Exons;
    $best->add_Exon($exon);
  }
  
  # remove translation
  $best->translation(undef);
  $best->recalculate_coordinates;
  map { $_->phase(-1); $_->end_phase(-1) } @{$best->get_all_Exons};
    
  return Bio::EnsEMBL::Gene->new(-transcripts => [$best]);
}


#########################################################
sub set_gene_stop_codon {
  my ($self, $gene) = @_;

  my $seq = uc($gene->slice->seq);
  
  my %stops = (TAA => 1,
               TGA => 1,
               TAG => 1);

  foreach my $t (@{$gene->get_all_Transcripts}) {
    my $reg_start = $t->end - 2;

    my $subseq = substr($seq, $t->end - 3, 3);

    if (not exists $stops{$subseq}) {
      my $next_seq = substr($seq, $t->end, 3);

      if (exists $stops{$next_seq}) {
        $self->adjust_transcript($t, 0, 3); 
      }
    }
  }
}

#########################################################
sub trim_gene_translation5 {
  my ($self, $gene) = @_;
  
  my $seq = uc($gene->slice->seq);
  
  my $codons_to_check = 2;
  
  foreach my $t (@{$gene->get_all_Transcripts}) {
    my $exon = $t->get_all_Exons->[0];
    my $spare = (3 - $exon->phase) % 3;
    my $adjust = 0;
    
    for(my $i=0; $i < $codons_to_check; $i++) {
      my $start = $exon->start + $spare + (3 * $i);
      my $codon = substr($seq, $start - 1, 3);
      if ($codon eq 'TAA' or 
          $codon eq 'TGA' or 
          $codon eq 'TAG') {
        $adjust = (3 * $i) + 3 + $spare;
      }
    }
    
    if ($adjust) {
      $t->translation->start($adjust);
    }
  }
}

#########################################################
sub trim_gene_translation3 {
  my ($self, $gene) = @_;
  
  my $seq = uc($gene->slice->seq);
  
  my $codons_to_check = 2;
  
  foreach my $t (@{$gene->get_all_Transcripts}) {
    my $exon = $t->get_all_Exons->[-1];
    my $spare = $exon->end_phase;
    my $adjust = 0;
    
    for(my $i=0; $i < $codons_to_check; $i++) {
      my $start = $exon->end - $spare - 2 - (3 * $i);
      my $codon = substr($seq, $start - 1, 3);
      if ($codon eq 'TAA' or 
          $codon eq 'TGA' or 
          $codon eq 'TAG') {
        $adjust = (3 * $i) + 3 + $spare;
      }
    }
    
    if ($adjust) {
      $t->translation->end($exon->length - $adjust);
    }
  }
}

#########################################################
sub adjust_gene_5 {
  my ($self, $gene, $motif) = @_;

  my $seq = uc($gene->slice->seq);

  # we will consider adding/removing a maximum of
  # 2 bps to the transcript start (to allow for
  # a partial codon at the start sometimes being
  # included in the translation, and sometimes not)

  my (@adjusted_trans);
  
  foreach my $t (@{$gene->get_all_Transcripts}) {
    my $reg_end = $t->start + 1;
    my $reg_start = $t->start - 2 - length($motif);
    my $subseq = substr($seq, $reg_start - 1, $reg_end - $reg_start + 1);

    # favour extension so take the first match
    if ((my $subseqoff = index($subseq, $motif)) >= 0) {
      my $newstart = $subseqoff + length($motif) + $reg_start;
      $self->adjust_transcript($t, $t->start - $newstart, 0);
    }
  }
}


#########################################################
sub adjust_gene_3 {
  my ($self, $gene, $motif) = @_;

  my $seq = uc($gene->slice->seq);

  # we will consider adding/removing a maximum of
  # 2 bps to the transcript end (to allow for
  # a partial codon at the end sometimes being
  # included in the translation, and sometimes not)

  foreach my $t (@{$gene->get_all_Transcripts}) {
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
  }
}

#########################################################
sub adjust_transcript {
  my ($self, $tran, $offset5, $offset3) = @_;

  # offsets in the following are how much to EXTEND transcript
  # by at each end (-> negative offsets are truncations)

  my @exons = @{$tran->get_all_Exons};

  if ($offset5) {
    $exons[0]->start($exons[0]->start - $offset5);
    $exons[0]->phase( ($exons[0]->phase - $offset5) % 3 );
    if ($tran->translation) {
      $tran->translation->end($exons[-1]->length);
    }
    $tran->recalculate_coordinates;
  }

  if ($offset3) {
    $exons[-1]->end($exons[-1]->end + $offset3);    
    $exons[-1]->end_phase( ($exons[-1]->end_phase + $offset3) % 3 );
    if ($tran->translation) {
      $tran->translation->end($exons[-1]->length);
    }
    $tran->recalculate_coordinates;
  }
}


sub gene_close_to_others {
  my ($self, $gene, $max_dist, @glists) = @_;

  # return true if gene is within 2Mb of at least
  # one other gene in each of the given lists

  foreach my $list (@glists) {
    my $member_close = 0;
    foreach my $g (@$list) {
      if (abs($g->start - $gene->start) <= $max_dist) {
        $member_close = 1;
        last;
      }
    }
    if (not $member_close) {
      return 0;
    }
  }

  return 1;
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


sub LV_OUTPUT_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_output_lv_biotype} = $val;
  }

  return $self->{_output_lv_biotype};
}


sub D_OUTPUT_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_output_d_biotype} = $val;
  }

  return $self->{_output_d_biotype};
}


sub J_OUTPUT_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_output_j_biotype} = $val;
  }

  return $self->{_output_j_biotype};
}


sub C_OUTPUT_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_output_c_biotype} = $val;
  }

  return $self->{_output_c_biotype};
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




1;
