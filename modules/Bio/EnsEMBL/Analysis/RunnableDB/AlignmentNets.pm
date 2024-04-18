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

Bio::EnsEMBL::Analysis::RunnableDB::AlignmentNets - 

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $genscan = Bio::EnsEMBL::Analysis::RunnableDB::AlignmentNets->new (
                                                    -db      => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $genscan->fetch_input();
  $genscan->run();
  $genscan->write_output(); #writes to DB


=head1 DESCRIPTION

Given an compara MethodLinkSpeciesSet identifer, and a reference genomic
slice identifer, fetches the GenomicAlignBlocks from the given compara
database, infers chains from the group identifiers, and then forms
an alignment net from the chains and writes the result
back to the database. 

This module (at least for now) relies heavily on Jim Kent\'s Axt tools.
=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::AlignmentNets;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Config::General;

use Bio::EnsEMBL::Analysis::RunnableDB::AlignmentFilter;
use Bio::EnsEMBL::Analysis::Config::AlignmentFilter;
use Bio::EnsEMBL::Analysis::Runnable::AlignmentNets;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::AlignmentFilter);

############################################################
sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
 
  $self->read_and_check_config($NET_CONFIG_BY_LOGIC);

  return $self;
}


############################################################
sub fetch_input {
  my( $self) = @_; 
  
  $self->throw("No input id") unless defined($self->input_id);
  
  my ($seq_name, $seq_start, $seq_end);
  if ($self->input_id =~ /^([^:]+):(\d+):(\d+)$/) {
    ($seq_name, $seq_start, $seq_end) = ($1, $2, $3);
  } elsif ($self->input_id =~ /(\S+)/) {
    $seq_name = $1;
  } else {
    throw("Input id could not be parsed: ", $self->input_id);
  }

  if ($self->QUERY_CORE_DB) {
    my $q_dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$self->QUERY_CORE_DB});
    my $sl = $q_dbh->get_SliceAdaptor->fetch_by_region('toplevel',
                                                       $seq_name);
    my @segs = @{$sl->project('seqlevel')};
    $self->query_seq_level_projection(\@segs);

    #foreach my $seg (@segs) {
    #  printf("FROM_START %d FROM_END %d TO_SLICE %s TO_START %d TO_END %d\n", 
    #         $seg->from_start, $seg->from_end, $seg->to_Slice->seq_region_name,
    #         $seg->to_Slice->start, $seg->to_Slice->end);
    #}
  }


  my $compara_dbh = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(%{$self->COMPARA_DB});

  my $query_species = $self->QUERY_SPECIES;
  my $target_species = $self->TARGET_SPECIES;
  
  my $q_gdb = $compara_dbh->get_GenomeDBAdaptor->fetch_by_name_assembly($query_species);
  my $t_gdb = $compara_dbh->get_GenomeDBAdaptor->fetch_by_name_assembly($target_species);

  throw("Could not get GenomeDB for '$query_species'") if not defined $q_gdb;
  throw("Could not get GenomeDB for '$target_species'") if not defined $t_gdb;
  
  ################################################################
  # get the compara data: MethodLinkSpeciesSet, reference DnaFrag, 
  # and GenomicAlignBlocks
  ################################################################
  my $mlss = $compara_dbh->get_MethodLinkSpeciesSetAdaptor
      ->fetch_by_method_link_type_GenomeDBs($self->INPUT_METHOD_LINK_TYPE,
                                            [$q_gdb, $t_gdb]);
  throw("No MethodLinkSpeciesSet for :\n" .
        $self->INPUT_METHOD_LINK_TYPE . "\n" . 
        $query_species . "\n" . 
        $target_species)
      if not $mlss;

  my $out_mlss = $compara_dbh->get_MethodLinkSpeciesSetAdaptor
      ->fetch_by_method_link_type_GenomeDBs($self->OUTPUT_METHOD_LINK_TYPE,
                                            [$q_gdb, $t_gdb]);

  throw("No MethodLinkSpeciesSet for :\n" .
        $self->OUTPUT_METHOD_LINK_TYPE . "\n" . 
        $query_species . "\n" . 
        $target_species)
      if not $out_mlss;

  ######## needed for output####################
  $self->output_MethodLinkSpeciesSet($out_mlss);
  
  my $ref_dnafrag = $compara_dbh->get_DnaFragAdaptor->fetch_by_GenomeDB_and_name($q_gdb,
                                                                                 $seq_name);

  my $gen_al_blocks = $compara_dbh->get_GenomicAlignBlockAdaptor
      ->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss, $ref_dnafrag, $seq_start, $seq_end);

  ###################################################################
  # get the target slices and bin the GenomicAlignBlocks by group id
  ###################################################################
  my (%features_by_group);
  
  foreach my $block (@$gen_al_blocks) {
    my ($qy_al) = $block->reference_genomic_align;
    my ($tg_al) = @{$block->get_all_non_reference_genomic_aligns};

    if (not exists($self->query_DnaFrag_hash->{$qy_al->dnafrag->name})) {
      ######### needed for output ######################################
      $self->query_DnaFrag_hash->{$qy_al->dnafrag->name} = $qy_al->dnafrag;
    }
    if (not exists($self->target_DnaFrag_hash->{$tg_al->dnafrag->name})) {
      ######### needed for output #######################################
      $self->target_DnaFrag_hash->{$tg_al->dnafrag->name} = $tg_al->dnafrag;
    }

    my $daf_cigar = $self->daf_cigar_from_compara_cigars($qy_al->cigar_line,
                                                         $tg_al->cigar_line);
    
    my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new
        (-seqname => $qy_al->dnafrag->name,
         -start    => $qy_al->dnafrag_start,
         -end      => $qy_al->dnafrag_end,
         -strand   => $qy_al->dnafrag_strand,
         -hseqname => $tg_al->dnafrag->name,
         -hstart   => $tg_al->dnafrag_start,
         -hend     => $tg_al->dnafrag_end,
         -hstrand  => $tg_al->dnafrag_strand,         
         -score    => $block->score,
         -cigar_string => $daf_cigar,
         -group_id => $block->group_id);

    push @{$features_by_group{$block->group_id}}, $daf;
  }

  $self->chains($self->sort_chains_by_max_block_score([values %features_by_group]));

  if ($self->METHOD eq 'STANDARD' or $self->METHOD eq 'SYNTENIC') {
    my $run = $self->make_runnable($self->METHOD eq 'SYNTENIC');
    $self->runnable($run);
  }
}

############################################################
sub run{
  my ($self) = @_;

  if (@{$self->runnable}) {
    $self->SUPER::run;
  } else {
    my $filtered_chains;
    
    if ($self->METHOD eq 'SIMPLE_HIGH') {
      $filtered_chains = $self->calculate_simple_high_net($self->chains);
    } elsif ($self->METHOD eq 'SIMPLE_LOW') {
      $filtered_chains = $self->calculate_simple_low_net($self->chains);
    } 
    
    if (defined $filtered_chains) {
      my $converted_out = $self->convert_output($filtered_chains);
      $self->output($converted_out);
    }
  }
}
    

###############################################################
sub calculate_simple_high_net {
  my ($self, $chains) = @_;

  my (@net_chains, @retained_blocks, %contigs_of_kept_blocks);

  # Junk chains that have extent overlap
  # with retained chains so far
  foreach my $c (@$chains) {
    my @b = sort { $a->start <=> $b->start } @$c; 
    my $c_st = $b[0]->start;
    my $c_en = $b[-1]->end;
    
    my $keep_chain = 1;
    foreach my $oc (@net_chains) {
      my @ob = sort { $a->start <=> $b->start } @$oc; 
      
      my $oc_st = $ob[0]->start;
      my $oc_en = $ob[-1]->end;
      
      if ($oc_st <= $c_en and $oc_en >= $c_st) {
        # overlap; junk
        $keep_chain = 0;
        last;
      }
    }
    
    if ($keep_chain) {
      my (%contigs_of_blocks, @split_blocks);

      foreach my $bl (@b) {
        my ($inside_seg, @overlap_segs);

        foreach my $seg (@{$self->query_seq_level_projection}) {
          if ($bl->start >= $seg->from_start and
              $bl->end    <= $seg->from_end) {
            $inside_seg = $seg;
            last;
          } elsif ($seg->from_start <= $bl->end and 
              $seg->from_end   >= $bl->start) {
            push @overlap_segs, $seg;
          } elsif ($seg->from_start > $bl->end) {
            last;
          }          
        }
        if (defined $inside_seg) {
          push @split_blocks, $bl;
          $contigs_of_blocks{$bl} = $inside_seg;
        } else {
          my @cut_blocks;
          foreach my $seg (@overlap_segs) {
            my ($reg_start, $reg_end) = ($bl->start, $bl->end);
            $reg_start = $seg->from_start if $seg->from_start > $reg_start;
            $reg_end   = $seg->from_end   if $seg->from_end   < $reg_end;

            my $cut_block = $bl->restrict_between_positions($reg_start, 
                                                            $reg_end, 
                                                            "SEQ");
            if (defined $cut_block) {
              push @cut_blocks, $cut_block;
              $contigs_of_blocks{$cut_block} = $seg;
            }
          }
          push @split_blocks, @cut_blocks;
        }
      }
      @b  = @split_blocks;

      my %new_block;
      map { $new_block{$_} = 1 } @b;
      
      my @new_list = (@retained_blocks, @b);
      @new_list = sort { $a->start <=> $b->start } @new_list;
      
    NEWBLOCK:
      for(my $i = 0; $i < @new_list; $i++) {
        my $bl = $new_list[$i];

        if (exists($new_block{$bl})) {

          my ($flank_left, $flank_right);

          if ($i > 0 and 
              not exists($new_block{$new_list[$i-1]})) {
            $flank_left = $new_list[$i-1];
          }
          if ($i < scalar(@new_list) - 1 and 
              not exists($new_block{$new_list[$i+1]})) {
            $flank_right = $new_list[$i+1];
          }
          
          if (defined $flank_left and
              $contigs_of_blocks{$bl}->to_Slice->seq_region_name eq
              $contigs_of_kept_blocks{$flank_left}->to_Slice->seq_region_name
              ) {
            $keep_chain = 0;
            last NEWBLOCK;
          }
          
          if (defined $flank_right and
              $contigs_of_blocks{$bl}->to_Slice->seq_region_name eq
              $contigs_of_kept_blocks{$flank_right}->to_Slice->seq_region_name
              ) {
            $keep_chain = 0;
            last NEWBLOCK;
          }
        }
      }

      if ($keep_chain) {
        foreach my $bid (keys %contigs_of_blocks) {
          $contigs_of_kept_blocks{$bid} = $contigs_of_blocks{$bid};
        }

        push @net_chains, \@b;
        push @retained_blocks, @b;
        @retained_blocks = sort { $a->start <=> $b->start } @retained_blocks;
      }
    }
  }
  
  return \@net_chains;
}


################################################################
sub calculate_simple_low_net {
  my ($self, $chains) = @_;


  my (@net_chains, @retained_blocks, %contigs_of_kept_blocks);
  
  foreach my $c (@$chains) {

    my @b = sort { $a->start <=> $b->start } @$c; 
    
    my $keep_chain = 1;
    BLOCK: foreach my $b (@b) {
      OTHER_BLOCK: foreach my $ob (@retained_blocks) {
        if ($ob->start <= $b->end and $ob->end >= $b->start) {
          $keep_chain = 0;
          last BLOCK;
        } elsif ($ob->start > $b->end) {
          last OTHER_BLOCK;
        }
      }
    }
    
    if ($keep_chain) {
      my (%contigs_of_blocks, @split_blocks);
      
      foreach my $bl (@b) {
        my ($inside_seg, @overlap_segs);

        foreach my $seg (@{$self->query_seq_level_projection}) {
          if ($bl->start >= $seg->from_start and
              $bl->end    <= $seg->from_end) {
            $inside_seg = $seg;
            last;
          } elsif ($seg->from_start <= $bl->end and 
              $seg->from_end   >= $bl->start) {
            push @overlap_segs, $seg;
          } elsif ($seg->from_start > $bl->end) {
            last;
          }          
        }
        if (defined $inside_seg) {
          push @split_blocks, $bl;
          $contigs_of_blocks{$bl} = $inside_seg;
        } else {
          my @cut_blocks;
          foreach my $seg (@overlap_segs) {
            my ($reg_start, $reg_end) = ($bl->start, $bl->end);
            $reg_start = $seg->from_start if $seg->from_start > $reg_start;
            $reg_end   = $seg->from_end   if $seg->from_end   < $reg_end;

            my $cut_block = $bl->restrict_between_positions($reg_start, 
                                                            $reg_end, 
                                                            "SEQ");
            if (defined $cut_block) {
              push @cut_blocks, $cut_block;
              $contigs_of_blocks{$cut_block} = $seg;
            }
          }
          push @split_blocks, @cut_blocks;
        }
      }
      @b = @split_blocks;

      my %new_block;
      map { $new_block{$_} = 1 } @b;

      my @new_list = (@retained_blocks, @b);
      @new_list = sort { $a->start <=> $b->start } @new_list;

    NEWBLOCK:
      for(my $i = 0; $i < @new_list; $i++) {
        my $bl = $new_list[$i];

        if (exists($new_block{$bl})) {

          my ($flank_left, $flank_right);

          if ($i > 0 and 
              not exists($new_block{$new_list[$i-1]})) {
            $flank_left = $new_list[$i-1];
          }
          if ($i < scalar(@new_list) - 1 and 
              not exists($new_block{$new_list[$i+1]})) {
            $flank_right = $new_list[$i+1];
          }
          
          if (defined $flank_left and
              #($flank_left->hseqname ne $bl->hseqname or
              # $flank_left->hstrand != $bl->hstrand or
              # ($bl->hstrand > 0 and $flank_left->hend >= $bl->hstart) or
              # ($bl->hstrand < 0 and $flank_left->hstart <= $bl->hend)) and
              $contigs_of_blocks{$bl}->to_Slice->seq_region_name eq
              $contigs_of_kept_blocks{$flank_left}->to_Slice->seq_region_name
              ) {
            $keep_chain = 0;
            last NEWBLOCK;
          }
          
          if (defined $flank_right and
              #($flank_right->hseqname ne $bl->hseqname or
              # $flank_right->hstrand != $bl->hstrand or
              # ($bl->hstrand > 0 and $flank_right->hstart <= $bl->hend) or
              # ($bl->hstrand < 0 and $flank_right->hend >= $bl->hstart)) and
              $contigs_of_blocks{$bl}->to_Slice->seq_region_name eq
              $contigs_of_kept_blocks{$flank_right}->to_Slice->seq_region_name
              ) {
            $keep_chain = 0;
            last NEWBLOCK;
          }
        }
      }
  
      if ($keep_chain) {  
        foreach my $bid (keys %contigs_of_blocks) {
          $contigs_of_kept_blocks{$bid} = $contigs_of_blocks{$bid};
        }
        push @net_chains, \@b;
        push @retained_blocks, @b;
        @retained_blocks = sort { $a->start <=> $b->start } @retained_blocks;
      }
    }
  }

  return \@net_chains;
}


############################################################
sub make_runnable {
  my ($self, $syntenic) = @_;

  my (%query_lengths, %target_lengths);

  foreach my $nm (keys %{$self->query_DnaFrag_hash}) {
    $query_lengths{$nm} = $self->query_DnaFrag_hash->{$nm}->length;
  }
  foreach my $nm (keys %{$self->target_DnaFrag_hash}) {
    $target_lengths{$nm} = $self->target_DnaFrag_hash->{$nm}->length;
  }
    
  my %parameters = (-analysis             => $self->analysis, 
                    -chains               => $self->chains,
                    -query_lengths        => \%query_lengths,
                    -target_lengths       => \%target_lengths,
                    -min_chain_score      => $self->MIN_CHAIN_SCORE,
                    -filter_non_syntenic  => $syntenic);
  
  $parameters{-chainNet} = $self->CHAIN_NET 
      ? $self->CHAIN_NET
      : $BIN_DIR . "/" . "chainNet";
  $parameters{-netSyntenic} = $self->NET_SYNTENIC 
      ? $self->NET_SYNTENIC 
      : $BIN_DIR . "/" . "netSyntenic";
  $parameters{-netFilter} = $self->NET_FILTER
      ? $self->NET_FILTER
      : $BIN_DIR . "/" . "netFilter";
    
  my $run = Bio::EnsEMBL::Analysis::Runnable::AlignmentNets->new(%parameters);
  return $run;

}


##############################################################
sub chains {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_chains} = $value;
  }
  return $self->{_chains};
}


sub query_seq_level_projection {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_query_proj_segs} = $val;
  }
  return $self->{_query_proj_segs};
}


####################################
# config variable holders
####################################

sub read_and_check_config {
  my ($self, $hash) = @_;

  $self->SUPER::read_and_check_config($hash);
 
  my $logic = $self->analysis->logic_name;

  foreach my $var (qw(INPUT_METHOD_LINK_TYPE
                      OUTPUT_METHOD_LINK_TYPE
                      QUERY_SPECIES
                      TARGET_SPECIES
                      COMPARA_DB
                      METHOD)) {

    throw("You must define $var in config for logic '$logic'" . 
          " or in the DEFAULT entry")
        if not $self->$var;
  }

  # check for sanity
  my %allowable_methods = 
      (
       STANDARD      => 1,
       SYNTENIC      => 1,
       SIMPLE_HIGH   => 1,
       SIMPLE_MEDIUM => 1,
       SIMPLE_LOW    => 1,
       );

  if ($self->METHOD and 
      not $allowable_methods{$self->METHOD}) {
    throw("You must set METHOD to one of the reserved names\n" .
          "See the .example file for these names");
  }

}

sub QUERY_CORE_DB {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_query_core_db} = $val;
  }
  
  return $self->{_query_core_db};
}


sub QUERY_SPECIES {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_query_species} = $val;
  }
  
  return $self->{_query_species};
}

sub TARGET_SPECIES {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_target_species} = $val;
  }
  
  return $self->{_target_species};
}


sub METHOD {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_primary_method} = $val;
  }
  
  return $self->{_primary_method};
}


sub CHAIN_NET {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_chain_net_prog} = $val;
  }
  
  return $self->{_chain_net_prog};
}


sub NET_SYNTENIC {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_net_syntenic_prog} = $val;
  }
  
  return $self->{_net_syntenic_prog};
}

sub NET_FILTER {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_net_filter_prog} = $val;
  }
  
  return $self->{_net_filter_prog};
}



1;
