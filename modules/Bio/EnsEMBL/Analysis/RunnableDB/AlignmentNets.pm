# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::AlignmentNets

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


=cut
package Bio::EnsEMBL::Analysis::RunnableDB::AlignmentNets;

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




=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   
    Returns :   nothing
    Args    :   none

=cut

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

  $self->GROUP_TYPE("chain") unless (defined $self->GROUP_TYPE);
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
         -cigar_string => $daf_cigar);

    my $group_id = $qy_al->genomic_align_group_id_by_type("chain");
    if ($group_id != $tg_al->genomic_align_group_id_by_type("chain")) {
      throw("GenomicAligns in a GenomicAlignBlock belong to different group");
    }
    
    push @{$features_by_group{$group_id}}, $daf;
  }

  $self->chains($self->sort_chains_by_max_block_score([values %features_by_group]));

  if (not $self->SIMPLE_NET) {
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
                      -filter_non_syntenic  => $self->FILTER_NON_SYNTENIC);
    
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
    $self->runnable($run);
  } 
}



sub calculate_simple_net {
  my ($self, $mode, $chains) = @_;

  my @net_chains;

  if ($mode =~ /HIGH/) {
    # VERY simple: just keep the best chain

    if (@$chains) {
      @net_chains = ($chains->[0]);
    }

  } elsif ($mode =~ /MEDIUM/) {
    # Slightly less simple. Junk chains that have extent overlap
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
        push @net_chains, $c;
      }
    }

  } elsif ($mode =~ /LOW/) {
    # not quite so simple. Junk chains that have BLOCK overlap
    # with retained chains so far
    my @retained_blocks;

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
        push @net_chains, $c;
        push @retained_blocks, @$c;
        @retained_blocks = sort { $a->start <=> $b->start } @retained_blocks;
      }
    }
  }

  return \@net_chains;
}

=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : cycles through all the runnables, calls run and pushes
  their output into the RunnableDBs output array
  Returntype: array ref
  Exceptions: none
  Example   : 

=cut

sub run{
  my ($self) = @_;
  if (not $self->SIMPLE_NET) {
    $self->SUPER::run();
  } else {
    my $filtered_chains = $self->calculate_simple_net($self->SIMPLE_NET,
                                                      $self->chains);
    my $converted_output = $self->convert_output($filtered_chains);
    $self->output($converted_output);
  }
}


=head2 chains

 Title   : chains
 Usage   :
 Function:
 Returns : 
 Args    :

=cut

sub chains {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_chains} = $value;
  }
  return $self->{_chains};
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
                      COMPARA_DB)) {

    throw("You must define $var in config for logic '$logic'" . 
          " or in the DEFAULT entry")
        if not $self->$var;
  }
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

sub SIMPLE_NET {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_simple_net} = $val;
  }

  return $self->{_simple_net};
}

sub FILTER_NON_SYNTENIC {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_filter_non_syntenic} = $val;
  }
  
  return $self->{_filter_non_syntenic};
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
