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

  if ($self->PRIMARY_METHOD eq 'STANDARD' or 
      $self->SECONDARY_METHOD eq 'STANDARD') {
    my $run = $self->make_runnable(0);
    $self->standard_runnable($run);
  }
  if ($self->PRIMARY_METHOD eq 'SYNTENIC' or 
      $self->SECONDARY_METHOD eq 'SYNTENIC') {
    my $run = $self->make_runnable(1);
    $self->syntenic_runnable($run);
  }
}

############################################################
sub run{
  my ($self) = @_;
  
  $self->run_method($self->PRIMARY_METHOD);
  if (not @{$self->output} and $self->SECONDARY_METHOD) {
    $self->run_method($self->SECONDARY_METHOD);
  }

}

#############################################################
sub run_method {
  my ($self, $meth) = @_;

  if ($meth eq 'STANDARD' or $meth eq 'SYNTENIC') {
    $self->runnable( $meth eq 'STANDARD' 
                     ? $self->standard_runnable 
                     : $self->syntenic_runnable ); 

    $self->SUPER::run;
  } else {
    my $filtered_chains;

    if ($meth eq 'SIMPLE_HIGH') {
      $filtered_chains = [$self->chains->[0]];
    } elsif ($meth eq 'SIMPLE_MEDIUM') {
      $filtered_chains = $self->calculate_simple_medium_net($self->chains);
    } elsif ($meth eq 'SIMPLE_LOW') {
      $filtered_chains = $self->calculate_simple_low_net($self->chains);
    } 

    if (defined $filtered_chains) {
      my $converted_out = $self->convert_output($filtered_chains);
      $self->output($converted_out);
    }
  }
}


###############################################################
sub calculate_simple_medium_net {
  my ($self, $chains) = @_;

  my (@net_chains);
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
  
  return \@net_chains;
}


################################################################
sub calculate_simple_low_net {
  my ($self, $chains) = @_;

  my (@net_chains, @retained_blocks);
  
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
sub runnable {
  my ($self, $runnable) = @_;

  my @runnables;

  if (defined $runnable) {
    $self->{_single_runnable} = $runnable;
  }
  if (exists $self->{_single_runnable}) {
    push @runnables, $self->{_single_runnable};
  }

  return \@runnables;
}

##############################################################
sub standard_runnable {
  my ($self, $val) = @_;
  
  if (defined $val) {
    $self->{_standard_runnable} = $val;
  }

  return $self->{_standard_runnable};
}

##############################################################
sub syntenic_runnable {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_syntenic_runnable} = $val;
  }

  return $self->{_syntenic_runnable};
}

##############################################################
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
                      COMPARA_DB
                      PRIMARY_METHOD)) {

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

  if ($self->PRIMARY_METHOD and 
      not $allowable_methods{$self->PRIMARY_METHOD}) {
    throw("You must set PRIMARY_METHOD to one of the reserved names\n" .
          "See the .example file for these names");
  }
  if ($self->SECONDARY_METHOD and 
      not $allowable_methods{$self->PRIMARY_METHOD}) {
    throw("You must set SECONDARY_METHOD to one of the reserved names\n" .
          "See the .example file for these names");
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


sub PRIMARY_METHOD {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_primary_method} = $val;
  }
  
  return $self->{_primary_method};
}

sub SECONDARY_METHOD {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_secondary_method} = $val;
  }
  
  return $self->{_secondary_method};
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
