=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::AlignmentChains - 

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $genscan = Bio::EnsEMBL::Analysis::RunnableDB::AlignmentChains->new (
                                                    -db      => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $genscan->fetch_input();
  $genscan->run();
  $genscan->write_output(); #writes to DB


=head1 DESCRIPTION

Given an compara MethodLinkSpeciesSet identifer, and a reference genomic
slice identifer, fetches the GenomicAlignBlocks from the given compara
database, forms them into sets of alignment chains, and writes the result
back to the database. 

This module (at least for now) relies heavily on Jim Kent\'s Axt tools.
=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::AlignmentChains;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Config::General;

use Bio::EnsEMBL::Analysis::RunnableDB::AlignmentFilter;
use Bio::EnsEMBL::Analysis::Config::AlignmentFilter;
use Bio::EnsEMBL::Analysis::Runnable::AlignmentChains;

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
 
  $self->read_and_check_config($CHAIN_CONFIG_BY_LOGIC);

  return $self;
}




=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my( $self) = @_; 
  
  throw("No input id") unless defined($self->input_id);

  my ($seq_name, $target_chunk_total, $target_chunk_id);
  if ($self->input_id =~ /^([^:]+):(\d+):(\d+)$/) {
    ($seq_name, $target_chunk_total, $target_chunk_id) = ($1, $2, $3);
  } elsif ($self->input_id =~ /(\S+)/) {
    ($seq_name, $target_chunk_total, $target_chunk_id) = ($1, 1, 1);
  } else {
    throw("Input id could not be parsed: ", $self->input_id);
  }
  
  $self->OUTPUT_GROUP_TYPE("chain") unless (defined $self->OUTPUT_GROUP_TYPE);
  my $q_dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$self->QUERY_CORE_DB});
  my $t_dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$self->TARGET_CORE_DB});

  my $compara_dbh = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(%{$self->COMPARA_DB});
  
  my $query_species = $q_dbh->get_MetaContainerAdaptor->get_Species->binomial;
  my $target_species = $t_dbh->get_MetaContainerAdaptor->get_Species->binomial;
  
  my $q_gdb = $compara_dbh->get_GenomeDBAdaptor->fetch_by_name_assembly($query_species);
  my $t_gdb = $compara_dbh->get_GenomeDBAdaptor->fetch_by_name_assembly($target_species);
  
  ################################################
  # check that the default assembly for the query and target agrees with that
  # for the method_link_species_set GenomeDBs
  #################################################

  my ($q_assembly_version, $t_assembly_version);
  eval {
    $q_assembly_version = $q_dbh->get_CoordSystemAdaptor->fetch_by_name
        ('toplevel',
         $q_gdb->assembly);

    $t_assembly_version = $t_dbh->get_CoordSystemAdaptor->fetch_by_name
        ('toplevel',
         $t_gdb->assembly);
  };
  $@ and do {
    throw("Had trouble fetching coord systems for ". 
          $q_gdb->assembly . " and " . 
          $t_gdb->assembly . " from core dbs: $@");
  };

  #############
  # query Slice
  #############

  # since compara only stores toplevel DNA fragments, we assume
  # supplied slice identifer corresponds to a top-level region.
  my $ref_slice = $q_dbh->get_SliceAdaptor->fetch_by_region('toplevel', 
                                                            $seq_name,
                                                            undef,
                                                            undef,
                                                            undef,
                                                            $q_assembly_version);

  throw("Could not fetch top level query slice for '$seq_name'") if not defined $ref_slice;    

  ################################################################
  # get the compara data: MethodLinkSpeciesSet, reference DnaFrag, 
  # and all GenomicAlignBlocks
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


  print STDERR "Fetching all GenomicAlignBlocks and sorting them by target...\n";

  my %blocks_by_target;
  foreach my $block (@{$compara_dbh->get_GenomicAlignBlockAdaptor
                           ->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss, $ref_dnafrag)}) {
        
    #my ($qy_al) = $block->reference_genomic_align;
    my ($tg_al) = @{$block->get_all_non_reference_genomic_aligns};
    my $tg_name = $tg_al->dnafrag->name;

    # the following awful hack releases the contained genomic_aligns;
    # since we wmay only be using a fraction of the total set of blocks,
    # the tome lost in having to requery for the relevant genomic_aligns
    # later is more than offset by a big memory saving
    $block->{'genomic_align_array'} = undef;

    push @{$blocks_by_target{$tg_name}}, $block;
  }

  print STDERR "Gathering sequence converting features for id group $target_chunk_id of $target_chunk_total...\n";

  ###################################################################
  # Fetch slices and make features for blocks involving only the
  # designated subset of target sequences
  ###################################################################
  my @all_target_names = sort keys %blocks_by_target;

  my (%target_slices, %features_by_target);

  for(my $i = $target_chunk_id - 1; 
      $i < @all_target_names; 
      $i += $target_chunk_total) {
   
    my $target_name = $all_target_names[$i];
   
    $target_slices{$target_name} = 
        $t_dbh->get_SliceAdaptor->fetch_by_region('toplevel', 
                                                  $target_name,
                                                  undef,
                                                  undef,
                                                  undef,
                                                  $t_assembly_version); 

    foreach my $target_block (@{$blocks_by_target{$target_name}}) {
      my ($qy_al) = $target_block->reference_genomic_align;
      my ($tg_al) = @{$target_block->get_all_non_reference_genomic_aligns};

      if (not exists($self->query_DnaFrag_hash->{$qy_al->dnafrag->name})) {
        ######### needed for output ######################################
        $self->query_DnaFrag_hash->{$qy_al->dnafrag->name} = $qy_al->dnafrag;
      }
      
      if (not exists($self->target_DnaFrag_hash->{$target_name})) {
        ######### needed for output #######################################
        $self->target_DnaFrag_hash->{$tg_al->dnafrag->name} = $tg_al->dnafrag;
      }
      
      my $daf_cigar = $self->daf_cigar_from_compara_cigars($qy_al->cigar_line,
                                                           $tg_al->cigar_line);

      if (defined $daf_cigar) {
        my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new
            (-seqname => $qy_al->dnafrag->name,
             -start    => $qy_al->dnafrag_start,
             -end      => $qy_al->dnafrag_end,
             -strand   => $qy_al->dnafrag_strand,
             -hseqname => $tg_al->dnafrag->name,
             -hstart   => $tg_al->dnafrag_start,
             -hend     => $tg_al->dnafrag_end,
             -hstrand   => $tg_al->dnafrag_strand,         
             -cigar_string => $daf_cigar);
        
        push @{$features_by_target{$tg_al->dnafrag->name}}, $daf;      
      }
    }
  }

  print STDERR "Sorting resulting features into batches...\n";
  ######################################################################
  # each runnable comprises blocks from a subset of the target sequences
  ######################################################################
  my $target_batches = $self->form_target_batches(\%target_slices);
  foreach my $targets (@$target_batches) {

    my (%these_target_slices, @features);
    foreach my $t_id (@$targets) {
      $these_target_slices{$t_id} =  $target_slices{$t_id};
      push @features, @{$features_by_target{$t_id}};
    }

    printf(STDERR "Making runnable with %d targets (%d features)\n", scalar(@$targets), scalar(@features));
      
    my %parameters = (-analysis             => $self->analysis, 
                      -query_slice          => $ref_slice,
                      -target_slices        => \%these_target_slices,
                      -query_nib_dir        => $self->QUERY_NIB_DIR,
                      -target_nib_dir       => $self->TARGET_NIB_DIR,
                      -min_chain_score      => $self->MIN_CHAIN_SCORE,
                      -features             => \@features);

    foreach my $program (qw(faToNib lavToAxt axtChain)) {
      $parameters{'-' . $program} = $BIN_DIR . "/" . $program;
    }
    
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::AlignmentChains->new(%parameters);    
    $self->runnable($runnable);
  }

}


###########################################

sub form_target_batches {
  my ($self, $t_slices) = @_;

  my @batches;
  my $batch_index = 0;
  my ($total_len, $total_count) = (0,0);

  foreach my $hname (keys %{$t_slices}) {
    push @{$batches[$batch_index]}, $hname;
    $total_count++;
    
    if ($total_count >= 1000) {
      $total_count = 0;      
      $batch_index++;
    }
  }

  return \@batches;
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
                      QUERY_CORE_DB
                      TARGET_CORE_DB
                      COMPARA_DB)) {

    throw("You must define $var in config for logic '$logic'" . 
          " or in the DEFAULT entry")
        if not $self->$var;
  }
}


sub QUERY_CORE_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->{_query_core_db} = $db;
  }

  return $self->{_query_core_db};
}


sub TARGET_CORE_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->{_target_core_db} = $db;
  }

  return $self->{_target_core_db};
}


sub QUERY_NIB_DIR {
  my ($self, $dir) = @_;

  if (defined $dir) {
    $self->{_query_nib_dir} = $dir;
  }

  return $self->{_query_nib_dir};
}


sub TARGET_NIB_DIR {
  my ($self, $dir) = @_;

  if (defined $dir) {
    $self->{_target_nib_dir} = $dir;
  }

  return $self->{_target_nib_dir};
}



1;
