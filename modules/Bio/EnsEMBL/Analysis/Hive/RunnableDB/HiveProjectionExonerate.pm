package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;

  my $input_ids = $self->param('iid');
  $self->param('exon_region_padding',100);
  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');

  # Define the source transcript and target transcript dbs
  my $source_transcript_dba = $self->hrdb_get_dba($self->param('source_db'));
  my $target_transcript_dba = $self->hrdb_get_dba($self->param('target_db'));
  $self->hrdb_set_con($source_transcript_dba,'source_transcript_db');
  $self->hrdb_set_con($target_transcript_dba,'target_transcript_db');

  # Define the compara db
  my $compara_dba = $self->hrdb_get_dba($self->param('compara_db'),undef,'Compara');
  $self->hrdb_set_con($compara_dba,'compara_db');

  # Get the genome db adpator
  my $genome_dba = $compara_dba->get_GenomeDBAdaptor;

  # Retrieve the production names for the query and target species
  my $source_species = $source_transcript_dba->get_MetaContainerAdaptor->get_production_name();
  my $target_species = $target_transcript_dba->get_MetaContainerAdaptor->get_production_name();

  my $source_genome_db = $genome_dba->fetch_by_core_DBAdaptor($source_transcript_dba);
  my $target_genome_db = $genome_dba->fetch_by_core_DBAdaptor($target_transcript_dba);

  ########
  # check that the default assembly for the query and target agrees
  # with that for the method_link_species_set GenomeDBs
  ########

  my $source_assembly = $source_genome_db->assembly;
  my $target_assembly = $target_genome_db->assembly;

  my ($source_assembly_version, $target_assembly_version);
  eval {
    $source_assembly_version = $source_transcript_dba->get_CoordSystemAdaptor->fetch_by_name('toplevel',$source_genome_db->assembly);
    $target_assembly_version = $target_transcript_dba->get_CoordSystemAdaptor->fetch_by_name('toplevel',$target_genome_db->assembly);
  };
  if($@) {
    $self->throw("Had trouble fetching coord systems for ".$source_genome_db->assembly . " and " .$target_genome_db->assembly . " from core dbs:\n".$@);
  }

  my $transcript_align_slices = {};
  my $mlss = $compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_by_method_link_type_GenomeDBs($self->param('method_link_type'),
                                                                                                [$source_genome_db,
                                                                                                $target_genome_db]);
  unless($mlss) {
    $self->throw("No MethodLinkSpeciesSet for :\n" .$self->param('method_link_type') . "\n" .$source_species . "\n" .$target_species);
  }

  foreach my $input_id (@$input_ids) {
    my $transcript = $source_transcript_dba->get_TranscriptAdaptor->fetch_by_dbID($input_id);
    say "Processing source transcript: ".$transcript->stable_id;
    $self->process_transcript($transcript,$compara_dba,$mlss,$source_genome_db);
  } #close foreach input_id

}

sub run {
  my ($self) = @_;

}


sub write_output {
  my ($self) = @_;

 die;

}

sub process_transcript {
  my ($self,$transcript,$compara_dba,$mlss,$source_genome_db) = @_;

  my $max_cluster_gap_length = $self->max_cluster_gap_length($transcript);
  say "Max align gap: ".$max_cluster_gap_length;
  my $all_target_genomic_aligns = [];
  my $exons = $transcript->get_all_Exons;
  foreach my $exon (@{$exons}) {
    my $exon_region_padding = $self->param('exon_region_padding');
    my $exon_padded_start = $exon->start - $exon_region_padding;
    my $exon_padded_end = $exon->end + $exon_region_padding;
    if($exon_padded_end > $exon->slice->length) {
      $exon_padded_end = $exon->slice->length;
    }
    if($exon_padded_start < 1) {
      $exon_padded_start = 1;
    }

    my $dna_fragments = $compara_dba->get_DnaFragAdaptor->fetch_by_GenomeDB_and_name($source_genome_db,$exon->slice->seq_region_name);
    my $genomic_align_block_adaptor = $compara_dba->get_GenomicAlignBlockAdaptor;
    my $genomic_align_block = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,$dna_fragments,$exon_padded_start,$exon_padded_end);

    foreach my $genomic_align_block (@$genomic_align_block) {
      my $source_genomic_align = $genomic_align_block->reference_genomic_align;
      my @target_genomic_aligns = @{$genomic_align_block->get_all_non_reference_genomic_aligns};
      foreach my $target_genomic_align (@target_genomic_aligns) {
        push(@{$all_target_genomic_aligns},$target_genomic_align);
      }
    }
  }

  my @sorted_target_genomic_aligns  = sort { $a->get_Slice->seq_region_name cmp $b->get_Slice->seq_region_name ||
                                             $a->get_Slice->start <=> $b->get_Slice->start
                                           } @{$all_target_genomic_aligns};

  unless(scalar(@sorted_target_genomic_aligns)) {
    say "No align blocks so skipping transcript";
    next;
  }

  my $transcript_slices = $self->make_cluster_slices(\@sorted_target_genomic_aligns,$max_cluster_gap_length);
  unless($transcript_slices) {
    $self->throw("The sorted align blocks were not converted into transcript slices, something went wrong");
  }

  foreach my $transcript_slice (@{$transcript_slices}) {
    say "Created transcript slice: ".$transcript_slice->name."\n";
  }
}

sub max_cluster_gap_length {
  my ($self,$transcript) = @_;

  # This sub will loop through the introns and decide on how long the max allowed value should be
  # This will be used to decide when to break up clusters on the same seq_region
  # The max value multiplier will change based on how long the biggest intron is
  # Have 1000 as a min default value to handle small breaks in single exon genes
  my $longest_intron = 50;
  my $introns = $transcript->get_all_Introns();
  unless(scalar(@{$introns}) > 0) {
    return($longest_intron);
  }

  foreach my $intron (@{$introns}) {
    if($intron->length > $longest_intron) {
      $longest_intron = $intron->length;
    }
  }

  if($longest_intron <= 10000) {
    return(int($longest_intron * 3));
  } elsif($longest_intron <= 50000) {
    return(int($longest_intron * 2));
  } else {
    return(int($longest_intron * 1.5))
  }
}

sub make_cluster_slices {
  my ($self,$genomic_aligns,$max_cluster_gap_length) = @_;

  my $slice_adaptor = $self->hrdb_get_con('target_dna_db')->get_SliceAdaptor;
  my $cluster_slices = [];
  my $previous_genomic_align = shift(@{$genomic_aligns});
  my $cluster_start = $previous_genomic_align->get_Slice->start();
  my $cluster_end = $previous_genomic_align->get_Slice->end();
  unless($previous_genomic_align) {
    return;
  }

  foreach my $current_genomic_align (@{$genomic_aligns}) {
    if($previous_genomic_align->get_Slice->seq_region_name ne $current_genomic_align->get_Slice->seq_region_name) {
      my $slice = $slice_adaptor->fetch_by_region($previous_genomic_align->get_Slice->coord_system_name,
                                                  $previous_genomic_align->get_Slice->seq_region_name, $cluster_start, $cluster_end);
      push(@{$cluster_slices},$slice);
      $cluster_start = $current_genomic_align->get_Slice->start();
      $cluster_end = $current_genomic_align->get_Slice->end();
      $previous_genomic_align = $current_genomic_align;
    } elsif($current_genomic_align->get_Slice->start - $previous_genomic_align->get_Slice->end > $max_cluster_gap_length) {
      my $slice = $slice_adaptor->fetch_by_region($previous_genomic_align->get_Slice->coord_system_name,
                                                  $previous_genomic_align->get_Slice->seq_region_name, $cluster_start, $cluster_end);
      push(@{$cluster_slices},$slice);
      $cluster_start = $current_genomic_align->get_Slice->start();
      $cluster_end = $current_genomic_align->get_Slice->end();
      $previous_genomic_align = $current_genomic_align;
    } else {
      $cluster_end = $current_genomic_align->get_Slice->end();
    }
  } # end foreach my $current_genomic_align

  # push the final slice
  my $slice = $slice_adaptor->fetch_by_region($previous_genomic_align->get_Slice->coord_system_name,
                                              $previous_genomic_align->get_Slice->seq_region_name, $cluster_start, $cluster_end);
  push(@{$cluster_slices},$slice);

  return($cluster_slices);
}

1;
