package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar;

use warnings ;
use strict;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;
use Bio::EnsEMBL::Analysis::Tools::ClusterFilter;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(replace_stops_with_introns);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;

  my $input_id = $self->param('iid');

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');

  # Define the source transcript and target transcript dbs
  my $source_transcript_dba = $self->hrdb_get_dba($self->param('source_db'),$source_dna_dba);
  my $target_transcript_dba = $self->hrdb_get_dba($self->param('target_db'),$target_dna_dba);
  $self->hrdb_set_con($source_transcript_dba,'source_transcript_db');
  $self->hrdb_set_con($target_transcript_dba,'target_transcript_db');

  # Define the compara db
  my $compara_dba = $self->hrdb_get_dba($self->param('compara_db'),undef,'Compara');
  $self->hrdb_set_con($compara_dba,'compara_db');

#  say "DUMPER:".Dumper($compara_dba);

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

  my $gene = $source_transcript_dba->get_GeneAdaptor->fetch_by_dbID($input_id);
  my @unique_translateable_exons = $self->get_unique_translateable_exons($gene);


  my $exon_align_blocks;

  #########
  # get the compara data: MethodLinkSpeciesSet, reference DnaFrag,
  # and all GenomicAlignBlocks
  #########
  my $mlss = $compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_by_method_link_type_GenomeDBs($self->param('method_link_type'),
                                                                                                [$source_genome_db,
                                                                                                 $target_genome_db]);
  unless($mlss) {
    $self->throw("No MethodLinkSpeciesSet for :\n" .$self->param('method_link_type') . "\n" .$source_species . "\n" .$target_species);
  }

  my $dna_fragment = $compara_dba->get_DnaFragAdaptor->fetch_by_GenomeDB_and_name($source_genome_db,$gene->slice->seq_region_name);
  my $genomic_align_block_adaptor = $compara_dba->get_GenomicAlignBlockAdaptor;
  my $exon_region_padding = $self->param('exon_region_padding');
  foreach my $exon (@unique_translateable_exons) {
    my $exon_padded_start = $exon->start - $exon_region_padding;
    if($exon_padded_start < 0) {
      $exon_padded_start = 0;
    }

    my $exon_padded_end = $exon->end + $exon_region_padding;
    if($exon_padded_end > $gene->slice->length) {
      $exon_padded_end = $gene->slice->length;
    }

    my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,
                                                                                                       $dna_fragment,
                                                                                                       $exon_padded_start,
                                                                                                       $exon_padded_end);

#    my (%chains, @chains);
#    foreach my $block (@{$genomic_align_blocks}) {
#      my $source_alignment_block = $block->reference_genomic_align;
#      my ($target_alignment_blocks) = @{$block->get_all_non_reference_genomic_aligns};

      # fetch the target slice for later reference
#      if (not exists $self->target_slices->{$target_alignment_blocks->dnafrag->name}) {
#        $self->target_slices->{$target_alignment_blocks->dnafrag->name} = $target_transcript_dba->get_SliceAdaptor->fetch_by_region('toplevel',
#                                                                                                                                    $target_alignment_blocks->dnafrag->name);
#      }

#      if ($block->reference_genomic_align->dnafrag_strand < 0) {
#        $block->reverse_complement;
#      }

#      push @{$chains{$block->group_id}}, $block;
#    }

#    foreach my $chain_id (keys %chains) {
#      push(@chains,[sort {$a->reference_genomic_align->dnafrag_start <=> $b->reference_genomic_align->dnafrag_start;} @{$chains{$chain_id}}]);
#    }

#    $exon_align_block_chains->{$exon->dbID} = \@chains;
#  }

    $exon_align_blocks->{$exon->dbID} = $genomic_align_blocks;
  }
  $self->param('_exons',\@unique_translateable_exons);
  $self->param('_exon_align_blocks',$exon_align_blocks);
}

sub run {
  my ($self) = @_;

  my $exons = $self->param('_exons');
  my $results = [];
  my $fail_count = 0;
  foreach my $exon (@{$exons}) {
    my $result = $self->project_exon($exon);
    if($result) {
      push(@{$results},$result);
    } else {
      say "Failed to project exon: ".$exon->stable_id;
      $fail_count++;
    }
  }

  say "Had a total of ".$fail_count."/".scalar(@{$exons})." failed exon projections";
}


sub write_output {

}

sub project_exon {
  my ($self,$exon) = @_;
  my $exon_align_blocks = $self->param('_exon_align_blocks')->{$exon->dbID};

  unless($exon_align_blocks) {
    return 0;
  }

  my $aln_map = $self->make_alignment_mapper($exon_align_blocks);

  foreach my $exon_align_block (@{$exon_align_blocks}) {
  my ($target_blocks) = @{$exon_align_block->get_all_non_reference_genomic_aligns};
    say "Target blocks: ".$target_blocks->dnafrag->slice->name;
  }


  say "Exon: ".$exon->slice->name.":".$exon->start.":".$exon->end;

  return(1);
}

sub target_slices {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_target_slices',$val);
  }

  unless($self->param_is_defined('_target_slices')) {
    $self->param('_target_slices',{});
  }

  return $self->param('_target_slices');
}

sub get_unique_translateable_exons {
  my ($self,$gene) = @_;

  my $translateable_exons = {};
  unless($gene->biotype eq 'protein_coding') {
      $self->input_job->autoflow(0);
      $self->complete_early('Gene does not have protein_coding biotype!');
  }

  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}) {
    unless($transcript->biotype eq 'protein_coding') {
      next;
    }

    my $exons = $transcript->get_all_translateable_Exons;
    foreach my $exon (@{$exons}) {
      my $dbid = $exon->dbID();
      unless($translateable_exons->{$dbid}) {
        $translateable_exons->{$dbid} = $exon;
      }
    }
  }

  return(values(%{$translateable_exons}));
}

sub make_alignment_mapper {
  my ($self,$gen_al_blocks) = @_;

  my $FROM_CS_NAME = 'chromosome';
  my $TO_CS_NAME   = 'scaffold';

  my $mapper = Bio::EnsEMBL::Mapper->new($FROM_CS_NAME,
                                         $TO_CS_NAME);

  say "FM2 ALIGN MAP: ".ref($gen_al_blocks);
  foreach my $bl (@$gen_al_blocks) {
    foreach my $ugbl (@{$bl->get_all_ungapped_GenomicAlignBlocks}) {
      my ($from_bl) = $ugbl->reference_genomic_align;
      my ($to_bl)   = @{$ugbl->get_all_non_reference_genomic_aligns};

      say "FM2 from_bl: ".$from_bl->dnafrag_start.":".$from_bl->dnafrag_end;
      say "FM2 to_bl: ".$to_bl->dnafrag_start.":".$to_bl->dnafrag_end;
      $mapper->add_map_coordinates($from_bl->dnafrag->name,
                                   $from_bl->dnafrag_start,
                                   $from_bl->dnafrag_end,
                                   $from_bl->dnafrag_strand * $to_bl->dnafrag_strand,
                                   $to_bl->dnafrag->name,
                                   $to_bl->dnafrag_start,
                                   $to_bl->dnafrag_end);
    }
  }

  return $mapper;
}

1;
