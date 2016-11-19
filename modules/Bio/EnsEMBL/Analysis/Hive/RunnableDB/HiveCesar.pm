package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
#use File::chdir;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;
use Bio::EnsEMBL::Analysis::Tools::ClusterFilter;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(replace_stops_with_introns calculate_exon_phases);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;

  unless(-e $self->param('output_path')) {
    system("mkdir -p ".$self->param('output_path'));
  }

  my $input_id = $self->param('iid');

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
  $self->parent_gene($gene);

  my @unique_translateable_exons = $self->get_unique_translateable_exons($gene);


  my $exon_align_slices;

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

    my $slice_adaptor = $source_dna_dba->get_SliceAdaptor();
    my $exon_slice = $slice_adaptor->fetch_by_region($exon->slice->coord_system_name, $exon->slice->seq_region_name, $exon_padded_start, $exon_padded_end);

    my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $exon_slice);
    my $exon_slices = [];
    foreach my $genomic_align_block (@{$genomic_align_blocks}) {
      my $restricted_gab = $genomic_align_block->restrict_between_reference_positions($exon_padded_start, $exon_padded_end);
      foreach my $genomic_align ( @{ $restricted_gab->get_all_non_reference_genomic_aligns() } ) {
        my $genomic_align_slice = $genomic_align->get_Slice();
        push(@{$exon_slices},$genomic_align_slice);
        say "GAS NAME: ".$genomic_align_slice->name;
        say "GAS SEQ: ".$genomic_align_slice->seq;
      }
    }

    $exon_align_slices->{$exon->dbID} = $exon_slices;
  }
  $self->param('_exons',\@unique_translateable_exons);
  $self->param('_exon_align_slices',$exon_align_slices);
}

sub run {
  my ($self) = @_;

  my $exons = $self->param('_exons');
  my $projected_exons = {};
  my $fail_count = 0;
  foreach my $exon (@{$exons}) {
    my $projected_exon = $self->project_exon($exon);
    if($projected_exon) {
      $projected_exons->{$exon->dbID} = $projected_exon;
    } else {
      say "Failed to project exon: ".$exon->stable_id;
      $fail_count++;
    }
  }

  my $transcripts = $self->build_transcripts($projected_exons);
  unless($transcripts) {
    
  }



  say "Had a total of ".$fail_count."/".scalar(@{$exons})." failed exon projections";
}


sub write_output {
  my ($self) = @_;

  my $gene_adaptor = $self->hrdb_get_con('target_transcript_db')->get_GeneAdaptor;
  my $slice_adaptor = $self->hrdb_get_con('target_transcript_db')->get_SliceAdaptor;

  my $genes = $self->output_genes();
  foreach my $gene (@{$genes}) {
    empty_Gene($gene);
    $gene_adaptor->store($gene);
  }

}

sub build_transcripts {
  my ($self,$projected_exons) = @_;

  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => 'cesar',
                                              -module => 'HiveCesar',
                                            );

  say "Building transcripts from projected exons";

  my $gene = $self->parent_gene;
  say "Source gene SID: ".$gene->stable_id;
  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}) {
    unless($transcript->biotype eq 'protein_coding') {
      next;
    }

    say "Source transcript SID: ".$transcript->stable_id;
    my $exons = $transcript->get_all_translateable_Exons();
    my $projected_exon_set = [];
    foreach my $exon (@{$exons}) {
      say "Checking for exon ".$exon->stable_id;
      if($projected_exons->{$exon->dbID}) {
        push(@{$projected_exon_set},$projected_exons->{$exon->dbID});
      }
    }

    unless(scalar(@{$projected_exon_set}) > 0) {
      next;
    }

    my $transcript_slice = ${$projected_exon_set}[0]->slice->seq_region_Slice();
    my $transcript_analysis = ${$projected_exon_set}[0]->analysis;
    my $projected_transcript = Bio::EnsEMBL::Transcript->new(-exons => $projected_exon_set,
                                                             -analysis  => $analysis,
                                                             -stable_id => $transcript->stable_id.".".$transcript->version,
                                                             -slice     => $transcript_slice);
    say "Transcript SID: ".$projected_transcript->stable_id;
    say "Transcript start: ".$transcript->seq_region_start;
    say "Transcript end: ".$transcript->seq_region_end;
    say "Transcript exon count: ".scalar(@{$projected_exon_set});

    my $start_exon = ${$projected_exon_set}[0];
    my $end_exon = ${$projected_exon_set}[$#{$projected_exon_set}];
    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($start_exon);
    $translation->start(1);
    $translation->end_Exon($end_exon);
    $translation->end($end_exon->length());
    $projected_transcript->translation($translation);

    # Set the phases
    calculate_exon_phases($projected_transcript, 0);

    say "Transcript translation:\n".$transcript->translation->seq;
    say "Projected transcript translation:\n".$projected_transcript->translation->seq;

    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->add_Transcript($projected_transcript);
    $gene->analysis($analysis);
    $self->output_genes($gene);
  }

}
sub project_exon {
  my ($self,$exon) = @_;
  my $exon_align_slices = $self->param('_exon_align_slices')->{$exon->dbID};

  unless($exon_align_slices) {
    return 0;
  }

  my $seq = $exon->seq->seq();
  my $phase = $exon->phase();
  my $end_phase = $exon->end_phase();

  my $phase_string;
  my $end_phase_string;
  my $start_coord;
  my $end_coord;

  say "Phase/End phase: ".$phase."/".$end_phase;
  say "S0: ".$seq;
  if($phase == 0 || $phase == -1) {
  } elsif($phase == 1) {
    my $split_codon = substr($seq,0,2);
    say "Split coding base start: ".lc($split_codon);
    $seq = lc($split_codon).substr($seq,2);
  } elsif($phase == 2) {
    my $split_codon = substr($seq,0,1);
    say "Split coding base start: ".lc($split_codon);
    $seq = lc($split_codon).substr($seq,1);
  } else {
    $self->throw("Unexpected phase found for exon ".$exon->stable_id." (".$exon->dbID()."): ".$phase);
  }

  say "S1: ".$seq;
  if($end_phase == 0 || $end_phase == -1) {
  } elsif($end_phase == 1) {
    my $split_codon = substr($seq,length($seq)-1);
    say "Split coding base end: ".lc($split_codon);
    $seq = substr($seq,0,length($seq)-1).lc($split_codon);
  } elsif($end_phase == 2) {
    my $split_codon = substr($seq,length($seq)-2);
    say "Split coding base end: ".lc($split_codon);
    $seq = substr($seq,0,length($seq)-2).lc($split_codon);
  } else {
    $self->throw("Unexpected end phase found for exon ".$exon->stable_id." (".$exon->dbID()."): ".$end_phase);
  }
  say "S2: ".$seq;
  my $rand = int(rand(10000));
  # Note as each accession will occur in only one file, there should be no problem using the first one
  my $outfile_path = $self->param('output_path')."/cesar_".$$."_".$exon->stable_id."_".$rand.".fasta";

  open(OUT,">".$outfile_path);
  say OUT ">".$exon->stable_id;
  say OUT $seq;

  foreach my $exon_align_slice (@{$exon_align_slices}) {
    say $exon->stable_id.": ".$exon_align_slice->name();
    say OUT ">".$exon_align_slice->name();
    say OUT $exon_align_slice->seq();
  }
  close OUT;

  chdir $self->param('cesar_path');
  my $extra_commands = "";
  if($exon->{'is_first_exon'}) {
   $extra_commands .= "--is_first_exon 1 ";
  }

  if($exon->{'is_last_exon'}) {
   $extra_commands .= "--is_last_exon 1 ";
  }

  my $cesar_command = "python2.7 CESAR/CESAR.py ".$outfile_path." ".$extra_commands."--clade human > ".$outfile_path.".ces";
  if(system($cesar_command)) {
    return(0);
  }

  my $projected_exon = $self->parse_exon($exon,$outfile_path.".ces");
  if($projected_exon) {
   return($projected_exon);
  } else {
    return(0);
  }
}

sub parse_exon {
  my ($self,$source_exon,$projected_outfile_path) = @_;
  open(IN,$projected_outfile_path);
  my @projection_array = <IN>;
  close IN;

  if(scalar(@projection_array) > 4) {
    $self->warning("Output file has more than one projection. We should put in code for this. Exon: ".$source_exon->stable_id);
  }

  my $slice_name = $projection_array[2];
  my $proj_seq = $projection_array[3];

  # Code needs to go here to tag frameshifts

  # For now just sub out the frameshifts
  $proj_seq =~ s/-//g;

  unless($proj_seq =~ /(^[atgcn]+)[ATGCN]+([atgcn]+$)/) {
    $self->throw("Couldn't match on the projected exon seq. Sequence: ".$proj_seq);
  }

  my $flank_length_start = length($1);
  my $flank_length_end = length($2);
#  >chromosome:Mmul_8.0.1:20:2231611:2231829:1
  unless($slice_name =~ /^>(.+\:.+\:.+\:)(.+)\:(.+)\:(.+)$/) {
    $self->throw("Couldn't parse the header to get the slice name. Header: ".$slice_name);
  }

  my $proj_exon_slice_name = $1;
  my $start_coord = $2 + $flank_length_start;
  my $end_coord = $3 - $flank_length_end;
  my $strand = $4;
  $proj_exon_slice_name .= join(":",($start_coord,$end_coord,$strand));
  my $slice_adaptor =  $self->hrdb_get_con('target_dna_db')->get_SliceAdaptor();
#  say "FM2 DUMPER: ".Dumper($slice_adaptor);

  my $proj_slice = $slice_adaptor->fetch_by_name($proj_exon_slice_name)->seq_region_Slice;
  unless($proj_slice) {
    $self->throw("Couldn't retrieve a slice for: ".$proj_exon_slice_name);
  }
#  say "Proj slice name: ".$proj_slice->name;
  my $proj_exon = new Bio::EnsEMBL::Exon(
      -START     => $start_coord,
      -END       => $end_coord,
      -STRAND    => $strand,
      -SLICE     => $proj_slice,
      -ANALYSIS  => $source_exon->analysis,
      -STABLE_ID => $source_exon->stable_id.".".$source_exon->version,
      -VERSION   => 1,
    );

  say "Proj exon slice: ".$proj_exon->slice->name;
  say "Proj exon seq: ".$proj_exon->seq->seq;

  return($proj_exon);

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
    for(my $i=0; $i<scalar(@{$exons}); $i++) { # my $exon (@{$exons}) {
      my $exon = ${$exons}[$i];
      if($i==0) {
        $exon->{'is_first_exon'} = 1;
      }

      if($i==scalar(@{$exons})-1) {
        $exon->{'is_last_exon'} = 1;
      }

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

sub parent_gene {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_parent_gene',$val);
  }

  return($self->param('_parent_gene'));
}

sub output_genes {
  my ($self,$val) = @_;
  unless($self->param('_output_genes')) {
    $self->param('_output_genes',[]);
  }

  if($val) {
    push(@{$self->param('_output_genes')},$val);
  }

  return($self->param('_output_genes'));
}

1;
