=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionMinimap

=head1 SYNOPSIS

To be used as part of an eHive pipeline config file:

{
  -logic_name => 'minimap_project_transcripts',
  -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionMinimap',
  -parameters => {
                   'logic_name' => 'minimap_projection',
                   'module' => 'Minimap2',
                   'source_dna_db' => $self->default_options()->{'projection_source_db'},
                   'target_dna_db' => $self->o('dna_db'),
                   'source_db' => $self->o('projection_source_db'),
                   'target_db' => $self->o('projection_db'),
                   'compara_db' => $self->o('projection_compara_db'),
                   'method_link_type' => $self->o('method_link_type'),
                   'minimap_path' => $self->o('minimap2_path'),
                   'paftools_path' => $self->o('paftools_path'),
                   'minimap_coverage' => 80,
                   'minimap_percent_id' => 60,
                 },
  -rc_name    => 'default',
  -hive_capacity => 900,
  -flow_into => {
                  -3 => ['failed_coding_jobs'],
                },
},

=head1 DESCRIPTION

HiveProjectionMinimap fetches the transcripts corresponding to the given array of transcript_id or a single transcript_id in source_db,
projects them based on the given Compara lastz alignment and the Minimap2 aligner
and builds single-transcript genes from these projections to be written to target_db while filtering out the specified
transcripts if they don't meet the threshold criteria set for the protein muscle alignments (if applicable)
by coverage and percentage identity.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionMinimap;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(attach_Analysis_to_Gene attach_Slice_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(attach_Slice_to_Transcript empty_Transcript);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (align_proteins write_seqfile write_sliceseq2fastafile);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Runnable::Minimap2;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

 Arg [1]    : None
 Description: It fetches the source db transcript and its corresponding target db slices from the compara db
              to make the Minimap2 runnables.
 Returntype : None
 Exceptions : Throws if the toplevel coordinate systems could not be fetched from the source and target transcript dbs or
              if the MethodLinkSpeciesSet could not be fetched from the compara db.

=cut

sub fetch_input {

  my($self) = @_;

  $self->create_analysis;
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
  $target_transcript_dba->dnadb($target_dna_dba);
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

  $self->param('_transcript_biotype',{});
  foreach my $input_id (@$input_ids) {

    my $transcript = $source_transcript_dba->get_TranscriptAdaptor->fetch_by_dbID($input_id);
    my $biotype = $transcript->biotype;
    my $stable_id = $transcript->stable_id.".".$transcript->version;
    $self->param('_transcript_biotype')->{$stable_id} = $biotype;
    say "Processing source transcript: ".$transcript->stable_id . "from: " . $source_genome_db->assembly . "\n";
    say "with input_id: $input_id and number of exons: " . scalar(@{$transcript->get_all_Exons()})  . "\n";     
    my $transcript_slices = $self->process_transcript($transcript,$compara_dba,$mlss,$source_genome_db,$source_transcript_dba);
    $self->make_runnables($transcript->spliced_seq(),$transcript_slices,$input_id,$target_transcript_dba);
    #$self->make_runnables($transcript->feature_Slice(), $transcript_slices, $input_id, $target_transcript_dba);
  } #close foreach input_id

}

=head2 run

 Arg [1]    : None
 Description: It runs each runnable previously created and stores their transcript ID attributes in the corresponding output attribute.
              If a source transcript has been projected to multiple regions, it selects the best one based on combined coverage and percent id.
 Returntype : None
 Exceptions : Throws if Minimap2 does not put 1 transcript in a projected gene, the coverage and percent ID cannot be fetched or are equal to 0,
              or there was no best transcript selected.

=cut

sub run {

  my ($self) = @_;

  $self->runnable_failed(0);
  foreach my $runnable (@{$self->runnable}) {
    eval {
      $runnable->run;
    };
    if ($@) {
      my $except = $@;
      $self->runnable_failed($runnable->{'_transcript_id'});
      $self->warning("Issue with running Minimap2, will dataflow input id on branch -3. Exception:\n".$except);
      $self->param('_branch_to_flow_on_fail', -3);
    } else {
      # Store original transcript id for realignment later. Should implement a cleaner solution at some point
      foreach my $output (@{$runnable->output}) {
        $output->{'_old_transcript_id'} = $runnable->{'_transcript_id'};
      }

      # If the transcript has been projected to multiple places then select the best one in terms of combined coverage and
      # percent identity but also any that fall within 5 percent of this value
      my $preliminary_genes = $runnable->output(); # note that minimap2 returns a reference to an array of genes
      my @preliminary_transcripts;
      foreach my $preliminary_gene (@{$preliminary_genes}) {
        my @transcripts = @{$preliminary_gene->get_all_Transcripts()};
        my $num_transcripts = scalar(@transcripts);
        if ($num_transcripts == 1) {
          my $preliminary_transcript = $transcripts[0];

          # minimap2 stores the percent id in the description column and the coverage in the version column
          # I need to store them in the transcript so they can be used in the select_best_transcripts sub
          $preliminary_transcript->version($preliminary_gene->version());
          $preliminary_transcript->description($preliminary_gene->description());

          # the source transcript id is stored in the _old_transcript_id key of the output gene
          $preliminary_transcript->{'_old_transcript_id'} = $preliminary_gene->{'_old_transcript_id'};

          push(@preliminary_transcripts,$transcripts[0]);
        } elsif ($num_transcripts > 1) {
          $self->throw("Minimap2 put more than 1 transcript in a gene for the source transcript ".$runnable->{'_transcript_id'});
        } elsif ($num_transcripts < 1) {
          $self->throw("Minimap2 put 0 transcripts in a gene for the source transcript ".$runnable->{'_transcript_id'});
        }
      }

      my $selected_transcripts;
      if (scalar(@preliminary_transcripts) < 1) {
        # consider it as a failed runnable as the transcript was not projected
        $self->runnable_failed($runnable->{'_transcript_id'});
        $self->warning("Issue with running Minimap2. Transcript ".$runnable->{'_transcript_id'}." did not get projected. Input id will be passed to branch -3.\n");
        $self->param('_branch_to_flow_on_fail', -3);
        next;
      } elsif (scalar(@preliminary_transcripts) == 1) {
        $selected_transcripts = \@preliminary_transcripts;
      } else {
        $selected_transcripts = $self->select_best_transcripts(\@preliminary_transcripts);
      }
      # foreach transcript check translation: 
	  # Returns the peptide translation of the exons as a Bio::Seq
      # foreach my $tran (@$selected_transcripts) {
      #  if ( $tran->translation() ) {
      #    my $pep = $tran->translate();
      #    print "DEBUG::Transcript ", $tran->stable_id(), " is protein_coding\n";
      #  } else {
      #    print "DEBUG::Transcript ", $tran->stable_id(), " is non-coding\n";
      #  }
      # }
      $self->output($selected_transcripts);
    }
  }

  return 1;
}

=head2 write_output

 Arg [1]    : None
 Description: It creates single-transcript genes for each projected gene in the output, sets biotypes and stores them in the target database.
 Returntype : None
 Exceptions : None

=cut

sub write_output {

  my ($self) = @_;
 
  my $adaptor = $self->hrdb_get_con('target_transcript_db')->get_GeneAdaptor;
  my $slice_adaptor = $self->hrdb_get_con('target_transcript_db')->get_SliceAdaptor;

  my @output = @{$self->output};
  my $analysis = $self->analysis;
  foreach my $transcript (@output){

    my $slice_id = $transcript->start_Exon->seqname;
    my $slice = $slice_adaptor->fetch_by_name($slice_id);
    my $biotype = $self->retrieve_biotype($transcript);
    $transcript->biotype($biotype);

    attach_Slice_to_Transcript($transcript, $slice);

    if ($self->filter_transcript($transcript)) {
      # consider it as a failed runnable as the transcript was filtered out
      $self->runnable_failed($transcript->{'_old_transcript_id'});
      $self->warning("The transcript was filtered out after projection due to cov, pid or stop codons. Transcript ".$transcript->{'_old_transcript_id'}." . Input id will be passed to branch -3.\n");
      $self->param('_branch_to_flow_on_fail', -3);
      next;
    }

    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->biotype($biotype);
    $gene->add_Transcript($transcript);
    $gene->slice($slice);
    attach_Analysis_to_Gene($gene, $analysis);
    $adaptor->store($gene);
  }
  my $output_hash = {};
  my $failure_branch_code = $self->param('_branch_to_flow_on_fail');
  my $failed_transcript_ids = $self->runnable_failed;
  if (scalar @$failed_transcript_ids ) {
    $output_hash->{'iid'} = $failed_transcript_ids;
    $self->dataflow_output_id($output_hash, $failure_branch_code);
  }

  $output_hash = {};
  $failure_branch_code = $self->param('_branch_to_flow_on_noalignblocks');
  $failed_transcript_ids = $self->transcripts_noalignblocks();
  if (scalar @$failed_transcript_ids ) {
    $output_hash->{'iid'} = $failed_transcript_ids;
    $self->dataflow_output_id($output_hash, $failure_branch_code);
  }

  return 1;
}

=head2 runnable_failed

 Arg [1]    : 
 Description: It updates the _runnable_failed array by emptying it if called without any parameter or by adding the element in the 'runnable_failed' parameter, which
              is a transcript id.
 Returntype : Array of transcript_id
 Exceptions : None

=cut

sub runnable_failed {
  my ($self,$runnable_failed) = @_;
  unless ($self->param_is_defined('_runnable_failed')) {
    $self->param('_runnable_failed',[]);
  }
  if ($runnable_failed) {
    push (@{$self->param('_runnable_failed')},$runnable_failed);
  }
  return ($self->param('_runnable_failed'));
}

=head2 transcripts_noalignblocks

 Arg [1]    : 
 Description: It updates the _transcripts_noalignblocks array by emptying it if called without any parameter or by adding the element in the 'transcripts_noalignblocks' parameter, which
              is a transcript id.
 Returntype : Array of transcript_id
 Exceptions : None

=cut

sub transcripts_noalignblocks {
  my ($self,$transcripts_noalignblocks) = @_;
  unless ($self->param_is_defined('_transcripts_noalignblocks')) {
    $self->param('_transcripts_noalignblocks',[]);
  }
  if ($transcripts_noalignblocks) {
    push (@{$self->param('_transcripts_noalignblocks')},$transcripts_noalignblocks);
  }
  return ($self->param('_transcripts_noalignblocks'));
}

=head2 process_transcript

 Arg [1]    : Bio::EnsEMBL::Transcript to be projected
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor corresponding to the compara db containing the align blocks between the projection source and projection target assemblies.
 Arg [3]    : Bio::EnsEMBL::Compara::DBSQL::MethodLinkSpeciesSetAdaptor corresponding to the compara db in Arg [2] 
 Arg [4]    : Bio::EnsEMBL::Compara::GenomeDB corresponding to the projection source database. 
 Arg [5]    : Bio::EnsEMBL::DBSQL::DBAdaptor corresponding to the projection source database.
 Description: It fetches all align blocks from the compara db for the given transcript and it creates slices by clustering them,
              which will be the projection target slices.
 Returntype : Arrayref of Slice (transcript slices)
 Exceptions : It throws if the align blocks could not be converted into transcript slices.

=cut

sub process_transcript {
  my ($self,$transcript,$compara_dba,$mlss,$source_genome_db,$source_transcript_dba) = @_;

  my $max_cluster_gap_length = $self->max_cluster_gap_length($transcript);
  say "Max align gap: ".$max_cluster_gap_length;
  my $all_target_genomic_aligns = [];
  my $exons = $transcript->get_all_Exons;
  my $exon_region_padding = $self->param('exon_region_padding');
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
    my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,$dna_fragments,$exon_padded_start,$exon_padded_end);

    foreach my $genomic_align_block (@$genomic_align_blocks) {
      push(@{$all_target_genomic_aligns},@{$genomic_align_block->get_all_non_reference_genomic_aligns});
    }
  }

  my $unique_target_genomic_aligns = $self->unique_genomic_aligns($all_target_genomic_aligns);
  my @sorted_target_genomic_aligns = sort { $a->get_Slice->seq_region_name cmp $b->get_Slice->seq_region_name ||
                                            $a->get_Slice->start <=> $b->get_Slice->start
                                          } @{$unique_target_genomic_aligns};

  unless(scalar(@sorted_target_genomic_aligns)) {
    #say "No align blocks so skipping transcript";
    # consider it as a failed runnable as the transcript was not projected
    $self->transcripts_noalignblocks($transcript->dbID());
    $self->warning("No align blocks so skipping transcript. Transcript ".$transcript->dbID()." . Input id will be passed to branch -4.\n");
    $self->param('_branch_to_flow_on_noalignblocks', -4);
    next;
  }

  my $transcript_slices = $self->make_cluster_slices(\@sorted_target_genomic_aligns,$max_cluster_gap_length);
  unless($transcript_slices) {
    $self->throw("The sorted align blocks were not converted into transcript slices, something went wrong");
  }

  foreach my $transcript_slice (@{$transcript_slices}) {
    say "Created transcript slice: ".$transcript_slice->name."\n";
  }
  return $transcript_slices;
}

=head2 make_runnables

 Arg [1]    : Bio::EnsEMBL::Slice corresponding to the slice where the source transcript lies on 
 Arg [2]    : Arrayref of Slice, projection target transcript slices
 Arg [3]    : Integer, transcript id of the source transcript to be projected 
 Arg [4]    : Bio::EnsEMBL::DBSQL::DBAdaptor corresponding to the projection target database.
 Description: It makes the Minimap2 runnables and stores the source transcript id as an attribute for each runnable object.
              It writes the source and target transcripts sequences into temporary FASTA files as required by the Minimap2 runnables. 
 Returntype : None
 Exceptions : 

=cut

sub make_runnables {
  my ($self,$transcript_seq,$transcript_slices,$input_id,$target_transcript_dba) = @_;
  my %parameters = %{$self->parameters_hash};
  my $source_sequence_fasta_file = $self->param('tmpdir')."/source_sequence_".$input_id;
  my $target_sequences_fasta_file = $self->param('tmpdir')."/target_sequences_".$input_id;

  my ($type,$assembly,$chrname,$target_genomic_start,$target_genomic_end,$step); 
  if (scalar(@{$transcript_slices}) == 1) {
  	my $transcript_slice = @{$transcript_slices}[0]; 
    ($type,$assembly,$chrname,$target_genomic_start,$target_genomic_end,$step) = split(":", $transcript_slice->name);
  } elsif (scalar(@{$transcript_slices}) > 1) {
  	
  } else {
    die("I expect one value only"); 
  }
  foreach my $transcript_slice (@{$transcript_slices}) {
    say "Created transcript slice: ".$transcript_slice->name."\n";
    ($type,$assembly,$chrname,$target_genomic_start,$target_genomic_end,$step) = split(":", $transcript_slice->name); 
  } 

  # dump the transcript sequence into a file which will be the input source for Minimap2
  #write_sliceseq2fastafile($transcript_seq,$source_sequence_fasta_file);
  #my $source_sequence_fasta_file = write_seqfile($transcript_seq);
  open(OUT,">".$source_sequence_fasta_file);
  say OUT ">".$input_id;
  say OUT $transcript_seq;
  close OUT;

  # dump transcript slices sequences into a file which will be the input target for Minimap2
  write_sliceseq2fastafile($transcript_slices,$target_sequences_fasta_file);

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis          => $self->analysis(),
       -program           => $self->param('minimap_path'),
       -paftools_path     => $self->param('paftools_path'),
#       -options           => $self->param('minimap2_options'),
       -genome_index      => $target_sequences_fasta_file,#$genome_index,
       -input_file        => $source_sequence_fasta_file, # input is different
       -database_adaptor  => $target_transcript_dba,
       -delete_input_file => 0, # only set this when creating ranged files, not when using the original input file
       -skip_introns_check       => 1,
       -add_offset               => $target_genomic_start - 1,
       -skip_compute_translation => 0, # If you turn that to 1, it will not return translation
       -percent_id_cutoff => $self->param('minimap_percent_id'),
       -coverage_cutoff => $self->param('minimap_coverage'),
       #-canonical_intron_cutoff => 0.8, this is the default one in the Minimap2 runnable
  );

  $runnable->{'_transcript_id'} = $input_id;
  $self->runnable($runnable);
}

=head2 max_cluster_gap_length

 Arg [1]    : Bio::EnsEMBL::Transcript 
 Description: This sub will loop through the introns and decide on how long the max allowed value should be.
              It will be used to decide when to break up clusters on the same seq_region.
              The max value multiplier will change based on how long the biggest intron is.
 Returntype : None
 Exceptions : None

=cut

sub max_cluster_gap_length {
  my ($self,$transcript) = @_;

  # Have 50 as a min default value to handle small breaks in single exon genes
  # kbillis increased this value to 200000 in order to avoid having small slices 
  # that will be hard to parse the results due to different offsets. 
  my $longest_intron = 200000; # 50
  my $introns = $transcript->get_all_Introns();
  unless(scalar(@{$introns}) > 0) {
    return($longest_intron);
  }

  foreach my $intron (@{$introns}) {
    if($intron->length > $longest_intron) {
      $longest_intron = $intron->length;
    }
  }

  # This is only an intial guess, using orthologs would allow more accurate estimations for these
  # values in future. It would make sense to have different values for different clades
  if($longest_intron <= 10000) {
    return(int($longest_intron * 2));
  } elsif($longest_intron <= 50000) {
    return(int($longest_intron * 1.5));
  } else {
    return(int($longest_intron * 1.25))
  }
}

=head2 make_cluster_slices

 Arg [1]    : ArrayRef of Bio::EnsEMBL::Compara::GenomicAlign, sorted by seq_region_name and start
 Arg [2]    : Integer, maximum cluster gap length  
 Description: It makes slices by clustering the GenomicAlign features.
 Returntype : ArrayRef of Bio::EnsEMBL::Slice
 Exceptions : None

=cut

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
      $previous_genomic_align = $current_genomic_align;
    }
  } # end foreach my $current_genomic_align

  # push the final slice
  my $slice = $slice_adaptor->fetch_by_region($previous_genomic_align->get_Slice->coord_system_name,
                                              $previous_genomic_align->get_Slice->seq_region_name, $cluster_start, $cluster_end);
  push(@{$cluster_slices},$slice);
  return($cluster_slices);
}


=head2 unique_genomic_aligns

 Arg [1]    : ArrayRef of Bio::EnsEMBL::Compara::GenomicAlign
 Description: It makes an array of unique genomic aligns by dbID.
 Returntype : ArrayRef of Bio::EnsEMBL::Compara::GenomicAlign
 Exceptions : None

=cut

sub unique_genomic_aligns {
  my ($self,$genomic_aligns) = @_;

  my $seen_id_hash = {};
  my $unique_genomic_aligns = [];
  foreach my $genomic_align (@{$genomic_aligns}) {
    unless($seen_id_hash->{$genomic_align->dbID}) {
      push(@{$unique_genomic_aligns},$genomic_align);
      $seen_id_hash->{$genomic_align->dbID} = 1;
    }
  }
  return($unique_genomic_aligns);
}

=head2 select_best_transcripts

 Arg [1]    : ArrayRef of Bio::EnsEMBL::Transcript
 Description: It makes an array of selected transcripts based on coverage and percentage id.
 Returntype : ArrayRef of Bio::EnsEMBL::Transcript
 Exceptions : It throws if the coverage and percentage id cannot be fetched or they are zero, or if there was no transcript selected.

=cut

sub select_best_transcripts {
  my ($self,$preliminary_transcripts) = @_;
  my $selected_transcripts = [];
  my $best_score = 0;
  foreach my $preliminary_transcript (@{$preliminary_transcripts}) {
    my $cov = $preliminary_transcript->version();     # this was propagated from minimap gene version, which is the coverage
    my $pid = $preliminary_transcript->description(); # this was propagated from minimap gene description, which is the percent id
    my $combined_score = $cov + $pid;
    if ($combined_score > $best_score) {
      $best_score = $combined_score;
    }
    $preliminary_transcript->{'combined_score'} = $combined_score;
  }

  unless($best_score) {
    $self->throw('Issue with calculating the best score from projection result set. This should not happen.');
  }
  # At this point we know the highest value in terms of combined percent identity and coverage. Now
  # we want all transcripts that fall within 5 percent of this value
  foreach my $preliminary_transcript (@{$preliminary_transcripts}) {
    my $combined_score = $preliminary_transcript->{'combined_score'};
    if($combined_score/$best_score >= 0.95) {
      push(@{$selected_transcripts},$preliminary_transcript);
    }
  }

  unless(scalar(@{$selected_transcripts})) {
    $self->throw("No transcripts selected in select_best_transcripts, something went wrong");
  }

  return($selected_transcripts);

}

=head2 filter_transcript

 Arg [1]    : ArrayRef of Bio::EnsEMBL::Transcript
 Description: It checks if the projected transcript in Arg [1] and its source transcript protein alignment meets the coverage and percentage id
              thresholds set as parameters.
 Returntype : Boolean: 1 if the transcript does not meet the minimum coverage and percentage id or it contains stop codons; 0 if the transcript meets the minimum
              coverage and percentage id or if the transcript does not translate.
 Exceptions : None

=cut

sub filter_transcript {
  my ($self,$transcript) = @_;

  my $transcript_coverage = $transcript->version();     # this was propagated from minimap gene version, which is the coverage
  my $transcript_identity = $transcript->description(); # this was propagated from minimap gene description, which is the percent id

  unless($transcript_identity >= $self->param_required('minimap_percent_id') && $transcript_coverage >= $self->param_required('minimap_coverage')) {
    print("Transcript failed coverage (".$self->param_required('minimap_coverage').") and/or percent id (".$self->param_required('minimap_percent_id').") filter, will not store. Transcript coverage and percent id are: ".$transcript_coverage." ".$transcript_identity."\n");
    return(1);
  }

  my $original_transcript = $self->hrdb_get_con('source_transcript_db')->get_TranscriptAdaptor->fetch_by_dbID($transcript->{'_old_transcript_id'});

  if ($transcript->translation() and $transcript->translation()->seq() and
      $original_transcript->translation() and $original_transcript->translation()->seq()) {

    if ($transcript->translation()->seq() =~ /\*/) {
      say "Transcript translation seq has stop, will not store";
      return 1;
    }

    my ($translation_coverage,$translation_identity) = align_proteins($original_transcript->translation->seq, $transcript->translation->seq);
    unless($translation_identity >= $self->param_required('minimap_percent_id') && $translation_coverage >= $self->param_required('minimap_coverage')) {
      print("Translation failed coverage (".$self->param_required('minimap_coverage').") and/or percent id (".$self->param_required('minimap_percent_id').") filter, will not store. Translation coverage and percent id are: ".$translation_coverage." ".$translation_identity."\n");
      return(1);
    }
  } # end if($transcript->translation)
  return(0);
}

=head2 retrieve_biotype

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: It fetches the biotype for the projected transcript in Arg [1] from its source transcript.
 Returntype : String, biotype
 Exceptions : It throws if the biotype could not be fetched.

=cut

sub retrieve_biotype {
  my ($self, $transcript) = @_;

  my $original_transcript = $self->hrdb_get_con('source_transcript_db')->get_TranscriptAdaptor->fetch_by_dbID($transcript->{'_old_transcript_id'});
  my $stable_id = $original_transcript->stable_id_version();
  my $biotype = $self->param('_transcript_biotype')->{$stable_id};
  unless($biotype) {
    $self->throw("Failed to retieve biotype for output transcript. Should have been set in fetch_input");
  }
  return($biotype);
}

1;
