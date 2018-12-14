=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the
EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar

=cut

=head1 DESCRIPTION

HiveCesar fetches the genes corresponding to the given array of gene_id or a single gene_id in source_db,
projects its exons based on the given Compara lastz alignment and the
CESAR2.0 aligner and builds single-transcript genes from these
projections to be written to target_db while filtering out the specified
transcripts by applying the filter in TRANSCRIPT_FILTER.

#not done, under review
#seqEdits are inserted in order to deal with alignment gaps

Selenocysteine attributes are inserted in order to deal
with seleno-like TGA stops which will be converted
into NNN triplets to make the best possible aligment.

=head1 OPTIONS

-iid                  gene_id or array of gene_id from the source_db corresponding to the
gene to be projected from the source_dna_db to the target_dna_db.
-output_path          Path where the output files will be stored.
-source_dna_db        Ensembl database containing the DNA sequences that
correspond to the input gene_id from the source_db.
-target_dna_db        Ensembl database containing the DNA sequences
corresponding to the target_db species where the gene_id will be
projected to.
-source_db            Ensembl database containing the genes whose
transcripts will be projected to the target species db target_db.
-target_db            Ensembl database containing the DNA sequences
corresponding to the target species the input gene_id will be projected to.
-compara_db           Compara database containing the lastz alignments
of the source and target species.
-method_link_type     Default to 'LASTZ_NET' so it works with the
Compara lastz alignments.
-cesar_path           Path to the directory containing the CESAR2.0
binary to be run (excluding the binary filename).
-canonical            If set to 1, then only the canonical transcript for each gene will be fetched from the source db.
-canonical_or_longest If set to 1, then only the canonical transcript for each gene will be projected. If the projection is not done successfully, the next transcript having the longest translation will be projected until there is a successful projection.
-common_slice         If set to 1, all the transcripts projected from the same gene will be put on the same slice (and gene) based on the most common seq region name and min and max coordinates covering them. The projected transcripts on the other slices will be discarded. If set to 0 (default), the projected transcripts will be used to make single-transcript genes.
-stops2introns        Number of stops within a translation which will be replaces with introns. Default 0.
-max_stops            Only the transcripts whose translations contain a number of stops equal to or less than max_stops will be stored.  Default 0 (translations with stops are not allowed by default).
-TRANSCRIPT_FILTER    Hash containing the parameters required to apply
to the projected transcript to exclude some of them. Default to
ExonerateTranscriptFilter pid,cov 50,50 although note that the actual
implementation of this filter allows pid,cov below 50,50 in some cases.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar;

use warnings;
use strict;
use feature 'say';
use Scalar::Util 'reftype';
use List::Util qw[min max];
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;
use Bio::EnsEMBL::Analysis::Tools::ClusterFilter;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils
qw(replace_stops_with_introns                                                                     
   calculate_exon_phases                                                                  
   set_alignment_supporting_features                                                                  
   features_overlap);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(align_proteins);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(replace_stops_with_introns);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
    return {
      iid => '',
      output_path => '',
      source_dna_db => '',
      target_dna_db => '',
      source_db => '',
      target_db => '',
      compara_db => '',
      method_link_type => 'LASTZ_NET',
      transcript_region_padding => 50,
      cesar_path => '',
      canonical => 0,
      canonical_or_longest => 0,
      common_slice => 0,
      stops2introns => 0,
      max_stops => 0,
      #TRANSCRIPT_FILTER => {
      #                       OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
      #                       PARAMETERS => {
      #                         -coverage => 50,
      #                         -percent_id => 50,
      #                       },
      #                     }
   }
}

sub fetch_input {
  my($self) = @_;

  unless(-e $self->param('output_path')) {
    system("mkdir -p ".$self->param('output_path'));
  }

  my @input_id = ();
  if (reftype($self->param('iid')) eq "ARRAY") {
    @input_id = @{$self->param('iid')};
  } else {
    # make single-element array
    @input_id = ($self->param('iid'));
  }

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

  # Get the genome db adaptor
  my $genome_dba = $compara_dba->get_GenomeDBAdaptor();

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
  if ($@) {
    $self->throw("Had trouble fetching coord systems for ".
                 $source_genome_db->assembly." and ".$target_genome_db->assembly.
                 " from core dbs:\n".$@);
  }
 
  #########
  # get the compara data: MethodLinkSpeciesSet, reference DnaFrag,
  # and all GenomicAlignBlocks
  #########
  my $mlss = $compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_by_method_link_type_GenomeDBs($self->param('method_link_type'),[$source_genome_db,$target_genome_db]);

  if (!($mlss)) {
    $self->throw("No MethodLinkSpeciesSet for :\n".$self->param('method_link_type')."\n".$source_species."\n".$target_species);
  }

  foreach my $ii (@input_id) {

    my $gene = $source_transcript_dba->get_GeneAdaptor->fetch_by_dbID($ii);
    my @unique_translateable_transcripts = $self->get_unique_translateable_transcripts($gene,$self->param('canonical'));
    my $transcript_align_slices;
    my $genomic_align_block_adaptor = $compara_dba->get_GenomicAlignBlockAdaptor();
    my $transcript_region_padding = $self->param('transcript_region_padding');

    foreach my $transcript (@unique_translateable_transcripts) {
     
      my $transcript_group_id_lengths = {};
      my $transcript_group_id_min_starts = {};
      my $transcript_group_id_max_ends = {};
      my $transcript_group_id_seq_region_names = {};
      my $transcript_group_id_seq_region_strands = {};
   
      my $transcript_padded_start = $transcript->start()-$transcript_region_padding;
      if ($transcript_padded_start < 0) {
        $transcript_padded_start = 0;
      }

      my $transcript_padded_end = $transcript->end()+$transcript_region_padding;
      if ($transcript_padded_end > $gene->slice()->length()) {
        $transcript_padded_end = $gene->slice()->length();
      }

      my $slice_adaptor = $source_dna_dba->get_SliceAdaptor();
      my $transcript_slice = $slice_adaptor->fetch_by_region($transcript->slice()->coord_system_name(),$transcript->slice()->seq_region_name(),$transcript_padded_start,$transcript_padded_end,$transcript->seq_region_strand());

      say "---transcript slice: ".$transcript_slice->coord_system_name()." ".$transcript_slice->name()."\n"."length of transcript slice seq: ".length($transcript_slice->seq());

      my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss,$transcript_slice);
      my $transcript_slices = [];

      foreach my $genomic_align_block (@{$genomic_align_blocks}) {
        my $gab = $genomic_align_block->restrict_between_reference_positions($transcript_padded_start,$transcript_padded_end);
        if ($gab) {
          my $gab_group_id = $gab->group_id();
          foreach my $genomic_align (@{$gab->get_all_non_reference_genomic_aligns()}) {
            my $genomic_align_slice = $genomic_align->get_Slice();
            $transcript_group_id_lengths->{$gab_group_id} += length($genomic_align_slice->seq());
            
            if (!($transcript_group_id_min_starts->{$gab_group_id})) {
              $transcript_group_id_min_starts->{$gab_group_id} = $genomic_align_slice->start();
            } else {
              $transcript_group_id_min_starts->{$gab_group_id} = min($transcript_group_id_min_starts->{$gab_group_id},
                                                                     $genomic_align_slice->start());
            }
            $transcript_group_id_max_ends->{$gab_group_id} = max($transcript_group_id_max_ends->{$gab_group_id},
                                                                 $genomic_align_slice->end());
            $transcript_group_id_seq_region_names->{$gab_group_id} = $genomic_align_slice->seq_region_name();
            $transcript_group_id_seq_region_strands->{$gab_group_id} = $genomic_align_slice->strand();

            say "GAS NAME: ".$genomic_align_slice->name();
            say "GAS START: ".$genomic_align_slice->start();
            say "GAS END: ".$genomic_align_slice->end();
            say "GAS SEQ length: ".length($genomic_align_slice->seq());         
          } 
        }
      }

      my @sorted_group_ids = sort {$$transcript_group_id_lengths{$a} <=> $$transcript_group_id_lengths{$b}} keys %{$transcript_group_id_lengths};
      my $longest_group_id = $sorted_group_ids[-1];
     
      if ($longest_group_id) {

        print "longest group is: ".$longest_group_id."\n";
        print "length: ".$transcript_group_id_lengths->{$longest_group_id}."\n";
      
        my $sa = $self->hrdb_get_con('target_dna_db')->get_SliceAdaptor();       
        my $target_transcript_slice = $sa->fetch_by_region(undef,
                                                  $transcript_group_id_seq_region_names->{$longest_group_id},
                                                  $transcript_group_id_min_starts->{$longest_group_id},
                                                  $transcript_group_id_max_ends->{$longest_group_id},
                                                  $transcript_group_id_seq_region_strands->{$longest_group_id});
        
        #if ($transcript->length() <= $target_transcript_slice->length()) {
          $transcript_align_slices->{$transcript->dbID()} = $target_transcript_slice;
        #} else {
        #  $self->warning("Skipping transcript ".$transcript->dbID()."( ".$transcript->stable_id()." ) because its length (sum of exons length) is greater than the target transcript slice length.");
        #}
      }
    }
    
    if (\@unique_translateable_transcripts and $transcript_align_slices and $gene) {
      $self->parent_genes($gene);
      $self->unique_translateable_transcripts(\@unique_translateable_transcripts);
      $self->transcript_align_slices($transcript_align_slices);
    } else {
      $self->warning("Gene ".$gene->dbID()."( ".$gene->stable_id()." ) does not have unique_translateable_transcripts or transcript_align_slices.");
    }

  } # foreach my $ii

  # check that each gene has a set of unique translateable transcripts and transcript align slices
  if (scalar(@{$self->parent_genes()}) != scalar(@{$self->unique_translateable_transcripts()}) or
      scalar(@{$self->parent_genes()}) != scalar(@{$self->transcript_align_slices()})) {
    $self->throw("Different number of elements in parent_genes, unique_translateable_transcripts and transcript_align_slices arrays.");
  }
}

sub run {
  my ($self) = @_;

  my $gene_index = 0;
  foreach my $gene (@{$self->parent_genes()}) {

    my $transcripts = @{$self->unique_translateable_transcripts()}[$gene_index];

    my @projected_transcripts = ();
    my $fail_count = 0;

    my $himem_required = 0;

    foreach my $transcript (@{$transcripts}) {

      my $projected_transcript = $self->project_transcript($transcript,$gene_index);
      if ($projected_transcript == -1) {
        # it will be retried in the himem analysis
        say "Failed to project transcript due to himem required: ".$transcript->stable_id();
        $fail_count++;
        $himem_required = 1;
        last;
      } elsif ($projected_transcript) {
        push(@projected_transcripts,$projected_transcript);
      } else {
        say "Failed to project transcript: ".$transcript->stable_id();
        $fail_count++;
      }
    }

    if (!$himem_required) {
      $self->build_gene(\@projected_transcripts,$gene_index,$self->param('canonical'),$self->param('canonical_or_longest'));
    }
    
    say "Had a total of ".$fail_count."/".scalar(@{$transcripts})." failed transcript projections for gene ".$gene->dbID();
    $gene_index++;
  }
}


sub write_output {
  my ($self) = @_;

  my $gene_adaptor = $self->hrdb_get_con('target_transcript_db')->get_GeneAdaptor;
  my $slice_adaptor = $self->hrdb_get_con('target_transcript_db')->get_SliceAdaptor;

  my $genes = $self->output_genes();
  foreach my $gene (@{$genes}) {
    my $transcript = @{$gene->get_all_Transcripts}[0]; # any transcript
    if (!($gene_adaptor->fetch_by_transcript_stable_id($transcript->stable_id()))) {
      say "Storing gene: ".$gene->start.":".$gene->end.":".$gene->strand." (g.start:g.end:g.strand). Transcript stable ID used to fetch gene: ".$transcript->stable_id();
      empty_Gene($gene);
      $gene->biotype('projection');
      $gene_adaptor->store($gene);
    } else {
      say "NOT storing gene because it has already been stored: ".$gene->start.":".$gene->end.":".$gene->strand."(g.start,g.end,g.strand). Transcript stable ID used to fetch gene: ".$transcript->stable_id();
    }
  }
}

sub build_gene {
  my ($self,$projected_transcripts,$gene_index,$canonical,$canonical_or_longest) = @_;

  if (scalar(@$projected_transcripts) > 0) {
    my $analysis = Bio::EnsEMBL::Analysis->new(
                                                -logic_name => 'cesar',
                                                -module => 'HiveCesar',
                                              );

    say "Building genes from projected transcripts";

    my $gene = @{$self->parent_genes}[$gene_index];
    say "Source gene SID: ".$gene->stable_id();

    my @projected_transcripts_for_gene = ();  
    if ($canonical_or_longest) {
   
      # sort transcripts by translation length
      my @projected_transcripts_sorted = sort {$b->translation()->length() <=> $a->translation()->length()} @$projected_transcripts;
    
      # give maximum priority to the canonical transcript by swapping it to the first element
      my $t_index = 0;
      my $num_transcripts = scalar(@projected_transcripts_sorted);
      my $canonical_t_sid = $gene->canonical_transcript()->stable_id();
    
      while ($t_index < $num_transcripts) {
        my $curr_t = $projected_transcripts_sorted[$t_index];
        if ($curr_t->stable_id() eq $canonical_t_sid) {
           # canonical transcript projection found
          unshift @projected_transcripts_for_gene,$curr_t;
        } else {
          push(@projected_transcripts_for_gene,$curr_t);
        }
        $t_index++;
      }

    } else {
      @projected_transcripts_for_gene = @$projected_transcripts;
    }

    my $projected_transcripts_on_common_slice;
    if ($self->param('common_slice')) {
      $projected_transcripts_on_common_slice = $self->set_common_slice(\@projected_transcripts_for_gene);
    }

    if ($projected_transcripts_on_common_slice) {
      @projected_transcripts_for_gene = @$projected_transcripts_on_common_slice;
    }

    #my $projected_gene = $gene->flush_Transcripts();
    my $projected_gene = Bio::EnsEMBL::Gene->new();
    $projected_gene->stable_id($gene->stable_id());

TRANSCRIPT: foreach my $projected_transcript (@projected_transcripts_for_gene) {
      # do not store transcripts containing stops
      if ($projected_transcript->translate()) {
        my $projected_transcript_translate_seq = $projected_transcript->translate()->seq();
        my $num_stops = $projected_transcript_translate_seq =~ s/\*/\*/g;
        if ($num_stops > $self->param('max_stops')) {
          say "The projected transcript has been filtered out because its translation contains more than the maximum number of stops (".$num_stops." stops, the maximum is ".$self->param('max_stops')." stops).";
        } else {
          # filter out transcripts below given pid and cov
          if ($self->TRANSCRIPT_FILTER) {
            if (scalar(@{$projected_transcript->get_all_supporting_features()}) > 0) {
              my $filtered_transcripts = $self->filter->filter_results([$projected_transcript]);
              if (scalar(@$filtered_transcripts) > 0) {
               
                if ($self->param('common_slice')) {
                  # only one projected gene per source gene is built
                  $projected_gene->add_Transcript($projected_transcript);
                  $projected_gene->analysis($analysis);
                  $self->output_genes($projected_gene);
                  if ($canonical_or_longest) {
                    # only one projected transcript per projected gene
                    last TRANSCRIPT;
                  }

                } else {
                  # multiple projected single-transcript genes per source gene is built
                  $self->output_single_transcript_gene($projected_transcript,$analysis);
                  if ($canonical_or_longest) {
                    # only one projected gene per source gene is built
                    last TRANSCRIPT;
                  }
                }
                
              } else {
                say "The projected transcript has been filtered out because its pid and cov are too low.";
              }
            }
          } else {
            if ($self->param('common_slice')) {
              # only one projected gene per source gene is built
              $projected_gene->add_Transcript($projected_transcript);
              $projected_gene->analysis($analysis);
              $self->output_genes($projected_gene);
              if ($canonical_or_longest) {
                # only one projected gene per source gene is built
                last TRANSCRIPT;
              }
            } else {
              # multiple projected single-transcript genes per source gene is built
              $self->output_single_transcript_gene($projected_transcript,$analysis);
              if ($canonical_or_longest) {
                # only one projected gene per source gene is built
                last TRANSCRIPT;
              }
            } # end else common_slice
          } # end else TRANSCRIPT_FILTER
        } # end else num_stops > 0
      } # end if project_transcript->translate
    } # end foreach TRANSCRIPT
  }
}

sub largest_value_mem {
# it returns the key containing the largest value in a given hash
  my $hash = shift;
  my ($key,@keys) = keys %$hash;
  my ($big,@vals) = values %$hash;

  for (0 .. $#keys) {
    if ($vals[$_] > $big) {
      $big = $vals[$_];
      $key = $keys[$_];
    }
  }
  $key
}

sub set_common_slice {
# it sets the same slice for all transcripts on the same seq region
# so they can be added to the same gene later
# based on the most common seq region name
# and minimum and maximum transcript coordinates for that seq region name
# It returns a new array containing the transcripts on the same slice only
# after having discarded the transcripts on other seq regions
  my ($self,$projected_transcripts) = @_;

  my @projected_transcripts_on_common_slice = ();
  my %common_regions;
  my $min = 9999999999999999;
  my $max = 0;
  my $sa = $self->hrdb_get_con('target_dna_db')->get_SliceAdaptor(); 

  foreach my $projected_transcript (@$projected_transcripts) {
    $common_regions{$projected_transcript->seq_region_name()} += 1;
  }

  my $most_common_seq_region_name = largest_value_mem(\%common_regions);

  foreach my $projected_transcript (@$projected_transcripts) {
    if ($projected_transcript->seq_region_name() eq $most_common_seq_region_name) {
      $min = min($projected_transcript->seq_region_start(),$min);
      $max = max($projected_transcript->seq_region_end(),$max);
      push(@projected_transcripts_on_common_slice,$projected_transcript);
    }
  }

  my $common_slice = $sa->fetch_by_region(undef,$most_common_seq_region_name,$min,$max);

  foreach my $projected_transcript (@$projected_transcripts) {
    $projected_transcript->slice($common_slice);
  }

  return \@projected_transcripts_on_common_slice;
}

sub project_transcript {
  my ($self,$transcript,$gene_index) = @_;

  my $transcript_align_slice = @{$self->transcript_align_slices()}[$gene_index]->{$transcript->dbID()};  

  if (!$transcript_align_slice) {
    $self->warning("transcript_align_slice is empty for transcript dbID ".$transcript->dbID()." (stable ID ".$transcript->stable_id_version()." ) gene index: ".$gene_index);
    return 0;
  }

  # set output filename
  my $rand = int(rand(10000));
  # Note as each accession will occur in only one file, there should be no problem using the first one
  my $outfile_path = $self->param('output_path')."/cesar_".$$."_".$transcript->stable_id()."_".$rand.".fasta";
  $self->files_to_delete($outfile_path);

  open(OUT,">".$outfile_path);
  my $exon_index = 0;
EXON:  foreach my $exon (@{$transcript->get_all_translateable_Exons()}) {
    my $seq = $exon->seq->seq();
    my $phase = $exon->phase();
    my $end_phase = $exon->end_phase();

    my $start_coord;
    my $end_coord;

    # Find 5' split codon and lower case bases
    if ($phase == 0 or $phase == -1) {
      ;
    } elsif($phase == 1) {
      my $split_codon = substr($seq,0,2);
      say "Split coding base start: ".lc($split_codon);
      $seq = lc($split_codon).substr($seq,2);
      $transcript->{$exon_index}->{'five_split_codon'} = $split_codon;
    } elsif($phase == 2) {
      my $split_codon = substr($seq,0,1);
      say "Split coding base start: ".lc($split_codon);
      $seq = lc($split_codon).substr($seq,1);
      $transcript->{$exon_index}->{'five_split_codon'} = $split_codon;
    } else {
      $self->throw("Unexpected phase found for exon ".$exon->stable_id." (".$exon->dbID()."): ".$phase);
    }

    # Find 3' split codon and lower case bases
    if ($end_phase == 0 or $end_phase == -1) {
      ;
    } elsif($end_phase == 1) {
      my $split_codon = substr($seq,length($seq)-1);
      say "Split coding base end: ".lc($split_codon);
      $seq = substr($seq,0,length($seq)-1).lc($split_codon);
      $transcript->{$exon_index}->{'three_split_codon'} = $split_codon;
    } elsif($end_phase == 2) {
      my $split_codon = substr($seq,length($seq)-2);
      say "Split coding base end: ".lc($split_codon);
      $seq = substr($seq,0,length($seq)-2).lc($split_codon);
      $transcript->{$exon_index}->{'three_split_codon'} = $split_codon;
    } else {
      $self->throw("Unexpected end phase found for exon ".$exon->stable_id." (".$exon->dbID()."): ".$end_phase);
    }

    # remove bases from the 3' end in case the sequence is not multiple of 3
    while (($seq =~ tr/ACGTN//)%3 != 0) {
      $seq = substr($seq,0,length($seq)-1);
      say("Removed last base because the end phase is -1 and the sequence is not multiple of 3.");
    }

    $exon_index++; # the exon index is only used to record the split codons for each exon

    # replace TGA stops/selenocysteines with NNN so CESAR2.0 makes it match with anything
    my $i_step = 1;
    for (my $i = 0; $i < length($seq); $i += $i_step) {
      my $base_1 = substr($seq,$i,1);
      if ($base_1 !~ /[acgtn]/) {
        # we have reached the first (upper case or -) base of the exon sequence
        $i_step = 3;
      }
      if ($i_step == 3) {
        my $base_2 = substr($seq,$i+1,1);
        my $base_3 = substr($seq,$i+2,1);
          if ($base_1 eq "T" and
              $base_2 eq "G" and
              $base_3 eq "A" and
              $i+$i_step < length($seq)) { # ignore the last stop codon
            # selenocysteine stop found needs to be replaced with cysteine
            $seq = substr($seq,0,$i)."NNN".substr($seq,$i+3);
            $exon->{'selenocysteine'} = $i; # create new exon attribute to store the start of the selenocysteine

            $self->warning("Potential selenocysteine/TGA stop codon found at position $i (including lower case flanks). Exon ".$exon->stable_id().". Sequence (including lower case flanks): $seq");
          } elsif ( ($base_1 eq "T" and
                     $base_2 eq "A" and
                     ($base_3 eq "A" or $base_3 eq "G")) and
                    $i+$i_step < length($seq) ) { # ignore the last stop codon
            $self->warning("Stop codon TAA or TAG found in reference. Exon ".$exon->stable_id()." skipped.");
            next EXON;
          }
       }
    }
    
    # replace any base different from A,C,G,T with N
    $seq =~ tr/ykwmsrdvhbxYKWMSRDVHBX/nnnnnnnnnnnNNNNNNNNNNN/;

    say OUT ">".$transcript->stable_id()."_".$exon->stable_id();
    say OUT $seq;
  }

  # CESAR2.0 requires references and queries to be separated by a line starting with '#'.
  # References are the exons (together with their reading frame) that you want to align to the query sequence.
  say OUT "#";

  say $transcript->stable_id().": ".$transcript_align_slice->name();
  say OUT ">".$transcript_align_slice->name();
    
  my $transcript_align_slice_seq = $transcript_align_slice->seq();
    
  # replace any base different from A,C,G,T with N
  $transcript_align_slice_seq =~ tr/ykwmsrdvhbxRYKWMSDVHBX/nnnnnnnnnnnNNNNNNNNNNN/;
  say OUT $transcript_align_slice_seq;

  close OUT;

  chdir $self->param('cesar_path');
  my $cesar_command = $self->param('cesar_path')."/cesar ".$outfile_path." --clade human ";
  if ($self->param('cesar_mem')) {
    $cesar_command .= "--max-memory ".$self->param('cesar_mem'); # set max mem in GB
  }

  say $cesar_command;

  my $cesar_output;
  $cesar_output = `$cesar_command 2>&1`;
  my $fces_name_tmp = $outfile_path.".ces.tmp";
  my $fces_name = $outfile_path.".ces";

  if ($cesar_output =~ /Your attempt requires ([0-9]+)\./) {
    # Parse error message from CESAR2.0 to retry the job according to the memory required:
    # CRITICAL src/Cesar.c:117 main():  The memory consumption is limited to 20.0000 GB by default. Your attempt requires 73.6664 GB. You can change the limit via --max-memory.
    my $output_hash = {};
    push(@{$output_hash->{'iid'}},@{$self->parent_genes()}[$gene_index]->dbID());

    if ($1 < 15) {
      $self->dataflow_output_id($output_hash,15);
    } elsif ($1 < 25) {
      $self->dataflow_output_id($output_hash,25);
    } elsif ($1 < 35) {
      $self->dataflow_output_id($output_hash,35);
    } else {
      $self->dataflow_output_id($output_hash,55);
    }
    
    $self->warning("cesar command FAILED and it will be passed to cesar_XXX: ".$cesar_command.". Gene ID: ".@{$self->parent_genes()}[$gene_index]->dbID()." CESAR output: ".$cesar_output."\n");
    say "projected exon will return -1";
    return (-1);
  } elsif ($cesar_output =~ /CRITICAL/) {
    $self->throw("cesar command FAILED: ".$cesar_command."\n");
  } else {
    open(FCES,'>',$fces_name_tmp) or die $!;
    print FCES $cesar_output;
    close(FCES);
    system("grep -v WARNING $fces_name_tmp > $fces_name"); # remove CESAR2.0 warnings
  }

  $self->files_to_delete($fces_name_tmp);
  $self->files_to_delete($fces_name);
  my $projected_transcript = $self->parse_transcript($transcript,$fces_name);
#  while(my $file_to_delete = shift(@{$self->files_to_delete})) {
#    system('rm '.$file_to_delete);
#  }

  if ($projected_transcript) {
    return ($projected_transcript);
  } else {
    return (0);
  }
}

sub parse_transcript {
  my ($self,$source_transcript,$projected_outfile_path) = @_;

  open(IN,$projected_outfile_path);
  my @projection_array = <IN>;
  close IN;
  
  # remove last line if blank and not corresponding to the last sequence
  if ($projection_array[-1] =~ /^\$/ and $projection_array[-3] =~ /^>/) {
    pop(@projection_array);
  }
  
  my $reference_exon_header = shift(@projection_array);
  my $source_seq =  shift(@projection_array);
  my $slice_name = shift(@projection_array);
  my $proj_seq = shift(@projection_array);

  # CESAR sometimes produces projected exon sequences with no actual exonic sequence like 
  # >reference
  # agACACATaa
  # >projected
  # tg------aa
  # which we are going to discard.
  my $num_proj_seq_exonic_bases = $proj_seq =~ tr/ACGTNYKWMSRDVHBX//;

  if (scalar(@projection_array) > 0) {
    $self->throw("Output file has more than one projection. The projection having fewer gaps will be chosen. Transcript: ".$source_transcript->stable_id());
  }

  chomp($source_seq);
  chomp($proj_seq);
  
  if ($slice_name !~ /^>(.+\:.+\:.+\:)(.+)\:(.+)\:(.+)$/) {
    $self->throw("Couldn't parse the header to get the slice name. Header: ".$slice_name);
  }

  my $proj_transcript_slice_name = $1;
  my $transcript_start_coord = $2;
  my $transcript_end_coord = $3;
  my $original_proj_transcript_strand = $4;
  my $strand = $original_proj_transcript_strand;

  $source_seq =~ /( *)([\-atgcnATGCN> ]+[-atgcnATGCN>]+)( *)/; # '>' means do not expect a splice site in the query because the intron has been deleted, annotate as one composite exon
  
  my $transcript_left_flank = $1;
  my $source_transcript_align = $2;
  my $transcript_right_flank = $3;
  my $exon_offset_from_start = 0;
  
  $exon_offset_from_start = length($transcript_left_flank);

  $proj_transcript_slice_name .= join(":",($transcript_start_coord,$transcript_end_coord,$original_proj_transcript_strand));
  my $slice_adaptor = $self->hrdb_get_con('target_dna_db')->get_SliceAdaptor();

  my $proj_transcript_slice = $slice_adaptor->fetch_by_name($proj_transcript_slice_name);

  if (!($proj_transcript_slice)) {
    $self->throw("Couldn't retrieve a slice for transcript: ".$proj_transcript_slice_name);
  }

  # I have to store the source sequence exons as they appear in the alignment file
  # so I can compare the split codons from the source and the projected sequence later on.
  my @source_exons = ();
  my $exon_index = 0;
  while ($source_seq =~ /([atgcn]*)([\-ATGCN>]+)([atgcn]*)/g) { # source exon sequences contain split codon bases as lower case bases

    # Find split codons in the alignment file as CESAR might have changed them compared to the ones in the fasta file
    $source_transcript->{$exon_index}->{'alignment_five_split_codon'} = $1;
    $source_transcript->{$exon_index}->{'alignment_three_split_codon'} = $3;

    push(@source_exons,$1.$2.$3);
    $exon_index++;
  }

  my @projected_exons = ();
  my $source_exon_index = 0;
  my $accum_proj_seq_gap_length = 0;

PROJSEQ: while ($proj_seq =~ /([\-ATGCN]+)/g) {
    
    my $proj_exon_sequence = $1;

    # @- and @+ are the start and end positions of the last match.
    # $-[0] and $+[0] are the entire pattern.
    # $-[N] and $+[N] are the $N submatches.
    my $exon_start = $-[0]+1; # +1 because exon coordinates start at 1 for Exon objects
    my $exon_end = $+[0]+1-1; # +1 because exon coordinates start at 1 for Exon objects
                          # -1 because $+[] gives the index of the character following the match, not the last character of the match.

    if ($proj_exon_sequence =~ /^\-+\-+$/) {
      # skip projected exon sequences which only contain '-'
      say "Skipping projected exon sequence because it only contains '-'. File: ".$projected_outfile_path;
      $source_exon_index = 0;
      next PROJSEQ;
    } elsif ($proj_exon_sequence =~ /(^\-*)[ATGCNatgcn]+(\-*$)/) {
      # ignore any number of '-' which would change the phase at the beginning and at the end of the sequence
      
      if (length($1) % 3 == 1) {
        # ignore 1 '-' at the beginning to keep the phase
        say "Ignoring 1 '-' at the beginning of the exon to keep the phase. $proj_transcript_slice_name File: ".$projected_outfile_path;
        $source_exons[$source_exon_index] = substr($source_exons[$source_exon_index],1);
        $proj_exon_sequence = substr($proj_exon_sequence,1);
        $exon_start += 1;
      } elsif (length($1) % 3 == 2) {
        # ignore 2 '-' at the beginning to keep the phase
        say "Ignoring 2 '-' at the beginning of the exon to keep the phase. $proj_transcript_slice_name File: ".$projected_outfile_path;
        $source_exons[$source_exon_index] = substr($source_exons[$source_exon_index],2);
        $proj_exon_sequence = substr($proj_exon_sequence,2);
        $exon_start += 2;
      }

      if (length($2) % 3 == 1) {
        # ignore 1 '-' at the end to keep the phase
        say "Ignoring 1 '-' at the end of the exon to keep the phase. $proj_transcript_slice_name File: ".$projected_outfile_path;
        $source_exons[$source_exon_index] = substr($source_exons[$source_exon_index],0,-1);
        $proj_exon_sequence = substr($proj_exon_sequence,0,-1);
        $exon_end -= 1;
      } elsif (length($2) % 3 == 2) {
        # ignore 2 '-' at the end to keep the phase
        say "Ignoring 2 '-' at the beginning of the exon to keep the phase. $proj_transcript_slice_name File: ".$projected_outfile_path;
        $source_exons[$source_exon_index] = substr($source_exons[$source_exon_index],0,-2);
        $proj_exon_sequence = substr($proj_exon_sequence,0,-2);
        $exon_end -= 2;
      }
    }

    $accum_proj_seq_gap_length = substr($proj_seq,0,$exon_start-1) =~ tr/\-//;
    $exon_start -= $accum_proj_seq_gap_length;
                          
    my $base_index_offset = $exon_start;
    my $original_exon_start = $exon_start;
    my $still_first = 1; # boolean to indicate whether we are still going through a gap before
                         # finding the first valid codon within an exon

    my $end_gap_length = 0;
    if ($proj_exon_sequence =~ /[ACGT]*(\-+)$/) {
      $end_gap_length = length($1);
    }

    $accum_proj_seq_gap_length = substr($proj_seq,0,$exon_end) =~ tr/\-//;
    $exon_end -= $accum_proj_seq_gap_length;

    # find the source exon sequence which corresponds to the current projected sequence
    my $source_exon_sequence = $source_exons[$source_exon_index];
    while (length($proj_exon_sequence) != length($source_exon_sequence) and
           $source_exon_index < @source_exons) {
      say "Source exon sequence and projected exon sequence lengths IN THE ALIGNMENT (including '-') do not match. Source exon index: ".$source_exon_index.". Skipping source exon. File: ".$projected_outfile_path;
      $source_exon_index++;
      $source_exon_sequence = $source_exons[$source_exon_index];
    }

    #say "proj exon sequence length: ".length($proj_exon_sequence);
    #say "source exon sequence length: ".length($source_exon_sequence);
   
    if (length($proj_exon_sequence) != length($source_exon_sequence)) {
      $self->warning("Source exon sequence not found for current projected sequence ".$proj_exon_sequence.". Projected exon not added to the projected transcript. File: ".$projected_outfile_path);
      $source_exon_index = 0;
      next PROJSEQ;
    } else {
      # source exon sequence found for current projected sequence
      my @source_exon_sequence_arr = split('',$source_exon_sequence);
      my @proj_exon_sequence_arr = split('',$proj_exon_sequence);
      #say "source exon sequence: ".$source_exon_sequence;
      #say "proj_exon_sequence exon sequence: ".$proj_exon_sequence;
      
      # get the split codon at the 3' end so we can use its length to know
      # when to stop looping through the codons
      my $source_split_codon_3 = "";
      my $source_split_codon_3_length = 0;
      if ($source_transcript->{$source_exon_index}) {
        if ($source_transcript->{$source_exon_index}->{'alignment_three_split_codon'}) {
          $source_split_codon_3 = $source_transcript->{$source_exon_index}->{'alignment_three_split_codon'};
          $source_split_codon_3_length = length($source_split_codon_3);
        }
      }

      my $current_codons_in_gap_count = 0;
      my $exon_made = 0;
      my $exon_made_dash_count = 0;
      my $split_codon_exon_made = 0;
      my $previous_exon_end = 0;
      my $current_proj_seq_gap_length = 0;
      my $base_index = 0;
      CODON: while ($base_index < @source_exon_sequence_arr-$end_gap_length) {

        if ($source_exon_sequence_arr[$base_index] eq ">") {
          # there is a fixed-length gap of arbitrary length of 19 bases like ">" in the source and "-" in the projected
          # sequence inserted by CESAR to represent a merged/composite exon.
          # The gap has to be skipped since it would change the phase otherwise.
          $base_index += 19;
          say "Merged/composite exon fixed-length gap of 19 '>' found and skipped.";
          next CODON;
        }

        if ($base_index == 0) {
          # first codon in the exon
          if ($source_transcript->{$source_exon_index}) { # if the source exon contained any split codon
            if ($source_transcript->{$source_exon_index}->{'alignment_five_split_codon'}) { # if the source exon contained a split codon at the 5' end
              my $source_split_codon_5_length = length($source_transcript->{$source_exon_index}->{'alignment_five_split_codon'});
              $base_index += $source_split_codon_5_length;
              
              my $next_codon_string = $proj_exon_sequence_arr[$base_index].$proj_exon_sequence_arr[$base_index+1].$proj_exon_sequence_arr[$base_index+2];
              my $dash_count = $next_codon_string =~ tr/\-//;
              if ($dash_count) {
                say "New exon made after the first codon because there is a gap after the split codon.";
                $split_codon_exon_made = 1;
                push(@projected_exons,
                     new Bio::EnsEMBL::Exon(-START     => $exon_start,#+$exon_offset_from_start,#+$proj_transcript_slice->start(),
                                            -END       => $exon_start-1+$source_split_codon_5_length,#-1,#+$exon_offset_from_start,#+$proj_transcript_slice->start(), # $+[] gives the index of the character following the match, not the last character of the match.
                                            -STRAND    => 1, # the proj_transcript_slice is already on the reverse strand
                                            -SLICE     => $proj_transcript_slice,
                                            -ANALYSIS  => $source_transcript->analysis(),
                                            -STABLE_ID => $source_transcript->stable_id_version(),
                                            -VERSION   => 1));
               
                $exon_start += $source_split_codon_5_length; # exon_start will be ready for next new exon within the current exon or to be greater than the end meaning no more exons should be made
              }
              $still_first = 0;
              next CODON;
            }
          }

          my $codon_string = $proj_exon_sequence_arr[$base_index].$proj_exon_sequence_arr[$base_index+1].$proj_exon_sequence_arr[$base_index+2];

          my $dash_count = $codon_string =~ tr/\-//;
          my $current_proj_seq_gap_length = substr($proj_exon_sequence,0,$base_index+1-1) =~ tr/\-//;

          if ($dash_count) {
            # skip the codon to keep the phase

            if ($current_codons_in_gap_count == 0) {
             
              if ($dash_count == 1 or $dash_count == 2) {
                ;
              }
            } else {
              # going through a gap in the projected sequence
              $exon_start += (3-$dash_count); # any base in the gap which is not forming a complete 3-base codon
                                              # will be skipped to be part of an intron
              $current_codons_in_gap_count++;
            }
          } else { # no more dashes so there is no gap
                    
            if ($current_codons_in_gap_count) {
              # the gap just finished 
              $current_codons_in_gap_count = 0;
              
              my $previous_codon_string = $proj_exon_sequence_arr[$base_index-3].$proj_exon_sequence_arr[$base_index-2].$proj_exon_sequence_arr[$base_index-1];
              my $previous_dash_count = $previous_codon_string =~ tr/\-//;
              
              if (!$exon_made and ($previous_dash_count == 1 or $previous_dash_count == 2)) {
                say "New exon made at the end of the gap because no exon was made and codon dash count is 1 or 2.";
                push(@projected_exons,
                     new Bio::EnsEMBL::Exon(-START     => $exon_start,
                                            -END       => $exon_end,
                                            -STRAND    => 1, # the proj_transcript_slice is already on the reverse strand
                                            -SLICE     => $proj_transcript_slice,
                                            -ANALYSIS  => $source_transcript->analysis(),
                                            -STABLE_ID => $source_transcript->stable_id_version(),
                                            -VERSION   => 1));
                $exon_start = $base_index_offset+$base_index-$current_proj_seq_gap_length;
                
                $original_exon_start = $exon_start;
              } elsif ($exon_made and ($previous_dash_count == 1 or $previous_dash_count == 2)) {
                say "Not making exon at the end of the gap because the exon was already made at the beginning.";           
              } else {
                # exon not made neither at the beginning nor at the end of the gap because it is multiple of 3
                $exon_end = $previous_exon_end;
                say "Not making exon neither at the beginning nor at the end of the gap because it is multiple of 3."; 
              }            
              $exon_made = 0;
              $exon_made_dash_count = 0;
            }
          }
          
        } elsif ($base_index+3 >= length($source_exon_sequence)-$end_gap_length) {
          # last codon in the exon

          $exon_end = $base_index_offset+$base_index;
          
          my $last_codon_length = length($source_exon_sequence)-$end_gap_length-$base_index;
          if ($last_codon_length % 3 != $source_split_codon_3_length) {
            $exon_end += $source_split_codon_3_length;
          } else {
            $exon_end += $last_codon_length-1; 
          }
          my $current_proj_seq_gap_length = substr($proj_exon_sequence,0,$base_index+1-1) =~ tr/\-//;
          $exon_end -= $current_proj_seq_gap_length;

        } else {
          # any other codon between the first one and the last one
          my $codon_string = $proj_exon_sequence_arr[$base_index].$proj_exon_sequence_arr[$base_index+1].$proj_exon_sequence_arr[$base_index+2];
          my $dash_count = $codon_string =~ tr/\-//;
          my $current_proj_seq_gap_length = substr($proj_exon_sequence,0,$base_index+1-1) =~ tr/\-//;
          if ($dash_count) {
            # skip the codon to keep the phase

            if ($current_codons_in_gap_count == 0) {
              # if this is the first codon in the gap and the number of '-' is not 3, then make a new exon

              # we keep what the exon_end should be in case we need to make a new exon
              # we won't need to make a new exon if all the gap codons have 3 '-'
              # we will make a new exon if the first codon has 1 or 2 '-' or
              # the last codon has 1 or 2 '-'
              $previous_exon_end = $exon_end;
              $exon_end = $base_index_offset+$base_index-1-$current_proj_seq_gap_length;
              
              if (($dash_count == 1 or $dash_count == 2) and
                  (!$still_first) and (!$split_codon_exon_made)) {
                say "New exon made at the beginning of the gap because dash count is 1 or 2.";
                push(@projected_exons,
                     new Bio::EnsEMBL::Exon(-START     => $exon_start,
                                            -END       => $exon_end,
                                            -STRAND    => 1, # the proj_transcript_slice is already on the reverse strand
                                            -SLICE     => $proj_transcript_slice,
                                            -ANALYSIS  => $source_transcript->analysis(),
                                            -STABLE_ID => $source_transcript->stable_id_version(),
                                            -VERSION   => 1));
                $exon_made = 1;
                $exon_made_dash_count = $dash_count;
                $exon_start = $exon_end+4-$dash_count;
                $original_exon_start = $exon_start;

              } elsif ($still_first) {
                $exon_start += (3-$dash_count); # any base in the gap which is not forming a complete 3-base codon
                                                # will be skipped to be part of an intron
              }
              $current_codons_in_gap_count++;

            } else {
              # going through a gap in the projected sequence
              $exon_start += (3-$dash_count); # any base in the gap which is not forming a complete 3-base codon
                                                # will be skipped to be part of an intron
                             
              $current_codons_in_gap_count++;
            }
          } else { # no more dashes so there is no gap
            
            $still_first = 0;
               
            if ($current_codons_in_gap_count) {
              # the gap just finished 
              $current_codons_in_gap_count = 0;
              
              my $previous_codon_string = $proj_exon_sequence_arr[$base_index-3].$proj_exon_sequence_arr[$base_index-2].$proj_exon_sequence_arr[$base_index-1];
              my $previous_dash_count = $previous_codon_string =~ tr/\-//;
              
              if (!$exon_made and ($previous_dash_count == 1 or $previous_dash_count == 2) and $still_first) {
                # exon_start <= exon_end means that we are skipping a gap at the beginning of the exon so we haven't reached
                # the first actual codon yet, so no need to make any new exon in this case

                say "New exon made at the end of the gap because no exon was made and dash count is 1 or 2.";
                push(@projected_exons,
                     new Bio::EnsEMBL::Exon(-START     => $exon_start,#+$exon_offset_from_start,#+$proj_transcript_slice->start(),
                                            -END       => $exon_end,#-1,#+$exon_offset_from_start,#+$proj_transcript_slice->start(), # $+[] gives the index of the character following the match, not the last character of the match.
                                            -STRAND    => 1, # the proj_transcript_slice is already on the reverse strand
                                            -SLICE     => $proj_transcript_slice,
                                            -ANALYSIS  => $source_transcript->analysis(),
                                            -STABLE_ID => $source_transcript->stable_id_version(),
                                            -VERSION   => 1));
                $exon_start = $base_index_offset+$base_index-$current_proj_seq_gap_length;               
                $original_exon_start = $exon_start;
              } elsif ($exon_made and ($previous_dash_count == 1 or $previous_dash_count == 2)) {
                say "Not making exon at the end of the gap because the exon was already made at the beginning.";           
              } else {
                # exon not made neither at the beginning nor at the end of the gap because it is multiple of 3
                $exon_end = $previous_exon_end;
                say "Not making exon neither at the beginning nor at the end of the gap because it is multiple of 3."; 
              }
              $exon_made = 0;
              $exon_made_dash_count = 0;
            }
          }
        }
        $base_index += 3;
      } # end while codon
    } # end else length($proj_exon_sequence) != length($source_exon_sequence

    if ($exon_start <= $exon_end) {
      push(@projected_exons,
        new Bio::EnsEMBL::Exon(-START     => $exon_start,
                               -END       => $exon_end,
                               -STRAND    => 1, # the proj_transcript_slice is already on the reverse strand
                               -SLICE     => $proj_transcript_slice,
                               -ANALYSIS  => $source_transcript->analysis(),
                               -STABLE_ID => $source_transcript->stable_id_version(),
                               -VERSION   => 1));
    } else {
      say "exon_start is greater than exon_end, not making final exon at the end. This should follow a single-codon exon formed by a split codon. Source transcript: ".$source_transcript->stable_id_version()." Projected exon sequence: ".$proj_seq;
    }

    $source_exon_index++;
  } # end while PROJSEQ

  if (scalar(@projected_exons) > 0) {

    my $projected_transcript = Bio::EnsEMBL::Transcript->new(-exons => \@projected_exons,
                                                             -analysis => $source_transcript->analysis(),
                                                             -stable_id => $source_transcript->stable_id_version(),
                                                             -strand => 1,
                                                             -slice => $proj_transcript_slice);

    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($projected_exons[0]);
    $translation->start(1);
    $translation->end_Exon($projected_exons[-1]);
    $translation->end($projected_exons[-1]->length());
    $projected_transcript->translation($translation);

    # Set the phases  
    calculate_exon_phases($projected_transcript,$source_transcript->translation()->start_Exon()->phase());

    if ($self->param('stops2introns') > 0 and $projected_transcript->translate() and $projected_transcript->translate()->seq()) {

      my $projected_transcript_translate_seq = $projected_transcript->translate()->seq();
      my $num_stops = $projected_transcript_translate_seq =~ s/\*/\*/g;
      say $projected_transcript->stable_id()." : number of stops before replacing up to ".$self->param('stops2introns')." stops with introns: ".$num_stops;
      
      my $projected_transcript_after_replaced_stops = replace_stops_with_introns($projected_transcript,$self->param('stops2introns'));
      # if there are more stops than the stops2introns value then 'replace_stops_with_introns'
      # returns 0 so we want to use the original projected_transcript
      if ($projected_transcript_after_replaced_stops and $projected_transcript_after_replaced_stops->translate()) {
        $projected_transcript = $projected_transcript_after_replaced_stops;
      }

      $projected_transcript_translate_seq = $projected_transcript->translate()->seq();
      $num_stops = $projected_transcript_translate_seq =~ s/\*/\*/g;
      say $projected_transcript->stable_id()." : number of stops after replacing up to ".$self->param('stops2introns')." stops with introns: ".$num_stops;
    }

    if ($projected_transcript) {
      # Set the exon and transcript supporting features
      if ($projected_transcript->translation()->seq()) { 
        set_alignment_supporting_features($projected_transcript,$source_transcript->translation()->seq(),$projected_transcript->translation()->seq());
      }

      say "Transcript translation:\n".$source_transcript->translation()->seq();
      say "Projected transcript translation:\n".$projected_transcript->translation()->seq();

      my ($coverage,$percent_id) = (0,0);
      if ($projected_transcript->translation()->seq()) {
        ($coverage,$percent_id) = align_proteins($source_transcript->translate()->seq(),$projected_transcript->translate()->seq());
      }
      $projected_transcript->source($coverage);
      $projected_transcript->biotype($percent_id);
      $projected_transcript->description("stable_id of source: ".$source_transcript->stable_id());

      # add a 'seq_edits' attribute to the proj_exon object
      # to store the seq edits that will be added to the transcript
      # when the transcript is built
      #my @seq_edits = make_seq_edits($source_seq,$proj_seq);
      #$proj_exon->{'seq_edits'} = \@seq_edits;

      return ($projected_transcript);
    } else { # replace_stops_with_introns returns zero if translation contains stop codon adjacent to gap
      $self->warning("Transcript with internal stop codon next to gapped sequence. Source transcript ".$source_transcript->stable_id()." skipped.");
      return 0;
    }
  } else {
    # no exons projected after parsing
    $self->warning("No exons projected after parsing. Source transcript ".$source_transcript->stable_id()." skipped.");
    return 0;
  }
}

sub make_seq_edits {
  # It returns an array of SeqEdit objects for the target sequence to make
  # the insertions for the alignment gaps between the source and target sequences
  # created for an alignment between two dna sequence in cesar output format ie string containing acgtACGT-.
  # A SeqEdit object is added to the array for each substring of any number of "-" not multiple of 3.
  # Inserted bases are taken from the source sequence.

  my ($source_seq,$target_seq) = @_;

  my @seq_edits = ();
  my $acumm_gap_length = 0;
 
  # count the number of lowercase bases before the start of the actual (uppercase) target sequence
  my $num_lowercase_left_flank = 0;
  my $target_seq_copy = $target_seq;
  if ($target_seq_copy =~ m/([acgtn]+)[ACGTN-]+/g) {
    $num_lowercase_left_flank = length($1);
  }

  while ($target_seq =~ /(\-+)/g) {
    $acumm_gap_length += length($1);
    my $start = pos($target_seq)+1-$acumm_gap_length-$num_lowercase_left_flank;
    my $end = $start-1;

    push(@seq_edits,Bio::EnsEMBL::SeqEdit->new(-CODE    => '_rna_edit',
                                               -NAME    => 'rna_edit',
                                               -DESCRIPTION    => 'Cesar alignment',
                                               -START   => $start,
                                               -END     => $end,
                                               -ALT_SEQ => substr($source_seq,pos($target_seq)-length($1),length($1))
                                              ));
  }
  return (@seq_edits);
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

sub get_unique_translateable_transcripts {
  my ($self,$gene,$canonical) = @_;

  my $translateable_transcripts = {};
  if ($gene->biotype ne 'protein_coding') {
    $self->input_job->autoflow(0);
    $self->complete_early('Gene does not have protein_coding biotype!');
  }

  my @transcripts = ();
  if ($canonical) {
    push(@transcripts,$gene->canonical_transcript());
  } else {
    @transcripts = @{$gene->get_all_Transcripts()};
  }
  foreach my $transcript (@transcripts) {
    if ($transcript->biotype eq 'protein_coding') {
      my $dbid = $transcript->dbID();
      $translateable_transcripts->{$dbid} = $transcript;
    }
  }
  return(values(%{$translateable_transcripts}));
}

sub make_alignment_mapper {
  my ($self,$gen_al_blocks) = @_;

  my $FROM_CS_NAME = 'chromosome';
  my $TO_CS_NAME   = 'scaffold';

  my $mapper = Bio::EnsEMBL::Mapper->new($FROM_CS_NAME,
                                         $TO_CS_NAME);

  foreach my $bl (@$gen_al_blocks) {
    foreach my $ugbl (@{$bl->get_all_ungapped_GenomicAlignBlocks}) {
      my ($from_bl) = $ugbl->reference_genomic_align;
      my ($to_bl)   = @{$ugbl->get_all_non_reference_genomic_aligns};

      $mapper->add_map_coordinates($from_bl->dnafrag->name,
                                   $from_bl->dnafrag_start,
                                   $from_bl->dnafrag_end,
                                   $from_bl->dnafrag_strand*$to_bl->dnafrag_strand,
                                   $to_bl->dnafrag->name,
                                   $to_bl->dnafrag_start,
                                   $to_bl->dnafrag_end);
    }
  }

  return $mapper;
}

sub parent_genes {
  my ($self,$val) = @_;
  if (!($self->param('_parent_genes'))) {
    $self->param('_parent_genes',[]);
  }

  if ($val) {
    push(@{$self->param('_parent_genes')},$val);
  }

  return $self->param('_parent_genes');
}

sub unique_translateable_transcripts {
  my ($self,$val) = @_;
  if (!($self->param('_unique_translateable_transcripts'))) {
    $self->param('_unique_translateable_transcripts',[]);
  }

  if ($val) {
    push(@{$self->param('_unique_translateable_transcripts')},$val);
  }

  return($self->param('_unique_translateable_transcripts'));
}

sub transcript_align_slices {
  my ($self,$val) = @_;
  if (!($self->param('_transcript_align_slices'))) {
    $self->param('_transcript_align_slices',[]);
  }

  if ($val) {
    push(@{$self->param('_transcript_align_slices')},$val);
  }

  return($self->param('_transcript_align_slices'));
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

sub files_to_delete {
  my ($self,$val) = @_;
  unless($self->param('_files_to_delete')) {
    $self->param('_files_to_delete',[]);
  }

  if($val) {
    push(@{$self->param('_files_to_delete')},$val);
  }

  return($self->param('_files_to_delete'));
}

sub remove_overlapping_exons {
# any exon overlapped by another longer exon is removed
# and it will not be part of the returned array reference of exons
  my ($exons) = shift;

  print("Removing overlapping projected exons... Before: ".scalar(@$exons)." exons.\n");

  my @discarded_exon_indexes = ();
  my $exon1_index = 0;
 
  foreach my $exon1 (@$exons) {
    my $exon2_index = 0;
    foreach my $exon2 (@$exons) {
      if ($exon1_index != $exon2_index and !($exon2_index ~~ @discarded_exon_indexes)) {
        if (features_overlap($exon1,$exon2)) {
          if ($exon1->length() <= $exon2->length()) {
            push(@discarded_exon_indexes,$exon1_index);
            last;
          }
        }
      }
      $exon2_index++;
    }
    $exon1_index++;
  }

  my $no_overlap_exons = [];
  my $exon_index = 0;
EXON: foreach my $exon (@$exons) {
    foreach my $discarded_exon_index (@discarded_exon_indexes) {
      if ($exon_index == $discarded_exon_index) {
        $exon_index++;
        next EXON;
      }
    }
    push(@{$no_overlap_exons},$exon);
    $exon_index++;
  }

  print("Removing overlapping projected exons... After: ".scalar(@{$no_overlap_exons})." exons.\n");

  return $no_overlap_exons;
}

####################################
# config variable holders
####################################
#
# transcript editing and filtering
#

sub TRANSCRIPT_FILTER {
   my ($self, $val) = @_;

  if (defined $val) {
    $self->param('TRANSCRIPT_FILTER',$val);
  }

  if ($self->param_is_defined('TRANSCRIPT_FILTER')) {
    return $self->param('TRANSCRIPT_FILTER');
  }
  else {
    return;
  }
}

sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->param('_runnable_filter',$val);
  }

  # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
  if ($self->TRANSCRIPT_FILTER and !$self->param_is_defined('_runnable_filter')) {
    if (not ref($self->TRANSCRIPT_FILTER) eq "HASH" or
        not exists($self->TRANSCRIPT_FILTER->{OBJECT}) or
        not exists($self->TRANSCRIPT_FILTER->{PARAMETERS})) {

      $self->throw("FILTER in config for '".$self->analysis->logic_name."' must be a hash ref with elements:\n" .
            "  OBJECT : qualified name of the filter module;\n" .
            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
    } else {
      $self->require_module($self->TRANSCRIPT_FILTER->{OBJECT});
     
$self->filter($self->TRANSCRIPT_FILTER->{OBJECT}->new(%{$self->TRANSCRIPT_FILTER->{PARAMETERS}}));
    }
  }
  if ($self->param_is_defined('_runnable_filter')) {
    return $self->param('_runnable_filter');
  }
  else {
    return;
  }
}

sub output_single_transcript_gene {
  my ($self,$projected_transcript,$analysis) = @_;

  my $gene = Bio::EnsEMBL::Gene->new();
  $gene->add_Transcript($projected_transcript);
  $gene->analysis($analysis);
  $self->output_genes($gene);
}

1;
