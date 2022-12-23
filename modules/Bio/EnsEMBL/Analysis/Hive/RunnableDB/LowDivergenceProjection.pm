=head1 LICENSE

Copyright [2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::LowDivergenceProjection

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::LowDivergenceProjection;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use File::Spec::Functions qw(tmpdir catfile);
use POSIX;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs);
use Bio::EnsEMBL::Analysis::Runnable::LowDivergenceProjection;
#use Bio::EnsEMBL::Analysis::Tools::Utilities qw (align_nucleotide_seqs);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene clone_Gene);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    no_projection => 0,
    coverage_cutoff => 0.98,
    perc_id_cutoff => 0.99,
    extended_length_variation_cutoff => 0.1,
    anchor_coverage_cutoff => 0.95,
    anchor_perc_id_cutoff => 0.99,
  }
}


sub fetch_input {
  my($self) = @_;

  $self->setup_fasta_db;
  $self->create_analysis;
  $self->param_required('genome_index');

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->say_with_header($self->is_fasta_file_ready);
  if ($self->is_fasta_file_ready) {
    if ($self->param_is_defined('source_dna_fasta')) {
      $self->say_with_header($self->is_fasta_file_ready);
      $self->say_with_header($self->param('source_dna_fasta'));
      $source_dna_dba->get_SequenceAdaptor->fasta($self->param_required('source_dna_fasta'));
    }
    if ($self->param_is_defined('genome_index')) {
      my ($fasta) = $self->param('genome_index') =~ /(.*).mmi/;
      $self->say_with_header($fasta);
      $target_dna_dba->get_SequenceAdaptor->fasta($fasta);
    }
  }

  # Define the source and target gene dbs
  my $source_gene_dba = $self->hrdb_get_dba($self->param('source_gene_db'), $source_dna_dba);
  my $target_gene_dba = $self->get_database_by_name('target_gene_db', $target_dna_dba);
  $self->hrdb_set_con($target_gene_dba,'target_gene_db');

#  my $input_genes = $self->fetch_input_genes_by_id([53671,56952,59930,59295,59266,59293],$source_gene_dba);
  my $input_genes = $self->fetch_input_genes_by_id($self->param('iid'),$source_gene_dba);

  my $sequence_adaptor = $source_dna_dba->get_SequenceAdaptor();


  say "Processing ".scalar(@$input_genes)." genes";


################
# Put in code to sort input genes (might just be pre sorted by api, but just in case)
# Then for each gene record two genes to the left and right. Store these in a gene based hash
# Each key in the hash should be a gene id that points at up to 4 other gene ids that are keys on a hash
# Later on, once the final mapped set has been created and sorted, it will look at the current gene id
# and get the gene ids to the left and right and then compare them to the entries in the hash created now
# If we have two or more matches then the we can be confident in the gene's location
################

  my $sorted_input_genes = [sort { $a->slice->name() cmp $b->slice->name() or
                                   $a->start() <=> $b->start() or
                                   $a->end() <=> $b->end() }  @{$input_genes}];


#  my $batched_input_genes = $self->batch_input_genes();

  my $gene_synteny_hash = {};
  my $gene_genomic_seqs_hash = {};
  my $parent_gene_id_hash = {};
  my $genomic_reads = [];

  for(my $i=0; $i<scalar(@$sorted_input_genes); $i++) {
    my $gene = ${$sorted_input_genes}[$i];
    my $gene_id = $gene->dbID();
    my $gene_stable_id = $gene->stable_id();
    my $gene_version = $gene->version();
    my $gene_biotype = $gene->biotype();

    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      my $transcript_id = $transcript->dbID();
      my $transcript_stable_id = $transcript->stable_id();
      my $transcript_version = $transcript->version();
      my $biotype = $transcript->get_Biotype();
      my $biotype_group = $biotype->biotype_group();
      my $is_canonical = $transcript->is_canonical();
      my $source = $transcript->source();

      $parent_gene_id_hash->{$transcript_id}->{'gene_id'} = $gene_id;
      $parent_gene_id_hash->{$transcript_id}->{'gene_stable_id'} = $gene_stable_id;
      $parent_gene_id_hash->{$transcript_id}->{'gene_version'} = $gene_version;
      $parent_gene_id_hash->{$transcript_id}->{'gene_biotype'} = $gene_biotype;
      $parent_gene_id_hash->{$transcript_id}->{'transcript_stable_id'} = $transcript_stable_id;
      $parent_gene_id_hash->{$transcript_id}->{'transcript_version'} = $transcript_version;
      $parent_gene_id_hash->{$transcript_id}->{'biotype_group'} = $biotype_group;
      $parent_gene_id_hash->{$transcript_id}->{'is_canonical'} = $is_canonical;
      $parent_gene_id_hash->{$transcript_id}->{'source'} = $source;
    }

    my $slice = $gene->slice();
    my $stable_id = $gene->stable_id.".".$gene->version;


    # There's a big issue in terms of small genes, even with a fair amount of padding. To counter this have a minimum
    # target region size of about 50kb. If the gene is bigger then the padding itself should be enough as it's likely
    # there are many neutral sites in the gene already
    my $min_padding = 1000;
    my $min_region_size = 50000;
    my $region_diff = $min_region_size - $gene->length;
    if($region_diff > 0) {
      $region_diff = ceil($region_diff/2);
    } else {
      $region_diff = 0;
      $min_padding = 5000;
    }

    my $region_start = $gene->seq_region_start - $min_padding - $region_diff;
    if($region_start < $slice->seq_region_start()) {
       $region_start = $slice->seq_region_start();
    }

    my $region_end = $gene->seq_region_end + $min_padding + $region_diff;
    if($region_end > $slice->seq_region_end()) {
      $region_end = $slice->seq_region_end();
    }

    my $strand = $gene->strand;
    my $genomic_seq = ${ $sequence_adaptor->fetch_by_Slice_start_end_strand($slice,$region_start,$region_end,$strand) };
    say "FERGAL DEBUG ORIG REGION S/E, GENE S/E: ".$region_start.'/'.$region_end.", ". $gene->seq_region_start."/".$gene->seq_region_end;

    my $fasta_record = ">".$gene_id."\n".$genomic_seq;
    push(@$genomic_reads, $fasta_record);
    $gene_genomic_seqs_hash->{$gene_id} = [$region_start,$region_end,$genomic_seq];
  } #close foreach sorted_input_gene

  $self->param('gene_synteny_hash',$gene_synteny_hash);
  my $input_file = $self->write_input_file($genomic_reads);
  $self->param('input_file',$input_file);

  my $analysis = $self->create_analysis;
  $analysis->logic_name("minimap2remap");
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::LowDivergenceProjection->new(
       -analysis          => $analysis,
       -program           => $self->param('minimap2_path'),
       -paftools_path     => $self->param('paftools_path'),
       -genome_index      => $self->param_required('genome_index'),
       -input_file        => $input_file,
       -source_adaptor    => $source_gene_dba,
       -target_adaptor    => $target_gene_dba,
       -parent_genes      => $sorted_input_genes,
       -parent_gene_ids   => $parent_gene_id_hash,
       -gene_synteny_hash => $gene_synteny_hash,
       -gene_genomic_seqs_hash => $gene_genomic_seqs_hash,
       -no_projection     => $self->param('no_projection'),
       -coverage_cutoff   => $self->param('coverage_cutoff'),
       -perc_id_cutoff    => $self->param('perc_id_cutoff'),
       -extended_length_variation_cutoff => $self->param('extended_length_variation_cutoff'),
       -anchor_coverage_cutoff => $self->param('anchor_coverage_cutoff'),
       -anchor_perc_id_cutoff => $self->param('anchor_perc_id_cutoff'),
  );
  $self->runnable($runnable);
  if ($self->param('disconnect_jobs')) {
    $target_gene_dba->dbc->disconnect_when_inactive(1);
    $source_gene_dba->dbc->disconnect_when_inactive(1);
  }
}


sub run {
  my ($self) = @_;
  $self->dbc->disconnect_when_inactive(1) if ($self->param('disconnect_jobs'));
  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    if ($self->can('filter_results')) {
      $self->output($self->filter_results($runnable->output));
    }
    else {
      $self->output($runnable->output);
    }
  }
  $self->dbc->disconnect_when_inactive(0);
  return $self->output;
}


sub write_output {
  my ($self) = @_;

  my $output_gene_adaptor = $self->hrdb_get_con('target_gene_db')->get_GeneAdaptor;
  $output_gene_adaptor->dbc->disconnect_when_inactive(0);
  foreach my $output_gene (@{$self->output}) {
    say "Final gene: ".$output_gene->stable_id();
    empty_Gene($output_gene);
    $self->set_parent_attribs($output_gene);
    $output_gene_adaptor->store($output_gene);
  }

  return 1;
}


sub set_parent_attribs {
  my ($self,$gene) = @_;

  my ($current_gene_attrib) = @{$gene->get_all_Attributes('proj_parent_g')};
  unless($current_gene_attrib and $current_gene_attrib->value() =~ /^ENS/) {
    my $parent_versioned_stable_id = $gene->{'parent_gene_versioned_stable_id'};
    unless($parent_versioned_stable_id) {
      $self->throw("Couldn't find a parent_gene_versioned_stable_id key on the output gene, something has gone wrong");
    }

    my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_g',-VALUE => $parent_versioned_stable_id);
    $gene->add_Attributes($parent_attribute);
  }

  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@$transcripts) {
    my ($current_transcript_attrib) = @{$transcript->get_all_Attributes('proj_parent_t')};
    unless($current_transcript_attrib and $current_transcript_attrib->value() =~ /^ENS/) {
      my $parent_versioned_stable_id = $transcript->{'parent_transcript_versioned_stable_id'};
      unless($parent_versioned_stable_id) {
        $self->throw("Couldn't find a parent_transcript_versioned_stable_id key on the output gene, something has gone wrong");
      }
      my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_t',-VALUE => $parent_versioned_stable_id);
      $transcript->add_Attributes($parent_attribute);
    }
  }
}

sub fetch_input_genes_by_id {
  my ($self,$gene_ids,$source_gene_dba) = @_;

  my $input_genes = [];
  my $source_gene_adaptor = $source_gene_dba->get_GeneAdaptor();
  foreach my $gene_id (@$gene_ids) {
    my $gene = $source_gene_adaptor->fetch_by_dbID($gene_id);
    unless($gene) {
      $self->throw("Could not retrieve gene from source gene db with dbID: ".$gene->dbID());
    }

    if($gene->seq_region_name() eq 'MT') {
      next;
    }

    my $is_readthrough = 0;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      if($is_readthrough) {
         last;
      }

      my $attributes = $transcript->get_all_Attributes();
      foreach my $attribute (@{$attributes}) {
        if($attribute->value eq 'readthrough') {
          $is_readthrough = 1;
          last;
        }
      } # foreach my $attribute
    } # End foreach my $transcript

    unless($is_readthrough) {
      push(@$input_genes,$gene);
    }
  }

  return($input_genes);
}

sub write_input_file {
  my ($self,$genomic_reads) = @_;

  my $output_file = $self->create_filename();
  open(OUT,">".$output_file);
  foreach my $genomic_read (@$genomic_reads) {
    say OUT $genomic_read;
  }
  close OUT;

  return($output_file);
}


sub create_filename{
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}


1;
