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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2Remap

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2Remap;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use File::Spec::Functions qw(tmpdir catfile);
use POSIX;
use List::Util qw(min max);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs);
use Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;

  if ($self->param('use_genome_flatfile')) {
    $self->setup_fasta_db;
  }
  $self->create_analysis;
  my $genome_index = $self->param_required('genome_index');

#  $self->param('region_padding',10000);

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');
  if ($self->param('use_genome_flatfile')) {
    if ($self->param_is_defined('source_dna_fasta')) {
      $source_dna_dba->get_SequenceAdaptor->fasta($self->param_required('source_dna_fasta'));
    }
    my $target_dna_fasta = $genome_index;
    $target_dna_fasta =~ s/(\.fa(sta)?).*$/$1/;
    $target_dna_dba->get_SequenceAdaptor->fasta($target_dna_fasta);
  }

  # Define the source and target gene dbs
  my $source_gene_dba = $self->hrdb_get_dba($self->param('source_gene_db'));
  my $target_gene_dba = $self->hrdb_get_dba($self->param('target_gene_db'));
  $source_gene_dba->dnadb($source_dna_dba);
  $target_gene_dba->dnadb($target_dna_dba);
  $self->hrdb_set_con($source_gene_dba,'source_gene_db');
  $self->hrdb_set_con($target_gene_dba,'target_gene_db');

# [19649] -> this is a gene/batch mapped to the -1 strand in the target, this e.g. has a 0 dist and no issues, so should be interesting
#  my $input_genes = $self->fetch_input_genes_by_id([59266,59269,59270,59293,59271,59274,59295,59276,59278],$source_gene_dba);
#my $input_genes = $self->fetch_input_genes_by_id([59266],$source_gene_dba);

  my $input_id_file = $self->param_required('input_id_file');
  unless(-e $input_id_file) {
    $self->throw("Did not find an input id file containing the list of source genes to be projected/mapped");
  }

  # TEST
#  my $test_slice = $target_gene_dba->get_SliceAdaptor->fetch_by_region('toplevel','17');
#  my $initial_projected_genes = $target_gene_dba->get_GeneAdaptor->fetch_all_by_Slice($test_slice);
#  my $initial_projected_genes = $target_gene_dba->get_GeneAdaptor->fetch_all();

#  my ($missing_genes,$problematic_transcripts) = $self->analyses_initial_projections($initial_projected_genes,$input_id_file);

  my $input_genes = $self->fetch_source_genes($input_id_file,$source_gene_dba);
#  my $test_slice = $source_gene_dba->get_SliceAdaptor->fetch_by_region('toplevel','7');
#  my $input_genes = $source_gene_dba->get_GeneAdaptor->fetch_all_by_Slice($test_slice);

  my $sequence_adaptor = $source_dna_dba->get_SequenceAdaptor();
  my $target_gene_adaptor = $target_gene_dba->get_GeneAdaptor();


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


  my $parent_gene_id_hash = {};

#  $self->set_parent_info($sorted_input_genes,$parent_gene_id_hash,$sequence_adaptor);


  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap->new(
       -analysis          => $self->analysis,
       -program           => $self->param('minimap2_path'),
       -paftools_path     => $self->param('paftools_path'),
       -genome_index      => $genome_index,
       -source_adaptor    => $source_gene_dba,
       -target_adaptor    => $target_gene_dba,
       -parent_genes      => $sorted_input_genes,
       -parent_gene_ids   => $parent_gene_id_hash,
       -no_projection     => $self->param('no_projection'),
  );
  $self->runnable($runnable);
}


sub run {
  my ($self) = @_;
  $self->dbc->disconnect_if_idle() if ($self->param('disconnect_jobs'));
  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    if ($self->can('filter_results')) {
      $self->output($self->filter_results($runnable->output));
    }
    else {
      $self->output($runnable->output);
    }
  }

  return $self->output;
}


sub write_output {
  my ($self) = @_;
  my $output_dba = $self->hrdb_get_con('target_gene_db');

  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  my $target_gene_dba = $self->hrdb_get_dba($self->param('target_gene_db'));
  $target_gene_dba->dnadb($target_dna_dba);

  my $output_gene_adaptor = $target_gene_dba->get_GeneAdaptor;
  my $output_genes = $self->output();
  foreach my $output_gene (@$output_genes) {
#    say "Final gene: ".$output_gene->stable_id()." ".$output_gene->seq_region_start.":".$output_gene->seq_region_end.":".$output_gene->seq_region_strand.":".$output_gene->seq_region_name;
    if($output_gene->{'to_remove'} and !($output_gene->{'to_write'})) {
      $output_gene_adaptor->remove($output_gene);
    } elsif($output_gene->{'to_write'} and !($output_gene->{'to_remove'})) {
      empty_Gene($output_gene);
      $output_gene_adaptor->store($output_gene);
    }
  }

  return 1;
}


sub fetch_source_genes {
  my ($self,$input_id_file,$source_gene_dba) = @_;

  my $source_genes = [];

  my @id_list;
  open(IN,$input_id_file) or $self->throw("Could not open $input_id_file");
  while(<IN>) {
    my $line = $_;
    my @eles = split("\t",$line);
    push(@id_list, $eles[0]);
  }
  close IN or $self->throw("Could not close $input_id_file");
  $source_genes = $source_gene_dba->get_GeneAdaptor->fetch_all_by_dbID_list(\@id_list);
  my %unique_ids = map {$_ => 1} @id_list;
  if (@$source_genes != scalar(keys %unique_ids)) {
    $self->throw("Fetched ".scalar(@$source_genes).' genes but expected '.scalar(@id_list));
  }

  return($source_genes);
}

sub set_parent_info {
  my ($self,$sorted_input_genes,$parent_gene_id_hash,$sequence_adaptor) = @_;

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
  }
}

sub calculate_target_midpoint {
  my ($self,$source_hit_start,$source_hit_end,$source_hit_midpoint,$target_genomic_start,$target_genomic_end,$target_strand,$target_genomic_length) = @_;

  # In this case the midpoint is within the hit boundaries, so no adjustment is needed. Note that the midpoint is only quite approximate as there
  # could be changes in length between the sounce and target
  my $target_midpoint;
  if($source_hit_midpoint >= $source_hit_start and $source_hit_midpoint <= $source_hit_end) {
    if($target_strand eq '+') {
      my $source_hit_offset = $source_hit_midpoint - $source_hit_start;
      # Could put in a target offset here in terms of the difference in the length of the source and target sets
      $target_midpoint = $target_genomic_start + $source_hit_offset;
      return($target_midpoint);
    } else {
      my $source_hit_offset = $source_hit_end - $source_hit_midpoint;
      $target_midpoint = $target_genomic_start + $source_hit_offset;
    }
  } else{
    # This means the source midpoint is outside of the hit, so we want to guess where it is, though the guess is going to be wrong since the implication is
    # that the reason the hit is incomplete is because there's some sort of gap
    if($target_strand eq '+') {
      if($source_hit_midpoint < $source_hit_start) {
        my $source_hit_offset = $source_hit_start - $source_hit_midpoint;
        $target_midpoint = $target_genomic_start - $source_hit_offset;
        if($target_midpoint < 1) {
          $target_midpoint = 1;
        }
        return($target_midpoint);
      } else {
        my $source_hit_offset = $source_hit_midpoint - $source_hit_end;
        $target_midpoint = $target_genomic_end + $source_hit_offset;
        if($target_midpoint > $target_genomic_length) {
          $target_midpoint = $target_genomic_length;
        }
        return($target_midpoint);
      }
    } else {
      if($source_hit_midpoint < $source_hit_start) {
        my $source_hit_offset = $source_hit_start - $source_hit_midpoint;
        $target_midpoint = $target_genomic_end + $source_hit_offset;
        if($target_midpoint > $target_genomic_length) {
          $target_midpoint = $target_genomic_length;
        }
        return($target_midpoint);
      } else {
        my $source_hit_offset = $source_hit_midpoint - $source_hit_end;
        $target_midpoint = $target_genomic_start - $source_hit_offset;
        if($target_midpoint < 1) {
          $target_midpoint = 1;
        }
        return($target_midpoint);
      }
    }
  } # End outer else
}


1;
