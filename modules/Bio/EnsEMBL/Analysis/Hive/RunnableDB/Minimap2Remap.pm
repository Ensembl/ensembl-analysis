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

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (align_nucleotide_seqs);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;
  $self->create_analysis;

  $self->param('region_padding',500);

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');

  # Define the source and target gene dbs
  my $source_gene_dba = $self->hrdb_get_dba($self->param('source_gene_db'));
  my $target_gene_dba = $self->hrdb_get_dba($self->param('target_gene_db'));
  $source_gene_dba->dnadb($source_dna_dba);
  $target_gene_dba->dnadb($target_dna_dba);
  $self->hrdb_set_con($source_gene_dba,'source_gene_db');
  $self->hrdb_set_con($target_gene_dba,'target_gene_db');


  my $slice_array = $self->param('iid');
  my $slice_adaptor = $source_gene_dba->get_SliceAdaptor();
  my $sequence_adaptor = $source_gene_dba->get_SequenceAdaptor();
  my $input_genes = [];
  foreach my $slice_name (@$slice_array) {
    my $slice = $slice_adaptor->fetch_by_name($slice_name);
    my $genes = $slice->get_all_Genes();
    my $filtered_genes = [];
    # Put the filtering in a subroutine at some point if it gets more complex
    # At the moment the filtering is just on readthroughs
    say "Got ".scalar(@$genes)." unfiltered genes";
    foreach my $gene (@$genes) {
      # Skip X genes on Y
      if($slice->seq_region_name eq 'Y' and $gene->seq_region_name ne 'Y') {
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
        push(@$filtered_genes,$gene);
      }
    } # End foreach my $gene
    push(@$input_genes,@$filtered_genes);
  }

############################################
# TEST

#  my $all_slices = $slice_adaptor->fetch_all('toplevel');
#  foreach my $slice (@$all_slices) {
#    my $genes = $slice->get_all_Genes();
#    push(@$input_genes,@$genes);
#  }

############################################

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

  my $gene_synteny_hash = {};

  my $parent_gene_id_hash = {};
  my $genomic_reads = [];
#  foreach my $gene (@$sorted_input_genes) {
  for(my $i=0; $i<scalar(@$sorted_input_genes); $i++) {
    my $gene = ${$sorted_input_genes}[$i];
    my $gene_id = $gene->dbID();
    my $gene_stable_id = $gene->stable_id();
    my $gene_version = $gene->version();
    my $gene_biotype = $gene->biotype();

    # Store the ids of up to four genes, two upstream and two downstream around the current gene
    if($i-1 >= 0 && ${$sorted_input_genes}[$i-1]->slice->name() eq $gene->slice->name()) {
      $gene_synteny_hash->{$gene_id}->{${$sorted_input_genes}[$i-1]->dbID()} = 1;
    }

    if($i-2 >= 0 && ${$sorted_input_genes}[$i-2]->slice->name() eq $gene->slice->name()) {
      $gene_synteny_hash->{$gene_id}->{${$sorted_input_genes}[$i-2]->dbID()} = 1;
    }

    if($i+1 < scalar(@$sorted_input_genes) && ${$sorted_input_genes}[$i+1]->slice->name() eq $gene->slice->name()) {
      $gene_synteny_hash->{$gene_id}->{${$sorted_input_genes}[$i+1]->dbID()} = 1;
    }

    if($i+2 < scalar(@$sorted_input_genes) && ${$sorted_input_genes}[$i+2]->slice->name() eq $gene->slice->name()) {
      $gene_synteny_hash->{$gene_id}->{${$sorted_input_genes}[$i+2]->dbID()} = 1;
    }


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


    my $region_start = $gene->seq_region_start - $self->param('region_padding');
    if($region_start <= 0) {
      $region_start = 1;
    }

    my $region_end = $gene->seq_region_end + $self->param('region_padding');
    if($region_end > $slice->length()) {
      $region_end = $slice->length();
    }

    my $strand = $gene->strand;
    my $genomic_seq = ${ $sequence_adaptor->fetch_by_Slice_start_end_strand($slice, $region_start, $region_end, $strand) };
    my $fasta_record = ">".$gene_id."\n".$genomic_seq;
    push(@$genomic_reads, $fasta_record);
  } #close foreach sorted_input_gene

  $self->param('gene_synteny_hash',$gene_synteny_hash);
  my $input_file = $self->write_input_file($genomic_reads);
  $self->param('input_file',$input_file);

  my $analysis = $self->create_analysis;
  $analysis->logic_name("minimap2remap");
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap->new(
       -analysis          => $analysis,
       -program           => $self->param('minimap2_path'),
       -paftools_path     => $self->param('paftools_path'),
       -genome_index      => $self->param_required('genome_index'),
       -input_file        => $input_file,
       -source_adaptor    => $source_gene_dba,
       -target_adaptor    => $target_gene_dba,
       -parent_gene_ids   => $parent_gene_id_hash,
       -gene_synteny_hash => $gene_synteny_hash,
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

  my $final_genes = $self->output();
  my $sorted_final_genes = [sort { $a->slice->name() cmp $b->slice->name() or
                                   $a->start() <=> $b->start() or
                                   $a->end() <=> $b->end() }  @{$final_genes}];

  my $gene_synteny_hash = $self->param('gene_synteny_hash');
  for(my $i=0; $i<scalar(@$sorted_final_genes); $i++) {
    my $gene = ${$sorted_final_genes}[$i];
    my $gene_id = $gene->{'parent_gene_id'};
    my $synteny_score = 0;
#    say "FERGAL SYN GID: ".$gene_id;
#    say "FERGAL SYN DMP: ".Dumper($gene_synteny_hash->{$gene_id});
#    say "FERGAL SYN SLC: ".$gene->slice->name();


    # Store the ids of up to four genes, two upstream and two downstream around the current gene
    if($i-1 >= 0 && ${$sorted_final_genes}[$i-1]->slice->name() eq $gene->slice->name() &&
       $gene_synteny_hash->{$gene_id}->{${$sorted_final_genes}[$i-1]->{'parent_gene_id'}}) {
#       my $flanking_gene_id = ${$sorted_final_genes}[$i-1]->{'parent_gene_id'};
#       say "FERGAL SYN FLK: ".$flanking_gene_id;
       $synteny_score++;
     }

    if($i-2 >= 0 && ${$sorted_final_genes}[$i-2]->slice->name() eq $gene->slice->name() &&
       $gene_synteny_hash->{$gene_id}->{${$sorted_final_genes}[$i-2]->{'parent_gene_id'}}) {
       $synteny_score++;
     }

    if($i+1 < scalar(@$sorted_final_genes) && ${$sorted_final_genes}[$i+1]->slice->name() eq $gene->slice->name() &&
       $gene_synteny_hash->{$gene_id}->{${$sorted_final_genes}[$i+1]->{'parent_gene_id'}}) {
       $synteny_score++;
     }

    if($i+2 < scalar(@$sorted_final_genes) && ${$sorted_final_genes}[$i+2]->slice->name() eq $gene->slice->name() &&
       $gene_synteny_hash->{$gene_id}->{${$sorted_final_genes}[$i+2]->{'parent_gene_id'}}) {
       $synteny_score++;
    }

    say "FERGAL SYNTENY SCORE: ".$synteny_score;
    if($synteny_score < 2) {
      $gene->description("weak synteny");
    }

  }

  return $self->output;
}



sub write_output {
  my ($self) = @_;
  my $output_dba = $self->hrdb_get_con('target_gene_db');
  my $output_gene_adaptor = $output_dba->get_GeneAdaptor;
  my $output_genes = $self->output();
  foreach my $output_gene (@$output_genes) {
    $output_gene_adaptor->store($output_gene);
  }

  return 1;
}


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
