=head1 LICENSE

 Copyright [2019-2020] EMBL-European Bioinformatics Institute

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

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Star

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::Star->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence. Star is a splice aware
aligner. It creates output files with the reads overlapping splice sites and the reads
aligning on the exons. Some reads are aligned multiple times in the genome.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::Minimap2;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name);
use File::Spec;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable');


=head2 new

 Arg [DECOMPRESS]           : String as a command like 'gzip -c -'
 Arg [EXPECTED_ATTRIBUTES]  : String specify the attribute expected for the output, see STAR manual
 Description                : Creates a  object to align reads to a genome using STAR
 Returntype                 : 
 Exceptions                 : Throws if WORKDIR does not exist
                              Throws if the genome has not been indexed

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($genome_index, $input_file, $paftools_path, $source_adaptor, $target_adaptor, $delete_input_file, $parent_gene_ids) = rearrange([qw (GENOME_INDEX INPUT_FILE PAFTOOLS_PATH SOURCE_ADAPTOR TARGET_ADAPTOR DELETE_INPUT_FILE PARENT_GENE_IDS)],@args);
  $self->genome_index($genome_index);
  $self->input_file($input_file);
  $self->paftools_path($paftools_path);
  $self->source_adaptor($source_adaptor);
  $self->target_adaptor($target_adaptor);
  $self->delete_input_file($delete_input_file);
  $self->parent_gene_ids($parent_gene_ids);
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run Star to align reads to an indexed genome. The resulting output file will be stored in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $leftover_genes = [];
  my $paf_file = $self->create_filename(undef,'paf');
  $self->files_to_delete($paf_file);

  my $genome_index  = $self->genome_index;
  my $input_file    = $self->input_file;
  my $options = $self->options;

  # run minimap2
  my $minimap2_command = $self->program." --cs --secondary=no -x map-pb ".$genome_index." ".$input_file." > ".$paf_file;
  $self->warning("Command:\n".$minimap2_command."\n");
  if(system($minimap2_command)) {
    $self->throw("Error running minimap2\nError code: $?\n");
  }

  my $paf_results = [];
  open(IN,$paf_file);
  while(<IN>) {
    chomp($_);
    push(@$paf_results,$_);
  }
  close IN;

  foreach my $paf_result (@$paf_results) {
    my @result_cols = split("\t",$paf_result);
    $self->process_results(\@result_cols);
  }
} # End run



sub process_results {
  my ($self, $paf_result) = @_;

  my $source_gene_id = ${$paf_result}[0];
  my $source_genomic_length = ${$paf_result}[1];
  my $source_genomic_start = ${$paf_result}[2];
  my $source_genomic_end = ${$paf_result}[3];
  my $same_strand = ${$paf_result}[4];
  my $target_genomic_name = ${$paf_result}[5];
  my $target_genomic_length = ${$paf_result}[6];
  my $target_genomic_start = ${$paf_result}[7];
  my $target_genomic_end = ${$paf_result}[8];
  my $matching_bases = ${$paf_result}[9];
  my $total_bases = ${$paf_result}[10];
  my $mapping_quality = ${$paf_result}[11];

  my $target_adaptor = $self->target_adaptor();
  my $target_slice_adaptor = $target_adaptor->get_SliceAdaptor();
  my $target_parent_slice = $target_slice_adaptor->fetch_by_region('toplevel',$target_genomic_name);
  my $target_strand = 1;
  my $target_sequence_adaptor = $target_adaptor->get_SequenceAdaptor;
  my $target_genomic_seq = ${ $target_sequence_adaptor->fetch_by_Slice_start_end_strand($target_parent_slice, $target_genomic_start, $target_genomic_end, $target_strand) };
  my $target_genomic_fasta = ">".$target_genomic_name."\n".$target_genomic_seq;
  my $target_genome_file = $self->write_input_file([$target_genomic_fasta]);
  my $target_genome_index = $target_genome_file.".mmi";
  my $target_index_command = $self->program()." -d ".$target_genome_index." ".$target_genome_file;
  my $index_result = system($target_index_command);
  if($index_result) {
    $self->throw('The minimap2 index command returned a non-zero exit code. Commandline used:\n'.$target_index_command);
  }

  # Covert the transcripts into fasta records
  my $source_adaptor = $self->source_adaptor();
  my $source_gene_adaptor = $source_adaptor->get_GeneAdaptor();
  my $source_gene = $source_gene_adaptor->fetch_by_dbID($source_gene_id);
  my $source_transcripts = $source_gene->get_all_Transcripts();

  my $source_transcript_fasta_seqs = [];
  foreach my $source_transcript (@$source_transcripts) {
    my $source_transcript_sequence = $source_transcript->seq->seq();
    my $fasta_record = ">".$source_transcript->dbID()."\n".$source_transcript_sequence;
    push(@$source_transcript_fasta_seqs,$fasta_record);
  }

  my $analysis = $self->analysis();
  my $program = $self->program();
  my $paftools = $self->paftools_path();
  my $source_input_file = $self->write_input_file($source_transcript_fasta_seqs);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis                 => $analysis,
       -program                  => $program,
       -paftools_path            => $paftools,
       -genome_index             => $target_genome_index,
       -input_file               => $source_input_file,
       -database_adaptor         => $target_adaptor,
       -skip_introns_check       => 1,
       -add_offset               => $target_genomic_start - 1,
       -skip_compute_translation => 1,
  );

  $runnable->run();
  my $output_genes = $runnable->output();
  my $parent_gene_ids = $self->parent_gene_ids();
  my $final_gene_hash = {};
  foreach my $gene (@$output_genes) {
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      unless($parent_gene_ids->{$transcript->stable_id()}) {
        $self->throw("The following mapped transcript stable id was not found in the initial list of dbIDs: ".$transcript->stable_id());
      }

      my $parent_gene_id = $parent_gene_ids->{$transcript->stable_id()}->{'gene_id'};
      my $biotype_group = $parent_gene_ids->{$transcript->stable_id()}->{'biotype_group'};

      unless($final_gene_hash->{$parent_gene_id}) {
        $final_gene_hash->{$parent_gene_id} = [];
      }

      if($biotype_group eq 'coding') {
        compute_translation($transcript);
      }

      push(@{$final_gene_hash->{$parent_gene_id}},$transcript);
    }
  } # End foreach my $gene

  my $final_genes = [];
  foreach my $gene_id (keys(%{$final_gene_hash})) {
    my $transcripts = $final_gene_hash->{$gene_id};
    my $gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);
    foreach my $transcript (@$transcripts) {
      $gene->add_Transcript($transcript);
    }
    push(@$final_genes, $gene);
  }

  $self->output($final_genes);
}


sub write_input_file {
  my ($self,$fasta_records) = @_;

  my $output_file = $self->create_filename();
  open(OUT,">".$output_file);
  foreach my $fasta_record (@$fasta_records) {
    say OUT $fasta_record;
  }
  close OUT;

  return($output_file);
}


sub create_filename{
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}



sub genome_index {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_genome_index} = $val;
  }

  return $self->{_genome_index};
}


sub input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_input_file} = $val;
  }

  return $self->{_input_file};
}


sub paftools_path {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_paftools_path} = $val;
  }

  return $self->{_paftools_path};
}


sub source_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_source_adaptor} = $val;
  }

  return $self->{_source_adaptor};
}


sub target_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_target_adaptor} = $val;
  }

  return $self->{_target_adaptor};
}


sub delete_input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_delete_input_file} = $val;
  }

  return $self->{_delete_input_file};
}


sub parent_gene_ids {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_parent_gene_ids} = $val;
  }

  return $self->{_parent_gene_ids};
}

1;
