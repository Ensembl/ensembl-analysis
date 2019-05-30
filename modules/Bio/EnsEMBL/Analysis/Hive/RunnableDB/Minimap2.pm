=head1 LICENSE

 Copyright [2019] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::Star

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Analysis::Runnable::Minimap2;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    input_id_type => 'fastq_range',
    canonical_only => 0,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for minimap2
 Returntype : None
 Exceptions : Throws if 'genome_file' does not exist
              Throws if 'input_file' does not exist

=cut

sub fetch_input {
  my ($self) = @_;

  my $program = $self->param('minimap2_path');
  my $paftools =  $self->param('paftools_path');

  unless($program) {
    $program = "minimap2";
  }

  unless($paftools) {
    $paftools = "paftools.js";
  }

  my $genome_file = $self->param_required('genome_file');
  unless(-e $genome_file) {
    $self->throw("Could not find the genome file. Path used:\n".$genome_file);
  }

  my $genome_index = $self->param_required('minimap2_genome_index');
  unless(-e $genome_index) {
    $self->throw("Could not find the genome index. Path used:\n".$genome_index);
  }

  my $master_input_file = $self->param('input_file');
  unless($self->param_required('input_id_type') eq "gene_id" || -e $master_input_file) {
    $self->throw("Could not find the input file. Path used:\n".$master_input_file);
  }

  my $runnable_input_file = "";
  my $delete_input_file = 0;
  if($self->param_required('input_id_type') eq "fastq_range") {
    my $input_range = shift(@{$self->param('iid')});
    my $range_start = $$input_range[0];
    my $range_end = $$input_range[1];

    say "Creating tmp input file based on index range: ".$range_start."..".$range_end;
    $runnable_input_file = $self->create_input_file($master_input_file,$range_start,$range_end);
    say "Finished creating input file";
    $delete_input_file = 1;
  } elsif($self->param('input_id_type') eq "gene_id") {
    my $source_gene_dba = $self->hrdb_get_dba($self->param_required('source_gene_db'));
    my $source_gene_dna_dba = $self->hrdb_get_dba($self->param_required('source_gene_dna_db'));
    $source_gene_dba->dnadb($source_gene_dna_dba);
    say "Creating tmp input file based on gene ids...";
    $runnable_input_file = $self->create_gene_input_file($self->param('iid'),$source_gene_dba);
    say "Finished creating input file";
    $delete_input_file = 1;
  } else {
    say "No input id type specified, will pass input file directly to runnable and will not delete";
    $runnable_input_file = $master_input_file;
  }

  my $target_dba = $self->hrdb_get_dba($self->param('target_db'));

  my $dna_dba;
  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $target_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($target_dba,'target_db');

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis          => $self->create_analysis,
       -program           => $program,
       -paftools_path     => $paftools,
       -genome_index      => $genome_index,
       -input_file        => $runnable_input_file,
       -database_adaptor  => $target_dba,
       -delete_input_file => $delete_input_file, # NB!! only set this when creating ranged files, not when using the original input file
    );

  $self->runnable($runnable);
}

=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $target_dba = $self->hrdb_get_con("target_db");
  my $gene_adaptor = $target_dba->get_GeneAdaptor();

  my $genes = $self->output();
  say "Total genes to output: ".scalar(@$genes);
  foreach my $gene (@$genes) {
    $gene_adaptor->store($gene);
  }

}


sub create_gene_input_file {
  my ($self,$gene_ids,$gene_dba) = @_;

  my $gene_adaptor = $gene_dba->get_GeneAdaptor();
  my $output_file = $self->create_filename();
  open(OUT,">".$output_file);
  foreach my $gene_id (@$gene_ids) {
    my $gene = $gene_adaptor->fetch_by_dbID($gene_id);
    unless($gene) {
      $self->throw("Found no gene in the gene source db for the following gene id: ".$gene_id);
    }


    my $transcripts = [];
    if($self->param('canonical_only')) {
      push(@$transcripts,$gene->canonical_transcript);
    } else {
      $transcripts = $gene->get_all_Transcripts;
    }

    foreach my $transcript (@$transcripts) {
      my $accession = "";
      if($transcript->stable_id) {
        $accession = $transcript->stable_id.".".$transcript->version;
       } else {
        $accession = $transcript->dbID;
      }

      unless($accession) {
        $self->throw("Could not find a transcript accession (stable id or dbID) for a transcript from the following gene_id: ".$gene_id);
      }

      my $seq = $transcript->seq->seq;
      say OUT ">".$accession;
      say OUT $seq;
    }
  }
  close OUT;
  return($output_file);
}


sub create_input_file {
  my ($self,$input_file,$start,$end) = @_;

  my $output_file = $self->create_filename();
  open(OUT,">".$output_file);
  my $seq_count = `wc -l $input_file | awk '{print \$1}'`;
  chomp($seq_count);
  $seq_count = $seq_count / 4;

  unless($seq_count > 0) {
    $self->throw("You have selected to generate a fasta range, but the fasta file doesn't have any headers. Path specified:\n".$input_file);
  }

  my $index = Bio::DB::HTS::Faidx->new($input_file);
  my @seq_ids = $index->get_all_sequence_ids();
  for(my $i=$start; $i<= $end; $i++) {
    my $seq_id = $seq_ids[$i];
    my $length = $index->length($seq_id);
    my $location  = $seq_id.":1-".$length;
    my $seq = $index->get_sequence_no_length($location);

    say OUT ">".$seq_id;
    say OUT $seq;
  }
  close OUT;
  return($output_file);
}


sub create_filename{
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}


1;
