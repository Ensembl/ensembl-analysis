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


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for minimap2
 Returntype : None
 Exceptions : Throws if 'genome_file' does not exist
              Throws if 'input_file' does not exist

=cut

sub fetch_input {
  my ($self) = @_;

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

  my $genome_file = $self->param_required('genome_file');
  unless(-e $genome_file) {
    $self->throw("Could not find the genome file. Path used:\n".$genome_file);
  }

  my $genome_index = $self->param_required('minimap2_genome_index');
  unless(-e $genome_index) {
    $self->throw("Could not find the genome index. Path used:\n".$genome_index);
  }

  my $input_file = $self->param('input_file');
  unless(-e $input_file) {
    $self->throw("Could not find the input file. Path used:\n".$input_file);
  }

  my $input_range = shift(@{$self->param('iid')});
  my $range_start = $$input_range[0];
  my $range_end = $$input_range[1];

  say "Creating tmp input file based on index range: ".$range_start."..".$range_end;
  my $ranged_input_file = $self->create_input_file($input_file,$range_start,$range_end);
  say "Finished creating input file";

#  say "FERGAL RANGED DEBUG: ".$ranged_input_file;
#  say "FERGAL RANGE DEBUG: ".$range_start."..".$range_end;
#  sleep(60);
#  $self->throw("DEBUG");
  my $program = $self->param('minimap2_path');
  my $paftools =  $self->param('paftools_path');

  unless($program) {
    $program = "minimap2";
  }

  unless($paftools) {
    $paftools = "paftools.js";
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis          => $self->create_analysis,
       -program           => $program,
       -paftools_path     => $paftools,
#       -options        => $self->param('minimap2_options'),
       -genome_index      => $genome_index,
       -input_file        => $ranged_input_file,
       -database_adaptor  => $target_dba,
       -delete_input_file => 1, # NB!! only set this when creating ranged files, not when using the original input file
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


#  my $counter = 0;
#  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($input_file);
#  while($parser->next() && $counter <= $end) {
#    if($counter < $start) {
#      $counter++;
#      next;
#    }
#    my $seq = $parser->getSequence();
#    if(length($seq < 200)) {
#      next;
#    }
#    my $header = $parser->getHeader();
#    say OUT ">".$header;
#    say OUT $seq;
#    $counter++;
#  }

  close OUT;
  return($output_file);
}


sub create_filename{
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}


1;
