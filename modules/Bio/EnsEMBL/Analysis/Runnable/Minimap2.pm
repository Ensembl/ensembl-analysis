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


package Bio::EnsEMBL::Analysis::Runnable::Minimap2;

use warnings;
use strict;
use feature 'say';

use File::Spec;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
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
  my ($genome_file, $input_file, $paftools_path, $database_adaptor) = rearrange([qw (GENOME_FILE INPUT_FILE PAFTOOLS_PATH DATABASE_ADAPTOR)],@args);
  $self->genome_file($genome_file);
  $self->input_file($input_file);
  $self->paftools_path($paftools_path);
  $self->database_adaptor($database_adaptor);
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

  my $output_file = $self->create_filename();
  $self->files_to_delete($output_file);

  my $genome_file   = $self->genome_file;
  my $input_file    = $self->input_file;
  my $paftools_path = $self->paftools_path;
  my $options       = $self->options;

  unless($paftools_path) {
    $self->throw("Paftools path was empty");
  }

  # run minimap2
  my $command = $self->program." -N 1 -ax splice -uf -C5 ".$genome_file." ".$input_file." | ".$paftools_path." splice2bed - > ".$output_file;
  $self->warning("Command:\n".$command."\n");
  if (system($command)) {
    $self->throw("Error running minimap2\nError code: $?\n");
  }

  $self->output($self->parse_results($output_file));

}


sub parse_results {
  my ($self,$output_file) = @_;

# 13  0   84793   ENST00000380152.7   1000    +   0   84793   0,128,255   27  194,106,249,109,50,41,115,50,112,1116,4932,96,70,428,182,188,171,355,156,145,122,199,164,139,245,147,2105,  0,948,3603,9602,10627,10768,11025,13969,15445,16798,20791,29084,31353,39387,40954,42268,47049,47705,54928,55482,61196,63843,64276,64533,79215,81424,82688,

  my $dba = $self->database_adaptor();
  my $slice_adaptor = $dba->get_SliceAdaptor();
  my $genes = [];

  unless(-e $output_file) {
    $self->throw("Output file does not exist. Path used:\n".$output_file);
  }

  open(IN,$output_file);
  while(<IN>) {
    my $line = $_;
    my @results = split("\t",$line);
    my $seq_region_name = $results[0];
    my $offset = $results[1];
    my $slice = $slice_adaptor->fetch_by_region('toplevel',$seq_region_name);
    my $hit_name = $results[3];
    my $strand = $results[5];
    if($strand eq '+') {
      $strand = 1;
    } elsif($strand eq '-') {
      $strand = -1;
    } else {
      $self->throw("Expected strand info to be + or -, found: ".$strand);
    }
    my $block_sizes = $results[10];
    my $block_starts = $results[11];

    my @block_sizes = split(",",$block_sizes);
    my @block_starts = split(",",$block_starts);

    my @exons = ();
    for(my $i=0; $i<scalar(@block_sizes); $i++) {
      my $block_start = $offset + $block_starts[$i] + 1; # We need to convert to 1-based
      my $block_end = $block_start + $block_sizes[$i] - 1;
      my $exon = $self->create_exon($slice,$block_start,$block_end,$strand);
      unless($exon) {
        $self->throw("Tried to create an exon and failed: ".$seq_region_name.", ".$block_start.", ".$block_end.", ".$strand);
      }
      push(@exons,$exon);
    }

    if($strand == -1) {
      @exons = reverse(@exons);
    }

    my $gene = $self->create_gene(\@exons,$slice);
    unless($self->filter_gene($gene)) {
#      say "Pushing gene: ".$gene->start."..".$gene->end;
      push(@$genes,$gene);
    }
  }
  close IN;

  return($genes);
}


sub create_exon {
  my ($self,$slice,$exon_start,$exon_end,$strand) = @_;

  my $exon = Bio::EnsEMBL::Exon->new(-start     => $exon_start,
                                     -end       => $exon_end,
                                     -strand    => $strand,
                                     -phase     => -1,
                                     -end_phase => -1,
                                     -analysis  => $self->analysis,
                                     -slice     => $slice);

#  say "Created exon: ".$slice->name." (".$exon_start."..".$exon_end.":".$strand.")";
  return($exon);
}


sub create_gene {
  my ($self,$exons,$slice) = @_;

  my $transcript = Bio::EnsEMBL::Transcript->new(-exons    => $exons,
                                                 -slice    => $slice,
                                                 -analysis => $self->analysis);

  compute_translation($transcript);
  my $gene = Bio::EnsEMBL::Gene->new(-slice    => $slice,
                                     -analysis => $self->analysis);

  $gene->add_Transcript($transcript);

  return($gene);
}


sub filter_gene {
  my ($self) = @_;
  return(0);
}


sub genome_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_genome_file} = $val;
  }

  return $self->{_genome_file};
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


sub database_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_database_adaptor} = $val;
  }

  return $self->{_database_adaptor};
}

1;
