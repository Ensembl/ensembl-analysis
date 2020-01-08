# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Bam2Genes

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $refine_genes = Bio::EnsEMBL::Analysis::RunnableDB::Bam2Genes->new (
          -db      => $db,
          -input_id   => $input_id
          -analysis   => $analysis );
  $refine_genes->fetch_input();
  $refine_genes->run();
  $refine_genes->write_output(); #writes to DB

=head1 DESCRIPTION

The module creates "proto-transcripts" based on the alignments of short reads.
It will first create blocks from overlapping reads which represent possible exons and
it will link these blocks by using pairing information if it is available or by
using a predefined "proto-transcript" length when it uses single reads.
The "proto-transcripts" will be stored in an Ensembl database.
If your genome file is using UCSC naming, ie. chr1, set 'wide_use_ucsc_naming' to 1

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::Bam2GenesDepth;

use warnings ;
use strict;
use feature 'say';

use Bio::DB::HTS;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Runnable::Bam2Genes;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(create_Exon);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    stranded_read => 0,
    _repeat_libs => ['dust', 'trf'],
    _repeats_allowed => 0.95,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: It will fetch the sequence based on the input_id, retrieve all alignments
               overlapping the region and create the "exon" blocks. Because the alignment
               is made on the whole genome, we need to provide the runnable with a full length
               slice. If your genome use UCSC style names (chr1,...), set 'wide_use_ucsc_naming' to 1
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
    my ($self) = @_;

    my $reference_db = $self->get_database_by_name('dna_db');
    my $slice_adaptor = $reference_db->get_SliceAdaptor;
    $self->hrdb_set_con($self->get_database_by_name('output_db'), 'output_db');

    my $slice = $self->fetch_sequence($self->input_id, $reference_db);
    if ($self->param('disconnect_jobs')) {
      $reference_db->dbc->disconnect_if_idle;
      $self->dbc->disconnect_if_idle;
    }
    my $sam = Bio::DB::HTS->new(
            -bam => $self->param('alignment_bam_file'),
            -expand_flags => 1,
            );
    $self->throw('Bam file ' . $self->param('alignment_bam_file') . "  not found \n") unless $sam;
    $self->create_analysis;

    my $seq_region_name = $slice->seq_region_name;
    my $region = $seq_region_name.':'.$slice->start.'-'.$slice->end;
    say "Getting depth: ";
    my $depth_array = $self->get_base_depth($self->param('alignment_bam_file'),$region);
    say "...got depth";

   $self->depth_array($depth_array);
   $self->slice($slice);
}


sub run {
  my ($self) = @_;
  my $rough_exons = $self->generate_rough_exons();
  $self->rough_exons($rough_exons);
}

=head2 write_output

 Arg [1]    : None
 Description: Write the proto transcripts in the database specified as 'output_db'.
 Returntype : None
 Exceptions : Throws if it cannot write all the genes in the database for any reason

=cut

sub write_output{
    my ($self) = @_;

    my $outdb = $self->hrdb_get_con('output_db');
    my $gene_adaptor = $outdb->get_GeneAdaptor;

    my $rough_exons = $self->rough_exons;

    $gene_adaptor->dbc->disconnect_when_inactive(0);
    foreach my $exon ( @{$rough_exons} ) {
      my $exon_gene = $self->make_gene($exon);
      $gene_adaptor->store($exon_gene);
    }
}


sub generate_rough_exons {
  my ($self) = @_;

  my $rough_exons = [];
  my $depth_array = $self->depth_array();
  my $slice = $self->slice();
  my @clusters = ();
  my $current_start = 0;
  my $current_end = 0;
  my $min_rough_exon_size = 20;
  my $min_depth = 50;
  my $in_potential_exon = 0;
  foreach(my $i=$slice->start; $i <= $slice->end; $i++) {
    my $res = $$depth_array[$i] ? $$depth_array[$i] : 0;

    if($res >= $min_depth && !$in_potential_exon) {
      $in_potential_exon = 1;
      $current_start = $i;
      $current_end = $i;
    } elsif($res >= $min_depth) {
      $current_end = $i;
    } elsif($res < $min_depth && ($current_end - $current_start + 1) >= $min_rough_exon_size) {
      say "FERGAL FOUND EXON: ".$current_start.":".$current_end;
      push(@clusters,[$current_start,$current_end]);
      $current_start = 0;
      $current_end = 0;
      $in_potential_exon = 0;
    } elsif($res < $min_depth && ($current_start && $current_end) && ($current_end - $current_start + 1) < $min_rough_exon_size) {
      say "FERGAL SKIPPING SMALL SIGNAL: ".$current_start.":".$current_end;
      $current_start = 0;
      $current_end = 0;
      $in_potential_exon = 0;
    } elsif($res < $min_depth && ($current_end - $current_start + 1) < $min_rough_exon_size) {
      $current_start = 0;
      $current_end = 0;
      $in_potential_exon = 0;
    } else {
      $self->throw("Shouldn't be logically able to get here!!!!");
    }
  }

  if(($current_start && $current_end) && ($current_end - $current_start + 1) >= $min_rough_exon_size) {
    say "FERGAL FOUND LAST EXON: ".$current_start.":".$current_end;
    push(@clusters,[$current_start,$current_end]);
  }

  my $full_slice = $self->slice->seq_region_Slice;
  foreach my $cluster (@clusters) {
    my ($start,$end) = @{$cluster};

    my $rough_exon = create_Exon(
                                  $start,
                                  $end,
                                  -1,
                                  -1,
                                  -1,
                                  $self->analysis,
                                  undef,
                                  undef,
                                  $full_slice,
                                );


    push(@$rough_exons,$rough_exon);
  }

  return($rough_exons);
}


sub get_base_depth {
  my ($self,$bam_file,$region) = @_;

  my $depth_array = [];
#  open(SAMTOOLS_DEPTH,"samtools depth -Q 30 -r ".$region." ".$bam_file." |");
  open(SAMTOOLS_DEPTH,"samtools depth -r ".$region." ".$bam_file." |");
  while(my $depth_line = <SAMTOOLS_DEPTH>) {
    chomp($depth_line);
    my ($region,$base,$cov) = split('\t',$depth_line);
    ${$depth_array}[$base] = $cov;
  }
  close SAMTOOLS_DEPTH;

  return($depth_array);
}


sub depth_array {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_depth_array',$val);
  }

  return($self->param('_depth_array'));
}


sub rough_exons {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_rough_exons',$val);
  }

  return($self->param('_rough_exons'));
}


sub slice {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_slice',$val);
  }

  return($self->param('_slice'));
}


sub make_gene {
  my ($self,$exon) = @_;

  my $transcript =  new Bio::EnsEMBL::Transcript(-EXONS => [$exon]);
  $transcript->analysis($self->analysis);

  my $gene = new Bio::EnsEMBL::Gene();
  $gene->add_Transcript($transcript);
  $gene->analysis($self->analysis);

  return($gene);
}

1;
