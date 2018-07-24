# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::RecalculateRoughModels;

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
    refine_gene_padding   => 50,
    min_intergenic_length => 500,
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

  my $dna_db = $self->get_database_by_name('dna_db');
  my $refine_db = $self->get_database_by_name('refine_db');
  my $rough_db = $self->get_database_by_name('rough_db');

  my $padding = $self->param_required('refine_gene_padding');
  my $refine_logic_name = $self->param_required('logic_name');

  my $min_intergenic_length = $self->param_required('min_intergenic_length');

  my $slice_name = $self->param_required('iid');
  my $slice = $dna_db->get_SliceAdaptor->fetch_by_name($slice_name);
  $self->query($slice);

  my $refine_genes = $refine_db->get_GeneAdaptor->fetch_all_by_Slice_constraint($slice,undef,$refine_logic_name);
  say "Found ".scalar(@$refine_genes)." genes for ".$refine_logic_name." on slice ".$slice->name;

  foreach my $gene (@$refine_genes) {
    say "RG: ".$gene->seq_region_start.":".$gene->seq_region_end.":".$gene->strand;
  }
  my $input_gene_coords = $self->calculate_input_gene_coords($refine_genes);
  $self->input_gene_coords($input_gene_coords);

}


sub run {
  my ($self) = @_;

  my $input_gene_coords = $self->input_gene_coords();
  say "Processing input genes to find intergenic regions...";
  my $output_ids = $self->find_intergenic_coords($input_gene_coords);
  say "...finished finding intergenic regions";

  foreach my $output_id (@$output_ids) {
    say "INTERGENIC SLICE: ".$output_id;
  }

  $self->throw("FERGAL TEST");
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


sub calculate_input_gene_coords {
  my ($self,$genes) = @_;

  my @sorted_start_end_coords = ();
  my @unsorted_start_end_coords = ();
  my $slice = $self->query;
  my $skip_duplicates = {};
  my $filter = {};

  foreach my $gene (@{$genes}) {
    my $start = $gene->seq_region_start;
    my $end = $gene->seq_region_end;
    if($skip_duplicates->{$start.":".$end}) {
      next;
    } else {
      push(@unsorted_start_end_coords,{'s' => $start, 'e' => $end});
      $skip_duplicates->{$start.":".$end} = 1;
    }
  }

  @sorted_start_end_coords = sort {
                                    $a->{'s'} <=> $b->{'s'} ||
                                    $a->{'e'} <=> $b->{'e'}
                                  } @unsorted_start_end_coords;

  return \@sorted_start_end_coords;
}


sub input_gene_coords {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_input_gene_coords',$val);
  }

  return($self->param('_input_gene_coords'));
}


sub find_intergenic_coords {
  my ($self,$input_gene_coords) = @_;

  my $slice = $self->query;
  my @intergenic_coords = ();
  my @slice_eles = split(':',$slice->name());

  for(my $i=0; $i<scalar(@{$input_gene_coords}); $i++) {
    my $right_start = ${$input_gene_coords}[$i]->{'s'};
    my $right_end = ${$input_gene_coords}[$i]->{'e'};

    if($i==0) {
      $slice_eles[3] = $slice->start;
      $slice_eles[4] = $right_start-1;
      push(@intergenic_coords,join(':',@slice_eles));
      next;
    }

    my $left_start = ${$input_gene_coords}[$i-1]->{'s'};
    my $left_end = ${$input_gene_coords}[$i-1]->{'e'};

    $slice_eles[3] = $left_end+1;
    $slice_eles[4] = $right_start-1;
    push(@intergenic_coords,join(':',@slice_eles));
  }

  $slice_eles[3] = ${$input_gene_coords}[$#{$input_gene_coords}]->{'e'} + 1;
  $slice_eles[4] = $slice->seq_region_end;
  push(@intergenic_coords,join(':',@slice_eles));

  return(\@intergenic_coords);
}


1;
