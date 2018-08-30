=head1 LICENSE

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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME


Bio::EnsEMBL::Analysis::RunnableDB::HiveRemoveDuplicateGenes - 

=head1 SYNOPSIS

Get lincRNA candidate genes and break them down into single-transcript genes.

=head1 DESCRIPTION

This module collects the genes of a region, checks for overlaps, collapses overlapping models (with genebuilder) 
and writes back the genebuilder "unique" results. This module is usefull for lincRNA pipeline, it collapses
all RNAseq models of the different tissues. This is memory efficient and avoids doing the same calculations 
many times. 

=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveDuplicateGenes;


use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::LincRNA qw(get_genes_of_biotypes_by_db_hash_ref) ;  
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(print_Gene attach_Slice_to_Gene empty_Gene);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input{
  my ($self) = @_;

  # set up config
  $self->create_analysis;
  my $dna_db = $self->get_database_by_name('dna_db');
  my $dba = $self->get_database_by_name('output_db', $dna_db);
  $self->hrdb_set_con($dba, 'output_db');

  $self->query($self->fetch_sequence($self->input_id, $dba));
  # get all candidates lincRNA 
  my $lincrna_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->RNA_DB);
  if (@$lincrna_genes) {
    $self->say_with_header('We have '.scalar(@$lincrna_genes).' multi_transcript lincRNA candidates from RNAseq stage.');
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
            -query => $self->query,
            -analysis => $self->analysis,
            -genes => $lincrna_genes,
            -output_biotype => $self->param_required('biotype_output'),
            -max_transcripts_per_cluster => 20,
            -min_short_intron_len => 1,
            -max_short_intron_len => 15,
            -blessed_biotypes => {},
           );
    $self->runnable($runnable);
  }
  else {
    $self->input_job->autoflow(0);
    $self->complete_early('There is no gene to process');
  }
}


sub run {
  my $self = shift; 
 
  print  "\nCOLLAPSING RNAseq models (removing duplications) BY GENEBUILDER. \n Afterwards feed lincRNA pipeline.\n";
  my $tmp = $self->runnable->[0]; 
  print  "\n 1. Running GeneBuilder for " . scalar( @{$tmp->input_genes} ). " RNAseq models...\n";  
  foreach my $runnable (@{$self->runnable}) {
    $runnable->run;
    $self->output($runnable->output);
    print "Number of collapsed models: " . scalar(@{$self->output}) . "\n"; 
  }
}

sub write_output{
  my ($self) = @_; 

  print  "\nWRITING RESULTS IN OUTPUT DB... " . "\n";
  my $lincrna_ga  = $self->hrdb_get_con('output_db')->get_GeneAdaptor;
  my $genes_to_write = $self->output;
  print  "***HAVE ". scalar(@$genes_to_write) ." GENE(S) TO WRITE IN TOTAL (INCLUDING VALID AND REJECTED lincRNAs).\n";

  my $logic_name_to_be = "lincRNA_noDuplication";
  ## I will create a new gene without translations and only one transcript to be stored under different analysis ##
  # Make an analysis object (used later for storing stuff in the db)
  my $analysis = Bio::EnsEMBL::Analysis->new(
                                         -logic_name => $logic_name_to_be,
                                         -displayable => 1
                                         );
  my $slice = $self->query;
  foreach my $gene (@$genes_to_write) {
    empty_Gene($gene);
    $gene->analysis($analysis);
    $lincrna_ga->store($gene);
  } 
}

sub RNA_DB {
  my ($self, $arg) = @_;
  if($arg){
    $self->param('RNA_DB', $arg);
  }
  return $self->param('RNA_DB');
}


1;
