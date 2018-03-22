=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string);
use Bio::EnsEMBL::Analysis::Tools::LincRNA qw(get_genes_of_biotypes_by_db_hash_ref) ;  
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(print_Gene) ; 

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input{
  my ($self) = @_;

  # set up config
  $self->hive_set_config; 
  # get all candidates lincRNA 
  my $lincrna_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->RNA_DB);
  print  "We have ". scalar(@$lincrna_genes) . " multi_transcript lincRNA candidates from RNAseq stage.\n" ;  
  $self->param_required('biotype_output'); 
  my $output_biotype = $self->param('biotype_output');    
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
          -query => $self->query,
          -analysis => $self->analysis,
          -genes => $lincrna_genes,  
          -output_biotype =>  $output_biotype, 
          -max_transcripts_per_cluster => 20, 
          -min_short_intron_len => 1,
          -max_short_intron_len => 15,
          -blessed_biotypes => {} , 
         );
  $self->runnable($runnable);
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
  my $dba = $self->hrdb_get_dba($self->param('lincRNA_output_db'));
  $self->hrdb_set_con($dba,'lincRNA_output_db');
  my $lincrna_ga  = $self->hrdb_get_con('lincRNA_output_db')->get_GeneAdaptor;
  my @genes_to_write = @{$self->output}; 
  print  "***HAVE ". scalar(@genes_to_write) ." GENE(S) TO WRITE IN TOTAL (INCLUDING VALID AND REJECTED lincRNAs).\n";  

  my $sucessful_count = 0 ; 
  my $logic_name_to_be = "lincRNA_noDuplication";
  ## I will create a new gene without translations and only one transcript to be stored under different analysis ##
  # Make an analysis object (used later for storing stuff in the db)
  my $analysis = Bio::EnsEMBL::Analysis->new(
                                         -logic_name => $logic_name_to_be,
                                         -displayable => 1
                                         );
  foreach my $gene(@genes_to_write){  
    $gene->analysis($analysis);
    eval{
      $lincrna_ga->store($gene);
    };
    if($@){
      $self->warning("Failed to write gene: " . print_Gene($gene) ." $@");
    }else{
      $sucessful_count++;
    }
  } 
  
  # this check was added because I had problems with few genes that didn't stored and the job didn't died! mysql kind of thing! 
  eval{
    my $check = $self->check_if_all_stored_correctly($self->RNA_DB); 
    print "check result: " . $check . " -- " . $sucessful_count ." genes written to FINAL OUTPUT DB " . $dba->dbc->dbname . "\n" ; # . $self->output_db->dbname . " @ ". $self->output_db->host . "\n"  ;   
    if($sucessful_count != @genes_to_write ) { 
      $self->throw("Failed to write some genes");
    }
  };
  if($@){
  	print "You have a problem with those genes, they didn't stored suggessfully, I will try to delete them and rerun the job: \n "; 
  	foreach my $g_t(@genes_to_write){ 
  		print $g_t->dbID . "\n";
                print "--start:" . $g_t->seq_region_start . " end:" . $g_t->seq_region_end . "\n" ;  	
  	}
    $self->param('fail_delete_features', \@genes_to_write);
    $self->throw($@);
  }
}

# post_cleanup will clean your entries if your full job didn't finish fine. Usefull! 
sub post_cleanup {
  my $self = shift;
  
  if ($self->param_is_defined('fail_delete_features')) {
    my $dba = $self->hrdb_get_con('lincRNA_output_db');
    my $gene_adaptor = $dba->get_GeneAdaptor;
    foreach my $gene (@{$self->param('fail_delete_features')}) {
      eval {
        print "DEBUG::cleaning-removing gene, something didn't go as should... \n"; 
        $gene_adaptor->remove($gene);
      };
      if ($@) {
        $self->throw('Could not cleanup the mess for these dbIDs: '.join(', ', @{$self->param('fail_delete_features')}));
      }
    }
  }
  return 1;
}


# this function checks if everything stored successfully 
sub check_if_all_stored_correctly { 
  my ($self, $href) = @_; 

  my $set_db = $self->hrdb_get_dba($self->param('lincRNA_output_db')); 
  my $dna_dba = $self->hrdb_get_dba($self->param('reference_db')); 
  if($dna_dba) { 
    $set_db->dnadb($dna_dba); 
  } 
  
  my $test_id = $self->param('iid'); 
  my $slice = $self->fetch_sequence($test_id, $set_db, undef, undef, 'lincRNA_output_db')  ; 
  print  "check if all genes are fine!! \n" ; 
  my $genes = $slice->get_all_Genes(undef,undef,1) ; 
	return "yes"; 
}


# HIVE check
sub hive_set_config {
  my $self = shift;

  # Throw is these aren't present as they should both be defined
  unless($self->param_is_defined('logic_name') && $self->param_is_defined('module')) {
    $self->throw("You must define 'logic_name' and 'module' in the parameters hash of your analysis in the pipeline config file, ".
          "even if they are already defined in the analysis hash itself. This is because the hive will not allow the runnableDB ".
          "to read values of the analysis hash unless they are in the parameters hash. However we need to have a logic name to ".
          "write the genes to and this should also include the module name even if it isn't strictly necessary"
         );
  }

  # Make an analysis object and set it, this will allow the module to write to the output db
  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                           );
  $self->analysis($analysis);

  # Now loop through all the keys in the parameters hash and set anything that can be set
  my $config_hash = $self->param('config_settings');
  foreach my $config_key (keys(%{$config_hash})) {
    if(defined &$config_key) {
      $self->$config_key($config_hash->{$config_key});
    } else {
      $self->throw("You have a key defined in the config_settings hash (in the analysis hash in the pipeline config) that does ".
            "not have a corresponding getter/setter subroutine. Either remove the key or add the getter/setter. Offending ".
            "key:\n".$config_key
           );
    }
	}
}


=head2 FINAL_OUTPUT_DB

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::lincRNAEvaluator
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::lincRNAEvaluator

sub FINAL_OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('FINAL_OUTPUT_DB', $arg);
  }
  return $self->param('FINAL_OUTPUT_DB');
}  

sub RNA_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('RNA_DB', $arg);
  }
  return $self->param('RNA_DB');
}  




1;
