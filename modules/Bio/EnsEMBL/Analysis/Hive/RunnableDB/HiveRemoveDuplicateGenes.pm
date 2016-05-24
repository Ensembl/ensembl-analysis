=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveDuplicateGenes;


use warnings ;
use vars qw(@ISA);
use strict;
use Data::Dumper;

use Bio::EnsEMBL::Hive::Utils ('destringify');
use Bio::EnsEMBL::Analysis; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

# @ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild Bio::EnsEMBL::Analysis::RunnableDB);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::HiveRemoveDuplicateGenes
  Function  : instatiates a lincRNAFinder object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::HiveRemoveDuplicateGenes 
  Exceptions: 
  Example   : 

=cut


sub fetch_input{
  my ($self) = @_;

  # set up config
  $self->hive_set_config;
  
print "DEBUG::fetch_input!!\n"; 
  # Get lincRNA candidate genes and break each of them down into single-transcript genes:
print  "DEBUG:: fetch_input:: dump the object START_FROM" .  "\n"; 
  # get all candidates lincRNA 
  my $lincrna_genes = $self->get_genes_of_biotypes_by_db_hash_ref($self->RNA_DB);
  ### ### foreach my $gene (@{$lincrna_genes}) {
  ### ###   print  $gene->display_id, "\t", $gene->adaptor->dbc->dbname, , "\t", $gene->adaptor->dbc->host, "\n";
  ### ### }

  print  "We have ". scalar(@$lincrna_genes) . " multi_transcript lincRNA candidates from RNAseq stage.\n" ;  
  print  "analysis : " . $self->analysis . "\n";
  $self->param_required('biotype_output'); 
  my $output_biotype = $self->param('biotype_output');   #"human_rnaseq"; # hard coded now. need to change... 
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
          -query => $self->query,
          -analysis => $self->analysis,
          -genes => $lincrna_genes, #  \@genes_for_build, 
          -output_biotype =>  $output_biotype, 
          -max_transcripts_per_cluster => 20, 
          -min_short_intron_len => 1,
          -max_short_intron_len => 15,
          -blessed_biotypes => {} , 
         );
  $self->runnable($runnable);
  # $self->param('linc_rna_genes') = @$lincrna_genes; 
}


sub run {
  my $self = shift; 
 
  print  "\nCOLLAPSING RNAseq models (removing duplications) BY GENEBUILDER. \n Afterwards feed lincRNA pipeline.\n";

  my $tmp = $self->runnable->[0]; 
  # print "number of elements: " . scalar(@tmp) . "\n";

  print  "\n 1. Running GeneBuilder for " . scalar( @{$tmp->input_genes} ). " unclustered lincRNAs...\n";  

  foreach my $runnable (@{$self->runnable}) {
    $runnable->run;
    $self->output($runnable->output);
    print "DEBUG::" . scalar(@{$self->output}) . "\n"; 
  }
  
  #
  # First genebuilder run with unclustered lincRNAs
  #
  # print "DEBUG::HIVElincRNaEvaluator::dumper:: self:: " . Dumper($self) . "\n"; 
  # print "DEBUG::HIVElincRNaEvaluator::dumper:: gb:: " . Dumper($gb) . "\n"; 

  
   # my @output_clustered; 
   
# print "DEBUG::HIVElincRNaEvaluator::dumper::" . Dumper($gb) . "\n"; 
   # push @output_clustered, @{$gb->output()} ;  
   # $self->output( \@output_clustered );  #  This is just to store the lincRNA genes so test_RunnableDB script can find them.

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
  foreach my $gene(@genes_to_write){  
    for ( @{$gene->get_all_Transcripts} ) {    
     #  print  "DEBUG::write $_ : ".$_->seq_region_start . " " . $_->biotype . " " . $_->translation()->seq() . " " . $_->translation()->start() . " gene: " . $gene->dbID . " " . $gene->biotype . " " . $_->strand . " " . $gene->strand ."\n"; 
      my $end_exon = $_->end_Exon;
      print  "DEBUG::write::end_exon:: $end_exon  \n";  
    } 
    my $logic_name_to_be = "lincRNA_set_test_3"; 
    # my @t = @{ $gene->get_all_Transcripts}; 
    # my $strand_to_use  =  $gene->strand;    
    # my $biotype_to_use = "empty_something_is_wrong"; 

    ## I will create a new gene without translations and only one transcript to be stored under different analysis ##
    # Make an analysis object (used later for storing stuff in the db)
    my $analysis = Bio::EnsEMBL::Analysis->new(
                                           -logic_name => $logic_name_to_be,
                                           -displayable => 1
                                           );    
    $gene->status(undef); 
    $gene->analysis($self->analysis);   
    eval{
      $lincrna_ga->store($gene);
    };
    if($@){
      warning("Failed to write gene ".id($gene)." ".coord_string($gene)." $@");
    }else{
      $sucessful_count++;
      print  "STORED LINCRNA GENE ".$gene->dbID. " "  . $gene->biotype  . " " .  $gene->strand . "\n";
    }
  } 
  
  # this check was added because I had problems with few genes. 
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
  	}
    $self->param('fail_delete_features', \@genes_to_write);
    $self->throw($@);
  }
}


sub post_cleanup {
  my $self = shift;
  
  if ($self->param_is_defined('fail_delete_features')) {
    my $dba = $self->hrdb_get_con('lincRNA_output_db');
    my $gene_adaptor = $dba->get_GeneAdaptor;
    foreach my $gene (@{$self->param('fail_delete_features')}) {
      eval {
         $gene_adaptor->remove($gene);
      };
      if ($@) {
        $self->throw('Could not cleanup the mess for these dbIDs: '.join(', ', @{$self->param('fail_delete_features')}));
      }
    }
  }
  return 1;
}




####### START KB15 ADD METHOD
sub get_genes_of_biotypes_by_db_hash_ref { 
  my ($self, $href) = @_;

  my %dbnames_2_biotypes = %$href ; 

  ### ### print  "DEBUG::get_genes_of_biotypes_by_db_hash_ref::Get genes " . scalar(keys %dbnames_2_biotypes) . "\n dumper:" . Dumper(%$href) . "\n"; 

  my @genes_to_fetch;  
  foreach my $db_hash_key ( keys %dbnames_2_biotypes )  {
    print  "DEBUG::get_genes_of_biotypes_by_db_hash_ref::1 $db_hash_key\n";  # <--- name of the database to use

    my @biotypes_to_fetch = @{$dbnames_2_biotypes{$db_hash_key}};  
   
    # my $set_db = $self->hrdb_get_dba($self->param($db_hash_key));
    my $set_db = $self->hrdb_get_dba($self->param('source_cdna_db'));
    my $dna_dba = $self->hrdb_get_dba($self->param('reference_db'));
    if($dna_dba) {
      $set_db->dnadb($dna_dba);
    }

    print  "DEBUG::get_genes_of_biotypes_by_db_hash_ref::2 $db_hash_key\n";
    my $test_id = $self->param('iid');
    # my $test_id = "chromosome:GRCh38:13:1:89625480:1"; 
    my $slice = $self->fetch_sequence($test_id, $set_db, undef, undef, $db_hash_key)  ;
print "-----> $test_id -- $set_db " .  $set_db->dbc->dbname  . " \n";
    # implementation of fetch_all_biotypes ....  
    my $fetch_all_biotypes_flag ; 
    foreach my $biotype  ( @biotypes_to_fetch ) {   
      if ($biotype=~m/fetch_all_biotypes/ ) {    
        $fetch_all_biotypes_flag = 1 ; 
      }
    }  
    if ( $fetch_all_biotypes_flag ) {  
         print  "fetching ALL biotypes for slice out of db $db_hash_key :\n" ; 
         my $genes = $slice->get_all_Genes(undef,undef,1) ; 
         push @genes_to_fetch, @$genes;  
         my %tmp ; 
         for ( @$genes ) {  
           $tmp{$_->biotype}++; 
         }  
         foreach ( keys %tmp ) {  
           print  "found $_ $tmp{$_}\n" ; 
         }  
         print  scalar(@genes_to_fetch) . " genes fetched in total\n" ; 
    } else { 
      foreach my $biotype  ( @biotypes_to_fetch ) { 
      	# $biotype = "best"; 
         my $genes = $slice->get_all_Genes_by_type($biotype,undef,1);
         if ( @$genes == 0 ) {
           warning("No genes of biotype $biotype found in $set_db\n");
         } 
         # if ( $self->verbose ) { 
         print  "$db_hash_key [ " . $set_db->dbc->dbname  . " ] Retrieved ".@$genes." of type ".$biotype."\n";
         print  "DEBUG::HiveLincRNA::get_genes_of_biotypes_by_db_hash_ref: " . $db_hash_key . " Retrieved ".@$genes." of biotype ".$biotype."\n";
         # }
         push @genes_to_fetch, @$genes;
      }  
    }  
  } 

  return \@genes_to_fetch;
}
####### END KB15 ADD METHOD

sub check_if_all_stored_correctly {
  my ($self, $href) = @_;

  my $set_db = $self->hrdb_get_dba($self->param('lincRNA_output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('reference_db'));
  if($dna_dba) {
    $set_db->dnadb($dna_dba);
  }
  
  my $test_id = $self->param('iid'); 
  # print "-----> $test_id  \n";
  my $slice = $self->fetch_sequence($test_id, $set_db, undef, undef, 'lincRNA_output_db')  ;
  # implementation of fetch_all_biotypes ....  
  print  "check if all genes are fine!! \n" ; 
  my $genes = $slice->get_all_Genes(undef,undef,1) ; 
	return "yes";
}



# HIVE check
sub hive_set_config {
  my $self = shift;

  # Throw is these aren't present as they should both be defined
  unless($self->param_is_defined('logic_name') && $self->param_is_defined('module')) {
    throw("You must define 'logic_name' and 'module' in the parameters hash of your analysis in the pipeline config file, ".
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
      throw("You have a key defined in the config_settings hash (in the analysis hash in the pipeline config) that does ".
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
