=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincAfterChecks - 

=head1 SYNOPSIS

Check the final lincRNA dataset and classify or suggests a classification for each model

=head1 DESCRIPTION


=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincAfterChecks;
        
use warnings;
use strict;
use File::Basename;
use Bio::EnsEMBL::Analysis::Tools::LincRNA qw(get_genes_of_biotypes) ;  
use Bio::EnsEMBL::Analysis;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my ($self) = @_;
  $self->create_analysis;
  $self->param_required('file_l'); 
  $self->param_required('file_b'); 
  $self->param_required('assembly_name');
  $self->param_required('update_database');

  # get genes from lincRNA database
  my $set_db  = $self->hrdb_get_dba( $self->param('output_db') );
  my $dna_dba = $self->hrdb_get_dba( $self->param('dna_db') );
  $set_db->dnadb($dna_dba);
  $self->hrdb_set_con($set_db, 'output_db');
  my $biotype_lincRNA =   $self->Final_BIOTYPE_TO_CHECK;  
  my $dbName_to_check = 'output_db';
  my $genes_to_check  = $self->get_genes_of_biotypes($biotype_lincRNA, $dbName_to_check);

  # get genes from core database
  my $core_dba = $self->hrdb_get_dba( $self->param('dna_db') ); 
  $core_dba->dnadb($dna_dba);
  $self->hrdb_set_con($core_dba, 'dna_db');
  my $biotype_core = 'fetch_all_biotypes'; 
  my $core_to_check = 'dna_db';
  my $genes_from_core   = $self->get_genes_of_biotypes($biotype_core, $core_to_check);
  
  print "### " . "fetch input :: biotype: $biotype_lincRNA has number of lincRNA genes to check: " . scalar(@$genes_to_check) . "\n";
  print "### " . "fetch input :: biotype: $biotype_core number of core genes to check: " . scalar(@$genes_from_core) . "\n";
	
  $self->param('lincRNA_genes', $genes_to_check);
  $self->param('core_genes', $genes_from_core);
  return 1;
}


sub run {
  my $self = shift;

  my $genes_to_process = $self->param('lincRNA_genes');
  my $standard_geneset = $self->param('core_genes');
  my $output_file_l = $self->param('file_l');
  my $output_file_b = $self->param('file_b');
  my $short_genes = $self->check_length_of_genes($output_file_l, $genes_to_process);
  my $biotypes_updates = $self->update_biotypes($output_file_b, $genes_to_process );
  
  return 1;
}


sub write_output {
  my $self = shift;
  return 1;
}



#### Functions ###

sub check_length_of_genes {
  my ( $self, $output_file, $genes_to_process) = @_;
  my @genes_to_fetch;
  my $count =0;
  my $how_many = 0;
  print "### my output file for length info about the genes is: $output_file \n";
  open (MYFILE, ">" , $output_file) or die "Couldn't open: $!"; 
  print MYFILE "## WARNING stable_id(orGeneDisplayID) gl_length exons_count \n";
  foreach my $g ( @{ $genes_to_process } ) {
  	my $gl_length = $g->length ;
    my $exons_count = check_number_of_exons($g); 
    if ( ( $gl_length < 200 ) or ($exons_count < 3) ) { 
      $how_many = $how_many+1; 
      my $some_id = ""; 
      if ($g->stable_id()) {
        $some_id = $g->stable_id(); 
      }else {
        $some_id = $g->display_id;
      } 
      print MYFILE "## WARNING $some_id $gl_length $exons_count \n";
      my $tmp_biotype = "putative_truncated_lincRNA";
      print MYFILE "UPDATE gene, transcript SET transcript.biotype=\'$tmp_biotype\' , gene.biotype=\'$tmp_biotype\' WHERE gene.gene_id = $some_id AND gene.gene_id = transcript.gene_id ;\n"; 
    }
  }
  close MYFILE;
  return $how_many;
}


sub check_number_of_exons {
  my ( $hgene) = shift;
  my @exons = @{ $hgene->get_all_Exons }; 
  my $how_many_exons = scalar(@exons);
  return $how_many_exons; 
}


sub update_biotypes {
  my ( $self, $output_file, $lincRNA_geneset) = @_;

  open (MYFILE_B, ">" , $output_file) or die "Couldn't open: $!"; 
  my $assembly_version = $self->param('assembly_name');
  my $core_db = $self->hrdb_get_con('dna_db');
  my $core_db_adaptor = $core_db->get_SliceAdaptor;
  foreach my $g_lincRNA ( @{ $lincRNA_geneset } ) {
  	my $protein_switch = 0; # this is in case there are many genes that overlap with the lincRNA
  	my $tmp_biotype = "empty"; 
  	my $g_lincRNA_id = $g_lincRNA->dbID; 
  	my $g_start = $g_lincRNA->seq_region_start;
    my $g_end = $g_lincRNA->seq_region_end;
    my $chr_name = $g_lincRNA->seq_region_name; 
    my $slice_to_check = $core_db_adaptor ->fetch_by_region( 'toplevel', $chr_name, $g_start, $g_end, undef, $assembly_version );
    # print "### I have " . scalar(@{ $slice_to_check->get_all_Genes }) . " genes \n"; 
    if (scalar(@{ $slice_to_check->get_all_Genes }) == 0 ) {
    	# print "### " . "it is not overlap anything and will be marked as lincRNA\n" ;
    	$tmp_biotype = "lincRNA";             
    } elsif (scalar(@{ $slice_to_check->get_all_Genes }) == 1 ) {
      foreach my $g_core (@{ $slice_to_check->get_all_Genes } ) {	  
        if ($g_core->biotype eq "protein_coding" ) { 
      	  if ($g_lincRNA->strand ne $g_core->strand) {
            # print "### " . $g_lincRNA->dbID . " " . $g_lincRNA->seq_region_start . " " . $g_lincRNA->seq_region_end . " " . $g_lincRNA->seq_region_name . " " . $g_lincRNA->biotype  . " gBiotype:" . $g_core->biotype  .  " linc_strand: " . $g_lincRNA->strand . " gStrand: " . $g_core->strand  .  " antisense\n";
            $tmp_biotype = "antisense"; 
          } else {
            # print "### " . $g_lincRNA->dbID . " " . $g_lincRNA->seq_region_start . " " . $g_lincRNA->seq_region_end . " " . $g_lincRNA->seq_region_name . " " . $g_lincRNA->biotype  . " gBiotype:" . $g_core->biotype  .  " linc_strand: " . $g_lincRNA->strand . " gStrand: " . $g_core->strand  .  " senseANDneedsAttention\n";      		
            $tmp_biotype = "sense"; 
      	  }
        } else {
        	# print "### " . "Biotype is not protein coding, but " . $g_core->biotype . " , so I will use a lincRNA biotype \n";
          $tmp_biotype = "lincRNA"; 
        }
      }
    } else {
      # in case there are more than 2 genes in that region. 
      my $count_protein_coding_overlaps = 0;
      foreach my $g_core (@{ $slice_to_check->get_all_Genes } ) {	  
        if ($g_core->biotype eq "protein_coding" ) {  
          $protein_switch = 1 +$protein_switch; # there is a protein_coding, if more than one, requires attention! 
      	  if ($g_lincRNA->strand ne $g_core->strand) {
            print "### " . $g_lincRNA->dbID . " " . $g_lincRNA->seq_region_start . " " . $g_lincRNA->seq_region_end . " " . $g_lincRNA->seq_region_name . " " . $g_lincRNA->biotype  . " gBiotype:" . $g_core->biotype  .  " linc_strand: " . $g_lincRNA->strand . " gStrand: " . $g_core->strand  .  " antisense\n";
            $tmp_biotype = "antisense";
          } else {
            print "### " . $g_lincRNA->dbID . " " . $g_lincRNA->seq_region_start . " " . $g_lincRNA->seq_region_end . " " . $g_lincRNA->seq_region_name . " " . $g_lincRNA->biotype  . " gBiotype:" . $g_core->biotype  .  " linc_strand: " . $g_lincRNA->strand . " gStrand: " . $g_core->strand  .  " senseANDneedsAttention\n";      		
            $tmp_biotype = "sense";
      	  }
        } 
      }
      if ($protein_switch == 0) { 
        print "### " . "Biotypes of all overlapping genes are not protein coding, so I will put lincRNA as biotype \n";
        $tmp_biotype = "lincRNA";
      } elsif ($protein_switch == 1) {
      	print "### " . "There is only one protein coding that overlaps and will use the previous, " . $tmp_biotype .  " biotype\n";
        $tmp_biotype = $tmp_biotype;
      } elsif ($protein_switch > 1) {
      	print "### " . "There are many protein coding that overlap and the region requires attention, will set the biotype to warning " . "\n";
        $tmp_biotype = "CHECK_ME_MANUALLY"; 
      } else {
      	print "### " . " Don't expect to be here! \n";
      }
    } # lincRNA gene loop done
    print MYFILE_B "UPDATE gene, transcript SET transcript.biotype=\'$tmp_biotype\' , gene.biotype=\'$tmp_biotype\' WHERE gene.gene_id = $g_lincRNA_id AND gene.gene_id = transcript.gene_id ;\n"; 
  } 
  close MYFILE_B; 

  # update the biotypes with the mysql queries... if the user is fine with it! 
  if ( $self->param('update_database') eq "yes" ) {
    my $db_lincRNA_out = $self->hrdb_get_dba( $self->param('output_db') );
    my $command_to_upload = 'mysql -h  ' . $db_lincRNA_out->dbc->host . ' -P ' . $db_lincRNA_out->dbc->port .' -u ' . $db_lincRNA_out->dbc->user . ' -p' . $db_lincRNA_out->dbc->password . '  ' . $db_lincRNA_out->dbc->dbname . ' < ' .  $output_file ; 
    print "I will apply the suggested changes to the database... " . $command_to_upload . "\n"; 
    $self->throw("can't update the biotypes check the file: ". $output_file)
      unless ((system ($command_to_upload)) == 0);
  } else {
  	print "Check the files and think if you need to apply the biotypes...\n"; 
  }
} 




#### Standard Functions ###


sub Final_BIOTYPE_TO_CHECK{
  my ($self, $arg) = @_;
  if($arg){
    $self->param('Final_BIOTYPE_TO_CHECK', $arg);
  }
  return $self->param('Final_BIOTYPE_TO_CHECK');
}  


1;



