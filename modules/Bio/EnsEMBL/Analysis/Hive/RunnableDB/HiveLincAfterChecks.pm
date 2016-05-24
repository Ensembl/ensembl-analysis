package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincAfterChecks;
        
use warnings;
use vars qw(@ISA);
use strict;
use Data::Dumper;

use feature 'say';
use File::Basename;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub fetch_input {
  my $self = shift;
  
  $self->param_required('file');
      	
	# my $biotype_tmp = $self->param_required('Final_BIOTYPE_TO_CHECK');
	my $biotype_tmp = 'lincRNA_withReg';
	print "pak\n";
	my $output_file = $self->param('file');
	# $self->param('file'); 
	my $short_genes = $self->check_length_of_genes($biotype_tmp, $output_file);
	# file_for_output

	 
  return 1;
}

sub run {
  my $self = shift;
  
  


  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}



sub check_length_of_genes {
	my ( $self, $href, $output_file ) = @_;
	my @genes_to_fetch;
	my $set_db  = $self->hrdb_get_dba( $self->param('regulation_reform_db') );
	my $dna_dba = $self->hrdb_get_dba( $self->param('reference_db') );
	if ($dna_dba) {
		$set_db->dnadb($dna_dba);
	}

  my $sa = $dna_dba->get_SliceAdaptor;
  my $count =0;
  my $how_many = 0;
  
  # my $test_id = $self->param('iid');  
  # print "-----> $test_id  \n";
  # my $slice = $self->fetch_sequence($test_id, $set_db, undef, undef, 'lincRNA_output_db')  ;
  
  open (MYFILE, ">" , $output_file) or die "Couldn't open: $!"; 
  foreach my $slice (@{$sa->fetch_all('toplevel')}) {
    foreach my $g ( @{ $slice->get_all_Genes } ) {
  	  next if ( $g->biotype ne "lincRNA" );  
      for my $t ( @{$g->get_all_Transcripts } ) {  
        my $stable_id = $g->display_id; 
        my $tl_length = $t->length ; 
        if ( $tl_length < 200 ) { 
          print MYFILE "WARNING:: $stable_id and $tl_length \n";
        }
      }
    }
  }
  close MYFILE;

return $how_many;
}



sub check_number_of_exons {
	my ( $self, $href ) = @_;
	my @genes_to_fetch;
	my $set_db  = $self->hrdb_get_dba( $self->param('regulation_reform_db') );
	my $dna_dba = $self->hrdb_get_dba( $self->param('reference_db') );
	if ($dna_dba) {
		$set_db->dnadb($dna_dba);
	}

  my $sa = $dna_dba->get_SliceAdaptor;
  my $count =0;
  my $how_many = 0;
  foreach my $slice (@{$sa->fetch_all('toplevel')}) {
    foreach my $g ( @{ $slice->get_all_Genes } ) {
  	  next if ( $g->biotype ne "lincRNA" ); 
      my $stable_id = $g->display_id; 
  	  my @exons = @{ $g->get_all_Exons }; 
  	  if (scalar(@exons) < 2){
        print "WARNING:: $stable_id and " . scalar(@exons) . " number of exons\n";
      }
    } 
  } 
return $how_many; 
}


1;



