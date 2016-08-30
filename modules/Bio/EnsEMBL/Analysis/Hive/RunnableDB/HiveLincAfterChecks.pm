package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincAfterChecks;
        
use warnings;
use vars qw(@ISA);
use strict;
use Data::Dumper;

use feature 'say';
use File::Basename;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Analysis::Tools::LincRNA qw(get_genes_of_biotypes) ;  
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub fetch_input {
	my ($self) = @_;
	$self->hive_set_config;
  $self->param_required('file');      	

	# get data
	my $set_db  = $self->hrdb_get_dba( $self->param('regulation_reform_db') );
	my $dna_dba = $self->hrdb_get_dba( $self->param('reference_db') );
  $set_db->dnadb($dna_dba);
  $self->hrdb_set_con($set_db, 'regulation_reform_db');

	my $biotype_tmp = 'lincRNA_withReg';
	my $dbName_tmp = 'regulation_reform_db';
	my $genes_to_check  = $self->get_genes_of_biotypes($biotype_tmp, $dbName_tmp);
	print "fetch input :: number of lincRNA genes to check: " . scalar(@$genes_to_check) . "\n";
  $self->param('genes', $genes_to_check);
  return 1;
}

sub run {
  my $self = shift;

  my $genes_to_process = $self->param('genes');
	my $output_file = $self->param('file');
	my $short_genes = $self->check_length_of_genes($output_file, $genes_to_process);
 
  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}



sub check_length_of_genes {
	my ( $self, $output_file, $genes_to_process) = @_;
	my @genes_to_fetch;
  my $count =0;
  my $how_many = 0;

  open (MYFILE, ">" , $output_file) or die "Couldn't open: $!"; 
  foreach my $g ( @{ $genes_to_process } ) {
      for my $t ( @{$g->get_all_Transcripts } ) {  
        my $tl_length = $t->length ;
        my $exons_count = check_number_of_exons($g); 
        if ( ( $tl_length < 200 ) or ($exons_count < 2) ) { 
          my $some_id = ""; 
          if ($g->stable_id()) {
            $some_id = $g->stable_id(); 
          }else {
      	    $some_id = $g->display_id;
          } 
          print MYFILE "WARNING:: $some_id :: $tl_length :: $exons_count \n";
        }
      }
  }
  close MYFILE;
  return $how_many;
}



sub check_number_of_exons {
	my ( $hgene) = shift;
  my @exons = @{ $hgene->get_all_Exons }; 
  my $how_many = scalar(@exons);
  return $how_many; 
}




#### Standard Functions ###

=head2 hive_set_config 

  Function  : loop through parameters and set them. All parameters of basic config (config_settings) will be accessible via param
  Returntype:  

=cut

sub hive_set_config {
	my $self = shift;

	# Throw is these aren't present as they should both be defined
	unless ( $self->param_is_defined('logic_name')
		&& $self->param_is_defined('module') )
	{
		warn("You must define 'logic_name' and 'module' in the parameters hash of your analysis in the pipeline config file, "
				. "even if they are already defined in the analysis hash itself. This is because the hive will not allow the runnableDB "
				. "to read values of the analysis hash unless they are in the parameters hash. However we need to have a logic name to "
				. "write the genes to and this should also include the module name even if it isn't strictly necessary"
		);
	}

  # Make an analysis object and set it, this will allow the module to write to the output db
	my $analysis = new Bio::EnsEMBL::Analysis(
		-logic_name => $self->param('logic_name'),
		-module     => $self->param('module'),
	);
	$self->analysis($analysis);

  # Now loop through all the keys in the parameters hash and set anything that can be set
	my $config_hash = $self->param('config_settings');
	foreach my $config_key ( keys( %{$config_hash} ) ) {
		if ( defined &$config_key ) {
			$self->$config_key( $config_hash->{$config_key} );
		}
		else {
			warn("You have a key defined in the config_settings hash (in the analysis hash in the pipeline config) that does "
					. "not have a corresponding getter/setter subroutine. Either remove the key or add the getter/setter. Offending "
					. "key:\n"
					. $config_key );
		}
	}
}



1;



