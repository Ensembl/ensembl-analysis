#This is to be a baseclass for the genebuild code to have various
#object like methods which genebuild modules want. GeneBuild code
#should ideally inherit from this


package Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Config::Databases qw(DATABASES DNA_DBNAME);
use Bio::EnsEMBL::Analysis::RunnableDB;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);

#first of all we just need methods to get databases from the
#database config


sub database_hash{
  my ($self, $name, $db) = @_;
  if(!$self->{'db_hash'}){
    $self->{'db_hash'} = {};
  }
  if($name && $db){
    $self->{'db_hash'}->{$name} = $db;
  }
  return $self->{'db_hash'};
}

=head2 get_dbadaptor

  Arg [1]   : String - key of database hash
  Arg [2]   : return a pipeline db adaptor flag
  Function  : Returns a Bio::EnsEMBL::DBSQL::DBAdaptor for a given hash key.
              or a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor if requested
              Requires proper configuration of
              Bio::EnsEMBL::Analysis::Config::Databases

  Returntype: Bio::EnsEMBL:DBSQL::DBAdaptor or Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Exceptions: throw if key can't be found in Databases.pm

=cut

sub get_dbadaptor{
  my ($self, $name, $use_pipeline_adaptor, $not_use_dna_database) = @_;
  my $hash = $self->database_hash;
  my $db;
  if(!$hash->{$name}){
    if (exists $DATABASES->{$name}) {
      my $constructor_args = $DATABASES->{$name};

      foreach my $arg ( qw ( -user -port -host -dbname) ) {
        unless ( $$constructor_args{$arg}){
          throw ("Database-connection-details not properly configured : Argument : $arg missing in Databases.pm for $name \n") ;
        }
      }
      if ( $use_pipeline_adaptor ) {
         if ( $use_pipeline_adaptor == 1 || $use_pipeline_adaptor eq "pipeline") {
            $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
							    %$constructor_args,
	  						   );
         } elsif ( $use_pipeline_adaptor == 2 || $use_pipeline_adaptor =~m/compara/ ) {
         	require Bio::EnsEMBL::Compara::DBSQL::DBAdaptor; 
            unless ( $$constructor_args{'-species'}){  
               throw("need species !\n"); 
            } 
            $db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
                                                            %$constructor_args,
                                                              );
         } elsif ( $use_pipeline_adaptor == 3 || $use_pipeline_adaptor =~m/functgenomics/i || $use_pipeline_adaptor =~m/funcgen/i) {
         	require Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
            $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
                                                            %$constructor_args,
                                                              );
         }
      } else {
	$db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
						  %$constructor_args,
						 );
      }
      if($name ne $DNA_DBNAME ){
        if (length($DNA_DBNAME) ne 0 ){
          if ( $not_use_dna_database ) {
            # if two different species are considered, the not_use_dna_database is set to 1 to avoid adding the second species to the first one
          } else { 
           # there's a little danger if we have multiple diffent species in our "Databases.pm" file. 
           # We need to avoid that the wrong dna db is attached, ie a mouse core with a human dna db. 
         

           my $dnadb = $self->get_dbadaptor($DNA_DBNAME);  

           # get species name for DNA_DBNAME and OTHER db 
           my $core_db_species = $db->get_MetaContainer->get_Species->binomial();
           my $dna_db_species = $dnadb->get_MetaContainer()->get_Species->binomial(); 


           # get default asm for DNA_DBNAME and OTHER db 
           my $core_db_asm = $db->get_MetaContainer->get_default_assembly();   
           my $dna_db_asm = $dnadb->get_MetaContainer->get_default_assembly();     

           unless ( $core_db_asm eq $dna_db_asm ) { 
               throw("you try to add  a DNA_DB with assembly $dna_db_asm to a core/cdna/otherfeatures DB with assembly $core_db_asm ...\n\t".
                    "that's incompatbile. try to not use any DNA_DATABASE name in Analysis/Config/Databases.pm\n" ); 
           }   

           unless ( $core_db_species eq $dna_db_species ) { 
               throw("you try to add  a DNA_DB with species $dna_db_species to a core database with speices : $core_db_species - this does not work\n\t".
                     "that's incompatbile. try to not use any DNA_DATABASE name in Analysis/Config/Databases.pm\n"); 
           }  

           $db->dnadb($dnadb); 
           
         }
        }else{
          warning("You haven't defined a DNA_DBNAME in Config/Databases.pm");
        }
      }
      $self->database_hash($name, $db);
    } else {
      throw("No entry in Config/Databases.pm hash for $name");
    }
  }else{
    $db = $hash->{$name};
  }
  throw("Unable to find a db with name ".$name) if not $db;
  return $db;
}



1;
