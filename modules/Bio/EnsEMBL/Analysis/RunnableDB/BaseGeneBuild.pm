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
  my ($self, $name, $use_pipeline_adaptor) = @_;
  my $hash = $self->database_hash;
  my $db;
  if(!$hash->{$name}){
    if (exists $DATABASES->{$name}) {
      my $constructor_args = $DATABASES->{$name}; 

      foreach my $arg ( qw ( -user -port -host -dbname) ) {  
        unless ( $$constructor_args{$arg}){ 
          throw ("Database-connection-details not properly configured : Arguemnt : $arg missing in Databases.pm\n") ; 
        }
      }
      if ( $use_pipeline_adaptor ) {
	$db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
							    %$constructor_args,
							   );
      } else {
	$db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
						  %$constructor_args,
						 );
      }
      if($name ne $DNA_DBNAME ){
        if (length($DNA_DBNAME) ne 0 ){
          my $dnadb = $self->get_dbadaptor($DNA_DBNAME);
          $db->dnadb($dnadb);
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
