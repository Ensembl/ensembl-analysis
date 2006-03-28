#This is to be a baseclass for the genebuild code to have various
#object like methods which genebuild modules want. GeneBuild code 
#should ideally inherit from this


package Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info); 

use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases qw(DATABASES);
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases;


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


sub get_dbadaptor{
  my ($self, $name) = @_;
  my $hash = $self->database_hash;
  my $db;
  if(!$hash->{$name}){
    my $constructor_args = $DATABASES->{$name};
    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                              %$constructor_args,
                                             );
    $self->database_hash($name, $db);
  }else{
    $db = $hash->{$name};
  }
  throw("Unable to find a db with name ".$name);
  return $db;
}



1;
