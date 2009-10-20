#This is to be a baseclass for the genebuild code to have various
#object like methods which genebuild modules want. GeneBuild code
#should ideally inherit from this


package Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
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
  my ($self, $name, $non_standard_db_adaptor , $not_use_dna_database) = @_; 

  my ($db, $try_to_attach_dna_db) ;  

  my $hash = $self->database_hash; 

  if(!$hash->{$name}){
    if (exists $DATABASES->{$name}) {
      my $constructor_args = $DATABASES->{$name};

      # check if we got all arguments 
      foreach my $arg ( qw ( -user -port -host -dbname) ) {
        unless ( $$constructor_args{$arg}){
          throw ("Database-connection-details not properly configured : Argument : $arg missing in Databases.pm for $name \n") ;
        }
      } 
      if ( defined $non_standard_db_adaptor) { 
      if ( $non_standard_db_adaptor =~m/1/ || $non_standard_db_adaptor eq "pipeline") { 
         require Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor; 
         $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( %$constructor_args,); 

      } elsif ( $non_standard_db_adaptor =~m/compara/ ) {
         require Bio::EnsEMBL::Compara::DBSQL::DBAdaptor; 
         unless ( $$constructor_args{'-species'}){  
             throw("need species !\n"); 
          } 
         $db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( %$constructor_args,); 

      } elsif ( $non_standard_db_adaptor =~m/functgenomics/i || $non_standard_db_adaptor =~m/funcgen/i) {
         require Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
         # funcgen adaptor needs species   
         if ( $$constructor_args{'-species'}) { 
            $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new( %$constructor_args,); 
         }else {  
           throw("if you require a Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor you need to provide a -speices flag \n" ) ; 
         } 
      }
      }  else {
         $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %$constructor_args); 
         # it's a core db so try to attach dna_db 
         $try_to_attach_dna_db = 1 ; 
     } 

     #  this bit is attaching a dna_db . 
       
     if($name ne $DNA_DBNAME  && length($DNA_DBNAME) ne 0 && $try_to_attach_dna_db ){
        if ( $not_use_dna_database ) {
          # if two different species are considered, the not_use_dna_database is set to 1 to avoid adding the second species to the first one
        } else { 
         # there's a little danger if we have multiple diffent species in our "Databases.pm" file. 
         # We need to avoid that the wrong dna db is attached, ie a mouse core with a human dna db. 

         my $dnadb = $self->get_dbadaptor($DNA_DBNAME);  

         # try to get default asm+ species name for OTHER db - does not work for comapra database 
         my $core_db_asm = $db->get_MetaContainer->get_default_assembly();    
         my $core_db_species = $db->get_MetaContainer->get_Species->binomial();  

        # get the same for dna-db 
         my $dna_db_asm = $dnadb->get_MetaContainer->get_default_assembly();      
         my $dna_db_species = $dnadb->get_MetaContainer()->get_Species->binomial(); 
   
         my $dna_db_and_core_db_are_compatible = 1 ; 
   
         unless ( $core_db_asm eq $dna_db_asm ) { 
           warning ("you try to add  a DNA_DB with assembly $dna_db_asm to a core/cdna/otherfeatures DB with assembly $core_db_asm ...\n\t".
                      "that's incompatbile. I will not add dna_database " . $dnadb->dbname . " to core " . $db->dbname . "\n"); 
           $dna_db_and_core_db_are_compatible = 0 ; 
         }    
   
         unless ( $core_db_species eq $dna_db_species ) { 
            warning ("you try to add  a DNA_DB with species $dna_db_species to a core database with speices : $core_db_species - this does not work\n\t".
                     "that's incompatbile. try to not use any DNA_DATABASE name in Analysis/Config/Databases.pm\n"); 
            $dna_db_and_core_db_are_compatible = 0 ; 
         }  
         if ( $dna_db_and_core_db_are_compatible ) {  
              $db->dnadb($dnadb);    
              print "adding dna_db " . $dnadb->dbname . " to " . $db->dbname . "\n" ;  
         } 
       } 
    }else{
          warning("You haven't defined a DNA_DBNAME in Config/Databases.pm ");
    }
   } else { 
    throw("No entry in Config/Databases.pm hash for $name");
   }
  $self->database_hash($name, $db); 
  } else{
    $db = $hash->{$name};
  }
  return $db;
}



# parses a config section like this : 
#   VALIDATION_DBS => {
#                      HUMAN_DB => ['protein_coding','processed_transcript'],
#                      MOUSE_DB => ['ncrna','cdna'], 
#                     }                                                                                                 },
# -> gets DBAdaptor for HUAM_DB ( connection details defined in Databaess.pm ) and then 
# gets gene_adaptor to fetch the genes of the biotypes specified.  Only genes on input_d are fetched 
# 

sub get_genes_of_biotypes_by_db_hash_ref { 
  my ($self,$href ) = @_;

  my %dbnames_2_biotypes = %$href ;  

  my @genes_to_fetch;  

  foreach my $db_hash_key ( keys %dbnames_2_biotypes )  {

    my @biotypes_to_fetch = @{$dbnames_2_biotypes{$db_hash_key}};  

    my $set_db = $self->get_dbadaptor($db_hash_key); 
    my $slice = $self->fetch_sequence($self->input_id, $set_db)  ;
    
    # implementation of fetch_all_biotypes ....  
    my $fetch_all_biotypes_flag ; 
    foreach my $biotype  ( @biotypes_to_fetch ) {   
      if ($biotype=~m/fetch_all_biotypes/ ) {    
        $fetch_all_biotypes_flag = 1 ; 
      }
    }  
    if ( $fetch_all_biotypes_flag ) {  
         print "fetching ALL biotypes for slice out of db $db_hash_key :\n" ; 
         my $genes = $slice->get_all_Genes(undef,undef,1) ; 
         push @genes_to_fetch, @$genes;  
         my %tmp ; 
         for ( @$genes ) {  
           $tmp{$_->biotype}++; 
         }  
         foreach ( keys %tmp ) {  
           print "found $_ $tmp{$_}\n" ; 
         } 
         print scalar(@genes_to_fetch) . " genees fetched \n" ; 
    } else { 
      foreach my $biotype  ( @biotypes_to_fetch ) {  
         my $genes = $slice->get_all_Genes_by_type($biotype,undef,1);
         if ( @$genes == 0 ) {
           warning("No genes of biotype $biotype found in $set_db\n");
         } 
         if ( $self->verbose ) { 
           print "$db_hash_key [ " .$set_db->dbname  . " ] Retrieved ".@$genes." of type ".$biotype."\n";
         }
         push @genes_to_fetch, @$genes;
      }  
    }  
  } 
  return \@genes_to_fetch;
}




1;
