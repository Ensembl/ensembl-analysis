
=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild - 

=head1 SYNOPSIS


=head1 DESCRIPTION

This is to be a baseclass for the genebuild code to have various
object like methods which genebuild modules want. GeneBuild code
should ideally inherit from this

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use strict;
use warnings;

use Data::Dumper;

use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Config::Databases; 
use Bio::EnsEMBL::Analysis::Config::General qw(BIN_DIR ANALYSIS_WORK_DIR);
use Bio::EnsEMBL::Analysis::RunnableDB;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

# first of all we just need methods to get databases from the
# database config


sub database_hash {
  my ($self, $name, $db) = @_;
  if(!$self->{'db_hash'}){
    $self->{'db_hash'} = {};
  }
  if($name && $db){
    $self->{'db_hash'}->{$name} = $db;
  }
  return $self->{'db_hash'};
}


=head2 select_random_db 

  Arg       : String - HashKey pointing to an entry in Databases.pm - 
              Either this key points to a key in the hash %DATABASES or in %DISTRIBUTED_DBS 

  Fuction   : The fuction reads the 2 hashes %DATABASES and %DISTRIBUTED_DBS which are exported by 
              Databases.pm. If the Argument is found in  %DATABASES, the name is returned. 
              if the Argument points to an entry in %DISTRIBUTED_DBS and %DISTRIBUTED_DBS{$arg} is 
              an array reference, an element is randomly picked out of this array and returned. 

              This function is basically used to spread the load over different db servers randomly. 
              with the $DISTRIBUTED_DBS array in Databases.pm 
 
  Returntype: String 

=cut


sub select_random_db {  
  my (  $name ) = @_; 
 
  my $tmp;    
  if (exists $DATABASES->{$name} && ref($DATABASES->{$name}) =~m/AREF/  ) { 
    $tmp = $DATABASES ; 
  } elsif ( exists  $DISTRIBUTED_DBS->{$name} && ref($DISTRIBUTED_DBS->{$name}) =~m/ARRAY/   ) {  
    $tmp = $DISTRIBUTED_DBS; 
  } 
  if ( defined $tmp ) { 
     my @array = @{ $tmp->{$name} };   
     my $randomIndex = rand(@array);  
     $name = $array[$randomIndex]; 
     print "Random database selected : $name \n"; 
  }
  return $name ; 
}

=head2 get_dbadaptor

  Arg [0]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [1]   : String - key of database hash
  Arg [2]   : return a non-standard adaptor [ valie values : 'pipeline' 'compara' 'functgenomics' or undef ] 
  Arg [3]   : flag to attch dna_db nor not 

  Function  : Returns a Bio::EnsEMBL::DBSQL::DBAdaptor for a given hash key.
              or a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor if requested
              Requires proper configuration of
              Bio::EnsEMBL::Analysis::Config::Databases

  Returntype: Bio::EnsEMBL:DBSQL::DBAdaptor or Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
              Returns undef if an entry for the database is found but is empty.
  Exceptions: throw if key can't be found in Databases.pm

=cut



sub get_dbadaptor {
  my ( $self, $name, $non_standard_db_adaptor, $not_use_dna_database ) = @_;

  my ( $db, $try_to_attach_dna_db );
  my $hash = $self->database_hash;

  $name = select_random_db($name);  

  if ( !$hash->{$name} ) {   # if we don't already have an entry for this ...
    if ( exists $DATABASES->{$name} ) { 

      my $constructor_args = $DATABASES->{$name};

      if ( scalar( keys( %{$constructor_args} ) ) == 0 ) {
        # The entry is empty.  Warn about this, but don't throw.
        # Return undef.
        warning(
             sprintf( "Empty entry for database '%s' in Databases.pm\n",
                      $name ) );
        return undef;
      }

      # check if we got all arguments
      foreach my $arg (qw ( -user -port -host -dbname)) {
        unless ( $$constructor_args{$arg} ) {
          throw(   "Database-connection-details not properly configured : "
                 . "Argument : $arg missing in Databases.pm for $name \n" );
        }
      } 


      if ( defined $non_standard_db_adaptor ) { # value of 
        if (    $non_standard_db_adaptor =~ m/1/ || $non_standard_db_adaptor eq "pipeline" ) {
          require Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
          $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( %$constructor_args );

        } elsif ( $non_standard_db_adaptor =~ m/compara/ ) {
          require Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
          unless ( $$constructor_args{'-species'} ) {
            throw("need species !\n");
          }
          $db =
            Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( %$constructor_args );

        } elsif (    $non_standard_db_adaptor =~ m/functgenomics/i
                  || $non_standard_db_adaptor =~ m/funcgen/i )
        {
          require Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
          # funcgen adaptor needs species
          if ( $$constructor_args{'-species'} ) {
            $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
                                                        %$constructor_args, );
          } else {
            throw( "if you require a Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor "
                   . "you need to provide a -speices flag.\n" );
          }
        }
        else {
          throw( "No matching non-standard adaptor could be found. If you want ".
                 "a standard adaptor you should pass in undef.");
        }

      } else {
        $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$constructor_args);
        # it's a core db so try to attach dna_db
        $try_to_attach_dna_db = 1;
      }

      # this bit is attaching a dna_db .
      if (    $name ne $DNA_DBNAME
           && defined($DNA_DBNAME)
           && $try_to_attach_dna_db )
      {
        if ($not_use_dna_database) {
          print "\nNot attaching a DNA_DB to $name\n";
          # if two different species are considered, the
          # not_use_dna_database is set to 1 to avoid adding the second species to
          # the first one

        } else {
          # there's a little danger if we have multiple diffent
          # species in our "Databases.pm" file. We need to avoid that the wrong
          # dna db is attached, ie a mouse core with a human dna db.

          print "\nAttaching DNA_DB $DNA_DBNAME to $name...\n"; 
          if ( length ( $DNA_DBNAME ) == 0 ) {  
            throw("you're using an empty string as dna_dbname in your Databases.pm config"); 
          } 
          my $dnadb = $self->get_dbadaptor($DNA_DBNAME);

          # try to get default asm+ species name for OTHER db - does not work
          # for comapra database
          my $core_db_asm = $db->get_CoordSystemAdaptor->get_default_version();
          my $core_db_species =
            $db->get_MetaContainer->get_common_name();

          # get the same for dna-db
          my $dna_db_asm = $dnadb->get_CoordSystemAdaptor->get_default_version();
          my $dna_db_species =
            $dnadb->get_MetaContainer()->get_common_name();

          my $dna_db_and_core_db_are_compatible = 1;

          unless ( $core_db_asm eq $dna_db_asm ) {
            warning( "you try to add  a DNA_DB with assembly $dna_db_asm to "
              . "a core/cdna/otherfeatures DB with assembly $core_db_asm ...\n\t"
              . "that's incompatbile. I will not add dna_database "
              . $dnadb->dbname . " to core " . $db->dbname . "\n" );

            $dna_db_and_core_db_are_compatible = 0;
          }

          unless ( $core_db_species eq $dna_db_species ) {
            warning( "you try to add a DNA_DB with species ".$dna_db_species." to "
                . "a core database with species: '" .$core_db_species . "' - this does not work. \n"
                . "Check that you are using the correct DNA_DB and that the species.common_name values in the meta tables match\n"
            );
            $dna_db_and_core_db_are_compatible = 0;
          }
          if ($dna_db_and_core_db_are_compatible) {
            $db->dnadb($dnadb);
            print "\nAttaching DNA_DB "
              . $dnadb->dbc->dbname . " to "
              . $db->dbc->dbname . "\n";
          }
        } ## end else [ if ($not_use_dna_database)
      } else {
        if ( $name eq $DNA_DBNAME ) {
          print "\nNot attaching DNA_DB to $name ...\n" ; 
        } else {
          warning("You haven't defined a DNA_DBNAME in Config/Databases.pm ");
        }
      }
    } else {
      throw("No entry in Config/Databases.pm hash for $name");
    }
    $self->database_hash( $name, $db );
  } else {
    $db = $hash->{$name};
  }
  return $db;
} ## end sub get_dbadaptor



# parses a config section like this : 
#   VALIDATION_DBS => {
#                      HUMAN_DB => ['protein_coding','processed_transcript'],
#                      MOUSE_DB => ['ncrna','cdna'], 
#                     }                                                                                                 },
# -> gets DBAdaptor for HUAM_DB ( connection details defined in Databaess.pm ) and then 
# gets gene_adaptor to fetch the genes of the biotypes specified.  Only genes on input_d are fetched 
# 

=head2 get_genes_of_biotypes_by_db_hash_ref

  Arg [0]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [1]   : Hashref. where the keys represent keys in the Databases-Hash, 
              with a reference to an array where the values represent gene-biotypes 
  Function  : Loops through the keys of the hash, builds a DB-connection to each DB specified by the key
              and fetches all genes of the biotypes specified 

  Returntype: Returns an Array-reference to Bio::EnsEMBL::Gene objects which are located on the slice

=cut


sub get_genes_of_biotypes_by_db_hash_ref { 
  my ($self, $href) = @_;

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
         print scalar(@genes_to_fetch) . " genes fetched in total\n" ; 
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

=head2 main_reference_db

    Function  : Return the value of MAIN_REFERENCE_DB from the Databases.pm file
                It point to the DB that should be THE reference
    Returntype: Returns a string

=cut

sub main_reference_db {
    my $self = shift;
    return $MAIN_REFERENCE_DB;
}

=head2 vital_tables

   Function  : Return the values of VITAL_TABLES from the Databases.pm file
               It is all the vital tables from the schema that should be
               used to create a new database with a reference
   Returntype: Returns a reference on an array of strings 

=cut

sub vital_tables {
    my $self = shift;
    return $VITAL_TABLES;
}


=head2 create_new_database

    Arg [0]   : the hash key of the reference database
    Arg [1]   : the hash key of the new database
    Function  : It create a new database using a reference, usually the MAIN_REFERENCE_DB,
                it import the values in the vital tables (VITAL_TABLES found in Databases.pm)
    ReturnType: None
    Exception : throw if the new database already exists

=cut

sub create_new_database {
    my ($self, $old_db_key, $new_db_key) = @_;

    if (!db_exists($new_db_key)) {
        my $template_dump = $self->dump_template_db($old_db_key);
        my %h_new_db = %{$$DATABASES{$new_db_key}};

        print "creating new database ".$h_new_db{-dbname}." \@ ".$h_new_db{-host}." ".$template_dump."\n"; 
        my $create_cmd = 'mysql -h '.$h_new_db{-host}.' -u '.$h_new_db{-user}.' -p'.$h_new_db{-pass}.' -P'.$h_new_db{-port}.' -e" create database  '.$h_new_db{-dbname}.'"';
        system($create_cmd); 
        $create_cmd = 'mysql -h '.$h_new_db{-host}.' -u '.$h_new_db{-user}.' -p'.$h_new_db{-pass}.' -P'.$h_new_db{-port}.' -D '.$h_new_db{-dbname}.' < '.$template_dump; 
        system($create_cmd);  
        return $template_dump;
    }
    else {
        throw($new_db_key." already exists!");
    }
}


=head2 dump_template_db

    Arg [0]   : the hash key of the database to dump
    Function  : Dump the template database given in parameters. A file is created,
                containing all the tables without data and all the vital tables
                (VITAL_TABLES from Databases.pm) and their data are append after
    Returntype: Returns the path to the file
    Exception : throw if the dumping of the table fails

=cut

sub dump_template_db {
   my ($self, $db_key) = @_ ;  

#    throw("Database parameters are not properly configured for $db_key") unless (check_db_params($db_key));
    my $dbn = $self->get_dbadaptor($db_key);
    my $dbname = $dbn->dbname();
    my $host = $dbn->host();
    my $port = $dbn->port();
    my $file = $db_key."_".$dbname."_".$host."_".$port."_".$$.".dump.sql"; 
    my $dir = $ANALYSIS_WORK_DIR;  

    if ( defined $ENV{BASE} ) {  
        if ( -e "$ENV{BASE}/sql" ) {  
          $dir = "$ENV{BASE}/sql";
        } 
    }   

    my $file_name = "$dir/$file" ; 
    my $cmd = -x $BIN_DIR.'/mysqldump' ? $BIN_DIR."/mysqldump" : "mysqldump";

    $cmd .= ' -h '.$dbn->host.' -u '.$dbn->username.' -P'.$dbn->port;  

    if ( $dbn->password ne '' ) { 
        if ( length($dbn->password) > 0) {  
         $cmd .= ' -p'.$dbn->password; 
        }
    } 
    # print dump create table statements first ...  
    my $cmd1 = "$cmd --no-data ".$dbn->dbname." > $file_name "; 
    info("CMD 1 : $cmd1\n");
    print "dumping data to $file_name\n" ; 
    system($cmd1);
    throw("Could not dump empty tables from ".$dbn->dbname." in file $file_name") if $@;

    # now dump some tables with data ...
    $cmd .= " --add-drop-table ".$dbn->dbname." ".join(' ', @{$self->vital_tables}) ." >> $file_name ";  
    info("CMD 2 : $cmd\n");
    system($cmd); 
    throw("Could not dump vital tables from ".$dbn->dbname." in file $file_name") if $@;
    return $file_name ; 
}

#sub check_db_params {
#  my ($db_key) = @_ ;
#  my %href = %{$$DATABASES{$db_key}}; 
#
#  for my $arg (qw( -dbname -user -host  -port ) ) { 
#    if ( !defined $href{$arg}  || length($href{$arg}) == 0 ) {  
#        if ( $arg !~ m/-dbname/) {  
#         warning(" Database parameters are not properly configured for $db_key - argument $arg s missing, SKIPPING creation of $db_key \n");   
#        }
#        return 0 ; 
#    } 
#  }  
#  return 1;
#} 


=head2 db_exists

    Arg [0]   : the hash key of a database
    Function  : Check if the database already exists
    Returntype: Returns 1 if the database exists, otherwise undef

=cut

sub db_exists {  
  my $db_key = shift;
  $db_key = shift if ($db_key =~ /Bio::EnsEMBL/);
  my %href = %{$$DATABASES{$db_key}};  
  $href{'-dbname'} = "mysql" ; 
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%href); 
  my $sth=$dba->prepare('show databases like "'.$$DATABASES{$db_key}{-dbname}.'"' )  ; 
  $sth->execute;  
  my $result = $sth->fetchall_arrayref(); 
  for my $arg ( @$result ) {    
     my $db_name_on_server = $$arg[0];   
     if ( $db_name_on_server eq $$DATABASES{$db_key}{-dbname} ) { 
       return 1;
     } 
  }   
  return undef ;
} 

1;
