
=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Utilities

- base class which exports utility methods which don't take Bio::XX objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::Utilities qw(shuffle);

  or

  use Bio::EnsEMBL::Analysis::Tools:Utilities

  to get all methods

=head1 DESCRIPTION

This is a class which exports Utility methods for genebuilding and
other gene manupulation purposes.

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::Utilities;

use strict;
use warnings;
use Exporter;
use Bio::EnsEMBL::Analysis::Tools::Stashes qw( package_stash ) ; # needed for read_config()
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);

@EXPORT = qw( shuffle
              parse_config
              parse_config_mini
              create_file_name
              write_seqfile
              merge_config_details
              get_input_arg
              get_db_adaptor_by_string read_config
              import_var
              get_database_connection_parameters_by_string ) ;






=head2 merge_config_details

  Arg [0]   : Array of Hashreferences
  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Function  : This func. merges the Configurations out differnt configuration-files into one Hash
  Returntype: Hashref.
  Exceptions: throws as this method should be implemented by any child
  Example   : merge_database_configs ($DATABASES, $EXONERATE2GENES, $TRANSCRIPT_COALESCER) ;

=cut


sub merge_config_details {
  my ($self,  @config_hashes )= @_ ;

  my %result ;

  # loop through all hrefs which are passed as input

  foreach my $config_file ( @config_hashes ) { 
    my %file = %$config_file ;
    foreach my $db_class ( keys %file ) {
      # process Exonerate2Genes.pm config (has section --> OUTDB)

      if ( exists ${$file{$db_class}}{OUTDB} ) {

       if ( ref(${$file{$db_class}}{OUTDB}) !~m/HASH/ && defined ${$file{$db_class}}{OUTDB}) { 
         # section in Exonerate2Genes is not a HREF which defines details for
         # database-connection - it's a hash-key pointing to Databases.pm

         my $href = get_database_connection_parameters_by_string(${$file{$db_class}}{OUTDB}) ;

         unless ( $href )  {
          print " $db_class  parameters are not defined in Databases.pm - skipping\n";
          next ;
         } else {
            #print "Used database : $$href{'-dbname'}\n" ;
            $result{$db_class}{db} = $href ;
         }
        }else {
          if ( defined ${$file{$db_class}}{OUTDB}
          && length(${$file{$db_class}}{OUTDB}{'-dbname'}) > 0  ) {
          # don't process undefiend OUT-DB's and
          # don't process defined OUT-DB's which have no name
            $result{$db_class}{db} = ${$file{$db_class}}{OUTDB} ;

          }else {
           next ;
          }
        }
      }

      # process /Conf/Databases.pm

      if (defined ( ${$file{$db_class}}{'-dbname'}) &&  length ( ${$file{$db_class}}{'-dbname'}) > 0 )  {
        # we process Databases.pm // parameteres for db-connection are ok
        $result{$db_class}{db} = \%{$file{$db_class}} ;

      } elsif (defined ( ${$file{$db_class}}{'-dbname'}) &&  length ( ${$file{$db_class}}{'-dbname'}) == 0  ) {
        next ;
      }

      # add / process data from other configs in format TranscriptCoalescer.pm
      # and attach data to main config hash

      for my $key (keys %{$file{$db_class}}) {
        $result{$db_class}{$key} = $file{$db_class}{$key};
      }
    }
  }
  return \%result ;
}


=head2 shuffle

  Arg [1]   : Reference to Array
  Function  : randomizes the order of an array
  Returntype: arrayref
  Exceptions: none
  Example   :

=cut

sub shuffle {
  my $tref = shift ;
  my $i = @$tref ;
  while ($i--) {
     my $j = int rand ($i+1);
     @$tref[$i,$j] = @$tref[$j,$i];
  }
  return $tref ;
}



=head2 get_input_arg

  Function  : waits for input from STDIN and returns '1' if input =~m/y/i
              and '0' if input matches /n/i.
  Returntype: 1 or 0
  Exceptions: none

=cut

sub get_input_arg {
  while (defined (my $line=<STDIN>)){
   chomp($line) ;
   if ( $line=~m/y/i){
      return 1 ;
   }elsif( $line =~m/n/i){
     return 0 ;
   }
   print "Wrong input - only answer 'y' or 'n'\n" ;
  }
}


sub parse_config{
  my ($obj, $var_hash, $label, $ignore_throw) = @_; 

  throw("Can't parse the ".$var_hash." hash for object ".$obj." if we are give no label") if(!$label); 

  my $DEFAULT_ENTRY_KEY = 'DEFAULT';
  if(!$var_hash || ref($var_hash) ne 'HASH'){
    my $err = "Must pass read_and_check_config a hashref with the config ".
      "in ";
    $err .= " not a ".$var_hash if($var_hash);
    $err .= " Utilities::read_and_and_check_config";
    throw($err);
  }

  my %check;
  foreach my $k (keys %$var_hash) {
    my $uc_key = uc($k);
    if (exists $check{$uc_key}) {
      throw("You have two entries in your config with the same name (ignoring case)\n");
    }
    $check{$uc_key} = $k;
  }
  # replace entries in config has with lower case versions.
  foreach my $k (keys %check) {
    my $old_k = $check{$k};
    my $entry = $var_hash->{$old_k};
    delete $var_hash->{$old_k};

    $var_hash->{$k} = $entry;
  }

  if (not exists($var_hash->{$DEFAULT_ENTRY_KEY})) {
    throw("You must define a $DEFAULT_ENTRY_KEY entry in your config");
  }

  my $default_entry = $var_hash->{$DEFAULT_ENTRY_KEY};
  # the following will fail if there are config variables that
  # do not have a corresponding method here
  foreach my $config_var (keys %{$default_entry}) {
    if ($obj->can($config_var)) {
      $obj->$config_var($default_entry->{$config_var});
    } else {
      throw("no method defined in Utilities for config variable '$config_var'");
    }
  }

  #########################################################
  # read values of config variables for this logic name into
  # instance variable, set by method
  #########################################################
  my $uc_logic = uc($label);
  if (exists $var_hash->{$uc_logic}) {
    # entry contains more specific values for the variables
    my $entry = $var_hash->{$uc_logic};

    foreach my $config_var (keys %{$entry}) {

      if ($obj->can($config_var)) {

        $obj->$config_var($entry->{$config_var});
      } else {
        throw("no method defined in Utilities for config variable '$config_var'");
      }
    }
  }else{ 
    if ( defined $ignore_throw && $ignore_throw== 1 ){ 
      warning("Your logic_name ".$uc_logic." doesn't appear in your config file hash - using default settings\n".  $var_hash);  
    }else{
      throw("Your logic_name ".$uc_logic." doesn't appear in your config file hash - using default settings\n".  $var_hash); 
    }
  }
}


sub parse_config_mini{
  my ($obj, $var_hash, ) = @_; 


  if(!$var_hash || ref($var_hash) ne 'HASH'){
    my $err = "Must pass read_and_check_config a hashref with the config ".
      "in ";
    $err .= " not a ".$var_hash if($var_hash);
    $err .= " Utilities::read_and_and_check_config";
    throw($err);
  }

  my %check;
  foreach my $k (keys %$var_hash) {
    my $uc_key = uc($k);
    if (exists $check{$uc_key}) {
      throw("You have two entries in your config with the same name (ignoring case)\n");
    }
    $check{$uc_key} = $k;
  }

  #########################################################
  # read values of config variables for this logic name into object methods which are defined in the class
  #########################################################

  foreach my $config_var (keys %{$var_hash}) { 
    if ($obj->can($config_var)) {
      $obj->$config_var($$var_hash{$config_var}); 
    } else {
      throw("no method defined in Utilities for config variable '$config_var'");
    }
  } 
}





=head2 create_file_name

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : string, stem of filename
  Arg [3]   : string, extension of filename
  Arg [4]   : directory file should live in
  Function  : create a filename containing the PID and a random number
  with the specified directory, stem and extension
  Returntype: string, filename
  Exceptions: throw if directory specifed doesnt exist
  Example   : my $queryfile = $self->create_filename('seq', 'fa');

=cut



sub create_file_name{
  my ($stem, $ext, $dir) = @_;
  if(!$dir){
    $dir = '/tmp';
  }
  $stem = '' if(!$stem);
  $ext = '' if(!$ext);
  throw($dir." doesn't exist SequenceUtils::create_filename")
    unless(-d $dir);
  my $num = int(rand(100000));
  my $file = $dir."/".$stem.".".$$.".".$num.".".$ext;
  while(-e $file){
    $num = int(rand(100000));
    $file = $dir."/".$stem.".".$$.".".$num.".".$ext;
  }
  return $file;
}



=head2 write_seq_file

  Arg [1]   : Bio::Seq
  Arg [2]   : string, filename
  Function  : This uses Bio::SeqIO to dump a sequence to a fasta file
  Returntype: string, filename
  Exceptions: throw if failed to write sequence
  Example   :

=cut

sub write_seqfile{
  my ($seq, $filename, $format) = @_;
  $format = 'fasta' if(!$format);
  my @seqs;
  if(ref($seq) eq "ARRAY"){
    @seqs = @$seq;
    throw("Seqs need to be Bio::PrimarySeqI object not a ".$seqs[0])
      unless($seqs[0]->isa('Bio::PrimarySeqI'));
  }else{
    throw("Need a Bio::PrimarySeqI object not a ".$seq)
      if(!$seq || !$seq->isa('Bio::PrimarySeqI'));
    @seqs = ($seq);
  }
  $filename = create_file_name('seq', 'fa', '/tmp')
    if(!$filename);
  my $seqout = Bio::SeqIO->new(
                               -file => ">".$filename,
                               -format => $format,
                              );
  foreach my $seq(@seqs){
    eval{
      $seqout->write_seq($seq);
    };
    if($@){
      throw("FAILED to write $seq to $filename SequenceUtils:write_seq_file $@");
    }
  }
  return $filename;
}



=head2 get_db_adaptor_by_string

  Arg [1]   : String
  Arg [2]   : verbose-flag
  Arg [3]   : return a pipeline db adaptor flag 
  Arg       : binary - -do_not_attach_dna_db 
  Function  : Returns a Bio::EnsEMBL::DBSQL::DBAdaptor for a given string.
              or a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor if requested
              Requires proper configuration of
              Bio::EnsEMBL::Analysis::Config::Databases

  Returntype: Bio::EnsEMBL:DBSQL::DBAdaptor or Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Exceptions: throw if string can't be found in Databases.pm

  Example : 

          get_db_adaptor_by_string("SOLEXA_DB" , 1, -do_not_attach_dna_db =>1  ) ;
=cut





sub get_db_adaptor_by_string {
   my ($string, $verbose, $use_pipeline_adaptor,@args) = @_ ;

   my ($no_dna_db) = rearrange( 'do_not_attach_dna_db', @args ); 

   #print "Fetching ".$string."\n";
   require "Bio/EnsEMBL/Analysis/Config/Databases.pm" ;
   no strict ;
   Bio::EnsEMBL::Analysis::Config::Databases->import("DATABASES");
   Bio::EnsEMBL::Analysis::Config::Databases->import("DNA_DBNAME");

   unless ( ${$DATABASES}{$string} ) {
     print "WARNING : Database parameters undefined for - skipping \n" ;
     return undef ;
   }

   if ( length(${$DATABASES}{$string}{'-dbname'}) == 0 ) {
     print "WARNING : You haven't defined a database-name in the Databases.pm config-file for $string\n" ;
     return undef ;
   }

   my $db;
   my $dnadb;
   if ( $use_pipeline_adaptor ) {
      require 'Bio/EnsEMBL/Pipeline/DBSQL/DBAdaptor.pm';
     $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor( %{ ${$DATABASES}{$string} } ) ;
   } else {
     $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ ${$DATABASES}{$string} } ) ;
   } 

   #print "Got ".$db."\n";
   if ( $verbose ) {
     my %tmp =  %{${$DATABASES}{$string}} ;
     print STDERR "Database : $tmp{'-dbname'} @ $tmp{'-host'} : $tmp{'-port'} AS $tmp{'-user'} - $tmp{'-pass'}\n" ;
   }

   if (! $no_dna_db ) { 
     if($string ne $DNA_DBNAME ){
       if (length($DNA_DBNAME) ne 0 ){
        my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ ${$DATABASES}{$DNA_DBNAME} } ) ;
        $db->dnadb($dnadb);
       }else{
        warning("You haven't defined a DNA_DBNAME in Config/Databases.pm");
      }
     }
   }
  use strict ;
  return $db;
}




=head2 get_database_connection_parameters_by_string

  Arg [1]   : String
  Function  : Returns a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor for a given string.
              Requires proper configuration of
              Bio::EnsEMBL::Analysis::Config::Databases

  Returntype: Hashref
  Exceptions: throw if string can't be found in Databases.pm

=cut

sub get_database_connection_parameters_by_string {
   my ($string) = @_ ;


   
   require "Bio/EnsEMBL/Analysis/Config/Databases.pm" ;
   no strict ;
   Bio::EnsEMBL::Analysis::Config::Databases->import("DATABASES");
   Bio::EnsEMBL::Analysis::Config::Databases->import("DNA_DBNAME");

   unless ( ${$DATABASES}{$string} ) {
     print "WARNING : Database parameters undefined - skipping \n" ; 
     print stack_trace_dump() ; 
     return undef ;
   }

   if ( length(${$DATABASES}{$string}{'-dbname'}) == 0 ) {
     print "You haven't defined a database-name in the Databases.pm config-file for $string\n" ;
     return undef ;
   }
   return ${$DATABASES}{$string} ;
}




=head2 read_config

  Arg [1]   : String
  Arg [2]   : optional Array-reference
  Function  : reads a configuration file in the ensembl perl module format
              in runtime and, accesses the variables ( either only the variables
              listend in $aref or all variables ) and returns a hash-reference to this.

  Returntype:
  Exceptions:
  Requires  : use Bio::EnsEMBL::Analysis::Tools::Stashes qw( package_stash ) ;

=cut

sub read_config {
   my ($module_name , $aref ) = @_ ;

   (my $module_path = $module_name )=~s/::/\//g;
   require "$module_path.pm" ;

   # get the names of the variables
   unless ($aref) {
     my ($config_href, $varname ) = @{package_stash("$module_name")};
     map { $module_name->import($_) } keys %$config_href ;

     return $config_href;
   }

   # import only variables specified in $aref
   no strict ;
   map { $module_name->import($_) } @$aref ;
   my %import ;
   map {$import{$_} = ${$_}} @$aref ;
   use strict ;
   return \%import;
}



=head2 import_var

  Arg [0]   : Array of Hashreferences
  Arg [1]   : optional Array-ref
  Function  : gets a hash-reference and adds them to the namespace of the module - they can
              be accesed in the package by using 'no strict;'
  Returntype: Hashref.
  Examples  : import_var ($href) ;
              import_var(read_config("Bio::EnsEMBL::Analysis::Config::Databases"));


=cut



sub import_var {
    my ($callpack) = caller(0); # Name of the calling package
   # my $pack = shift; # Need to move package off @_
    my $vars_to_import = shift ;

   # $vars_to_import = $pack unless $vars_to_import ;
    # Get list of variables supplied, or else all

    my @vars = @_ ? @_ : keys(%{$vars_to_import});

    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;

    foreach (@vars) {
        if (defined ${$vars_to_import}{ $_ }) {
            no strict 'refs';
            # Exporter does a similar job to the following
            # statement, but for function names, not
            # scalar variables:
            *{"${callpack}::$_"} = \${$vars_to_import}{ $_ };
        } else {
            die "Error: Config: $_ not known\n";
        }
    }
}
1;
