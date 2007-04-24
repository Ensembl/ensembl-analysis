
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
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);

@EXPORT = qw( shuffle parse_config create_file_name write_seqfile merge_config_details
              get_input_arg get_db_adaptor_by_string );






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
        if ( defined ${$file{$db_class}}{OUTDB} && length(${$file{$db_class}}{OUTDB}{'-dbname'}) > 0  ) {
        # don't process undefiend OUT-DB's and  don't process defiened OUT-DB's which have no name
           #print "-dbname "  .${$file{$db_class}}{OUTDB}{'-dbname'}. "\n\n\n" ;

          $result{$db_class}{db} = ${$file{$db_class}}{OUTDB} ;

        }else {
         next ;
        }
      }

      # process /Conf/GeneBuild/Databases.pm 

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
  my ($obj, $var_hash, $label) = @_;

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
    throw("Your logic_name ".$uc_logic." doesn't appear in your config file hash ".
          $var_hash);
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
  Function  : Returns a Bio::EnsEMBL::DBSQL::DBAdaptor for a given string.
              Requires proper configuration of 
              Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases 
 
  Returntype: Bio::EnsEMBL::Analysis::DBSQL::DBAdaptor 
  Exceptions: throw if string can't be found in Databases.pm 

=cut

sub get_db_adaptor_by_string {
   my ($string) = @_ ;


   require "Bio/EnsEMBL/Analysis/Config/GeneBuild/Databases.pm" ;
   no strict ;
   Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases->import("DATABASES");
   Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases->import("DNA_DBNAME");

   my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ ${$DATABASES}{$string} } ) ;
   my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ ${$DATABASES}{$DNA_DBNAME} } ) ;

   if($string ne $DNA_DBNAME ){
     if (length($DNA_DBNAME) ne 0 ){
        $db->dnadb($dnadb);
     }else{
        warning("You haven't defined a DNA_DBNAME in Config/Databases.pm");
     }
   }
  use strict ;
  throw("No entry in Config/GeneBuild/Databases.pm hash for $string") unless $db ;
  return $db;
}


1;
