# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

please send any questions to http://lists.ensembl.org/mailman/listinfo/dev

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::Utilities;

use strict;
use warnings;
use feature 'say';

use Exporter qw(import);
use File::Spec::Functions qw(catfile tmpdir);
use File::Which;
use File::Temp;
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Analysis::Tools::Stashes qw( package_stash ) ; # needed for read_config()
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our @EXPORT_OK = qw(
              shuffle
              parse_config
              parse_config_mini
              parse_config_value
              create_file_name
              write_seqfile
              merge_config_details
              get_evidence_set
              convert_prediction_transcripts_to_genes
              get_input_arg
              get_db_adaptor_by_string read_config
              import_var
              is_canonical_splice
              get_analysis_settings
              get_database_connection_parameters_by_string
              run_command
              send_email
              hrdb_get_dba
              convert_to_ucsc_name
              align_proteins
              align_proteins_with_alignment
              align_nucleotide_seqs
              map_cds_location
              locate_executable
              first_upper_case
              execute_with_wait
              execute_with_timer
              parse_timer
              is_slice_name
              get_database_from_registry
              get_biotype_groups
              get_feature_name
              create_production_directory
              );






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



=head2 _get_evidence_set ($logic_name_or_biotype)

  Name     : get_evidence_set( $logic_name_or_biotype )
  Arg      : String 
  Func     : returns the name of the evidence_set of a genee / PredictionTranscript 
  Returnval: String describing evidence_set_name

=cut

sub get_evidence_set {
  my ($self, $logic_name_or_biotype) = @_ ;

  my %ev_sets = %{ $self->{evidence_sets} } ;
  my $result_set_name ;
  for my $set_name (keys %ev_sets){
    my @logic_names = @{$ev_sets{$set_name}} ;
    for my $ln (@logic_names ) {
       if ($ln eq $logic_name_or_biotype) {
         $result_set_name = $set_name ;
      }
    }
  }
  return $result_set_name ;
}

=head2 convert_prediction_transcripts_to_genes 

  Arg [0]   : reference to an array of Bio::EnsEMBL::PredictionTranscript-objects
  Arg [1]   : String describing the logic_name  
  Arg [2]   : String describing name of the evidence-set 
  Function  : Creates Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended-Objects from a 
              set of Bio::EnsEMBL::PredictionTranscript-Objects and adds further information to these Transcripts.
              The PredictionExons are re-blessed to Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended-objects
  Returntype: Ref. to arrray of  Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended-objects
  Example   : 

=cut


sub convert_prediction_transcripts_to_genes {
  my ($pt,$logic_name_becomes_biotype,$ev_set_name ) = @_ ;
  my @new_genes ;
  for my $pt (@$pt) {
    # conversion 
    my $gene_from_pt = Bio::EnsEMBL::Gene->new(
                       -start => $pt->start ,
                       -end => $pt->end ,
                       -strand => $pt->strand ,
                       -slice =>$pt->slice ,
                       -biotype => $logic_name_becomes_biotype,
                       -analysis=>$pt->analysis,
                       ) ;

    my $new_tr = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended->new(
                    -BIOTYPE => $logic_name_becomes_biotype ,
                    -ANALYSIS => $pt->analysis ,
                 ) ;

    my @pt_exons  = @{$pt->get_all_Exons} ;

    for (my $i=0 ; $i<scalar(@pt_exons) ; $i++) {

      # converting Bio::EnsEMBL::PredictionExon into ExonExtened (ISA Bio::EnsEMBL::Exon)  
      my $pte =$pt_exons[$i] ;
      bless $pte,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended" ;
      $pte->biotype($logic_name_becomes_biotype) ;
      $pte->ev_set($ev_set_name) ;
      $pte->end_phase(0);
      $pte->phase(0);
      $pte->next_exon($pt_exons[$i+1]) ;
      $pte->prev_exon($pt_exons[$i-1]) ;
      $pte->transcript($new_tr) ;
      $pte->analysis($pt->analysis) ;
    } ;

    #
    # Extending the Bio::EnsEMBL::Transcript object by ev_set methods 
    #
    for (@pt_exons) {
      $new_tr->add_Exon($_);
    }

    $gene_from_pt->add_Transcript($new_tr) ;

    push @new_genes , $gene_from_pt ;
  }
  return \@new_genes ;
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

sub parse_config_value{
  my ($obj, $var_hash, $label, $values_to_get) = @_; 

  throw("Can't parse the ".$var_hash." hash for object ".$obj." if we are give no label") if(!$label); 

  my $DEFAULT_ENTRY_KEY = 'DEFAULT';
  if(!$var_hash || ref($var_hash) ne 'HASH'){
    my $err = "Must pass read_and_check_config a hashref with the config ".
      "in ";
    $err .= " not a ".$var_hash if($var_hash);
    $err .= " Utilities::read_and_and_check_config_value";
    throw($err);
  }

  if (not exists($var_hash->{$DEFAULT_ENTRY_KEY})) {
    throw("You must define a $DEFAULT_ENTRY_KEY entry in your config");
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

  my $default_entry = $var_hash->{$DEFAULT_ENTRY_KEY};
  # the following will fail if there are config variables that
  # do not have a corresponding method here
  foreach my $config_var (@$values_to_get) {
    throw("$config_var does not exist in your config file for $label\n") unless (exists $default_entry->{$config_var});
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
      throw("Your logic_name ".$uc_logic." doesn't appear in your config file hash - using default settings\n".  $var_hash); 
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

  Arg [1]   : string, stem of filename
  Arg [2]   : string, extension of filename
  Arg [3]   : directory file should live in
  Function  : create a filename using File::Temp
  with the specified directory, stem and extension
  Returntype: File::Temp object
  Exceptions: throw if directory specifed doesnt exist
  Example   : my $queryfile = create_file_name('seq', 'fa');

=cut



sub create_file_name{
  my ($stem, $ext, $dir, $no_clean) = @_;

  my $random = 'XXXXX';
  my %params = (DIR => tmpdir);
  if ($dir) {
    if (-d $dir) {
      $params{DIR} = $dir;
    }
    else {
      throw(__PACKAGE__."::create_file_name: $dir doesn't exist");
    }
  }
  $params{TEMPLATE} = $stem.'_'.$random if ($stem);
  $params{SUFFIX} = '.'.$ext if ($ext);
  if($no_clean) {
    $params{UNLINK} = 0;
  }
  my $fh = File::Temp->new(%params);
  return $fh;
}



=head2 write_seqfile

  Arg [1]   : Bio::Seq
  Arg [2]   : string, filename
  Function  : This uses Bio::SeqIO to dump a sequence to a fasta file
  Returntype: string, filename
  Exceptions: throw if failed to write sequence
  Example   :

=cut

sub write_seqfile{
  my ($seq, $filename, $format, $no_clean) = @_;

  my %params = ( -format => $format || 'fasta');
  my $seqs;
  if(ref($seq) eq "ARRAY"){
    throw("Seqs need to be Bio::PrimarySeqI object not a ".$seq->[0])
      unless($seq->[0]->isa('Bio::PrimarySeqI'));
    $seqs = $seq
  }else{
    throw("Need a Bio::PrimarySeqI object not a ".$seq)
      if(!$seq || !$seq->isa('Bio::PrimarySeqI'));
    $seqs = [$seq];
  }
  if($filename) {
    if (ref($filename) eq 'File::Temp') {
      $params{-fh} = $filename;
    }
    else {
      $params{-file} = $filename;
    }
  }
  else {
    $params{-fh} = create_file_name('seq', 'fa', undef, $no_clean);
    $filename = $params{-fh};
  }

  my $seqout = Bio::SeqIO->new(%params);
  foreach my $seq(@$seqs){
    eval{
      throw("Problem writing to $filename") unless ($seqout->write_seq($seq));
    };
    if($@){
      throw("FAILED to write $seq to $filename ".__PACKAGE__.":write_seqfile $@");
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

# this method copied from sw4's ensembl-analysis/modules/Bio/EnsEMBL/Analysis/RunnableDB/ExonerateSolexaTranscript.pm
sub is_canonical_splice {
  my ($intron, $slice_adaptor, $slice) = @_;

  my $prev_Exon = $intron->prev_Exon;
  my $next_Exon = $intron->next_Exon;

  my $canonical;
  my $donor;
  my $acceptor;

  if ( $prev_Exon && $next_Exon) {
    my $donor_splice;
    my $acceptor_splice;
    if ($prev_Exon->strand  == 1 ) {
      # we are working on the forward strand
      $donor_splice = $slice_adaptor->fetch_by_region('toplevel',
                                                        $slice->seq_region_name,
                                                        $prev_Exon->end+1,
                                                        $prev_Exon->end+2,
                                                        $prev_Exon->strand
                                                       );
      $acceptor_splice = $slice_adaptor->fetch_by_region('toplevel',
                                                         $slice->seq_region_name,
                                                         $next_Exon->start-2,
                                                         $next_Exon->start-1,
                                                         $prev_Exon->strand
                                                        );
    } else {
      # we are working on the reverse strand
      $donor_splice = $slice_adaptor->fetch_by_region('toplevel',
                                                        $slice->seq_region_name,
                                                        $prev_Exon->start-2,
                                                        $prev_Exon->start-1,
                                                        $prev_Exon->strand
                                                       );
      $acceptor_splice = $slice_adaptor->fetch_by_region('toplevel',
                                                         $slice->seq_region_name,
                                                         $next_Exon->end+1,
                                                         $next_Exon->end+2,
                                                         $prev_Exon->strand
                                                        );
    }
    if ( $donor_splice->seq eq 'NN' && $acceptor_splice->seq eq 'NN' ) {
      warn("Cannot find dna sequence for prev_Exon " . $prev_Exon->stable_id . " next_Exon ".$next_Exon->stable_id.
           " this is used in detetcting non canonical splices\n");
    } else {

      #print "Splice type " . $acceptor_splice->seq ."- ".  $donor_splice->seq ."\n";
      #is it GTAG?
      $donor = $donor_splice->seq;
      $acceptor = $acceptor_splice->seq;
      if ( $acceptor_splice->seq eq 'AG' && ($donor_splice->seq eq 'GT' || $donor_splice->seq eq 'GC')) {
        # these combinations are canonical according to the HAVANA
        # annotation guidelines published on 18 June 2010
        $canonical = 1;
      } elsif ($donor_splice->seq eq 'AT' && $acceptor_splice->seq eq 'AC') {
        # this rare splcie site is canonical according to the HAVANA
        # annotation guidelines published on 18 June 2010
        $canonical = 1;
        print STDERR "Found AT-AC\n";
      } else {
        $canonical = 0;
      }
    }

  } else {
    throw("Cannot find previous or next exon");
  }
  return ($canonical, $donor, $acceptor);
}

=head2 run_command

  Arg [0]   : string containing the command to run
  Arg [1]   : optional description of the command to run that will be printed to the standard output if defined
  Arg [2]   : optional expected command numeric return value. If it does not match the command numeric return value, a exception will be thrown.
  Function  : It runs a command and returns its return value including an optional checking of the returned value.
  Returntype: string
  Examples  : run_command("wc -l stable_ids.txt","Counting stable IDs...",19831);

=cut

sub run_command {
  my ($command,$name,$expected_result) = @_;

  print($name) if ($name);
  print("\nRunning command:\n$command\n");

  my $result = `$command`;
  if ($?) {
    throw("Command FAILED: `$command`");
  }

  if (defined($expected_result)) {
    if (int($result) != $expected_result) {
      throw("Command: $command\nResult: $result\nExpected result: $expected_result\n");
    }
    else {
      print ("\nResult and expected result match: $expected_result\n");
    }
  }
  return $result;
}

=head2 send_email

  Arg [0]   : email to
  Arg [1]   : email from
  Arg [2]   : email subject
  Arg [3]   : email body
  
  Function  : It sends an email by using the 'sendmail' command.
  Returntype: n/a
  Examples  : send_email("pollo@granja.es","hen@farm.co.uk","Farm issues","My body is mine.");
  
=cut

sub send_email {
  my ($to,$from,$subject,$body) = @_;
  
  open(my $sendmail_fh, '|-', "sendmail '$to'");
  print $sendmail_fh "Subject: $subject\n";
  print $sendmail_fh "From: $from\n";
  print $sendmail_fh "\n";
  print $sendmail_fh "$body\n";
  close $sendmail_fh;
}

=head2 align_proteins

  Arg [0]   : source protein sequence
  Arg [1]   : target protein sequence

  Function  : It aligns the source protein sequence to the target protein sequence to
              calculate the coverage and the percent identity of the source against the target.
  Returntype: List containing (coverage,percent_identity) i.e. (82.7%,91.22%)
  Examples  : align_proteins("ADCDA","ADCTM");

=cut

sub align_proteins {
  my ($source_protein_seq,$target_protein_seq) = @_;

  my (undef,undef,$coverage,$percent_id) = align_proteins_with_alignment($source_protein_seq,$target_protein_seq);

  return ($coverage,$percent_id);
}

=head2 align_proteins_with_alignment

  Arg [0]   : source protein sequence
  Arg [1]   : target protein sequence

  Function  : It aligns the source protein sequence to the target protein sequence to
              calculate the coverage and the percent identity of the source against the target.
  Returntype: List containing (aligned_source_protein_seq,aligned_target_protein_seq,coverage,percent_identity) i.e. (82.7%,91.22%)
  Examples  : align_proteins_with_alignment("ADCDA","ADCTM");

=cut

sub align_proteins_with_alignment {
  my ($source_protein_seq,$target_protein_seq) = @_;

  my $align_input_file = "/tmp/align_".$$.".fa";
  my $align_output_file = "/tmp/align_".$$.".aln";

  open(INPUT,">".$align_input_file);
  say INPUT ">query";
  say INPUT $source_protein_seq;
  say INPUT ">target";
  say INPUT $target_protein_seq;
  close INPUT;

  my $align_program_path = 'mafft';

  my $cmd = $align_program_path." --amino ".$align_input_file." > ".$align_output_file;
  my $result = system($cmd);

  if ($result) {
    throw("Got a non-zero exit code from alignment. Command line used:\n".$cmd);
  }

#  my $align_program_path = 'muscle';

#  my $cmd = $align_program_path." -in ".$align_input_file." -out ".$align_output_file;
#  my $result = system($cmd);

#  if ($result) {
#    throw("Got a non-zero exit code from alignment. Command line used:\n".$cmd);
#  }

  my $file = "";
  open(ALIGN,$align_output_file);
  while (<ALIGN>) {
    $file .= $_;
  }
  close ALIGN;

  if ($file !~ /\>.+\n(([^>]+\n)+)\>.+\n(([^>]+\n)+)/) {
    warning("Could not parse the alignment file for the alignment sequences. Alignment file: ".$align_output_file);
    return (undef,undef,0,0);
  }

  my $aligned_source_protein_seq = $1;
  my $aligned_target_protein_seq = $3;

  $aligned_source_protein_seq =~ s/\n//g;
  $aligned_target_protein_seq =~ s/\n//g;

  `rm $align_input_file`;
  `rm $align_output_file`;

  # Work out coverage
  my $coverage;
  my $temp = $aligned_target_protein_seq;
  my $projected_gap_count = $temp =~ s/\-//g;
  my $ungapped_source_protein_seq = $aligned_source_protein_seq;
  $ungapped_source_protein_seq  =~ s/\-//g;

  if (length($ungapped_source_protein_seq) == 0) {
    $coverage = 0;
  } else {
    $coverage = 100-(($projected_gap_count/length($ungapped_source_protein_seq))*100);
  }

  # Work out percent identity
  my $match_count = 0;
  my $aligned_positions = 0;
  for (my $j = 0; $j < length($aligned_source_protein_seq); $j++) {
    my $char_query = substr($aligned_source_protein_seq,$j,1);
    my $char_target = substr($aligned_target_protein_seq,$j,1);
    if ($char_query eq '-' || $char_target  eq '-') {
      next;
    }
    if ($char_query eq $char_target) {
      $match_count++;
    }
    $aligned_positions++;
  }

  if ($aligned_positions <= 0) {
    throw("Pairwise alignment between the query protein sequence and the target protein sequence shows zero aligned positions. Something has gone wrong.");
  }
  my $percent_id = ($match_count/$aligned_positions)*100;
  $coverage = sprintf "%.2f", $coverage;
  $percent_id = sprintf "%.2f", $percent_id;
  return ($aligned_source_protein_seq,$aligned_target_protein_seq,$coverage,$percent_id);
}



=head2 align_nucleotide_seqs

  Arg [0]   : source sequence
  Arg [1]   : target sequence

  Function  : It aligns the source sequence to the target sequence to
              calculate the coverage and the percent identity of the source against the target.
  Returntype: List containing (coverage,percent_identity) i.e. (82.7%,91.22%)
  Examples  : align_nucleotide_seqs("ATTTA","ATCTA");

=cut

sub align_nucleotide_seqs {
  my ($source_protein_seq,$target_protein_seq) = @_;

  my $align_input_file = "/tmp/align_".$$.".fa";
  my $align_output_file = "/tmp/align_".$$.".aln";

  open(INPUT,">".$align_input_file);
  say INPUT ">query";
  say INPUT $source_protein_seq;
  say INPUT ">target";
  say INPUT $target_protein_seq;
  close INPUT;

  my $align_program_path = 'mafft';

  my $cmd = $align_program_path." --nuc ".$align_input_file." > ".$align_output_file;
  my $result = system($cmd);

  if ($result) {
    throw("Got a non-zero exit code from alignment. Command line used:\n".$cmd);
  }

  my $file = "";
  open(ALIGN,$align_output_file);
  while (<ALIGN>) {
    $file .= $_;
  }
  close ALIGN;

  if ($file !~ /\>.+\n(([^>]+\n)+)\>.+\n(([^>]+\n)+)/) {
    warning("Could not parse the alignment file for the alignment sequences. Alignment file: ".$align_output_file);
    return (undef,undef,0,0);
  }

  my $aligned_source_seq = $1;
  my $aligned_target_seq = $3;

  $aligned_source_seq =~ s/\n//g;
  $aligned_target_seq =~ s/\n//g;

  `rm $align_input_file`;
  `rm $align_output_file`;

  # Work out coverage
  my $coverage;
  my $temp = $aligned_target_seq;
  my $projected_gap_count = $temp =~ s/\-//g;
  my $ungapped_source_seq = $aligned_source_seq;
  $ungapped_source_seq  =~ s/\-//g;

  if (length($ungapped_source_seq) == 0) {
    $coverage = 0;
  } else {
    $coverage = 100-(($projected_gap_count/length($ungapped_source_seq))*100);
  }

  # Work out percent identity
  my $match_count = 0;
  my $aligned_positions = 0;
  for (my $j = 0; $j < length($aligned_source_seq); $j++) {
    my $char_query = substr($aligned_source_seq,$j,1);
    my $char_target = substr($aligned_target_seq,$j,1);
    if ($char_query eq '-' || $char_target  eq '-') {
      next;
    }
    if ($char_query eq $char_target) {
      $match_count++;
    }
    $aligned_positions++;
  }

  if ($aligned_positions <= 0) {
    throw("Pairwise alignment between the query sequence and the target sequence shows zero aligned positions. Something has gone wrong.");
  }
  my $percent_id = ($match_count/$aligned_positions)*100;
  $coverage = sprintf "%.2f", $coverage;
  $percent_id = sprintf "%.2f", $percent_id;
  return ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq);
}


=head2 map_cds_location

  Arg [0]   : source transcript with CDS
  Arg [1]   : target transcript without a CDS
  Arg [2]   : source transcript alignment seq from pairwise aligment of both transcripts
  Arg [3]   : target transcript alignment seq from pairwise aligment of both transcripts
  Function  : It takes in a source transcript with an annotated CDS and then tries to find that CDS in
              the target based on a pairwise alignement of both sequences. This works in conjunction
              with the output of align_nucleotide_seqs, which returns the aligned sequences from the
              pairwise alignment
  Returntype: None
  Examples  : align_nucleotide_seqs($source_transcript,$target_transcript,$aligned_source_seq,$aligned_target_seq);

=cut

sub map_cds_location {
  my ($source_transcript,$target_transcript,$aligned_source_seq,$aligned_target_seq) = @_;

  my $source_seq = $source_transcript->seq->seq();
  my $target_seq = $target_transcript->seq->seq();

  unless($source_transcript->translation()) {
    throw("Source transcript does not have a translation, so can't map a CDS");
  }

  my $source_cds_start = $source_transcript->cdna_coding_start();
  my $source_cds_end = $source_transcript->cdna_coding_end();

  my @align_source = split("",$aligned_source_seq);
  my @align_target = split("",$aligned_target_seq);

  my ($source_align_start_index,$target_start_pos) = find_base_alignment_index($source_cds_start,$aligned_source_seq,$aligned_target_seq);

  my $target_codon = substr($target_seq,$target_start_pos - 1, 3);
  say "Target codon: ".$target_codon;
  if($target_codon eq 'ATG') {
    my $target_stop_pos = find_stop_pos($target_seq,$target_start_pos);
    my $target_stop_codon = substr($target_seq,$target_stop_pos - 3, 3);
    say "Target stop codon: ".$target_stop_codon;
#    my @coords = $target_transcript->cdna2genomic($target_start_pos,$target_start_pos);
#    my $transcript_mapper = $target_transcript->get_TranscriptMapper();
#
#    say "Genomic start/end: ".${$coords}[0]."..".${$coords}[1];
#    use Data::Dumper;
#    say "FERGAL DUMPER: ".Dumper(@coords);
    my ($target_start_exon,$start_offset) = find_base_exon_position($target_transcript,$target_start_pos);
#    say "FERGAL T1";
    my ($target_end_exon,$end_offset) = find_base_exon_position($target_transcript,$target_stop_pos);
    # Once we have all of these it's enough to attach a translation to the target
#    say "FERGAL T2";
    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($target_start_exon);
    $translation->start($start_offset);
    $translation->end_Exon($target_end_exon);
    $translation->end($end_offset);
    $target_transcript->translation($translation);

    # Set the phases
    calculate_exon_phases($target_transcript, 0);
 #   say "Translation: ".$target_transcript->translation->seq();
 #   throw("FERGAL DEBUG");
  }

}


=head2 find_base_alignment_index

  Arg [0]   : the position of the base in the unaligned source sequence
  Arg [1]   : source transcript alignment seq from pairwise aligment of both transcripts
  Arg [2]   : target transcript alignment seq from pairwise aligment of both transcripts
  Function  : It takes in the position of a base in the source transcript (for example $transcript->cdna_coding_start)
              along with the two aligned sequences from the source and target transcripts and then tries to figure out
              the location of the base in the target transcript and the index of the base in the alignment
  Returntype: source_align_pos_index, target_pos
  Examples  : my ($source_align_pos_index,$target_pos) = find_base_alignment_index($source_base_index,$aligned_source_seq,$aligned_target_seq);

=cut

sub find_base_alignment_index {
  my ($source_base_index,$aligned_source_seq,$aligned_target_seq) = @_;

  my @align_source = split("",$aligned_source_seq);
  my @align_target = split("",$aligned_target_seq);

  my $source_align_pos = 0;
  my $source_align_pos_index = 0;
  my $target_pos = 0;
  # Loop and find the cds start index in the alignment
  for(my $i=0; $i<scalar(@align_source); $i++) {
    my $align_char_source = $align_source[$i];
    my $align_char_target = $align_target[$i];
    if($align_char_source ne '-') {
      $source_align_pos++;
    }

    if($align_char_target ne '-') {
      $target_pos++;
    }

    if($source_align_pos == $source_base_index) {
      $source_align_pos_index = $i;
      last;
    }
  }

  say "Source base align pos: ".$source_align_pos_index;
  say "Target base seq pos: ".$target_pos;
  return($source_align_pos_index,$target_pos);
}


=head2 find_stop_pos

  Arg [0]   : the target sequence
  Arg [1]   : the position of the orf start
  Function  : takes in the orf start pos and then
  Returntype: stop index
  Examples  : my $stop_index = find_stop_pos($target_seq,$target_start_pos);

=cut

sub find_stop_pos {
  my ($target_seq,$target_start_pos) = @_;

  my @target_seq_array = split("",$target_seq);

  # This is the second codon index
  my $second_codon_index = $target_start_pos - 1 + 3;
  my $last_codon_index = 0;
  for(my $i=$second_codon_index; $i<scalar(@target_seq_array) - 2; $i += 3) {
    my $codon = $target_seq_array[$i].$target_seq_array[$i+1].$target_seq_array[$i+2];
    $last_codon_index = $i+2;
    if($codon eq 'TAA' or $codon eq 'TAG' or $codon eq 'TGA') {
      last;
    }
  }

  # Add 1 since the position in seq coords is 1-based
  return($last_codon_index + 1);
}


=head2 find_base_exon_position

  Arg [0]   : target_transcript
  Arg [1]   : target_pos
  Function  : takes the target transcript and a base position in it and figures out where it is in
              exons
  Returntype: exon, offset in exon
  Examples  : my ($target_end_exon,$end_offset) = find_base_exon_position($target_transcript,$target_stop_pos);

=cut

#sub find_base_exon_position {
#  my ($transcript,$pos) = @_;

#  my $exons = $transcript->get_all_Exons();
#  my $running_count = 0;
#  my $selected_exon;
#  my $offset;
#  foreach my $exon (@$exons) {
#    $running_count += $exon->length();
#    if($running_count >= $pos) {
#      $selected_exon = $exon;
#      my $start_base_count = $running_count - $exon->length() + 1;
#      my $offset = $pos - $start_base_count + 1;
#    }
#  }

  # Stop examples
  # pos = 12
  #  123456789012345
  #  AAATTTAAATGATTT
  # rc = 15
  # sbc = 15 - 15  + 1 = 1
  # offset = 12 - 1 + 1 = 12
  #
  # pos = 21
  # 123456789 012345678901234
  # TTTAAATTT AAACCCAAATGAAAA
  # loop 1 (exon 1):
  # rc = 9
  # loop 2 (exon 2):
  # rc = 9 + 15 = 24
  # sbc = 24 - 15 + 1 = 10
  # offset = 21 - 10 + 1 = 12

#  unless($selected_exon and $offset) {
#    throw("Couldn't fine the position of the base in the exon set, something went wrong");
#  }

#  return($selected_exon,$offset);
#}


sub find_base_exon_position {
  my ($transcript,$pos) = @_;

  my $exons = $transcript->get_all_Exons();
  my @genomic_coords = $transcript->cdna2genomic($pos,$pos);
  my $pos_genomic = $genomic_coords[0]->start();
  say "Pos genomic: ".$pos_genomic;
  my $selected_exon;
  foreach my $exon (@$exons) {
    say "Exon: ".$exon->seq_region_start()."..".$exon->seq_region_end();
    if($pos_genomic >= $exon->seq_region_start and $pos_genomic <= $exon->seq_region_end) {
      $selected_exon = $exon;
      last;
    }
  }

  unless($selected_exon) {
    throw("Couldn't find the position of the base in the exon set, something went wrong");
  }

  my $offset = 0;
  if($selected_exon->strand() == 1) {
    $offset = $pos_genomic - $selected_exon->seq_region_start() + 1;
  } else {
    $offset = $selected_exon->seq_region_end() - $pos_genomic + 1;
  }

  return($selected_exon,$offset);
}


=head2 calculate_exon_phases

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : Given a transcript, calculates and sets
    exon phases according to translation and given
    start phase

=cut

sub calculate_exon_phases {
  my ($transcript, $start_phase) = @_;

  foreach my $e (@{$transcript->get_all_Exons}) {
    $e->phase(-1);
    $e->end_phase(-1);
  }

  if ($transcript->translation) {
    my $tr = $transcript->translation;

    my @exons = @{$transcript->get_all_Exons};

    while($exons[0] != $tr->start_Exon) {
      shift @exons;
    }
    while($exons[-1] != $tr->end_Exon) {
      pop @exons;
    }

    # set phase of for first coding exon
    my $cds_len = $exons[0]->length;
    if ($tr->start == 1) {
      $exons[0]->phase($start_phase);
      if ($start_phase > 0) {
        $cds_len += $start_phase;
      }
    } else {
      $cds_len -= ($tr->start - 1);
    }
    $exons[0]->end_phase($cds_len % 3);
    # set phase for internal coding exons
    for(my $i=1; $i < @exons; $i++) {
      $exons[$i]->phase($exons[$i-1]->end_phase);
      $exons[$i]->end_phase(($exons[$i]->length + $exons[$i]->phase) % 3);
    }

    # set phase for last coding exon
    if ($exons[-1]->length > $tr->end) {
      $exons[-1]->end_phase(-1);
    }
  }
}


=head2 hrdb_get_dba

 Arg [1]    : Hashref $connection_info, containing the connection details for the database:
              -host, -user, -dbname, -port [, -pass, -dna_db,...]
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor object (optional), the database will have the dna
 Arg [3]    : String $alternative_class (optional), Allowed class are Variation, Compara, Funcgen
 Example    : hrdb_get_dba->($self->param('target_db'));
 Description: It creates a object based on the information contained in $connection_info.
              If the hasref contains -dna_db or if the second argument is populated, it will
              try to attach the DNA database
 Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
 Exceptions : Throws if it cannot connect to the database.
              Throws if $connection_info is not a hashref
              Throws if $dna_db is not a Bio::EnsEMBL::DBSQL::DBAdaptor object

=cut

sub hrdb_get_dba {
  my ($connection_info, $dna_db, $alternative_class) = @_;

  my $dba;
  my %params;
  if(ref($connection_info) eq 'HASH') {
    my $module_name = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
    if ($alternative_class) {
      if ($alternative_class =~ /::/) {
        $module_name = $alternative_class;
      }
      else {
        $module_name = 'Bio::EnsEMBL::'.$alternative_class.'::DBSQL::DBAdaptor';
      }
      eval "use $module_name";
      if ($@) {
        throw("Cannot find module $module_name");
      }
    }
    eval {
      $dba = $module_name->new(%$connection_info, %params);
    };

    if($@) {
      throw("Error while setting up database connection:\n".$@);
    }
    if (!$dba->isa($module_name)) {
      warning("Hardcore blessing $dba into $module_name");
      weaken($dba);
      bless $dba, $module_name;
    }
  } else {
    throw("DB connection info passed in was not a hash:\n".$connection_info);
  }

  if (defined $dna_db) {
      if ($dna_db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
          $dba->dnadb($dna_db);
      }
      else {
          throw(ref($dna_db)." is not a Bio::EnsEMBL::DBSQL::DBAdaptor\n");
      }
  }

  return $dba;
}

=head2 convert_to_ucsc_name

 Arg [1]    : String $ensembl_name, an Ensembl seq_region name
 Arg [2]    : (optional) Bio::EnsEMBL::Slice Object slice, the slice you want to get the UCSC name from
 Example    : convert_to_ucsc_name($slice->seq_region_name, $slice);
 Description: It returns the UCSC name of the region by fetching the UCSC name from the seq_region_synonym table.
              If a Bio::EnsEMBL::Slice object is not provided or if it cannot find the synonym, it returns the
              Ensembl name prefixed with 'chr'.
 Returntype : String
 Exceptions : None


=cut

sub convert_to_ucsc_name {
    my ($ensembl_name, $slice) = @_;

    my $ucsc_name = 'chr'.$ensembl_name;
    if ($slice) {
        my $ucsc_synonyms = $slice->get_all_synonyms('UCSC');
        if (scalar(@$ucsc_synonyms)) {
            $ucsc_name = $ucsc_synonyms->[0]->name;
        }
    }
    return $ucsc_name;
}


=head2 get_analysis_settings

 Arg [1]    : String $class which represents the class of the config file you want to use
 Arg [2]    : String $key, the key value of the hash you want to retrieve
 Arg [3]    : Hashref $additional_data (optional), additional values to add or to overwrite the default/original values
 Arg [4]    : String $data_type (optional), to specify the type of data which will be return, default is hashref. Values are:
                ARRAY
                HASH
 Example    : my $config_hash = get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic', 'BlastGenscanPep');
 Description: Retrieve Blast, Exonerate,... configuration hash which are similar for most of the analyses like running raw computes
 Returntype : Reference, it will return an empty hashref/arrayref if the Config file does not exists or is not in PERL5LIB
 Exceptions : None

=cut

sub get_analysis_settings {
    my ($class, $key, $additional_data, $data_type) = @_;

    eval "use $class";
    if ($@) {
        if (defined $data_type and $data_type eq 'ARRAY') {
          return [];
        }
        else {
          return {};
        }
    }
    my $config = $class->new();
    if (defined $data_type and $data_type eq 'ARRAY') {
      return $config->get_array_config_settings($key, $additional_data);
    }
    else {
      return $config->get_config_settings($key, $additional_data);
    }
}


=head2 locate_executable

 Arg [1]    : String $executable, executable to locate
 Arg [2]    : String $bindir (optional), directory to look for the executable
 Description: First checks if Arg[1] is executable, if not checks if the name
              concatenated with Arg[2] is executable, if not then uses File::Which
              to search in PATH
 Returntype : String absolute path to the executable
 Exceptions : Throws if Arg[1] is not defined
              Throws if it couldn't find the path to the executable

=cut

sub locate_executable {
  my ($name, $bindir) = @_;

  my $path;
  if ($name) {
    if (-x $name) {
      $path = $name;
    }
    elsif ($bindir && -x catfile($bindir, $name)) {
      $path = catfile($bindir, $name);
    }
    else {
      $path = File::Which::which($name);
    }
    throw('Could not find the absolute path for '.$name) unless ($path);
  }
  else {
    throw("Must pass locate_executable a name if the program is to be located");
  }
  return $path;
}


=head2 first_upper_case

 Arg [1]    : String $string
 Description: Set the first letter of the string to upper case
 Returntype : String
 Exceptions : None

=cut

sub first_upper_case {
  my ($string) = @_;

  $string =~ s/^(\w)/\U$1/;
  return $string;
}


=head2 execute_with_wait

 Arg [1]    : String $cmd, the command to run
 Arg [2]    : String $msg (optional), message to the user if something goes wrong
 Arg [3]    : Int $wait (optional), how long we should wait before throwing an exception
              Default is 30 seconds
 Description: Execute a command and wait for Arg[3] seconds before throwing an execption
              This is really useful in Hive while using LSF as processes are killed in random order
              Without the wait a LSF signal like TERM_MEMLIMIT might not be caught before the exit
              code from the executable. Hive would not use the -1, -2 branches
 Returntype : Int 1 if successfull
 Exceptions : Throws after a Arg[3] seconds wait if the execution failed

=cut

sub execute_with_wait {
  my ($cmd, $failed_msg, $wait) = @_;

  if (system($cmd)) {
    sleep($wait || 30);
    throw($failed_msg || 'Failed to run with code: '.$?."\n".$cmd);
  }
  return 1;
}


=head2 execute_with_timer

 Arg [1]    : String $cmd, the command to run
 Arg [2]    : Int $timer, how long we should wait before killing your job
 Description: Execute a system command and kill it if it doesn't finish in time
              You can either specify the time in seconds as digits only or
              you can use M and H to specify hours and/or minutes, without white spaces.
 Returntype : Int 1 if successfull
 Exceptions : Throws after your jobs did not finish in time


=cut

sub execute_with_timer {
  my ($cmd, $timer) = @_;

  my $realtimer = parse_timer($timer);
  my $remaining_time = 0;

  # As seen on perldoc: http://perldoc.perl.org/5.14.2/functions/alarm.html
  eval {
    local $SIG{ALRM} = sub {die("alarm\n")};
    alarm $realtimer;
    execute_with_wait($cmd);
    $remaining_time = alarm 0;
  };
  if ($@) {
    if ($@ eq "alarm\n") {
      throw("Your job $cmd was still running after your timer: $realtimer\n");
    }
    else {
      throw($@);
    }
  }
  return $remaining_time;
}


=head2 parse_timer

 Arg [1]    : Int $timer, how long we should wait before killing your job
 Description: Parse a timer value of the form 2h30m or 1H50M or 3600 (seconds)
              You can either specify the time in seconds as digits only or

              you can use M and H to specify hours and/or minutes, without white spaces.
 Returntype : Value of the timer convered to seconds
 Exceptions : Throws if Arg[1] is not set or incorrect or 0

=cut

sub parse_timer {
  my ($timer) = @_;

  my $realtimer = 0;
  if ($timer) {
    if ($timer =~ tr/mhMH/MHMH/) {
      if ($timer =~ s/^(\d+)H//) {
        $realtimer = $1*3600;
      }
      if ($timer =~ s/^(\d+)M//) {
        $realtimer += $1*60;
      }
    }
    if ($timer) {
      if ($timer =~ /^\d*\s*$/) {
        $realtimer += $timer;
      }
      else {
        throw("This is what is left of your timer: $timer and this is what I could calculate: $realtimer\n");
      }
    }
  }
  else {
    throw("Something is wrong with your timer: $timer\n");
  }

  return($realtimer);
}



=head2 is_slice_name

 Arg [1]    : String, string to check
 Example    : is_slice_name($input_id);
 Description: Return 1 if the string given is an Ensembl slice name
 Returntype : Boolean
 Exceptions : None

=cut

sub is_slice_name {
    my ($string) = @_;

    return $string =~ /^[^:]+:[^:]+:[^:]+:\d+:\d+:(1|-1)$/;
}


=head2 get_database_from_registry

 Arg [1]    : String species name
 Arg [2]    : String type of the database, default is 'Core'
 Description: Use the Bio::EnsEMBL::Registry to get a database for the species
              Arg[1] of type Arg[2].
              If the databases are not released yet, it uses the databases
              from the previous release.
              This could be dangerous but needed in some cases.
              This method is for the test of modules
 Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
 Exceptions : None

=cut

sub get_database_from_registry {
  my ($species, $type) = @_;

  $type = 'Core' unless ($type);
  my $registry = 'Bio::EnsEMBL::Registry';
  my %hash = (
    -host => 'ensembldb.ensembl.org',
    -port => 3306,
    -user => 'anonymous',
  );
  if ($species) {
    $hash{'-species'} = $species;
  }
  $registry->load_registry_from_db(
    %hash,
  );
  my $db;
  eval {
    $db = $registry->get_DBAdaptor($species, $type);
  };
  if ($@) {
# Because the branching happens earlier I need to put this piece of code
# to make sure that we connect to the latest release.
    $registry->load_registry_from_db(
      %hash,
      -db_version => Bio::EnsEMBL::ApiVersion->software_version-1,
    );
    $db = $registry->get_DBAdaptor($species, $type);
  }
  return $db;
}


=head2 get_biotype_groups

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor, your database should have a biotype table
 Arg [2]    : String, database type, default to core
 Description: Retrieve all biotypes for a certain database type from the biotype table
              which is synchronised with the ensembl_production database.
 Returntype : Hashref, key is biotype, value is biotype_group
 Exceptions : Throws if Arg[1] is not a Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub get_biotype_groups {
  my ($db, $db_type) = @_;

  $db_type = 'core' unless ($db_type);
  throw('Bio::EnsEMBL::DBSQL::DBAdaptor needed, not '.ref($db)) unless (ref($db) eq 'Bio::EnsEMBL::DBSQL::DBAdaptor');
  my %biotype2group;

  # list all biotypes
  # and tag them by the group they belong to
  my $biotype_adaptor = $db->get_BiotypeAdaptor;
  foreach my $biotype (@{$biotype_adaptor->fetch_all}) {
    if ($biotype->{db_type} =~ /$db_type/ and $biotype->{db_type} !~ /$db_type\w+/) {
      $biotype2group{$biotype->name} = $biotype->biotype_group;
    }
  }

  return \%biotype2group;
}

=head2 get_feature_name

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: This method is mostly for the Ensembl RefSeq comparison scripts.
              It will find the stable id for the feature. If it is an Ensembl
              feature, it will return the stable id. If it is a RefSeq model
              if will try to get the display_xref then the EntrezGene id. If it
              cannot find anything it return the stable id which will be the id
              in the GFF file
 Returntype : String
 Exceptions : None

=cut

sub get_feature_name {
  my ($feature) = @_;

  my $name = $feature->stable_id;
  if (!($name =~ /^ENS\w+/ or $name =~ /^[NX][MR]_\d+/)) {
    if ($feature->display_xref) {
      $name = $feature->display_xref->display_id;
    }
    if (!$name) {
      my $dbentries = $feature->get_all_DBEntries('EntrezGene');
      if (@$dbentries) {
        $name = $dbentries->[0]->display_id;
      }
      if (!$name) {
        $name = $feature->stable_id;
      }
    }
  }
  return $name;
}


=head2 create_production_directory

 Arg [1]    : String, path of the directory to create
 Arg [2]    : Boolean (optiona), true to stripe the directory
 Arg [3]    : Int (optional), Permissions to set, need to start with 0
 Description: Create a directory in the filesystem. By default, it sets the
              permissions to 2775 (rwxrwsrx). The permissions can be given
              as Arg[3]. It will warn if the directory already exists but it
              will change the permissions.
              It is also possible to stripe the directory if it is created
              in LFS.
 Returntype : None
 Exceptions : Throws if the path cannot be created
              Throws is striping fails

=cut

sub create_production_directory {
  my ($path, $do_lsf_stripe, $mode) = @_;

  my $default_mode = 02775;
  $mode ||= $default_mode;
  if (-d $path) {
    warning("'$path' exists, only changing permissions");
  }
  else {
    mkdir $path;
  }
  chmod $mode, $path;
  if ($do_lsf_stripe) {
    if (!system("lfs getstripe $path")) {
      if(system("lfs setstripe  -c -1 $path")) {
        throw("Failed to stripe '$path'");
      }
    }
  }
}

1;
