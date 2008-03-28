
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ImportArrays;


=head1 SYNOPSIS

my $affy = 
  Bio::EnsEMBL::Analysis::RunnableDB::ImportArrays->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
    -database   => $EST_GENOMIC,
    -query_seqs => \@sequences,
  );

$affy->fetch_input();
$affy->run();
$affy->output();
$affy->write_output(); #writes to DB and a big fasta file

=head1 DESCRIPTION

This object imports array fasta files of a common array format and collapses redundant 
probes into unique records based on their probeset name and sequence identity.  The 
probes supplied are redundant - a single probe (characterised by a unique sequence) may 
occur in different positions of different arrays (i.e. for AFFY the chip coord is the 
probe name). Parsing of the fasta file is facilitated by dynamic configuration of regular 
expressions and array parameters based on the array format and array being imported. To 
enable multiple instances of ImportArrays to run at the same time using different 
configurations, the correct config type is read from an ImportArrays.config file, 
using the input_id as the key.

A non-redundant array format specific output file is written, using the probe dbIDs in 
the fasta headers.  This can then be cat'd with out no-redundant array format fasta file 
for use in the mapping carried out by ProbeAlign and ProbeTranscriptAlign.


Note that probes are defined as redundant when they share the same
- sequence and 
- probeset 


To do:

- Write run methods for other import formats e.g. CSV?
- Add recover and rollback functionality so we can skip Arrays which have been fully imported and rollback ones where we fell over half way through.


=head1 AUTHOR

This module was written by Nathan Johnson, based on the CollapseAffy/Oligo code.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>


=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ImportArrays;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeSet;
use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Analysis::Config::ImportArrays;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);

############################################################
sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config;
  
  return $self;
}

sub fetch_input {
  my ($self) = @_;

  my $logic = $self->analysis->logic_name;

  my ($query_file, $chunk_number, $chunk_total);

  my $query = $self->QUERYSEQS;

  if ( -e $query and -d $query ) {
    throw "I need to have all affy probes input in one big file\n";
  } elsif ( -e $query and -s $query ) {
    #store this for reference later
    $self->query_file($query);
  } else {
    throw("'$query' refers to something that could not be made sense of\n");
  }


  #Why do we bother setting this?  Why don't we just use the config methods directly?
  $self->non_redundant_probe_seqs($self->NON_REDUNDANT_PROBE_SEQS);
}


sub run{
  my $self = shift;

  #If we are dealing with fasta then just the one run method should do
  #If not then we need to write a different parse method for the input format
  #Need INPUT_FORMAT in config?

  my $method = 'run_'.$self->get_INPUT_FORMAT;

  $self->$method;
}


sub run_FASTA{
  my ($self) = @_;
 
  my (%probes_by_sequence, %array_chips, %probe_attrs, $probeset);
  my ($current_array_chip, $existing_probe, $current_sequence, $sequence_fragment);

  #These get_CONFIGFIELD methods use the ARRAY_FORMAT to retrieve the correct config
  my $header_regex = $self->get_IIDREGEXP;
  my %field_order = %{$self->get_IFIELDORDER};

  #This is done so we can dynamically assign fields using their match order
  #Don't actually use 4 and 5 here, but here just in case config regexps are extended
  my @match_refs = (\$1, \$2, \$3, \$4, \$5);
  #or is there a special array to hold these?
  
  open( PROBES, "<".$self->query_file);
 
  while(<PROBES>){
    chomp;
   
    if(/$header_regex$/){
	
      if($current_sequence){

        if(! $current_array_chip){
          throw ("Have sequence $current_sequence but no current array chip!\n");
        }

		#You can only place the probe into the hash at this point, because
        #it's only at this point that you know you have the full sequence
		$probeset = (exists $probe_attrs{'-probeset'}) ? $probe_attrs{'-probeset'} : 'NO_PROBE_SET';
        $existing_probe = $probes_by_sequence{$probeset}{$current_sequence};

        if(! $existing_probe){

          $existing_probe = $self->create_new_probe(
													$current_array_chip, 
													#$probeset, 
													#$probe_name, 
													\%probe_attrs,
													length($current_sequence),
												   );

          $probes_by_sequence{$probeset}{$current_sequence} = $existing_probe;
        }
		else{
        
		  $self->add_array_to_existing_probe(
											 $existing_probe, 
											 $current_array_chip, 
											 #$probeset,
											 #$probe_name ,
											 \%probe_attrs,
											);
        }

        $current_sequence = undef;
      }


	  #set probe attrs using field order config
	  #Do not use map as this is slow
	  foreach my $field(keys %field_order){
		$probe_attrs{$field} = ${$match_refs[$field_order{$field}]};
	  }

	  #Set current array chip
      $current_array_chip = $array_chips{$probe_attrs{-array_chip}};

	  #Create new array chip
      if(! $current_array_chip){
        $current_array_chip = $self->create_new_array_chip($probe_attrs{-array_chip});
        $array_chips{$probe_attrs{-array_chip}} = $current_array_chip;
      }

    }
	else{#build sequence
      $sequence_fragment = $_;

      if($current_sequence){
        $current_sequence = $current_sequence . $sequence_fragment;
      }
	  else{
        $current_sequence = $sequence_fragment;
      }
    }
  }

  #deal with last entry
  $array_chips{$probe_attrs{-array_chip}} = $current_array_chip;
  $probeset = (exists $probe_attrs{'-probeset'}) ? $probe_attrs{'-probeset'} : 'NO_PROBE_SET';
  $existing_probe = $probes_by_sequence{$probeset}{$current_sequence};

  if(! $existing_probe){
    $existing_probe =  $self->create_new_probe(
											   $current_array_chip, 
											   #$probeset, 
											   #$probe_name, 
											   \%probe_attrs,
											   length($current_sequence),
											  );

    $probes_by_sequence{$probeset}{$current_sequence} = $existing_probe;
  }
  else{
    $self->add_arraychip_to_existing_probe(
									   $existing_probe, 
									   $current_array_chip, 
									   #$probeset,
									   #$probe_name,
									   \%probe_attrs,
									  );
  }  

  $self->probes(\%probes_by_sequence);

  return;
}

sub add_array_chip_to_existing_probe{#_probe_set{
  my ($self, $probe, $array_chip, $probeset, $probename) = @_;

  my @all_probenames = @{$probe->get_all_probenames};
      
  if(!($probe->probeset eq $probeset)){
    throw (
      "Inconsistency: have found a probe ".$array_chip->name.":$probeset:$probename with ".
      "identical sequence but different probeset to another one: ".$probe->probeset."\n"
    );
  }

  $probe->add_array_chip_probename($array_chip->dbID, $probename, $array_chip->get_Array);
}

#can we replace these var with $_[n] for speed?
sub create_new_probe {
  my($self, $array_chip, $probe_hash, $length) = @_;
  
  #can we set a template hash with the array and array chip id in?
  $probe_hash->{'-length'} = $length;
  $probe_hash->{'-class'}  = 'EXPERIMENTAL';#Will they all be experimental?
  $probe_hash->{'-array_chip_id'} = $array_chip->dbID;
  $probe_hash->{'-array'}  = $array_chip->get_Array;

  return Bio::EnsEMBL::Funcgen::Probe->new(%{$probe_hash});

  #This causes problem as the params names are not the same as the attr fields
  #return Bio::EnsEMBL::Funcgen::Probe->_new_fast(%{$probe_attrs});
}

sub create_new_array_chip {
  my($self, $array_name) = @_;


  
  #Let's make the assumption that related arrays will have the same format header
  #Therefore we don't have to do any heuristics to get the parse regex and hence the name of the array
  #This would assume that array names are unique across vendors
  #Very probably, and can just comment out duplicate names on running
  #Need to check this in validate? Will two keys of the same value be returned, pointing to noe set of data?
  #Should also check that all arrays in arrays.list are present in config
  #Or can we just do this in a text file which we parse and populate the config with?
  #This way we don't have to keep editing the code
  #And we can test for duplications in names easily
  #Would need to check in the config file separately for each run
  #Can we automate this, with a comment which is the pid or something to link it to the the log?
  #Or maybe we can just append the config to the log
  #Implement Helper.pm!
  #Is this not reinventing the wheel slightly?
  
  #We also need to alter config to set arrays which are to be imported together
  #We can validate each one, by just taking the first record and testing the parse regex works


  my $array_adaptor = $self->outdb->get_ArrayAdaptor;
  my $array_params = $self->get_ARRAY_PARAMS_by_array_name($array_name);

  if($array_adaptor->fetch_by_name_vendor($array_name, $array_params->{'-vendor'})){

	#Could test for status here and rollback or skip?
	#return no array_chip if it is already fully import

	throw("Array $array_name already exists in the DB, please remove from DB before attempting import");
  }

  my $array = Bio::EnsEMBL::Funcgen::Array->new(%{$array_params});
  ($array) = @{$self->outdb->get_ArrayAdaptor->store($array)};

  my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new
	(
	 -name => $array->name,
	 -design_id => $array->name,#This is a placeholder as we don't have the real ArrayChip details
	 -array_id  => $array->dbID,
	);


  ($array_chip) = @{$self->outdb->get_ArrayChipAdaptor->store($array_chip)};

  return $array_chip;
}
############################################################

sub write_output {
  my ( $self, @output ) = @_;

  my $outdb = $self->outdb;
  my $outfile = $self->non_redundant_probe_seqs;
  my $probe_adaptor = $outdb->get_ProbeAdaptor;
  my $probeset_adaptor = $outdb->get_ProbeSetAdaptor;
  
  open (OUTFILE, ">".$outfile) || throw("Failed to open ouput file:\t".$outfile);
    
  #Store all probes
  foreach my $probeset(keys %{$self->probes}){
	
	my %probes = %{$self->probes->{$probeset}};

	if($probeset ne 'NO_PROBE_SET'){
	  $probeset = Bio::EnsEMBL::Funcgen::ProbeSet->new
		(
		 -name => $probeset,
		 -size => scalar(values %probes),
		 #-array_chip => #do we need to add array_chip_id to the probeset table?
		 #-family => ?,
		);

	  ($probeset) = @{$probeset_adaptor->store($probeset)};
	}
	else{
	  undef $probeset;
	}

	foreach my $sequence(keys %probes){
	  my $probe = $probes{$sequence};
	  $probe->probeset($probeset);
	  
	  ($probe) = @{$probe_adaptor->store($probe)};

	  print OUTFILE ">".$probe->dbID."\n".$sequence."\n";
	}
  }

  close(OUTFILE);
}

############################################################

sub outdb {
  my ($self) = @_;

  my ($outdb, $dnadb);

  #Where are db and dnadb methods coming from?
  #Is this defaulting to the pipeline DB?

  if(! defined $self->{'_outdb'}){

	if( my $dnadb_params =  $self->DNADB ){
	  #Need to test if they have been defined
	  #Just test dbname
	  #All blank env vars are defined as the null string '';
	  
	  if($dnadb_params->{-dbname}){
		$dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{$dnadb_params});
	  }
	}  #else will use default ensembldb dnadb, this may fail for oddly named dbs
	
	
	if ( my $outdb_params = $self->OUTDB ) {
	  #We need to set the dnadb here if we are dealing with a pre release DB, i.e. dnadb is not on ensembldb
	  $outdb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new( %{ $outdb_params }, -dnadb => $dnadb);
	} 
	else {
	  #This is the pipeline DBAdaptor?
	  $outdb = $self->db;
	  $outdb->dnadb($dnadb) if $dnadb;
	}

	$self->{'_outdb'} = $outdb;
  }

  return $self->{'_outdb'};
}

############################################################
#
# The QUERYSEQS config variable.
#
sub query_file {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_query_file'} = $value;
  }

  if ( exists( $self->{'_query_file'} ) ) {
    return $self->{'_query_file'};
  } else {
    return undef;
  }
}

#############################################################


sub non_redundant_probe_seqs {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_non_redudant_probe_seqs'} = $value;
  }

  if ( exists( $self->{'_non_redudant_probe_seqs'} ) ) {
    return $self->{'_non_redudant_probe_seqs'};
  } else {
    return undef;
  }
}


#############################################################

sub probes {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_probes'} = $value;
  }

  if ( exists( $self->{'_probes'} ) ) {
    return $self->{'_probes'};
  } else {
    return undef;
  }
}


#############################################################
# Declare and set up config variables
#############################################################

sub read_and_check_config {
  my $self = shift;

  #This global variable is defined in the ImportArrays config
  $self->SUPER::read_and_check_config($ARRAY_CONFIG);

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
  foreach my $config_var (
    qw(
	   OUTDB
	   DB_HOME
	   IIDREGEXP
	   ARRAY_FORMAT_FILE
	   INPUT_FORMAT
	   IFIELDORDER
	  )
  ){
    if ( ! defined $self->$config_var ){
      throw("You must define $config_var in config for logic '$logic'");
    }
  }

  
  open(ARRAY_FORMAT_CONFIG, '<'.$self->ARRAY_FORMAT_FILE) || throw("Could not open ARRAY_FORMAT_FILE:\t"..$self->ARRAY_FORMAT_FILE);

  my @config = <ARRAY_FORMAT_CONFIG>;
  my $input_id = $self->input_id;
  my ($array_config) = grep(/^$input_id/, @config);
  my (undef, $array_format, $raw_fasta) = split/\s+/, $array_config;
  $self->ARRAY_FORMAT($array_format);
  $self->QUERYSEQS($raw_fasta);
  $self->NON_REDUNDANT_PROBE_SEQS($self->DB_HOME."/arrays_nr.${array_format}.fasta");
 
  #NON_REDUNDANT_PROBE_SEQS is now a format specific file which will need cating after we have finished importing
  #This should be built from the DB_HOME and the ARRAY_FORMAT.
  
  #ARRAY_FORMAT QUERYSEQS(RAW_FASTA) are retrieved from the IportArrays.config file
  #Key is input_id


  #now check for ARRAY_FORMAT defined config
  #we could set this here too
  foreach my $config_var ('IIDREGEXP', 'IFIELDORDER', 'INPUT_FORMAT'){
	
	if ( ! defined $self->{"_CONFIG_${config_var}"}{$self->ARRAY_FORMAT}){
	  throw('You must define a '.$self->ARRAY_FORMAT." entry in the $config_var config for logic '$logic'");
	}
  }


  #Why does ouput DB not have to be defined? We are actually checking for it above
  #What about DNADB?
  
  #remove?
  # output db does not have to be defined, but if it is, it should be a hash
  if ( ref( $self->OUTDB ) ne "HASH" ) {
    throw("OUTDB in config for '$logic' must be a hash ref of db connection pars.");
  }

  #Test for DNADB here and warn?

}

### Dynamic config

sub ARRAY_FORMAT {
  my ( $self, $value ) = @_;

  $self->{'_ARRAY_FORMAT'} = $value if defined $value;
  return $self->{'_ARRAY_FORMAT'};
}

sub QUERYSEQS{
  my ( $self, $value ) = @_;

  $self->{'_QUERYSEQS'} = $value if defined $value;
  return $self->{'_QUERYSEQS'};
}


# This is now format specific file, which will need cating after the individual Imports have completed

sub NON_REDUNDANT_PROBE_SEQS {
  my ( $self, $value ) = @_;

  $self->{'_NON_REDUNDANT_PROBE_SEQS'} = $value if defined $value;
  return $self->{'_NON_REDUNDANT_PROBE_SEQS'};
}



### Config from ImportArrays

#Why are we prefixing everything with config?


sub INPUT_FORMAT{
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_INPUT_FORMAT'} = $value  if defined $value ;
  
  return $self->{'_CONFIG_INPUT_FORMAT'};
}

sub get_INPUT_FORMAT{
  my $self = shift;

  return $self->{'_CONFIG_INPUT_FORMAT'}{$self->ARRAY_FORMAT};
}


sub ARRAY_FORMAT_FILE {
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_ARRAY_FORMAT_FILE'} = $value  if defined $value ;
  
  return $self->{'_CONFIG_ARRAY_FORMAT_FILE'};
}


sub ARRAY_PARAMS {
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_ARRAY_PARAMS'} = $value  if defined $value;
  return $self->{'_CONFIG_ARRAY_PARAMS'};
}


#Don't really need this unless we see probe records not in consecutive array blocks
#This accessor will slow down import
#Can we implement this as a local cache in the caller?
#sub get_array_by_name{
#  my ($self, $array_name)
#
#}

sub get_ARRAY_PARAMS_by_array_name {
  my ( $self, $array_name ) = @_;

  #Should cache the 
  #Need to test for array name here

  if(! exists $self->{'_CONFIG_ARRAY_PARAMS'}{$array_name}){
	throw("No ARRAY_PARAMS config available for $array_name.  You must add this to the ImportArrays config before importing");
  }


  return $self->{'_CONFIG_ARRAY_PARAMS'}{$array_name};
}


sub IIDREGEXP {
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_IIDREGEXP'} = $value  if defined $value;
  return $self->{'_CONFIG_IIDREGEXP'};
}

sub get_IIDREGEXP{
  my $self = shift;

  return $self->{'_CONFIG_IIDREGEXP'}{$self->ARRAY_FORMAT};
}

sub get_IFIELDORDER{
  my $self = shift;

  return $self->{'_CONFIG_IFIELDORDER'}{$self->ARRAY_FORMAT};
}



sub IFIELDORDER {
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_IFIELDORDER'} = $value if defined $value;
  return $self->{'_CONFIG_IFIELDORDER'};
}

sub OUTDB {
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_OUTDB'} = $value  if defined $value;
  return $self->{'_CONFIG_OUTDB'};
}


sub DNADB {
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_DNADB'} = $value  if defined $value;
  return $self->{'_CONFIG_DNADB'};
}

sub DB_HOME {
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_DB_HOME'} = $value  if defined $value;
  return $self->{'_CONFIG_DB_HOME'};
}


###############################################
###     end of config
###############################################

1;
