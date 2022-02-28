# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ImportArrays;


=head1 SYNOPSIS

my $importer = 
  Bio::EnsEMBL::Analysis::RunnableDB::ImportArrays->new(

  );

$importer->fetch_input();
$importer->run();
$importer->output();
$importer->write_output(); #writes to DB and a big fasta file

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



=head1 AUTHOR

This module was written by Nathan Johnson, based on the CollapseAffy/Oligo code.

=head1 CONTACT

Post general queries to B<http://lists.ensembl.org/mailman/listinfo/dev>


=cut

# TO DO
# Write run methods for other import formats e.g. CSV?
# Add IMPORTED status for each array_chip?
# We can't do the same rollback with the probes on these array_chip
# as they may exist on another array.
# Write fasta file to db fasta cache?
# Add mode to recreate nr_fasta if already imported?

package Bio::EnsEMBL::Analysis::RunnableDB::ImportArrays;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeSet;
use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;

use Bio::EnsEMBL::Analysis::Config::ImportArrays;

use vars qw(@ISA);

$|=1;

#How do we get this to print warnings or to STDERR, it's getting caught somewhere!!!
#


@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);

############################################################

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config;
  
  #Set up and check DBs before we do the mapping
  $self->outdb->dbc->db_handle;

  #array names hash so we can print out a list of the actual names
  #i.e. not the names in the filename which can be formated differently.
  $self->{'_array_names'} = {}; #ArrayChip values

  #Have to set Helper here as it doesn't use the same param names
  #This is for checking the import status of the arrays and rolling back if required
  #This should enable failed imports to be picked

  $self->{'_helper'} = Bio::EnsEMBL::Funcgen::Utils::Helper->new( no_log => 1 );


  my $db = $self->db;

  return $self;
}


sub helper{
  return $_[0]->{'_helper'};
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
  #$self->non_redundant_probe_seqs($self->NON_REDUNDANT_PROBE_SEQS);
}


sub run{
  my $self = shift;

  #If we are dealing with fasta then just the one run method should do
  #If not then we need to write a different parse method for the input format
  #Need INPUT_FORMAT in config?

  my $method = 'run_'.$self->INPUT_FORMAT;#remove get from this method

  $self->$method;
}


#The idea is that we could handle other formats here.
#This overlaps with the idea of Parsers, so need to merge these at some point
#Fastas are not design files, parsers handle design files.

sub run_FASTA{
  my ($self) = @_;
 
  my (%probes_by_sequence, %array_chips, %probe_attrs, $probe_set);
  my ($current_array_chip, $existing_probe, $current_sequence, $sequence_fragment);
  
  #These get_CONFIGFIELD methods use the ARRAY_FORMAT to retrieve the correct config
  my $header_regex = $self->get_IIDREGEXP;
  my %valid_fields = 
    (
     -probe_set   => undef,
     -name        => undef,
     -array       => undef,
     -array_chip  => undef,
     -description => undef,
    );
  my %field_order = %{$self->get_IFIELDORDER};

  #Validate field
  foreach my $config_field (keys(%field_order)) {

    if (! exists $valid_fields{$config_field}) {
      throw("Found invalid field on ImportArrays.pm config:\t$config_field\n".
            "IFIELDORDER must only contain keys:\t".join("\t", keys %valid_fields));
    }
  }

  #This is done so we can dynamically assign fields using their match order
  #Don't actually use 4 and 5 here, but here just in case config regexps are extended
  my @match_refs = (\$1, \$2, \$3, \$4, \$5);
  #or is there a special array to hold these?

  #  warn "We need to count and report the number of probes imported, else throw if none imported\n".
  #"or if it doesn't match half the wc of the input file";
  
  open( PROBES, "<".$self->query_file);
 
  #todo count these in an array specific way
  my $cnt     = 0;
  my $nr_cnt  = 0;
  my $skipped_reps = 0;
  my $redundant_seq = 0;

  while (<PROBES>) {
    chomp;
	
    if (/$header_regex/) {     #This places header values into @match_refs
      $cnt++;

      #Annoyingly some probes have IDs of '0', but have differing probe seqs (and x;y; coords)
      #This is prevalent on AFFY_ST arrays e.g.
      #>probe:HuGene-1_0-st-v1:0;780:205; TranscriptClusterID=8180325; Sense; ProbeSetType=rescue->FLmRNA->unmapped
      #GGGCTGTCGCACACTGCACAGTTTG
      #>probe:HuGene-1_0-st-v1:0;956:935; TranscriptClusterID=8180325; Sense; ProbeSetType=rescue->FLmRNA->unmapped
      #CATGTCATGCGAAGGCAAGTCTGAA


      if ($current_sequence) {
		
        if (! $current_array_chip) {
          throw ("Have sequence $current_sequence but no current array chip!\n");
        }

        #You can only place the probe into the hash at this point, because
        #it's only at this point that you know you have the full sequence
        #Set probeset here as we delete from hash when creating probe


        $probe_set = (exists $probe_attrs{'-probe_set'}) ?  $probe_attrs{'-probe_set'} : 'NO_PROBESET';
        $existing_probe = $probes_by_sequence{$probe_set}{$current_sequence};
        #As this is keyed on the sequence it will not return probes with the same ID on the same 
        #array but with different sequence. (regardless of the probe_set_id)
        #Hence these, will be stored as separate probes, which will look identical apart from their probe_features

        if (! $existing_probe) {
          $nr_cnt++;
          #warn "Creating new probe:\t".$probe_attrs{'-name'}."\n";

          $existing_probe = $self->create_new_probe(
                                                    $current_array_chip, 
                                                    \%probe_attrs,
                                                    length($current_sequence),
                                                   );

          $probes_by_sequence{$probe_set}{$current_sequence} = $existing_probe;
        } else {
          #Sanity check here that it is not a technical replicate
          #as these will break the primary key of the probe table
          #due to the fact we don't model spatial differences here (e.g. x/y coords)
          #This is done implicitly for Affy arrays as they integrate the x/y coords into the probe names


          my $existing_name = $existing_probe->get_probename($probe_attrs{'-array'});

          if((defined $existing_name) && 
             ($existing_name eq $probe_attrs{'-name'}) ){
          
            #Must have found either a technical replicate with an identical name
            #Currently the model only allows one unique name per probeset per array

            #Check $current_sequence matches here also?
            warn "Skipping identical technical replicate probe:\t".
              $probe_attrs{'-name'}.' '.$current_sequence."\n";
            $skipped_reps ++;  
          }
          else {
           
            if(defined ($existing_name)){
              warn "Found replicate probe for ($existing_name) replicate:\t".
                $probe_attrs{'-name'}."\n";
              $redundant_seq ++;
            }
            
            $existing_probe->add_array_chip_probename($current_array_chip->dbID, 
                                                      $probe_attrs{'-name'}, 
                                                      $current_array_chip->get_Array);
          }
        }

        $current_sequence = undef;
      }


      #Set probe attrs using field order config
      foreach my $field (keys %field_order) {

        #Skip empty string as we prefer NULLs in the DB
        if(${$match_refs[$field_order{$field}]} ne ''){

          $probe_attrs{$field} = ${$match_refs[$field_order{$field}]};
        }
      }

      #Set current array chip
      $current_array_chip = $array_chips{$probe_attrs{-array_chip}};

      #Create new array chip
      if (! $current_array_chip) {		
        $current_array_chip = $self->create_new_array_chip($probe_attrs{-array}, $probe_attrs{-array_chip});
        $array_chips{$probe_attrs{-array_chip}} = $current_array_chip;
      }

    } elsif (/^[NSWRYKMBDHVATGCU]+$/i) { #build sequence
      $sequence_fragment = $_;
  
      if ($current_sequence) {
        $current_sequence = $current_sequence.$sequence_fragment;
      } else {
        $current_sequence = $sequence_fragment;
      }
    } else {
      throw('Found header which does not match '.$self->INPUT_FORMAT." regex($header_regex):\n$_");
    }
  }

  #deal with last entry
  $array_chips{$probe_attrs{-array_chip}} = $current_array_chip;
  $probe_set = (exists $probe_attrs{'-probe_set'}) ? $probe_attrs{'-probe_set'} : 'NO_PROBESET';
  $existing_probe = $probes_by_sequence{$probe_set}{$current_sequence};

  if (! $existing_probe) {
    $nr_cnt++;
    $existing_probe =  $self->create_new_probe(
                                               $current_array_chip, 
                                               \%probe_attrs,
                                               length($current_sequence),
                                              );

    $probes_by_sequence{$probe_set}{$current_sequence} = $existing_probe;
  } else {

    my $existing_name = $existing_probe->get_probename($probe_attrs{'-array'});

    if($existing_name && 
       ($existing_name eq $probe_attrs{'-name'})){
      
      warn "Skipping identical technical replicate probe:\t".
        $probe_attrs{'-name'}.' '.$current_sequence."\n";
      $skipped_reps ++;  
    }
    else{

      if($existing_name){
        warn "Found replicate probe for ($existing_name) replicate:\t".
                $probe_attrs{'-name'}."\n";
              $redundant_seq ++;
      }
      
      $existing_probe->add_array_chip_probename($current_array_chip->dbID, 
                                                $probe_attrs{'-name'}, 
                                                $current_array_chip->get_Array);
    }
  }  

  print "Seen $cnt fasta records\n";
  print "Created $nr_cnt probes\n";
  print "Skipped $skipped_reps identical replicates\n";
  print "Found $redundant_seq replicates with non-unique names\n";
  

  throw('No probes stored! Maybe you need to tweak your IIDREGEXP config for this format') if($nr_cnt == 0);
  $self->probes(\%probes_by_sequence);
  return;
}



#can we replace these var with $_[n] for speed?
sub create_new_probe {
  my($self, $array_chip, $probe_hash, $length) = @_;
  
  #can we set a template hash with the array and array chip id in?
  $probe_hash->{'-length'} = $length;
  $probe_hash->{'-class'}  = 'EXPERIMENTAL';#Will they all be experimental?
  $probe_hash->{'-array_chip_id'} = $array_chip->dbID;
  $probe_hash->{'-array'}  = $array_chip->get_Array;
  #remove probeset as this is not ProbeSet yet
  delete $probe_hash->{-probe_set};
 
  return Bio::EnsEMBL::Funcgen::Probe->new(%{$probe_hash});
}


sub create_new_array_chip {
  my($self, $array_name, $design_id) = @_;

  if(! ($array_name && $design_id)){
    throw('Need to pass an Array name and an ArrayChip design ID');
  }

  
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
  my $achip_adaptor = $self->outdb->get_ArrayChipAdaptor;
  my $array_params = $self->get_ARRAY_PARAMS_by_array_name($array_name);

  if(exists $array_params->{skip_config}){
    $self->{skip_config}{$array_name} = $array_params->{skip_config};
    delete $array_params->{skip_config};
  }

  my $array = $array_adaptor->fetch_by_name_vendor($array_name, $array_params->{'-vendor'});
  my $array_chip;

  if(! defined $array){
    $array = Bio::EnsEMBL::Funcgen::Array->new(%{$array_params});
    ($array) = @{$self->outdb->get_ArrayAdaptor->store($array)};
  }
  
  
  if($array_chip = $achip_adaptor->fetch_by_array_design_ids($array->dbID, $design_id)){

	if($array_chip->has_status('IMPORTED')){
	  throw("$array_name ArrayChip has already been IMPORTED. Please rollback_ArrayChip or recreate your arrays_nr.".$self->ARRAY_FORMAT.'.fasta file for alignment');
	}
	else{#Rollback
	  $self->helper->rollback_ArrayChips([$array_chip], 'probe');#, 'force');
	  #should we force roll back here?
	  #If not we may have to manually remove probe2transcript xrefs first
	  #Forcing here will mean that we are silently deleting the probe2transcript xrefs
	  #before we have a chance to back up
	  #should we env var this?
	}
  }
  else{
	#override design_id name default if we actually have a design_id in the config
	if($array_params->{'-design_id'}){
	  $design_id = $array_params->{'-design_id'};
	}

	$array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new
	  (
	   -name => $array->name,
	   -design_id => $design_id,#This is a placeholder as we may have the real ArrayChip details
	   -array_id  => $array->dbID,
	  );
	
	($array_chip) = @{$self->outdb->get_ArrayChipAdaptor->store($array_chip)};
  }

  if(! exists $self->{'_array_names'}->{$array_name}){
	$self->{'_array_names'}->{$array->name} = $array_chip;
	$self->{'_arrays'}->{$array->name} = $array;
  }



  return $array_chip;
}
############################################################

sub write_output {
  my ( $self, @output ) = @_;

  my $outdb = $self->outdb;
  my $outfile = $self->NON_REDUNDANT_PROBE_SEQS;
  my $probe_adaptor = $outdb->get_ProbeAdaptor;
  my $probeset_adaptor = $outdb->get_ProbeSetAdaptor;
  
  open (OUTFILE, ">".$outfile) || throw("Failed to open ouput file:\t".$outfile);
    
  #Store all probes
  foreach my $probeset (keys %{$self->probes}) {	
    my %probes = %{$self->probes->{$probeset}};

    if ($probeset ne 'NO_PROBESET') {
      $probeset = Bio::EnsEMBL::Funcgen::ProbeSet->new
        (
         -name => $probeset,
         -size => scalar(values %probes),
         #-array_chip => #do we need to add array_chip_id to the probeset table?
         #-family => ?,
        );
      ($probeset) = @{$probeset_adaptor->store($probeset)};
    }
  
    foreach my $sequence (keys %probes) {
      
      my $probe = $probes{$sequence};
      $probe->probeset($probeset) if $probeset ne 'NO_PROBESET';
      
      ($probe) = @{$probe_adaptor->store($probe)};
      
      print OUTFILE ">".$probe->dbID."\n".$sequence."\n";
    }
  }

  close(OUTFILE);


  #No we record imported as we need the nr_fasta dump?
  #Or are we going to allow this to run if already import
  #and just create the nr_fasta by querying the db for probe dbID
  #and using the sequence from the input file?


  #And now the array names and states
  $outfile = $self->NAMES_FILE;
  open (OUTFILE, ">".$outfile) || throw("Failed to open ouput file:\t".$outfile);


  print "Setting IMPORTED & DISPLAYABLE states for:\n";

  #MART_DISPLAYABLE & DISPLAYABLE should be on Array level
  #IMPORTED should be on ArrayChip
  #This all looks good, but everything for AGILENT re-imported is on array_chip?!!

  foreach my $aname(keys %{$self->{'_array_names'}}){
    print "\t".$aname."\n";
    print OUTFILE $aname."\n";
    $self->{'_array_names'}->{$aname}->add_status('IMPORTED');
    $self->{'_arrays'}->{$aname}->add_status('DISPLAYABLE');#This should be on Array?!
    #$self->{'_arrays'}->{$aname}->add_status('MART_DISPLAYABLE');#Now done in probe2transcript

    $self->{'_array_names'}->{$aname}->adaptor->store_states($self->{'_array_names'}->{$aname});	
    $self->{'_arrays'}->{$aname}->adaptor->store_states($self->{'_arrays'}->{$aname});	
  }


  close(OUTFILE);

  return;
}

############################################################

#problem here is that the super new is callign this before we get a chance to define it here
#rename method
#is db is super pipeline db?

sub outdb {
  my ($self) = @_;

  #Don't need DNADB for collapse
  #But need to define as will try and find a default on ensembldb which may not exist?

  my ($outdb);



  #Where are db and dnadb methods coming from?
  #Is this defaulting to the pipeline DB?
 
  if(! defined $self->{'_outdb'}){

	#This DNADB testing is a work around to avoid having to edit Config::ProbeAlign
	my $dnadb;

	#not defined as an empty env var will give the defined null string
	if($self->DNADB->{-dbname}){
	  $dnadb =  new Bio::EnsEMBL::DBSQL::DBAdaptor(%{ $self->DNADB });
	}

	$self->{'_outdb'} = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
	  (
	   %{ $self->OUTDB }, 
	   -dnadb => $dnadb,
	  );
   
	if(! $self->DNADB->{-dbname}){
	  print "WARNING: Using default DNADB ". $self->{'_outdb'}->dnadb->dbname."\n";
	}
	
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

  my ($array_format) = split/\:/, $self->input_id;
  $self->ARRAY_FORMAT($array_format);#Do we even need this now?


  # check that compulsory options have values
  #removed  ARRAY_FORMAT_FILE
  foreach my $config_var (
						  qw(
							 OUTDB
							 OUTPUT_DIR
							 IIDREGEXP
							 INPUT_FORMAT
							 IFIELDORDER
							 IIDREGEXP
							 IFIELDORDER
							 INPUT_FORMAT
							)
						 ){

	if ( ! defined $self->$config_var ){
      throw("You must define $config_var in config for logic '$logic'");
    }
  }

  
  $self->QUERYSEQS($self->OUTPUT_DIR."/arrays.${array_format}.fasta");
  $self->NON_REDUNDANT_PROBE_SEQS($self->OUTPUT_DIR."/arrays_nr.${array_format}.fasta");
  $self->NAMES_FILE($self->OUTPUT_DIR."/arrays.${array_format}.names");
 


  #Now defined as separate logic name config in ImportArrays config
  #now check for ARRAY_FORMAT defined config
  #we could set this here too
  #foreach my $config_var ('IIDREGEXP', 'IFIELDORDER', 'INPUT_FORMAT'){
	
#	if ( ! defined $self->{"_CONFIG_${config_var}"}{$self->ARRAY_FORMAT}){
#	  throw('You must define a '.$self->ARRAY_FORMAT." entry in the $config_var config for logic '$logic'");
#	}
#  }


  #Why does ouput DB not have to be defined? We are actually checking for it above
  #What about DNADB?
  
  #remove?
  # output db does not have to be defined, but if it is, it should be a hash
  if ( ref( $self->OUTDB ) ne "HASH"  || ! defined $self->OUTDB->{-dbname}) {
    throw("OUTDB in config for '$logic' must be a hash ref of db connection pars.");
  }


}

### Dynamic config
#Now built from input id

sub ARRAY_FORMAT {
  my ( $self, $value ) = @_;

  $self->{'_ARRAY_FORMAT'} = $value if defined $value;
  return $self->{'_ARRAY_FORMAT'};
}

sub NON_REDUNDANT_PROBE_SEQS {
  my ( $self, $value ) = @_;

  $self->{'_NON_REDUNDANT_PROBE_SEQS'} = $value if defined $value;
  return $self->{'_NON_REDUNDANT_PROBE_SEQS'};
}


sub QUERYSEQS{
  my ( $self, $value ) = @_;

  $self->{'_QUERYSEQS'} = $value if defined $value;
  return $self->{'_QUERYSEQS'};
}

sub NAMES_FILE{
  my ( $self, $value ) = @_;

  $self->{'_NAMES_FILE'} = $value if defined $value;
  return $self->{'_NAMES_FILE'};
}

### Config from ImportArrays ##

sub INPUT_FORMAT{
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_INPUT_FORMAT'} = $value  if defined $value ;
  
  return $self->{'_CONFIG_INPUT_FORMAT'};
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

  return $self->{'_CONFIG_IIDREGEXP'};
}

sub get_IFIELDORDER{
  my $self = shift;

  return $self->{'_CONFIG_IFIELDORDER'};
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

sub OUTPUT_DIR {
  my ( $self, $value ) = @_;

  $self->{'_CONFIG_OUTPUT_DIR'} = $value  if defined $value;
  return $self->{'_CONFIG_OUTPUT_DIR'};
}

#Only used to define dnadb if not on ensembldb
#else dnadb autoguessing will fail
sub DNADB {
  my ( $self, $value ) = @_;


  $self->{'_CONFIG_DNADB'} = $value  if defined $value;
  return $self->{'_CONFIG_DNADB'};
}



###############################################
###     end of config
###############################################

1;
