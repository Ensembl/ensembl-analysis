=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::CollapseAffyProbes - 

=head1 SYNOPSIS

my $affy = 
  Bio::EnsEMBL::Analysis::RunnableDB::CollapseAffyProbes->new(
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

This object runs the first step in the process for mapping affymetrix probes to a genome.
The probes supplied are redundant - a single probe (characterised by a unique sequence) 
may occur in different positions of different arrays. SO we will first collect all 
redundant sequences, find a non-redudant subset, and then simultaneously 
-- write the probes to our db, and
-- write an output file of non-redundant sequences to a flat file, so that
exonerate can map them against the genome in the next step.

Note that probes are defined as redundant when they share the same
- sequence and 
- probeset 

=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::CollapseAffyProbes;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::AffyProbe;
use Bio::EnsEMBL::AffyArray;

use Bio::EnsEMBL::Analysis::Config::CollapseAffyProbes;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);

############################################################
sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($AFFY_CONFIG);
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

  $self->non_redundant_probe_seqs($self->NON_REDUNDANT_PROBE_SEQS);
}

sub run {
  my ($self) = @_;
  my %probes_by_sequence;
  my %affy_arrays;
  
  open( PROBES, "<".$self->query_file);

  my $current_sequence = undef;
  my $sequence_fragment = undef;

  #
  # Dont forget to keep a set of arrays for writing as well!
  #
  my $array_name = undef; 
  my $probeset = undef;
  my $probe_name = undef;
  my $current_affy_array = undef;
  my $existing_probe = undef;
  

  while(<PROBES>){

    chomp;
    if(/^>probe:(\S+):(\S+):(\S+:\S+;).*$/){

      if($current_sequence){
        if(!$current_affy_array){
          throw ("Have sequence $current_sequence but no current array !\n");
        }

        #
        #You can only place the probe into the hash at this point, because
        #it's only at this point that you know you have the full sequence

        $existing_probe = $probes_by_sequence{$probeset.":---:".$current_sequence};

        if(!$existing_probe){
          $existing_probe = 
            $self->create_new_probe(
              $current_affy_array, 
              $probeset, 
              $probe_name, 
              $current_sequence
            );

          $probes_by_sequence{$probeset.":---:".$current_sequence} = $existing_probe;
        }else{
          $self->add_array_to_existing_probe(
            $existing_probe, 
            $current_affy_array, 
            $probeset,
            $probe_name 
          );
        }

        $current_sequence = undef;
      }

      $array_name = $1;
      $probeset = $2;
      $probe_name = $3;

      $current_affy_array = $affy_arrays{$array_name};
      if(!$current_affy_array){
        $current_affy_array = $self->create_new_array($array_name);
        $affy_arrays{$array_name} = $current_affy_array;
      }
    }else{
      $sequence_fragment = $_;
      if($current_sequence){
        $current_sequence = $current_sequence . $sequence_fragment;
      }else{
        $current_sequence = $sequence_fragment;
      }
    }
  }

  $affy_arrays{$array_name} = $current_affy_array;
  $existing_probe = $probes_by_sequence{$probeset.":---:".$current_sequence};

  if(!$existing_probe){
    $existing_probe = 
      $self->create_new_probe(
        $current_affy_array, 
        $probeset, 
        $probe_name, 
        $current_sequence
      );

    $probes_by_sequence{$probeset.":---:".$current_sequence} = $existing_probe;
  }else{
    $self->add_array_to_existing_probe(
      $existing_probe, 
      $current_affy_array, 
      $probeset,
      $probe_name 
    );
  }  
  $self->affy_probes(\%probes_by_sequence);
  $self->affy_arrays(\%affy_arrays);
}

sub add_array_to_existing_probe{
  my ($self, $probe, $array, $probeset, $probename) = @_;
  my @all_probenames = @{$probe->get_all_probenames};

      
  if(!($probe->probeset eq $probeset)){
    throw (
      "Inconsistency: have found a probe ".$array->name.":$probeset:$probename with ".
      "identical sequence but different probeset to another one: ".$probe->probeset."\n"
    );
  }
  $probe->add_Array_probename($array, $probename);
}

sub create_new_probe {
  my($self, $affy_array, $probeset, $probe_name, $current_sequence) = @_;

  my $affy_probe =
     new Bio::EnsEMBL::AffyProbe(
       -probeset => $probeset,
       -name => $probe_name,
       -array => $affy_array
     );

  return $affy_probe;
}

sub create_new_array {
  my($self, $array_name) = @_;
  my $affy_array =
     new Bio::EnsEMBL::AffyArray(
       -name => $array_name,
       -setsize => 0,
     );
}
############################################################

sub write_output {
  my ( $self, @output ) = @_;

  my $outdb        = $self->get_output_db;
  my $outfile = $self->non_redundant_probe_seqs;
  open (OUTFILE, ">".$outfile);
  my $affy_array_adaptor = $outdb->get_AffyArrayAdaptor;
  my $affy_probe_adaptor = $outdb->get_AffyProbeAdaptor;

  foreach my $affy_array (values %{$self->affy_arrays}){
    my $existing_affy_array = $affy_array_adaptor->fetch_by_name($affy_array->name);
    if(!$affy_array->dbID){
      eval{ $affy_array_adaptor->store($affy_array) };
      if ($@) {
        $self->throw("Unable to store affy array!\n $@");
      }
    }
  }

  foreach my $probeset_sequence_key (keys %{$self->affy_probes}){
    my $affy_probe = $self->affy_probes->{$probeset_sequence_key};
    if(!$affy_probe->dbID){
      eval{ $affy_probe_adaptor->store($affy_probe) };
      if ($@) {
        $self->throw("Unable to store affy probe!\n $@");
      }
    }
    
    #Now that you've stored it, you have the probe's dbID, so you can write out the NEW
    #fasta file, with that dbID as a header...
    $probeset_sequence_key =~ /.*:---:(.*)/;
    my $sequence = $1;
    print OUTFILE ">".$affy_probe->dbID."\n";
    print OUTFILE $sequence."\n";
  }
  close(OUTFILE) or throw("Failed top close ".$outfile);
}

############################################################

sub clean_affy_features {
  my ( $self, @affy_features ) = @_;

  my $slice_adaptor = $self->db->get_SliceAdaptor;

  my %genome_slices;

  foreach my $affy_feature (@affy_features) {

    $affy_feature->analysis( $self->analysis );

    # get the slice based on the seqname stamped on in the runnable
    my $slice_id = $affy_feature->seqname;

    if ( not exists $genome_slices{$slice_id} ) {
      # assumes genome seqs were named in the Ensembl API Slice naming
      # convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
      $genome_slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
    }
    my $slice = $genome_slices{$slice_id};

    $affy_feature->slice($slice);

    #This is a hack: a probe might sit on two arrays, but I don't 
    # know how to resolve this ambiguity now, or whether it's significant
    my $array_name = $affy_feature->probe->get_all_AffyArrays->[0]->name;
    my $probe_name = $array_name.":".$affy_feature->probeset.":".$affy_feature->probe->get_all_probenames->[0];
    my $real_probe = $self->affy_probes->{$probe_name};
    
    if(!$real_probe){
      throw "Inconsistency! I can't find an affy probe corresponding to $probe_name\n";
    }

    #print "setting real probe: ".$real_probe->get_all_AffyArrays->[0]->name.":".
    $affy_feature->probe($real_probe);
  }
  
}

############################################################

sub populate_affy_arrays_and_probes{
  my($self, @args) = @_;
  
  my $query_file = $self->query_file;

  my $outdb        = $self->get_output_db;
  my $affy_array_adaptor = $outdb->get_AffyArrayAdaptor;
  my $affy_probe_adaptor = $outdb->get_AffyProbeAdaptor;
  
  #
  #make a pass through each fasta header in the query file,
  #and sift out the array and probe id.
  #If they don't exist in our %affy_array and %affy_probe hashes,
  # store an object for them.
  open(PROBES, "<".$self->query_file);

  while(<PROBES>){
    chomp;
    /^>probe:(\S+):(\S+):(\S+:\S+;).*$/;
    my $array_name = $1;
    my $probeset = $2;
    my $probe_name = $3;
    
    if (!$array_name){
      throw "array name could not be deduced from: ".$_."\n";
    }
    
    if( !$probeset){
      throw "probeset could not be deduced from: ".$_."\n";
    }
    
    if(!$probe_name){
      throw "probename could not be deduced from: ".$_."\n";
    }
    
    my $affy_array = $affy_array_adaptor->fetch_by_name($array_name);
    if(!$affy_array){
      $affy_array = 
        new Bio::EnsEMBL::AffyArray(
          -name => $array_name
        );
    }
    
    $self->affy_arrays->{$array_name} = $affy_array;

    my $affy_probe = 
      $affy_probe_adaptor->fetch_by_array_probeset_probe(
        $array_name,
        $probeset,
        $probe_name
      );

    if(!$affy_probe){
      $affy_probe = 
        new Bio::EnsEMBL::AffyProbe(
          -probeset => $probeset,
          -name => $probe_name,
          -array => $affy_array
        );
    }

    $self->affy_probes->{$array_name.":".$probeset.":".$probe_name} = $affy_probe;
  }
  
}

############################################################

sub get_output_db {
  my ($self) = @_;

  my $outdb;

  if ( $self->OUTDB ) {
    $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $self->OUTDB }, -dnadb => $self->db );
  } else {
    $outdb = $self->db;
  }

  return $outdb;
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

sub affy_arrays {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_affy_arrays'} = $value;
  }

  if ( exists( $self->{'_affy_arrays'} ) ) {
    return $self->{'_affy_arrays'};
  } else {
    return undef;
  }
}

#############################################################

sub affy_probes {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_affy_probes'} = $value;
  }

  if ( exists( $self->{'_affy_probes'} ) ) {
    return $self->{'_affy_probes'};
  } else {
    return undef;
  }
}


#############################################################
# Declare and set up config variables
#############################################################

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($AFFY_CONFIG);

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
  foreach my $config_var (
    qw(
      QUERYSEQS
      NON_REDUNDANT_PROBE_SEQS
      OUTDB
    )
  ){
    if ( not defined $self->$config_var ){
      throw("You must define $config_var in config for logic '$logic'");
    }
  }

  # output db does not have to be defined, but if it is, it should be a hash
  if ( $self->OUTDB and ref( $self->OUTDB ) ne "HASH" ) {
    throw("OUTDB in config for '$logic' must be a hash ref of db connection pars.");
  }
}

sub QUERYSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_QUERYSEQS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_QUERYSEQS'} ) ) {
    return $self->{'_CONFIG_QUERYSEQS'};
  } else {
    return undef;
  }
}

sub NON_REDUNDANT_PROBE_SEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_NON_REDUNDANT_PROBE_SEQS'} = $value;
  }

  if ( exists( $self->{'_NON_REDUNDANT_PROBE_SEQS'} ) ) {
    return $self->{'_NON_REDUNDANT_PROBE_SEQS'};
  } else {
    return undef;
  }
}

sub OUTDB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_OUTDB'} = $value;
  }

  if ( exists( $self->{'_CONFIG_OUTDB'} ) ) {
    return $self->{'_CONFIG_OUTDB'};
  } else {
    return undef;
  }
}

###############################################
###     end of config
###############################################

1;
