#
# package Bio::EnsEMBL::Analysis::Config::ImportArrays
# 
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Config::ImportArrays

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis:Config::ImportAarrys;

=head1 DESCRIPTION

This contains the configuration for importing arrays from flat files.
It is entirely dependant on the arrays.env environment which can be used 
to set up and run the pipeline in an easy and interactive way. This contains 
all possible configurations which will then be set dynamically by the RunnableDB
for each instance using the input_id as a key into a separate ImportArrays.conf 
file, listed here as ARRAY_FORMAT_FILE.


The layout of the configuration is a set of hashes,
each one keyed by logic name. There is also a DEFAULT hash,
which is used as the default for all logic names (this
was the configuration pattern stolen from Exonerate2Genes,
although in this case it's very unlikely you will need to have
different configs by logic name).

=head1 CONTACT

=cut


package Bio::EnsEMBL::Analysis::Config::ImportArrays;

use strict;
use vars qw( %Config );

# Hash containing config info
# -- one hashnode per logic name, with a 'DEFAULT' logic name provided
#

%Config = 
  (
   ARRAY_CONFIG => 
   {
	DEFAULT => 
	{
	 #These are now defined dynamically or via the ImportArrays.conf file
	 # All input probes must be kept in one huge (possibly redundant) fasta file
	 #QUERYSEQS            => $ENV{'RAW_FASTA'},
	 # The output of this module writes a set of affy probes into the OUTDB.affy_probe table,
	 # and also writes the nonredundant probes into this fasta file,
	 # with the fasta headers keyed with the affy probes' internal id. 
	 #NON_REDUNDANT_PROBE_SEQS => $ENV{'NR_FASTA'},
	 
	 # DB containing all affy_arrays, affy_probes and (next step) affy_features
	 OUTDB => {
			   -dbname => $ENV{'DB_NAME'},
			   -host   => $ENV{'DB_HOST'},
			   -port   => $ENV{'DB_PORT'},
			   -user   => $ENV{'DB_USER'},
			   -pass   => $ENV{'DB_PASS'},
			  },

	 #Optional, must define if dnadb is not on ensembldb
	 DNADB => {
			   -dbname => $ENV{'DNADB_NAME'},
			   -host   => $ENV{'DNADB_HOST'},
			   -port   => $ENV{'DNADB_PORT'},
			   -user   => $ENV{'DNADB_USER'},
			   -pass   => $ENV{'DNADB_USER'},
			  },
	 
	 #The ImportArray.conf file for this instance of the RunnableDB
	 ARRAY_FORMAT_FILE => $ENV{'ARRAY_FORMAT_FILE'},

	 #Used for building the format specific NR fasta file
	 DB_HOME           => $ENV{'DB_HOME'},


	 #This defines how to parse the file headers
	 IIDREGEXP => {
				   AFFY      => '^>probe:(\S+):(\S+):(\S+:\S+;).*$',
				   #AFFY_ST => 
				   #ILLUMINA  =>

				   #if(/^>probe:([^:]+):([^:]+):([0-9:]+;).*$/){#hacked affy,  one o
				   #  if(/^>probe:(\S+):(\S+).*$/){   #NATH hack to get non probe_set arrays to work   
				  },

	 #We also need a has to define the input field order
	 #This will be used to set the relevant hash values
	 IFIELDORDER => {
					 AFFY => {
							  -name       => 2,
							  -array_chip => 0,
							  -probeset   => 1,#This is not using a class yet? How will this impact on the collapsing?
							 }

					},

	 #This is used to store Arrays
	 ARRAY_PARAMS => {
					  'MG-U74Cv2' => {
									  -name => 'MG-U74Cv2',
									  -vendor => 'AFFY',
									  #-setsize => undef,
									  -format  => 'EXPRESSION',#? UTR?
									  -type    => 'OLIGO',
									  #-description => '',
									  
									 },

					 },

	 #Used to call the correct run method
	 INPUT_FORMAT => {
					  AFFY        => 'FASTA',
					  AFFY_ST     => 'FASTA',
					  #ILLUMINA
					  #ILLUMINA_V1
					  #ILLUMINA_V2
					  #CODELINK
					  #AGILENT
					  #?
					 },
	 
	},


	IMPORTARRAYS => 
	{
	 OUTDB => {
			   -dbname => $ENV{'DB_NAME'},
			   -host   => $ENV{'DB_HOST'},
			   -port   => $ENV{'DB_PORT'},
			   -user   => $ENV{'DB_USER'},
			   -pass   => $ENV{'DB_PASS'},
			  },
	 
	 DNADB => {
			   -dbname => $ENV{'DNADB_NAME'},
			   -host   => $ENV{'DNADB_HOST'},
			   -port   => $ENV{'DNADB_PORT'},
			   -user   => $ENV{'DNADB_USER'},
			   -pass   => $ENV{'DNADB_USER'},
			  },

	 ARRAY_FORMAT_FILE => $ENV{'ARRAY_FORMAT_FILE'},
	 DB_HOME           => $ENV{'DB_HOME'},

	 IIDREGEXP => {
				   AFFY      => '^>probe:(\S+):(\S+):(\S+:\S+;).*$',
				   #AFFY_ST => 
				   #ILLUMINA  =>

				   #if(/^>probe:([^:]+):([^:]+):([0-9:]+;).*$/){#hacked affy,  one o
				   #  if(/^>probe:(\S+):(\S+).*$/){   #NATH hack to get non probe_set arrays to work   
				  },

	 IFIELDORDER => {
					 AFFY => {
							  -name       => 2,
							  -array_chip => 0,
							  -probeset   => 1,
							 }

					},

	 ARRAY_PARAMS => {
					  'MG-U74Cv2' => {
									  -name => 'MG-U74Cv2',
									  -vendor => 'AFFY',
									  #-setsize => undef,
									  -format  => 'EXPRESSION',#? UTR?
									  -type    => 'OLIGO',
									  #-description => '',
									  
									 },

					 },

	 INPUT_FORMAT => {
					  AFFY        => 'FASTA',
					  AFFY_ST     => 'FASTA',
					  #ILLUMINA
					  #ILLUMINA_V1
					  #ILLUMINA_V2
					  #CODELINK
					  #AGILENT
					  #?
					 },

	},
   }
  );

sub import {
  my ($callpack) = caller(0); # Name of the calling package
  my $pack = shift; # Need to move package off @_

  # Get list of variables supplied, or else everything
  my @vars = @_ ? @_ : keys( %Config );
  return unless @vars;
  
  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $Config{$_} ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$Config{ $_ };
	} else {
	    die "Error: Config: $_ not known\n";
	}
    }
}

1;
