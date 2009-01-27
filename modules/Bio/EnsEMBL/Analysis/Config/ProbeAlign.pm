#
# package Bio::EnsEMBL::Analysis::Config::Funcgen::ProbeAlign?
# 
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Config::Funcgen::ProbeAlign?

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::Funcgen::ProbeAlign;

=head1 DESCRIPTION

This contains the configuration for step 2 of the 
process which maps probes to a Genome. This step
is an alignment of probes (dna) against a genome (dna)
using exonerate. So this config looks very similar to that
of any other exonerate-driving config.

The layout of the configuration is a set of hashes,
each one keyed by logic name. There is also a DEFAULT hash,
which is used as the default for all logic names.

There are genomic and transcript based logic names and config
hashes for each discrete format of array.


=head1 CONTACT

=cut


package Bio::EnsEMBL::Analysis::Config::ProbeAlign;

use strict;
use vars qw( %Config );



%Config =
  (

   #This entire hash is exported as the global $PROBE_CONFIG var
   #each key will be exported as $PROBE_CONFIG->{'_CONFIG_'.$key}
   #Dependant on logic name of RunnableDB


   PROBE_CONFIG => 
   {
	DEFAULT => 
	{
	 
	 # path to softmasked, dusted genomic sequence or transcript seqs on the farm 
	 # 	 #'/data/blastdb/Ensembl/Rmacaque/MMUL_2/genome/softmasked_dusted.fa',
	 #/data/blastdb/Ensembl/Human/NCBI35/softmasked_dusted/', #allowed to be a dir.
	 TARGETSEQS         =>  $ENV{'GENOMICSEQS'},#or $ENV{'TRANSCRIPTSEQS'}
	 QUERYTYPE           => 'dna',
	 
	 # must be a single file containing all (non-redundant) probes indexed by affy_probe_id
	 # QUERYSEQS refers to the value of the parameter NON_REDUNDANT_PROBE_SEQS
	 # in the config-file ensembl-analysis/Config/CollapseAffyProbes.pm	 
	 #QUERYSEQS           => $ENV{'NR_FASTA'},
	 #Removed this now as we want to run different analyses at the same time so we have to hardcode below
	 

	 
	 # must supply one, since the queryseqs MUST be a single file
	 #InputIDREGEXP this is used to infer chunk number from the headers of a single fasta file
	 #Therefore we cannot have mixed type in the same file, must be in a different array set
	 #If not related, or reformated prior to Import if they are related
	 IIDREGEXP           => '(\d+):(\d+)',
	 
	 #DNADB is not essential, but we need this if we are going to define a DNADB not on ensembldb 
	 #e.g. new release on staging
	 #Add species and group here?
	 DNADB => {
			   -dbname => $ENV{'DNADB_NAME'},
			   -host   => $ENV{'DNADB_HOST'},
			   -port   => $ENV{'DNADB_PORT'},
			   -user   => $ENV{'DNADB_USER'},
			   -pass   => $ENV{'DNADB_PASS'},
			   
			  },
	 
	 OUTDB => {
			   -dbname  => $ENV{'DB_NAME'},
			   -host    => $ENV{'DB_HOST'},
			   -port    => $ENV{'DB_PORT'},
			   -user    => $ENV{'DB_USER'},
			   -pass    => $ENV{'DB_PASS'},
			   -species => $ENV{'SPECIES'},#required for auto generation fo DNADB
			  },
	 
	 
	 #25 mers
	 OPTIONS             => ' --bestn 100 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 25 --dnawordlimit 11 ',
	 #50 mers
	 #OPTIONS             => ' --bestn 100 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 50 --dnawordthreshold 11 ',				  			      
	 # if the number of hits reaches or exceeds the figure below, we reject 
	 # all hits to that probe
	 HIT_SATURATION_LEVEL => 100,#2 for unique tiling arrays mappings

	 MAX_MISMATCHES => 1,#No way to dynamically test method prerequisite config vars without setting a FILTER_METHOD hash?
	 #This would not be bullet proof a the hash could then be edited
	 
	 #Would need to add another method to Runnable::ExonerateProbe to change this
	 #FILTER_METHOD => 'filter_mismatches',
	 	 
	 #Or we can pass a code ref here to allow easy extension without editing Runnables etc.
	 #Can't name code ref subs!
	 #This is used in ExonerateProbe.pm
	 #e.g.
	 #FILTER_METHOD => sub {
	 #  my ($self, $query_match_length, $q_length, $score) = @_;
	 #  my $mismatch;
	 #  my $full_score = $q_length * 5;
	   
	 #  if($query_match_length == $q_length){
		 
	#	 if($score == $full_score){
	#	   $mismatch = 0;
	#	 }
	#   }
	   
	 #  if(! defined $mismatch){
		 
	#	 my $max_mismatch = $self->allowed_mismatches;
		 
	#	 for my $i(1..$max_mismatch){
		   
	#	   my $mismatch_length = $q_length - $i;
	#	   my $mismatch_score = $mismatch_length * 5;
		   
		#   if($query_match_length == $q_length){
			 
		#	 if ($score == ($full_score - ($i*9))) {
		#	   $mismatch = $i;
		#	 }
		#   }
		#   elsif($query_match_length == $mismatch_length){
		#	 $mismatch = $i if ($score == $mismatch_score);
		#   }
		# }
	  # } 
	  # return $mismatch;
	 #},

	},#end of DEFAULT


	#IIDREGEXP, DNADB, OUTDB and QUERYSEQS and QUERYTYPE should be same for all these configs


	#Need to add  ILLUMINA_PROBE_ALIGN,  ILLUMINA_PROBE_TRANSCRIPT_ALIGN, CODELINK, AGILENT etc
	 
	 
	#There is no point in using a % threshold here as bestn value will most likely cause only high quality
	#hits to be returned.  However there is a possiblity with longer sequences that we may get a duff alignment.
	#We could set -percent to a conservative 95%, but this entirely depends on the length and number of allowed
	#mismatches.

	#Define QUERYSEQS here instead of in Runnable to prevent convolution 
	#of RunnableDB and environment, 
	#i.e. we can still run ProbeAlign so long as we change this config
	#The only downfall is lack of validation of QUERYSEQS files
	#WARNING CHECK YOUR QUERYSEQS!

	AFFY_UTR_PROBEALIGN => 
	{
	 MAX_MISMATCHES       => 1,
	 TARGETSEQS           => $ENV{'GENOMICSEQS'},
	 QUERYSEQS            => $ENV{'WORK_DIR'}.'/arrays_nr.AFFY_UTR.fasta',
	 #at least 25 mers allowing 1bp mismatch
	 #this will still work in the worst case where the mismatch is at the centre of a 25bp probe
	 
	 OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 12 --seedrepeat 2 --dnahspthreshold 118 --dnawordlimit 0',
	 

	 #OPTIONS              => ' --bestn 101 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 14 --dnawordlimit 11 ',
	 HIT_SATURATION_LEVEL => 100,
	},
	
	#Essentially same as AFFY but with different NR_FASTA
	AFFY_ST_PROBEALIGN => 
	{
	 TARGETSEQS           => $ENV{'GENOMICSEQS'},
	 QUERYSEQS            => $ENV{'WORK_DIR'}.'/arrays_nr.AFFY_ST.fasta',
	 #25 mers 
	 #OPTIONS              => ' --bestn 101 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 14 --dnawordlimit 11 ',
	 OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 12 --seedrepeat 2 --dnahspthreshold 118 --dnawordlimit 0',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES       => 1,
	},

	
	NIMBLEGEN_PROBEALIGN => 
	{
	 TARGETSEQS         =>  $ENV{'GENOMICSEQS'},
	 
	 #Need to define this dynamically based on oligo length (40-60mers)
	 #50mers
	 #Can we up the dnaword limit 10 50 here if we only want unique matches?
	 OPTIONS => ' --bestn 2 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 50 --dnawordlimit 11 ',
	 HIT_SATURATION_LEVEL => 2, #We only want unique mappings for tiling probes
	 MAX_MISMATCHES => 0, #Unique mappings for tiling probes
	},		  


	#ILLUMINA_WG are 51mers. These settings allow for at least 1bp mismatch
	ILLUMINA_WG_PROBEALIGN => 
	{
	 TARGETSEQS         => $ENV{'GENOMICSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.ILLUMINA_WG.fasta',
	 #Need to define this dynamically based on oligo length (40-60mers)
	 #50mers
	 OPTIONS => ' --bestn 101 --dnahspthreshold 246 --fsmmemory 256 --dnawordlen 25 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1, #Unique mappings for tiling probes
	},		  

	ILLUMINA_WG_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.ILLUMINA_WG.fasta',
	 #Need to define this dynamically based on oligo length (40-60mers)
	 #50mers
	 OPTIONS => ' --bestn 101 --dnahspthreshold 246 --fsmmemory 256 --dnawordlen 25 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1, #Unique mappings for tiling probes
	},		  
	

	AFFY_UTR_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.AFFY_UTR.fasta',
	 #25 mers 
	 OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 12 --seedrepeat 2 --dnahspthreshold 118 --dnawordlimit 0',
	 #OPTIONS => ' --bestn 101 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 14 --dnawordlimit 11 ',
	 #HIT_SATURATION_LEVEL => 100,#I don't think we want this for the transcript mappings
	 #Defaults to 100 anyway, but not used
	 #FILTER_METHOD => 'filter_mismatches',#Would need to add another method to Runnable::Exonerate
	 MAX_MISMATCHES => 1,
	},


	#Essentially same as AFFY but with different NR_FASTA
	AFFY_ST_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         =>  $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS      => $ENV{'WORK_DIR'}.'/arrays_nr.AFFY_ST.fasta',
	 #25 mers 
	 #OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 12 --seedrepeat 2 --dnahspthreshold 118 --dnawordlimit 0',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 14 --dnawordlimit 11 ',
	 #HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1,
	}, 
   }
  );



sub import {
  my ($callpack) = caller(0);	# Name of the calling package
  my $pack = shift;				# Need to move package off @_

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
