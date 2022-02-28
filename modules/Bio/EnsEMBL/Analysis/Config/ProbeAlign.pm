
=head1 NAME

Bio::EnsEMBL::Analysis::Config::ProbeAlign

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::ProbeAlign;

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


package Bio::EnsEMBL::Analysis::Config::ProbeAlign;

use warnings ;
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
			   -dbname          => $ENV{'DNADB_NAME'},
			   -host            => $ENV{'DNADB_HOST'},
			   -port            => $ENV{'DNADB_PORT'},
			   -user            => $ENV{'DNADB_USER'},
			   -pass            => $ENV{'DNADB_PASS'},
			   -species         => $ENV{'SPECIES'},
			   -multispecies_db => $ENV{'DNADB_MULTISPECIES_DB'},
			   -species_id      => $ENV{'DNADB_SPECIES_ID'}
			  },
	 
	 OUTDB => {
			   -dbname          => $ENV{'DB_NAME'},
			   -host            => $ENV{'DB_HOST'},
			   -port            => $ENV{'DB_PORT'},
			   -user            => $ENV{'DB_USER'},
			   -pass            => $ENV{'DB_PASS'},
			   -species         => $ENV{'SPECIES'},#required for auto generation fo DNADB
			   -multispecies_db => $ENV{'MULTISPECIES_DB'},
			   -species_id      => $ENV{'SPECIES_ID'}
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

	#Remember: Exonerate scoring is +5 for match -4 for mismatch


	AFFY_UTR_PROBEALIGN => 
	{
	 MAX_MISMATCHES       => 1,
	 TARGETSEQS           => $ENV{'GENOMICSEQS'},
	 QUERYSEQS            => $ENV{'WORK_DIR'}.'/arrays_nr.AFFY_UTR.fasta',
	 #at least 25 mers allowing 1bp mismatch
	 #this will still work in the worst case where the mismatch is at the centre of a 25bp probe
	 
	 OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 12 --seedrepeat 2 --dnahspthreshold 118 --dnawordlimit 0',
	 
	 #Perfect matches only for Jing
	 #MAX_MISMATCHES       => 0,
	 #OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 25 --seedrepeat 1 --dnahspthreshold 125 --dnawordlimit 0',
	 

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


	#ILLUMINA_WG are 50mers. These settings allow for at least 1bp mismatch
	ILLUMINA_WG_PROBEALIGN => 
	{
	 TARGETSEQS         => $ENV{'GENOMICSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.ILLUMINA_WG.fasta',
	 #50mers
	 OPTIONS => ' --bestn 101 --dnahspthreshold 246 --fsmmemory 256 --dnawordlen 25 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1, 
	},		  

	ILLUMINA_WG_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.ILLUMINA_WG.fasta',
	 #50mers
	 OPTIONS => ' --bestn 101 --dnahspthreshold 246 --fsmmemory 256 --dnawordlen 25 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1, 
	},		  
	

	#ILLUMINA_INFINIUM are 50mers. These settings allow for at least 1bp mismatch
	ILLUMINA_INFINIUM_PROBEALIGN => 
	{
	 TARGETSEQS         => $ENV{'GENOMICSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.ILLUMINA_INFINIUM.fasta',
	 #50mers
	 OPTIONS => ' --bestn 101 --dnahspthreshold 246 --fsmmemory 256 --dnawordlen 25 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,#This may need changing?
	 MAX_MISMATCHES => 1, 
	},		  



	#CODELINK are 30mers
	CODELINK_PROBEALIGN => 
	{
	 TARGETSEQS         => $ENV{'GENOMICSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.CODELINK.fasta',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 141 --fsmmemory 256 --dnawordlen 15 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1, 
	},		  

	CODELINK_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.CODELINK.fasta',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 141 --fsmmemory 256 --dnawordlen 15 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1, 
	},	


	#PHALANX are 60mers
	PHALANX_PROBEALIGN => 
	{
	 TARGETSEQS         => $ENV{'GENOMICSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.PHALANX.fasta',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 30 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1, 
	},		  

	PHALANX_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.PHALANX.fasta',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 30 --seedrepeat 2 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1, 
	},	


	AFFY_UTR_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.AFFY_UTR.fasta',
	 #25 mers 
	 OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 12 --seedrepeat 2 --dnahspthreshold 118 --dnawordlimit 0',
	 MAX_MISMATCHES       => 1,
	 #Perfect matches only for Jing
	 #MAX_MISMATCHES       => 0,
	 #OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 25 --seedrepeat 1 --dnahspthreshold 125 --dnawordlimit 0',
	 

	 #OPTIONS => ' --bestn 101 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 14 --dnawordlimit 11 ',
	 #HIT_SATURATION_LEVEL => 100,#I don't think we want this for the transcript mappings
	 #Defaults to 100 anyway, but not used
	 #FILTER_METHOD => 'filter_mismatches',#Would need to add another method to Runnable::Exonerate
	 #MAX_MISMATCHES => 1,
	},


	#Essentially same as AFFY but with different NR_FASTA
	AFFY_ST_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         =>  $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS      => $ENV{'WORK_DIR'}.'/arrays_nr.AFFY_ST.fasta',
	 #25 mers 
	 OPTIONS             => ' --bestn 101 --fsmmemory 256 --dnawordlen 12 --seedrepeat 2 --dnahspthreshold 118 --dnawordlimit 0',
	 #OPTIONS => ' --bestn 101 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 14 --dnawordlimit 11 ',
	 #HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 1,
	},




	#AGILENT 60 mers
	#Min length 45bp
	AGILENT_PROBEALIGN => 
	{
	 TARGETSEQS         => $ENV{'GENOMICSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.AGILENT.fasta',
	 #OPTIONS => ' --bestn 101 --dnahspthreshold 216 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ', #ORIG
	 HIT_SATURATION_LEVEL => 100,
	 #MAX_MISMATCHES => 1,

	 #Danio Zv7 params
	 #Do we need a way of setting these in env so we don't have to edit here?
	 #Can we do a $ENV{'PARAM'} || ref self default?
	 #These would have to be set before this hash by reading ini into hash?
	 OPTIONS => ' --bestn 101 --dnahspthreshold 216 --fsmmemory 256 --dnawordlen 15 --seedrepeat 4 --dnawordlimit 0 ',
	 MAX_MISMATCHES => 3,

	},		  

	AGILENT_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.AGILENT.fasta',
	 #OPTIONS => ' --bestn 101 --dnahspthreshold 216 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ', #ORIG
	 HIT_SATURATION_LEVEL => 100,
	 #MAX_MISMATCHES => 1,

	 #Danio Zv7 params
	 #Do we need a way of setting these in env so we don't have to edit here?
	 #Can we do a $ENV{'PARAM'} || ref self default?
	 #These would have to be set before this hash by reading ini into hash?
	 OPTIONS => ' --bestn 101 --dnahspthreshold 216 --fsmmemory 256 --dnawordlen 15 --seedrepeat 4 --dnawordlimit 0 ',
	 MAX_MISMATCHES => 3,
	 

	},	


        CATMA_PROBEALIGN =>
        {
         TARGETSEQS         => $ENV{'GENOMICSEQS'},
         QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.CATMA.fasta',
         OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 59 --seedrepeat 2 --dnawordlimit 0 ', #ORIG
         HIT_SATURATION_LEVEL => 100,
         MAX_MISMATCHES => 1,


        },

        CATMA_PROBETRANSCRIPTALIGN =>
        {
         TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
         QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.CATMA.fasta',
         OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 59 --seedrepeat 2 --dnawordlimit 0 ', #ORIG
         HIT_SATURATION_LEVEL => 100,
         MAX_MISMATCHES => 1,

        },
	
     NSF_PROBEALIGN =>
        {
         TARGETSEQS         => $ENV{'GENOMICSEQS'},
         QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.NSF.fasta',
         OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ', #ORIG
         HIT_SATURATION_LEVEL => 100,
         MAX_MISMATCHES => 1,


        },

        NSF_PROBETRANSCRIPTALIGN =>
        {
         TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
         QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.NSF.fasta',
         OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ', #ORIG
         HIT_SATURATION_LEVEL => 100,
         MAX_MISMATCHES => 1,

        },




	#LEIDEN 50 mers
	#Actually some are 50 some are 60.
	#3 mismatches for Zv7 due to low quality assembly
	LEIDEN_PROBEALIGN => 
	{
	 TARGETSEQS         => $ENV{'GENOMICSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.LEIDEN.fasta',
	 #OPTIONS => ' --bestn 101 --dnahspthreshold 241 --fsmmemory 256 --dnawordlen 25 --seedrepeat 2 --dnawordlimit 0 ',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 223 --fsmmemory 256 --dnawordlen 13 --seedrepeat 4 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 3, 
	},		  

	LEIDEN_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.LEIDEN.fasta',
	 #OPTIONS => ' --bestn 101 --dnahspthreshold 241 --fsmmemory 256 --dnawordlen 25 --seedrepeat 2 --dnawordlimit 0 ',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 223 --fsmmemory 256 --dnawordlen 13 --seedrepeat 4 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES =>3, 
	},

	
	#STEMPLE_LAB_SANGER 65 mers
	#3 mismatches  due to low quality assembly
	

	STEMPLE_LAB_SANGER_PROBEALIGN => 
	{
	 TARGETSEQS         => $ENV{'GENOMICSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.STEMPLE_LAB_SANGER.fasta',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 298 --fsmmemory 256 --dnawordlen 17 --seedrepeat 4 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES => 3, 
	},		  

	STEMPLE_LAB_SANGER_PROBETRANSCRIPTALIGN => 
	{
	 TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
	 QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.STEMPLE_LAB_SANGER.fasta',
	 OPTIONS => ' --bestn 101 --dnahspthreshold 298 --fsmmemory 256 --dnawordlen 17 --seedrepeat 4 --dnawordlimit 0 ',
	 HIT_SATURATION_LEVEL => 100,
	 MAX_MISMATCHES =>3, 

	},
    

    #WUSTL Custom arrays (only used for C.elegans AFAIK) 60 mers
    WUSTL_PROBEALIGN => 
    {
     TARGETSEQS         => $ENV{'GENOMICSEQS'},
     QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.WUSTL.fasta',
     OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ', 
     HIT_SATURATION_LEVEL => 100,
     MAX_MISMATCHES => 1,
     
    },                
    
    WUSTL_PROBETRANSCRIPTALIGN => 
    {
     TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
     QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.WUSTL.fasta',
     OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ',
     HIT_SATURATION_LEVEL => 100,
     MAX_MISMATCHES => 1,    
     
    },


    #UCSF Custom arrays (only used for C.elegans) 50-70 mers
    UCSF_PROBEALIGN => 
    {
     TARGETSEQS         => $ENV{'GENOMICSEQS'},
     QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.UCSF.fasta',
     OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ', 
     HIT_SATURATION_LEVEL => 100,
     MAX_MISMATCHES => 1,
    },                
    
    UCSF_PROBETRANSCRIPTALIGN => 
    {
     TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
     QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.UCSF.fasta',
     OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ',
     HIT_SATURATION_LEVEL => 100,
     MAX_MISMATCHES => 1,    
    },

    #SLRi Custom arrays (only used for C.elegans) 50-70 mers
    SLRI_PROBEALIGN => 
    {
     TARGETSEQS         => $ENV{'GENOMICSEQS'},
     QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.SLRI.fasta',
     OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ', 
     HIT_SATURATION_LEVEL => 100,
     MAX_MISMATCHES => 1,
    },                
    
    SLRI_PROBETRANSCRIPTALIGN => 
    {
     TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
     QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.SLRI.fasta',
     OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ',
     HIT_SATURATION_LEVEL => 100,
     MAX_MISMATCHES => 1,    
    },   


    # NIMBLEGen modENCODE arrays (only used for C.elegans) 60 mers
    NIMBLEGEN_MODENCODE_PROBEALIGN => 
    {
     TARGETSEQS         => $ENV{'GENOMICSEQS'},
     QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.NIMBLEGEN_MODENCODE.fasta',
     OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ', 
     HIT_SATURATION_LEVEL => 100,
     MAX_MISMATCHES => 1,
     
    },                
    
    NIMBLEGEN_MODENCODE_PROBETRANSCRIPTALIGN => 
    {
     TARGETSEQS         => $ENV{'TRANSCRIPTSEQS'},
     QUERYSEQS          => $ENV{'WORK_DIR'}.'/arrays_nr.NIMBLEGEN_MODENCODE.fasta',
     OPTIONS => ' --bestn 101 --dnahspthreshold 291 --fsmmemory 256 --dnawordlen 22 --seedrepeat 2 --dnawordlimit 0 ',
     HIT_SATURATION_LEVEL => 100,
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
