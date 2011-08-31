#!/usr/local/ensembl/bin/perl

use setup_rnaseq_pipeline_config;
use vars qw(%Config);
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use strict;
use Getopt::Long;


my $verbose;
my $file;
my $dbname = $RNASEQCONFIG->{DB};
my $analysisconfigdir = $RNASEQCONFIG->{ANALYSISCONFIG_DIR};
my $pipelineconfigdir =  $RNASEQCONFIG->{PIPELINECONFIG_DIR};
my $delimiter = $RNASEQCONFIG->{DELIMITER};
my $summaryfile  = $RNASEQCONFIG->{SUMMARY};
my $input_dir = $RNASEQCONFIG->{INPUT_DIR};
my $output_dir = $RNASEQCONFIG->{OUTPUT_DIR};
my $all_paired = $RNASEQCONFIG->{ALL_PAIRED};
my $regex = $RNASEQCONFIG->{PAIRING_REGEX};
my $write;

my $usage = "perl setup_rnaseq_pipeline.pl
-verbose    $verbose,
-write      $write, write config to directory
Need to fill in the config in the setup_rnaseq_pipeline_config.pm module.
";

$| = 1;
&GetOptions(
	    'verbose!'     => \$verbose,
	    'write!'       => \$write,
	   );

die($usage) unless ($dbname && $analysisconfigdir && $pipelineconfigdir && $delimiter && $summaryfile && $input_dir && $output_dir );
throw("Cannot find input directory $input_dir\n") unless -e $input_dir;
throw("Cannot find output directory $output_dir\n") unless -e $output_dir;
# get database hash
my %database_hash;
my %databases = %{$DATABASES};
for my $category (keys %databases ) {
  if(!$databases{$category}{-host} || !$databases{$category}{-dbname}){
    print STDERR "Skipping category ".$category." no dbinfo ".
      "defined\n";
    next;
  }
  print STDERR "\n$category: $databases{$category}{-host} $databases{$category}{-dbname} :\n--------------------------------------\n" ;
  my %constructor_args = %{$databases{$category}};
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     %constructor_args,
    );
  $database_hash{$category} = $dba;
}

my $dba = $database_hash{$dbname};
my $sa = $dba->get_SliceAdaptor;
# also need the pipeline adaptor
my %constructor_args = %{$databases{$dbname}};
my $pipelinea = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  (
   %constructor_args,
  );

# parse the summary file

open ( FILE, $summaryfile) or die("Cannot open summary file $summaryfile\n");
my $line = 0;
my @rows;
my %map;
while (<FILE>){
  my %data;
  chomp;
  next if $_ =~ /^#/;
  my @cells = split($delimiter,$_);
  if ( $line == 0 ) {
    # header row
    print STDERR "Got these headers\n";
    for (  my $i = 0 ; $i < $#cells ; $i++ ) {
      my $header =  $cells[$i];
      print STDERR $i+1 . ") $header\n";
    }
    print STDERR "Please assign some / all columns to the some / all of the following categories multiple values can be separted with commas:\n";
    print STDERR "Tag\tDescription\n";
    push(@{$map{"ID"}}, @{assign_categories(\@cells,1,$RNASEQCONFIG->{ID},"ID")});
    push(@{$map{"SM"}}, @{assign_categories(\@cells,1,$RNASEQCONFIG->{SM},"SM")});
    push(@{$map{"LB"}}, @{assign_categories(\@cells,0,$RNASEQCONFIG->{LB},"LB")});
    push(@{$map{"DS"}}, @{assign_categories(\@cells,0,$RNASEQCONFIG->{DS},"DS")});
    push(@{$map{"PU"}}, @{assign_categories(\@cells,0,$RNASEQCONFIG->{PU},"PU")});
    push(@{$map{"CN"}}, @{assign_categories(\@cells,0,$RNASEQCONFIG->{CN},"CN")});
    push(@{$map{"ST"}}, @{assign_categories(\@cells,0,$RNASEQCONFIG->{ST},"ST")});
    push(@{$map{"PL"}}, @{assign_categories(\@cells,0,$RNASEQCONFIG->{PL},"PL")});
    push(@{$map{"FILE"}}, @{assign_categories(\@cells,1,$RNASEQCONFIG->{FILE},"FILE")});
    push(@{$map{"LENGTH"}}, @{assign_categories(\@cells,0,$RNASEQCONFIG->{LENGTH},"LENGTH")});
    push(@{$map{"PAIRED"}}, @{assign_categories(\@cells,0,$RNASEQCONFIG->{PAIRED},"PAIRED")});
    print STDERR "Sample data:\n";
    $line++;
    next;
  } else {
     foreach my $key( keys %map ) {
       foreach my $col ( @{$map{$key}} ) {
	 unless ( $key eq "FILE" ) {
	   $data{$key} .= $cells[$col-1] ." ";
	 } else {
	   if ( $data{$key} ) {
	     $data{$key} .= "/" . $cells[$col-1];
	   } else {
	     $data{$key} = $cells[$col-1];
	   }
	 }
	 # no room for whitespace in the paired or length flag
	 $data{$key} =~ s/\s+//g if $key eq 'PAIRED' or $key eq 'LENGTH';
       }
     }
   }
  # add the path to the file name
 # $data{FILE} = $input_dir."/".$data{FILE};
  if ( $line == 1 ){
    foreach my $key ( keys %data ) {
      print STDERR "$key - " . $data{$key} ."\n";
    }
    print STDERR "Continue?(y/n)";
    my $reply = <>;
    chomp $reply;
    exit unless $reply eq "y" or $reply eq "Y";
  }
  push @rows,\%data;
  $line++;
}
print STDERR "Processed $line lines of data \n";

print  STDERR "Using " .$dba->dbc->dbname ."@" .$dba->dbc->host."  as pipeline db\n";

print STDERR "Creating analyses...\n";
# get the pipeline adaptor
# need to delete this hash ref in order 
# to call the pipeline version


# make and store the analyses
my $pipeline_analysis = $pipelinea->get_AnalysisAdaptor;
my $ra = $pipelinea->get_RuleAdaptor;
my $sic = $pipelinea->get_StateInfoContainer;
my %pairs;
# loop through the rows and create the analyses, rules and input_ids
 $line = 0;
foreach my $row ( @rows ) {
  $line++;
  #print  "ROW $line\t";
  #analyses
  $row->{ID} =~ s/ //g;
  my $ln = $row->{ID};

  my $submit_bwa    = new Bio::EnsEMBL::Pipeline::Analysis(
							   -logic_name      => "submit_".$ln."_bwa",
							   -input_id_type   => 'BWA'. $ln,
							  );
  my $bwa_wait    = new Bio::EnsEMBL::Pipeline::Analysis(
							   -logic_name      => "bwa_".$ln."_wait",
							   -module          => "Accumulator",
							   -input_id_type   => 'ACCUMULATOR',
							  );  

  my $submit_bwa2bam    = new Bio::EnsEMBL::Pipeline::Analysis(
							   -logic_name      => "submit_".$ln."_bwa2bam",
							   -input_id_type   => 'BWA2BAM'.$ln,
							  );
  my $skip=0;
  # catch paired analyses and store the file names
  if ( $all_paired || $row->{PAIRED} ) {
    # need to pair up the ids using the regexp
    if ( $row->{FILE} =~ /$regex/ ) {
      $pairs{$1."-".$3}->{$2} = $row->{FILE};
      $pairs{$1."-".$3}->{ANALYSIS} = $submit_bwa2bam;
      $skip = 1;
      unless ( $1 && $2 && $3 ) {
	throw("need 3 parts to the filename start - read num -  file extension. I have $1 $2 $3\n");
      }
      
    } else {
      throw("Cannot pair analyses of type " .$row->{ID} ." with filename " . $row->{FILE} . " and regex $regex\n");
    }
  }
  my $bwa    = new Bio::EnsEMBL::Pipeline::Analysis(
						    -logic_name      => "bwa_" .$ln,
						    -program         => "bwa",
						    -program_file    => "/software/solexa/bin/bwa",
						    -module          => "BWA",
						    -description     => $row->{DS},
						    -display_label   => $row->{ID},
						    -displayable     => '1',
						    -input_id_type   => 'BWA'.$ln,
						   );
  
  my $bwa2bam = new Bio::EnsEMBL::Pipeline::Analysis(
						     -logic_name      => "bwa2bam_" . $ln,
						     -program         => "bwa",
						     -program_file    => "/software/solexa/bin/bwa",
						     -module          => "BWA2BAM",
						     -description     => $row->{DS},
						     -display_label   => $row->{ID},
						     -displayable     => '1',
						     -input_id_type   => 'BWA2BAM'.$ln,
						    );
 
  my $bwa_rule = Bio::EnsEMBL::Pipeline::Rule->new
    (
     -goalanalysis => $bwa
    );
  $bwa_rule->add_condition($submit_bwa->logic_name);

  my $bwa_wait_rule = Bio::EnsEMBL::Pipeline::Rule->new
    (
     -goalanalysis => $bwa_wait
    );
  $bwa_wait_rule->add_condition($bwa->logic_name);

  my $bwa2bam_rule = Bio::EnsEMBL::Pipeline::Rule->new
    (
     -goalanalysis => $bwa2bam
    );
  $bwa2bam_rule->add_condition($bwa_wait->logic_name);
  $bwa2bam_rule->add_condition($submit_bwa2bam->logic_name);
  
#  print  "\n";
  # store the analyses
  if ($write) {
  #  print "Storing\n";
    $pipeline_analysis->store($submit_bwa);
    $pipeline_analysis->store($bwa_wait);
    $pipeline_analysis->store($submit_bwa2bam);
    $pipeline_analysis->store($bwa);
    $pipeline_analysis->store($bwa2bam);
    $ra->store($bwa_rule);
    $ra->store($bwa_wait_rule);
    $ra->store($bwa2bam_rule);
  }
  # input_ids

  if ( $write ) {
    $sic->store_input_id_analysis($row->{FILE},$submit_bwa,"dummy");
    $sic->store_input_id_analysis($row->{FILE},$submit_bwa2bam,"dummy") unless $skip;
  }
}
 
if ( $write ) {
  # store paired input ids
  foreach my $key ( keys %pairs ) {
    my $iid;
    throw("Cannot parse file names using regex $regex for lanes $key\n")
      unless scalar(keys  %{$pairs{$key}} == 3 ) ;
    foreach my $key2 ( keys %{$pairs{$key}} ) {
      next if $key2 eq 'ANALYSIS';
      $iid .= $pairs{$key}->{$key2} .":";
    }
    $iid =~ s/:$//;
    $sic->store_input_id_analysis($iid,$pairs{$key}->{ANALYSIS},"dummy") ;
  }
}



# check config directory
  

# write the extra information header files
foreach my $row ( @rows ) {
  open(HEAD,">$output_dir/" . $row->{ID} ."_header.txt") or die("Cannot open  $output_dir/" . $row->{ID} ."_header.txt for writing\n");
  print HEAD "\@RG\tID:" . $row->{ID} ."\tPU:" . $row->{PU} . "\tSM:" . $row->{SM} ."\t";
  print HEAD "LB:" . $row->{LB} ."\tDS:" . $row->{DS} . "\tCN:" . $row->{CN} ."\t";
  print HEAD "ST:" . $row->{ST} ."\tPL:" . $row->{PL}  ."\n";
}

print STDERR "Have these config files to modify:\n$analysisconfigdir/Genebuild/BWA.pm\n$pipelineconfigdir/BatchQueue.pm - backing them up\n";

system ("mv $pipelineconfigdir/BatchQueue.pm $pipelineconfigdir/BatchQueue.pm_bk") if -e "$pipelineconfigdir/BatchQueue.pm";
system ("mv $analysisconfigdir/GeneBuild/BWA.pm $analysisconfigdir/GeneBuild/BWA.pm_bk") if -e "$analysisconfigdir/GeneBuild/BWA.pm";

# write config
open(BATCHQUEUE,">$pipelineconfigdir/BatchQueue.pm") or die("Cannot open " .$pipelineconfigdir."/BatchQueue.pm  for writing\n");
open(BWA,">$analysisconfigdir/GeneBuild/BWA.pm") or die("Cannot open " .$analysisconfigdir."/GeneBuild/BWA.pm  for writing\n");

# write headers
my $string =  bwa_header();
print BWA $string;
# batchqueue
$string =  batchqueue_header();
print BATCHQUEUE $string;
#write logic name specific config
my %seen;
foreach my $row ( @rows ) {
  my $length = $RNASEQCONFIG->{READ_LENGTH};
  $length =  $row->{LENGTH} if $row->{LENGTH};
  next if $seen{$row->{ID}};
  # bwa
  print BWA "           'bwa_" . $row->{ID} ."' => {
                        INDIR   => \"". $input_dir  ."\",
                        OUTDIR  => \"". $output_dir ."\",
                        OPTIONS => \"-n " . ($length / 2). " -i " . $length  ."\",
                      },\n";
  
  # bwa2bam
  print BWA "           'bwa2bam_" . $row->{ID} ."' => {
                        HEADER  => \"$output_dir/" . $row->{ID} ."_header.txt\",
                        INDIR   => \"". $input_dir  ."\",
                        OUTDIR  => \"". $output_dir ."\",\n";
  print BWA "                        PAIRED => 1,\n" if $row->{PAIRED};
  print BWA "                 },\n";
  # batchqueue
  print BATCHQUEUE "       {
             logic_name => 'bwa_" . $row->{ID} ."',
             output_dir => '".$output_dir ."/".$row->{ID}."_pipeline',
       },
{
             logic_name => 'bwa_" . $row->{ID} ."_wait',
             output_dir => '".$output_dir ."/".$row->{ID}."_pipeline',
       },
       {
             logic_name => 'bwa2bam_" . $row->{ID} ."',
             output_dir => '".$output_dir ."/".$row->{ID}."_pipeline',
       },\n";
  
  $seen{$row->{ID}} = 1;
}

print BWA "     }\n";
print BATCHQUEUE "     ]\n";
print BWA tail();
print BATCHQUEUE tail();


sub assign_categories {
  my ($array,$required,$answer,$category) = @_;
  my @answers = split(/,/,$answer);
  if ( $required ){
    throw("Answer required for category $category\n")
      if scalar(@answers == 0 ) ;
  }
  my $chosen = [];
  print "Selections for category $category $answer:\n";
  foreach my $ans ( @answers ) {
    throw("Selection $ans not recognised\n") 
      unless $ans =~ /\d+/;
    throw("Selection $ans out of range unless\n") 
      unless $array->[$ans-1];
    print "$ans " . $array->[$ans-1] ."\n";
    push @$chosen,$ans;
  }
  return $chosen;
}


sub batchqueue_header {
  my $str = '# EnsEMBL module for Bio::EnsEMBL::Pipeline::Config::BatchQueue;
#
# You may distribute this module under the same terms as perl itself


=head1 NAME

    Bio::EnsEMBL::Pipeline::Config::BatchQueue

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
    use Bio::EnsEMBL::Pipeline::Config::BatchQueue qw();

=head1 DESCRIPTION

    Configuration for pipeline batch queues. Specifies per-analysis
    resources and configuration, e.g. so that certain jobs are run
    only on certain nodes.

    It imports and sets a number of standard global variables into
    the calling package. Without arguments all the standard variables
    are set, and with a list, only those variables whose names are
    provided are set. The module will die if a variable which doesn\'t
    appear in its C<%Config> hash is asked to be set.

    The variables can also be references to arrays or hashes.

    Edit C<%Config> to add or alter variables.

    All the variables are in capitals, so that they resemble
    environment variables.

    To run a job only on a certain host, you have to add specific
    resource-requirements. This can be useful if you have special
    memory-requirements, for example if you like to run the job only
    on linux 64bit machines or if you want to run the job only on a
    specific host group. The commands bmgroup and lsinfo show you
    information about certain host-types / host-groups.

    Here are some example resource-statements / sub_args statements:

        sub_args => \'-m bc_hosts\',              # only use hosts of host-group \'bc_hosts\' (see bmgroup)
        sub_args => \'-m bc1_1\',                 # only use hosts of host-group \'bc1_1\'

        resource => \'select[type==X86_64]\',     # use Linux 64 bit machines only
        resource => \'select[model==IBMBC2800]\', # only run on IBMBC2800 hosts

        resource => \'alpha\',                    # only run on DEC alpha
        resource => \'linux\',                    # run on any machine capable of running 32-bit X86 Linux apps

=head2 Database throttling

    This runs a job on a linux host, throttles ecs4:3350 to not have
    more than 300 active connections, 10 connections per job in the
    duration of the first 10 minutes when the job is running (means 30
    hosts * 10 connections = 300 connections):

        resource =>\'select[linux && ecs4my3350 <=300] rusage[ecs4my3350=10:duration=10]\',

    Running on \'linux\' hosts with not more than 200 active connections
    for myia64f and myia64g, 10 connections per job to each
    db-instance for the first 10 minutes:

        resource =>\'select[linux && myia64f <=200 && myia64g <=200] rusage[myia64f=10:myia64g=10:duration=10]\',

    Running on hosts of model \'IBMBC2800\' hosts with not more than
    200 active connections to myia64f, 10 connections per job for the
    first 10 minutes:

        resource =>\'select[model==IBMBC2800 && myia64f<=200] rusage[myia64f=10:duration=10]\',

    Running on hosts of host_group bc_hosts with not more than 200
    active connections to myia64f, 10 connections per job for the
    first 10 minutes:

        resource =>\'select[myia64f<=200] rusage[myia64f=10:duration=10]\',
        sub_args =>\'-m bc_hosts\'

=head1 LICENSE

    Copyright (c) 1999-2009 The European Bioinformatics Institute and
    Genome Research Limited.  All rights reserved.

    This software is distributed under a modified Apache license.
    For license details, please see

      http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

    Please email comments or questions to the public Ensembl
    developers list at <ensembl-dev@ebi.ac.uk>.

    Questions may also be sent to the Ensembl help desk at
    <helpdesk@ensembl.org>.

=cut


package Bio::EnsEMBL::Pipeline::Config::BatchQueue;

use strict;
use vars qw(%Config);

%Config = (

  # Depending on the job-submission-system you\'re using,
  # use LSF, you can also use \'Local\'.
  #
  # For mor info look into:
  # /ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/BatchSubmission
  QUEUE_MANAGER => \'LSF\', # use "SGE_GridEngine_v6" for running in ensembl cloud evironment

  DEFAULT_BATCH_SIZE  => 1,
  DEFAULT_RETRIES     => 3,
  DEFAULT_BATCH_QUEUE => \'long\',    # Put in the queue of your choice, eg. \'normal\'
  DEFAULT_RESOURCE    => \'select[mem>'.($RNASEQCONFIG->{MEMORY} / 1000) .'] rusage[mem='.($RNASEQCONFIG->{MEMORY} / 1000) .']\',
  DEFAULT_SUB_ARGS    => \'-M'.$RNASEQCONFIG->{MEMORY}.'\',
  DEFAULT_OUTPUT_DIR  => \'\',
  DEFAULT_CLEANUP     => \'no\',
  DEFAULT_VERBOSITY   => \'WARNING\',

  # The two variables below are to overcome a bug in LSF. 
  # We\'re currently running the pre-exec with a different perl. lsf currently unsets the LD_LIBRARY_PATH 
  # which we need for certain 64bit libraries in pre-exec commands. (more info see LSF_LD_SECURITY variable ) 
    
  DEFAULT_LSF_PRE_EXEC_PERL =>\'/usr/local/ensembl32/bin/perl\', # ONLY use 32bit perl for lsf -pre-exec jobs
  DEFAULT_LSF_PERL =>\'/usr/local/ensembl32/bin/perl\', # ONLY use ensembl64/bin/perl for memory jobs > 4 gb
                                                     # SANGER farm : don\'t forget to source source /software/intel_cce_80/bin/iccvars.csh for big mem jobs 
                                                     #
  # At this number of jobs RuleManager will sleep for a certain period
  # of time if you effectively want this never to run set the value to
  # very high ie 100000 for a certain period of time this is important
  # for queue managers which cannot cope with large numbers of pending
  # jobs (e.g. early LSF versions and SGE)
  JOB_LIMIT           => 10000,

  JOB_STATUSES_TO_COUNT => [\'PEND\'],    # These are the jobs which will be 
                                        # counted. valid statuses for
                                        # this array are RUN, PEND, SSUSP, EXIT, DONE ; use \'qw\' for Sun Grid Engine
  MARK_AWOL_JOBS => 1,
  MAX_JOB_SLEEP  => 3600,   # The maximun time to sleep for when job limit
                            # reached
  MIN_JOB_SLEEP => 120, # The minimum time to sleep for when job limit reached
  SLEEP_PER_JOB => 30,  # The amount of time to sleep per job when job limit
                        # reached

  DEFAULT_RUNNABLEDB_PATH => \'Bio/EnsEMBL/Analysis/RunnableDB\',

  DEFAULT_RUNNER         => \'\',
  DEFAULT_RETRY_QUEUE    => \'long\',
  DEFAULT_RETRY_SUB_ARGS => \'\',
  DEFAULT_RETRY_RESOURCE => \'\',

  QUEUE_CONFIG => [
';
  return $str;
}

sub batchqueue_tail {
  my $str = " ";
  return $str;
}

sub bwa_header {
  my $str = '# package Bio::EnsEMBL::Analysis::Config::GeneBuild::BWA
# 
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Config::GeneBuild::BWA

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::GeneBuild::BWA

=head1 DESCRIPTION

This contains the specific configuraton for 
Bio::EnsEMBL::Analysis::RunnableDB::BWA and 
Bio::EnsEMBL::Analysis::RunnableDB::BWA2BAM

=head1 CONTACT

=cut


package Bio::EnsEMBL::Analysis::Config::GeneBuild::BWA;

use strict;
use vars qw( %Config );

%Config = (
  BWA_CONFIG_BY_LOGIC =>  {
            DEFAULT =>  {

	      # base path to the fastq
	      INDIR => "/path/to/my/input",	
	      
	      # path to the output directory
	      OUTDIR => "/path/to/my/output",	
			      
	      # path to dumped genome file used for the alignment
	      # it will make an index for it if one does not already exist
	      GENOMEFILE => "'. $RNASEQCONFIG->{GENOME} .'",
	      
	      # alignment options
	      OPTIONS => "-n 20 -i 75",
	      
	      # options for BWA2BAM
	      #####################
	          
	      # are the reads paired end? (1/0)
	      PAIRED => "0",
       
       	      # parameters for sampe ( BWA paired alignment processing ) 
              SAMPE_OPTIONS => "-A -a 200000",
	      
	      #paramteres for samse ( BWA unpaired alignment processing)
	      SAMSE_OPTIONS => "",
	      
	      # path to the samtools binaries
	      SAMTOOLS_PATH => "' . $RNASEQCONFIG->{SAMTOOLS} .'",

              # optional header with additional information describing the sample
              HEADER => "/path/to/my/SAM/header/file.txt",
            },	
';
  return $str;
}

sub tail {
  my $str = ');

sub import {
  my ($callpack) = caller(0); # Name of the calling package
  my $pack = shift; # Need to move package off @_

  # Get list of variables supplied, or else everything
  my @vars = @_ ? @_ : keys( %Config );
  return unless @vars;
  
  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(\' \', map { \'$\'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $Config{$_} ) {
            no strict \'refs\';
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
';
  return $str;
}
