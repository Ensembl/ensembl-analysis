#!/usr/bin/env perl

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

use warnings ;
use setup_rnaseq_pipeline_config;
use vars qw(%Config);
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Analysis::Tools::ConfigWriter;


my $verbose = 0;
my $file;
my $dbname = $RNASEQCONFIG->{DB};
my $configdir = $RNASEQCONFIG->{CONFIG_DIR};
my $analysisconfigdir = $configdir.'/Bio/EnsEMBL/Analysis/Config';
my $pipelineconfigdir =  $configdir.'/Bio/EnsEMBL/Pipeline/Config';
my $delimiter = $RNASEQCONFIG->{DELIMITER};
my $summaryfile  = $RNASEQCONFIG->{SUMMARY};
my $input_dir = $RNASEQCONFIG->{INPUT_DIR};
my $output_dir = $RNASEQCONFIG->{OUTPUT_DIR};
my $all_paired = $RNASEQCONFIG->{ALL_PAIRED};
my $regex = $RNASEQCONFIG->{PAIRING_REGEX};
my $merge_dir = $RNASEQCONFIG->{MERGE_DIR};
$merge_dir = $RNASEQCONFIG->{OUTPUT_DIR} unless $merge_dir;
my $use_gsnap = $RNASEQCONFIG->{USE_GSNAP};
my $gsnap_path = $RNASEQCONFIG->{GSNAP_PATH};
my $rgt = $RNASEQCONFIG->{READ_GROUP_TAG};
my $queue_manager = $RNASEQCONFIG->{BATCHQUEUE_MANAGER} || 'LSF';
my $default_lsf_pre_exec_perl = $RNASEQCONFIG->{BATCHQUEUE_DEFAULT_LSF_PRE_EXEC_PERL} || '/software/ensembl/central/bin/perl';
my $default_lsf_perl = $RNASEQCONFIG->{BATCHQUEUE_DEFAULT_LSF_PERL} || '/software/ensembl/central/bin/perl';
my $splice_aligner = $RNASEQCONFIG->{SPLICING_ALIGNER} || '/software/ensembl/genebuild/usrlocalensemblbin/exonerate-0.9.0';

$rgt = 'ID' unless $rgt;
my %id_groups;
my %ids_by_tissue;
my %tissue_by_id;
my $stage = "Initialization";
my $bwa_analysis_written = 0;
my $slice_batches = 1;
my $rough_batches = 1;
my $check = 0 ;
my $rough_load;
my $refine_load;
my $ref_load;
my $blast_load;
my $update_analyses = 0;
my $force_stage = 0;
my $jdi;
my $use_existing;
my $update_config = 0;
my $bwa_run = 0;
my $introns_run = 0;
my $refine_run = 0;
my $blast_run = 0;

my $usage = "perl setup_rnaseq_pipeline.pl
-verbose    $verbose,
-check      $check, print out which columns are used for which RG tag
-update_analyses $update_analyses, only write the analyses - do not alter the config,
-jdi        Just Do it, skip the continue(Y/n) prompt.
Need to fill in the config in the setup_rnaseq_pipeline_config.pm module.
-stage      $force_stage, Force the pipeline to start from a particular stage -
            could be dangerous unless your pipeline has finished the previous stages but useful
        if you have run part of the pipeline and want to change config..
	    The stages are:
	    1	bwa_complete
	    2	bam2genes complete
	    3	bam2introns complete
	    Enter 1,2 or 3 to force the pipeline to start from the stage of choice.
            You must run the steps sequentially though so to go to stage 2 first run
            the script with no stage, then with -stage 1 and finally with -stage 2 to
            ensure all the analyses get written
";

$| = 1;

GetOptions( 'verbose!'         => \$verbose,
            'check!'           => \$check,
            'stage:s'          => \$force_stage,
            'jdi!'            => \$jdi,
            'update_analyses!' => \$update_analyses,
            'use_existing!'    => \$use_existing, );

die($usage) unless ($dbname && $analysisconfigdir &&
$pipelineconfigdir && $delimiter && $summaryfile && $input_dir &&
$output_dir );

my %database_hash;
my $rows;
if ($update_config) {
    $rows = read_pipeline_config();
    check_directories();
    update_config($bwa_run);
    update_config($introns_run);
    update_config($refine_run);
    update_config($blast_run);
}
elsif ($force_stage) {
    $rows = read_pipeline_config();
    check_directories();
    check_databases(\%database_hash);

}
else {
    $rows = read_pipeline_config();
    check_directories();
    check_databases(\%database_hash);
}




# also need the pipeline
my $pipelinea = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( -dbconn => $database_hash{REFERENCE_DB}->dbc );
# parse the summary file
# get the pipeline adaptor
# need to delete this hash ref in order
# to call the pipeline version
my $pipeline_analysis = $pipelinea->get_AnalysisAdaptor;
my $ra                = $pipelinea->get_RuleAdaptor;
my $sic               = $pipelinea->get_StateInfoContainer;
my $sa = $pipelinea->get_SliceAdaptor;
my $stored_ids;


# make the database load variables
$rough_load =  $database_hash{$RNASEQCONFIG->{ROUGHDB}}->dbc->host;
$rough_load =~ s/gene//;
$refine_load =  $database_hash{$RNASEQCONFIG->{REFINEDDB}}->dbc->host;
$refine_load =~ s/gene//;
$ref_load =  $database_hash{$RNASEQCONFIG->{DB}}->dbc->host;
$ref_load =~ s/gene//;
$blast_load =  $database_hash{$RNASEQCONFIG->{BLASTDB}}->dbc->host;
$blast_load =~ s/gene//;


print "Stage = $stage\n";


if ( $stage eq "Initialization" ) {
  print STDERR 'Processed ', scalar(@$rows), "lines of data \n";
  print STDERR "Using " . $pipelinea->dbc->dbname . "@" . $pipelinea->dbc->host . "  as pipeline db\n";
  print STDERR "Creating analyses...\n";
}

# figure out the relationship between the tissues and the ids
foreach my $key1 ( keys %id_groups ) {
  foreach my $key2 (  keys %{$id_groups{$key1}} ) {
    $key1 =~ s/ //g;
    $key2 =~ s/ //g;
    $tissue_by_id{$key2} = $key1;
    $ids_by_tissue{$key1} .= $key2.'","';
  }
}



my %pairs;
# loop through the rows and create the analyses, rules and input_ids
my $line = 0;
foreach my $row (@$rows) {
  $line++;
  #analyses
  $row->{ID} =~ s/ //g;
  my $ln = $row->{ID};
  # trim off trailing commas
  $ids_by_tissue{$tissue_by_id{$ln}} =~ s/"\,"$//;
  my $submit_bwa =
    new Bio::EnsEMBL::Pipeline::Analysis(
                                      -logic_name => "submit_" . $ln . "_bwa",
                                      -input_id_type => 'BWA' . $ln, );
  my $bwa_wait =
    new Bio::EnsEMBL::Pipeline::Analysis(
                                        -logic_name => "bwa_" . $ln . "_wait",
                                        -module     => "Accumulator",
                                        -input_id_type => 'ACCUMULATOR', );
  my $submit_gsnap =
    new Bio::EnsEMBL::Pipeline::Analysis(
                                    -logic_name => "submit_" . $ln . "_gsnap",
                                    -input_id_type => 'GSNAP' . $ln, );

  my $submit_bwa2bam =
    new Bio::EnsEMBL::Pipeline::Analysis(
                                  -logic_name => "submit_" . $ln . "_bwa2bam",
                                  -input_id_type => 'BWA2BAM' . $ln, );
  my $refine =
    new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "refine_" . $tissue_by_id{$ln},
                                          -input_id_type => 'CHROMOSOME',
                                          -module => 'RefineSolexaGenes', );
  my $skip = 0;
  # catch paired analyses and store the file names
  if ( $all_paired || $row->{PAIRED} ) {
    # need to pair up the ids using the regexp
    if ( $row->{FILE} =~ /$regex/ ) {
      $pairs{ $1 . "-" . $3 }->{$2} = $row->{FILE};
      unless ($use_gsnap) {
        $pairs{ $1 . "-" . $3 }->{ANALYSIS} = $submit_bwa2bam;
      } else {
        $pairs{ $1 . "-" . $3 }->{ANALYSIS} = $submit_gsnap;
      }
      $skip = 1;
      unless ( $1 && $2 && $3 ) {
        throw(   "need 3 parts to the filename start - read num -  "
               . "file extension. I have $1 $2 $3\n" );
      }

    } else {
      throw(   "Cannot pair analyses of type " . $row->{ID}
             . " with filename " . $row->{FILE}
             . " and regex $regex\n" );
    }
  }
  my $bwa =
    new Bio::EnsEMBL::Pipeline::Analysis(
                                  -logic_name   => "bwa_" . $ln,
                                  -program      => "bwa",
                                  -program_file => $RNASEQCONFIG->{BWA_PATH},
                                  -module       => "BWA",
                                  -description  => $row->{DS},
                                  -display_label => $row->{ID},
                                  -displayable   => '1',
                                  -input_id_type => 'BWA' . $ln, );

  my $bwa2bam =
    new Bio::EnsEMBL::Pipeline::Analysis(
                                  -logic_name   => "bwa2bam_" . $ln,
                                  -program      => "bwa",
                                  -program_file => $RNASEQCONFIG->{BWA_PATH},
                                  -module       => "BWA2BAM",
                                  -description  => $row->{DS},
                                  -display_label => $row->{ID},
                                  -displayable   => '1',
                                  -input_id_type => 'BWA2BAM' . $ln, );
  my $gsnap =
    new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "gsnap_" . $ln,
                                          -program       => "gsnap",
                                          -program_file  => $gsnap_path,
                                          -module        => "Gsnap",
                                          -description   => $row->{DS},
                                          -display_label => $row->{ID},
                                          -displayable   => '1',
                                          -input_id_type => 'GSNAP' . $ln, );

  my $bwa_rule = Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $bwa );
  $bwa_rule->add_condition( $submit_bwa->logic_name );

  my $bwa_wait_rule =
    Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $bwa_wait );
  $bwa_wait_rule->add_condition( $bwa->logic_name );

  my $gsnap_rule =
    Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $gsnap );
  $gsnap_rule->add_condition( $submit_gsnap->logic_name );

  my $bwa2bam_rule =
    Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $bwa2bam );
  $bwa2bam_rule->add_condition( $bwa_wait->logic_name );
  $bwa2bam_rule->add_condition( $submit_bwa2bam->logic_name );
  my $refine_rule =
    Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $refine );
  $refine_rule->add_condition("submit_chromosome");
  $refine_rule->add_condition("sam2bam_wait");
  # store the analyses
  unless ($use_gsnap) {
    $pipeline_analysis->store($submit_bwa);
    $pipeline_analysis->store($bwa_wait);
    $pipeline_analysis->store($submit_bwa2bam);
    $pipeline_analysis->store($bwa);
    $pipeline_analysis->store($bwa2bam);
    $ra->store($bwa_rule)      if check_rule($bwa_rule);
    $ra->store($bwa_wait_rule) if check_rule($bwa_wait_rule);
    $ra->store($bwa2bam_rule)  if check_rule($bwa2bam_rule);
    # input_ids
    # dont store duplicate ids
    $sic->store_input_id_analysis( $row->{FILE}, $submit_bwa, "dummy" )
      unless $stored_ids->{ $row->{FILE} }->{ $submit_bwa->logic_name };
    $sic->store_input_id_analysis( $row->{FILE}, $submit_bwa2bam, "dummy" )
      unless ( $skip
          || $stored_ids->{ $row->{FILE} }->{ $submit_bwa2bam->logic_name } );
  } else {
    $pipeline_analysis->store($submit_gsnap);
    $pipeline_analysis->store($gsnap);
    $ra->store($gsnap_rule) if check_rule($gsnap_rule);
    # input_ids
    # dont store duplicate ids
    $sic->store_input_id_analysis( $row->{FILE}, $submit_gsnap, "dummy" )
      unless ( $skip
            || $stored_ids->{ $row->{FILE} }->{ $submit_gsnap->logic_name } );
  }
  $pipeline_analysis->store($refine)
    if $RNASEQCONFIG->{SINGLE_TISSUE} && $stage eq 'configured';
  if ( $RNASEQCONFIG->{SINGLE_TISSUE} && $stage eq 'configured' ) {
    $ra->store($refine_rule) if check_rule($refine_rule);
  }
} ## end foreach my $row (@rows)


# store paired input ids
foreach my $key ( keys %pairs ) {
  my $iid;
  throw("Cannot parse file names using regex $regex for lanes $key\n")
    unless scalar( keys %{ $pairs{$key} } == 3 );
  foreach my $key2 ( sort keys %{ $pairs{$key} } ) {
    next if $key2 eq 'ANALYSIS';
    $iid .= $pairs{$key}->{$key2} . ":";
    #print "KEy $key2 $key $iid\n";
    #print $pairs{$key}->{ANALYSIS}->logic_name."\n";
  }
  $iid =~ s/:$//;
  $sic->store_input_id_analysis( $iid, $pairs{$key}->{ANALYSIS}, "dummy" )
    unless $stored_ids->{$iid}->{ $pairs{$key}->{ANALYSIS}->logic_name };
}


# RNASeq pipeline anaysis that are run on all lanes
my $submit_chromosome =
  new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "submit_chromosome",
                                        -input_id_type => 'CHROMOSOME', );
my $bam2genes =
  new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "bam2genes",
                                        -input_id_type => 'CHROMOSOME',
                                        -module        => "Bam2Genes", );

my $submit_bam2introns =
  new Bio::EnsEMBL::Pipeline::Analysis( -logic_name => "submit_bam2introns",
                                        -input_id_type => 'STABLEID', );
my $bam2introns =
  new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "bam2introns",
                                        -input_id_type => 'STABLEID',
                                        -program_file => $splice_aligner,
                                        -module        => 'Bam2Introns', );
my $bam2introns_wait =
  new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "bam2introns_wait",
                                        -module        => "Accumulator",
                                        -input_id_type => 'ACCUMULATOR', );

my $submit_sam2bam =
  new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "submit_sam2bam",
                                        -input_id_type => 'GENOME', );
my $sam2bam =
  new Bio::EnsEMBL::Pipeline::Analysis(
                                   -logic_name   => "sam2bam",
                                   -program_file => $RNASEQCONFIG->{SAMTOOLS},
                                   -input_id_type => 'GENOME',
                                   -module        => 'Sam2Bam', );
my $sam2bam_wait =
  new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "sam2bam_wait",
                                        -module        => "Accumulator",
                                        -input_id_type => 'ACCUMULATOR', );

my $refine_all =
  new Bio::EnsEMBL::Pipeline::Analysis( -logic_name    => "refine_all",
                                        -input_id_type => 'CHROMOSOME',
                                        -module        => 'RefineSolexaGenes',
  );
my $rnaseq_blast;
if ($RNASEQCONFIG->{BLASTP} eq 'ncbi') {
  $rnaseq_blast = new Bio::EnsEMBL::Pipeline::Analysis(
                                  -logic_name    => "rnaseqblast",
                                  -input_id_type => 'CHROMOSOME',
                                  -module        => 'BlastRNASeqPep',
                                  -parameters => '-p blastp -a 1 -A 40 -F F',
                                  -program_file => 'blastall',
                                  -program      => 'blastall',
                                  -db_file      => $RNASEQCONFIG->{UNIPROTDB},
  );
}
elsif ($RNASEQCONFIG->{BLASTP} eq 'wu') {
  $rnaseq_blast = new Bio::EnsEMBL::Pipeline::Analysis(
                                  -logic_name    => "rnaseqblast",
                                  -input_id_type => 'CHROMOSOME',
                                  -module        => 'BlastRNASeqPep',
                                  -parameters => '-cpus 1 -hitdist 40',
                                  -program_file => 'wublastp',
                                  -program      => 'wublastp',
                                  -db_file      => $RNASEQCONFIG->{UNIPROTDB},
  );
}
else {
  die("I don't know your blastp program".$RNASEQCONFIG->{BLASTP});
}


my $bam2genes_rule =
  Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $bam2genes );
$bam2genes_rule->add_condition("submit_chromosome");

my $bam2introns_rule =
  Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $bam2introns );
$bam2introns_rule->add_condition("submit_bam2introns");

my $bam2introns_wait_rule =
  Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $bam2introns_wait );
$bam2introns_wait_rule->add_condition("bam2introns");

my $sam2bam_rule =
  Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $sam2bam );
$sam2bam_rule->add_condition("submit_sam2bam");
$sam2bam_rule->add_condition("bam2introns_wait");

my $sam2bam_wait_rule =
  Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $sam2bam_wait );
$sam2bam_wait_rule->add_condition("sam2bam");

my $refine_all_rule =
  Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $refine_all );
$refine_all_rule->add_condition("submit_chromosome");
unless ($use_gsnap) {
  $refine_all_rule->add_condition("sam2bam_wait");
} else {
  $refine_all_rule->add_condition("bam2genes");
}

my $rnaseqblast_rule =
  Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $rnaseq_blast );

$rnaseqblast_rule->add_condition("submit_chromosome");
$rnaseqblast_rule->add_condition("refine_all");

# store the analyses

print "Stage $stage\n";
$pipeline_analysis->store($submit_chromosome);
$pipeline_analysis->store($bam2genes);

unless ($use_gsnap) {
  $pipeline_analysis->store($submit_bam2introns);
  $pipeline_analysis->store($bam2introns);
  $pipeline_analysis->store($bam2introns_wait);
  $pipeline_analysis->store($submit_sam2bam);
  $pipeline_analysis->store($sam2bam);
  $pipeline_analysis->store($sam2bam_wait);
}

$pipeline_analysis->store($refine_all);
$pipeline_analysis->store($rnaseq_blast);

  $ra->store($bam2genes_rule) if check_rule($bam2genes_rule);
    $ra->store($bam2introns_rule)      if check_rule($bam2introns_rule);
    $ra->store($bam2introns_wait_rule) if check_rule($bam2introns_wait_rule);
    $ra->store($sam2bam_rule)      if check_rule($sam2bam_rule);
    $ra->store($sam2bam_wait_rule) if check_rule($sam2bam_wait_rule);
  $ra->store($refine_all_rule)  if check_rule($refine_all_rule);
  $ra->store($rnaseqblast_rule) if check_rule($rnaseqblast_rule);
  # need to add a dummy input id for submit_sam2bam
    $sic->store_input_id_analysis( 'dummy', $submit_sam2bam, "dummy" )
      unless $stored_ids->{'dummy'}->{'submit_sam2bam'};

exit if $update_analyses;

# check config directory
# write the extra information header files
open( ALL, ">$output_dir/all_headers.txt" ) or die("Cannot open  $output_dir/all_headers.txt for writing\n");


my %seen_it;
foreach my $row (@$rows) {
    next if $seen_it{ $row->{ID} };
    $seen_it{ $row->{ID} } = 1;
    open( HEAD, ">$output_dir/" . $row->{ID} . "_header.txt" )
        or die( "Cannot open  $output_dir/" . $row->{ID} . "_header.txt for writing\n" );
    my $read_group = "\@RG\tID:" . $row->{ID};
    $read_group .= "\tPU:" . $row->{PU} if (exists $row->{PU});
    $read_group .= "\tSM:" . lc($row->{SM}) if (exists $row->{SM});
    $read_group .= "\tLB:" . $row->{LB} if (exists $row->{LB});
    $read_group .= "\tDS:" . $row->{DS} if (exists $row->{DS});
    $read_group .= "\tCN:" . $row->{CN} if (exists $row->{CN});
    # Might need to change ST to DT but we will need to change the previous configs
    $read_group .= "\tDT:" . $row->{ST} if (exists $row->{ST});
    if (exists $row->{PL})
    {
        my ($field_pl) = $row->{PL} =~ /([^+]+)/;
        $read_group .= "\tPL:" . $field_pl ;
    }
    $read_group .= "\n" ;
    print HEAD $read_group;
    print ALL $read_group;
    close(HEAD) || die( "Cannot close $output_dir/" . $row->{ID} . "_header.txt for writing\n" );
} ## end foreach my $row (@rows)
close(ALL) || die("Cannot close $output_dir/all_headers.txt for writing\n");





my $is_example = 1;
if ($stage eq "Initialization" || $check) {
  $is_example = 1;
  print STDERR "Have these config files to modify:\n
    $analysisconfigdir/GeneBuild/Bam2Genes.pm\n";
  if ( $use_gsnap ) {
      print STDERR "$analysisconfigdir/GeneBuild/Gsnap.pm\n";
  } else {
      print STDERR "$analysisconfigdir/GeneBuild/BWA.pm\n
        $analysisconfigdir/GeneBuild/Bam2Introns.pm\n
        $analysisconfigdir/GeneBuild/Sam2Bam.pm\n";
  }
  print STDERR "$analysisconfigdir/GeneBuild/RefineSolexaGenes.pm\n
    $analysisconfigdir/GeneBuild/BlastRNASeqPep.pm\n
    $pipelineconfigdir/BatchQueue.pm\n
     - backing them up\n" ;
}
$is_example = 0 if ($force_stage and $use_existing);

my $bq_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
    -modulename => 'Bio::EnsEMBL::Pipeline::Config::BatchQueue',
    -moduledir  => $configdir,
    -is_example => $use_existing ? 0 : 1);
$bq_cfg->root_value('DEFAULT_LSF_PERL', $default_lsf_perl);
$bq_cfg->root_value('DEFAULT_OUTPUT_DIR', $output_dir);
my $blast_general_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
    -modulename => 'Bio::EnsEMBL::Analysis::Config::Blast',
    -moduledir  => $configdir,
    -is_example => $use_existing ? 0 : 1);
$blast_general_cfg->write_config(1);
my $blast_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
    -modulename => 'Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastRNASeqPep',
    -moduledir  => $configdir,
    -is_example => $is_example);
    $blast_cfg->default_value('OUTPUT_DB', $RNASEQCONFIG->{BLASTDB});
    $blast_cfg->default_value('MODEL_DB', $RNASEQCONFIG->{REFINEDDB});
    $blast_cfg->default_value('LOGICNAME', 'refine_all');
    $blast_cfg->default_value('INDEX', $RNASEQCONFIG->{UNIPROTINDEX});
$blast_cfg->write_config(1);
my $bam2genes_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
    -modulename => 'Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Genes',
    -moduledir  => $configdir,
    -is_example => $is_example);
$bam2genes_cfg->add_analysis_to_config('bam2genes',
    {'OUTPUT_DB' => $RNASEQCONFIG->{ROUGHDB},
    'ALIGNMENT_BAM_FILE' => $RNASEQCONFIG->{MERGE_DIR}.'/merged.bam'},
    );
$bam2genes_cfg->default_value('ALIGNMENT_BAM_FILE', $RNASEQCONFIG->{MERGE_DIR}.'/merged.bam');
my $refine_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
    -modulename => 'Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes',
    -moduledir  => $configdir,
    -is_example => $is_example);
$refine_cfg->default_value('MODEL_DB', $RNASEQCONFIG->{ROUGHDB});
my $introns_bam_files = $refine_cfg->default_value('INTRON_BAM_FILES');
$introns_bam_files->[0]->{DEPTH} = 0;
$introns_bam_files->[0]->{GROUPNAME} = [];
$introns_bam_files->[0]->{FILE} = $RNASEQCONFIG->{MERGE_DIR}.'/merged.bam';
my $gsnap_cfg;
my $bwa_cfg;
my $bam2introns_cfg;
my $sam2bam_cfg;
if ($use_gsnap) {
    $introns_bam_files->[0]->{MIXED_BAM} = 1;
    $gsnap_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
        -modulename => 'Bio::EnsEMBL::Analysis::Config::GeneBuild::Gsnap',
        -moduledir  => $configdir,
        -is_example => $is_example);
    $gsnap_cfg->default_value('GENOMENAME', $RNASEQCONFIG->{GENOME_FILE});
    $gsnap_cfg->default_value('GENOMEDIR', $RNASEQCONFIG->{GENOME_DIR});
}
else {
    $introns_bam_files->[0]->{MIXED_BAM} = 0;
    $bwa_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
        -modulename => 'Bio::EnsEMBL::Analysis::Config::GeneBuild::BWA',
        -moduledir  => $configdir,
        -is_example => $is_example);
    $bwa_cfg->default_value('GENOMEFILE', $RNASEQCONFIG->{GENOME_DIR}.'/'.$RNASEQCONFIG->{GENOME_FILE});
    $bwa_cfg->default_value('SAMTOOLS_PATH', $RNASEQCONFIG->{SAMTOOLS});
    $bam2introns_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
        -modulename => 'Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Introns',
        -moduledir  => $configdir,
        -is_example => $is_example);
    $bam2introns_cfg->default_value('TRANSDB', $RNASEQCONFIG->{ROUGHDB});
    $bam2introns_cfg->add_analysis_to_config('bam2introns', {
        'OUT_SAM_DIR' => $RNASEQCONFIG->{MERGE_DIR}.'/SAM',
        'BAM_FILE' => $RNASEQCONFIG->{MERGE_DIR}.'/merged.bam'}
        );
    $sam2bam_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new(
        -modulename => 'Bio::EnsEMBL::Analysis::Config::GeneBuild::Sam2Bam',
        -moduledir  => $configdir,
        -is_example => $is_example);
    $sam2bam_cfg->default_value('HEADERFILE', $RNASEQCONFIG->{OUTPUT_DIR}.'/all_headers.txt');
    $sam2bam_cfg->default_value('GENOMEFILE', $RNASEQCONFIG->{GENOME_DIR}.'/'.$RNASEQCONFIG->{GENOME_FILE});
    $sam2bam_cfg->add_analysis_to_config('sam2bam', {
        'SAM_DIR' => $RNASEQCONFIG->{MERGE_DIR}.'/SAM',
        'BAMFILE' => $RNASEQCONFIG->{MERGE_DIR}.'/introns.bam'}
    );
}
$refine_cfg->add_analysis_to_config('refine_all', {
    'OUTPUT_DB' => $RNASEQCONFIG->{REFINEDDB},
    'INTRON_BAM_FILES' => $introns_bam_files,}
    );
$bq_cfg->add_analysis_to_batchqueue('refine_all',
  { output_dir => $output_dir .'/refine_all_pipeline',
    batch_size => $slice_batches,
    memory     => [ '1GB', '2GB', '5GB', '15GB','30GB' ],
    resource   => 'rusage[myens_'.$ref_load.'tok=25, myens_'.$refine_load.'tok=25, myens_'.$refine_load.'tok=25]]',
    retry_resource   => 'rusage[myens_'.$ref_load.'tok=25, myens_'.$refine_load.'tok=25, myens_'.$refine_load.'tok=25]]',
  });
$bq_cfg->add_analysis_to_batchqueue('bam2genes',
  { output_dir => $output_dir .'/bam2genes_pipeline',
    batch_size => $slice_batches,
    memory     => [ '1GB', '2GB', '5GB', '15GB','30GB' ],
    resource   => 'rusage[myens_'.$ref_load.'tok=25, myens_'.$rough_load.'tok=25]',
    retry_resource   => 'rusage[myens_'.$ref_load.'tok=25, myens_'.$rough_load.'tok=25]',
  });
$bq_cfg->add_analysis_to_batchqueue('bam2genes_wait',
  { output_dir => $output_dir .'/bam2genes_wait_pipeline',
    batch_size => 1,
  });
$bq_cfg->add_analysis_to_batchqueue('bam2introns',
  { output_dir => $output_dir .'/bam2introns_pipeline',
    batch_size => $rough_batches,
    memory     => [ '2GB', '5GB', '15GB', '20GB', '30GB' ],
    resource   => 'rusage[myens_'.$ref_load.'tok=25]',
    retry_resource   => 'rusage[myens_'.$ref_load.'tok=25]',
  });
$bq_cfg->add_analysis_to_batchqueue('bam2introns_wait',
  { output_dir => $output_dir .'/bam2introns_wait_pipeline',
  });
$bq_cfg->add_analysis_to_batchqueue('sam2bam',
  { output_dir => $output_dir .'/sam2bam_pipeline',
    batch_size => 1,
    memory     => [ '2GB', '5GB', '10GB', '20GB', '30GB' ],
    resource   => 'rusage[myens_'.$ref_load.'tok=25]',
    retry_resource   => 'rusage[myens_'.$ref_load.'tok=25]',
  });
$bq_cfg->add_analysis_to_batchqueue('sam2bam_wait',
  { output_dir => $output_dir .'/sam2bam_wait_pipeline',
    batch_size => 1,
  });
$bq_cfg->add_analysis_to_batchqueue('rnaseqblast',
  { output_dir => $output_dir .'/rnaseqblast_pipeline',
    batch_size => $slice_batches,
    memory     => [ '1GB', '1GB', '2GB' ],
    resource   => 'rusage[myens_'.$ref_load.'tok=25, myens_'.$blast_load.'tok=25]',
    retry_resource   => 'rusage[myens_'.$ref_load.'tok=25, myens_'.$blast_load.'tok=25]',
  });

my %seen;
foreach my $row ( @$rows ) {
  my $length = $RNASEQCONFIG->{READ_LENGTH};
  $length = $row->{LENGTH} if $row->{LENGTH};
  next if $seen{$row->{ID}};
  my $bam_arrayref = [
              FILE => $RNASEQCONFIG->{MERGE_DIR}.'/merged.bam',
              MIXED_BAM => 1,
              DEPTH => 0,
              GROUPNAME => [$ids_by_tissue{$tissue_by_id{$row->{ID}}}],
              ];
  if ($use_gsnap) {
    my %gsnap_hash = (
                      INDIR  => $input_dir,
                      OUTDIR => $output_dir,
                      HEADER => $output_dir.'/'.$row->{ID}.'_header.txt',
                      PAIRED => defined $row->{PAIRED} ? 1 : 0,
                     );
    $gsnap_cfg->add_analysis_to_config('gsnap_'. $row->{ID}, \%gsnap_hash);
    $bq_cfg->add_analysis_to_batchqueue('gsnap_' . $row->{ID},
        { output_dir => $output_dir .'/'.$row->{ID}.'_pipeline',
          memory     => [ '10GB', '20GB', '30GB' ],
          queue      => 'long',
          resource   => 'rusage[myens_'.$ref_load.'tok=25]',
          retry_resource   => 'rusage[myens_'.$ref_load.'tok=25]',
        });
  }
  else {
    my %bwa_hash = (
                      INDIR  => $input_dir,
                      OUTDIR => $output_dir,
                      OPTIONS => '-n ' . int($length / 2). ' -i ' . $length,
                     );
    $bwa_hash{OPTIONS} .= ' -t '.$RNASEQCONFIG->{USE_THREADS} if ( $RNASEQCONFIG->{USE_THREADS});
    $bwa_cfg->add_analysis_to_config('bwa_'. $row->{ID}, \%bwa_hash);
    my %bwa2bam_hash = (
                      INDIR  => $input_dir,
                      OUTDIR => $output_dir,
                      HEADER => $output_dir.'/' . $row->{ID} .'_header.txt',
                      PAIRED => defined $row->{PAIRED} ? 1 : 0,
                     );
    $bwa_cfg->add_analysis_to_config('bwa2bam_'. $row->{ID}, \%bwa2bam_hash);
    my %bq_bwa_cfg = (
        output_dir => $output_dir .'/'.$row->{ID}.'_pipeline',
        memory     => [ '5GB', '10GB', '20GB', '30GB' ],
        resource   => 'rusage[myens_'.$ref_load.'tok=25]',
        retry_resource   => 'rusage[myens_'.$ref_load.'tok=25]',
        );
    if ($RNASEQCONFIG->{USE_THREADS}) {
        $bq_bwa_cfg{resource} .= ' span[hosts=1]';
        $bq_bwa_cfg{retry_resource} .= ' span[hosts=1]';
        $bq_bwa_cfg{sub_args} .= ' -n '.$RNASEQCONFIG->{USE_THREADS};
        $bq_bwa_cfg{retry_sub_args} .= ' -n '.$RNASEQCONFIG->{USE_THREADS};
    }
    $bq_cfg->add_analysis_to_batchqueue('bwa_' . $row->{ID}, \%bq_bwa_cfg);
    $bq_cfg->add_analysis_to_batchqueue('bwa_' . $row->{ID} . '_wait',
        { output_dir => $output_dir .'/'.$row->{ID}.'_pipeline',
          resource   => 'rusage[myens_'.$ref_load.'tok=25]',
        });
    $bq_cfg->add_analysis_to_batchqueue('bwa2bam_' . $row->{ID}, \%bq_bwa_cfg);
  }
  $refine_cfg->value_by_logic_name('refine_'. $tissue_by_id{$row->{ID}}, 'SINGLE_EXON_MODEL', 'single_exon_'.$tissue_by_id{$row->{ID}});
  if ($RNASEQCONFIG->{SINGLE_TISSUE} && !$seen{$tissue_by_id{$row->{ID}}}) {
      $refine_cfg->add_analysis_to_config('refine_'. $tissue_by_id{$row->{ID}}, {
          INTRON_BAM_FILES  => $bam_arrayref,
          SINGLE_EXON_MODEL => 'single_exon_'.$tissue_by_id{$row->{ID}},
          });
      $bq_cfg->add_analysis_to_batchqueue('refine_' . $tissue_by_id{$row->{ID}},
          { output_dir => $output_dir .'/refine_'.$tissue_by_id{$row->{ID}}.'_pipeline',
            memory     => [ '1GB', '2GB', '5GB', '15GB','30GB' ],
            resource   => 'rusage[myens_'.$ref_load.'tok=25, myens_'.$refine_load.'tok=25], myens_'.$rough_load.'tok=25]',
            retry_resource   => 'rusage[myens_'.$ref_load.'tok=25, myens_'.$refine_load.'tok=25], myens_'.$rough_load.'tok=25]',
          });
  }
  $seen{$tissue_by_id{$row->{ID}}} = 1;
  $seen{$row->{ID}} = 1;
}

$bq_cfg->write_config(1);
$bam2genes_cfg->write_config(1);
$refine_cfg->write_config(1);

if ($use_gsnap) {
    $gsnap_cfg->write_config(1);
}
else {
    $bwa_cfg->write_config(1);
    $bam2introns_cfg->write_config(1);
    $sam2bam_cfg->write_config(1);
}

make_bam2genes_input_ids($sic, $pipeline_analysis, $sa);

sub assign_categories {
  my ($array,$required,$answer,$category) = @_;
  my @answers = split(/,/,$answer);
  if ( $required ){
    throw("Answer required for category $category\n")
      if scalar(@answers == 0 ) ;
  }
  my $chosen = [];
  print "Selections for category $category $answer:\n"  if ($stage eq "Initialization") || $check;
  foreach my $ans ( @answers ) {
    throw("Selection $ans not recognised\n")
      unless $ans =~ /\d+/;
    throw("Selection $ans out of range unless\n")
      unless $array->[$ans-1];
    print "$ans " . $array->[$ans-1] ."\n"  if ($stage eq "Initialization") || $check;
    push @$chosen,$ans ;
  }
  return $chosen;
}





sub check_rule {
  my ($rule) = @_;
  #  print "CHECK $check \n";
  # see if you can fetch it first
  my $check = $ra->fetch_by_goal($rule->goalAnalysis);
  unless ( $check ) {
    return 1;
  }
  return;
}

sub check_directories {
    my $throw_msg = '';
    $throw_msg .= "Cannot find input directory $input_dir\n" unless (-e $input_dir);
    $throw_msg .= "Cannot find output directory $output_dir\n" unless (-e $output_dir);
    $throw_msg .= "Cannot find merge directory $merge_dir\n" unless (-e $merge_dir);
    system( "mkdir -p " . $merge_dir . "/SAM" ) unless -e $merge_dir . "/SAM";
    throw($throw_msg) if ($throw_msg);
}

sub check_databases {
  my $databases = shift;
  foreach my $config_dbname ('REFERENCE_DB', $RNASEQCONFIG->{ROUGHDB}, $RNASEQCONFIG->{REFINEDDB}, $RNASEQCONFIG->{BLASTDB}) {
      if (exists $DATABASES->{$config_dbname}) {
          my %constructor_args = %{ $DATABASES->{$config_dbname} };
          $databases->{$config_dbname} = new Bio::EnsEMBL::DBSQL::DBAdaptor( %constructor_args, );
      }
      else {
          throw("Database $config_dbname is not defined")
      }
  }

  my $dba = $databases->{REFERENCE_DB};
  # test that repeat_feature, repeat_consensus, meta_cord and dna are populated in the ref db
  my $throw_msg = '';
  my $dbname = $dba->dbname;
  foreach my $table ('repeat_feature', 'repeat_consensus', 'dna', 'meta_coord') {
      my $sth = $dba->dbc->prepare( "SELECT COUNT(*) FROM $table" );
      $sth->execute;
      $throw_msg .= "Db $dbname has unpopulated $table table\n" unless ($sth->fetchrow > 0);
  }
  throw($throw_msg) if ($throw_msg);
}

sub read_pipeline_config {

    my $line = 0;
    my @rows;
    my %map;
    open( FILE, $summaryfile ) or die("Cannot open summary file $summaryfile\n");
    while (<FILE>) {
        next if ($_ =~ /^#/) or ($_ eq '');
        chomp;
        my %data;
        my @cells = split( $delimiter, $_ );
        if ( $line == 0 ) {
            # header row
            print STDOUT "Got these headers\n";
            for ( my $i = 0 ; $i <= $#cells ; $i++ ) {
                my $header = $cells[$i];
                print STDOUT $i + 1 . ") $header\n";
            }

            print STDOUT "Please assign some / all columns to the some / all ",
                "of the following categories multiple values can be separated with commas:\n",
                "Tag\tDescription\n";

            push( @{ $map{"ID"} }, @{ assign_categories( \@cells, 1, $RNASEQCONFIG->{ID}, "ID" ) } );
            push( @{ $map{"SM"} }, @{ assign_categories( \@cells, 1, $RNASEQCONFIG->{SM}, "SM" ) } );
            push( @{ $map{"LB"} }, @{ assign_categories( \@cells, 0, $RNASEQCONFIG->{LB}, "LB" ) } );
            push( @{ $map{"DS"} }, @{ assign_categories( \@cells, 0, $RNASEQCONFIG->{DS}, "DS" ) } );
            push( @{ $map{"PU"} }, @{ assign_categories( \@cells, 0, $RNASEQCONFIG->{PU}, "PU" ) } );
            push( @{ $map{"CN"} }, @{ assign_categories( \@cells, 0, $RNASEQCONFIG->{CN}, "CN" ) } );
            push( @{ $map{"ST"} }, @{ assign_categories( \@cells, 0, $RNASEQCONFIG->{ST}, "ST" ) } );
            push( @{ $map{"PL"} }, @{ assign_categories( \@cells, 0, $RNASEQCONFIG->{PL}, "PL" ) } );
            push( @{ $map{"FILE"} },   @{assign_categories(  \@cells, 1, $RNASEQCONFIG->{FILE},   "FILE"  ) });
            push( @{ $map{"LENGTH"} }, @{ assign_categories( \@cells, 0, $RNASEQCONFIG->{LENGTH}, "LENGTH") } );
            push( @{ $map{"PAIRED"} }, @{ assign_categories( \@cells, 0, $RNASEQCONFIG->{PAIRED}, "PAIRED") } );

            print STDOUT "Sample data:\n";
            $line++;
            next;

        }
        else {
            foreach my $key ( keys %map ) {
                foreach my $col ( @{ $map{$key} } ) {
                    if ( $key eq "FILE" ) {
                        if ( exists $data{$key} ) {
                            $data{$key} .= "/" . $cells[ $col - 1 ];
                        } else {
                            $data{$key} = $cells[ $col - 1 ];
                        }
                    }
                    elsif ($key eq 'PAIRED' or $key eq 'LENGTH') {
                        # no room for whitespace in the paired or length flag
                        $data{$key} .= $cells[ $col - 1 ];
                        $data{$key} =~ s/\s+//g;

                    }
                    else {
                        $data{$key} .= $cells[ $col - 1 ] . " ";
                    }
                    # group together IDs by tissue if specified
                    if ( $key eq $rgt ) {
                        $id_groups{$data{$rgt}}->{$data{'ID'}}++;
                    }
                }
                $data{$key} =~ s/\s*$// if (exists $data{$key});

            }
			if ( $line == 1 ) {
				foreach my $key ( keys %data ) {
					print STDOUT "$key - " . $data{$key} . "\n";
				}
				if (!$jdi) {
					print STDOUT "Continue?(y/n)";
					my $reply = <>;
					chomp $reply;
					exit unless $reply eq "y" or $reply eq "Y";
				}
				$line = 2;
			}
        }
        push @rows, \%data;
    } ## end while (<FILE>)
    close(FILE) or die("Cannot close summary file $summaryfile\n");

    return \@rows;
}

sub check_pipeline_state {
    my ($pipeline_analysis, $sic, $pipeline_state) = @_;

    my $submit_bwa2bam_count     = 0;
    my $bwa2bam_count            = 0;
    my $submit_bam2genes_count   = 0;
    my $bam2genes_count          = 0;
    my $submit_chromosome_count  = 0;
    my $submit_bam2introns_count = 0;
    my $bam2introns_count        = 0;
    my $gsnap_count              = 0;
    my $submit_gsnap_count       = 0;
    foreach my $analysis ( @{ $pipeline_analysis->fetch_all } ) {
        if ( $analysis->logic_name =~ /submit_.+_bwa2bam/ ) {
            $submit_bwa2bam_count += @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) };
        }
        elsif ( $analysis->logic_name =~ /bwa2bam_.+/ ) {
            $bwa2bam_count += @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) };
        }
        elsif ( $analysis->logic_name =~ /submit_chromosome/ ) {
            $submit_chromosome_count += @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) };
        }
        elsif ( $analysis->logic_name =~ /bam2genes/ ) {
            $bam2genes_count += @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) };
        }
        elsif ( $analysis->logic_name =~ /submit_.*bam2introns/ ) {
            $submit_bam2introns_count += @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) };
        }
        elsif ( $analysis->logic_name =~ /bam2introns/ ) {
            $bam2introns_count += @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) };
        }
        elsif ( $analysis->logic_name =~ /submit_.+_gsnap/ ) {
            $submit_gsnap_count += @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) };
        }
        elsif ( $analysis->logic_name =~ /gsnap_.+/ ) {
            $gsnap_count += @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) };
        }
#        foreach my $id ( @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) } ) {
#            $stored_ids->{$id}->{ $analysis->logic_name } = 1;
#        }
    }
    $pipeline_state->{bwa2bam_complete} = 1 if ($submit_bwa2bam_count == $bwa2bam_count);
    $pipeline_state->{bam2genes_complete} = 1 if ($submit_bam2genes_count == $bam2genes_count);
    if ($use_gsnap) {
        $pipeline_state->{introns_complete} = 1 if ($pipeline_state->{bwa2bam_complete} and $submit_gsnap_count == $gsnap_count);
    }
    else {
        $pipeline_state->{introns_complete} = 1 if ($pipeline_state->{bwa2bam_complete} and $submit_bam2genes_count == $bam2genes_count and $submit_bam2introns_count == $bam2introns_count);
    }
#    $pipeline_state->{refine_complete} = 1 if ($pipeline_state->{introns_complete} and
}


sub make_bam2genes_input_ids {
    my ($sic, $pipeline_analysis, $sa) = @_;

    # fetch dummy analysis
    my $analysis = $pipeline_analysis->fetch_by_logic_name("submit_chromosome");
    print "Writing input ids for bam2genes\n";
    # add the submit chromosome input ids to start the next phase of the analysis
    foreach my $slice ( @{ $sa->fetch_all('toplevel') } ) {
        # we dont build on the mitochondrial sequences
        next if ( $slice->seq_region_name eq 'MT' );
        if ( $RNASEQCONFIG->{ALL_PAIRED} ) {
            # run on chromosomes if you have pairing to help separate out the models
            $sic->store_input_id_analysis( $slice->name, $analysis, "dummy" ) unless $stored_ids->{ $slice->name }->{"submit_chromosome"};
            $stored_ids->{$slice->name}->{"submit_chromosome"} = 1;
        }
        else {
            # otherwise run on slices
            my @iid_sections = split( /:/, $slice->name );
            for ( my $i = 1; $i <= $slice->length; $i += $RNASEQCONFIG->{SLICE_LENGTH} ) {
                my $end = $i + $RNASEQCONFIG->{SLICE_LENGTH} - 1;
                $end = $slice->length if $i + $RNASEQCONFIG->{SLICE_LENGTH} - 1 > $slice->length;
                my $id = $iid_sections[0].":".$iid_sections[1].":".$iid_sections[2] .":"."$i:$end:1";
                # print "ID $id \n";
                $sic->store_input_id_analysis( $id, $analysis, "dummy" ) unless (exists $stored_ids->{$id}->{"submit_chromosome"});
                $stored_ids->{$id}->{"submit_chromosome"} = 1;
            }
        }
    }
}

