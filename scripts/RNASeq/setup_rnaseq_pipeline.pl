#!/usr/local/ensembl/bin/perl
# 
# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/scripts/RNASeq/setup_rnaseq_pipeline.pl,v $
# $Revision: 1.23 $
#

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
my $merge_dir = $RNASEQCONFIG->{MERGE_DIR};
$merge_dir = $RNASEQCONFIG->{OUTPUT_DIR} unless $merge_dir;
my $use_gsnap = $RNASEQCONFIG->{USE_GSNAP};
my $gsnap_path = $RNASEQCONFIG->{GSNAP_PATH};
my $rgt = $RNASEQCONFIG->{READ_GROUP_TAG};
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
my $update_analyses;
my $force_stage;

my $usage = "perl setup_rnaseq_pipeline.pl
-verbose    $verbose,
-check      $check, print out which columns are used for which RG tag
-update_analyses $update_analyses, only write the analyses - do not alter the config,
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
            'update_analyses!' => \$update_analyses, );

die($usage) unless ($dbname && $analysisconfigdir &&
$pipelineconfigdir && $delimiter && $summaryfile && $input_dir &&
$output_dir );

if ($force_stage) {
  die($usage)
    unless ( $force_stage == 1 || $force_stage == 2 || $force_stage == 3 );
  $force_stage = 'bwa_complete'       if $force_stage == 1;
  $force_stage = 'gsnap_complete'     if $force_stage == 1 && $use_gsnap;
  $force_stage = 'bam2genes complete' if $force_stage == 2;
  $force_stage = 'configured'         if $force_stage == 3;
  print "Starting from stage $force_stage\n";
}

$stage = $force_stage if $force_stage;
throw("Cannot find input directory $input_dir\n")   unless -e $input_dir;
throw("Cannot find output directory $output_dir\n") unless -e $output_dir;
throw("Cannot find merge directory $merge_dir\n")   unless -e $merge_dir;
system( "mkdir -p " . $merge_dir . "/SAM" ) unless -e $merge_dir . "/SAM";

# get database hash
my %database_hash;
my %databases = %{$DATABASES};
for my $category ( keys %databases ) {
  if ( !$databases{$category}{-host} || !$databases{$category}{-dbname} ) {
    print STDERR "Skipping category "
      . $category
      . " no dbinfo "
      . "defined\n";
    next;
  }
  print STDERR "\n$category: $databases{$category}{-host} "
      . "$databases{$category}{-dbname} :\n--------------------------------------\n";

  my %constructor_args = %{ $databases{$category} };
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor( %constructor_args, );
  $database_hash{$category} = $dba;
}

# test that dbs all exist
my $dba = $database_hash{$dbname};
throw("Db $dbname not found in Databases.pm\n") unless $dba;
throw( "Db " . $RNASEQCONFIG->{ROUGHDB} . " not found in Databases.pm\n" )
  unless $database_hash{ $RNASEQCONFIG->{ROUGHDB} };
throw( "Db " . $RNASEQCONFIG->{REFINEDDB} . " not found in Databases.pm\n" )
  unless $database_hash{ $RNASEQCONFIG->{REFINEDDB} };
throw( "Db " . $RNASEQCONFIG->{BLASTDB} . " not found in Databases.pm\n" )
  unless $database_hash{ $RNASEQCONFIG->{BLASTDB} };

my $sa = $dba->get_SliceAdaptor;
# also need the pipeline adaptor
my %constructor_args = %{$databases{$dbname}};
my $pipelinea =
  new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor( %constructor_args, );
# parse the summary file
# get the pipeline adaptor
# need to delete this hash ref in order 
# to call the pipeline version
my $pipeline_analysis = $pipelinea->get_AnalysisAdaptor;
my $ra                = $pipelinea->get_RuleAdaptor;
my $sic               = $pipelinea->get_StateInfoContainer;
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

my $submit_bwa2bam_count    = 0;
my $bwa2bam_count           = 0;
my $submit_bam2genes_count  = 0;
my $bam2genes_count         = 0;
my $submit_chromosome_count = 0;
my $bam2introns_count       = 0;
my $gsnap_count             = 0;
my $submit_gsnap_count      = 0;
my $rough_count             = 0;
my @files;
# determine the state of the db 
foreach my $analysis ( @{ $pipeline_analysis->fetch_all } ) {
  # print $analysis->logic_name ."\n";
  foreach my $id ( @{ $sic->list_input_ids_by_analysis( $analysis->dbID ) } )
  {
    $stored_ids->{$id}->{ $analysis->logic_name } = 1;
    if ( $analysis->logic_name =~ /submit_.+_bwa2bam/ ) {
      $submit_bwa2bam_count++;
    }
    if ( $analysis->logic_name =~ /bwa2bam_.+/ ) {
      $bwa2bam_count++;
      # remove the file extension
      my @name = split( /\./, $id );
      push @files, shift(@name);
    }
    if ( $analysis->logic_name =~ /submit_chromosome/ ) {
      $submit_chromosome_count++;
    }
    if ( $analysis->logic_name =~ /bam2genes/ ) {
      $bam2genes_count++;
    }
    if ( $analysis->logic_name =~ /bam2introns/ ) {
      $bam2introns_count++
        unless $analysis->logic_name eq 'submit_bam2introns';
    }
    if ( $analysis->logic_name =~ /submit_.+_gsnap/ ) {
      $submit_gsnap_count++;
    }
    if ( $analysis->logic_name =~ /gsnap_.+/ ) {
      $gsnap_count++;
      # remove the : from paired reads
      my @name = split( /\:/, $id );
      push @files, shift(@name);
    }
  } ## end foreach my $id ( @{ $sic->list_input_ids_by_analysis...
} ## end foreach my $analysis ( @{ $pipeline_analysis...

unless ($update_analyses) {
  unless ($use_gsnap) {
    if (    ( $bwa2bam_count >= 1 && $bwa2bam_count == $submit_bwa2bam_count )
         or ( $force_stage eq "bwa_complete" ) )
    {
      $stage = "bwa_complete";
      unless ($check) {
        print "\n\nBWA finished successfully\n------------------------\n";
        # print out the merge command
        print "Before running the next stage of the analysis "
            . "you will need to merge the individual bam files, "
            . "the following commands will do this.\n";
        print "If you wish to remove any of the lanes from the analysis "
            . "at this point just delete them from the merge command.\n\n";
        print "#MERGE\nbsub -o $output_dir" . "/merge.out "
            . "-e $output_dir" . "/merge.err "
            . "-R 'select[mem>5000] rusage[mem=5000]' -M5000000 ";
        print $RNASEQCONFIG->{SAMTOOLS}
            . " merge $merge_dir" . "/merged_unsorted.bam ";
        foreach my $file (@files) {
          print "$output_dir/$file" . "_sorted.bam ";
        }
        print "\n\n";
        print "#SORT\nbsub -o $output_dir" . "/sort.out "
          . "-e $output_dir" . "/sort.err "
          . "-R 'select[mem>5000] rusage[mem=5000]' -M5000000 ";
        print $RNASEQCONFIG->{SAMTOOLS} . " sort $merge_dir"
          . "/merged_unsorted.bam  $merge_dir" . "/merged \n\n";
        print "#INDEX\nbsub -o $output_dir" . "/index.out "
          . "-e $output_dir" . "/index.err "
          . "-R 'select[mem>5000] rusage[mem=5000]' -M5000000 ";
        print $RNASEQCONFIG->{SAMTOOLS} . " index $merge_dir" . "/merged.bam\n\n";
      } ## end unless ($check)
    } ## end if ( ( $bwa2bam_count ...
  } else {
    if (    ( $gsnap_count >= 1 && $gsnap_count == $submit_gsnap_count )
         or ( $force_stage eq "bwa_complete" ) )
    {
      $stage = "gsnap_complete";
      unless ($check) {
        print "\n\nGSNAP finished successfully\n------------------------\n";
        # print out the merge command
        print "Before running the next stage of the analysis "
            . "you will need to merge the individual bam files, "
            . "the following commands will do this.\n";
        print "If you wish to remove any of the lanes from the analysis "
            . "at this point just delete them from the merge command.\n\n";
        print "#MERGE\nbsub -o $output_dir" . "/merge.out "
          . "-e $output_dir" . "/merge.err "
          . "-R 'select[mem>5000] rusage[mem=5000]' -M5000000 ";
        print $RNASEQCONFIG->{SAMTOOLS} . " merge $merge_dir" . "/merged_unsorted.bam ";
        foreach my $file (@files) {
          print "$output_dir/$file" . "_sorted.bam ";
        }
        print "\n\n";
        print "#SORT\nbsub -o $output_dir" . "/sort.out "
          . "-e $output_dir" . "/sort.err "
          . "-R 'select[mem>5000] rusage[mem=5000]' -M5000000 ";
        print $RNASEQCONFIG->{SAMTOOLS}
          . " sort $merge_dir" . "/merged_unsorted.bam  $merge_dir" . "/merged \n\n";
        print "#INDEX\nbsub -o $output_dir" . "/index.out "
          . "-e $output_dir" . "/index.err "
          . "-R 'select[mem>5000] rusage[mem=5000]' -M5000000 ";
        print $RNASEQCONFIG->{SAMTOOLS} . " index $merge_dir" . "/merged.bam\n\n";
      } ## end unless ($check)
    } ## end if ( ( $gsnap_count >=...
  } ## end else
  if ( $bam2genes_count > 0 ) {
    $stage = "bam2genes";
  }
  unless ( $stage eq "Initialization" ) {
    my $slice_count = 0;
    if ( $pipeline_analysis->fetch_by_logic_name("submit_chromosome") ) {
      my $slice_count =
        scalar(@{ $sic->list_input_ids_by_analysis($pipeline_analysis->fetch_by_logic_name("submit_chromosome")->dbID)});
    }
    if ( $pipeline_analysis->fetch_by_logic_name("submit_bam2introns") ) {
      $rough_count =
        scalar(@{ $sic->list_input_ids_by_analysis($pipeline_analysis->fetch_by_logic_name("submit_bam2introns")->dbID)});
    }
    if ( (    $submit_chromosome_count > 0
           && $submit_chromosome_count == $bam2genes_count
           && $bam2introns_count == 0 )
         or ( $force_stage eq "bam2genes complete" ) )
    {
      $stage = 'bam2genes complete';
      my $analysis;
      unless ($use_gsnap) {
        $analysis =
          $pipeline_analysis->fetch_by_logic_name("submit_bam2introns");
        unless ($analysis) {
          throw("submit_bam2introns analysis not found\n");
        }
      }
  # assign stable ids to the models and make input ids for the bam2introns run
      my $rough_db = $database_hash{ $RNASEQCONFIG->{ROUGHDB} };
      my $sql = "update gene set stable_id =  concat( 'RNASEQG',lpad(gene_id,11,'0'))";
      my $sth = $rough_db->dbc->prepare($sql);
      $sth->execute();
      $sql = "update transcript set stable_id =  concat( 'RNASEQT',lpad(gene_id,11,'0'))";
      $sth = $rough_db->dbc->prepare($sql);
      $sth->execute();
      # get the input ids
      $sql = "select stable_id from gene";
      $sth = $rough_db->dbc->prepare($sql);
      $sth->execute();
      unless ($use_gsnap) {
        my @ids;
        while ( my @array = $sth->fetchrow_array ) {
          push @ids, $array[0];
        }
        foreach my $id (@ids) {
          $sic->store_input_id_analysis( $id, $analysis, "dummy" )
            unless $stored_ids->{$id}->{"submit_bam2introns"};
        }
        print "Stored " . scalar(@ids) . " submit_bam2introns ids\n";
        $rough_count = scalar(@ids);
      }
    } else {
      if ( $stage eq 'bwa_complete' or $stage eq 'gsnap_complete' ) {
        my $count = 0;
        # fetch dummy analysis
        my $analysis = $pipeline_analysis->fetch_by_logic_name("submit_chromosome");
        if ($analysis) {
          print "Writing input ids for bam2genes\n";
          # add the submit chromosome input ids to start the next phase of the analysis
          foreach my $slice ( @{ $sa->fetch_all('toplevel') } ) {
            # we dont build on the mitochondrial sequences
            next if ( $slice->seq_region_name eq 'MT' );
            unless ( $RNASEQCONFIG->{SLICE_LENGTH} ) {
              $count++;
              # run on chromosomes if you have pairing to help separate out the models
              $sic->store_input_id_analysis( $slice->name, $analysis, "dummy" )
                unless $stored_ids->{ $slice->name }->{"submit_chromosome"};
            } else {
              # otherwise run on slices
              my @iid_sections = split( /:/, $slice->name );
              for ( my $i = 1 ;
                    $i <= $slice->length ;
                    $i += $RNASEQCONFIG->{SLICE_LENGTH} )
              {
                my $end = $i + $RNASEQCONFIG->{SLICE_LENGTH} - 1;
                $end = $slice->length
                  if $i + $RNASEQCONFIG->{SLICE_LENGTH} - 1 > $slice->length;
                my $id =
                    $iid_sections[0] . ":"
                  . $iid_sections[1] . ":"
                  . $iid_sections[2] . ":"
                  . "$i:$end:1";
                # print "ID $id \n";
                $count++;
                $sic->store_input_id_analysis( $id, $analysis, "dummy" )
                  unless $stored_ids->{$id}->{"submit_chromosome"};
              }
            }
          } ## end foreach my $slice ( @{ $sa->fetch_all...
        } else {
          print "Cannot find submit_chromosome analysis - "
              . "run the script with the -write option to refresh the analyses\n";
        }
        $slice_count = $count;
      } ## end if ( $stage eq 'bwa_complete'...
    } ## end else [ if ( ( $submit_chromosome_count...

    # calculate batch sizes for slice jobs
    $slice_batches = int( $slice_count/100 );
    $slice_batches = 1 unless $slice_batches;
    $slice_batches = 500 if $slice_batches > 500;
    $rough_batches = int( $rough_count/100 );
    $rough_batches = 1 unless $rough_batches;
    $rough_batches = 500 if $rough_batches > 500;
    print "Got $slice_count slices using a batch size of $slice_batches \n"
      if ($stage eq 'bwa_complete') or ($stage eq 'gsnap_complete');
    print "Got $rough_count (batches of) rough models using a batch size of $rough_batches \n";

    if (    ( $bam2introns_count > 0 )
         or ( $force_stage eq "configured" )
         or ( $use_gsnap && $stage eq 'bam2genes complete' ) )
    {
      $stage = 'configured';
    }
  } ## end unless ( $stage eq "Initialization")
} ## end unless ($update_analyses)

print "Stage = $stage\n";

open( FILE, $summaryfile ) or die("Cannot open summary file $summaryfile\n");
my $line = 0;
my @rows;
my %map;
while (<FILE>) {
  my %data;
  chomp;
  next if $_ =~ /^#/;
  my @cells = split( $delimiter, $_ );
  if ( $line == 0 ) {
    # header row
    print STDERR "Got these headers\n"
      if ( $stage eq "Initialization" ) || $check;
    for ( my $i = 0 ; $i <= $#cells ; $i++ ) {
      my $header = $cells[$i];
      print STDERR $i + 1 . ") $header\n"
        if ( $stage eq "Initialization" ) || $check;
    }
    
    print STDERR "Please assign some / all columns to the some / all "
        . "of the following categories multiple values can be separted with commas:\n"
      if ( $stage eq "Initialization" ) || $check;
    print STDERR "Tag\tDescription\n" if ( $stage eq "Initialization" ) || $check;

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

    print STDERR "Sample data:\n" if ( $stage eq "Initialization" ) || $check;
    $line++;
    next;

  } else {
    foreach my $key ( keys %map ) {
      foreach my $col ( @{ $map{$key} } ) {
        unless ( $key eq "FILE" ) {
          $data{$key} .= $cells[ $col - 1 ] . " ";
        } else {
          if ( $data{$key} ) {
            $data{$key} .= "/" . $cells[ $col - 1 ];
          } else {
            $data{$key} = $cells[ $col - 1 ];
          }
        }
        # no room for whitespace in the paired or length flag
        $data{$key} =~ s/\s+//g if $key eq 'PAIRED' or $key eq 'LENGTH';
	print "$key  - " . $cells[ $col-1] ."\n";
	# group together IDs by tissue if specified
	if ( $key eq $rgt ) {
	  $id_groups{$data{$rgt}}->{$data{'ID'}}++;
	}
      }
    }
  }
  # add the path to the file name
  # $data{FILE} = $input_dir."/".$data{FILE};
  if ( $stage eq "Initialization" ) {
    if ( $line == 1 ) {
      foreach my $key ( keys %data ) {
        print STDERR "$key - " . $data{$key} . "\n";
      }
      print STDERR "Continue?(y/n)";
      my $reply = <>;
      chomp $reply;
      exit unless $reply eq "y" or $reply eq "Y";
    }
  }
  push @rows, \%data;
  $line++;
} ## end while (<FILE>)
if ( $stage eq "Initialization" ) {
  print STDERR "Processed $line lines of data \n";
  print STDERR "Using " . $dba->dbc->dbname . "@" . $dba->dbc->host . "  as pipeline db\n";
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
 $line = 0;
foreach my $row (@rows) {
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
                                  -program_file => "/software/solexa/bin/bwa",
                                  -module       => "BWA",
                                  -description  => $row->{DS},
                                  -display_label => $row->{ID},
                                  -displayable   => '1',
                                  -input_id_type => 'BWA' . $ln, );

  my $bwa2bam =
    new Bio::EnsEMBL::Pipeline::Analysis(
                                  -logic_name   => "bwa2bam_" . $ln,
                                  -program      => "bwa",
                                  -program_file => "/software/solexa/bin/bwa",
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
  foreach my $key2 ( keys %{ $pairs{$key} } ) {
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
my $rnaseq_blast =
  new Bio::EnsEMBL::Pipeline::Analysis(
                                  -logic_name    => "rnaseqblast",
                                  -input_id_type => 'CHROMOSOME',
                                  -module        => 'BlastRNASeqPep',
                                  -paramteres => '-cpus => 1, -hitdist => 40',
                                  -program_file => 'wublastp',
                                  -program      => 'wublastp',
                                  -db_file      => $RNASEQCONFIG->{UNIPROTDB},
  );

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
$pipeline_analysis->store($bam2genes) if $stage eq 'bwa_complete' or  $stage eq 'gsnap_complete';

unless ($use_gsnap) {
  $pipeline_analysis->store($submit_bam2introns);
  $pipeline_analysis->store($bam2introns)       if $stage eq 'bam2genes complete';
  $pipeline_analysis->store($bam2introns_wait)  if $stage eq 'bam2genes complete';
  $pipeline_analysis->store($submit_sam2bam)    if $stage eq 'configured' or $stage eq 'bam2genes complete';
  $pipeline_analysis->store($sam2bam)           if $stage eq 'configured';
  $pipeline_analysis->store($sam2bam_wait)      if $stage eq 'configured';
}

$pipeline_analysis->store($refine_all)   if $stage eq 'configured';
$pipeline_analysis->store($rnaseq_blast) if $stage eq 'configured';

if ( $stage eq 'bwa_complete' or $stage eq 'gsnap_complete' ) {
  $ra->store($bam2genes_rule) if check_rule($bam2genes_rule);
}
unless ($use_gsnap) {
  if ( $stage eq 'bam2genes complete' ) {
    $ra->store($bam2introns_rule)      if check_rule($bam2introns_rule);
    $ra->store($bam2introns_wait_rule) if check_rule($bam2introns_wait_rule);
  }
}
if ( $stage eq 'configured' ) {
  unless ($use_gsnap) {
    $ra->store($sam2bam_rule)      if check_rule($sam2bam_rule);
    $ra->store($sam2bam_wait_rule) if check_rule($sam2bam_wait_rule);
  }
  $ra->store($refine_all_rule)  if check_rule($refine_all_rule);
  $ra->store($rnaseqblast_rule) if check_rule($rnaseqblast_rule);
}
unless ($use_gsnap) {
  # need to add a dummy input id for submit_sam2bam
  if ( $stage eq 'configured' or $stage eq 'bam2genes complete' ) {
    $sic->store_input_id_analysis( 'dummy', $submit_sam2bam, "dummy" )
      unless $stored_ids->{'dummy'}->{'submit_sam2bam'};
  }
}

exit if $update_analyses;

# check config directory
# write the extra information header files
open( ALL, ">$output_dir/all_headers.txt" )
  or die("Cannot open  $output_dir/all_headers.txt for writing\n");

my %seen_it;
foreach my $row (@rows) {
  next if $seen_it{ $row->{ID} };
  $seen_it{ $row->{ID} } = 1;
  open( HEAD, ">$output_dir/" . $row->{ID} . "_header.txt" )
    or die( "Cannot open  $output_dir/" . $row->{ID} . "_header.txt for writing\n" );

  print HEAD "\@RG\tID:" . $row->{ID} . "\tPU:" . $row->{PU} . "\tSM:" . $row->{SM} . "\t";
  print HEAD "LB:"       . $row->{LB} . "\tDS:" . $row->{DS} . "\tCN:" . $row->{CN} . "\t";
  print HEAD "ST:"       . $row->{ST} . "\tPL:" . $row->{PL} . "\n";
  print ALL "\@RG\tID:"  . $row->{ID} . "\tPU:" . $row->{PU} . "\tSM:" . $row->{SM} . "\t";
  print ALL "LB:"        . $row->{LB} . "\tDS:" . $row->{DS} . "\tCN:" . $row->{CN} . "\t";
  print ALL "ST:"        . $row->{ST} . "\tPL:" . $row->{PL} . "\n";
} ## end foreach my $row (@rows)

if ($stage eq "Initialization" || $check) {
  unless ( $use_gsnap ) {
    print STDERR "Have these config files to modify:
$analysisconfigdir/GeneBuild/BWA.pm
$analysisconfigdir/GeneBuild/Bam2Genes.pm
$analysisconfigdir/GeneBuild/Bam2Introns.pm
$analysisconfigdir/GeneBuild/Sam2Bam.pm
$analysisconfigdir/GeneBuild/RefineSolexaGenes.pm
$analysisconfigdir/GeneBuild/BlastRNASeqPep.pm
$pipelineconfigdir/BatchQueue.pm
 - backing them up\n" ;
  } else {
    print STDERR "Have these config files to modify:
$analysisconfigdir/GeneBuild/Gsnap.pm
$analysisconfigdir/GeneBuild/Bam2Genes.pm
$analysisconfigdir/GeneBuild/RefineSolexaGenes.pm
$analysisconfigdir/GeneBuild/BlastRNASeqPep.pm
$pipelineconfigdir/BatchQueue.pm
 - backing them up\n" ;
  }
}

system("mv $pipelineconfigdir/BatchQueue.pm $pipelineconfigdir/BatchQueue.pm_bk")
  if -e "$pipelineconfigdir/BatchQueue.pm" &
    !-e "$pipelineconfigdir/BatchQueue.pm_bk";

system("mv $analysisconfigdir/GeneBuild/BWA.pm $analysisconfigdir/GeneBuild/BWA.pm_bk")
  if -e "$analysisconfigdir/GeneBuild/BWA.pm" &
    !-e "$analysisconfigdir/BWA.pm_bk";

system("mv $analysisconfigdir/GeneBuild/Bam2Genes.pm $analysisconfigdir/GeneBuild/Bam2Genes.pm_bk")
  if -e "$analysisconfigdir/GeneBuild/Bam2Genes.pm" &
    !-e "$analysisconfigdir/GeneBuild/Bam2Genes.pm_bk";

system("mv $analysisconfigdir/GeneBuild/Bam2Introns.pm $analysisconfigdir/GeneBuild/Bam2Introns.pm_bk")
  if -e "$analysisconfigdir/GeneBuild/Bam2Introns.pm" &
    !-e "$analysisconfigdir/GeneBuild/Bam2Introns.pm_bk";

system("mv $analysisconfigdir/GeneBuild/Sam2Bam.pm $analysisconfigdir/GeneBuild/Sam2Bam.pm_bk")
  if -e "$analysisconfigdir/GeneBuild/Sam2Bam.pm" &
    !-e "$analysisconfigdir/GeneBuild/Sam2Bam.pm_bk";

system("mv $analysisconfigdir/GeneBuild/RefineSolexaGenes.pm $analysisconfigdir/GeneBuild/RefineSolexaGenes.pm_bk")
  if -e "$analysisconfigdir/GeneBuild/RefineSolexaGenes.pm" &
    !-e "$analysisconfigdir/GeneBuild/RefineSolexaGenes.pm_bk";

system("mv $analysisconfigdir/GeneBuild/BlastRNASeqPep.pm $analysisconfigdir/GeneBuild/BlastRNASeqPep.pm_bk")
  if -e "$analysisconfigdir/GeneBuild/BlastRNASeqPep.pm" &
    !-e "$analysisconfigdir/GeneBuild/BlastRNASeqPep.pm_bk";

system("mv $analysisconfigdir/GeneBuild/Gsnap.pm $analysisconfigdir/GeneBuild/Gsnap.pm_bk")
  if -e "$analysisconfigdir/GeneBuild/Gsnap.pm" &
    !-e "$analysisconfigdir/GeneBuild/Gsnap.pm_bk";

system("mv $analysisconfigdir/Blast.pm $analysisconfigdir/Blast.pm_bk")
  if -e "$analysisconfigdir/Blast.pm" & !-e "$analysisconfigdir/Blast.pm_bk";

# write config
open( BATCHQUEUE, ">$pipelineconfigdir/BatchQueue.pm" )
  or die( "Cannot open " . $pipelineconfigdir . "/BatchQueue.pm  for writing\n" );
open( BAM2GENES, ">$analysisconfigdir/GeneBuild/Bam2Genes.pm" )
  or die( "Cannot open " . $analysisconfigdir . "/GeneBuild/Bam2Genes.pm  for writing\n" );
open( REFINE, ">$analysisconfigdir/GeneBuild/RefineSolexaGenes.pm" )
  or die( "Cannot open " . $analysisconfigdir . "/GeneBuild/RefineSolexaGenes.pm  for writing\n" );
open( BLAST, ">$analysisconfigdir/GeneBuild/BlastRNASeqPep.pm" )
  or die( "Cannot open " . $analysisconfigdir . "/GeneBuild/BlastRNASeqPep.pm  for writing\n" );
open( BLAST2, ">$analysisconfigdir/Blast.pm" )
  or die( "Cannot open " . $analysisconfigdir . "/Blast.pm  for writing\n" );
unless ($use_gsnap) {
  open( BWA, ">$analysisconfigdir/GeneBuild/BWA.pm" )
    or die( "Cannot open " . $analysisconfigdir . "/GeneBuild/BWA.pm  for writing\n" );
  open( BAM2INTRONS, ">$analysisconfigdir/GeneBuild/Bam2Introns.pm" )
    or die( "Cannot open " . $analysisconfigdir . "/GeneBuild/Bam2Introns.pm  for writing\n" );
  open( SAM2BAM, ">$analysisconfigdir/GeneBuild/Sam2Bam.pm" )
    or die( "Cannot open " . $analysisconfigdir . "/GeneBuild/Sam2Bam.pm  for writing\n" );
} else {
  open( GSNAP, ">$analysisconfigdir/GeneBuild/Gsnap.pm" )
    or die( "Cannot open " . $analysisconfigdir . "/GeneBuild/Gsnap.pm  for writing\n" );
}

# write headers
my $string;
unless ( $use_gsnap ) {
  $string =  bwa_header();
  print BWA $string;
  $string =  BAM2INTRONS_header();
  print BAM2INTRONS $string;
  print BAM2INTRONS tail();
  $string =  SAM2BAM_header();
  print SAM2BAM $string;
  print SAM2BAM tail();
} else {
  $string =  GSNAP_header();
  print GSNAP $string;
}

# batchqueue
$string =  batchqueue_header();
print BATCHQUEUE $string;
$string =  BAM2GENES_header();
print BAM2GENES $string;
print BAM2GENES tail();
$string =  REFINE_header();
print REFINE $string;
$string =  BLAST_header();
print BLAST $string;
print BLAST tail();
$string =  blast2_header();
print BLAST2 $string;
print BLAST2 tail();

#write logic name specific config
my %seen;
foreach my $row ( @rows ) {
  my $length = $RNASEQCONFIG->{READ_LENGTH};
  $length = $row->{LENGTH} if $row->{LENGTH};
  next if $seen{$row->{ID}};
  unless ( $use_gsnap ) {
    # bwa
    print BWA "           'bwa_" . $row->{ID} ."' => {
                        INDIR   => \"". $input_dir  ."\",
                        OUTDIR  => \"". $output_dir ."\",
                        OPTIONS => \"-n " . int($length / 2). " -i " . $length  ."\",
                      },\n";

    # bwa2bam
    print BWA "           'bwa2bam_" . $row->{ID} ."' => {
                        HEADER  => \"$output_dir/" . $row->{ID} ."_header.txt\",
                        INDIR   => \"". $input_dir  ."\",
                        OUTDIR  => \"". $output_dir ."\",\n";
    print BWA "                        PAIRED => 1,\n" if $row->{PAIRED};
    print BWA "                 },\n";
  } else {
    print GSNAP "           'gsnap_" . $row->{ID} ."' => {
                        INDIR   => \"". $input_dir  ."\",
                        OUTDIR  => \"". $output_dir ."\",
                        HEADER  => \"$output_dir/" . $row->{ID} ."_header.txt\",\n";
    print GSNAP "                        PAIRED => 1,\n" if $row->{PAIRED};
    print GSNAP "                 },\n";
  }
  # make the single tissue config for refine
  my $str;
  unless ( $use_gsnap ) {
$str .= '  			INTRON_BAM_FILES => [
				{
		    		# location of the bam file
				FILE => "'.$RNASEQCONFIG->{MERGE_DIR}.'/introns.bam",
				# does the bam file(s) contain a mixture of spliced and unspliced reads such as you might get using Gsnap?
				MIXED_BAM => "0",
				# only take introns above this depth
				DEPTH => "0",
				# return only reads from the specified read groups
	               		GROUPNAME => ["' .$ids_by_tissue{$tissue_by_id{$row->{ID}}}  . '"],
				},
			],
';
} else {
$str .= '  			INTRON_BAM_FILES => [
				{
		    		# location of the bam file
				FILE => "'.$RNASEQCONFIG->{MERGE_DIR}.'/merged.bam",
				# does the bam file(s) contain a mixture of spliced and unspliced reads such as you might get using Gsnap?
				MIXED_BAM => "1",
				# only take introns above this depth
				DEPTH => "0",
				# return only reads from the specified read groups
	               		GROUPNAME => ["' .$ids_by_tissue{$tissue_by_id{$row->{ID}}} . '"],
				},
			], 
';}
  print REFINE "             'refine_" . $tissue_by_id{$row->{ID}} ."' => {
$str
			SINGLE_EXON_MODEL => ''},\n" if  $RNASEQCONFIG->{SINGLE_TISSUE} &! $seen{$tissue_by_id{$row->{ID}}};  
												

  # batchqueue
  print BATCHQUEUE "       {
             logic_name => 'bwa_" . $row->{ID} ."',
             output_dir => '".$output_dir ."/".$row->{ID}."_pipeline',
             memory    => [ '5GB', '10GB', '20GB','30GB' ],
             resource  => 'select[myens_".$ref_load."tok>800] ' .  'rusage[myens_".$ref_load."tok=25]',
       },
       {
             logic_name => 'bwa_" . $row->{ID} ."_wait',
             output_dir => '".$output_dir ."/".$row->{ID}."_pipeline',
             resource  => 'select[myens_".$ref_load."tok>800] ' .  'rusage[myens_".$ref_load."tok=25]',
       },
       {
             logic_name => 'bwa2bam_" . $row->{ID} ."',
             output_dir => '".$output_dir ."/".$row->{ID}."_pipeline',
             memory    => [ '2GB', '5GB', '10GB', '20GB', '30GB' ],
             resource  => 'select[myens_".$ref_load."tok>800] ' .  'rusage[myens_".$ref_load."tok=25]',
       },
       {
             logic_name => 'gsnap_" . $row->{ID} ."',
             output_dir => '".$output_dir ."/".$row->{ID}."_pipeline',
             memory    => [ '10GB', '20GB', '30GB' ],
             queue     => 'long',
             resource  => 'select[myens_".$ref_load."tok>800] ' .  'rusage[myens_".$ref_load."tok=25]',
       },
";


 print BATCHQUEUE "       {
             logic_name => 'refine_" . $tissue_by_id{$row->{ID}} ."',
             output_dir => '".$output_dir ."/refine_".$row->{ID}."_pipeline',
             batch_size => ".$slice_batches.",
             lsf_perl  => '/usr/local/bin/perl',
             memory    => [ '1GB', '2GB', '5GB', '15GB','30GB' ],
             resource  => 'select[myens_".$ref_load."tok>800] && select[myens_".$refine_load."tok>800] ' .  'rusage[myens_".$ref_load."tok=25:myens_".$refine_load."tok=25]',
       },
" if  $RNASEQCONFIG->{SINGLE_TISSUE} &! $seen{$tissue_by_id{$row->{ID}}};  ;
# only write each tissue once
$seen{$tissue_by_id{$row->{ID}}} = 1;

  $seen{$row->{ID}} = 1;
}
unless ( $use_gsnap ) {
print BWA "     }\n";
print BWA tail();
} else {
print GSNAP "     }\n";
print GSNAP tail();
}
print BATCHQUEUE "     ]\n";
print REFINE "     }\n";

print BATCHQUEUE tail();
print REFINE tail();

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
    developers list at <dev@ensembl.org>.

    Questions may also be sent to the Ensembl help desk at
    <helpdesk@ensembl.org>.

=cut


package Bio::EnsEMBL::Pipeline::Config::BatchQueue;

use strict;
use vars qw(%Config);

%Config = (

  # Depending on the job-submission-system you"re using,
  # use LSF, you can also use "Local".
  #
  # For more info look into:
  # /ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/BatchSubmission
  QUEUE_MANAGER => "LSF", # use "SGE_GridEngine_v6" for running in ensembl cloud evironment

  DEFAULT_BATCH_SIZE  => 1,
  DEFAULT_RETRIES     => 5,
  DEFAULT_BATCH_QUEUE => "normal",    # Put in the queue of your choice, eg. "normal"
  DEFAULT_RESOURCE    => "",
  DEFAULT_SUB_ARGS    => "",
  DEFAULT_OUTPUT_DIR  => "",
  DEFAULT_CLEANUP     => "no",
  DEFAULT_VERBOSITY   => "WARNING",

  # The two variables below are to overcome a bug in LSF. 
  # We"re currently running the pre-exec with a different perl. lsf currently unsets the LD_LIBRARY_PATH 
  # which we need for certain 64bit libraries in pre-exec commands. (more info see LSF_LD_SECURITY variable ) 

  DEFAULT_LSF_PRE_EXEC_PERL =>"/usr/local/ensembl32/bin/perl", # ONLY use 32bit perl for lsf -pre-exec jobs
  DEFAULT_LSF_PERL =>"/usr/local/ensembl32/bin/perl", # ONLY use ensembl64/bin/perl for memory jobs > 4 gb
                                                     # SANGER farm : don"t forget to source source /software/intel_cce_80/bin/iccvars.csh for big mem jobs 
                                                     #
  # At this number of jobs RuleManager will sleep for a certain period
  # of time if you effectively want this never to run set the value to
  # very high ie 100000 for a certain period of time this is important
  # for queue managers which cannot cope with large numbers of pending
  # jobs (e.g. early LSF versions and SGE)
  JOB_LIMIT           => 10000,

  JOB_STATUSES_TO_COUNT => ["PEND"],    # These are the jobs which will be 
                                        # counted. valid statuses for
                                        # this array are RUN, PEND, SSUSP, EXIT, DONE ; use "qw" for Sun Grid Engine
  MARK_AWOL_JOBS => 1,
  MAX_JOB_SLEEP  => 3600,   # The maximun time to sleep for when job limit
                            # reached
  MIN_JOB_SLEEP => 120, # The minimum time to sleep for when job limit reached
  SLEEP_PER_JOB => 30,  # The amount of time to sleep per job when job limit
                        # reached

  DEFAULT_RUNNABLEDB_PATH => "Bio/EnsEMBL/Analysis/RunnableDB",

  DEFAULT_RUNNER         => "",
  DEFAULT_RETRY_QUEUE    => "long",
  DEFAULT_RETRY_SUB_ARGS    => "",
  DEFAULT_RETRY_RESOURCE => "",

  QUEUE_CONFIG => [
       {
             logic_name => "refine_all",
             output_dir => "'.$output_dir .'/refine_all_pipeline",
             batch_size => '.$slice_batches.',
             lsf_perl  => "/usr/local/bin/perl",
             memory    => [ "2GB", "5GB", "10GB", "20GB","30GB" ],
             resource  => \'select[myens_'.$ref_load.'tok>800] && select[myens_'.$refine_load.'tok>800] \' .  \'rusage[myens_'.$ref_load.'tok=25:myens_'.$refine_load.'tok=25]\',
       },
       {
             logic_name => "bam2genes",
             output_dir => "'.$output_dir .'/bam2genes_pipeline",
             batch_size => '.$slice_batches.',
             lsf_perl  => "/usr/local/bin/perl",
             memory    => [ "2GB", "5GB", "10GB", "20GB","30GB" ],
             resource  => \'select[myens_'.$ref_load.'tok>800] && select[myens_'.$rough_load.'tok>800] \' .  \'rusage[myens_'.$ref_load.'tok=25:myens_'.$rough_load.'tok=25]\',
       },
       {
             logic_name => "bam2genes_wait",
             output_dir => "'.$output_dir .'/bam2genes_wait_pipeline",
       },
       {
             logic_name => "bam2introns",
             output_dir => "'.$output_dir .'/bam2introns_pipeline",
             batch_size => '.$rough_batches.',
             lsf_perl  => "/usr/local/bin/perl",
             memory    => [ "2GB", "5GB", "10GB", "20GB","30GB" ],
             resource  => \'select[myens_'.$ref_load.'tok>800] \' .  \'rusage[myens_'.$ref_load.'tok=25]\',
       },
       {
             logic_name => "bam2introns_wait",
             output_dir => "'.$output_dir .'/bam2introns_wait_pipeline",
       },
       {
             logic_name => "sam2bam",
             output_dir => "'.$output_dir .'/sam2bam_pipeline",
             memory    => [ "2GB", "5GB", "10GB", "20GB","30GB" ],
             resource  => \'select[myens_'.$ref_load.'tok>800] \' .  \'rusage[myens_'.$ref_load.'tok=25]\',
       },
       {
             logic_name => "sam2bam_wait",
             output_dir => "'.$output_dir .'/sam2bam_wait_pipeline",
       },
       {
             logic_name => "rnaseqblast",
             output_dir => "'.$output_dir .'/rnaseqblast_pipeline",
             batch_size => '.$slice_batches.',
             memory    => [ "2GB", "5GB", "10GB", "20GB","30GB" ],
             resource  => \'select[myens_'.$ref_load.'tok>800] && select[myens_'.$blast_load.'tok>800] \' .  \'rusage[myens_'.$ref_load.'tok=25:myens_'.$blast_load.'tok=25]\',
       },
';
  return $str;
}

sub blast2_header {
  my $str = '# Ensembl module for Bio::EnsEMBL::Analysis::Config::Blast
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Config::Blast

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Config::Blast;
  
  use Bio::EnsEMBL::Analysis::Config::Blast qw(BLAST_CONTIG);
=head1 CONTACT

Post questions to the Ensembl development list: dev@ensembl.org

=cut

package Bio::EnsEMBL::Analysis::Config::Blast;

use strict;
use vars qw(%Config);


%Config = (
  BLASTDB => "/data/blastdb/Ensembl/Uniprot/", 
  BLAST_CONFIG => {
    DEFAULT => {
      BLAST_PARSER => "Bio::EnsEMBL::Analysis::Tools::BPliteWrapper",
      PARSER_PARAMS => {
        -regex => \'^(\w+)\',
        -query_type => undef,
        -database_type => undef,
      },
      BLAST_FILTER => undef,
      FILTER_PARAMS => {},
      BLAST_PARAMS => {
        -unknown_error_string => "FAILED",
        -type => "wu",
      }
    },
    rnaseqblast =>
    {
     BLAST_PARSER => "Bio::EnsEMBL::Analysis::Tools::FilterBPlite",
     PARSER_PARAMS => {
                       -regex => \'^(\w+\W\d+)\',
                       -query_type => "pep",
                       -database_type => "pep",
                       -threshold_type => "PVALUE",
                       -threshold => 0.01,
                      },
     BLAST_FILTER => "Bio::EnsEMBL::Analysis::Tools::FeatureFilter",
     FILTER_PARAMS => {
                       -min_score => 200,
                       -prune => 1,
                      },
     },

  },
  BLAST_AB_INITIO_LOGICNAME => ["genscan"],
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
# Cared for by EnsEMBL (dev@ensembl.org)
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
          GENOMEFILE => "'. $RNASEQCONFIG->{GENOME_DIR}.'/'.$RNASEQCONFIG->{GENOME_FILE} .'",

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

sub BAM2GENES_header {
  my $str = '# Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Genes;
#
# Cared for by EnsEMBL (dev@ensembl.org)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code


package Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Genes;

use strict;
use vars qw( %Config );

# Hash containing config info
%Config = (

BAM2GENES_CONFIG_BY_LOGIC =>
     {
      DEFAULT => {
                  # databases are defined as hash keys from Bio::EnsEMBL::Analysis::Config::Databases
                  OUTPUT_DB    => "'.$RNASEQCONFIG->{ROUGHDB}.'",

          # location of sorted and indexed bam file containing genomic alignments
              ALIGNMENT_BAM_FILE => "'.$RNASEQCONFIG->{MERGE_DIR}.'/merged.bam",

          # logic_name for repeats will be used to merge exons separated by repeats
          # leave blank if you don"t want to fill in the gaps ( shouldn"t be needed with BWA anyway )
                  REPEAT_LN  => "",

                  # options for filtering out small gene models
                  MIN_LENGTH => 300,
                  MIN_EXONS  =>   1,

          # are you using paired reads?
          PAIRED => '.$RNASEQCONFIG->{ALL_PAIRED}.',
               # If using paired reads we need to remove the 1 or 2 tag from the end of the paired read names so
          # we can pair them by name, this regex will usually work but if you have
          # reads with differently structured headers you may need to tweak it
          PAIRING_REGEX => \'_\d\',

          # if we are using unpaired reads we use the max intron length to
          # join clusters of reads into transcripts
          MAX_INTRON_LENGTH => 200000,

          # genes with a span < MIN_SPAN will also be considered single exon
          MIN_SINGLE_EXON_LENGTH => 1000,

                  # "span" = genomic extent / cdna length
                  MIN_SPAN   =>   1.5,
                  },
          bam2genes =>{},
     }
';
  return $str;
}

sub BAM2INTRONS_header {
  my $str = '# package Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Introns
#
# Cared for by EnsEMBL (dev@ensembl.org)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Introns

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Introns

=head1 DESCRIPTION

This contains the specific configuraton for 
Bio::EnsEMBL::Analysis::RunnableDB::Bam2Introns and 
Bio::EnsEMBL::Analysis::RunnableDB::Bam2IntronsTranscript

=head1 CONTACT

=cut


package Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Introns;

use strict;
use vars qw( %Config );

%Config = (
  BAM2INTRONS_CONFIG_BY_LOGIC => {
    DEFAULT => {
      ##############################################
      # Write out the alignments in SAM / BAM format
      # specify the path to an output directory here
      # files will be created as input_id.sam
      OUT_SAM_DIR => "'.$RNASEQCONFIG->{MERGE_DIR}.'/SAM",

      # dont allow more then X % missmatches ie a number of 6% = 2
      # missmatches on a 35 bp read and 4 missmatches on a 75 bp read
      # etc..
      MISSMATCH => 6,

      # Database to fetch the trancripts from
      TRANSDB => "'.$RNASEQCONFIG->{ROUGHDB}.'",

      # Loaction of BAM file containg the genomic alignments\
      BAM_FILE => "'.$RNASEQCONFIG->{MERGE_DIR}.'/merged.bam",

      # Exonerate word length, smaller = more accurate takes longer
      WORD_LENGTH => "10",

      # Exonerate saturate threshold - smaller = quicker but you lose
      # significant numbers of alignments where there is high depth
      # of reads, best to leave it undefined but if jobs are taking
      # forever try a value like 1000
      SATURATE_THRESHOLD => "10000",

      # repeat masks the transcript sequences - quick but you might
      # miss something
      MASK => "1",

      # Filter the alignments based on percent_id and coverage
      PERCENT_ID => 97,
      COVERAGE   => 90,

      # use the full genomic sequence rather than just the transcript sequence
      FULLSEQ => 1,

      # maximum (genomic) length roughmodel for using fullseq, any
      # larger and it will switch to transcript seq to save CPU and
      # mem
      MAX_TRANSCRIPT => 1000000,

      # number of reads to align in each batch
      BATCH_SIZE => 10000,
    },
    bam2introns => {},
  }
';
  return $str;
}

sub GSNAP_header { 
  my $str = '# package Bio::EnsEMBL::Analysis::Config::GeneBuild::Gsnap
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Config::GeneBuild::Gsnap

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::GeneBuild::Gsnap

=head1 DESCRIPTION

This contains the specific configuraton for 
Bio::EnsEMBL::Analysis::RunnableDB::Gsnap 

=head1 CONTACT

=cut


package Bio::EnsEMBL::Analysis::Config::GeneBuild::Gsnap;

use strict;
use vars qw( %Config );

%Config = (
  GSNAP_CONFIG_BY_LOGIC =>  {
      DEFAULT =>  {

          # base path to the fastq
          INDIR => "/path/to/my/input",

          # path to the output directory
          OUTDIR => "/path/to/my/output",

          # Nmme given to the indexed genome when using gmap build
          GENOMENAME => "'.$RNASEQCONFIG->{GENOME_FILE}.'",
          # Directory containing the genome files
          GENOMEDIR  => "' .$RNASEQCONFIG->{GENOME_DIR}. '",
          # alignment options ( just for example )
          OPTIONS => "",

          # are the reads paired end? (1/0)
          PAIRED => "0",

          # path to the samtools binaries
          SAMTOOLS_PATH => "/software/solexa/bin/samtools",

          # optional header with additional information describing the sample
          HEADER => "",
    },';
  return $str;
}

sub SAM2BAM_header {
  my $str = '# package Bio::EnsEMBL::Analysis::Config::GeneBuild::Sam2Bam
#
# Cared for by EnsEMBL (dev@ensembl.org)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code


package Bio::EnsEMBL::Analysis::Config::GeneBuild::Sam2Bam;

use strict;
use vars qw( %Config );

%Config = (
  SAM2BAM_CONFIG_BY_LOGIC => {
    DEFAULT => {
      # directory containg the sam file(s)
      SAM_DIR => "'.$RNASEQCONFIG->{MERGE_DIR}.'/SAM",

      # path to the bam file to produce as output
      BAMFILE => "'.$RNASEQCONFIG->{MERGE_DIR}.'/introns.bam",

      # regex to identify which SAM files to merge
      REGEX => ".sam",

      # file containing all the readgroup headers used in the alignments (optional)
      HEADERFILE => "'.$RNASEQCONFIG->{OUTPUT_DIR}.'/all_headers.txt",

      # path to dumped genome file used for the alignment
      # it will make an index for it if one does not already exist
      GENOMEFILE =>
        "'.$RNASEQCONFIG->{GENOME_DIR}.'/'.$RNASEQCONFIG->{GENOME_FILE}.'",
    },
    sam2bam => {},
  }
';
  return $str;
}

sub REFINE_header {
  my $str = '# Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes;
#
# Cared for by EnsEMBL (dev@ensembl.org)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code


package Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes;

use strict;
use vars qw( %Config );

# Hash containing config info
%Config = (

         REFINESOLEXAGENES_CONFIG_BY_LOGIC => {
            DEFAULT => {
                # databases are defined as hash keys from Bio::EnsEMBL::Analysis::Config::Databases
                        OUTPUT_DB => "'.$RNASEQCONFIG->{REFINEDDB}.'",
            # intron db is not needed if you are using BAM files
            INTRON_DB => "",
            MODEL_DB  => "'.$RNASEQCONFIG->{ROUGHDB}.'",

            # Using bam file to fetch intron features overrides the INTRON_DB\n
';

  unless ( $use_gsnap ) {
$str .= '  	    INTRON_BAM_FILES => [
			{
		    	# location of the bam file
			FILE => "'.$RNASEQCONFIG->{MERGE_DIR}.'/introns.bam",
			# does the bam file(s) contain a mixture of spliced and unspliced reads such as you might get using Gsnap?
			MIXED_BAM => "0",
			# only take introns above this depth
			DEPTH => "0",
			# return only reads from the specified read groups
	               	GROUPNAME => [],
			},
		],
';
} else {
$str .= '  	    INTRON_BAM_FILES => [
			{
		    	# location of the bam file
			FILE => "'.$RNASEQCONFIG->{MERGE_DIR}.'/merged.bam",
			# does the bam file(s) contain a mixture of spliced and unspliced reads such as you might get using Gsnap?
			MIXED_BAM => "1",
			# only take introns above this depth
			DEPTH => "0",
			# return only reads from the specified read groups
	               	GROUPNAME => [],
			},
		  ], 
';
}

$str .= '            # write the intron features into the OUTPUT_DB along with the models
            WRITE_INTRONS => 1,

            # maximum number of times to loop when building all possible paths through the transcript
            MAX_RECURSIONS => 100000,

            # analysis logic_name for the dna_align_features to fetch from the INTRON_DB
            # If left blank all features will be fetched
            LOGICNAME => [],

            # logic name of the gene models to fetch
            MODEL_LN  => "",

            # penalty for removing a retined intron
            RETAINED_INTRON_PENALTY => 2,

            # minimum size for an intron
            MIN_INTRON_SIZE  => 30,
            MAX_INTRON_SIZE  => 200000,

            # biotype to give to single exon models if left blank single exons are ignored
            SINGLE_EXON_MODEL => "single_exon",

            # minimum single exon size (bp)
            MIN_SINGLE_EXON => 1000,

            # minimum percentage of single exon length that is coding
            SINGLE_EXON_CDS => 66,

            # Intron with most support determines the splice sites for an internal exon
            # lower scoring introns with different splice sites are rejected

            STRICT_INTERNAL_SPLICE_SITES => 1,

            # In some species alternate splice sites for end exons seem to be common
            STRICT_INTERNAL_END_EXON_SPLICE_SITES => 1,

            # biotypes to give gene models if left blank these models will not get written to the output database
            # best score - model with most supporting intron features
            BEST_SCORE => "best",
            # all other possible models
            OTHER_ISOFORMS => "",
                        # max number of other models to make - blank = all
            OTHER_NUM      => "10",

                        # max number of other models to process - blank = all
            MAX_NUM      => "1000",

            # biotype to label bad models ( otherwise they are not written )
            BAD_MODELS     => "",

            # do you want to trim UTR
            TRIM_UTR => 1,
            # config for trimming UTR
            MAX_3PRIME_EXONS => 2,
            MAX_3PRIME_LENGTH => 5000,
            MAX_5PRIME_EXONS => 3,
            MAX_5PRIME_LENGTH => 1000,
            # % of average intron score that a UTR intron must have
            REJECT_INTRON_CUTOFF => 5,
                   },
              refine_all => {},
';
  return $str;
}

sub BLAST_header {
  my $str = '# Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastRNASeqPep;
#
# Cared for by EnsEMBL (dev@ensembl.org)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code


package Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastRNASeqPep;

use strict;
use vars qw( %Config );

# Hash containing config info
%Config = (

         BLASTRNASEQPEP_CONFIG_BY_LOGIC =>
           {
            DEFAULT => {
                # databases are defined as hash keys from Bio::EnsEMBL::Analysis::Config::Databases
                        OUTPUT_DB => "'.$RNASEQCONFIG->{BLASTDB}.'",
            MODEL_DB  => "'.$RNASEQCONFIG->{REFINEDDB}.'",

            # If left blank all refined genes will be fetched
            LOGICNAME => "refine_all",

            # path to index to fetch the sequence of the blast hit to calculate % coverage
            INDEX => "'.$RNASEQCONFIG->{UNIPROTINDEX}.'",
                 },
            rnaseqblast => {},
         }
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

sub check_rule {
  my ($rule) = @_;
  # see if you can fetch it first
  my $check = $ra->fetch_by_goal($rule->goalAnalysis);
#  print "CHECK $check \n";
  unless ( $check ) {
    return 1;
  }
  return;
}
