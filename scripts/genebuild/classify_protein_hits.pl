# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::UniProtKB;

use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Utils::Exception qw ( throw warning ) ;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw ( get_input_arg ) ;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Analysis::Config::GeneBuild::ClassifyProteins;
$|=1;

my $inserted_new_analysis = 0;

#  SETTING THE DEFAULT VALUES FOR DB-CONNECTION ETC
my %opt = (
           dbname => '',
           dbhost => '',
           dbuser => 'ensro',
           dbpass => '',
           dbport => '3306' ,
           outfile => "protein_classification.sql",
           filename_acc_no_desc => "acc_no_desc_found.txt" ,
           acc_input_file => 'acc_input_file.txt' ,
           output_dir => '.' ,
          );

my $script_call ="perl " .$0." ".join(" ",@ARGV);

my $default_paf_logic_name = 'uniprot' ;
my $paf_logic_name = '';
&GetOptions(
             \%opt ,
            'dbhost|host|h=s',
            'dbname|db|D=s',
            'dbuser|user|u=s',
            'dbpass|pass|p=s',
            'dbport|port|P=i',
            'config=s',
            'delete' ,      # avoid interactive mode and overwrite files
            'test!',        # only fetch the first 5000 entries
            'output_dir=s', # dir where you want the output written to
            'outfile=s',    # name of sql output file
            'paf_logic_name|paf_logic_names=s',   # name of protein_align_feature logic_name to fetch ( usually Uniprot )
                                                  # or use comma-separeated list uniprot_1,uniprot_pe2,uniprot_pe_3
            'verbose!',
            'acc_input_file=s',
            'filename_acc_no_desc=s',
            'rerun!',
            'run_multiple_jobs!',
            'no_mfetch_sprot_dat=s',   # filepath to sprot dat file, use when mfetch is not available
            'no_mfetch_trembl_dat=s',  # filepath to trembl dat file, use when mfetch is not available
            'website!', # Use it only if you have few number of protein to check otherwise you will hammer UniProt
           ) ;

if ( ($opt{no_mfetch_sprot_dat} and !$opt{no_mfetch_trembl_dat})
     or
     (!$opt{no_mfetch_sprot_dat} and $opt{no_mfetch_trembl_dat})
   ) {
  throw("If you want to run on no_mfetch mode, you must provide two filepaths. One of the following arguments is missing: no_mfetch_sprot_dat, no_mfetch_trembl_dat\n");
}

if ($opt{no_mfetch_sprot_dat} and $opt{no_mfetch_trembl_dat}) {

  if ($opt{no_mfetch_sprot_dat} !~ /uniprot_sprot/) {
    throw("The no_mfetch_sprot_dat filename must contain the string 'uniprot_sprot'\n");
  }

  if (!(-e $opt{no_mfetch_sprot_dat})) {
    throw("File $opt{no_mfetch_sprot_dat} does not exist.\n");
  }
  if (!(-e $opt{no_mfetch_trembl_dat})) {
    throw("File $opt{no_mfetch_trembl_dat} does not exist.\n");
  }
}

if ( $opt{rerun} ) {

#   $opt{filename_acc_no_desc}.="_".$opt{acc_input_file};
#   $opt{outfile}.="_".$opt{acc_input_file};
} else {
  $opt{outfile} = $opt{output_dir}."/".$opt{outfile};
  $opt{acc_input_file} = $opt{output_dir}."/".$opt{acc_input_file};
  $opt{filename_acc_no_desc} = $opt{output_dir}."/".$opt{filename_acc_no_desc};
}

unless ( $opt{paf_logic_name} ) {
 print "\nYou did not specify any paf_logic_name .. I will use default logic_name : $default_paf_logic_name \n\n".
       "\tYou can specify a logic_name using option : -paf_logic_name uniprot\n" ;
 $paf_logic_name = $default_paf_logic_name ;
}else {
  print "using logic_name(s) : $opt{paf_logic_name}\n" ;
  $paf_logic_name = $opt{paf_logic_name} ;
}

print "\n\n ACC    - acc-numbers which could'nt be checked are written to file $opt{filename_acc_no_desc}\n" ;
print " MEMORY - consider using  bsub -q normal -M6000 -R\"select[mem>6000] rusage[mem=6000]\"\n" ;
print " FILE   - sql-statements will be written to $opt{outfile}\n" ;
print " QUICK  - use -run_multiple_jobs flag to split the run into 10 jobs which you can run in bsub\n\n" ;
print "REMEMBER : \nYou have to use the 64bit perl to run jobs which use over 3 gigs of memory\n";
print "SQL-output will be written to $opt{outfile} \n" ;

if (-e $opt{outfile} ) {
     if ( $opt{delete}  ) {
       system("rm $opt{outfile}") ;
     } else {
       print "\noutput file $opt{outfile} exists - please delete it first or rename it, or use -delete option  \n" ;
       exit(0);
     }
}



if ( -e $opt{filename_acc_no_desc}) {
      if ( $opt{delete}) {
        print "$opt{filename_acc_no_desc} will be over-written\n" ;
      }else {
        throw("file $opt{filename_acc_no_desc} exist - please use -filename_acc_no_desc <NEW_FILENAME> and supply new filename \n".
              " or use -delete option to delete file and avoid interactive mode");
      }
  }


check_phylum() ;
open SQL, (">$opt{outfile}") || die " can't write to file $opt{outfile}\n" ;

#unless ( $opt{config} ) {
#   print "you have to supply a config file to run this script .. \n" ;
#   exit(0);
#}
#
#if ( -e $opt{config} ) {
#   use lib "$opt{config}"  ;
#   print "config read\n" ;
#} else {
#   print "the config file : $opt{config} does not exist ...... EXITING ...\n" ;

#   exit(0);
#}



print "\n\nI will use the following logic_names out of the config to classify your proteins : \n" ;
for ( keys %$PROTEIN_CLASSIFICATION) {
  print " - $_\n" ;
}


if ( $opt{rerun} ) {
  print " reading data from file $opt{acc_input_file} \n" ;
} else {
  print "\n\nFetching data out of protein align feature table \n" ;
}
my $dba = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
                          -dbname => $opt{dbname},
                          -host => $opt{dbhost},
                          -user => $opt{dbuser},
                          -pass => $opt{dbpass},
                          -port => $opt{dbport},
                         );


# test if input_id_analysis table exists ( this is 'cause we create  input_id_type_analysis entries
my $sth = $dba->prepare("show tables like 'input_id_type_analysis' " ) ;
$sth->execute();
my $result = $sth->fetchall_arrayref();
 for my $row (@$result) {
   my ($table_name)  = @$row ;
   unless ( $table_name =~/input_id_type_analysis/ ) {
      throw( " make sure that your database contains the pipeline-tables as analysis will be uploaded \n");
   }
}


my %paf;
my %paf_hash;

# data fetch - quick but dirty ...

my %analysis_hash;

my $statement =
  "select a.analysis_id, a.logic_name " .
  "from analysis a, protein_align_feature pf " .
  "where pf.analysis_id = a.analysis_id group by logic_name ";
$sth = $dba->prepare($statement);
$sth->execute();
$result = $sth->fetchall_arrayref();

print
  "These analysis have entries in the protein_align_feature-table : \n";
for my $row (@$result) {
  my ( $analysis_id, $logic_name ) = @$row;
  print "analysis_id = $analysis_id, logic_name = $logic_name\n";
  $analysis_hash{$logic_name} = $analysis_id;
}


if ( $opt{rerun} ) {
     print "rerunning -reading data from $opt{acc_input_file} \n" ;
     open (PAF,$opt{acc_input_file}) || die " Can't read file $opt{acc_input_file}\n" ;
     for ( <PAF> ) {
       my ( $logic_name, $hit_name ) = split ;
       ${$paf_hash{$logic_name}}{$hit_name} = 1 ;
     }
    close(PAF);
} else {
    print "\nFetching data.....\n" ;

     my $paf_logic_name_string  ;
     my @pf = split /\,/,$paf_logic_name ;
     for ( @pf ) {
       $paf_logic_name_string .= "'$_',";
     }
    $paf_logic_name_string =~s/,$//;
    $statement = "select hit_name, logic_name from analysis a,  protein_align_feature pf ".
               " where pf.analysis_id = a.analysis_id and logic_name in ($paf_logic_name_string) ";

    if ( $opt{test} ) {
      $statement .= " limit 5000 " ;
      print " TEST RUN - limiting data to 5000 entries ! \n" ;
    }
    print "SQL $statement \n" ;
    $sth = $dba->prepare($statement) ;
    $sth->execute();

    my $result = $sth->fetchall_arrayref();
    for my $row (@$result) {
      my ($hit_name, $logic_name ) = @$row ;
      ${$paf_hash{$logic_name}}{$hit_name} = 1 ;
   }
   my $line_counter = 0 ;
   print "writing fetched data to file $opt{acc_input_file}\n" ;
   open (PAF,">$opt{acc_input_file}") || die " can't write to file $opt{acc_input_file}\n" ;
     for my $logic_name ( keys %paf_hash ) {
        my %acc = %{ $paf_hash{$logic_name} } ;
        for my $hn ( keys %acc ) {
          $line_counter++;
          print PAF "$logic_name\t$hn\n" ;
        }
     }
   close PAF;


   if ( $opt{run_multiple_jobs} ) {
     open(OUTBS,">cmds_to_bsub.sh");
     my $nr_of_jobs = 10 ;
     my $line_count = sprintf("%.0f", $line_counter / $nr_of_jobs ) ;

     my $cmd = " split -d -l $line_count $opt{acc_input_file} $opt{acc_input_file}.acc_" ;
     system($cmd) ;

     my $files = sprintf("%.0f" , ($line_counter / $line_count)+0.5 );
     print "files : $files written \n"  ;

     print "Here are the command lines to submit the multiple jobs :\n" ;

     $script_call=~s/-?-run_multiple_jobs//g;
     for ( my $nr=0;$nr<$files; $nr++ ) {
        my $file_name = sprintf("%02d",$nr);
        my $cmd =  "bsub -M 500 -R 'select[mem>500] rusage[mem=500]' -oo ".$opt{output_dir}."/acc_$nr.log -J cfp_$nr $script_call -rerun -acc_input_file $opt{acc_input_file}.acc_$file_name ";
        $cmd .= " -outfile ".$opt{outfile}.".acc_$nr";
        $cmd .= " -filename_acc_no_desc $opt{filename_acc_no_desc}_$nr\n" ;
        #$cmd =~s/-filename_acc_no_desc\s+(\w+)\s+/acc_no_desc_found_$nr\.txt/;
        print  $cmd ;
        print OUTBS $cmd;
     }

     close OUTBS;  
     print "\n\nPlease ensure that your bsubs (-dbuser and -dbpass options) allow write access to the database\n\n";
     exit(0);
     # split written file $opt{acc_input_file} into multiple files
     # which are in different directories
     # run perl classifyProteins.pl with different output_dir options
   }
}


print "fetching *done*\n" ;
print "Entries found in protein_align_feature table of $opt{dbname} @ $opt{dbhost} : \n" ;

my @hit_ids_investigate;
for my $lg( keys %paf_hash ) {
  print "Analysis : $lg   hits_found:  " . scalar(keys %{$paf_hash{$lg}}) . "\n";
  push @hit_ids_investigate, keys %{$paf_hash{$lg}};
}

print "\n\nHave " . scalar(@hit_ids_investigate) . " different ID's to fetch ...\n" ;


# fetch analysis for input_id types


  my @fields_to_fetch  = qw ( Taxon acc org pe crd des ) ;
  my $obj;
  if ($opt{website}) {
    print "looking on www.uniprot.org to fetch description etc ..\n" ;
    $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::UniProtKB->new();
    $obj->use_website(1);
  } elsif ($opt{no_mfetch_sprot_dat} or $opt{no_mfetch_trembl_dat}) {
    $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::UniProtKB->new();
    $obj->uniprot_file([$opt{no_mfetch_sprot_dat}, $opt{no_mfetch_trembl_dat}]);
  } else {
    print "running mfetch to fetch description etc ..\n" ;
    $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch->new();
  }
  $obj->verbose(1) if $opt{verbose};

  # mfetch -f "Taxon acc org pe crd des" -d uniprot -i "acc:F6M3M6"
     # AC   F6M3M6;
     # DT   27-JUL-2011, integrated into UniProtKB/TrEMBL.
     # DT   27-JUL-2011, sequence version 1.
     # DT   29-MAY-2013, entry version 10.
     # DE   SubName: Full=Hypoxia-inducible factor 2 alpha;
     # DE   Flags: Fragment;
     # OS   Lepisosteus oculatus (Spotted gar).
     # OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
     # OC   Actinopterygii; Neopterygii; Semionotiformes; Lepisosteidae;
     # OC   Lepisosteus.
     # PE   2: Evidence at transcript level;

  my ($all_entries, $not_found) =  @{$obj->get_Entry_Fields(\@hit_ids_investigate,\@fields_to_fetch) } ;
  
  # print now try to  fetch with wildcards     - second round
  if (@$not_found) {
    print "\nNO description found by mfetch for ".scalar(@$not_found)." entries !!!! - we try with wildcards now ...  \n" ;
    for ( @$not_found ) {
      print  "$_ not_found\n" ;
    }
  }

  my @not_found_2 ;
  for my $acc ( @$not_found ) {
     my $nacc = $acc ;
     $nacc=~s/\..*//g;
     $nacc=~s/\-.*//g;
     print " fetching $acc ---->  $nacc\n" if $opt{verbose} ;

     my $search_acc = " -d uniprot -i \"acc:$nacc\%\"";

     my ($refetch_entries, $refetch_not_found ) = $obj->get_Entry_Fields($search_acc, \@fields_to_fetch);

     if ( scalar(@$refetch_not_found) > 0 ) {
        push @not_found_2, $acc ;
     }

     print "adding re-fetched values to entries{$acc}\n" if $opt{verbose} ;
     foreach my $k ( keys %$refetch_entries ) {
         my %tmp = %{$refetch_entries->{$k}};
         for my $field ( keys %tmp  ) {
            $all_entries->{$acc}{$field}  = $tmp{$field} ;
         }
     }
  }


  open (A,">$opt{filename_acc_no_desc}") || die "Can't write to file \n" ;
  for ( @not_found_2) {
     print A "entry_really_not_found_by_mfetch_for  $_\n" ;
  }
  close(A);
  print " descriptions fetched ........\n" ;

#  now let's classify this
 my @logic_names = @{get_logic_names_for_classification()} ;


#
# check if the analysis for classification already exist first ... if not create it and upload it .
#

print " checking for analysis in db \n" ;


for my $lg  ( @{get_logic_names_for_classification()}) {
  my $ana =  get_analysis($dba, $lg);
  my $dbID ;
  unless ( $ana ) {
    warning( "Analysis : $lg mentioned in Config-file does not exist....\n" .
              " I will create an analysis object for you and upload it in your db\n");

     my $upload_ana = Bio::EnsEMBL::Pipeline::Analysis->new( -logic_name => $lg,
                                                   -input_id_type => "classifyProteins" ,
                                                   -module => "BlastGenscanPep",
                                                 ) ;
     $dba->get_AnalysisAdaptor->store($upload_ana) ;
     $inserted_new_analysis += 1;
     warning("Added analysis $lg");
     # now get analysis_id ....
     $ana =  get_analysis($dba, $lg);
     $dbID = $ana->dbID ;
     print "Analysis $lg with analysis_id $dbID created successfully\n" ;
  } else {
    print "analysis $lg already exist :dbID ---> " . $ana->dbID . " \n" ;
  }
  $analysis_hash{$ana->logic_name}=$ana->dbID;
}

for my $lg ( @logic_names ) {
  my $analysis_id = $analysis_hash{$lg} ;

  my %filtered_hits;
  print "\n------  logic name : $lg [$analysis_id] ------\n" ;

  #
  # get PE levels, fragment vs non-fragment and swissprot vs trembl
  #
  my @pe_levels = @{get_pe_levels_for_classification($lg)} ;
  my @uniprot_dbs = @{get_uniprot_dbs_for_classification($lg)} ;
  my @fragment = @{get_fragments_for_classification($lg)} ; 

  my %filters;
  $filters{'PE_LEVEL'} = get_pe_levels_for_classification($lg);
  $filters{'UNIPROT_DB'} = get_uniprot_dbs_for_classification($lg);
  $filters{'FRAGMENT'} = get_fragments_for_classification($lg);

  # Looks like the logic has to work from bottom up here, hence
  # filter by fragment first
  # filter by PE level second
  # filter by uniprot db third
  if ( scalar(@pe_levels) > 0 && scalar(@uniprot_dbs) > 0 && scalar(@fragment) > 0 ) {
    print "PE-Levels wanted: " . join (" " , @pe_levels ) . "\n" ;
    print "UniProt dbs wanted: " . join (" " , @uniprot_dbs ) . "\n" ;
    print "filtering " . scalar(keys %$all_entries) . " by Fragment, PE LEVEL and UNIPROT_DB ......\n";
    %filtered_hits = %{get_protein_hits_by_fragment_pe_level_and_uniprot_db($all_entries, \%filters )} ;
    print "PE-filter and UNIPROT_DB-filter result: " . scalar(keys(%filtered_hits)) . " left\n" ;
    
  } elsif ( scalar(@pe_levels) > 0 && scalar(@uniprot_dbs) > 0 ) {
    print "PE-Levels wanted: " . join (" " , @pe_levels ) . "\n" ;
    print "UniProt dbs wanted: " . join (" " , @uniprot_dbs ) . "\n" ;
    print "filtering " . scalar(keys %$all_entries) . " by PE LEVEL and UNIPROT_DB ......\n";
    %filtered_hits = %{get_protein_hits_by_pe_level_and_uniprot_db($all_entries, \%filters )} ;
    print "PE-filter and UNIPROT_DB-filter result: " . scalar(keys(%filtered_hits)) . " left\n" ;

  } elsif ( scalar(@pe_levels) > 0 ) {
    print "PE-Levels wanted: " . join (" " , @pe_levels ) . "\n" ;
    print "filtering " . scalar(keys %$all_entries) . " by PE LEVEL ......\n";
    %filtered_hits = %{get_protein_hits_by_pe_level_and_uniprot_db($all_entries, \%filters )} ;
    print "PE-filter result: " . scalar(keys(%filtered_hits)) . " left\n" ;
    print "skipping .......: UNIPROT_DB-filter\n" ;

  } elsif ( scalar(@uniprot_dbs) > 0 ) {
    print "UniProt dbs wanted: " . join (" " , @uniprot_dbs ) . "\n" ;
    print "filtering " . scalar(keys %$all_entries) . " by UNIPROT_DB ......\n";
    %filtered_hits = %{get_protein_hits_by_pe_level_and_uniprot_db($all_entries, \%filters )} ;
    print "UNIPROT_DB-filter result: " . scalar(keys(%filtered_hits)) . " left\n" ;
    print "skipping .......: PE-filter\n" ;

  } else {
    print "skipping .......: UNIPROT-filter and PE-filter\n" ;
    %filtered_hits = %$all_entries ;
  }


  #
  # get proteins for phylum
  #

  my @wanted_phylum = @{get_phylum_for_classification($lg)} ;
  my @excluded_phylum = @{get_phylum_exclusions_for_classification($lg)};

  if ( scalar(@wanted_phylum) > 0 || scalar(@excluded_phylum) > 0 )  {

    print "phylum wanted   : " . join ( " " , @wanted_phylum ) . "\n" ;
    print "phylum excluded : " . join ( " " , @excluded_phylum) . "\n" ;
    print "filtering " . scalar(keys %filtered_hits) . " by phylum ......\n";
    %filtered_hits = %{get_and_exclude_protein_hits_by_phylum(\%filtered_hits, \@wanted_phylum, \@excluded_phylum)} ;
    print "phylum-filter   : " . scalar(keys(%filtered_hits)) . " left\n" ;
  } else {
    print "skipping .......: phylum_filter\n" ;
  }

  #
  # get organisms
  #
  print "try to get organism : $lg \n" ;

  my @wanted_org = @{get_organism_for_classification($lg)} ;
  print "done\n" ;
  my @excluded_org = @{get_organism_excluded_classification($lg)} ;

  if ( scalar(@wanted_org) > 0 || scalar(@excluded_org) > 0 )  {
    print "org wanted      : " . join ( " " , @wanted_org ) . "\n" ;
    print "org excluded    : " . join ( " " , @excluded_org) . "\n" ;
    print "filtering " . scalar(keys %filtered_hits) . " by organism ......\n";
    %filtered_hits = %{get_and_exclude_protein_hits_by_organism(\%filtered_hits, \@wanted_org, \@excluded_org)} ;
    print "organism_filter : " . scalar(keys(%filtered_hits)) . " left\n" ;
  } else {
    print "skipping .......: org-filter\n" ;
  }

  # create a new hash with the keys swapped so that
  # it is easier to fetch the analysis for a given hit_name
  my %paf_hit_log;
  foreach my $logic_name (keys %paf_hash) {
    foreach my $hit_name (keys %{$paf_hash{$logic_name}}) {
      my ($hit_name_no_version,$version) = split(/\./,$hit_name);
      $paf_hit_log{$hit_name}{$logic_name} = 1;
      $paf_hit_log{$hit_name_no_version}{$logic_name} = 1; # add hit name without version for no mfetch option
    }
  }

  foreach my $hit_name ( keys %filtered_hits ) {
    # if we have specified multiple logic names,
    # we'll get each hit_name multiple times (if it was present in each logic name)
    # and we want to create a SQL update command for each logic name
    my $orig_analysis_id;
    foreach my $logic_name (keys %{$paf_hit_log{$hit_name}}) {
      if ($paf_hit_log{$hit_name}{$logic_name} == 1) {
        $paf_hit_log{$hit_name}{$logic_name} = 0; # keep track of the already created SQL updates
        $orig_analysis_id = $analysis_hash{$logic_name};
        last;
      }
    }
    if (!$orig_analysis_id) {
      throw("Couldn't find original analysis id for $hit_name\n");
    }

    if ($opt{no_mfetch_sprot_dat} and $opt{no_mfetch_trembl_dat}) {
      printf( SQL "update protein_align_feature " .
        "set analysis_id = %d " .
        "where hit_name like '%s%' and analysis_id = %d;\n",
        $analysis_id, $hit_name, $orig_analysis_id );
    }
    else {
      printf( SQL "update protein_align_feature " .
        "set analysis_id = %d " .
        "where hit_name = '%s' and analysis_id = %d;\n",
        $analysis_id, $hit_name, $orig_analysis_id );
    }
  }
}    # end logic names
close (SQL);


if ($inserted_new_analysis >= 1) {
  warning("Script ended successfully but please be aware that one or more new analyses ".
          "have been added to this database, ".$opt{dbname}." at ".$opt{dbhost}.
          ". Please synchronise the analysis tables across all databases where necessary.");
}






sub get_and_exclude_protein_hits_by_field {
  my ($entries,$field,$expr_to_match,$expr_to_exclude) = @_ ;
  my %result ;

  unless ( $expr_to_exclude ) {
    $expr_to_exclude = [];
  }
  my $expression_to_match = 0 ;
  if ( scalar(@$expr_to_match) > 0 ) {
    $expression_to_match = 1
  }

  ACC: for my $acc ( keys %$entries ) {
     my $attr_to_test = ${$entries}{$acc}{$field} ;

     # this is the routine if you like to match AND exclude
     #print "processing $acc\n" ;

     if ( $expression_to_match ) {
       for my $string_match ( @$expr_to_match ) {
#        print "testing : $string_match  - $attr_to_test " ;
         if ($attr_to_test=~m/$string_match/ ) {
           print "MATCH  $attr_to_test $string_match $acc \n" if $opt{verbose} ;
            # we have a match
            if ( scalar(@$expr_to_exclude) > 0 ) {
              my $exc = 0;
              for my $exclude ( @$expr_to_exclude ) {
                if ( $attr_to_test =~m/$exclude/ ) {
                  $exc = 1;
                  print "excluding : $exclude - $attr_to_test \n" if $opt{verbose} ;
                  # we need to exclude this as it's listed to be excluded .......
                }# else {
                #   $result{$acc}= ${$entries}{$acc};
                #}
              }
              $result{$acc}= ${$entries}{$acc} unless $exc;
            } else {
              $result{$acc}= ${$entries}{$acc};
            }
            next ACC;
         } else {
           if ( $opt{verbose} ) {
             print " sorry - no match. the ACC $acc does not match the ".
             " $field-field \"$string_match\" criteria [ $attr_to_test =~ m/$string_match/ ] :-(\n" ;
           }

         }
       }
     } else {
        #print "we only have to exclude stuff, not to match stuff - for field $field .... \n" ;
         my $exc = 0;
         for my $exclude ( @$expr_to_exclude ) {
            if ( $attr_to_test =~m/$exclude/ ) {
                $exc = 1;
                #print "exclude matches : $exclude $attr_to_test \n" ;
                # we need to exclude this as it's listed to be excluded .......
             } #else {
                  #print "no match ----storing \n" ;
                #  $result{$acc}= ${$entries}{$acc};
             #}
        }
        $result{$acc}= ${$entries}{$acc} unless $exc;
     }
  }   # end of big loop .....
  return \%result;
}
 

sub get_protein_hits_by_fragment_pe_level_and_uniprot_db {
  my ($entries,$filter) = @_ ;
  my $de_result = {};
  my $frag_string = ['Fragment'];
  if (@{$filter->{'FRAGMENT'}}[0]) {
    # I want a fragment
    $de_result = get_and_exclude_protein_hits_by_field($entries,"DE",$frag_string) ;          
  } else { 
    # I don't want a fragment
    $de_result = get_and_exclude_protein_hits_by_field($entries,"DE", [], $frag_string) ;          
  }
  my $pe_result = get_and_exclude_protein_hits_by_field($de_result,"PE",$filter->{'PE_LEVEL'}) ;
  return  get_and_exclude_protein_hits_by_field($pe_result,"DT",$filter->{'UNIPROT_DB'}) ;
}

sub get_protein_hits_by_pe_level_and_uniprot_db {
  my ($entries,$filter) = @_ ;
  my $pe_result = get_and_exclude_protein_hits_by_field($entries,"PE",$filter->{'PE_LEVEL'}) ;
  return  get_and_exclude_protein_hits_by_field($pe_result,"DT",$filter->{'UNIPROT_DB'}) ;
}

sub get_protein_hits_by_pe_level {
  my ($entries,$pe_levels) = @_ ;
  return  get_and_exclude_protein_hits_by_field($entries,"PE",$pe_levels) ;
}

sub get_and_exclude_protein_hits_by_phylum {
  my ($entries ,$wanted_phylum,$excluded_phylum) = @_ ;
  return  get_and_exclude_protein_hits_by_field($entries,"OC",$wanted_phylum,$excluded_phylum) ;
}

sub get_and_exclude_protein_hits_by_organism{
  my ($entries ,$wanted_org, $excluded_org ) = @_ ;
  return  get_and_exclude_protein_hits_by_field($entries,"OS",$wanted_org,$excluded_org) ;
}

sub get_and_exclude_protein_hits_by_database{
  my ($entries ,$uniprot_dbs, ) = @_ ;
  return  get_and_exclude_protein_hits_by_field($entries,"DT",$uniprot_dbs) ;
}

sub get_logic_names_for_classification {
   return [keys %$PROTEIN_CLASSIFICATION];
}

sub get_pe_levels_for_classification {
  my ($logic_name) = @_ ;
  return ${$PROTEIN_CLASSIFICATION}{$logic_name}{"PE_LEVELS"} ;
}
sub get_phylum_for_classification {
  my ($logic_name) = @_ ;
  return ${$PROTEIN_CLASSIFICATION}{$logic_name}{"PHYLUM"};
}

sub get_phylum_exclusions_for_classification {
  my ($logic_name) = @_ ;
  return ${$PROTEIN_CLASSIFICATION}{$logic_name}{"EXCLUDED_PHYLUM"};
}

sub get_organism_for_classification {
  my ($logic_name) = @_ ;
  return ${$PROTEIN_CLASSIFICATION}{$logic_name}{"ORGANISM"};
}

sub get_organism_excluded_classification {
  my ($logic_name) = @_ ;
  return ${$PROTEIN_CLASSIFICATION}{$logic_name}{"EXCLUDED_ORGANISM"};
}

sub get_uniprot_dbs_for_classification {
  my ($logic_name) = @_ ;
  return ${$PROTEIN_CLASSIFICATION}{$logic_name}{"UNIPROT_DBS"} ;
}

sub get_fragments_for_classification {
  my ($logic_name) = @_ ; 
  return ${$PROTEIN_CLASSIFICATION}{$logic_name}{"FRAGMENT"} ;  
}                                                                                


sub get_analysis {
  my ( $dba, $logic_name ) = @_;
  my $ana =
    $dba->get_AnalysisAdaptor()->fetch_by_logic_name($logic_name);
  if ($ana) {
    print "fetched analysis $logic_name --> " . $ana->dbID . "\n";
    return $ana;
  }
  warning(
     "Analysis with logic_name : $logic_name could not be found in : " .
       $dba->dbname . " \@ " . $dba->host . "\n" );
  return undef;
}

sub check_phylum {
  my %phylum = ( 'Acanthomorpha'     => 1,
                 'Acanthopterygii'   => 1,
                 'Aconoidasida'      => 1,
                 'Actiniaria'        => 1,
                 'Actinopterygii'    => 1,
                 'Adrianichthyidae'  => 1,
                 'Aedes'             => 1,
                 'Agaricales'        => 1,
                 'Agaricomycetes'    => 1,
                 'Agaricomycetidae'  => 1,
                 'Agaricomycotina'   => 1,
                 'Alveolata'         => 1,
                 'Amniota'           => 1,
                 'Amphibia'          => 1,
                 'Anopheles'         => 1,
                 'Anophelinae'       => 1,
                 'Anthozoa'          => 1,
                 'Anura'             => 1,
                 'Apicomplexa'       => 1,
                 'Archosauria'       => 1,
                 'Arthropoda'        => 1,
                 'Ascomycota'        => 1,
                 'Atherinomorpha'    => 1,
                 'Aves'              => 1,
                 'Basidiomycota'     => 1,
                 'Batoidea'          => 1,
                 'Batrachia'         => 1,
                 'Beloniformes'      => 1,
                 'Bos'               => 1,
                 'Bovidae'           => 1,
                 'Bovinae'           => 1,
                 'Brachycera'        => 1,
                 'Callitrichinae'    => 1,
                 'Candida'           => 1,
                 'Carnivora'         => 1,
                 'Catarrhini'        => 1,
                 'Cebidae'           => 1,
                 'Cetartiodactyla'   => 1,
                 'Choanoflagellida'  => 1,
                 'Chondrichthyes'    => 1,
                 'Chordata'          => 1,
                 'Cnidaria'          => 1,
                 'Codonosigidae'     => 1,
                 'Coelurosauria'     => 1,
                 'Coprinopsis'       => 1,
                 'Craniata'          => 1,
                 'Crocodylidae'      => 1,
                 'Culex'             => 1,
                 'Culicidae'         => 1,
                 'Culicinae'         => 1,
                 'Culicini'          => 1,
                 'Culicoidea'        => 1,
                 'Cyprinidae'        => 1,
                 'Cypriniformes'     => 1,
                 'Danio'             => 1,
                 'Dikarya'           => 1,
                 'Dinosauria'        => 1,
                 'Diptera'           => 1,
                 'Discopyge'         => 1,
                 'Drosophila'        => 1,
                 'Drosophilidae'     => 1,
                 'Edwardsiidae'      => 1,
                 'Elasmobranchii'    => 1,
                 'Embryophyta'       => 1,
                 'Endopterygota'     => 1,
                 'Ephydroidea'       => 1,
                 'Euarchontoglires'  => 1,
                 'Eukaryota'         => 1,
                 'Euteleostei'       => 1,
                 'Euteleostomi'      => 1,
                 'Eutheria'          => 1,
                 'Fungi'             => 1,
                 'Galliformes'       => 1,
                 'Gallus'            => 1,
                 'Glires'            => 1,
                 'Haplorrhini'       => 1,
                 'Hexacorallia'      => 1,
                 'Hexapoda'          => 1,
                 'Hominidae'         => 1,
                 'Homo'              => 1,
                 'Hypnosqualea'      => 1,
                 'Insecta'           => 1,
                 'Laurasiatheria'    => 1,
                 'Lepidosauria'      => 1,
                 'Magnoliophyta'     => 1,
                 'Mammalia'          => 1,
                 'Mesobatrachia'     => 1,
                 'Metazoa'           => 1,
                 'Monosiga'          => 1,
                 'Muridae'           => 1,
                 'Murinae'           => 1,
                 'Muroidea'          => 1,
                 'Mus'               => 1,
                 'Muscomorpha'       => 1,
                 'Narcinidae'        => 1,
                 'Narcinoidei'       => 1,
                 'Nematocera'        => 1,
                 'Nematostella'      => 1,
                 'Neognathae'        => 1,
                 'Neoptera'          => 1,
                 'Neopterygii'       => 1,
                 'Neoteleostei'      => 1,
                 'Oryzias'           => 1,
                 'Oryziinae'         => 1,
                 'Ostariophysi'      => 1,
                 'Pecora'            => 1,
                 'Percomorpha'       => 1,
                 'Phasianidae'       => 1,
                 'Phasianinae'       => 1,
                 'Pipidae'           => 1,
                 'Pipoidea'          => 1,
                 'Piroplasmida'      => 1,
                 'Platyrrhini'       => 1,
                 'Pongo'             => 1,
                 'Primates'          => 1,
                 'Pterygota'         => 1,
                 'Rodentia'          => 1,
                 'Saccharomycetales' => 1,
                 'Saccharomycetes'   => 1,
                 'Saccharomycotina'  => 1,
                 'Saguinus'          => 1,
                 'Saurischia'        => 1,
                 'Sauropsida'        => 1,
                 'Sciurognathi'      => 1,
                 'Silurana'          => 1,
                 'Sophophora'        => 1,
                 'Spermatophyta'     => 1,
                 'Squalea'           => 1,
                 'Stegomyia'         => 1,
                 'Streptophyta'      => 1,
                 'Teleostei'         => 1,
                 'Testudines'        => 1,
                 'Tetradontoidea'    => 1,
                 'Tetraodon'         => 1,
                 'Tetraodontidae'    => 1,
                 'Tetraodontiformes' => 1,
                 'Theileria'         => 1,
                 'Theileriidae'      => 1,
                 'Theropoda'         => 1,
                 'Torpediniformes'   => 1,
                 'Tracheophyta'      => 1,
                 'Vertebrata'        => 1,
                 'Viridiplantae'     => 1,
                 'Vitaceae'          => 1,
                 'Vitales'           => 1,
                 'Vitis'             => 1,
                 'Xenopodinae'       => 1,
                 'Xenopus'           => 1,
                 'Xiphophorus'       => 1 );

  for my $k ( keys %$PROTEIN_CLASSIFICATION ) {
    my %conf = %{ ${$PROTEIN_CLASSIFICATION}{$k} };

    for my $p ( @{ $conf{PHYLUM} }, @{ $conf{EXCLUDED_PHYLUM} } ) {
      # print "k=$k, p=$p\n";
      if ( length($p) > 0 && !exists( $phylum{$p} ) ) {
        throw( "The phylum : $p for [$k] " .
              "can't be found in the phylum dictionary.\n" .
              "Make sure you've spelled it correctly " .
              "- if the dictionary is incomplete, please add it...\n" );
      }
    }

  }
} ## end sub check_phylum
