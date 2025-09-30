#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

=head1 NAME

  clone_analysis.pl

=head1 DESCRIPTION

  Given an input file this script will clone a specified analysis in
  the analysis table to a set of new logic names. The input id type
  will be cloned and by default the rules for the original analysis 
  will also be copied. In addition it will clone the analysis config 
  to the correpsonding logic names in any config files specified. It
  is possible to request that no rules are to be copied and also to
  request that the copy only occurs in the specified config files and
  not in the db.

=head1 OPTIONS

  -dbhost|host|h     Host of database to clone analyses in
  -dbport|port|P     Port of db
  -dbname|db|D       Name of db
  -dbuser|user|u     User to connect as (ensadmin if cloning in db tables)
  -dbpass|pass|p     Password required for cloning in db tables
  -file              The path to the input file
  -just_config       Use this to only clone the configs
  -no_rules          Do not clone any rules, overridden by -not_pipeline_db flag
  -not_pipeline_db   Only clone entry in analysis table/config files (no rules or
                     input_id_type entered to db)

=head1 INPUT FILE STRUCTURE

  The input file should consist of one or two lines. The first line
  should start with the logic name of the analysis you want to clone.
  This should be followed by the new logic names you want to clone to,
  with each new logic name separated by a space.

  The second line is optional and is for cloning the entry for the analysis
  in any config files listed. This line should start with 'CONFIG' and
  then a space separating each config. The path to the config should not
  be supplied, the script will search out the configs in whatever paths to
  'ensembl-config' you have in your PERL5LIB.

  Below are examples of input files:

  my_logic_name_to_clone my_new_logic_name_1 my_new_logic_name_2 my_new_logic_name_3
  CONFIG Blast.pm BatchQueue.pm

  Or if no config files just:

  my_logic_name_to_clone my_new_logic_name_1 my_new_logic_name_2 my_new_logic_name_3

  Real world example, clones the info for analysis rnaseqblast to spleen_blast,
  liver_blast, lung_blast in relevant db tables and Blast.pm, BlastRNASeqPep.pm
  and BatchQueue.pm:

  rnaseqblast spleen_blast liver_blast lung_blast
  CONFIG Blast.pm BlastRNASeqPep.pm BatchQueue.pm


=head1 EXAMPLES

  perl clone_analysis.pl -h my_host -u ensadmin -p **** \
    -D my_db -file my_input_file.txt

  The above would clone the analysis listed in the input file in the
  analysis, input_id_type_analysis, rule_goal and rule_conditions
  tables along with any configs listed.

  perl clone_analysis.pl -h my_host -u ensadmin -p **** \
  -D my_db -file my_input_file.txt -not_pipeline_db

  The above would only clone the entry in the analysis table and
  any configs listed in the input file. No rules or input id types
  would be inserted

  perl clone_analysis.pl -h my_host -u ensadmin -p **** \
  -D my_db -file my_input_file.txt -no_rules
  
  The above would clone any analyses listed in the input file in
  the analysis and input_id_type_analysis tables along with any 
  configs listed. Not rules would be added to the rule tables.
  
  perl clone_analysis.pl -h my_host -u ensadmin -p **** \
  -D my_db -file my_input_file.txt -just_config
  
  The above would not enter any data into the database, but would
  clone the analysis in any config files listed in the input file.


=head1 NOTES

  The script will (should) backup your config files in their respective
  directories. However, it's it's always a good idea to backup/commit
  in these situations. Nothing will be backed up in terms of db tables
  so do your own backups for that (though it is relatively simple to
  delete things from the relevant tables if something goes wrong).

  There is/was an issue with cloning entries in old copies of BatchQueue.pm, 
  this is not due to this script but to the rnaseq pipeline setup script
  malforming entries in BatchQueue.pm (the resource line was floating
  outside of the normal key-value pair structure). This has been fixed
  in that script, so shouldn't be a probably for recently created ones.

=cut


use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::Pipeline::AnalysisCreation;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::ConfigWriter;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use feature 'say';

sub usage {
  exec( 'perldoc', $0 );
  exit;
}

my $dbhost;
my $dbuser='ensadmin';
my $dbpass;
my $dbport= 3306;
my $dbname;
my $file;
my $help;
my $not_pipeline = 0;
my %to_clone = ();
my %analysis_logic_names = ();
my @config_list = ();
my %confirmed_config_paths;
my $no_rules = 0;
my $just_config = 0;
if ( !GetOptions( 
                  'host|dbhost|h:s'     => \$dbhost,
                  'dbname|db|D:s'       => \$dbname,
                  'user|dbuser|u:s'     => \$dbuser,
                  'pass|dbpass|p:s'     => \$dbpass,
                  'port|dbport|P:s'     => \$dbport,
                  'not_pipeline_db!'    => \$not_pipeline,
                  'no_rules!'           => \$no_rules,
                  'just_config!'        => \$just_config,
                  'file=s'              => \$file,
                  'help!'               => \$help, 
                ) || $help )
{
  usage();
}

if(!($dbhost) || !($dbuser) || !($dbname)) 
{
  print STDERR
    "need to pass in database arguments for script to work\n";
  print STDERR
    "-dbhost $dbhost -dbuser $dbuser -dbpass $dbpass -dbname" .
    " $dbname -dbport $dbport\n";
  usage();
}



if(!$file) 
{
  print STDERR
    "You need to pass a file name which either represents a " .
    "analysis config file to be read and stored in the database or " .
    "written based on the analysis objects in the database\n";
  usage();
}

my $db;

unless($not_pipeline) 
{
  $db =
    new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor( -host   => $dbhost,
                                                  -user   => $dbuser,
                                                  -pass   => $dbpass,
                                                  -dbname => $dbname,
                                                  -port   => $dbport, );
}

else 
{
  $db =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dbhost,
                                        -user   => $dbuser,
                                        -pass   => $dbpass,
                                        -dbname => $dbname,
                                        -port   => $dbport, );
}

# Load input file
open(IN,"<$file");
my @to_clone = <IN>;
close IN;

# Parse input file
process_input();

# Read analyses from db
my $analyses = read_db($db);
$analyses = [ sort { $a->dbID() <=> $b->dbID() } @{$analyses} ];

# Pull out all logic names
get_analysis_logic_names();


# Loop through all logic names
# If a logic name matches one on the clone list then loop through all the
# new logic names to create and clone them if they don't already exist.
# If there is a conflict, skip. If the just_config flag is in use skip.
foreach my $analysis (@{$analyses})
{
    
  if($to_clone{$analysis->logic_name()})
  {

    for(my $i=0; $i< scalar @{$to_clone{$analysis->logic_name()}}; $i++)
    {    
      my $no_conflict = check_name_conflict(${$to_clone{$analysis->logic_name()}}[$i]); 
      unless($no_conflict && !$just_config)
      {
        next;
      }

      say "Cloning analysis: ".$analysis->logic_name." to ".${$to_clone{$analysis->logic_name()}}[$i]; 
      clone_analysis($analysis,${$to_clone{$analysis->logic_name()}}[$i]);              

    } 
     
  }

}

# If a config list was provided clone the analysis in those files
if(@config_list)
{
  clone_analysis_config();
}

# Else if these is no conifg list but the just_config flag is in use warn the user that this makes no sense.
elsif($just_config)
{
  say "Warning: you requested to clone just the config files but no config list was parsed from the input file.";
  say "Nothing to clone."
}


=head2 process_input

  Function  : Reads the loaded input file from the to_clone array and places an
              array ref of the new logic names into toe to_clone hash with the
              logic name of the analysis to clone as a key
  Returntype: 
  Exceptions:
  Example   : 

=cut

sub process_input
{
  for(my $i=0; $i<scalar @to_clone; $i++)
  {
    my $current_line = $to_clone[$i];
    chomp $current_line;

    my @temp = split(' ',$current_line);
    my $clone_key = shift @temp;
    $to_clone{$clone_key} = \@temp;
    if($to_clone[$i+1] && $to_clone[$i+1] =~ s/^CONFIG //)
    {
      chomp $to_clone[$i+1];
      say "Config line: ".$to_clone[$i+1]."\n";
      @config_list = split(' ',$to_clone[$i+1]);            
      $i++; 
    }

  }

}


=head2 get_analysis_logic_names

  Function  : Reads the analyses array and pulls out the logic names
              from each analysis. Assigns each logic name as a key in
              the analysis_logic_names hash
  Returntype: 
  Exceptions:
  Example   : 

=cut

sub get_analysis_logic_names
{
  foreach my $analysis (@{$analyses})
  {
    $analysis_logic_names{$analysis->logic_name()} = 1;
  }

}


=head2 check_name_conflict

  Arg [1]   : Scalar logic name to check
  Function  : Takes in a logic name scalar and checks if it's already in the
              analysis_logic_names hash. If it is then this logic name already
              exists in the analysis table and thus is a conflict.
  Returntype: Scalar, false if there is conflict, true if there isn't conflict
  Exceptions:
  Example   : check_name_conflict("rnaseqblast")

=cut

sub check_name_conflict
{
  my ($new_logic_name) = @_;

  if($analysis_logic_names{$new_logic_name} && !$just_config)
  {
    say "Warning: logic name ".$new_logic_name." already exists in the analysis table, will skip";
    return 0;
  }

  else
  {
    return 1;
  }

}


=head2 clone_analysis

  Arg [1]   : Scalar, logic name to clone from 
  Arg [2]   : Scalar, new logic name to clone to
  Function  : Takes in a logic name to clone and the new logic name and then attemps to
              make a clone in various tables. If the db is a pipeline db it will attempt
              by default to copy the entry in the analysis table, the input_id_type and
              the rules. If it is not a pipeline db then it will only copy the entry in
              the analysis table. If the no_rules flag is set then it will not copy the
              rules even if the db is a pipeline db (it will copy the input_id_type).

  Returntype: 
  Exceptions: Throws if it is a pipeline db and an input id type can't be found.
  Example   : clone_analysis("rnaseqblast","spleen_blast")

=cut

sub clone_analysis
{
    my ($analysis_to_clone,$new_logic_name) = @_; 
    my $new_analysis;
    say "Adding analysis: ".$new_logic_name;  

    # If it's a pipeline db make a Pipeline::Analysis object, with an input id type
    unless($not_pipeline)
    { 
       $new_analysis =  new Bio::EnsEMBL::Pipeline::Analysis(
                                                              -logic_name      => $new_logic_name,
                                                              -db              => $analysis_to_clone->db(),
                                                              -db_version      => $analysis_to_clone->db_version(),
                                                              -db_file         => $analysis_to_clone->db_file(),
                                                              -program         => $analysis_to_clone->program(),
                                                              -program_version => $analysis_to_clone->program_version(),
                                                              -program_file    => $analysis_to_clone->program_file(),
                                                              -gff_source      => $analysis_to_clone->gff_source(),
                                                              -gff_feature     => $analysis_to_clone->gff_feature(),
                                                              -module          => $analysis_to_clone->module(),
                                                              -module_version  => $analysis_to_clone->module_version(),
                                                              -parameters      => $analysis_to_clone->parameters(),
                                                              -created         => $analysis_to_clone->created(),
                                                              -input_id_type   => $analysis_to_clone->input_id_type(),
                                                            );

      my $input_id_type = $new_analysis->input_id_type();
 
      if(!$input_id_type)
      {
        throw("Could not find a corresponding input id type! If you do not need one then run this script with the -pipeline flag set to 0");
      }

       say "Added input id type: ".$input_id_type;

    }

    # If it's not a pipeline db make a normal analysis object, without an input id type
    else
    {
      $new_analysis =  new Bio::EnsEMBL::Analysis(
                                                              -logic_name      => $new_logic_name,
                                                              -db              => $analysis_to_clone->db(),
                                                              -db_version      => $analysis_to_clone->db_version(),
                                                              -db_file         => $analysis_to_clone->db_file(),
                                                              -program         => $analysis_to_clone->program(),
                                                              -program_version => $analysis_to_clone->program_version(),
                                                              -program_file    => $analysis_to_clone->program_file(),
                                                              -gff_source      => $analysis_to_clone->gff_source(),
                                                              -gff_feature     => $analysis_to_clone->gff_feature(),
                                                              -module          => $analysis_to_clone->module(),
                                                              -module_version  => $analysis_to_clone->module_version(),
                                                              -parameters      => $analysis_to_clone->parameters(),
                                                              -created         => $analysis_to_clone->created(),
                                                            );
    }

    $db->get_AnalysisAdaptor->store($new_analysis);
    say "Analysis inserted in db\n";
  
    # As long as it's a pipeline db and the -no_rules flag isn't set then copy and store the rules
    unless($not_pipeline || $no_rules)
    {
      my $rule_adaptor = $db->get_RuleAdaptor;
      my $old_rule = $rule_adaptor->fetch_by_goal($analysis_to_clone);
      my $new_rule = Bio::EnsEMBL::Pipeline::Rule->new(
                                                        -goalanalysis => $new_analysis,
                                                        -adaptor => $db,
                                                      );

      my @old_conditions = @{$old_rule->list_conditions};
      foreach(@old_conditions)
      { 
        $new_rule->add_condition($_);
        say "Added condition: ".$_; 
      } 

      $db->get_RuleAdaptor->store($new_rule);
      say "Rule(s) inserted in db\n" ;
    }

    elsif($not_pipeline)
    {
      say "No input id type or rules were added as the -not_pipeline_db flag was set";
    }

    elsif($no_rules)
    {
      say "No rules were added to the db as the -no_rules flag is in use"; 
    }

}


=head2 clone_analysis_config

  Function  : Loops through the config list array and attempts to find any config
              files listed in it in the ensembl-config dir on the PERL5LIB paths.
              Once it finds a config it will use ConfigWriter to backup the config
              to the dir it's in and write out the cloned analyses. If the module 
              is BatchQueue.pm a different ConfigWriter subroutine is used as the
              structure of BatchQueue is different to a regular analysis config.

  Returntype: 
  Exceptions: Throws if it can't find the module in the PERL5LIB
  Example   : 

=cut

sub clone_analysis_config
{
  foreach my $module_name (@config_list)
  {

    my $path = find_config_path($module_name);
    unless($path)
    {
      throw("Couldn't find module ".$module_name." in your PERL5LIB. Check the PERL5LIB and module name.".
            " Make sure to have the .pm on the end of the module name in the input file.");
    }
    my $config_writer = new Bio::EnsEMBL::Analysis::Tools::ConfigWriter(
                                                                                     -MODULENAME => $path,
                                                                       );
    say "Copying entries in module: ".$module_name;
    say "Found module ".$module_name." on path:";
    say $path;
    
    foreach my $logic_name_to_copy (keys(%to_clone))
    { 
     
      foreach my $new_logic_name (@{$to_clone{$logic_name_to_copy}})
      { 
        
        if($module_name eq "BatchQueue.pm")
        {          
          $config_writer->copy_analysis_from_batchqueue($logic_name_to_copy,$new_logic_name);
        }

        else
        {
          $config_writer->copy_analysis_from_config($logic_name_to_copy,$new_logic_name);
        }

        say "Copied config info for ".$logic_name_to_copy." to ".$new_logic_name;
  
      }
  
    }
 
    $config_writer->write_config(1);
    print "\n";  
  }

}

=head2 find_config_path

  Arg [1]   : Scalar, module name with .pm extension
  Function  : Loops through @INC to attempt to find the module along the standard config
              paths. Once the path is found it is returned.

  Returntype: Scalar, path to module. NULL if module not found.
  Exceptions: 
  Example   : find_config_path("BatchQueue.pm")

=cut

sub find_config_path
{
  my ($module) = @_;
  my $module_no_pm = $module;
  $module_no_pm =~ s/\.pm//;

  my $return_path;

  if($confirmed_config_paths{$module})
  {
    return($confirmed_config_paths{$module});
  }

  else
  {
    foreach my $path (@INC)
    {

      if($path =~ /ensembl\-config/)
      {

        if(-e $path.'/Bio/EnsEMBL/Analysis/Config/GeneBuild/'.$module)
        {
          $confirmed_config_paths{$module} = 'Bio::EnsEMBL::Analysis::Config::GeneBuild::'.$module_no_pm; 
        }

        elsif(-e $path.'/Bio/EnsEMBL/Analysis/Config/'.$module)
        {
          $confirmed_config_paths{$module} = 'Bio::EnsEMBL::Analysis::Config::'.$module_no_pm;           
        }

        elsif(-e $path.'/Bio/EnsEMBL/Pipeline/Config/'.$module)
        {
          $confirmed_config_paths{$module} = 'Bio::EnsEMBL::Pipeline::Config::'.$module_no_pm;
        }

      }
 
    } # End foreach

    return($confirmed_config_paths{$module});
 
  }

}

