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

# Script to parse gff3 files for clones. At the moment based off parsing NCBI clone db files
# but should work with any correctly formatted file. It is assumed currently that there is
# order in the file, i.e:
# clone_insert line
# clone_insert_start line
# clone_insert_end line
# The reason is that the clone name from the clone_insert line needs to be used with the start/end lines
# Might change this in future
use warnings;
use strict;
use feature 'say';
use Getopt::Long;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::Parser::GFF3;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;

my $host;
my $port;
my $user;
my $pass;
my $dbname;
my $coord_system_name; # Coord system to get the slices from (toplevel usually)

my $prod_host;
my $prod_ruser = 'ensro';
my $prod_wuser;
my $prod_pass;
my $prod_port;
my $prod_dbname;
my $prod_backup_path; # where to backup the prod db to
my $prod_create_id = 106182; # id of genebuilder in the PRODUCTION DB (so not genebuilder id)
my $prod_name; # species.production_name, if not passed in this is gotten from the production db

my $species_name; # E.g 'mouse'
my $file; # The path the file with the clones
my $clone_state_file; # Path to state file. Get from ftp site e.g for mouse download ftp.ncbi.nih.gov/repository/clone/reports/Mus_musculus/clone_acstate_10090.out
my $set_code = "_clones"; # Will become 'rp24_clones' as code in misc_set
my $set_name = " clones"; # Will become 'RP24 mouse clones' as name in misc_set
my $set_desc = " library of clones for "; # Will become 'RP24 library of clones for mouse' as desc in misc_set
my $library_abbrev; # Clone library abbreviation, e.g. 'RP24'
my $clone_info_hash = {}; # Will store a variety of clone related info to pass stuff neatly into subroutine calls
my $set_code_char_lim = 25; # Current table character limit for the misc_set table code columm
my $write = 0; # For testing I guess
my $prod_write = 0;

GetOptions(
            'host=s'              => \$host,
            'port=s'              => \$port,
            'user=s'              => \$user,
            'pass=s'              => \$pass,
            'dbname=s'            => \$dbname,
            'coord_system_name'   => \$coord_system_name,
            'species_name=s'      => \$species_name,
            'prod_wuser'          => \$prod_wuser,
            'prod_pass'           => \$prod_pass,
            'prod_port'           => \$prod_port,
            'prod_dbname'         => \$prod_dbname,
            'prod_backup_path'    => \$prod_backup_path,
            'prod_create_id'      => \$prod_create_id,
            'species_prod_name'   => \$prod_name,
            'clone_file=s'        => \$file,
            'clone_state_file=s'  => \$clone_state_file,
            'set_code=s'          => \$set_code,
            'set_name=s'          => \$set_name,
            'set_desc=s'          => \$set_desc,
            'library_abbrev=s'    => \$library_abbrev,
            'ref_write!'          => \$write,
            'prod_write!'         => \$prod_write,
          );

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
           -host => $host,
           -user => $user,
           -port => $port,
           -pass => $pass,
           -dbname => $dbname
         );

my $prod_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                -host => $prod_host,
                -user => $prod_wuser,
                -port => $prod_port,
                -pass => $prod_pass,
                -dbname => $prod_dbname
              );


my $slice_adaptor = $db->get_SliceAdaptor;
my $misc_feature_adaptor = $db->get_MiscFeatureAdaptor();
my $simple_feature_adaptor = $db->get_SimpleFeatureAdaptor();
my $misc_set_adaptor_prod = $prod_db->get_MiscSetAdaptor();
my $misc_set_adaptor_core = $db->get_MiscSetAdaptor();

# The library abbreviation is usually the first set of characters before a dot in the file name
# So if the user does not define a library abbreviation then the script will try and pull it from
# the file name. If it still can't get it then throw
unless($library_abbrev) {
  say "FILE: ".$file;
  $file =~ /^([^\.]+)/;
  $library_abbrev = $1;
  unless($library_abbrev) {
  throw("No library abbreviation was provided and the script failed to parse the library abbreviation ".
        "from the file name. It should follow standard NCBI format, e.g:\n".
        "bMQ.GCF_000001635.22.103.unique_concordant.gff\n\nor\n\n".
        "bMQ.GCF_000001635.22.103.unique_discordant.gff\n\nor\n\n".
        "bMQ.GCF_000001635.22.103.multiple.gff");
  }
}

# Code for misc_set entry. Will take the form 'mouse_rp24_clones'. The char limit on the table is currently 25
# so the code will throw if it is too long
$set_code = $library_abbrev.$set_code;
if(length($set_code) > $set_code_char_lim) {
  throw("The code generated for the misc set was ".length($set_code)." characters long, which is greater than ".
        "the current limit for the table, which is ".$set_code_char_lim." characters. The generated misc set code was:\n".$set_code);
}

# Similar to the above, but there shouldn't be any space issues
$set_name = $library_abbrev." ".$species_name.$set_name;
$set_desc = $library_abbrev.$set_desc.$species_name;

# Build a misc set object based on the above info
my $clone_set = create_clone_set($set_code,$set_name,$set_desc,$misc_set_adaptor_prod);

# This hash ($clone_info_hash) is just an easy way of putting a lot of the data together
# so that the arguments list for calling the subroutines are smaller
set_clone_info_hash($clone_info_hash,$file,$species_name,$library_abbrev,$clone_set);

# If a clone state file is provided then load it and put the states into a hash so they
# can be added as attributes
my $clone_states;
if($clone_state_file) {
  $clone_states = load_clone_states($clone_state_file);
} else {
  warning("You have not provided a clone state file, so no states will be added for the clones.\n".
          "This file is usually available on the ftp site. For examples the mouse combined file is:\n".
          "https://ftp.ncbi.nih.gov/repository/clone/reports/Mus_musculus/clone_acstate_10090.out");
}

# This is a hash that will hold the parent inserts for the ends based on the parent id
# It is used to put the clone name of the parent into the display label of the simple feature
my $parent_clones;

# This will be the array that holds both the MiscFeature and SimpleFeature objects. The contents
# of this array will be stored at the end
my $clone_array;

# Open gff3 file and loop through it
my $gff_file = Bio::EnsEMBL::IO::Parser::GFF3->open($file);
while($gff_file->next) {
  my $type = $gff_file->get_type;
  unless($type) {
    next;
  }

  if($type eq 'clone_insert' || $type eq 'clone_insert_start' || $type eq 'clone_insert_end') {

    my $seq_region_name = $gff_file->get_seqname;
    $seq_region_name =~ s/^ref\|(.+)\|$/$1/;

    my $clone_name = $gff_file->get_attribute_by_name('Name');

    # This section of code checks to see if the line is a child line, if it is it changes $clone_name
    # to the value of the parent in the parent_clones hash. Otherwise the line itself is a parent and
    # the details are added to the parent hash. Note that as is it assume file order, but there is a
    # warning to at least point out when this happens
    my $clone_parent = $gff_file->get_attribute_by_name('Parent');
    if($clone_parent) {
      my $child_clone_name = $clone_name;
      $clone_name = $parent_clones->{$clone_parent};
      unless($clone_name) {
       warning("Failed to find the name of parent for: ".$child_clone_name);
      }
    } else {
      my $clone_id = $gff_file->get_attribute_by_name('ID');
      $parent_clones->{$clone_id} = $clone_name;
    }

    # Make a slice for the seq region name. Sometimes we might not have these seq regions, as they might
    # correspond to non-reference strains. For now just skip over with a warning. Will implement a count
    # at some point so that we can check for any abnormal behaviour
    my $slice = $slice_adaptor->fetch_by_region($coord_system_name,$seq_region_name);
    unless($slice) {
      warning("No slice found in the db for the following: ".$seq_region_name);
      next;
    }

    my $start = $gff_file->get_start;
    my $end  = $gff_file->get_end;
    my $strand = $gff_file->get_strand;
    unless(defined($start) && defined($end) && defined($strand)) {
      throw("Could not parse out start/end/strand for the following: ".$seq_region_name);
    }

    # Make a clone feature
    my $clone_feature = create_clone_feature($clone_info_hash,
                                             $type,
                                             $start,
                                             $end,
                                             $strand,
                                             $slice,
                                             $clone_name);

    # If this is a misc feature then add all the associated attributes
    if(ref($clone_feature) eq "Bio::EnsEMBL::MiscFeature") {
      add_clone_attribs($clone_feature,$clone_info_hash,$clone_name,$clone_states);
    }

    push(@{$clone_array},$clone_feature);

  } else {
    throw("Unrecognised clone type: ".$type);
  }
}

# Store all the features in the db as either misc or simple features
if($write) {
  if($prod_write) {
    backup_production_db($prod_backup_path);
    update_production_db($misc_set_adaptor_prod,$clone_info_hash,$prod_create_id,$prod_name);
  }
  store_clone_features($clone_array,$misc_feature_adaptor,$simple_feature_adaptor,$clone_info_hash);
}

exit;


sub backup_production_db {
  my ($prod_backup_path) = @_;
  my $timestamp = time();
  my $return_value = system('mysqldump -u'.$prod_ruser.' -h'.$prod_host.' '.$prod_dbname.' > '.
                            $prod_backup_path.'/prod_db_backup.'.$timestamp.'.sql');

  unless($return_value == 0 && -e $prod_backup_path.'/prod_db_backup.'.$timestamp.'.sql') {
    throw("Failed to backup the production db. Tried to back them up to the following dir:\n".$prod_backup_path."\n".
          "Potential name of failed backup file:\n".'prod_db_backup.'.$timestamp.'.sql');
  }
}

sub update_production_db {
  my ($misc_set_adaptor_prod,$clone_info_hash,$prod_create_id,$prod_name) = @_;
  # Store the misc_set
  $misc_set_adaptor_prod->store($clone_info_hash->{'misc_set'});

  my $logic_name = $clone_info_hash->{'analysis'}->logic_name;
  my $sth_select = $prod_db->dbc->prepare('SELECT logic_name from analysis_description where logic_name = ?');
  $sth_select->bind_param(1,$logic_name);
  $sth_select->execute();
  my ($result_row) = $sth_select->fetchrow_array();
  unless($result_row) {
    my $sth_insert = $prod_db->dbc->prepare('INSERT INTO analysis_description '.
                                            '(logic_name,description,display_label,db_version,is_current,created_by,created_at,default_web_data_id,default_displayable) '.
                                            'VALUES(?,?,?,1,1,?,NOW(),NULL,1)');
    $sth_insert->bind_param(1, $logic_name);
    $sth_insert->bind_param(2, $clone_info_hash->{'analysis'}->description);
    $sth_insert->bind_param(3, $clone_info_hash->{'analysis'}->display_label);
    $sth_insert->bind_param(4, $prod_create_id);
    $sth_insert->execute();
  } else {
    warning("The logic_name '".$logic_name."' was already present in the production databases so it was not added again");
  }

  unless($prod_name) {
    $sth_select = $db->dbc->prepare("SELECT meta_value from meta where meta_key='species.production_name'");
    $sth_select->execute();
    ($result_row) = $sth_select->fetchrow_array();
    unless($result_row) {
    throw("Failed to find the species production name from the meta table and no value was set by user. SQL used:\n".
          "SELECT meta_value from meta where meta_key='species.production_name'\nYou must either set the prod_name ".
          "variable or have the species.production_name key set in the meta table of the core db");
    }
    $prod_name = $result_row;
  }

  $sth_select = $prod_db->dbc->prepare('SELECT species_id from species where production_name = ?');
  $sth_select->bind_param(1,$prod_name);
  $sth_select->execute();
  ($result_row) = $sth_select->fetchrow_array();
  unless($result_row) {
    throw("Failed to find the species id from the species table in the production db. The table was queried with the following:\n".
          "SELECT species_id from species where production_name = '".$prod_name."';");
  }

  my $species_id = $result_row;

  $sth_select = $prod_db->dbc->prepare('SELECT analysis_description_id from analysis_description where logic_name = ?');
  $sth_select->bind_param(1,$logic_name);
  $sth_select->execute();
  ($result_row) = $sth_select->fetchrow_array();

  unless($result_row) {
    throw("Failed to retrieve an analysis_description_id for '".$logic_name."'. The table was queried with the following:\n".
          "SELECT analysis_description_id from analysis_description where logic_name = '".$logic_name."';");
  }

  my $analysis_description_id = $result_row;
  my $sth_insert = $prod_db->dbc->prepare('INSERT INTO analysis_web_data '.
                                          '(analysis_description_id,web_data_id,species_id,db_type,displayable,created_by,created_at) '.
                                          'VALUES(?,NULL,?,"core",1,?,NOW())');
  $sth_insert->bind_param(1, $analysis_description_id);
  $sth_insert->bind_param(2, $species_id);
  $sth_insert->bind_param(3, $prod_create_id);
  $sth_insert->execute();

}

# This subroutine sets a bunch of different general info related to the clone file
# into the clone info hash. It just makes it easier to pass this hash in later subroutine
# calls than passing all the individual arguments
sub set_clone_info_hash {
  my ($clone_info_hash,$file,$species_name,$library_abbrev,$clone_set) = @_;

  $clone_info_hash->{'library_abbrev'} = $library_abbrev;

  # Set concordance. Can be unique_concordant, unique_discordant or multiple. If this info is not
  # present then there is just a warning at the moment, the info itself is currently unused anyway,
  # ideally it will get added as a new attribute since the current ones in the webcode for MiscFeature
  # don't fit too well
  unless($file =~ /\.([^\.]+)\.gff$/) {
    warning("Error could not parse the file name for concordance info. It should follow standard NCBI format, e.g:\n".
            "bMQ.GCF_000001635.22.103.unique_concordant.gff\nor\n".
            "bMQ.GCF_000001635.22.103.unique_discordant.gff\nor\n".
            "bMQ.GCF_000001635.22.103.multiple.gff\n".
            "As this is not a requirement the script will continue without the info");
  } else {
    my @concordance = split('_',$1);
    if(scalar(@concordance) == 2) {
      $clone_info_hash->{'unique'} = 'yes';
      $clone_info_hash->{'concordance'} = $concordance[1];
    } else {
      $clone_info_hash->{'unique'} = 'no';
    }
  }

  # Set the logic_name. Here if you start off with mouse as the logic_name and the library is bMQ, you
  # get 'mouse_bMQ_clones'. Once this is done make an analysis object and add to hash
  unless($species_name) {
    throw("You should provide a generic logic name for the BAC ends as they are simple features. ".
          "Using the species name (e.g 'mouse') is enough as the library will be appended automatically");
  }

  my $logic_name = $species_name."_".$library_abbrev."_clones";
  $logic_name = lc($logic_name);
#  my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
#  if (!defined $analysis) {
  my $analysis = Bio::EnsEMBL::Analysis->new(
                   -logic_name => $logic_name,
                   -description => $library_abbrev." library clone ends for ".$species_name,
                   -display_label => $library_abbrev." ".$species_name." clone ends",
                   -displayable => 1,
                   -gff_source => "NCBI",
                   -gff_feature => "clone insert start/end",
                   -db_file => $file,
                   -db => "NCBI clone database",
                 );
#  }
  $clone_info_hash->{'analysis'} = $analysis;

  # Add clone_set
  $clone_info_hash->{'misc_set'} = $clone_set;

  return($clone_info_hash);
}

# This loads the clone states file. The second column in the file is the name, while the sixth
# column is the state. Remember that opening the file in emacs could lead to the removal of tabs,
# which is what the file is split on and would cause the script to throw
sub load_clone_states {
  my ($clone_state_file) = @_;

  my $clone_states = {};
  unless(open(IN,$clone_state_file)) {
    throw("You provided a clone state file but it could not be opened successfully:\n".$clone_state_file);
  }

  while(<IN>) {
    my $line = $_;
    if($line =~ /^\#/) {
      next;
    }

    my @line_ele = split("\t",$line);
    unless(scalar(@line_ele) == 10) {
     warning("Parsed a line in the clone state file that had more/less than 10 columns (skipping):\n".$line);
     next;
    }

    my $clone_name = $line_ele[1];
    my $clone_state = $line_ele[5];
    $clone_states->{$clone_name} = $clone_state;
  }
  close IN;

  return($clone_states);
}

# This creates a MiscSet object for the clones
sub create_clone_set {
  my ($set_code,$set_name,$set_desc,$misc_set_adaptor_prod) = @_;
  unless($set_code) {
    throw("Failed to create a set for the clones, you must at least have a set_code defined!");
  }

  my $clone_set = $misc_set_adaptor_prod->fetch_by_code($set_code);
  if($clone_set) {
    return $clone_set;
  } else {
    # If the set didn't exist then create one
    $clone_set = Bio::EnsEMBL::MiscSet->new(
                   -CODE        => $set_code,
                   -NAME        => $set_name,
                   -DESCRIPTION => $set_desc,
                 );
    return $clone_set;
  }
}

# This creates either a MiscFeature of a SimpleFeature depending on whether the line refers
# to the insert (MiscFeature) or the ends (SimpleFeature)
sub create_clone_feature {
  my ($clone_info_hash,$type,$start,$end,$strand,$slice,$clone_name) = @_;

  if($type eq 'clone_insert') {
    my $clone_feature = Bio::EnsEMBL::MiscFeature->new(
                          -START  => $start,
                          -END    => $end,
                          -STRAND => $strand,
                          -SLICE  => $slice,
                        );
#    $clone_feature->add_MiscSet($clone_info_hash->{'misc_set'});
    return $clone_feature;
  } else {
    $type =~ /([^\_]+)$/;
    my $subtype = $1;
    my $clone_feature = Bio::EnsEMBL::SimpleFeature->new(
                          -ANALYSIS => $clone_info_hash->{'analysis'},
                          -START  => $start,
                          -END    => $end,
                          -STRAND => $strand,
                          -SLICE  => $slice,
                          -DISPLAY_LABEL => $clone_name." (".$subtype.")",
                        );
    return $clone_feature;
  }
}

# This adds attributes to the MiscFeatures, therefore stored in misc_attrib. The codes themselves were used
# based on the fact that ensembl-webcode has a MiscFeature module that has an array for the clones. As such
# the codes used were based on this array, even if it is perhaps not the ideal setup (for example the label
# "Library:" on the zmenu is under the code "clone_name". In future we should update the table to make more
# sense than what's currently there
sub add_clone_attribs {
  my ($clone_feature,$clone_info_hash,$clone_name,$clone_states) = @_;

  my $attrib_name = Bio::EnsEMBL::Attribute->new(
                        -CODE  => 'name',
                        -VALUE => $clone_name,
                    );

  my $attrib_library = Bio::EnsEMBL::Attribute->new(
                         -CODE  => 'clone_name',
                         -VALUE => $clone_info_hash->{'library_abbrev'},
                       );

  $clone_feature->add_Attribute($attrib_name);
  $clone_feature->add_Attribute($attrib_library);

  if($clone_states->{$clone_name}) {
    my $attrib_state = Bio::EnsEMBL::Attribute->new(
                         -CODE  => 'state',
                         -VALUE => $clone_states->{$clone_name},
                       );
    $clone_feature->add_Attribute($attrib_state);
  }

  if($clone_name =~ /([^\-]+)$/) {
    my $well_name = $1;
    my $attrib_well = Bio::EnsEMBL::Attribute->new(
                         -CODE  => 'well_name',
                         -VALUE => $well_name,
                       );
    $clone_feature->add_Attribute($attrib_well);
  }

  return $clone_feature;
}

# Store sub to store either a MiscFeature or a SimpleFeature
sub store_clone_features {
  my ($clone_array_ref,$misc_feature_adaptor,$simple_feature_adaptor,$clone_info_hash) = @_;

  my $clone_feature;
  my $misc_set_code = $clone_info_hash->{'misc_set'}->code;
  my $misc_set_prod = $misc_set_adaptor_prod->fetch_by_code($misc_set_code);
  my $misc_set_prod_id = $misc_set_prod->dbID;
  $misc_set_adaptor_core->store($misc_set_prod);

  my $sth_update = $db->dbc->prepare('UPDATE misc_set set misc_set_id = ? where code = ?');
  $sth_update->bind_param(1, $misc_set_prod_id);
  $sth_update->bind_param(2, $misc_set_code);
  $sth_update->execute();

  my $misc_set_core = $misc_set_adaptor_core->fetch_by_dbID($misc_set_prod_id);

  while($clone_feature = pop(@{$clone_array_ref})) {
    if(ref($clone_feature) eq "Bio::EnsEMBL::MiscFeature") {
     $clone_feature->add_MiscSet($misc_set_core);
     $misc_feature_adaptor->store($clone_feature);
    } elsif(ref($clone_feature) eq "Bio::EnsEMBL::SimpleFeature") {
      $simple_feature_adaptor->store($clone_feature);
    } else {
      throw("Unexpected datastructure found on store. Should be either:\n".
            "\nBio::EnsEMBL::MiscFeature\n".
            "\nor\n".
            "Bio::EnsEMBL::SimpleFeature\n".
            "\nFound:\n".
            ref($clone_feature));
    }
  }

  my $sth_select = $prod_db->dbc->prepare('SELECT misc_set_id from misc_set where code = ?');
  $sth_select->bind_param(1, $misc_set_code);
  $sth_select->execute();
  my ($prod_misc_set_id) = $sth_select->fetchrow_array();

}
