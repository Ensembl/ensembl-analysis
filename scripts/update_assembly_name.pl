#!/usr/bin/env perl
#
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


# Update coord_system and meta tables with new assembly name

#

use strict;
use warnings;

use Getopt::Long;
use feature 'say';
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::CliHelper;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

#my $assembly_registry_host = $ENV{GBS1};
#my $assembly_registry_port = $ENV{GBP1};

#my $assembly_registry = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
#  -host    => $assembly_registry_host,
 # -port    => $assembly_registry_port,
 # -user    => 'ensro',
 # -dbname  => 'gb_assembly_registry');

my ($self) = @_;
#my $registry_assembly_name = $assembly_registry->fetch_assembly_name_by_gca($self->o('assembly_accession'));

if ( scalar(@ARGV) == 0 ) {
        usage();
        exit 0;
}

my $dbname = '';
my $host = '';
my $user = '';
my $port = '';
my $pass = '';
my $driver = '';
my $assembly_accession  = '';
my $assembly_name = '';
my $working_directory = '';
my $registry_host = '';
my $registry_port = '';
my $registry_db = '';

GetOptions('dbname:s' => \$dbname,
           'host:s'  => \$host,
           'user:s' => \$user,
           'port:s' => \$port,
           'pass:s' => \$pass,
           'driver:s' => \$driver,
           'assembly_accession:s' => \$assembly_accession,
           'assembly_name:s' => \$assembly_name,
           'working_dir:s' => \$working_directory,
           'registry_host:s' => \$registry_host,
           'registry_port:s' => \$registry_port,
           'registry_db:s' => \$registry_db,
          );

my $assembly_registry = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
  -host    => $registry_host,
  -port    => $registry_port,
  -user    => 'ensro',
  -dbname  => 'gb_assembly_registry');

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -dbname => $dbname,
  -host   => $host,
  -port   => $port,
  -user   => $user,
  -pass   => $pass,
  -driver => $driver,
);

my $registry_assembly_name = $assembly_registry->fetch_assembly_name_by_gca($assembly_accession);

#say "registry name is $registry_assembly_name and core name is $assembly_name";
if ($registry_assembly_name eq $assembly_name){
   say "nothing to update";
}
else{
  my $sth = $dba->dbc->prepare("UPDATE coord_system set version =?");
  $sth->bind_param(1,$registry_assembly_name);
  if ($sth->execute){
  }
  else{
    $self->throw("Could not update coord_system table");
  }
  my $sth1 = $dba->dbc->prepare("UPDATE meta set meta_value =? where meta_key = 'assembly.default'");
  $sth1->bind_param(1,$registry_assembly_name);
  if ($sth1->execute){
  }
  else{
    $self->throw("Could not update meta key assembly.default in meta table");
  }
  my $sth2 = $dba->dbc->prepare("UPDATE meta set meta_value =? where meta_key = 'assembly.name'");
  $sth2->bind_param(1,$registry_assembly_name);
  if ($sth2->execute){
  }
  else{
    $self->throw("Could not update meta_key assembly name in meta table");
  }
  &rename_bam_files($working_directory,$assembly_name,$registry_assembly_name);
}

sub rename_bam_files{
  my($dir, $prev_name, $new_name) = @_;
  $new_name =~ s/ /_/g;
  opendir(DH, $dir) || die "Can not open $dir: $!";
  my @files=readdir DH;
  close(DH);

  my $oldname;
  foreach(@files){
     $oldname=$_;
     s/$prev_name/$new_name/; # change $_ to new pattern
     next if(-e "$dir/$_");
     if(! rename "$dir/$oldname", "$dir/$_"){
        die("Could not rename $oldname to $_: $!");
     } else {
         print "File $oldname renamed to $_\n";
     }
   }
}

#------------------------------------------------------------------------------

sub usage {
	print <<EOF; exit(0);

Update coord_system and meta table with new assembly name. That is, if after pipeline started, the assembly name was changed.

Usage: perl $0 <options>

  --host|dbhost     Database host to connect to.

  --port|dbport     Database port to connect to (default is 3306).

  --user|dbuser     Database username.

  --pass|dbpass     Password for user.

  --driver|mysqldriver
  
  --assembly_accession|gca

  --assembly_name|name

  --working_dir|working_directory

  --help            This message.

  --verbose         Prints out the name of the database
                    which is going to be patched.




EOF
} ## end sub usage
