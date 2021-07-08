#!/usr/bin/env perl
#
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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


# Update assembly table with annotation status once genbuild is complete pre-handover

#

use strict;
use warnings;

use Getopt::Long;
use feature 'say';
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::CliHelper;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my ($self) = @_;

if ( scalar(@ARGV) == 0 ) {
        usage();
        exit 0;
}

my $user = '';
my $pass = '';
my $driver = '';
my $assembly_accession  = '';
my $registry_host = '';
my $registry_port = '';
my $registry_db = '';

GetOptions('user:s' => \$user,
           'pass:s' => \$pass,
           'driver:s' => \$driver,
           'assembly_accession:s' => \$assembly_accession,
           'registry_host:s' => \$registry_host,
           'registry_port:s' => \$registry_port,
           'registry_db:s' => \$registry_db,
          );

my $registry_dba = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
  -host    => $registry_host,
  -port    => $registry_port,
  -user    => $user,
  -dbname  => $registry_db,
  -pass    => $pass,
  -driver  => $driver,);

my $registry_assembly_id = $registry_dba->fetch_assembly_id_by_gca($assembly_accession);
my $sth = $registry_dba->dbc->prepare("UPDATE assembly set annotated_status =? where assembly_id=?");
$sth->bind_param(1,'completed');
$sth->bind_param(2,$registry_assembly_id);
if ($sth->execute){
}
else{
 $self->throw("Could not update annotation status");
}

#------------------------------------------------------------------------------

sub usage {
	print <<EOF; exit(0);

Update assembly table with annotation status once annotation is completed.

Usage: perl $0 <options>

  --user|dbuser     Database username.

  --pass|dbpass     Password for user.

  --driver|mysqldriver
  
  --assembly_accession|gca

  --help            This message.

  --verbose         Prints out the name of the database
                    which is going to be patched.




EOF
} ## end sub usage
