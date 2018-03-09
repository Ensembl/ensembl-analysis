#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
use strict;
use warnings;

use Cwd;
use File::Spec::Functions qw(catdir updir catfile splitpath);
use Test::More;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_database_from_registry);

use_ok('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase');

my %cloned_tables = map {$_ => 1} @{Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase->vital_tables()};

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor('core');
my $test_dbname = $multi->create_db_name;
my %target_db = (
  -dbname => $test_dbname,
  -host   => $db->dbc->host,
  -port   => $db->dbc->port,
  -user   => $db->dbc->user,
  -pass   => $db->dbc->pass,
  -driver => $db->dbc->driver,
);

my $enscode_root_dir;
if (exists $ENV{ENSCODE}) {
  $enscode_root_dir = $ENV{ENSCODE};
}
foreach my $path ($enscode_root_dir, cwd(), catdir(cwd(), updir()), catdir(cwd(), updir(), updir(), updir())) {
  if (-d catdir($path, 'ensembl')) {
    $enscode_root_dir = $path;
    last;
  }
}
my @split_dir = splitpath(cwd());
while ($split_dir[-1] ne 'ensembl-analysis') {
  pop(@split_dir);
}
my $analysis_dir = catdir(@split_dir);

###
# Testing create_type core_only
###
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
	{
    create_type => 'core_only',
    target_db => \%target_db,
    enscode_root_dir => $enscode_root_dir,
	},
);
$db->dbc->do("use $test_dbname");
my $sth = $db->dbc->prepare('SHOW TABLES');
$sth->execute();
cmp_ok(@{$sth->fetchall_arrayref}, '>=', 73, 'Checking all tables loaded');

###
# Testing create_type backup
###
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
	{
    create_type => 'backup',
    source_db => \%target_db,
    output_path => cwd(),
    backup_name => $test_dbname.'.sql',
	},
);
my $gzip = catfile(cwd(), $test_dbname.'.sql.gz');
ok(-e $gzip, "Testing backup file exists");
ok(-s $gzip, 'Testing backup file is not null');
my @stat = stat($gzip);
cmp_ok($stat[7], '>=', 8000, 'Checking file has data');
$db->dbc->do("DROP DATABASE $test_dbname");

###
# Testing create_type clone
###
#my $core_db = get_database_from_registry('sus_scrofa', 'Core');
my %source_db = (
  -dbname => $db->dbc->dbname,
  -host   => $db->dbc->host,
  -port   => $db->dbc->port,
  -user   => $db->dbc->user,
  -pass   => $db->dbc->pass,
  -driver => $db->dbc->driver,
);
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
	{
    create_type => 'clone',
    source_db => \%source_db,
    target_db => \%target_db,
	},
);
$db->dbc->do("use $test_dbname");
$sth = $db->dbc->prepare('SHOW TABLES');
$sth->execute();
my $tables = $sth->fetchall_arrayref;
cmp_ok(@$tables, '>=', 73, 'Checking all tables loaded');
my $cloned_count = 0;
my $empty_count = 0;
foreach my $table (@$tables) {
  my $query = 'SELECT COUNT(*) FROM '.$table->[0];
  $sth = $db->dbc->prepare($query);
  $sth->execute();
  my $data = $sth->fetchall_arrayref;
  if (exists $cloned_tables{$table->[0]}) {
    if ($data->[0]->[0] > 0) {
      ++$cloned_count;
    }
    else {
      $db->dbc->do('use '.$db->dbc->dbname);
      my $old_sth = $db->dbc->prepare($query);
      $sth->execute();
      my $old_data = $sth->fetchall_arrayref;
      if ($old_data->[0]->[0] == $data->[0]->[0]) {
        ++$cloned_count;
      }
      $db->dbc->do("use $test_dbname");
    }
  }
  elsif ($data->[0]->[0] == 0) {
    ++$empty_count;
  }
}
cmp_ok($cloned_count, '==', scalar(keys %cloned_tables), 'Checking vital tables are populated');
cmp_ok($empty_count+$cloned_count, '==', scalar(@$tables), 'Checking all other tables are empty');
$sth = $db->dbc->prepare('INSERT INTO protein_align_feature (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, hit_start, hit_end, analysis_id) VALUES (20, 1, 45, 1, "A1B2C3.1", 1, 15, 1)');
$sth->execute();
$sth = $db->dbc->prepare('SELECT protein_align_feature_id FROM protein_align_feature WHERE hit_name = "A1B2C3.1"');
$sth->execute();
my $data = $sth->fetchall_arrayref;
cmp_ok($data->[0]->[0], '==', 1, 'Checking autoincrement is reset to lower value');
$db->dbc->do("DROP DATABASE $test_dbname");

###
# Testing create_type dna_db
###
%source_db = (
  -dbname => $db->dbc->dbname,
  -host   => $db->dbc->host,
  -port   => $db->dbc->port,
  -user   => $db->dbc->user,
  -pass   => $db->dbc->pass,
  -driver => $db->dbc->driver,
);
standaloneJob(
	'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
	{
    create_type => 'dna_db',
    source_db => \%source_db,
    target_db => \%target_db,
	},
);
$db->dbc->do("use $test_dbname");
$sth = $db->dbc->prepare('SHOW TABLES');
$sth->execute();
$tables = $sth->fetchall_arrayref;
cmp_ok(@$tables, '>=', 73, 'Checking all tables loaded');
my $dna_count = 0;
$empty_count = 0;
my %dna_tables = (
  %cloned_tables,
  dna => 1,
  repeat_feature => 1,
  repeat_consensus => 1,
);
foreach my $table (@$tables) {
  my $query = 'SELECT COUNT(*) FROM '.$table->[0];
  $sth = $db->dbc->prepare($query);
  $sth->execute();
  my $data = $sth->fetchall_arrayref;
  if (exists $dna_tables{$table->[0]}) {
    if ($data->[0]->[0] > 0) {
      ++$dna_count;
    }
    else {
      $db->dbc->do('use '.$db->dbc->dbname);
      my $old_sth = $db->dbc->prepare($query);
      $sth->execute();
      my $old_data = $sth->fetchall_arrayref;
      if ($old_data->[0]->[0] == $data->[0]->[0]) {
        ++$dna_count;
      }
      $db->dbc->do("use $test_dbname");
    }
  }
  elsif ($data->[0]->[0] == 0) {
    ++$empty_count;
  }
}
cmp_ok($dna_count, '==', scalar(keys %dna_tables), 'Checking vital tables are populated');
cmp_ok($empty_count+$dna_count, '==', scalar(@$tables), 'Checking all other tables are empty');
$sth = $db->dbc->prepare('INSERT INTO protein_align_feature (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, hit_start, hit_end, analysis_id) VALUES (20, 1, 45, 1, "A1B2C3.1", 1, 15, 1)');
$sth->execute();
$sth = $db->dbc->prepare('SELECT protein_align_feature_id FROM protein_align_feature WHERE hit_name = "A1B2C3.1"');
$sth->execute();
$data = $sth->fetchall_arrayref;
cmp_ok($data->[0]->[0], '==', 1, 'Checking autoincrement is reset to lower value');

done_testing();

cleaning($db, $test_dbname);

sub cleaning {
  my ($db, $test_dbname) = @_;

  $db->dbc->do("DROP DATABASE $test_dbname");
  unlink $test_dbname.'.sql.gz';
}
