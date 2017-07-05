#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetMetaAndSeqRegionSynonym;

use strict;
use warnings;
use feature 'say';
use Time::Piece;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
                chromosomes_present => 0,
                has_mitochondria => 0,
 Returntype : Hashref, the default parameter for this module
 Exceptions : None

=cut

sub param_defaults {
  my $self = shift;

  return {
    %{$self->SUPER::param_defaults},
    chromosomes_present => 0,
    has_mitochondria => 0,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Check that 'target_db' and 'enscode_root_dir' are set
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db flag not passed into parameters hash. The target db to load the assembly info ".
                 "into must be passed in with write access");
  }

  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_root_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  return 1;
}


=head2 run

 Arg [1]    : None
 Description: backup the meta table, add information in the meta table and add seq_region_synonyms
 Returntype : None
 Exceptions : None

=cut

sub run {
  my $self = shift;

  say "Loading meta information seq region synonyms into reference db\n";
  my $target_db = $self->param_required('target_db');
  my $enscode_dir = $self->param_required('enscode_root_dir');
  my $primary_assembly_dir_name = $self->param_required('primary_assembly_dir_name');
  my $path_to_files = $self->param_required('output_path')."/".$primary_assembly_dir_name;
  my $meta_keys = $self->param_required('meta_key_list');
  say "\nBacking up meta and seq_region tables...";
  $self->backup_tables($path_to_files,$target_db);
  say "\nBackup of tables complete\n";

  say "Setting meta information in meta table...\n";
  $self->set_meta($target_db,$meta_keys,$path_to_files);
  say "\nMeta table insertions complete\n";

  say "Setting seq region synonyms...\n";
  $self->set_seq_region_synonyms($target_db,$path_to_files);
  say "\nSeq region synonyms inserted\n";

  say "\nFinished updating meta table and setting seq region synonyms";
  return 1;
}


=head2 write_output

 Arg [1]    : None
 Description: It will get the job related input ids, add 'mt_accession' and 'chromosomes_present'
              if needed and flow the data on branch #1
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

# This code should get the upstream input and add 'mt_accession' and 'chromosomes_present' if needed
  my $job_params = eval($self->input_job->input_id);
  if ($self->param_is_defined('mt_accession')) {
    $job_params->{mt_accession} = $self->param('mt_accession');
  }
  $self->dataflow_output_id($job_params, 1);
  return 1;
}


=head2 backup_tables

 Arg [1]    : String $path_to_files, directory
 Arg [2]    : Hashref $target_db, connection details
 Description: It will backup the meta table in the directory specified
              in Arg[1] in the format table_name.time.sql
 Returntype : None
 Exceptions : None

=cut

sub backup_tables {
  my ($self,$path_to_files,$target_db) = @_;

  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  for my $table ('seq_region','meta') {
    my $backup_file = $path_to_files."/".$table.".".time().".sql";
    my $cmd = "mysqldump".
              " -h".$dbhost.
              " -P".$dbport.
              " -u".$dbuser.
              " -p".$dbpass.
              " ".$dbname.
              " ".$table.
              " > ".$backup_file;
    my $return = system($cmd);
    if($return) {
      $self->throw("mysqldump to backup ".$table." table failed. Commandline used:\n".$cmd);
    } else {
      say $table." table backed up in the following location:\n".$backup_file;
    }
  }
}


=head2 set_meta

 Arg [1]    : Hashref $target_db, DB connection detail
 Arg [2]    : Int $genebuilder_id, numerical id of the person doing the annotation
 Arg [3]    : String $path_to_files, directory containing the file assembly_report.txt
 Description: Parse the assembly_report.txt file to get information about the assemlby date
              the taxon id, the assembly accession, the assembly level (chromosome, scaffold,...)
              and the presence of a mitochondria annotated by RefSeq
 Returntype : None
 Exceptions : None

=cut

sub set_meta {
  my ($self,$target_db,$meta_keys,$path_to_files) = @_;

  my $target_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$target_db});
  my $meta_adaptor = $target_dba->get_MetaContainerAdaptor;

  foreach my $meta_key (keys(%{$meta_keys})) {
    my $meta_value = $meta_keys->{$meta_key};
    $meta_adaptor->store_key_value($meta_key, $meta_value);
    say "Inserted into meta:\n".$meta_key." => ".$meta_value;
  }

  my $date = localtime->strftime('%Y-%m');
  $meta_adaptor->store_key_value('genebuild.start_date', $date."-Ensembl");

  unless(-e $path_to_files."/assembly_report.txt") {
    $self->throw("Could not find the assembly_report.txt file. Path checked:\n".$path_to_files."/assembly_report.txt");
  }

  open(IN,$path_to_files."/assembly_report.txt");
  while (my $line = <IN>) {
    # This check seems odd, might mess up if run on plants
    if($line !~ /^#/ or $self->param('has_mitochondria')) {
      if ($line =~ /(NC_\S+)\s+non-nuclear/) {
        $self->param('mt_accession', $1);
      }
    } elsif($line =~ /^#\s*Date:\s*(\d+)-(\d+)-\d+/) {
      $meta_adaptor->store_key_value('assembly.date', $1.'-'.$2);
      say "Inserted into meta:\nassembly.date => ".$1.'-'.$2;
   } elsif ($line =~ /^#\s*Assembly level:\s*[Cc]hromosome/) {
      $self->param('chromosomes_present', 1);
    } elsif ($line =~ /##\s*GCF_\S+\s*non-nuclear/) {
      $self->param('has_mitochondria', 1);
    }
  }

  close IN;
}


=head2 set_seq_region_synonyms

 Arg [1]    : Hashref $target_db, DB connection details
 Arg [2]    : String $path_to_files, directory name
 Description: Set the INSDC seq_region_syonyms for the chromosomes if present and
              the submitter's synonyms for the contigs and scaffolds
 Returntype : None
 Exceptions : None

=cut

sub set_seq_region_synonyms {
  my ($self,$target_db,$path_to_files) = @_;

  my $target_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$target_db});

  if($self->param('chromosomes_present')) {
    unless(-e $path_to_files."/chr2acc") {
      $self->throw("Could not find chr2acc file. No chromosome synonyms loaded. Expected location:\n".$path_to_files."/chr2acc");
    }

    open(IN,$path_to_files."/chr2acc");
    my $sth_select = $target_dba->dbc->prepare('SELECT sr.seq_region_id FROM seq_region sr, coord_system cs WHERE cs.coord_system_id = sr.coord_system_id AND sr.name = ? AND cs.rank = 1');
    my $sth_insdc = $target_dba->dbc->prepare('SELECT external_db_id FROM external_db WHERE db_name = "INSDC"');
    $sth_insdc->execute();
    my ($insdc_db_id) = $sth_insdc->fetchrow_array;

    my $sth_insert = $target_dba->dbc->prepare('INSERT INTO seq_region_synonym (seq_region_id, synonym, external_db_id) VALUES(?, ?, ?)');
    my $sth_update = $target_dba->dbc->prepare('UPDATE seq_region set name = ? WHERE seq_region_id = ?');
    my $insert_count = 0;
    while(my $line = <IN>) {
      if($line =~ /^#/) {
        next;
      }
      my ($synonym, $seq_region_name) = $line =~ /(\w+)\s+(\S+)/;
      $sth_select->bind_param(1, $seq_region_name);
      $sth_select->execute();
      my ($seq_region_id) = $sth_select->fetchrow_array();
      $sth_insert->bind_param(1, $seq_region_id);
      $sth_insert->bind_param(2, $seq_region_name);
      $sth_insert->bind_param(3, $insdc_db_id);
      $sth_insert->execute();
      $sth_update->bind_param(1, $synonym);
      $sth_update->bind_param(2, $seq_region_id);
      $sth_update->execute();
      $insert_count++;
    }
    close(IN);

    if($insert_count == 0) {
      $self->throw("The insert/update count after parsing chr2acc was 0, this is probably wrong. File used:\n".$path_to_files."/chr2acc");
    }

    say "\nInserted into seq_region_synonym and updated seq_region based on chr2acc. Total inserts/updates: ".$insert_count;
  } else {
    say "The chromosomes_present parameter was not set to 1 in the config, so assuming there are no chromosomes";
  }

  foreach my $file ('component_localID2acc', 'scaffold_localID2acc') {
    unless(-e $path_to_files."/".$file) {
      $self->throw("Could not find ".$file." file. No synonyms loaded. Expected location:\n".$path_to_files."/".$file);
    }

    open(IN,$path_to_files."/".$file);
    my $sth_select = $target_dba->dbc->prepare('SELECT seq_region_id FROM seq_region WHERE name = ?');
    my $sth_insert = $target_dba->dbc->prepare('INSERT INTO seq_region_synonym (seq_region_id, synonym) VALUES(?, ?)');
    my $insert_count = 0;
    while (my $line = <IN>) {
      if ($line =~ /^#/) {
        next;
      }

      my ($synonym, $seq_region_name) = $line =~ /(\S+)\s+(\S+)/;
      if($synonym eq 'na') {
        $synonym = $seq_region_name;
      }
      $sth_select->bind_param(1, $seq_region_name);
      $sth_select->execute();
      my ($seq_region_id) = $sth_select->fetchrow_array();
      $sth_insert->bind_param(1, $seq_region_id);
      $sth_insert->bind_param(2, $synonym);
      $sth_insert->execute();
      $insert_count++;
    }

    if($insert_count == 0) {
      $self->throw("The insert/update count after parsing ".$file." was 0, this is probably wrong. File used:\n".$path_to_files."/".$file);
    }

    say "\nInserted into seq_region_synonym and updated seq_region based on ".$file.". Total inserts/updates: ".$insert_count;
    say "You will need to update the external_db_id for the synonyms of scaffold or contigs!\n";
  }

}

1;
