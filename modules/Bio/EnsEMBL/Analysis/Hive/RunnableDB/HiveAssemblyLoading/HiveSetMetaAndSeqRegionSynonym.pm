#!/usr/bin/env perl

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetMetaAndSeqRegionSynonym;

use strict;
use warnings;
use feature 'say';
use Time::Piece;

use File::Spec::Functions qw (catfile catdir);
use Bio::EnsEMBL::Hive::Utils qw(destringify);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
                chromosomes_present => 0,
                has_mitochondria => 0,
                chromosome_level => 'chromosome',
                scaffold_level => 'scaffold',
                contig_level => 'contig',
 Returntype : Hashref, the default parameter for this module
 Exceptions : None

=cut

sub param_defaults {
  my $self = shift;

  return {
    %{$self->SUPER::param_defaults},
    chromosomes_present => 0,
    has_mitochondria => 0,
    chromosome_level => 'chromosome',
    scaffold_level => 'scaffold',
    contig_level => 'contig',
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
  my $target_db = $self->get_database_by_name('target_db');
  my $enscode_dir = $self->param_required('enscode_root_dir');
  my $path_to_files = catdir($self->param('output_path'), $self->param_required('primary_assembly_dir_name'));
  my $meta_key_list = $self->param_required('meta_key_list');

  say "\nBacking up meta and seq_region tables...";
  $self->backup_tables($path_to_files,$target_db);
  say "\nBackup of tables complete\n";

  say "Setting meta information in meta table...\n";
  $self->set_meta($target_db,$meta_key_list,$path_to_files);
  say "\nMeta table insertions complete\n";

  if(-e catdir($path_to_files), 'AGP') {
    say "Setting seq region synonyms...\n";
    $self->set_seq_region_synonyms($target_db,$path_to_files);
    say "\nSeq region synonyms inserted\n";
  } else {
    $self->warning("Could not find an AGP dir, so assuming assembly is single level. Will not load synonyms");
  }

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
  my $job_params = destringify($self->input_job->input_id);
  if ($self->param_is_defined('mt_accession')) {
    $job_params->{mt_accession} = $self->param('mt_accession');
  }
  $self->input_job->input_id($job_params);
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

  my $dbhost = $target_db->dbc->host;
  my $dbport = $target_db->dbc->port;
  my $dbuser = $target_db->dbc->user;
  my $dbpass = $target_db->dbc->password;
  my $dbname = $target_db->dbc->dbname;

  for my $table ('seq_region','meta') {
    my $backup_file = catfile($path_to_files, $table.'.'.time().'.sql');
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
  my ($self,$target_dba,$meta_keys,$path_to_files) = @_;

  my $meta_adaptor = $target_dba->get_MetaContainerAdaptor;

  foreach my $meta_key (keys(%{$meta_keys})) {
    my $meta_values = ref($meta_keys->{$meta_key}) ? $meta_keys->{$meta_key} : [$meta_keys->{$meta_key}];
    foreach my $value (@$meta_values) {
      $meta_adaptor->store_key_value($meta_key, $value);
      say "Inserted into meta:\n".$meta_key." => ".$value;
    }
  }

  my $date = localtime->strftime('%Y-%m-Ensembl');
  $meta_adaptor->store_key_value('genebuild.start_date', $date);

  unless(-e $path_to_files."/assembly_report.txt") {
    $self->throw("Could not find the assembly_report.txt file. Path checked:\n".$path_to_files."/assembly_report.txt");
  }

  open(IN, $path_to_files.'/assembly_report.txt') || $self->throw("Could not open $path_to_files/assembly_report.txt");
  while (my $line = <IN>) {
    # This check seems odd, might mess up if run on plants
    if($line !~ /^#/ or $self->param('has_mitochondria')) {
      if ($line =~ /(NC_\S+)\s+non-nuclear/) {
        $self->param('mt_accession', $1);
      }
    } elsif($line =~ /^#\s*Date:\s*(\d+)-(\d+)-\d+/) {
      $meta_adaptor->store_key_value('assembly.date', sprintf("%d-%02d", $1, $2));
      say "Inserted into meta:\nassembly.date => ".sprintf("%d-%02d", $1, $2);
   } elsif ($line =~ /^#\s*Assembly level:\s*[Cc]hromosome/) {
      $self->param('chromosomes_present', 1);
    } elsif ($line =~ /##\s*GCF_\S+\s*non-nuclear/) {
      $self->param('has_mitochondria', 1);
    }
  }

  close(IN) || $self->throw("Could not close $path_to_files/assembly_report.txt");
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
  my ($self, $target_dba ,$path_to_files) = @_;

  my $sa = $target_dba->get_SliceAdaptor;
  my $insdc_id = $target_dba->get_DBEntryAdaptor->get_external_db_id('INSDC');
  my %coord_system = (
    chr2acc => $self->param('chromosome_level'),
    scaffold_localID2acc => $self->param('scaffold_level'),
    component_localID2acc => $self->param('contig_level')
  );

  my @files = ('component_localID2acc', 'scaffold_localID2acc');
  if($self->param('chromosomes_present')) {
    push(@files,'chr2acc');
  }

  my $sth = $target_dba->dbc->prepare('UPDATE seq_region set name=? where name=?');
  foreach my $filename (@files) {
    my $na_count = 0;
    my $line_count = 0;
    my $file = catfile($path_to_files, $filename);
    if (-e $file) {
      my $insert_count = 0;
      open(IN, $file) || $self->throw("Could not open $file");
      while (my $line = <IN>) {
        if ($line =~ /^#/) {
          next;
        }

        $line_count++;

        my ($synonym, $seq_region_name) = $line =~ /(\S+)\s+(\S+)/;
        if($synonym eq 'na') {
          $na_count++;
          $self->warning($seq_region_name." does not have a synonym, skipping...\n");
          next;
        }
        # Make sure we don't do a fuzzy search
        my $slice = $sa->fetch_by_region($coord_system{$filename}, $seq_region_name, undef, undef, undef, undef, 1);
        if ($slice) {
          if ($filename eq 'chr2acc') {
            # In this file the synonym is what we want as sequence name
            $slice->add_synonym($seq_region_name, $insdc_id);
            $sa->update($slice);
            $sth->bind_param(1,$synonym);
            $sth->bind_param(2,$seq_region_name);
            $sth->execute;
          }
          else {
            $slice->add_synonym($synonym);
            $sa->update($slice);
          }
          $insert_count++;
        }
        else {
          $self->warning("Could not find $seq_region_name in the database to add synonym $synonym");
        }
      }
      close(IN) || $self->throw("Could not close $file");
      if($insert_count == 0 && ($na_count != $line_count)) {
        $self->warning("The insert/update count after parsing ".$file." was 0, this is probably wrong. File used:\n".$file);
      }
      say "\nInserted into seq_region_synonym and updated seq_region based on ".$file.". Total inserts/updates: ".$insert_count;
      say "You will need to update the external_db_id for the synonyms of scaffold or contigs!\n";
    }
    else {
      $self->warning("Could not find $filename file. No synonyms loaded. Expected location:\n".$file)
        unless ($filename eq 'chr2acc' and $self->param('chromosomes_present'));
    }
  }
}

1;
