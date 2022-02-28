=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadRefSeqSynonyms

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadRefSeqSynonyms;

use strict;
use warnings;

use Net::FTP;
use File::Fetch;
use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters for the module, these should not change
               ncbi_ftp_host => 'ftp.ncbi.nlm.nih.gov',
               ncbi_ftp_user => 'anonymous',
               ncbi_ftp_passwd => undef,
               external_db => 'RefSeq_genomic',
               min_synonym_load_threshold => 1, # 1 means all, 0.X means a fraction of them
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    ncbi_ftp_host => 'ftp.ncbi.nlm.nih.gov',
    ncbi_ftp_user => 'anonymous',
    ncbi_ftp_passwd => undef,
    external_db => 'RefSeq_genomic',
    min_synonym_load_threshold => 1, # 1 means all, 0.X means a fraction of them
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Use the assembly_refseq_accession parameters to retrieve all
              possible synonyms between the GenBank and the RefSeq assembly.
              It uses the assembly_report.txt, chr2acc (if it exists),
              scaffold_localID2acc and component_localID2acc. It looks in
              the "assembly root directory", %s_assembly_structure,
              Primary_Assembly and assembled_chromosomes (if it exists).
              It store the data in a hashref in '' where the key is the
              coord_system and the data is an arrayref of [seq_region_name, synonym]
 Returntype : None
 Exceptions : Throws if 'assembly_name' is not set or found in the meta table
              Throws if 'external_db_id' is not set or found in the database
              Throws if 'assembly_refseq_accession' is not set
              Throws if any problem is encounter while downloading data from
                the RefSeq FTP

=cut

sub fetch_input {
  my ($self) = @_;

  $self->create_analysis;
  my $db = $self->get_database_by_name('target_db');
  $self->hrdb_set_con($db, 'target_db');
  my $assembly_name;
  if ($self->param_is_defined('assembly_name')) {
    $assembly_name = $self->param('assembly_name');
  }
  else {
    $assembly_name = $db->get_MetaContainer->single_value_by_key('assembly.name');
  }
  $self->throw('Could not get assembly_name') unless ($assembly_name);
  my $primary_assembly_cs = $db->get_CoordSystemAdaptor->fetch_by_name('primary_assembly', $assembly_name);
  if ($primary_assembly_cs) {
# if the coord_system is 'primary_assembly' the refseq synonyms should have been imported during 'process_assembly_info', however, there are cases where the refseq synonyms have been updated but this is not reflected in the NCBI assembly report (i.e. can be found in GCF but not in GCA): in this case download the GCF assembly report and load synonyms from there 
    my $assembly_refseq_accession = $self->param_required('assembly_refseq_accession');
    my $refseq_db_id = $db->get_DBEntryAdaptor->get_external_db_id($self->param('external_db'), undef, 1);
    if ($refseq_db_id) {
      $self->param('external_db_id', $refseq_db_id);
    }
    else {
      $self->throw('Could not fetch the dbID for '.$self->param('external_db'));
    }
    my $refseq_ftp_url = $self->param('url');
    my $client = File::Fetch->new(uri => $refseq_ftp_url) || $self->throw('Could not create a fetcher for '.$refseq_ftp_url);
    $self->param('options', [to => $self->param('output_dir')]);

    my $file = $client->fetch(($self->param_is_defined('options') ? @{$self->param('options')}: undef));

    my $data_counter = 0;
    my %data;

    if ($file) {
      open(FH, "$file") || $self->throw("Could not open $file");
      while(my $line = <FH>) {
        next if ($line =~ /^#|^\s*$/);
        $line =~ s/\R$//;
        my @line = split("\t", $line);
        my $cs = 'primary_assembly';
	if ($line[6] ne "na"){
	  push(@{$data{$cs}}, [$line[4], $line[6]]);
	  ++$data_counter;
	}
      }
      close(FH) || $self->throw("Could not close $file");
      unlink($file);
    }
    $self->param('synonym_count', $data_counter);
    $self->param('synonyms', \%data);

  }
  else {
    my $refseq_db_id = $db->get_DBEntryAdaptor->get_external_db_id($self->param('external_db'), undef, 1);
    if ($refseq_db_id) {
      $self->param('external_db_id', $refseq_db_id);
    }
    else {
      $self->throw('Could not fetch the dbID for '.$self->param('external_db'));
    }
    my $ncbi_assembly_accession = $self->param_required('assembly_refseq_accession');
    my $ftpclient = Net::FTP->new($self->param('ncbi_ftp_host')) or $self->throw("Can't open ".$self->param('ncbi_ftp_host'));
    $ftpclient->login($self->param('ncbi_ftp_user'), $self->param('ncbi_ftp_passwd')) or $self->throw("Can't log ".$self->param('ncbi_ftp_user').' in');

    my $counter = 0;
    my $data_counter = 0;
    my %data;
    my %files = (
      chromosome => 'chr2acc',
      scaffold => 'scaffold_localID2acc',
      contig => 'component_localID2acc',
    );
    my %found_files = (
      chromosome => 0,
      scaffold => 0,
      contig => 0,
    );
    foreach my $dir ('genomes',
                     'all',
                     substr($ncbi_assembly_accession, 0, 3),
                     substr($ncbi_assembly_accession, 4, 3),
                     substr($ncbi_assembly_accession, 7, 3),
                     substr($ncbi_assembly_accession, 10, 3),
                     sprintf("%s_%s", $ncbi_assembly_accession, $assembly_name)) {
      $ftpclient->cwd($dir) || $self->throw("Could not got into $dir");
    }
    my $fh = $ftpclient->get(sprintf("%s_%s_assembly_report.txt", $ncbi_assembly_accession, $assembly_name));
    if ($fh) {
      open(FH, "$fh") || $self->throw("Could not open $fh");
      while(my $line = <FH>) {
        next if ($line =~ /^#|^\s*$/);
        $line =~ s/\R$//;
        my @line = split("\t", $line);
        my $cs = 'scaffold';
        $cs = 'chromosome' if ($line[1] eq 'assembled-molecule');
        push(@{$data{$cs}}, [$line[4], $line[6]]);
      }
      close(FH) || $self->throw("Could not close $fh");
      ++$counter;
      unlink($fh);
    }
    foreach my $dir ('.',
                     sprintf("%s_%s_assembly_structure", $ncbi_assembly_accession, $assembly_name),
                     'Primary_Assembly',
                     'assembled_chromosomes') {
      if ($ftpclient->cwd($dir)) {
        foreach my $file (keys %files) {
          my $fh = $ftpclient->get($files{$file});
          if ($fh) {
            open(FH, "$fh") || $self->throw("Could not open $fh");
            while(my $line = <FH>) {
              next if ($line =~ /^#|^\s*$/);
              $line =~ s/\R$//;
              my @synonyms = split("\t", $line);
              if ($synonyms[0] and $synonyms[0] ne 'na' and $synonyms[1] and $synonyms[1] ne 'na') {
                push(@{$data{$file}}, \@synonyms);
                ++$data_counter;
              }
              else {
                $self->say_with_header('No correct synonyms '.$line);
              }
            }
            close(FH) || $self->throw("Could not close $fh");
            ++$counter;
            unlink($fh);
            ++$found_files{$file};
          }
        }
      }
    }
    my $found_file_count = 0;
    foreach my $file (keys %found_files) {
      $self->throw("You have a problem with $file as it has been found ${$files{$file}} times")
        if ($found_files{$file} > 1);
      $found_file_count += $found_files{$file};
    }
    $self->throw("No files were found") unless ($found_file_count);
    $self->throw("No accession was found") unless ($counter);
    $self->param('synonym_count', $data_counter);
    $self->param('synonyms', \%data);
  }
}


=head2 run

 Arg [1]    : None
 Description: Does nothing
 Returntype : None
 Exceptions : None

=cut

sub run {
# There is nothing to do know, let's write the data
  return 1;
}


=head2 write_output

 Arg [1]    : None
 Description: Write the synonyms for all the seq_regions
 Returntype : None
 Exceptions : Throws if the number of synonyms is lower than 2/3 of the expectation

=cut

sub write_output {
  my ($self) = @_;

  # The real way to do it would be to only add synonyms in run then only update the slice in write_output
  # It may use too much memory in some cases
  my $external_db_id = $self->param('external_db_id');
  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  my $data = $self->param('synonyms');
  my $counter = 0;
  foreach my $coord_system (keys %$data) {
    foreach my $synonyms (@{$data->{$coord_system}}) {
      my $slice = $slice_adaptor->fetch_by_region($coord_system, $synonyms->[1]);
      if (!$slice) {
        $slice = $slice_adaptor->fetch_by_region($coord_system, $synonyms->[0]);
        if ($slice) {
          $slice->add_synonym($synonyms->[1], $external_db_id);
          $slice->adaptor->update($slice);
          $self->say_with_header('Added '.$synonyms->[1].' to '.$synonyms->[0]);
        }
        else {
          $self->warning('The slice '.$synonyms->[0].' with synonym '.$synonyms->[1].' does not exist on '.$coord_system);
        }
      }
# We also count it if it already exists
      ++$counter;
    }
  }
  if ($counter ne $self->param('synonym_count')) {
    if (($counter/$self->param('synonym_count')) <= $self->param('min_synonym_load_threshold')) {
      $self->throw("You are missing some synonyms: $counter instead of ".$self->param('synonym_count'));
    }
    else {
      $self->warning("There is a discrepancy between synonym counts: $counter instead of ".$self->param('synonym_count'));
    }
  }
}

1;
