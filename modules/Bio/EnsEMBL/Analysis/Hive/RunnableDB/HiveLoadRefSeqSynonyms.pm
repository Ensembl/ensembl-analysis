=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
use File::Temp;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    ncbi_ftp_host => 'ftp.ncbi.nlm.nih.gov',
    ncbi_ftp_user => 'anonymous',
    ncbi_ftp_passwd => undef,
    external_db => 'RefSeq_genomic',
  }
}
sub fetch_input {
  my ($self) = @_;

  $self->create_analysis;
  $self->hrdb_set_con($self->get_database_by_name('target_db'), 'target_db');
  my $assembly_name;
  if ($self->param_is_defined('assembly_name')) {
    $assembly_name = $self->param('assembly_name');
  }
  else {
    $assembly_name = $self->hrdb_get_con('target_db')->get_MetaContainer->single_value_by_key('assembly.name');
  }
  $self->throw('Could not get assembly_name') unless ($assembly_name);
  my $refseq_db_id = $self->hrdb_get_con('target_db')->get_DBEntryAdaptor->get_external_db_id($self->param('external_db'));
  if ($refseq_db_id) {
    $self->param('external_db_id', $refseq_db_id);
  }
  else {
    $self->throw('Could not fetch the dbID for '.$self->param('external_db'));
  }
  my $ncbi_assembly_accession = $self->param_required('assembly_refseq_accession');
  my $uri = sprintf("genomes/all/%s/%03d/%03d/%03d/%s_%s/%s_%s_assembly_structure/Primary_Assembly", substr($ncbi_assembly_accession, 0, 3), substr($ncbi_assembly_accession, 4, 3), substr($ncbi_assembly_accession, 7, 3), substr($ncbi_assembly_accession, 10, 3), $ncbi_assembly_accession, $assembly_name, $ncbi_assembly_accession, $assembly_name);
  my $ftpclient = Net::FTP->new($self->param('ncbi_ftp_host')) or $self->throw("Can't open ".$self->param('ncbi_ftp_host'));
  $ftpclient->login($self->param('ncbi_ftp_user'), $self->param('ncbi_ftp_passwd')) or $self->throw("Can't log ".$self->param('ncbi_ftp_user').' in');
  
  my $counter = 0;
  my $data_counter = 0;
# The order is important here because fetch_by_region returns the highest rank so when we have duplicated name for any reasons we will get the correct ones
  my %data;
  my %files = (
    chromosome => 'assembled_chromosomes/chr2acc',
    scaffold => 'scaffold_localID2acc',
    contig => 'component_localID2acc',
  );
  foreach my $file (keys %files) {
#    my $fh = File::Temp->new();
#    if ($ftpclient->get("$uri/".$files{$file}, $fh)) {
    my $fh = $ftpclient->get("$uri/".$files{$file});
    if ($fh) {
      open(FH, "$fh") || $self->throw("Could not open $fh");
      while(<FH>) {
        next if (/^#|^\s*$/);
        push(@{$data{$file}}, [split("\t", $_)]);
        ++$data_counter;
      }
      close(FH) || $self->throw("Could not close $fh");
      ++$counter;
      unlink($fh);
    }
  }
  $self->throw("No accession was found") unless ($counter);
  $self->param('synonym_count', $data_counter);
  $self->param('synonyms', \%data);
}

sub run {
# There is nothing to do know, let's write the data
  return 1;
}

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
      if (!$synonyms->[0] or !$synonyms->[1]) {
        print STDERR $synonyms->[0], ' ', $synonyms->[1], "\n";
      }
      my $slice = $slice_adaptor->fetch_by_region($coord_system, $synonyms->[1]);
      if ($slice) {
# It already exists so we can count it 
        ++$counter;
      }
      else {
        $slice = $slice_adaptor->fetch_by_region($coord_system, $synonyms->[0]);
        $self->throw('The slice '.$synonyms->[0].' with synonym '.$synonyms->[1].' does not exist')
          unless ($slice);
        $slice->add_synonym($synonyms->[1], $external_db_id);
        $slice->adaptor->update($slice);
        ++$counter;
      }
    }
  }
  $self->throw("You are missing some synonyms: $counter instead of ".$self->param('synonym_count'))
    unless ($counter eq $self->param('synonym_count'));
}
1;

