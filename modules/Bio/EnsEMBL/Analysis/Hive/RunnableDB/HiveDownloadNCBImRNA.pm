=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadNCBImRNA

=head1 SYNOPSIS


=head1 DESCRIPTION

Module to download sequences from NCBI using esearch and efetch

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadNCBImRNA;

use strict;
use warnings;

use LWP::UserAgent;
use File::Spec::Functions;
use File::Basename qw(dirname);
use File::Path qw(make_path);

use Bio::SeqIO;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 param_defaults

 Arg [1]    : None
 Description:
 Returntype :
 Exceptions :

=cut

sub param_defaults {
  my $self = shift;

  return {
    %{$self->SUPER::param_defaults()},
    base_url => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils',
    search_url => 'esearch.fcgi',
    fetch_url => 'efetch.fcgi',
    ncbidb => 'nucleotide',
    filetype => 'fasta',
    filemode => 'text',
    batch_size => 5000,
    http_proxy => undef,
    _input_id_name => 'query',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Search for data against NCBI using esearch based on the query parameter
 Returntype : None
 Exceptions : Throws if 'output_file' is not set

=cut

sub fetch_input {
  my $self = shift;

  $self->param_required('output_file');
  my $url = $self->param('base_url').'/'.$self->param('search_url').'?db='.$self->param('ncbidb').'&term='.$self->input_id.'&usehistory=y';

  if (!-d dirname($self->param('output_file')) ) {
    make_path(dirname($self->param('output_file')));
  }

  my $ua = LWP::UserAgent->new();
  $self->param('connector', $ua);
  if ($self->param_is_defined('http_proxy')) {
    $url =~ /^(\w+)/;
    $ua->proxy(['http', 'ftp', 'https'], $self->param('http_proxy'));
    $self->warning($1."\n");
  }
  else {
    $ua->env_proxy;
  }
  $self->warning($url);
  my $response = $ua->get($url);
  $self->throw('Could not retrieve data for '.$self->input_id."\n".$response->status_line) unless ($response->is_success);

  my $output = $response->decoded_content;
  $self->warning($output);
  $output =~ /<Count>(\d+)<\/Count>/m;
  if ($1 > 0) {
    $self->param('count', $1);
    $self->param('webenv', $1) if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
    $self->param('querykey', $1) if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
  }
  else {
    $self->input_job->autoflow(0);
    $self->complete_early('No sequences found for: '.$self->input_id);
  }
}


=head2 run

 Arg [1]    : None
 Description: Do nothing as storing the fetched data then write using write_output
              would use memory for nothing...
 Returntype : None
 Exceptions : None

=cut

sub run {
  return 1;
}


=head2 write_output

 Arg [1]    : None
 Description: Fetch the data from the search history using efetch and store the data in
              the file 'output_file'. It flows the filename of the written file on
              channel 2, 'iid'.
 Returntype : None
 Exceptions : Throws if it cannot write the file

=cut

sub write_output {
  my $self = shift;
  my $ua = $self->param('connector');
  my $fetch_url = $self->param('base_url').'/'.$self->param('fetch_url').'?db='.$self->param('ncbidb').'&WebEnv='.$self->param('webenv').'&query_key='.$self->param('querykey');
  my $datatype = '&retmax='.$self->param('batch_size').'&rettype='.$self->param('filetype').'&retmode='.$self->param('filemode');
  open(WH, '>'.$self->param('output_file')) || $self->throw('Could not open '.$self->param('output_file'));
  for (my $start = 0; $start < $self->param('count'); $start += $self->param('batch_size')) {
    my $response = $ua->get($fetch_url.'&retstart='.$start.$datatype);
    while($response->code == 502) {
      $response = $ua->get($fetch_url.'&retstart='.$start.$datatype);
    }
    if ($response->is_success) {
      print WH $response->decoded_content;
    }
    else {
      $self->throw('Could not retrieve data for '.$fetch_url.'&retstart='.$start.$datatype."\n".$response->status_line);
    }
  }
  close(WH) || $self->throw('Could not close '.$self->param('output_file'));

  sleep(45);
  # check the number of sequences that you have saved: 
  my $count = 0;
  my $format = $self->param('filetype');
  $format = 'genbank' if ($format eq 'gb');
  my $parser = Bio::SeqIO->new( -format => $format, -file => $self->param('output_file'));
  while ($parser->next_seq) {
    ++$count;
  }
  if ($count!=$self->param('count')) {
    $self->throw('Sequences in output file are ' . $count . '- while you should have: ' . $self->param('count'));
  }
  $self->dataflow_output_id({iid => $self->param('output_file')}, $self->param('_branch_to_flow_to'));
}

1;
