=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData;

use strict;
use warnings;

use Digest::MD5;
use File::Spec::Functions qw(splitpath file_name_is_absolute catfile);
use File::Path qw(make_path);
use File::Basename qw(basename);
use File::Fetch;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;

use Bio::EnsEMBL::Analysis::Runnable::Aspera;
use Bio::EnsEMBL::Analysis::Runnable::Samtools;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
               use_perl => 1,
               aspera_user => 'era-fasp',
               aspera_host => 'ftp.sra.ebi.ac.uk',
               uncompress => 1,
               create_faidx => 0,
               samtools => 'samtools',
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    use_perl => 1,
    aspera_user => 'era-fasp',
    aspera_host => 'ftp.sra.ebi.ac.uk',
    uncompress => 1,
    create_faidx => 0,
    samtools => 'samtools',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Check the different parameters and check the url. It will select the client to use
              based on the 'download_method' parameter.
 Returntype : None
 Exceptions : Throws if using Perl 5.24
              Throws if 'download_method' is not set
              Throws if 'output_dir' is not set
              Throws if 'url' is not set

=cut

sub fetch_input {
  my ($self) = @_;

  if ($] =~ '^5.024') {
    $self->throw("Perl 5.24 doesn't work with this module. If you manage to make it work, please submit a pull-request");
  }
  my $download_method = $self->param_required('download_method');
  my $output_dir = $self->param_required('output_dir');
  if (!-d $output_dir) {
    make_path($output_dir);
  }

  my $client;
  my $url = $self->param_required('url');
  if ($download_method eq 'aspera') {
    $client = Bio::EnsEMBL::Analysis::Runnable::Aspera->new();
    $client->options($self->param('commandline_params'));
    if (file_name_is_absolute($url)) {
      $url = $self->param('aspera_user').'@'.$self->param('aspera_host').':'.$url;
    }
    elsif ($url =~ /^[^@]+\//) {
      $url = $self->param('aspera_user').'@'.$url;
      $url =~ s/\//:\// if ($url !~ /:\//);
    }
    $client->source($url);
    $client->target($output_dir);
  }
  else {
    $url = "$download_method://$url" unless ($url =~ '^\w+://');
    $client = File::Fetch->new(uri => $url) || $self->throw('Could not create a fetcher for '.$url);
    $self->param('options', [to => $output_dir]);
  }
  $self->param('client', $client);
}


=head2 run

 Arg [1]    : None
 Description: Retrieve the file passed in parameters, check the integrity of the file when a md5sum
              has been provided. Decompress the file is 'uncompress' is set to 1.
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $client = $self->param('client');
  my $file;
  if ($self->param_is_defined('md5sum') and -e catfile($self->param_required('output_dir'), basename($self->param_required('url')))) {
    $file = catfile($self->param_required('output_dir'), basename($self->param_required('url')));
  }
  else {
    $file = $client->fetch(($self->param_is_defined('options') ? @{$self->param('options')}: undef));
  }
  if ($file) {
    $self->check_file($file);
    $file = $self->uncompress($file) if ($self->param('uncompress'));
    if ($self->param('create_faidx')) {
      my $samtools = Bio::EnsEMBL::Analysis::Runnable::Samtools->new(
                     -program => $self->param('samtools'),
                     );
      $samtools->index_genome($file);
    }
    $self->output([$file]);
  }
  else {
    $self->throw($client->error());
  }
}


=head2 uncompress

 Arg [1]    : String, path to file
 Description: Decompress the file
 Returntype : String, path to decompressed file
 Exceptions : Throws if the decrompression failed

=cut

sub uncompress {
  my ($self, $file) = @_;

  my ($output) = $file =~ /^(\S+)\.\w+$/;
  anyuncompress $file => $output
    or $self->throw("anyuncompress failed: $AnyUncompressError");
  unlink $file;
  return $output;
}


=head2 check_file

 Arg [1]    : String, path to file
 Description: If 'md5sum' is defined, check the md5 sum of the downloaded file
 Returntype : None
 Exceptions : Throws if the file doesn't exist
              Throws if the md5 sum differs for the expected checksum

=cut

sub check_file {
  my ($self, $file) = @_;

  $self->throw('The file "'.($file ||  '').'" does not exist') unless ($file and -e $file);
  if ($self->param_is_defined('md5sum')) {
    my $digest = Digest::MD5->new();
    open(my $fh, $file) || $self->throw('Could not open '.$file);
    binmode $fh;
    my $md5sum = $digest->addfile($fh)->hexdigest;
    close($fh) || $self->throw('Could not close '.$file);
    $self->throw('MD5 not as expected '.$self->param('md5sum').' but '.$md5sum)
      unless ($md5sum eq $self->param('md5sum'));
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Provide the name of the file downloaded on branch '_branch_to_flow_to'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my @output_ids;
  foreach my $file (@{$self->output}) {
    push(@output_ids, {filename => $file});
  }
  $self->dataflow_output_id(\@output_ids, $self->param('_branch_to_flow_to'));
}


=head2 pre_cleanup

 Arg [1]    : None
 Description: Remove the current file if it exists to avoid problems
 Returntype : None
 Exceptions : None

=cut

sub pre_cleanup {
  my ($self) = @_;

  my $filepath = catfile($self->param_required('output_dir'), basename($self->param_required('url')));
  if (-e $filepath) {
    unlink $filepath;
  }
}

1;
