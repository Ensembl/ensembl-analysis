=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

use Bio::EnsEMBL::Analysis::Runnable::Aspera;
use File::Fetch;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    use_perl => 1,
    aspera_user => 'era-fasp',
    aspera_host => 'ftp.sra.ebi.ac.uk',
    uncompress => 1,
  }
}


sub fetch_input {
  my ($self) = @_;

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
    $client = File::Fetch->new(uri => $url) || $self->throw('Could not create a fetcher for '.$url);
    $self->param('options', [to => $output_dir]);
  }
  $self->param('client', $client);
}

sub run {
  my ($self) = @_;

  my $client = $self->param('client');
  my $file = $client->fetch(($self->param_is_defined('options') ? @{$self->param('options')}: undef));
  $self->check_file($file);
  $file = $self->uncompress($file) if ($self->param('uncompress'));
  $self->output([$file]);
}


sub uncompress {
  my ($self, $file) = @_;

  my ($output) = $file =~ /^(\S+)\.\w+$/;
  anyuncompress $file => $output
    or $self->throw("anyuncompress failed: $AnyUncompressError");
  unlink $file;
  return $output;
}

sub check_file {
  my ($self, $file) = @_;

  $self->throw("The file $file does not exist") unless (-e $file);
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

sub write_output {
  my ($self) = @_;

  my @output_ids;
  foreach my $file (@{$self->output}) {
    push(@output_ids, {filename => $file});
  }
  $self->dataflow_output_id(\@output_ids, $self->param('_branch_to_flow_to'));
}

1;
