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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFTPURLs

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFTPURLs;

use strict;
use warnings;

use Net::FTP;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');


=head2 param_defaults

 Arg [1]    : None
 Description: Set the default column_names to 'url'
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    column_names => ['url'],
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Use a base url and a list of file names, can have '*'
              to generate a list of url to download from a FTP
 Returntype : None
 Exceptions : Throws if 'base_url' or 'file_list' is not set
              Throws if a file does not exist

=cut

sub fetch_input {
  my ($self) = @_;

  my $base_url = $self->param_required('base_url');
  my ($type, $host, $path) = $base_url =~ /(\w+:\/\/)([^\/]+)\/(.*)/;
  $path =~ s/\/$//;
  my $client = Net::FTP->new($host);
  $client->login;
  $client->cwd($path) || $self->throw('could not go in '.$path);
  my @iids;
  foreach my $filename (@{$self->param_required('file_list')}) {
    my $url = $base_url;
    if ($filename =~ /(.*)\/([^\/]+)$/) {
      my $fpath = $1;
      $filename = $2;
      $client->cwd($fpath) || $self->throw('Could not open '.$fpath);
      $fpath =~ s/^\///;
      $url .= '/'.$fpath;
    }
    if ($filename =~ s/\*/.*/g) {
      foreach my $ftp_file ($client->ls) {
        if ($ftp_file =~ /$filename/) {
          push(@iids, $url.'/'.$ftp_file);
        }
      }
    }
    else {
      $client->size($filename) || $self->throw("File $filename does not exist in ".$client->pwd);
      push(@iids, $url.'/'.$filename);
    }
    $client->cwd;
    $client->cwd($path);
  }
  $client->quit;

  $self->param('inputlist', \@iids);
}

1;
