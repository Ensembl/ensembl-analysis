#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions qw(catfile file_name_is_absolute);
use File::Path qw(make_path);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use LWP::UserAgent;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    base_url => 'https://rest.uniprot.org/uniprotkb/stream?query=',
    format => 'fasta',
    max_attempts => 4,
    retry_delay => 5,
    timeout => 300,
    min_total_sequences => 1,
    min_expected_fraction => 0.95,
  }
}

sub fetch_input {
  my $self = shift;

  if ($self->param_is_defined('query_url')) {
    my $url = $self->param('query_url');
    if ($url !~ /http|ftp/) {
      $url = $self->param('base_url').$url;
    }
    my $filename = $self->param_required('file_name');
    if (!file_name_is_absolute($filename)) {
      $filename = catfile($self->param_required('dest_dir'), $filename);
    }
    $self->param('query_url', [{url => $url, file_name => $filename}]);
  }
  else {
    if($self->param('multi_query_download')) {
      my @urls;
      foreach my $hash (values %{$self->param('multi_query_download')}) {
        push(@urls, $self->build_query($hash));
      }
      $self->param('query_url', \@urls);
    }
    else {
      my %hash = (
        dest_dir => $self->param_required('dest_dir'),
        file_name => $self->param_required('file_name'),
        pe_level => $self->param_required('pe_level'),
        format => $self->param('format'),
      );
      $hash{compress} = $self->param('compress') if ($self->param_is_defined('compress'));
      $hash{taxon_group} = $self->param('taxon_group') if ($self->param_is_defined('taxon_group'));
      $hash{taxon_id} = $self->param('taxon_id') if ($self->param_is_defined('taxon_id'));
      $hash{exclude_id} = $self->param('exclude_id') if ($self->param_is_defined('exclude_id'));
      $hash{exclude_group} = $self->param('exclude_group') if ($self->param_is_defined('exclude_group'));
      $hash{compress} = $self->param('compress') if ($self->param_is_defined('compress'));
      $hash{mito} = $self->param('mito') if ($self->param_is_defined('mito'));
      $hash{fragment} = $self->param('fragment') if ($self->param_is_defined('fragment'));
      $hash{isoforms} = $self->param('isoforms') if ($self->param_is_defined('isoforms'));
      $self->param('query_url', [$self->build_query(\%hash)]);
    }
  }
  return 1;
}

sub run {
  my $self = shift;

  my @iids;
  foreach my $query (@{$self->param_required('query_url')}) {
    my $query_url = $query->{url};
    my $filename = $query->{file_name};
    say "Downloading:\n$query_url\n";

    my $agent = LWP::UserAgent->new(
      agent   => 'libwww-perl',
      timeout => $self->param('timeout'),
    );
    $agent->env_proxy;

    my $expected_count = $self->_get_expected_count($agent, $query_url);
    say "UniProt preflight count: $expected_count sequences" if defined $expected_count;

    my ($response, $last_error);
    ATTEMPT: for my $attempt (1 .. $self->param('max_attempts')) {
      $response = eval { $agent->get($query_url) };
      $last_error = $@;
      last ATTEMPT if ($response && $response->is_success() &&
                       ($filename =~ /\.gz$/ || $self->_fasta_sequence_count($response->content())));

      if ($response && $response->status_line() =~ /^404\s/) {
        last ATTEMPT;
      }

      my $retry_after = $response ? $response->header('Retry-After') : undef;
      my $status = $response ? $response->status_line() : ($last_error || 'no response');
      if ($attempt < $self->param('max_attempts')) {
        my $delay = $retry_after && $retry_after =~ /^\d+$/
          ? $retry_after
          : $self->param('retry_delay') * (2 ** ($attempt - 1));
        say STDERR "Download attempt $attempt failed ($status); retrying in ${delay}s";
        sleep $delay;
      }
    }

    if ($response && $response->is_success()) {
      my $content = $response->content();
      my $sequence_count = $filename =~ /\.gz$/ ? 0 : $self->_fasta_sequence_count($content);
      $self->throw("Downloaded response for '$filename' is not a non-empty FASTA file")
        unless $filename =~ /\.gz$/ || $sequence_count;

      my $tmp_filename = $filename . '.tmp.' . $$;
      if ($filename =~ /\.gz$/) {
        $tmp_filename .= '.gz';
      }
      open(my $fh, '>', $tmp_filename) || $self->throw("Could not open file $tmp_filename\n");
      binmode $fh;
      print $fh $content;
      close($fh) || $self->throw("Could not close the file '$tmp_filename'");

      if ($tmp_filename =~ s/\.gz$//) {
        gunzip "$tmp_filename.gz" => $tmp_filename
          or $self->throw("gunzip failed for '$filename': $GunzipError");
        unlink "$tmp_filename.gz"
          or $self->throw("Could not remove temporary compressed file '$tmp_filename.gz': $!");
        $filename =~ s/\.gz$//;
        $sequence_count = $self->_fasta_sequence_count_file($tmp_filename);
        $self->throw("Downloaded response for '$filename' is not a non-empty FASTA file")
          unless $sequence_count;
      }
      rename($tmp_filename, $filename)
        or $self->throw("Could not atomically move '$tmp_filename' to '$filename': $!");
      push(@iids, $filename);
      $self->param('_downloaded_sequence_count', ($self->param('_downloaded_sequence_count') || 0) + $sequence_count);
      if (defined $expected_count && $sequence_count < $expected_count * $self->param('min_expected_fraction')) {
        $self->throw("Downloaded $sequence_count protein sequences for '$filename'; UniProt preflight reported $expected_count (minimum allowed is ".int($expected_count * $self->param('min_expected_fraction')).")");
      }
    } elsif ($response && $response->status_line() =~ /404 Not Found/) {
      $self->warning('Failed, got '.$response->status_line().' for '.$response->request()->uri()." . It's likely that this sequence has been deleted. Fine. \n");
    } else {
      my $status = $response ? $response->status_line() : ($last_error || 'no response');
      $self->throw("Failed to download '$query_url' after ".$self->param('max_attempts')." attempts: $status\n");
    }
  }
  my $min_total_sequences = $self->param('min_total_sequences');
  if (defined $min_total_sequences && ($self->param('_downloaded_sequence_count') || 0) < $min_total_sequences) {
    $self->throw("Downloaded only ".($self->param('_downloaded_sequence_count') || 0)." protein sequences; expected at least $min_total_sequences");
  }
  $self->output(\@iids);
  say "Finished downloading UniProt files";
  return 1;
}

sub _get_expected_count {
  my ($self, $agent, $stream_url) = @_;

  # The stream endpoint returns the FASTA body, whereas search with size=1
  # returns the total matching-record count in X-Total-Results.
  my $count_url = $stream_url;
  $count_url =~ s{/uniprotkb/stream\?}{/uniprotkb/search?};
  $count_url =~ s/&compress=[^&]+//;
  $count_url =~ s/&format=[^&]+/&format=tsv&size=1/;
  $count_url .= '&format=tsv&size=1' unless $count_url =~ /[&?]size=/;

  my ($response, $last_error);
  for my $attempt (1 .. $self->param('max_attempts')) {
    $response = eval { $agent->get($count_url) };
    $last_error = $@;
    if ($response && $response->is_success()) {
      my $count = $response->header('X-Total-Results');
      return $count if defined $count && $count =~ /^\d+$/;
    }
    last if $response && $response->status_line() =~ /^404\s/;
    if ($attempt < $self->param('max_attempts')) {
      my $delay = $response && $response->header('Retry-After');
      $delay = $self->param('retry_delay') * (2 ** ($attempt - 1))
        unless defined $delay && $delay =~ /^\d+$/;
      say STDERR "UniProt preflight attempt $attempt failed; retrying in ${delay}s";
      sleep $delay;
    }
  }
  my $status = $response ? $response->status_line() : ($last_error || 'no response');
  $self->throw("UniProt preflight count failed for '$stream_url' after ".$self->param('max_attempts')." attempts: $status");
}

sub _fasta_sequence_count {
  my ($self, $content) = @_;
  return 0 unless defined $content && $content =~ /^\s*>/;
  return scalar grep { /^>/ } split(/\n/, $content);
}

sub _fasta_sequence_count_file {
  my ($self, $filename) = @_;
  open(my $fh, '<', $filename) || $self->throw("Could not open downloaded FASTA '$filename'");
  local $/;
  my $content = <$fh>;
  close($fh) || $self->throw("Could not close downloaded FASTA '$filename'");
  return $self->_fasta_sequence_count($content);
}

sub write_output {
  my $self = shift;

  my @iids = map {{iid => $_, iid_type => 'filename'}} @{$self->output};
  $self->dataflow_output_id(\@iids, $self->param('_branch_to_flow_to'));
}


sub build_query {
  my ($self,$query_params) = @_;
  my $taxon_id = $query_params->{'taxon_id'};
  my $taxon_group = $query_params->{'taxon_group'};
  my $exclude_id = $query_params->{'exclude_id'};
  my $exclude_group = $query_params->{'exclude_group'};
  my $dest_dir = $query_params->{'dest_dir'};
  my $file_name = $query_params->{'file_name'};
  my $pe_level = $query_params->{'pe_level'};
  my $pe_string = "(";
  my $taxonomy_string = "";
  my $exclude_string = "";
  my $compress = "yes";
  my $fragment_string = "+AND+fragment%3Afalse";
  my $mito = "+NOT+organelle%3Amitochondrion";
  my $format = $self->param('format');

  if(exists($query_params->{'compress'})) {
    if($query_params->{'compress'} eq '0') {
      $compress = "no";
    }
  }

  if(exists($query_params->{'mito'})) {
    if($query_params->{'mito'}) {
      $mito = undef;
    }
  }

  if(exists($query_params->{'fragment'})) {
    if($query_params->{'fragment'}) {
      $fragment_string = "+AND+fragment%3Atrue";
    }
  }

  if(exists($query_params->{'format'})) {
    $format = $query_params->{'format'};
  }

  my $full_query = $self->param('base_url');

  # Must have file_name, pe_level, dest_dir and either taxon_id or taxonomy
  unless($file_name && $dest_dir && ($taxon_id || $taxon_group) && $pe_level) {
    $self->throw("Must define the following keys:\nfile_name\ntaxon_id or taxonomy\ndest_dir\npe_level");
  }

  my @pe_array = @{$pe_level};
  unless(scalar(@pe_array)) {
    $self->throw("Not PE levels found in value of pe_levels key. Format should be an array ref: ['1,2']");
  }

  foreach my $pe_level (@pe_array) {
    unless($pe_level =~ /\d+/) {
     $self->throw("Could not parse a PE level from the following: ".$pe_level);
    }

    my $parsed_pe_level = $&;
    unless($parsed_pe_level >= 1 && $parsed_pe_level <= 5) {
     $self->throw("Parsed PE level is outside the normal range of 1-5: ".$parsed_pe_level);
   }
   $pe_string .= 'existence%3A'.$pe_level.'+OR+';
  }

  $pe_string =~ s/\+OR\+$/\)/;

  # NOTE this bit of the code with taxonomy and exclude is shit and needs to be upgraded
  if($taxon_id) {
    $taxonomy_string = '+AND+taxonomy_id%3A'.$taxon_id;
  } elsif($taxon_group) {
    $taxonomy_string = '+AND+taxonomy_id%3A'.$taxon_group;
  }

#+NOT+taxonomy%3A%22
  if($exclude_id) {
    my @exclusion_array = @{$exclude_id};
    foreach my $id_to_exclude (@exclusion_array) {
      $exclude_string .= '+NOT+taxonomy_id%3A'.$id_to_exclude;
    }
  }

  $compress .= '&include=yes' if (exists $query_params->{isoforms} and $query_params->{isoforms});
  $full_query .= $pe_string.$taxonomy_string.$exclude_string.$fragment_string.$mito.'&compress='.$compress.'&format='.$format;
  if (!-d $query_params->{dest_dir}) {
    make_path($query_params->{dest_dir});
  }
  my $filename = catfile($query_params->{dest_dir}, $query_params->{file_name});
  $filename .= '.gz' if ($compress eq 'yes');;

  return {url => $full_query, file_name => $filename};
}

1;
