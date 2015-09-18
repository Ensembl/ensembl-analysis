#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my $self = shift;

  # This is a bit lazy, originally this was designed to work with a single query
  # but now I also want it to work with multiple files. I'm going to put that into
  # the if and the else will retain the original single file logic
  if($self->param('multi_query_download')) {
    $self->multi_query_download();
  }

  else {
    my $query_url = $self->build_query();
    say "Downloading:\n".$query_url."\n";

    my $query_exit_code;
    $query_exit_code = system($query_url);
    unless($query_exit_code == 0) {
      throw("The wget query ended in an non-zero exit code:\n".$query_exit_code);
    }

    if($query_url =~ /\.gz$/) {
      my $file_path = $self->param('dest_dir')."/".$self->param('file_name');
      my $gunzip_command = "gunzip ".$file_path;
      my $gunzip_exit_code;
      $gunzip_exit_code = system($gunzip_command);
      unless($gunzip_exit_code == 0) {
        throw("gunzip on file ended in an non-zero exit code:\n".$gunzip_exit_code);
      }
    }
  }
  say "Finished downloading UniProt files";
  return 1;
}

sub write_output {
  my $self = shift;

  if($self->param('multi_query_download')) {
    my $file_hash = $self->output_hash;
    my $output_hash = {};
    $output_hash->{'iid'} = $file_hash;
    $self->dataflow_output_id($output_hash,1);
  } else {
    my $file_path = $self->param('dest_dir')."/".$self->param('file_name');
    my $output_hash = {};
    $output_hash->{'iid'} = $file_path;
    $self->dataflow_output_id($output_hash,1);
  }
  return 1;
}

sub output_hash {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_output_hash',$val);
  }

  return($self->param('_output_hash'));

}

sub multi_query_download {
  my ($self) = @_;

  my $output_hash = {};
  my $multi_query_hash = $self->param('multi_query_download');
  unless(ref($multi_query_hash) eq 'HASH') {
    $self->throw("You're trying to download multiple queries, but haven't passed in a hash ref. You ".
                 "need to pass in hash of hashes, where each internal hash has the relevant query params");
  }

  foreach my $query_key (keys(%$multi_query_hash)) {
    say "Building query based on params for: ".$query_key;
    my $query_params = $multi_query_hash->{$query_key};

    my $query_url = $self->build_query($query_params);
    say "Downloading:\n".$query_url."\n";

    my $query_exit_code;
    $query_exit_code = system($query_url);
    unless($query_exit_code == 0) {
      throw("The wget query ended in an non-zero exit code:\n".$query_exit_code);
    }

    my $file_path = $query_params->{'dest_dir'}."/".$query_params->{'file_name'};

    if($query_url =~ /\.gz$/) {
      my $gunzip_command = "gunzip ".$file_path.".gz";
      my $gunzip_exit_code;
      $gunzip_exit_code = system($gunzip_command);
      unless($gunzip_exit_code == 0) {
        $self->throw("gunzip on file ended in an non-zero exit code:\n".$gunzip_exit_code);
      }
    }

    $output_hash->{$query_key} = $file_path;
  }

  $self->output_hash($output_hash);

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
  my $fragment_string = "+AND+fragment:no";
  my $mito = "+NOT+organelle%3Amitochondrion";
  my $format = "fasta";

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
      $fragment_string = undef;
    }
  }

  if(exists($query_params->{'format'})) {
    $format = $query_params->{'format'};
  }

  # http://www.uniprot.org/uniprot/?query=existence%3A%22evidence+at+protein+level%22+OR+existence%3A%22evidence+at+transcript+level%22+AND+taxonomy%3A%22Mammalia+%5B40674%5D%22+NOT+taxonomy%3A%22Primates+%5B9443%5D%22&sort=score

  my $full_query = "wget -q -O - \"http://www.uniprot.org/uniprot/?query=";
  my %pe_code = (
                  '1' => 'evidence+at+protein+level',
                  '2' => 'evidence+at+transcript+level',
                  '3' => 'inferred+from+homology',
                  '4' => 'predicted',
                  '5' => 'uncertain',
                );

  # Must have file_name, pe_level, dest_dir and either taxon_id or taxonomy
  unless($file_name && $dest_dir && ($taxon_id || $taxon_group) && $pe_level) {
    $self->throw("Must define the following keys:\nfile_name\ntaxon_id or taxonomy\ndest_dir\npe_level");
  }

  my @pe_array = @{$pe_level};
  unless(scalar(@pe_array)) {
    die "Not PE levels found in value of pe_levels key. Format should be an array ref: ['1,2']";
  }

  foreach my $pe_level (@pe_array) {
    unless($pe_level =~ /\d+/) {
     die "Could not parse a PE level from the following: ".$pe_level;
    }

    my $parsed_pe_level = $&;
    unless($parsed_pe_level >= 1 && $parsed_pe_level <= 5) {
     die "Parsed PE level is outside the normal range of 1-5: ".$parsed_pe_level;
   }
   $pe_string .= 'existence%3A%22'.$pe_code{$pe_level}.'%22+OR+';
  }

  $pe_string =~ s/\+OR\+$/\)/;

  # NOTE this bit of the code with taxonomy and exclude is shit and needs to be upgraded
  if($taxon_id) {
    $taxonomy_string = '+AND+taxonomy%3A+'.$taxon_id;
  } elsif($taxon_group) {
    $taxonomy_string = '+AND+taxonomy%3A'.$taxon_group;
  }

#+NOT+taxonomy%3A%22
  if($exclude_id) {
    $exclude_string = '+NOT+taxonomy%3A+'.$exclude_id;
  } elsif($exclude_group) {
    $exclude_string = '+NOT+taxonomy%3A'.$exclude_group;
  }

  $full_query .= $pe_string.$taxonomy_string.$exclude_string.$fragment_string.$mito."&compress=".$compress."&format=".$format.
                 "\" > ".$dest_dir."/".$file_name;
  if($compress eq 'yes') {
    $full_query .= ".gz";
  }

  unless(-e $dest_dir) {
    `mkdir -p $dest_dir`;
  }

  return($full_query);

}

1;
