#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny qw(path);
use Data::Dumper;

sub filter_rfam {
  my ($cm, $rfam_acc) = @_;

  my @filtered;
  my @cm_models = split(/\/\/\n/, $cm);
  
  my %rfam_acc = map { $_ => 1 } @$rfam_acc;
  foreach my $cm_model (@cm_models) {
    $cm_model =~ m/(RF\d+)/ ;
    my $rfam = $1;
    print $rfam . "\t";
    if (exists($rfam_acc{$rfam})) {
      push @filtered, $cm_model;
    } else {
      # print("Rfam model $rfam removed by filtering.\n");
    }
  }

  return join("//\n", @filtered)."//\n";
}

my $rfam_cm_file       = $ARGV[0]; #$self->param_required('rfam_cm_file');
my $rfam_accessions   = $ARGV[1]; #$self->param_required('rfam_accessions');
my $working_dir        = $ARGV[2]; #$self->param_required("working_dir");

my $cm_path = path($rfam_cm_file);
my $ra_path = path($rfam_accessions);
my $output = path($working_dir . "/Rfam.cm");

my $cm = $cm_path->slurp;
my $ra = $ra_path->slurp;

my @accessions = split(/\n/, $ra);

$cm = filter_rfam($cm, \@accessions);
$output->spew($cm);


