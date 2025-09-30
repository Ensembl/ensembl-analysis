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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveFindContigAccessions;

use strict;
use warnings;
use feature 'say';


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  return 1;
  unless($self->param('wgs_id') && $self->param('output_path') && $self->param('contigs_source')) {
    $self->throw("Must pass in the following parameters:\n".
          "wgs_id e.g AAEX for Dog".
          "output_path e.g /path/to/work/dir\n".
          "contigs_source e.g. 'ENA' or 'NCBI'");
  }

}

sub run {
  my $self = shift;

  my $primary_assembly_dir_name = $self->param('primary_assembly_dir_name');

  my $agp_files_path = $self->param('output_path')."/".$primary_assembly_dir_name."/AGP";

  unless(-e $agp_files_path) {
    $self->warning("Could not find an agp file path, therefore assuming the assembly is single level");
    return;
  }

  my $contig_accessions = $self->find_contig_accessions($agp_files_path);

  open(OUT,">".$agp_files_path."/contigs.txt");
  foreach my $contig_accession (keys(%{$contig_accessions})) {
    say OUT $contig_accession;
  }
  close OUT;

}

sub find_contig_accessions {
  my ($self,$agp_files_path) = @_;

#  my $contig_prefix_list = {};
  my $contig_accessions = {};
  my @agp_files = glob($agp_files_path."/*.agp");
  foreach my $agp_file (@agp_files) {

    unless($agp_file =~ /\.comp\.agp/ || $agp_file =~ /\.placed\.scaf\.agp/ ||
           $agp_file =~ /\.unlocalized\.scaf\.agp/ || $agp_file =~ /unplaced\.scaf\.agp/) {
     say "Skipping file: ".$agp_file;
     next;
    }

    say "File: ".$agp_file;
    open(IN,$agp_file);
    my @agp_file = <IN>;
    close IN;
    foreach my $row (@agp_file) {
      if($row =~ /^\#/) {
        next;
      }

      my @columns = split("\t",$row);

      unless($columns[5]) {
        next;
      }

      if($columns[4] eq 'N' || $columns[4] eq 'U') {
        next;
      }

      my $contig_accession = $columns[5];
      if($contig_accession =~ /KQ/) {
        say "FM2 Row: ".$row;
      }
#      my $contig_prefix = substr($contig_name,0,5);
      unless($contig_accessions->{$contig_accession}) {
        $contig_accessions->{$contig_accession} = 1;
      }

#      unless($contig_prefix_list->{$contig_prefix}) {
#        $contig_prefix_list->{$contig_prefix} = 1;
#      }
    }
  }

#  my $contig_accessions = [keys(%{$contig_names})];

  return $contig_accessions;
}


sub write_output {
  my $self = shift;

  return 1;
}

1;
