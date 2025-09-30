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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadAssembly;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions qw(catfile catdir);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadSeqRegions');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    convert_chromosome_names => 1,
    assembly_name => undef,
  }
}

sub fetch_input {
  my $self = shift;

  $self->hrdb_set_con($self->hrdb_get_dba($self->param_required('target_db')), 'target_db');

  my $path_to_files = $self->param_required('output_path');
  if ($self->param_is_defined('primary_assembly_dir_name')) {
    $path_to_files = catdir($self->param('output_path'), $self->param_required('primary_assembly_dir_name'), 'AGP');
  }
  if (!-d $path_to_files) {
    $self->complete_early("The AGP directory was not found. Assuming the assembly is single level. Path checked:\n".$path_to_files);
  }
  $self->param('input_dir', $path_to_files);
  if ($self->param_is_defined('chr2acc_file')) {
    my %acc_to_name;
    open(NF, $self->param('chr2acc_file')) or $self->throw("Can't open ".$self->param('chr2acc_file')." ".$!);
    while(<NF>){
      chomp;
      my @name_values = split(/\s+/,$_);
      $acc_to_name{$name_values[1]} = $name_values[0];
    }
    close(NF) || $self->throw('Could not open '.$self->param('chr2acc_file'));
    $self->param('acc_to_name', \%acc_to_name);
  }
}

sub run {
  my $self = shift;

  my $chromo_present = 0;

  my %assembled = (
    $self->param('scaffold_contig') => 'scaffold',
    $self->param('chromosome_contig') => 'chromosome',
    $self->param('chromosome_scaffold')  => 'chromosome',
  );
  my %component = (
    $self->param('scaffold_contig') => 'contig',
    $self->param('chromosome_contig') => 'contig',
    $self->param('chromosome_scaffold')  => 'scaffold',
  );
  my $db = $self->hrdb_get_con('target_db');
  my $csa = $db->get_CoordSystemAdaptor;
  my $sa = $db->get_SliceAdaptor;
  my $convert_chrom_names = $self->param('convert_chromosome_names');
  my @assembly_mapping;
  my $acc_to_name = $self->param_is_defined('acc_to_name') ? $self->param('acc_to_name') : {};

# DO NOT CHANGE THE ORDER as the code in write_output relies on the file to be in the chr2contig, chr2scaffold, scaffold2contig order
  foreach my $filename ($self->param('chromosome_contig'), $self->param('chromosome_scaffold'), $self->param('scaffold_contig')) {
    my $file = catfile($self->param('input_dir'), $filename);
    if (-e $file) {
      ++$chromo_present if ($assembled{$filename} eq 'chromosome');

      my $assembled_cs = $csa->fetch_by_name($assembled{$filename},
                                             $self->param('assembly_name'));
      my $component_cs = $csa->fetch_by_name($component{$filename},
                                             $component{$filename} eq 'contig' ? undef : $self->param('assembly_name'));

      my %assembled_ids;
      my %component_ids;
      my $mapping_delimiter = '|';
      my @rows;

      open(FH, $file) or $self->throw("Can't open $file");
      while(<FH>) {
        next if /^\#/;
        chomp;

        #CM000686.1  28033929  28215811  243 F AC007965.3  1 181883  +
        #CM000686.1  28215812  28388853  244 F AC006991.3  1 173042  +
        #cb25.fpc4250	119836	151061	13	W	c004100191.Contig2	1	31226	+
        #cb25.fpc4250	151062	152023	14	N	962	fragment	yes
       say "FERGAL: ".$_;

       my ($a_name_init, $a_start, $a_end, $ordinal, $type, $c_name, $c_start,
            $c_end, $ori) = split;

        if($type and ($type ne 'N' && $type ne 'U')) {
          my $strand = 1;
          $strand = -1 if ($ori eq '-');

          my $a_name = $a_name_init;
          $a_name = $acc_to_name->{$a_name_init} if (exists $acc_to_name->{$a_name_init});

          my ($a_id, $c_id);
          if (!$assembled_ids{$a_name}) {
            if($convert_chrom_names and ($a_name =~ /^chr(\S+)/)){
              $a_name = $1;
            }
            my $a_piece = $sa->fetch_by_region($assembled_cs->name, $a_name,
                                               undef, undef, undef,
                                               $assembled_cs->version);
            $self->throw($a_name." doesn't seem to exist in the database\n")
              unless ($a_piece);
            if($a_piece->length < $a_end){
              $self->throw($a_name." apparent length ".$a_end. " is longer than ".
                    "the length in the current database ".
                    $a_piece->length."\n");
            }
            $a_id = $sa->get_seq_region_id($a_piece);
            $assembled_ids{$a_name} = $a_id;
          }
          else {
            $a_id = $assembled_ids{$a_name};
          }
          if($component_ids{$c_name}){
            $self->warning("You are already using ".$c_name." in another place ".
                 "in your assembly are you sure you want to\n");
            $mapping_delimiter = '#';
            $c_id = $component_ids{$c_name};
          }
          else {
            my $c_piece = $sa->fetch_by_region($component_cs->name, $c_name,
                                               undef, undef, undef,
                                               $component_cs->version);
            $self->throw($c_name." doesn't seem to exist in the database\n")
              unless ($c_piece);
            if($c_piece->length < $c_end){
              $self->throw($c_name." apparent length ".$c_end. " is longer than ".
                    "the length in the current database ".
                    $c_piece->length."\n");
            }
            $c_id = $sa->get_seq_region_id($c_piece);
            $component_ids{$c_name} = $c_id;
          }
          push(@rows, [$a_id, $a_start, $a_end, $c_id, $c_start, $c_end, $strand]);
        }
      }
      close(FH) || $self->throw("Could not close $file");

      my $mapping_string = $assembled_cs->name;
      $mapping_string .= ":".$assembled_cs->version if($assembled_cs->version);
      $mapping_string .= $mapping_delimiter.$component_cs->name;
      $mapping_string .= ":".$component_cs->version if($component_cs->version);
      push(@assembly_mapping, [$mapping_string, \@rows]);
    }
  }

  # Have to ask about this scenario and whether it warrants a warning or a throw. It may not
  # actually be needed. But if it is then I need to make it more useful
  if ($chromo_present) {
    if (@assembly_mapping != 3) {
      $self->warning(@assembly_mapping.' assembly mapping found instead of 3, we will need to duplicate contigs to add the missing scaffolds and create the mappings.');
    }
  }
  else { # if there is no chromosome level
    if (@assembly_mapping != 1) {
      $self->trhow(@assembly_mapping.' assembly mapping found instead of 2: '.join(' ', @assembly_mapping));
    }
  }
  $self->output(\@assembly_mapping);
}

sub write_output {
  my $self = shift;

  my $db = $self->hrdb_get_con('target_db');
  my $assembly_mapping = $self->output;
  if (@$assembly_mapping < 3 and $assembly_mapping->[0]->[0] =~ /chromosome/) {
    my $csa = $db->get_CoordSystemAdaptor;
    my $slice_adaptor = $db->get_SliceAdaptor;
    my $scaffold_cs = $csa->fetch_by_name('scaffold');
    if (!$scaffold_cs) {
      my $chr_cs = $csa->fetch_by_name('chromosome');
      $self->throw('Could not find the chromosome coordinate system, something is wrong')
        unless ($chr_cs);
      my $contig_cs = $csa->fetch_by_name('contig');
      if ($contig_cs->rank == 2) {
        $csa->remove($contig_cs);
        $contig_cs->rank(3);
        $csa->store($contig_cs);
        foreach my $slice (@{$slice_adaptor->fetch_all('contig')}) {
          $slice->coord_system($contig_cs);
          $slice_adaptor->update($slice);
        }
      }
      $scaffold_cs = Bio::EnsEMBL::CoordSystem->new(
        -name => 'scaffold',
        -version => $chr_cs->version,
        -default => 1,
        -sequence_level => 0,
        -rank => 2,
      );
      $csa->store($scaffold_cs);
    }
    my $mapping_string = $assembly_mapping->[0]->[0];
    $mapping_string =~ s/contig/scaffold/;
    $assembly_mapping->[2]->[0] = $mapping_string.':'.$scaffold_cs->version;
    if (!defined $assembly_mapping->[1]) {
      $mapping_string = $assembly_mapping->[0]->[0];
      $mapping_string =~ s/chromosome/scaffold/;
      $assembly_mapping->[1]->[0] = $mapping_string;
    }

    foreach my $row (@{$assembly_mapping->[0]->[1]}) {
      my $contig = $slice_adaptor->fetch_by_seq_region_id($row->[3]);
      my $scaffold = $slice_adaptor->fetch_by_region('scaffold', $contig->seq_region_name);
      if (!$scaffold) {
        $scaffold = Bio::EnsEMBL::Slice->new(
          -seq_region_name => $contig->seq_region_name,
          -coord_system => $scaffold_cs,
          -start => $contig->start,
          -end => $contig->end,
        );
        $slice_adaptor->store($scaffold);
      }
      push(@{$assembly_mapping->[2]->[1]}, [$row->[0], $row->[1], $row->[2], $scaffold->get_seq_region_id, $row->[4], $row->[5], 1]);
      push(@{$assembly_mapping->[1]->[1]}, [$scaffold->get_seq_region_id, $row->[1], $row->[2], $row->[3], $row->[4], $row->[5], $row->[6]]);
    }
  }
  my $sql = 'INSERT INTO assembly (asm_seq_region_id, asm_start, asm_end, cmp_seq_region_id, cmp_start, cmp_end, ori) VALUES(?, ?, ?, ?, ?, ?, ?)';
  my $sth = $db->dbc->prepare($sql);
  my $mc = $db->get_MetaContainer();
  foreach my $mappings (@{$assembly_mapping}) {
    foreach my $row (@{$mappings->[1]}) {
      $sth->bind_param(1, $row->[0]);
      $sth->bind_param(2, $row->[1]);
      $sth->bind_param(3, $row->[2]);
      $sth->bind_param(4, $row->[3]);
      $sth->bind_param(5, $row->[4]);
      $sth->bind_param(6, $row->[5]);
      $sth->bind_param(7, $row->[6]);
      $sth->execute();
    }
    $mc->store_key_value('assembly.mapping', $mappings->[0]);
  }
  $sth->finish();
}

1;
