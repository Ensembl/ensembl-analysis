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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadAssemblyComponents

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadAssemblyComponents;

use strict;
use warnings;

use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Attribute;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    chromosome_level => 'chromosome',
    scaffold_level => 'scaffold',
    contig_level => 'contig',
  }
}


sub fetch_input {
  my ($self) = @_;

  my $file = $self->param_required('filename');
  $self->throw("$file does not exist") unless (-e $file);
  my $db = $self->get_database_by_name('target_db');
  $self->hrdb_set_con($db, 'target_db');
  my %coord_systems;
  my $coord_system_adaptor = $db->get_CoordSystemAdaptor;
  my $contig_rank = 3;
  my $cs = $coord_system_adaptor->fetch_by_name($self->param('chromosome_level'), $self->param_required('assembly_name'));
  if ($cs) {
    $coord_systems{$self->param('chromosome_level')} = $cs;
  }
  else {
    --$contig_rank;
  }
  $cs = $coord_system_adaptor->fetch_by_name($self->param('scaffold_level'), $self->param_required('assembly_name'));
  if ($cs) {
    $coord_systems{$self->param('scaffold_level')} = $cs;
  }
  else {
    $coord_systems{$self->param('scaffold_level')} = Bio::EnsEMBL::CoordSystem->new(
      -name => $self->param('scaffold_level'),
      -version => $self->param('assembly_name'),
      -default => 1,
      -rank => $contig_rank-1,
    );
  }
  my $sequence_level_cs = $coord_system_adaptor->fetch_sequence_level;
  if ($sequence_level_cs) {
    $coord_systems{$sequence_level_cs->name} = $sequence_level_cs;
    $self->param('contig_level', $sequence_level_cs->name);
  }
  else {
    $coord_systems{$self->param('contig_level')} = Bio::EnsEMBL::CoordSystem->new(
      -name => $self->param('contig_level'),
      -default => 1,
      -rank => $contig_rank,
      -sequence_level => 1,
    );
  }
  $self->param('coord_systems', \%coord_systems);
  if (!$self->param_is_defined('wgs_project')) {
    $self->param('wgs_project', undef);
  }
}

sub run {
  my ($self) = @_;

  my $coord_systems = $self->param('coord_systems');
  my $asm_coordinate_system = $self->param('filename') =~ /^chr/ ?
    $coord_systems->{$self->param('chromosome_level')} : $coord_systems->{$self->param('scaffold_level')};
  my $cmp_coordinate_system = $self->param('filename') =~ /comp/ ?
    $coord_systems->{$self->param('contig_level')} : $coord_systems->{$self->param('scaffold_level')};
  my %cmp_slices;
  my $multiple_mapping = 0;
  my @output;
  my @iids;
  my $wgs_project = $self->param('wgs_project');
  open(RH, $self->param('filename')) || $self->throw('Could not open AGP file');
  while (my $line = <RH>) {
    $line =~ s/\R$//;
    if ($line =~ /^#/) {
      if ($line =~ /(\w+)-from-(\w+)/) {
        if ($1 =~ /chromosome/i) {
          $asm_coordinate_system = $coord_systems->{$self->param('chromosome_level')};
        }
        else {
          $asm_coordinate_system = $coord_systems->{$self->param('scaffold_level')};
        }
        if ($2 =~ /components/i) {
          $cmp_coordinate_system = $coord_systems->{$self->param('contig_level')};
        }
        else {
          $cmp_coordinate_system = $coord_systems->{$self->param('scaffold_level')};
        }
      }
    }
    else {
      my @data = split("\t", $line);
      next if ($data[4] eq 'N');
      my $slice;
      my $cmp_slice;
      if (exists $cmp_slices{$data[5]}) {
        $multiple_mapping = 1;
      }
      push(@output, [$data[0], $data[1], $data[2], $data[5], $data[6], $data[7], ($data[8] eq '-' ? -1 : 1)]);
      push(@iids, $data[5]) unless ($wgs_project and $data[5] =~ /^$wgs_project/);
    }
  }
  close(RH) || $self->throw('Could not close AGP file');
  $self->param('multiple_mapping', $multiple_mapping);
  $self->param('asm_coordinate_system', $asm_coordinate_system);
  $self->param('cmp_coordinate_system', $cmp_coordinate_system);
  $self->output(\@output);
  $self->param('contigs_to_download', \@iids) if (@iids);
}

sub write_output {
  my ($self) = @_;

  my $db = $self->hrdb_get_con('target_db');
  my $asm_coordinate_system = $self->param('asm_coordinate_system');
  if (!$asm_coordinate_system->adaptor) {
    eval {
      $db->get_CoordSystemAdaptor->store($asm_coordinate_system);
    };
    if ($@) {
      $asm_coordinate_system = $db->get_CoordSystemAdaptor->fetch_by_name($asm_coordinate_system->name, $asm_coordinate_system->version);
    }
  }
  my $cmp_coordinate_system = $self->param('cmp_coordinate_system');
  if (!$cmp_coordinate_system->adaptor) {
    eval {
      $db->get_CoordSystemAdaptor->store($cmp_coordinate_system);
    };
    if ($@) {
      $cmp_coordinate_system = $db->get_CoordSystemAdaptor->fetch_by_name($cmp_coordinate_system->name, $cmp_coordinate_system->version);
    }
  }
  my $meta_adaptor = $db->get_MetaContainerAdaptor;
  my $string = '|';
  if ($self->param('multiple_mapping')) {
    $string = '#';
  }
  $meta_adaptor->store_key_value('assembly.mapping',
    $asm_coordinate_system->name.':'.$asm_coordinate_system->version.$string.
    $cmp_coordinate_system->name.($cmp_coordinate_system->verion ? ':'.$cmp_coordinate_system->version : '')
  );
  my %asm_slices;
  my %cmp_slices;
  my $sth = $db->dbc->prepare('INSERT INTO assembly VALUES (?, ?, ?, ?, ?, ?, ?)');
  foreach my $data (@{$self->output}) {
    my $slice;
    my $cmp_slice;
    if (exists $asm_slices{$data->[0]}) {
      $slice = $asm_slices{$data->[0]};
      if ($slice->end < $data->[2]) {
        $slice->end($data->[2]);
        $slice->adaptor->update($slice);
      }
    }
    else {
      $slice = $self->get_Slice($data->[0], $data->[2], $asm_coordinate_system);
      $asm_slices{$data->[0]} = $slice;
    }
    if (exists $cmp_slices{$data->[3]}) {
      $cmp_slice = $cmp_slices{$data->[3]};
      if ($slice->end < $data->[5]) {
        $cmp_slice->end($data->[5]);
        $slice->adaptor->update($cmp_slice);
      }
    }
    else {
      $cmp_slice = $self->get_Slice($data->[3], $data->[5], $cmp_coordinate_system);
      $cmp_slices{$data->[3]} = $slice;
    }
    $sth->bind_param(1, $slice->get_seq_region_id);
    $sth->bind_param(2, $data->[1]);
    $sth->bind_param(3, $data->[2]);
    $sth->bind_param(4, $cmp_slice->get_seq_region_id);
    $sth->bind_param(5, $data->[4]);
    $sth->bind_param(6, $data->[5]);
    $sth->bind_param(7, $data->[6]);
    $sth->execute();
  }
  $sth->finish();
  $self->dataflow_output_id({accessions => $self->param('contigs_to_download')}, $self->param('_branch_to_flow_to'))
    if ($self->param_is_defined('contigs_to_download'));
}


sub get_Slice {
  my ($self, $name, $end, $cs) = @_;
  
  my $slice_adaptor = $cs->adaptor->db->get_SliceAdaptor;
  my $slice = $slice_adaptor->fetch_by_region($cs->name, $name);
  if (!$slice) {
    $slice = Bio::EnsEMBL::Slice->new(
      -seq_region_name => $name,
      -coord_system => $cs,
      -start => 1,
      -end => $end,
    );
    eval {
      if ($cs->is_sequence_level) {
        my $seq = 'N'x$slice->length;
        $slice_adaptor->store($slice, \$seq);
      }
      else {
        $slice_adaptor->store($slice);
      }
    };
    if ($@) {
      $slice = $slice_adaptor->fetch_by_region($cs->name, $name);
    }
  }
  return $slice;
}

1;
