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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveIndexFlatFile

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveIndexFlatFile;

use strict;
use warnings;

use File::Spec::Functions qw (splitpath);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    batch_size => 1,
    primary_namespace => 'ACC',
    primary_pattern => '>(\S+)',
    primary_start_pattern => '^>',
    format => 'fasta',
  }
}


sub fetch_input {
  my ($self) = @_;

  my (undef, $index_dir, $database_name) = splitpath($self->param_required('seqfetcher_index'));
  if ($self->param_required('input_files')) {
    $self->throw('You need to specify an array ref') unless (ref($self->param('input_files')) eq 'ARRAY');
  }
  my $indexer;
  if ($self->param_is_defined('program')) {
    my $class = $self->require_module('Bio::EnsEMBL::Analysis::Runnable::Indicate');
    $indexer = $class->new;
  }
  else {
    my $class = $self->require_module('Bio::DB::Flat::BinarySearch');
    $indexer = $class->new(
      -directory         => $index_dir,
      -format            => $self->param('format'),
      -primary_pattern   => $self->param('primary_pattern'),
      -primary_namespace => $self->param('primary_namespace'),
      -start_pattern     => $self->param('primary_start_pattern'),
      -dbname            => $database_name,
      -write_flag        => 1,
    );
    if ($self->param_is_defined('secondary_namespaces')) {
      if (ref($self->param('secondary_namespaces')) eq 'ARRAY') {
        $indexer->secondary_namespaces($self->param('secondary_namespaces'));
      }
      else {
        $self->throw('You should provide an arrayref of string');
      }
      if (ref($self->param_required('secondary_patterns')) eq 'HASH') {
        $indexer->secondary_patterns($self->param('secondary_patterns'));
      }
      else {
        $self->throw('You should provide an hashref where the key is a secondary namespace and the value is a regex');
      }
    }
  }
  $self->param('indexer', $indexer);
}


sub run {
  my ($self) = @_;

  return 1;
}


sub write_output {
  my ($self) = @_;

  my $indexer = $self->param('indexer');
  $indexer->build_index(@{$self->param('input_files')});
  my @iids;
  open(FH, $indexer->primary_index_file) || $self->throw("Could not open index file");
  my $line = substr(<FH>, 4);
  close(FH) || $self->throw("Could not close index file");
  while ($line =~ /\s*(\S+)\t\s*\d+\t\s*\d+\t\s*\d+/g) {
    push(@iids, $1);
  }
  my @output_ids;
  my $seqfetcher_index = $self->param('seqfetcher_index');
  if ($self->param('batch_size') > 1) {
    my $batch_size = $self->param('batch_size');
    my $index = 0;
    push(@output_ids, {seqfetcher_index => $seqfetcher_index});
    foreach my $iid (@iids) {
      if (++$index > $batch_size) {
        push(@output_ids, {iid => [$iid], seqfetcher_index => $seqfetcher_index});
        $index = 0;
      }
      else {
        push(@{$output_ids[-1]->{iid}}, $iid);
      }
    }
  }
  else {
    foreach my $iid (@iids) {
      push(@output_ids, {iid => $_, seqfetcher_index => $seqfetcher_index});
    }
  }
  $self->dataflow_output_id(\@output_ids, $self->param('_branch_to_flow_to'));
}

1;
