=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatCoverage

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatCoverage;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::RepeatCoverage qw(get_genes get_repeat_blocks get_masked_stats);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command send_email);

=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
               coord_system_name => 'toplevel',
               coord_system_version => undef,
               include_non_reference => 0,
               include_duplicates => 0,
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    coord_system_name => 'toplevel',
    coord_system_version => undef,
    include_non_reference => 0,
    include_duplicates => 0,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Connect to 'source_db' database to retrieve all slices. It will check
              that 'repeat_logic_names' is an array ref and warn if one of the logic
              name is "trf".
              It uses 'dna_db' if present but it should not be needed.
              If 'gene_db' is set, retrieve the slice adaptor so it can fetch genes
              later
 Returntype : None
 Exceptions : Throws if 'source_db' is not set
              Throws if 'repeat_logic_names' is not set or not an array ref

=cut

sub fetch_input {
  my ($self) = @_;

  my $db = $self->get_database_by_name('source_db');
  if ($self->param_is_defined('dna_db')) {
    $db->dnadb($self->get_database_by_name('dna_db'));
  }
  if ($self->param_is_defined('gene_db')) {
    my $genedb = $self->get_database_by_name('gene_db');
    $self->param('gene_slice_adaptor', $genedb->get_SliceAdaptor);
  }
  my $repeats = $self->param_required('repeat_logic_names');
  $self->throw('"repeat_logic_names" should be an arrayref') unless (ref($repeats) eq 'ARRAY');
  foreach my $repeattype (@$repeats) {
    $self->warning('We usually do not use trf when calculating the coverage') if ($repeattype =~ /^trf$/i);
  }

  my $seq_regions = $db->get_SliceAdaptor->fetch_all(
                      $self->param('coord_system_name'),
                      $self->param('coord_system_version'),
                      $self->param('include_non_reference'),
                      $self->param('include_duplicates')
                    );
  $self->param('seq_regions', $seq_regions);
}


=head2 run

 Arg [1]    : None
 Description: Fetch all repeats and genes for each toplevel region
              and calculate the number of bases repeat masked and the
              number of genes which have at least one exon overlapping
              a repeat
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $total_masked = 0;
  my $total_genes = 0;
  my $total_overgene = 0;
  my $total_bases = 0;
  my $genesa;
  if ($self->param_is_defined('gene_slice_adaptor')) {
    $genesa = $self->param('gene_slice_adaptor');
  }
  my $repeattypes = $self->param('repeat_logic_names');
  foreach my $chr (@{$self->param('seq_regions')}) {
    my $genes;
    if ($genesa) {
      $genes = get_genes($genesa->fetch_by_name($chr->name), $repeattypes);
      $total_genes += scalar(@$genes);
    } else {
      print "NO genedb specified so no genes fetched\n";
    }

    print "Fetching repeats\n";

    my @repeats;
    foreach my $repeattype (@$repeattypes) {
      push(@repeats, @{$chr->get_all_RepeatFeatures($repeattype)});
    }
    print "Done fetching repeats (fetched " . scalar(@repeats) .")\n";
    my $repeat_blocks = get_repeat_blocks(\@repeats);

    my ($nmasked, $novergene) = get_masked_stats($repeat_blocks, $genes);

    $total_overgene += $novergene;
    $total_masked += $nmasked;
    $total_bases += $chr->length;
  }
  $self->param('bases', $total_bases);
  $self->param('masked', $total_masked);
  my $ratio = ($self->param('masked') / $self->param('bases')) * 100  ;
  $self->param('repeats_ratio', $ratio);
  if ($total_genes) {
    $self->param('genes', $total_genes);
    $self->param('genes_overlapped', $total_overgene);
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Create the message and send an email the specified address
              if 'email' is set. Other wise it write the result as a warning.
              It dataflows the repeats ratio on branch '_branch_to_flow_to'
              with the key 'repeat_mask_coverage'. If there is something extreme, 
              job will fail. 
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $db = $self->get_database_by_name('source_db');
  my $dbc = $db->dbc;
  my $dbname = $dbc->dbname;
  my $msg = $dbname."\nAnalyses:";
  my $repeats = $self->param_required('repeat_logic_names');
  foreach my $repeattype (@$repeats) {
    $msg .= " ".$repeattype;
  }

  $msg .= sprintf("\nTotal bases = %d\nTotal masked = %d\t( %.2f%% masked)\n", $self->param('bases'), $self->param('masked'), $self->param('repeats_ratio'));
  if ($self->param_is_defined('genes')) {
    $msg .= "Total genes = ".$self->param('genes')."\nTotal genes overlapped = ".$self->param('genes_overlapped')."\n";
  }
  if ($self->param_is_defined('email')) {
    send_email($self->param('email'), $self->param('email'), 'AUTOMATED REPORT - Repeat coverage stats', $msg);
  }
  else {
    $self->warning($msg)
  }
  
  $self->dataflow_output_id({repeat_mask_coverage => $self->param('repeats_ratio')}, $self->param('_branch_to_flow_to'));
}

1;
