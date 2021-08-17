=head1 LICENSE

 Copyright [2021] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::Star

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::GenerateSlices->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

A test module for generating slices outside of HiveSubmitAnalysis

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::GenerateSlices;

use warnings;
use strict;
use feature 'say';
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(hrdb_get_dba);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _branch_to_flow_to => 2,
    threads => 1,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Parameters for fetching slices
 Returntype : None
 Exceptions : Throws if 'filename' does not exist
              Throws if 'fastqpair' does not exist
              Throws if 'wide_short_read_aligner' is not defined

=cut

sub fetch_input {
  my ($self) = @_;

  my $dba = hrdb_get_dba($self->param('reference_db'));
  my $slice_adaptor = $dba->get_SliceAdaptor();
  my $slices = $slice_adaptor->fetch_all('toplevel');
  $self->param('slices',$slices);
}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $slices = $self->param('slices');
  foreach my $slice (@$slices) {
    if($slice->name =~ /\:MT\:/ or not $slice->name =~ /^chromosome/) {
      next;
    }

    say "Output slice name: ".$slice->name;
    $self->dataflow_output_id([{'iid' => [$slice->name],'core_db' => $self->param('core_db'),'genome_index' => $self->param('genome_index')}], $self->param('_branch_to_flow_to'));
  }
}

1;
