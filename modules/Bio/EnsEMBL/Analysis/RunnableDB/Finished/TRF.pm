# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Finished::TRF

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $trf = Bio::EnsEMBL::Analysis::RunnableDB::TRF->new(
    -dbobj      => $db,
    -input_id   => $input_id
    -analysis   => $analysis
);
$trf->fetch_input();
$trf->run();
$trf->output();
$trf->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::TRF to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction of
parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is required
for database access.

=head1 CONTACT

Refactored by Sindhu K. Pillai B<sp1@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::TRF;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB::Finished;
use Bio::EnsEMBL::Analysis::Runnable::TRF;
use Bio::EnsEMBL::Analysis::Config::General;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Finished);

=head2 fetch_input

=cut

sub fetch_input {

    my( $self) = @_;
    throw("No input id") unless defined($self->input_id);
    my $sliceid  = $self->input_id;
    my $sa = $self->db->get_SliceAdaptor();
    my $slice   = $sa->fetch_by_name($sliceid);
    $slice->{'seq'}=$slice->seq();
    my %parameters      = %{$self->parameters_hash};
    $parameters{-trf}   = $self->analysis->program_file || undef;
    $parameters{-query} = $slice->get_repeatmasked_seq($ANALYSIS_RM_TRF,$SOFT_MASKING) or throw("Unable to fetch slice");
    $parameters{-analysis} = $self->analysis;
    my $runnable = new Bio::EnsEMBL::Analysis::Runnable::TRF(%parameters);
    $self->runnable($runnable);
    return 1;

}


sub write_output{

  my ($self) = @_;
  my @features   = @{$self->output()->[0]};
  my $repeat_f_a = $self->db->get_RepeatFeatureAdaptor();
  my $slice;
  eval {
    $slice = $self->db->get_SliceAdaptor->fetch_by_name($self->input_id);
  };
  if ($@) {
    print STDERR "Slice not found, skipping writing output to db: $@\n";
  }
  foreach my $f(@features){
    $f->analysis($self->analysis);
    $f->slice($slice);
    $repeat_f_a->store($f);
  }
}


1;
