=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Accumulator - 

=head1 SYNOPSIS

  my $accumulator = Bio::EnsEMBL::Analysis::RunnableDB::Accumulator->
  new(
      -input_id => 'ACCUMULATOR',
      -db => $db,
      -analysis => $analysis,
     );
  $accumulator->fetch_input;
  $accumulator->run;
  $accumulator->write_output;

=head1 DESCRIPTION

This is a simple place holder module to allow the accumulator wait for all
stages in the pipeline to work. It does nothing just

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::Accumulator;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Dummy method to comply to the interface
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    throw("No input id") unless defined($self->input_id);

    return 1;

}

sub run {
    my ($self) = @_;
    print "Dummy RunnableDB - no runnable to run\n";

}

sub write_output {
    my ($self) = @_;

    print "Dummy RunnableDB - no output to write\n";

    return 1;
}

1;
