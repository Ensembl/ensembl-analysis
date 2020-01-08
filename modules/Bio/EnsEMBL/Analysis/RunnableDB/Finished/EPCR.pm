# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

=head1 NAME

EPCR.pm 

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Finished::EPCR->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;


=head1 DESCRIPTION

The Finished version of EPCR.

=head1 CONTACT

anacode@sanger.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::EPCR;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB::EPCR;
use Bio::EnsEMBL::Analysis::Runnable::Finished::EPCR;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::EPCR);

sub fetch_input{
  my ($self) = @_;
  my %parameters = %{$self->parameters_hash};
  if($self->analysis->db_file){
    $parameters{'-STS_FILE'} = $self->analysis->db_file 
      unless($parameters{'-STS_FILE'});
  }
  if(!$parameters{'-STS_FILE'}){
    my $sts = $self->db->get_MarkerAdaptor->fetch_all;
    throw("No markers in ".$self->db->dbname) unless(@$sts);
    $parameters{'-STS_FEATURES'} = $sts;
  }
  my $slice = $self->fetch_sequence;
  $self->query($slice);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::EPCR->new
    (
     -query => $slice,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     %parameters
    );
  $self->runnable($runnable);
}

1;
