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
=pod 

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Signalp

=head1 SYNOPSIS

  my $signalp = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Signalp->new ( -db      => $db,
    	  	                                                            -input_id   => $input_id,
                                                                            -analysis   => $analysis,
                                                                          );
  $signalp->fetch_input;  # gets sequence from DB
  $signalp->run;
  $signalp->write_output; # writes features to to DB

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Signalp;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Signalp;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);

sub fetch_input {
  my ($self, @args) = @_;

  $self->SUPER::fetch_input(@args);
  
  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Signalp->new(-query     => $self->query,
                                                                              -analysis  => $self->analysis);
  $self->runnable($run);
}


1;
