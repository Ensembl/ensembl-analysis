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

  Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::PrositePattern

=head1 SYNOPSIS

  my $tmhmm = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::PrositePattern->
    new ( 
    -db      => $db,
    -input_id   => $input_id,
    -analysis   => $analysis)
    );
  $tmhmm->fetch_input;  # gets sequence from DB
  $tmhmm->run;
  $tmhmm->write_output; # writes features to to DB

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::PrositePattern;

use warnings ;
use strict;
use vars qw(@ISA);


use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositePattern;
use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self, @args) = @_;
  
  $self->SUPER::fetch_input(@args);

  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositePattern->
      new(-query     => $self->query,
          -analysis  => $self->analysis,
          %{$self->parameters_hash}
          );
  $self->runnable($run);
}

1;
