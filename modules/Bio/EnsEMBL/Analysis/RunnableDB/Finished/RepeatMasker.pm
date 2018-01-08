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


=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Finished::RepeatMasker

=head1 SYNOPSIS

  my $repeat_masker = Bio::EnsEMBL::Analysis::RunnableDB::Finished::RepeatMasker->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $repeat_masker->fetch_input;
  $repeat_masker->run;
  $repeat_masker->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable RepeatMasker which wraps the program RepeatMasker

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the repeat_feature and repeat_consensus tables

=head1 CONTACT

Post questions to : anacode@sanger.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::RepeatMasker;

use strict;
use warnings;
use Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker;
use Bio::EnsEMBL::Analysis::Runnable::Finished::RepeatMasker;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker);



=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Finished::RepeatMasker
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   :

=cut



sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence;
  $self->query($slice);
  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::RepeatMasker->new
    (
     -query => $self->query,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     %parameters,
    );
  $self->runnable($runnable);
  return 1;
}

1;
