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

Bio::EnsEMBL::Analysis::RunnableDB::Genefinder - 

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Genefinder->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;


=head1 DESCRIPTION

fetches sequence data from database an instantiates and runs the
genefinder runnable


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Genefinder;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB::Genscan;
use Bio::EnsEMBL::Analysis::Runnable::Genefinder;
use Bio::EnsEMBL::Analysis::Config::General;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Genscan);



=head2 runnable_path

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genefinder
  Function  : return the runnable path
  Returntype: string
  Exceptions: 
  Example   : my $runnable = $self->runnable_path->new
                               (
                                -query    => $self->query,
                                -program  => $self->analysis->program_file,
                                -analysis => $self->analysis,
                                %parameters,
                               );

=cut


sub runnable_path{
  my ($self) = @_;
  return "Bio::EnsEMBL::Analysis::Runnable::Genefinder";
}




=head2 standard_args

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genefinder
  Function  : to return a hash which contains the 
  standard constructor args for the genscan runnable
  Returntype: hashref
  Exceptions: none
  Example   : 

=cut

sub standard_args{
  my ($self) = @_;
  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }
  return {
          -query => $self->query,
          -program => $self->analysis->program_file,
          -analysis => $self->analysis,
          %parameters,
         };
}

1;
