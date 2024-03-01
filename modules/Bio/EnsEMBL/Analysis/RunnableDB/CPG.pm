=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::CPG - 

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::CPG->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable CPG which wraps the program CPG

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the simple_feature table 

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::CPG;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::CPG;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::CPG
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence;
  $self->query($slice);
  if(!$self->analysis->program_file){
    $self->analysis->program_file('cpg');
  }
  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::CPG->new
    (
     -query => $self->query,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     %parameters,
    );
  $self->runnable($runnable);
  return 1;
}



=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::CPG
  Function  : get simple feature adaptor
  Returntype: Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor
  Exceptions: none
  Example   : 

=cut

sub get_adaptor{
  my ($self) = @_;
  return $self->db->get_SimpleFeatureAdaptor;
}


1;
