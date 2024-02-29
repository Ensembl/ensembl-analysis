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

Bio::EnsEMBL::Analysis::RunnableDB::EPCR - 

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::EPCR->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::EPCR;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::EPCR;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::EPCR
  Function  : fetch data out of database and create runnable
  Returntype: none
  Exceptions: throws if no markers in the database when no -STS_FILE
  is defined
  Example   : 

=cut

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
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::EPCR->new
    (
     -query => $slice,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     %parameters
    );
  $self->runnable($runnable);
}


=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::EPCR
  Function  : get marker feature adaptor
  Returntype: Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor
  Exceptions: none
  Example   : 

=cut

sub get_adaptor{
  my ($self) = @_;
  return $self->db->get_MarkerFeatureAdaptor;
}
