=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::CheckTranscriptomicRegistry;

use strict;
use warnings;
use POSIX;
use List::Util qw( min max );
use feature 'say';
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB'); 

sub fetch_input{
  my ($self) = @_;
}

sub run{
  my ($self) = @_;
  my $status = 0;
  my $taxon_id = $self->param('taxon_id');
  my $sth = $self->data_dbc()->prepare("SELECT * from short_read_data WHERE species_id =?");
  $sth->bind_param(1,$taxon_id);
  if ($sth->execute){
   while (my @result = $sth->fetchrow_array()){
    $status = 1;
   }
  }
  else{
   $self->throw("Could not determine the transcriptomic data status for this species $taxon_id" );
  }
  $self->param('transcriptomic_status',$status);
	 
} 

sub write_output{
	 my ($self) = @_;
	
}
1;
