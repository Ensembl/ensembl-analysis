=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam->
  new(
      -analysis => $analysis,
      -db => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

Modified blast runnable for specific use with RFAMSEQ.
Use for running BLASTN of genomic vs RFAMSEQ prior to 
ncRNA analysis using Infernal.
Slice size seems best around 200k

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB::Blast;
use Bio::EnsEMBL::Analysis::Runnable::BlastRfam;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::Databases qw(DATABASES DNA_DBNAME);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Blast Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Function  : fetch sequence out of database, instantiate the filter, 
  parser and finally the blast runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_; 
 
  # attach dna-db to db
  my $dna_db = $self->get_dbadaptor($DNA_DBNAME); 
  $self->db->dnadb($dna_db); 

  my $slice = $self->fetch_sequence($self->input_id, $self->db,['Dust']);
  $self->query($slice);
  my %blast = %{$self->BLAST_PARAMS};
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }
  my $seq = $self->query;
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastRfam->new
    (
     -query    => $seq,
     -program  => $self->analysis->program_file,
     -parser   => $parser,
     -filter   => $filter,
     -database => $self->analysis->db_file,
     -analysis => $self->analysis,
     %blast,
    );
  $self->runnable($runnable);
  return 1;
}


