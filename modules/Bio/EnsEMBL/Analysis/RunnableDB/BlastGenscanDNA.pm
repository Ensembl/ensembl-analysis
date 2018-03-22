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

Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanDNA - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanDNA->
  new(
      -analysis => $analysis,
      -db => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

 This module inherits from the Blast runnable and instantiates 
 BlastTranscriptDNA passing in prediction transcript

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanDNA;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA;
use Bio::EnsEMBL::Analysis::RunnableDB::Blast;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Blast);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanDNA 
  Function  : fetch sequence and prediction transcripts of database, 
  read config files instantiate the filter, parser and finally the blast 
  runnables
  Returntype: none
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence($self->input_id, $self->db);
  $self->query($slice);
  my %blast = %{$self->BLAST_PARAMS};
  my $logic_names = $BLAST_AB_INITIO_LOGICNAME;
  my @pts ;
  my $pta = $self->db->get_PredictionTranscriptAdaptor;
  $logic_names = ['Genscan'] if(scalar(@$logic_names) == 0);
  foreach my $logic_name (@$logic_names) {
    my $pt = $pta->fetch_all_by_Slice($self->query, $logic_name);
    push @pts, @$pt ;
  }
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }
  foreach my $t(@pts){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA->
      new(
          -transcript => $t,
          -query => $self->query,
          -program => $self->analysis->program_file,
          -parser => $parser,
          -filter => $filter,
          -database => $self->analysis->db_file,
          -analysis => $self->analysis,
          %blast,
         );
    $self->runnable($runnable);
  }
}


