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

Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep->
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
 BlastTranscriptPep passing in prediction transcript

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;
use Bio::EnsEMBL::Analysis::RunnableDB::Blast;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Blast);

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep
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
  my $pta = $self->db->get_PredictionTranscriptAdaptor;
  my $logic_names = $BLAST_AB_INITIO_LOGICNAME ;
  if ( !ref($logic_names) || scalar(@$logic_names) == 0 ) {
    $logic_names = ['Genscan'];
  }
  my @pts ;
  foreach my $logic_name (@$logic_names) {
    my $pt = $pta->fetch_all_by_Slice($self->query, $logic_name);
    push @pts, @$pt ;
  }
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }  

  # submit blast module to use via analysis_parameters column of analysis table 
  my $options_string ;  
  my %options = %{$self->PARSER_PARAMS};  

  if ( $blast{-type}=~m/ncbi/ ) { 
    if ( $options{-query_type}=~m/pep/ ) {  
      if ( $options{-database_type}=~m/pep/ ) { 
           $options_string = '-p blastp' ;  
      } elsif ( $options{-database_type}=~m/dna/ ) { 
         $options_string = '-p tblastn' ; 
      } 
    }   
  
    if ( $options{-query_type}=~m/dna/ ) {  
      if ( $options{-database_type}=~m/dna/ ) {  
           $options_string = '-p blastn' ;  
      }elsif ( $options{-database_type}=~m/pep/ ) {   
           $options_string = '-p blastx' ;  
      }
    }   
  }


  foreach my $t(@pts){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep->
      new(
          -transcript => $t,
          -query => $self->query,
          -program => $self->analysis->program_file,
          -parser => $parser,
          -filter => $filter,
          -database => $self->analysis->db_file,
          -analysis => $self->analysis,   
          -options => $options_string, 
          %blast,
         );
    $self->runnable($runnable);
  }
}




