
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;


=head1 DESCRIPTION

fetches sequence data from database and instantiates and runs the
augustus runnable, this inherits from the Genscan runnableDB and as such doesnt
implement much itself

=head1 CONTACT

Post questions to : anacode-people@sanger.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAugustus;

use strict;
use warnings;
use feature 'say';
use Bio::Seq;
use Bio::EnsEMBL::Utils::IO::FASTASerializer;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Analysis::Runnable::Finished::AugustusGene;


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(empty_Transcript);

sub fetch_input{
  my ($self) = @_;

  my $repeat_masking = $self->param('repeat_masking_logic_names');
  my $soft_masking = $self->param('soft_matching');
  my $write_dir = $self->param('write_dir');

#  my $dba;
  my $dna_db = $self->get_database_by_name('dna_db');
  my $dba = $self->get_database_by_name('target_db', $dna_db);
  $self->hrdb_set_con($dba,'target_db');

  # if($self->param('use_genome_flatfile') ) {
  #   say "Ingoring dna table and using fasta file for sequence fetching";
  #   unless($self->param_required('genome_file') && -e $self->param('genome_file') ) {
  #     $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
  #   } 
  #   setup_fasta(
  #                -FASTA => $self->param_required('genome_file'),
  #              );

  # } 

  # elsif($self->param('dna_db')) {
  #   say "Attaching dna db to target";
  #   my $dna_db = $self->get_database_by_name('dna_db');
  #   $dba = $self->get_database_by_name('target_db', $dna_db);
  #   $self->hrdb_set_con($dba,'target_db');
  # } 

  # else {
  #   say "Assuming the target db has dna";
  # }
  

  my $input_id = $self->param('iid');
  print $input_id,"\n";
  exit 1;
  my $slice = $self->fetch_sequence($input_id, $dba, $repeat_masking, $soft_masking);
  $self->query($slice);

  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                              -program_file => $self->param('augustus_path'),
                                            );

  $self->param('analysis',$analysis);

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::AugustusGene->new
                 (
                    -query => $self->query,
                    -program => $analysis->program_file,
                    -analysis => $self->param('analysis'),
                    -species => $self->param('species'),
                    -workdir => $write_dir,
                    -database_adaptor => $dba,
                 );

  $self->runnable($runnable);
  return 1;
}

=head2 runnable_path

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus
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
  my ($self);
  return "Bio::EnsEMBL::Analysis::Runnable::Finished::AugustusGene";
}

sub run {
  my ($self) = @_;

  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    $self->output($runnable->output);
  }

  return 1;
}

sub write_output{
  my ($self) = @_;

  my $adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;

  foreach my $transcript (@{$self->output}) {

    my $analysis = $self->analysis;
    $transcript->analysis($analysis);

    my $gene = Bio::EnsEMBL::Gene->new();

    $gene->analysis($transcript->analysis);
    $gene->biotype($gene->analysis->logic_name);
    $gene->add_Transcript($transcript);

    $adaptor->store($gene);

   }
  return 1;
}

sub get_biotype {
  my ($self,$biotype_hash) = @_;
  if($biotype_hash) {
    $self->param('_biotype_hash',$biotype_hash);
  }
  return($self->param('_biotype_hash'));
}

1;
