=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::TRF - 

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::RunnableDB::TRF->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnable->fetch_input;
  $runnable->run;
  $runnable->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable TRF which wraps the program TRF

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the repeat_feature and repeat_consensus tables 

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRF;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::TRF;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::TRF
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   :

=cut



sub fetch_input{
  my ($self) = @_;

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  my $rfa = $dba->get_RepeatFeatureAdaptor();
  $self->get_adaptor($rfa);

  if($self->param('use_genome_flatfile')) {
    say "Ingoring dna table and using fasta file for sequence fetching";
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } elsif($self->param('dna_db')) {
    say "Attaching dna db to target";
    my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $dba->dnadb($dna_dba);
  } else {
    say "Assuming the target db has dna";
  }

#  if($self->param_is_defined('dna_db')) {
#    say "Attaching dna_db to output db adaptor";
#    my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
#    $dba->dnadb($dna_dba);
#  } else {
#    say "No dna_db param defined, so assuming target_db has dna";
#  }

  $self->hrdb_set_con($dba,'target_db');

  my $slice_array = $self->param('iid');
  unless(ref($slice_array) eq "ARRAY") {
    $self->throw("Expected an input id to be an array reference. Type found: ".ref($slice_array));
  }

  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                              -program_file => $self->param('trf_path'),
                                              -parameters => $self->param('commandline_params'),
                                            );

  $self->analysis($analysis);

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  foreach my $slice_name (@{$slice_array}) {
    my $slice = $self->fetch_sequence($slice_name,$dba);
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::TRF->new
    (
     -query     => $slice,
     -program   => $analysis->program_file,
     -analysis  => $analysis,
     %parameters,
    );
    $self->runnable($runnable);
  }

  return 1;
}


=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::TRF
  Arg [2]   : Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
  Function  : get/set repeatfeature adaptor
  Returntype: Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
  Exceptions: none
  Example   :

=cut

sub get_adaptor {
  my ($self,$rfa) = @_;
  if($rfa) {
    $self->param('_rfa',$rfa);
  }

  return($self->param('_rfa'));
}

1;
