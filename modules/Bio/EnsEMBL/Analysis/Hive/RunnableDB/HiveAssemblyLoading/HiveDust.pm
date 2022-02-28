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

Bio::EnsEMBL::Analysis::RunnableDB::Dust -

=head1 SYNOPSIS

  my $dust = Bio::EnsEMBL::Analysis::RunnableDB::Dust->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $dust->fetch_input;
  $dust->run;
  $dust->write_output;


=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable Dust which wraps the program tcdust

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the repeat_feature and repeat_consensus tables 

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDust;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Dust
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  my $rfa = $dba->get_RepeatFeatureAdaptor;
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
                                              -program_file => $self->param('dust_path'),
                                              -parameters => $self->param('commandline_params'),
                                            );

  $self->analysis($analysis);

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  my $runnable_class = 'Bio::EnsEMBL::Analysis::Runnable::Dust';
  if ($self->param('dust_path') =~ /dustmasker/) {
    $runnable_class = 'Bio::EnsEMBL::Analysis::Runnable::DustMasker';
  }
  $self->require_module($runnable_class);

  foreach my $slice_name (@{$slice_array}) {
    my $slice = $self->fetch_sequence($slice_name,$dba);
    my $runnable = $runnable_class->new
                   (
                     -query => $slice,
                     -program => $analysis->program_file,
                     -analysis => $analysis,
                      %parameters,
                   );
    $self->runnable($runnable);
  }

  return 1;
}


=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : calculate which hits are mostly gap, and set analysis and 
  slice on each wanted feature and store it
  Returntype: 1
  Exceptions: none
  Example   :

=cut


#NOTE: convert features is an important method see its docs for more info

sub write_output{
  my ($self) = @_;
  my $rf_adaptor = $self->get_adaptor;
  my $analysis = $self->analysis();

  my @transformed;
  foreach my $feature (@{$self->output}) {
    if ( $feature->length > 50000 ) {
      my $transformed_features = $self->convert_feature($feature);
      push(@transformed, @$transformed_features) if($transformed_features);
    }else{
      push(@transformed, $feature);
    }
  }
  foreach my $feature (@transformed) {
    $feature->analysis($self->analysis);
    $feature->slice($self->query) if(!$feature->slice);
    $self->feature_factory->validate($feature);
    eval{
      $rf_adaptor->store($feature);
    };
    if ($@){
      $self->throw("RunnableDB:store failed, failed to write ".$feature." to ".
            "the database $@");
    }
  }

}


=head2 convert_feature

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Dust
  Arg [2]   : Bio::EnsEMBL::RepeatFeature
  Function  : This method takes a repeat feature projects it onto the 
  sequence level coord_system calculates what percentage of gaps make up
  the low complexity sequence. If the sequence is more than 25% low 
  complexity it is throw away. Only features which are more than 50k are
  checked for speed reasons. Once this has been calculated the features
  are converted by to the coordinate system they came from but they may now
  be in more pieces. The main reason for doing this was because features 
  which lied across long gaps were causing problems for the core api but
  it also means we dont store a lot of features which simply mask out Ns
  Returntype: Bio::EnsEMBL::RepeatFeature
  Exceptions: throws if it cant transform the feature
  Example   : 

=cut



sub convert_feature{
  my ($self, $rf) = @_;

  print "Converting ".$rf->start." ".$rf->end." ".
    $rf->slice->seq_region_name."\n";
  my $ff = $self->feature_factory;
  my $projections = $rf->project('seqlevel');
  my @converted;
  my $feature_length = $rf->length;
  my $projected_length = 0;
 PROJECT:foreach my $projection(@$projections){
    $projected_length += ($projection->from_end -
                          $projection->from_start) +1;
  }
  my $percentage = 100;
  if($projected_length != 0){
    $percentage = ($projected_length / $feature_length)*100;
  }
  if($percentage <= 75){
    return;
  }
 REPEAT:foreach my $projection(@$projections){
    my $start = 1;
    my $end = $projection->to_Slice->length;
    my $slice = $projection->to_Slice;
    my $rc = $ff->create_repeat_consensus('dust', 'dust', 'simple', 'N');
    my $rf = $ff->create_repeat_feature($start, $end, 0, 0, $start,
                                        $end, $rc, $slice->name,
                                        $slice);
    my $transformed = $rf->transform($rf->slice->coord_system->name,
                                     $rf->slice->coord_system->version);
    if(!$transformed){
      $self->throw("Failed to transform ".$rf." ".$rf->start." ".$rf->end."  ".$rf->seq_region_name." skipping \n");
    }
    push(@converted, $transformed);
  }
  return \@converted;
}


=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Dust
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
