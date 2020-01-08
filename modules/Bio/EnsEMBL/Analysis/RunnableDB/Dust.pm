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

package Bio::EnsEMBL::Analysis::RunnableDB::Dust;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::Dust;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Dust
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence;
  $self->query($slice);
  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Dust->new
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

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Dust
  Function  : get repeatfeature adaptor
  Returntype: Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
  Exceptions: none
  Example   : 

=cut


sub get_adaptor{
  my ($self) = @_;
  return $self->db->get_RepeatFeatureAdaptor;
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
  my $adaptor = $self->get_adaptor;
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
    #print "Have feature ".$feature."\n";
    $feature->analysis($self->analysis);
    $feature->slice($self->query) if(!$feature->slice);
    $self->feature_factory->validate($feature);
    eval{
      $adaptor->store($feature);
    };
    if ($@){
      throw("RunnableDB:store failed, failed to write ".$feature." to ".
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
    my $transformed = $rf->transform($self->query->coord_system->name,
                                     $self->query->coord_system->version);
    if(!$transformed){
      throw("Failed to transform ".$rf." ".$rf->start." ".
            $rf->end."  ".$rf->seq_region_name." skipping \n");
      
    }
    push(@converted, $transformed);
  }
  return \@converted;
}


1;
