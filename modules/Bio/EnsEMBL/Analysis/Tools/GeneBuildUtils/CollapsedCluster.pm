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
=pod 

=head1 NAME Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster

=head1 SYNOPSIS

my $collapsed_cluster = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster->new();

some examples:

add exons

$collapsed_cluster->add_exon($exon);

make and store introns

my $intron = $collapsed_cluster->make_intron_from_exons($exon1,$exon2);

$intron->ev_set($exon1->ev_set);

$intron->biotype('intron');

$collapsed_cluster->add_intron($intron);

count the number of identical exons that have been collapsed

$num = $collapsed_cluster->get_exon_count($exon);

given an intron what are the exons that join to it?

@exons = @{$collapsed_cluster->next_exons_by_intron($intron)};

=head1 DESCRIPTION

This module collapses down a redundant set of exons or introns contained within
Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended objects into a non redundant
set. It stores supporting features and provides a number of convenience methods for
assessing the cluster.

=head1 CONTACT

Post questions to the EnsEMBL developer list: <http://lists.ensembl.org/mailman/listinfo/dev>

=cut

package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Feature ;
use Bio::EnsEMBL::Exon ;
use Bio::EnsEMBL::Utils::Exception qw(throw warning );

@ISA = qw(Bio::EnsEMBL::Exon Bio::EnsEMBL::Feature);

=head2 new

   Name      : new
   Arg[0]    : none
   Function  : Cretes the ColapsedCluster object
   Returnval : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster

=cut

sub new  {
  my($class) = shift;
  if( ref $class ) {
      $class = ref $class;
  }
  my $self = $class->SUPER::new(@_);

  $self->{_exonhash} = undef;
  $self->{_endexonhash} = undef;
  $self->{_intronhash} = undef;
  $self->{_exonstarts} = undef;
  $self->{_intronstarts} = undef;
  $self->{_exoncluster} = undef;
  $self->{_evidencecount} = undef;
  $self->{_featurecount} = undef;
  $self->{_score} = undef;
  $self->{_endexon} = undef;
  $self->{_exontranscript} = undef;

  return $self ;
}

=head2 add_intron

   Name       : add_intron
   Arg[0]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended
   Function   : stores a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object representing an intron
   Exceptions : Throws unless it recieves a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut


sub add_intron{
  my ($self,$feature) = @_;
  if (defined($feature)){
    $self->throw("Intron must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    $self->_add_feature($feature,'intron');
  }
  return;
}

=head2 add_exon

   Name       : add_exon
   Arg[0]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended
   Function   : stores a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object representing an exon
   Exceptions : Throws unless it recieves a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut

sub add_exon{
  my ($self,$feature,$cluster) = @_;
  if (defined($feature)) {
    $self->throw("exon must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    $self->_add_feature($feature,'exon');
  }
  if (defined($cluster)){
    $self->_add_exoncluster($feature,$cluster);
  }
  return;
}

=head2 _add_exoncluster

   Name       : _add_exoncluster
   Arg[0]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster
   Function   : Private method, stores an exon cluster object by association with an exon
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] is a Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster
   Returnval  : none

=cut

sub _add_exoncluster{
  my ($self,$feature,$cluster) = @_;
  if (defined($feature)) {
    $self->throw("exon must be a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    if (defined $cluster){
      $self->throw("Cluster must be a Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster not a ".ref($cluster)."\n")
	unless $cluster->isa("Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster");
    }
    my $key = $self->_get_key_from_feature($feature);
    $self->{_exoncluster}{$key} = $cluster;
  }
  return;
}

=head2 add_feature

   Name       : add_feature
   Arg[0]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : (optional) string
   Function   : Stores an exon or intron feature in a hash keyed on its genomic coordinates
              : If an identical feature has already been stored it transfers the supporting features
              : Stores how many identical features have been stored so far.
              : Stores end exons in a seperate hash to aid idenifying the longest end exon in the stack
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] is 'exon' or 'intron'
   Returnval  : none

=cut

sub _add_feature{
  my ($self,$feature,$type) = @_;
  unless ($type eq 'exon' or $type eq 'intron'){
    $self->throw("Feature must be defined as in intron or an exon");
  }
  if (defined($feature)) {
    my $key = $self->_get_key_from_feature($feature);
    # store associated features along with the feature
    $self->add_ev_by_feature($feature);
    # exons are more coplictated because of supporting features
    # and end exons
    if ($type eq 'exon'){
      # single exons are terminal but not 3 or 5 prime
      $self->add_exon_count($feature);
      if ($feature->is_terminal_exon && $feature->number_exons > 1 ){
	my $terminal_key = $self->_get_key_from_terminal_feature($feature);
	$self->{'_endexoncount'}{$terminal_key}++;
	$self->add_ev_by_terminal_feature($feature);	
	# only store the terminal exon if it is the longest
	# otherwise just transfer the supporting features
	if (my $existing_end_exon =  $self->get_end_exon($feature)){
	  # push supporting features onto the feature if it exists
	  if ($existing_end_exon->length < $feature->length){
	    $self->{'_endexonhash'}{$terminal_key} = $feature;
	  }
	} else {
	  # end exon hasnt been added before 
	  $self->{'_endexonhash'}{$terminal_key} = $feature;
	}
      }
      if (my $existing_exon =  $self->get_exon($feature)){
	my @sfs = @{$existing_exon->get_all_supporting_features};
	$feature->add_supporting_features(@sfs);
	$self->{'_exonhash'}{$key} = $feature;
      } else {
	# its a new exon add it for the first time
	$self->{'_exonhash'}{$key} = $feature;	
      }
      # no point in adding get exon by start if its the first exon
      $self->add_exon_by_start($feature) unless $key =~ /start/;
      $self->add_exon_by_end($feature) unless $key =~ /end/;
      # want to store *all* the transcripts that share this exon
      my $transcript = $feature->transcript;
      $self->transcript_by_feature($feature,$transcript);
    }
    # introns are nice and easy
    if ($type eq 'intron'){
      $self->add_intron_count($feature);
      $self->{'_intronhash'}{$key} = $feature;
    }
  }
  return ;
}


=head2 _get_key_from_feature

   Name       : _get_key_from_feature
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Private method for creating the hash key from an objects genomic cordinates
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws if no key is produced
   Returnval  : string $key

=cut

sub _get_key_from_feature{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    if ($feature->start && $feature->end && $feature->strand && $feature->seq_region_name){
      my $key = $feature->seq_region_name.':'.$feature->start.':'.$feature->end.':'.$feature->strand;
      $self->throw("Key not found for $feature\n") unless $key;
      return $key;
    }
  }
  return;
}

=head2 _get_key_from_terminal_feature

   Name       : _get_key_from_terminal_feature
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Private method for creating the hash key from an end exons genomic cordinates
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] 
   Returnval  : string $key

=cut

sub _get_key_from_terminal_feature{
  my ($self,$feature) = @_;
  if (defined($feature)) {    
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key;
    return 0 unless $feature->is_terminal_exon;
    if (($feature->is_5prim_exon && $feature->strand == 1) or
	($feature->is_3prim_exon && $feature->strand == -1) ){
      $key =  $feature->seq_region_name.':start:'.$feature->end.':'.$feature->strand;
    }
    if (($feature->is_3prim_exon && $feature->strand == 1) or
	($feature->is_5prim_exon && $feature->strand == -1) ){
      $key =  $feature->seq_region_name.':'.$feature->start.':end:'.$feature->strand;
    }
    return $key;
  }
  return;
}

=head2 add_ev_by_terminal_feature

   Name       : add_ev_by_terminal_feature
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Stores the number of est and similarity exons that were collapsed to make this end exon
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] 
   Returnval  : none

=cut

sub add_ev_by_terminal_feature{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_terminal_feature($feature);
    my $evidence = $feature->ev_set;
    $self->throw("No evidence set found for $feature\n") unless $evidence;
    $self->{_terminalevidencecount}{$key}{$evidence}++;
  }
  return;
}

=head2 ev_count_by_terminal_feature

   Name       : ev_count_by_terminal_feature
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : String name of evidence set 'est' or 'simgw'
   Function   : Returns the number of est and similarity exons that were collapsed to make this end exon
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Scalar 

=cut

sub ev_count_by_terminal_feature{
  my ($self,$feature,$ev_set) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_terminal_feature($feature);
    if (defined($ev_set)) {   
      if ($self->{_terminalevidencecount}{$key}{$ev_set}){
	return $self->{_terminalevidencecount}{$key}{$ev_set} ;
      } else {
	return 0;
      }
    }
  }
  return 0;
}

=head2 add_ev_by_feature

   Name       : add_ev_by_feature
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : String name of evidence set 'est' or 'simgw'
   Function   : Stores the number of est and similarity exons that were collapsed to make this feature
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] 
   Returnval  : none

=cut

sub add_ev_by_feature{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    my $evidence = $feature->ev_set;
    $self->throw("No evidence set found for $feature\n") unless $evidence;
    $self->{_evidencecount}{$key}{$evidence}++;
  }
  return;
}

=head2 ev_count_by_feature

   Name       : ev_count_by_feature
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : String name of evidence set 'est' or 'simgw'
   Function   : Returns the number of est and similarity exons that were collapsed to make this feature
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] 
   Returnval  : scalar

=cut

sub ev_count_by_feature{
  my ($self,$feature,$ev_set) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
   my $key = $self->_get_key_from_feature($feature);
    if (defined($ev_set)) {
      if ($self->{_evidencecount}{$key}{$ev_set}){
	return $self->{_evidencecount}{$key}{$ev_set};
      } else {
	return 0;
      }
    }
  }
  return 0;
}

=head2 add_exon_by_end

   Name       : add_exon_by_end
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Stores exons by end position
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut

sub add_exon_by_end{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    my $start = $feature->end.":".$feature->strand;
    $self->{_exonends}{$start}{$key} = $feature;
  }
  return;
}

=head2 add_exon_by_start

   Name       : add_exon_by_start
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : (optional) 
   Function   : Stores exons by start position
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut

sub add_exon_by_start{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    my $start = $feature->start.":".$feature->strand;
    $self->{_exonstarts}{$start}{$key} = $feature;
  }
  return;
}

=head2 last_exons_by_intron

   Name       : last_exons_by_intron
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Given an intron object it returns an array ref of exons that join onto the start of the intron
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Array ref of Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended objects 

=cut

sub last_exons_by_intron{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $start = $feature->start.":".$feature->strand;
    my @exons = values %{$self->{_exonends}{$start}};
    return \@exons;
  }
  return;
}

=head2 next_exons_by_intron

   Name       : next_exons_by_intron
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Given an intron object it returns an array ref of exons that join onto the end of the intron
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Array ref of Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended objects 

=cut

sub next_exons_by_intron{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $start = $feature->end.":".$feature->strand;
    my @exons = values %{$self->{_exonstarts}{$start}};
    return \@exons;
  }
  return;
}

=head2 exon_score

   Name       : exon_score
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : (optional) scalar score
   Function   : Get / set a score to assign to an exon
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Scalar score

=cut

sub exon_score{
  my ($self,$feature,$score) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    if  (defined($score)) {
      $self->{_exonscore}{$key} = $score;
    }
    return $self->{_exonscore}{$key};
  }
  return;
}

=head2 intron_score

   Name       : intron_score
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : (optional) scalar score
   Function   : Get / set a score to assign to an intron
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Scalar score

=cut

sub intron_score{
  my ($self,$feature,$score) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    if  (defined($score)) {
      $self->{_intronscore}{$key} = $score;
    }
    return $self->{_intronscore}{$key};
  }
  return;
}

=head2 add_intron_count

   Name       : add_intron_count
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Adds up the number of introns that were collapsed to make this feature
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut

sub add_intron_count{
  my ($self,$feature) = @_;
  if (defined($feature)) {    
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    $self->{_introncount}{$key}++;
  }
  return;
}


=head2 add_intron_weight

   Name       : add_intron_weight
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Adds extra weight to RNAseq introns that represent a collapsed stack of introns
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut

sub add_intron_weight{
  my ($self,$feature,$num) = @_;
  # we will already have a weight of 1 for this feature
  $num--;
  if (defined($feature)) {    
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    $self->{_introncount}{$key}+= $num;
  }
  return;
}

=head2 get_intron_count

   Name       : get_intron_count
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Returns the number of introns that were collapsed to make this feature
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut

sub get_intron_count{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    return $self->{_introncount}{$key};
  }
  return;
}

=head2 add_exon_count

   Name       : add_exon_count
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Adds up the number of exons that were collapsed to make this feature
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut

sub add_exon_count{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    $self->{_exoncount}{$key}++;
  }
  return;
}

=head2 get_intron_count

   Name       : get_exon_count
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Returns the number of exons that were collapsed to make this feature
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : none

=cut

sub get_exon_count{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    return $self->{_exoncount}{$key};
  }
  return;
}

=head2 get_exon

   Name       : get_exon
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Takes an exon object and returns the associated collapsed exon with all associated features
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended 

=cut

sub get_exon{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    return $self->{_exonhash}{$key};
  }
  return;
}

=head2 make_intron_from_exons

   Name       : make_intron_from_exons
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Given 2 exons it makes an intron object that joins them together
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended

=cut

sub make_intron_from_exons{
  my ($self,$exon1,$exon2) = @_;
  if (defined($exon1) && defined($exon2)) {
    $self->throw("exon1 must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($exon1)."\n")
      unless $exon1->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    $self->throw("exon2 must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($exon2)."\n")
      unless $exon2->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my @array = ( $exon1, $exon2 );
    # sort them my start
    @array = sort { $a->start <=> $b->start } @array;
    my $intron = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended->new
      (
       -start   => $array[0]->end,
       -end     => $array[1]->start,
       -strand  => $exon1->strand,
       -slice   => $exon1->slice,
       -analysis=> $exon1->transcript->analysis,
      );
    return $intron;
  }
  return;
}

=head2 make_intron_from_rnaseq

   Name       : make_intron_from_rnaseq
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Given 2 rnaseq it makes an intron object that joins them together
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended

=cut

sub make_intron_from_rnaseq{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature");
    my $intron = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended->new
      (
       -start   => $feature->start -1,
       -end     => $feature->end +1,
       -strand  => $feature->strand,
       -slice   => $feature->slice,
       -analysis=> $feature->analysis,
      );
    return $intron;
  }
  return;
}

=head2 get_end_exon

   Name       : get_end_exon
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Given an exon object it returns the longest end exon that shares the
              : internal exon boundary
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended

=cut

sub get_end_exon{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    return 0 unless $feature->is_terminal_exon;
    my $key = $self->_get_key_from_terminal_feature($feature);
    return $self->{_endexonhash}{$key};
  }
  return;
}

=head2 count_end_exon

   Name       : count_end_exon
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Given an exon object it returns the number of end exons that share the
              : internal exon boundary
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Scalar exon count

=cut

sub count_end_exon{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    return 0 unless $feature->is_terminal_exon;
    my $key = $self->_get_key_from_terminal_feature($feature);
    return $self->{_endexoncount}{$key};
  }
  return;
}

=head2 get_all_exons

   Name       : get_all_exons
   Function   : returns a non redundant set of collapsed exons sorted by start position
   Returnval  : Array ref Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended

=cut

sub get_all_exons{ 
  my ($self) = @_;
  my @exons = values %{$self->{'_exonhash'}};
  @exons = sort {$a->start <=> $b->start} @exons;
  return \@exons;
}

=head2 get_all_introns

   Name       : get_all_introns
   Function   : returns a non redundant set of collapsed introns sorted by start position
   Returnval  : Array ref Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended

=cut

sub get_all_introns{ 
  my ($self) = @_;
  my @introns =values %{$self->{_intronhash}};
  @introns  = sort {$a->start <=> $b->start} @introns;
  return \@introns;
}

=head2 contains_exon

   Name       : contains_exon
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Boolean test to see if the collapsed cluster contains the exon in question
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Scalar 1 or 0

=cut

sub contains_exon{ 
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    if ( $self->{_exonhash}{$key} ){
      return 1;
    }
    else {
      return 0;
    }
  }
  return;
}

=head2 contains_intron

   Name       : contains_intron
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Boolean test to see if the collapsed cluster contains the intron in question
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Scalar 1 or 0

=cut

sub contains_intron{ 
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    if ( $self->{_intronhash}{$key} ){
      return 1;
    }
    else {
      return 0;
    }
  }
  return;
}

=head2 transcript_by_feature

   Name       : transcript_by_feature
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : (optional) Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended
   Function   : Get / Set  a transcript object by reference to an exon or intron
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
              : Throws unless ARG[1] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended
   Returnval  : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended

=cut

sub transcript_by_feature{
  my ($self,$feature,$trans) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    if  (defined($trans)) {
    $self->throw("Transcript must be a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended not a ".ref($trans)."\n")
      unless $trans->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended");
      push (@{$self->{_exontranscript}{$key}}, $trans);
    }
    return $self->{_exontranscript}{$key};
  }
  return;
}

=head2 exon_info

   Name       : exon_info
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Prints out the key used to store the exon
   Returnval  : String key

=cut

sub exon_info {
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    print "$key";
  }
}

=head2 get_exon_cluster

   Name       : get_exon_cluster
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Function   : Retrieves the exon cluster object associated with the feature provided
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster

=cut

sub get_exon_cluster{
  my ($self,$feature) = @_;
  if (defined($feature)) {
    $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($feature)."\n")
      unless $feature->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
    my $key = $self->_get_key_from_feature($feature);
    return $self->{_exoncluster}{$key};
  }
  return;
}


1;

