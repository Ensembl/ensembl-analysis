=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

IntronCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more introns which has been clustered 

=head1 CONTACT

ba1@sanger.ac.uk

=head1 METHODS

=cut


# Let the code begin ...

package Bio::EnsEMBL::Analysis::Tools::Algorithms::IntronCluster;

use warnings ;
use strict;

use Bio::EnsEMBL::Intron;
use Bio::EnsEMBL::Utils::Exception qw(throw warning );

=head2 new

  Arg[1]      :
  Example     : my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::IntronCluster->new() ;
  Description : Create new IntronCluster
  Return type : Bio::EnsEMBL::Analysis::Tools::Algorithms::IntronCluster
  Exceptions  : If arg passed in

=cut

sub new {
  my ($class,$whatever)=@_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);

  if ($whatever){
    throw( "Can't pass an object to new() method. Use put_Introns() to include Bio::EnsEMBL::Intron in cluster");
  }

  $self->{_cached_start}  = undef;
  $self->{_cached_end}    = undef;
  $self->{_cached_strand} = undef;
  $self->{v} = 0 ; # verbosity 
  return $self;
}

=head2 put_Introns

  Arg[1]      : arrayref of Bio::EnsEMBL::Intron 
  Arg[2]      : Bio::EnsEMBL::Transcript
  Example     : $intron_cluster->put_Introns([$intron],$trans); 
  Description : add an intron to the cluster
  Return type : None
  Exceptions  : Type-sets not set

=cut

sub put_Introns {
  my ($self, $new_introns, $transcript, $ignore_strand)= @_;
  if ( !defined( $self->{'_types_sets'} ) ){
    throw( "Cluster lacks references to intron-types, unable to put the intron");
  }

 INTRON:
  foreach my $intron (@$new_introns){
    throw("undef for intron") if (!$intron);

    $self->_add_transcript_reference($intron,$transcript);

    my $intron_biotype = $transcript->biotype;

    foreach my $set_name ( keys %{$self->{'_types_sets'}}) {

      my $set = $self->{'_types_sets'}{$set_name};
      foreach my $type ( @{$set} ){

        if ($intron_biotype eq $type) {
          push ( @{ $self->{'_intron_sets'}{$set_name} }, $intron );

          if (defined($self->{_cached_start})) {
            if ($intron->start != $self->{_cached_start}) {
              throw("Failed putting intron: start");
            } 
          }
          if (defined($self->{_cached_end})) {
            if ($intron->end != $self->{_cached_end}) {
              throw("Failed putting intron: end");
            }
          }
          if (!$ignore_strand) {
            if (defined($self->{_cached_strand})) {
              if ($intron->strand != $self->{_cached_strand}) {
                throw("Failed putting intron: strand");
              }
            }
          }
          next INTRON; 
        }
      }
    }
    throw("Failed putting intron of type " . $intron_biotype . "\n");
  }
}

=head2 get_Introns

  Arg[1]      : None
  Example     : foreach my $intron (@{ $self->get_Introns} ) { 
  Description : Gets all introns in an intron cluster
  Return type : Bio::EnsEMBL::Intron

=cut

sub get_Introns {
  my $self = shift @_;

  my @introns;
  if (!defined( $self->{'_intron_sets'} ) ) {
    $self->warning("The intron array you try to retrieve is empty");
    @introns = ();
  }

  foreach my $set_name (keys %{$self->{'_intron_sets'}}) {
    push( @introns, @{ $self->{'_intron_sets'}{$set_name} } );
  }

  return \@introns;
}

=head2 strand

  Arg[1]      : None
  Example     : $strand = $intron_cluster->strand
  Description : Gets strand of intron cluster
  Return type : Int

=cut

sub strand{
  my $self = shift;

  if (!defined($self->{_cached_strand})) {
    my @introns = @{$self->get_Introns};
    unless (@introns){
      $self->warning("cannot retrieve the strand in a cluster with no introns");
    }
    my $strand;
    foreach my $intron (@introns){
      # looking at only 1 should be enough
      $strand = $intron->prev_Exon->strand if (!$strand);
      my $tmp_strand = $intron->prev_Exon->strand;
      if ($tmp_strand != $strand) {
        throw("introns not on the same strand");
      }
    }
    $self->{_cached_strand} = $strand;
  }
  return $self->{_cached_strand};
}

=head2 get_Introns_by_Set

  Arg[1]      : String (set name) 
  Example     : my @introns = $intron_cluster->get_Introns_by_Set($setname); 
  Description : Gets all introns in the intron cluster
                belonging to the set
  Return type : arrayref of Bio::EnsEMBL::Intron

=cut

sub get_Introns_by_Set() {
  my ($self,$set) = @_;

  unless ($set){
    throw( "must provide a set");
  }

  my @selected_introns;
  if ($self->{v}){
    for (keys %{ $self->{_intron_sets} } ) {
      print " i know the following sets : $_\n" ;
    }
  }
  if (!defined($self->{'_intron_sets'}{$set})) {
    # throw("No introns of set name $set");
    warning("No introns of set name $set in cluster");
  }else{
     push @selected_introns, @{$self->{'_intron_sets'}{$set}};
  }
  return \@selected_introns;
}                 

=head2 start

  Arg[1]      : None
  Example     : $start = $intron_cluster->start
  Description : Gets start position of intron cluster
                (The smallest number, regardless of strand)
  Return type : Int

=cut

sub start{
  my ($self) = @_;

  if (!defined($self->{_cached_start})) {
    my $start;

    foreach my $intron (@{$self->get_Introns}) {
      my $this_start = $intron->start;
      unless ( $start ){
        $start = $this_start;
      }
      if ( $this_start < $start ){
        $start = $this_start;
      }
    }
    $self->{_cached_start} = $start;
  }
  return $self->{_cached_start};
}
      
=head2 end

  Arg[1]      : None
  Example     : $end = $intron_cluster->end
  Description : Gets end position of intron cluster
                (The largest number, regardless of strand)
  Return type : Int

=cut

sub end{
  my ($self) = @_;

  if (!defined($self->{_cached_end})) {
    my $end;

    foreach my $intron (@{$self->get_Introns}) {
      my $this_end = $intron->end;
      unless ( $end ){
        $end = $this_end;
      }
      if ( $this_end > $end ){
        $end = $this_end;
      }
    }
    $self->{_cached_end} = $end;
  }
  return $self->{_cached_end};
}

=head2 get_transcripts_having_Intron_in_IntronCluster

  Arg[1]      : Bio::EnsEMBL::Intron
  Example     : @transcripts = @{$clust->get_transcripts_having_Intron_in_IntronCluster($intron)}; 
  Description : Gets all transcripts in the intron cluster
                that have this intron
  Return type : Arrayref of Bio::EnsEMBL::Transcripts

=cut

sub get_transcripts_having_Intron_in_IntronCluster {
  my ($self,$intron) = @_;

  my @transcript_array;
  my %transhash =  $self->each_transcripts_introns;

  TRANS: foreach my $trans_id (keys %transhash) {
    foreach my $intron_to_test (@{$transhash{$trans_id}{'introns'}}) {
      #print STDERR "transcript $trans_id, intron start ".$intron_to_test->start.", intron end ".$intron_to_test->end.", self start ".$self->start.", self end ".$self->end."\n";
      if($intron_to_test->start == $self->start &&  $intron_to_test->end == $self->end){
        push @transcript_array, $transhash{$trans_id}{'transcript'};
         next TRANS;
      }
    }
  }
  return \@transcript_array;
}

=head2 each_transcripts_introns

  Arg[1]      : None
  Example     : my %transhash =  $self->each_transcripts_introns;
  Description : Get a hash of introns keyed on transcript unique 
                refernce key (dbname_dbID)
  Return type : Hash

=cut

sub each_transcripts_introns {
  my $self = shift;

  return %{$self->{_transcripthash}};
}

=head2 _add_transcript_reference

  Arg[1]      : Bio::EnsEMBL::Intron
  Arg[2]      : Bio::EnsEMBL::Transcript
  Example     : $self->_add_transcript_reference($intron,$transcript);
  Description : Create an internal hash of introns keyed on transcript 
                unique reference key
  Return type : None

=cut

sub _add_transcript_reference {
  my ($self,$intron,$transcript) = @_;
  #if there's not already a reference to transcript stored make an arrayref
  if (!$self->contains_transcript($transcript)) {
    $self->{_transcripthash}{$transcript->adaptor->dbc->dbname."_".$transcript->dbID}{'introns'} = [];
    $self->{_transcripthash}{$transcript->adaptor->dbc->dbname."_".$transcript->dbID}{'transcript'} = $transcript;
  }
  # store introns of transcript (key: transcript)
  push @{$self->{_transcripthash}{$transcript->adaptor->dbc->dbname."_".$transcript->dbID}{'introns'}}, $intron;
}

=head2 contains_transcript

  Arg[1]      : Bio::EnsEMBL::Transcript
  Example     : if (!$self->contains_transcript($transcript)) {
  Description : Looks to see if this transcript exists
                in the internal transcript/intron hash
  Return type : Boolean

=cut
 
sub contains_transcript {
  my ($self,$transcript) = @_;
  return  (exists $self->{_transcripthash}{$transcript->adaptor->dbc->dbname."_".$transcript->dbID});
}

1;
