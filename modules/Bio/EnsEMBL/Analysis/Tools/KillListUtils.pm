# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
package Bio::EnsEMBL::KillList::KillList;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::KillList::AnalysisLite;
use Bio::EnsEMBL::KillList::Comment;
use Bio::EnsEMBL::KillList::Filter;
use Bio::EnsEMBL::KillList::KillObject;
use Bio::EnsEMBL::KillList::Reason;
use Bio::EnsEMBL::KillList::Species;
use Bio::EnsEMBL::KillList::Team;
use Bio::EnsEMBL::KillList::User;
use Bio::EnsEMBL::KillList::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars qw (@ISA);
@ISA = qw();


=head2 new

  Arg[1]      :  
  Example     :
  Description :
  Return type : 
  Exceptions  :
  Caller      :
  Status      :

=cut
sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($type, $db_params, $filter_params) =
          rearrange([qw(TYPE
                        KILL_LIST_DB
                        FILTER_PARAMS
                        )],@args);

  if (!defined($type)) {throw " ERROR: need to set -type\n";}
  $self->type          ( $type );
  $self->db_params     ( $db_params ) if ( defined $db_params );
  $self->filter_params ( $filter_params ) if ( defined $filter_params );

  return $self; # success - we hope!
}

=head2 type

  Arg[1]      :
  Example     :
  Description :
  Return type :
  Exceptions  :
  Caller      :
  Status      :

=cut
sub type {
  my $self = shift;
  $self->{'type'} = shift if ( @_ );
  return $self->{'type'};
}

=head2 db_params

  Arg[1]      :
  Example     :
  Description :
  Return type :
  Exceptions  :
  Caller      :
  Status      :

=cut
sub db_params {
  my ($self, $db_params) = @_;
  if ( !$self->{'db_params'} ) {
    $self->{'db_params'} = {};
  }
  if ( $db_params ) {
    throw("Must pass KillList:db_params a hashref not a ".$db_params)
      unless(ref($db_params) eq 'HASH');
    $self->{'db_params'} = $db_params; 
  }
  return $self->{'db_params'};
}

=head2 filter_params

  Arg[1]      :
  Example     :
  Description :
  Return type :
  Exceptions  :
  Caller      :
  Status      :

=cut
sub filter_params {
  my ($self, $filter_params) = @_;
  if ( !$self->{'filter_params'} ) {
    $self->{'filter_params'} = {};
  }
  if ( $filter_params ) {
    throw("Must pass KillList:filter_params a hashref not a ".$filter_params)
      unless(ref($filter_params) eq 'HASH');
    $self->{'filter_params'} = $filter_params;
  }
  return $self->{'filter_params'};
}

=head2 read_and_check_config

  Arg[1]      :
  Example     :
  Description :
  Return type :
  Exceptions  :
  Caller      :
  Status      :

=cut
sub read_and_check_config {
  my ($self, $var_hash) = @_;

  parse_config($self, $var_hash, $self->type);
}

=head2 make_filter

  Arg[1]      :
  Example     :
  Description :
  Return type :
  Exceptions  :
  Caller      :
  Status      :

=cut
sub make_filter {
  my ($self, $hash) = @_;
  my %filter = %$hash;
  my $filter = $self->kill_list_filter->new(
                                        %filter
                                       );
  return $filter;
}

=head2 get_kill_list

  Arg[1]      :  $db, the kill_list database 
  Arg[2]      :  $filter_options - an array
                 (maybe better as a hash or Filter obj)  
                 specifying how to filter the kill_list.
  Example     :  my %kill_list = %{Bio::EnsEMBL::Analysis::Tools::KillListUtils::get_kill_list($kill_list_db, $filter_options)};
  Description :
  Return type :  Kill_list hash where keys = accession
                 and values = kill_objects 
  Exceptions  :
  Caller      :
  Status      :

=cut
sub get_kill_list {
  my ($self) = @_;

  #connect to the db
  my $db;
  if ($self->db_params) {
    if (!$self->db_params->dbname || !$self->db_params->dbhost ||
        !$self->db_params->dbuser || !$self->db_params->dbport ||
        !$self->db_params->dbpass ) {
      throw("Some db parameters are not set");
    }
    $db = Bio::EnsEMBL::KillList::DBSQL::DBAdaptor->new(
        '-dbname' => $self->db_params->dbname,
        '-host'   => $self->db_params->dbhost,
        '-user'   => $self->db_params->dbuser,
        '-port'   => $self->db_params->dbport,
        '-pass'   => $self->db_params->dbpass,
    );
  } else {
    throw("Must set db connection parameters in config");
  }

  #get the kill_list filter
  my $filter_params;
  if($self->filter_params){
    $filter_params = $self->filter_params;
  } else {
    throw("Must set filter parameters in config"); 
  }

  #now get all the filter parameters out of the hash
  my $user_id = $filter_params->{user_id} if (exists( $filter_params->{user_id} ));
  my $source_species = $filter_params->{source_species} if (exists( $filter_params->{source_species} ));
  my $date = $filter_params->{date} if (exists( $filter_params->{date} ));
  my $status = $filter_params->{status} if (exists( $filter_params->{status} ));
  my $mol_type = $filter_params->{mol_type} if (exists( $filter_params->{mol_type} ));
  my $reason_ids = $filter_params->{reason_ids} if (exists( $filter_params->{reason_ids} ));
  my $analysis_ids = $filter_params->{analysis_ids} if (exists( $filter_params->{analysis_ids} ));
  my $species_ids = $filter_params->{species_ids} if (exists( $filter_params->{species_ids} ));
  my $external_db_ids = $filter_params->{external_db_ids} if (exists( $filter_params->{external_db_ids} ));
  
  my (@reasons, @analyses, @species, @external_db_ids);

  #make necessary objects
  my $user = $db->get_UserAdaptor->fetch_by_dbID($user_id) if ($user_id);
  my $source_spp = $db->get_SpeciesAdaptor->fetch_by_dbID($source_species) if ($source_species);
  foreach my $id (@$reason_ids) {
    push @reasons, $db->get_ReasonAdaptor->fetch_by_dbID($id);
  }
  foreach my $id (@$analysis_ids) {
    push @analyses, $db->get_AnalysisLiteAdaptor->fetch_by_dbID($id);
  }
  foreach my $id (@$species_ids) {
    push @species, $db->get_SpeciesAdaptor->fetch_by_dbID($id);
  }

  my $filter = Bio::EnsEMBL::KillList::Filter->new(
               -user                    => $user,
               -from_source_species     => $source_spp,
               -before_date             => $date,
               -having_status           => $status,
               -only_mol_type           => $mol_type,
               -reasons                 => \@reasons,
               -for_analyses            => \@analyses,
               -for_species             => \@species,
               -for_external_db_ids     => $external_db_ids,
              );

  # an array of kill_objects
  my $kill_adaptor = $db->get_KillObjectAdaptor;
  my $kill_objects = $kill_adaptor->fetch_KillObjects_by_Filter($filter);

  my %kill_object_hash;
  foreach my $ko (@{$kill_objects}) {
    $kill_object_hash{$ko->accession} = $ko;
    print STDERR "got ".$ko->dbID." accesssion ".$ko->accession."\n";
  }
  return \%kill_object_hash;
}

=head2 update_meta_table

  Arg[1]      :  $db, the reference database to be updated
  Example     :
  Description :  Deletes any previous similar entries in the
                 meta table and adds a new entry, with 
                 meta_key = 'kill_list' and meta_value = now().
  Return type :  None
  Exceptions  :  None
  Caller      :
  Status      :

=cut
sub update_meta_table {
  my ($db) = @_;

  my $sth = $db->dbc->prepare(
            "DELETE FROM meta ".
            "WHERE meta_key = 'kill_list'");
  $sth->execute;
  $sth->finish;

  $sth = $db->dbc->prepare(
            "INSERT INTO meta ".
            "(meta_key, meta_value) ".
            "VALUES ('kill_list', now())");
  $sth->execute;
  $sth->finish;

  return;
}

1;
