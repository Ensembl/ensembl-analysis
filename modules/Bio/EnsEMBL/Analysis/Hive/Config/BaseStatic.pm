=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 DESCRIPTION

This base config should never be called directly. You should create configuration file
which inherits from this module. You should not call this module or the inheriting module
but you should use get_analysis_settings from Bio::EnsEMBL::Analysis::Tools::Utils

=head1 METHODS

  get_config_settings: Populate and return a hashref with the parameters for analyses like Blast
  _master_config: contains all possible parameters

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw (throw);


=head2 new

 Arg [1]    : None
 Description: Generic method to create an object
 Returntype : Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic
 Exceptions : None

=cut

sub new {
  my ($class) = @_;
  my $self = bless({}, $class);

  return $self;
}


=head2 _master_config

 Arg [1]    : String
 Description: Hashes containing config for programs such as blast, exonerate.
              Each inheriting class has to implement this method to be able to
              get the right parameters
 Returntype : Hashref
 Exceptions : Throws if the method has not been implemented in the inheriting class

=cut

sub _master_config {
  my ($self, $key) = @_;

  throw('You need to implement this method. It should return a hashref');
}


=head2 get_config_settings

  Arg [1]    : String config_group, name of a config found in _master_config
  Arg [2]    : String config_logic_name
  Description: This is a helper sub created to access parameters that historically were held in separate configs in the
               old pipeline. These are now stored in the _master_config sub below this one. In the analyses hashes
               earlier in the config sets of these param can be requested and stored in the config_settings hash which
               is them passed in as a parameter to the analysis. The converted analysis modules have code to take the key
               value pairs from the config_settings hash and assign the values to the getter/setter sub associated with the
               key.
  Returntype : Hash ref containing the parameters
  Exceptions : None


=cut

sub get_config_settings {
  my ($self, $config_group, $additional_hash) = @_;

  # Return a ref to the master hash for the group using the group name
  my $config = $self->_master_config('default');
  if ($config_group and $config_group !~ /^#.*#$/) {
    my $config_group_hash = $self->_master_config($config_group);
    throw("You have asked for a group name in _master_config of ".ref($self)." that doesn't exist. Group name:\n".$config_group)
      unless($config_group_hash);

    # All values in the first hash will be overwritten by the values in the second hash.
    # Unshared key/values will be added
    _add_keys($config, $config_group_hash);
  }
  _add_keys($config, $additional_hash) if ($additional_hash);

  if ($config) {
    return $config;
  }
  else {
    throw('Something went wrong, $config does not exist for '.$config_group.' in '.ref($self));
  }
}


=head2 _add_keys

 Arg [1]    : Hashref the Hash to write to
 Arg [2]    : Hashref the hash which data will overwrite data from Arg[1]
 Description: Add key from Arg[2] to Arg[1] and overwrite if the key is present.
              It does a deep copy when the value is a hashref
              If you have hashref in an arrayref, the arrayref in Arg[1] will be
              overwritten by the arrayref in Arg[2]. Maybe your config should not
              be that complicated.
 Returntype : None
 Exceptions : Throws if the value for a same key has a different data type

=cut

sub _add_keys {
  my ($hash1, $hash2) = @_;

  foreach my $key (keys %$hash2) {
    if (exists $hash1->{$key}) {
      if (ref($hash2->{$key}) eq ref($hash2->{$key})) {
        if (ref($hash2->{$key}) eq 'HASH') {
          _add_keys($hash1->{$key}, $hash2->{$key});
        }
        else {
          $hash1->{$key} = $hash2->{$key};
        }
      }
      else {
        throw("Discrepencies in type for $key");
      }
    }
    else {
      $hash1->{$key} = $hash2->{$key};
    }
  }
}


=head2 get_array_config_settings

 Arg [1]    : String $config_group, key for a hash containing an arrayref
 Arg [2]    : Arrayref $additional_array (optional), optional data to add to the array
 Description: Method to get an arrayref from a Static.pm config
 Returntype : Arrayref
 Exceptions : Throws if Arg[1] cannot be found in the config

=cut

sub get_array_config_settings {
  my ($self, $config_group, $additional_array) = @_;

  my $config = $self->_master_config('default');
  if ($config_group and $config_group !~ /^#.*#$/) {
    my $config_group_array = $self->_master_config($config_group);
    throw("You have asked for a group name in _master_config ".ref($self)." that doesn't exist. Group name:\n".$config_group)
      unless($config_group_array);
    push(@$config, @$config_group_array);
  }
  push(@$config, @$additional_array) if ($additional_array);
  return $config;
}
1;
