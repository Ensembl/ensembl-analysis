=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

This base config should be used by all pipeline created for the ensembl-annotation pipelines

=head1 METHODS

  default_options: returns the default options from HiveGeneric_conf and it adds pipeline_db,
    dna_db and use_tokens. The inheriting class needs to specify; pipe_dbname, pipe_db_server,
    port, user, password, reference_dbname, reference_db_server, user_r, dna_dbname, dna_db_server

  lsf_resource_builder: returns the parameters string for LSF meadow_type

  get_config_settings: Populate and return a hashref with the parameters for analyses like Blast

  _master_config_settings: contains all possible parameters
  _add_keys: Add keys from one hash to another

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

=head2 default_options

 Description: It returns a hashref containing the default options for HiveGeneric_conf
              and three DB connection hash: pipeline_db, reference_db and dna_db
 Returntype : Hashref
 Exceptions : None


=cut

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },

#        At the moment, we want to use tokens
        use_tokens => 1,
        'pipeline_db' => {
            -dbname => $self->o('pipe_dbname'),
            -host   => $self->o('pipe_db_server'),
            -port   => $self->o('port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('hive_driver'),
        },

        'dna_db' => {
            -dbname => $self->o('dna_dbname'),
            -host   => $self->o('dna_db_server'),
            -port   => $self->o('port'),
            -user   => $self->o('user_r'),
        },
    };
}


=head2 lsf_resource_builder

 Arg [1]    : String $queue, name of the queue to submit to, default is 'normal'
 Arg [2]    : Integer $mem, memory required in MB
 Arg [3]    : Arrayref String, list of server you will connect to during the analysis
 Arg [4]    : Arrayref Integer, list of tokens value for each server, default is 10
 Arg [5]    : Integer $num_threads, number of cores, use only if you ask for multiple cores
 Arg [6]    : String $extra_requirements, any other parameters you want to give to LSF option -R
 Example    : '1GB' => { LSF => $self->lsf_resource_builder('normal', 1000, [$self->default_options->{'pipe_db_server'}])},
              '3GB_multithread' => { LSF => $self->lsf_resource_builder('long', 3000, [$self->default_options->{'pipe_db_server'}], undef, 3)},
 Description: It will return the LSF requirement parameters you require based on the queue, the memory, the database servers, the number
              of CPUs. It uses options -q, -n, -M and -R. If you need any other other options you will need to add it to the returned string.
              If you want multiple cores from multiple CPUs, you will have to edit the returned string.
              If a server appears more than once, it will use the highest token value
              The command below will create a resource string for a job in the long queue, asking for 3GB of memory and it will reserve 10 token
              for server1, 20 for server2 and 10 for server3.
              $self->lsf_resource_builder('long', 3000, ['server1', 'server2', 'server3'], [, 20])
 Returntype : String
 Exceptions : None


=cut

sub lsf_resource_builder {
    my ($self, $queue, $memory, $servers, $tokens, $threads, $extra_requirements) = @_;

    my $lsf_requirement = '-q '.($queue || 'normal');
    my @lsf_rusage;
    my @lsf_select;
    $extra_requirements = '' unless (defined $extra_requirements);
    if (defined $memory) {
        $lsf_requirement .= ' -M '.$memory;
        push(@lsf_rusage, 'mem='.$memory);
        push(@lsf_select, 'mem>'.$memory);
    }
    if ($self->default_options->{use_tokens} and defined $servers) {
        my $i = 0;
        my %seen;
        foreach my $server (@$servers) {
            if (! exists $seen{$server}) {
            	  if ($server =~ /ens-livemirror/) {
                    push(@lsf_rusage, 'myens_livemirrortok='.($tokens->[$i] || 10));
                    $seen{$server} = $i;
            	  } else {
                    my ($server_id) = $server =~ /(\d+)$/;
                    push(@lsf_rusage, 'myens_build'.$server_id.'tok='.($tokens->[$i] || 10));
                    $seen{$server} = $i;
            	  }
            }
            else {
                my ($token) = $lsf_rusage[$seen{$server}] =~ /tok=(\d+)/;
                if (defined $token and ($tokens->[$i] || 10) > $token) {
                    $token = $tokens->[$i];
                    $lsf_rusage[$seen{$server}] =~ s/tok=\d+/tok=$token/;
                }
            }
            $i++;
        }
    }
    if (defined $threads) {
        $lsf_requirement .= ' -n '.$threads;
        $extra_requirements .= ' span[hosts=1]';
    }
    return $lsf_requirement.' -R"select['.join(', ', @lsf_select).'] rusage['.join(', ', @lsf_rusage).'] '.$extra_requirements.'"';
}

=head2 _master_config_settings

  Arg [1]    : String $config_group, key of one of the hash
  Example    : $_master_config_settings->{HiveBlastGenscanPep};
  Description: Contains all possible variations of parameters for some analyses
  Returntype : Hash ref with parameters
  Exceptions : None


=cut

sub _master_config_settings {

  my ($self,$config_group) = @_;
  my $master_config_settings = {

  HiveBlast => {
    DEFAULT => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::BPliteWrapper',
      PARSER_PARAMS => {
        -regex => '^(\w+)',
        -query_type => undef,
        -database_type => undef,
      },
      BLAST_FILTER => undef,
      FILTER_PARAMS => {},
      BLAST_PARAMS => {
        -unknown_error_string => 'FAILED',
        -type => 'wu',
      }
    },

    HiveBlastGenscanPep => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
      PARSER_PARAMS => {
                         -regex => '^(\w+\W\d+)',
                         -query_type => 'pep',
                         -database_type => 'pep',
                         -threshold_type => 'PVALUE',
                         -threshold => 0.01,
                       },
      BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
      FILTER_PARAMS => {
                         -min_score => 200,
                         -prune => 1,
                       },
    },

    HiveBlastGenscanVertRNA => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
      PARSER_PARAMS => {
                         -regex => '^(\w+\W\d+)',
                         -query_type => 'pep',
                         -database_type => 'dna',
                         -threshold_type => 'PVALUE',
                         -threshold => 0.001,
                       },
      BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
      FILTER_PARAMS => {
                         -prune => 1,
                       },
    },

    HiveBlastGenscanUnigene => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
      PARSER_PARAMS => {
                         -regex => '\/ug\=([\w\.]+)',
                         -query_type => 'pep',
                         -database_type => 'dna',
                         -threshold_type => 'PVALUE',
                         -threshold => 0.001,
                       },
      BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
      FILTER_PARAMS => {
                         -prune => 1,
                       },
      },
    },
  };

  return($master_config_settings->{$config_group});

}

=head2 get_config_settings

  Arg [1]    : String config_group, name of a config found in _master_config_settings
  Arg [2]    : String config_logic_name
  Example    : config_settings => $self->get_config_settings($config_group, $config_logic_name),
  Description: This is a helper sub created to access parameters that historically were held in separate configs in the
               old pipeline. These are now stored in the _master_config_settings sub below this one. In the analyses hashes
               earlier in the config sets of these param can be requested and stored in the config_settings hash which
               is them passed in as a parameter to the analysis. The converted analysis modules have code to take the key
               value pairs from the config_settings hash and assign the values to the getter/setter sub associated with the
               key.
  Returntype : Hash ref containing the parameters
  Exceptions : None


=cut

sub get_config_settings {


   # Shift in the group name (a hash that has a collection of logic name hashes and a default hash)
   # Shift in the logic name of the specific analysis
   my $self = shift;
   my $config_group = shift;
   my $config_logic_name = shift;

   # And additional hash keys will be stored in here
   my @additional_configs = @_;

   # Return a ref to the master hash for the group using the group name
   my $config_group_hash = $self->_master_config_settings($config_group);
   unless(defined($config_group_hash)) {
     die "You have asked for a group name in _master_config_settings that doesn't exist. Group name:\n".$config_group;
   }
   # Final hash is the hash reference that gets returned. It is important to note that the keys added have
   # priority based on the call to this subroutine, with priority from left to right. Keys assigned to
   # $config_logic_name will have most priority, then keys in any additional hashes, then keys from the
   # default hash. A default hash key will never override a $config_logic_name key
   my $final_hash;

   # Add keys from the logic name hash
   my $config_logic_name_hash = $config_group_hash->{$config_logic_name};
   unless(defined($config_logic_name_hash)) {
     die "You have asked for a logic name hash that doesn't exist in the group you specified.\n".
         "Group name:\n".$config_group."\nLogic name:\n".$config_logic_name;
   }

   $final_hash = $self->_add_keys($config_logic_name_hash,$final_hash);

   # Add keys from any additional hashes passed in, keys that are already present will not be overriden
   foreach my $additional_hash (@additional_configs) {
     my $config_additional_hash = $config_group_hash->{$additional_hash};
     $final_hash = $self->_add_keys($config_additional_hash,$final_hash);
   }

   # Default is always loaded and has the lowest key value priority
   my $config_default_hash = $config_group_hash->{'Default'};
   $final_hash = $self->_add_keys($config_default_hash,$final_hash);

   return($final_hash);
}

=head2 _add_keys

  Arg [1]    : Hash ref $hash_to_add, hashref to add to the final hashref
  Arg [2]    : Hash ref $final_hash, final hashref that will be return by get_config_settings
  Example    : $elf->_add_keys($hash_to_add, $final_hash);
  Description: Add keys from one hashref to another hashref
  Returntype : Hash ref
  Exceptions : None


=cut

sub _add_keys {
  my ($self,$hash_to_add,$final_hash) = @_;

  foreach my $key (keys(%$hash_to_add)) {
    unless(exists($final_hash->{$key})) {
      $final_hash->{$key} = $hash_to_add->{$key};
    }
  }

  return($final_hash);
}

1;
