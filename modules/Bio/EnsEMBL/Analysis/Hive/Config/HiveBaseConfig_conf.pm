=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

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
                my ($server_id) = $server =~ /(\d+)$/;
                push(@lsf_rusage, 'myens_build'.$server_id.'tok='.($tokens->[$i] || 10));
                $seen{$server} = $i;
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

1;
