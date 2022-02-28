=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Hive_primate_genebuilding_single_conf;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },

##########################################################################
#                                                                        #
# CHANGE STUFF HERE                                                      #
#                                                                        #
##########################################################################
        'pipeline_name'              => '',

        'user_r'                     => '',
        'user_w'                     => '',
        'password'                   => '',
        'port'                       => '',

        'pipe_db_name'                => '',
        'reference_db_name'           => '',
        'dna_db_name'                 => '',
        'genblast_output_db_name'     => '',
        'genewise_output_db_name'     => '',
        'exonerate_output_db_name'    => '',


        'pipe_db_server'             => '',
        'reference_db_server'        => '',
        'dna_db_server'              => '',
        'genblast_output_db_server'  => '',
        'genewise_output_db_server'  => '',
        'exonerate_output_db_server' => '',
        'killlist_db_server'         => '',
        'output_path'                => '',
        'genome_file'             => '',

        'repeat_masking_logic_names' => [],



##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

        'seqs_per_chunk'             => 25,
        'fastasplit_random_path'     => '/software/ensembl/bin/fastasplit_random',


        'genblast_path'              => 'genblast',
        'genblast_eval_primates'     => '1e-20',
        'genblast_eval_other'        => '1e-20',
        'genblast_cov_primates'      => '0.8',
        'genblast_cov_other'         => '0.5',

        'exonerate_path'             => '/software/ensembl/genebuild/usrlocalensemblbin/exonerate-0.9.0',
        'exonerate_pid_primates'     => '80',
        'exonerate_pid_other'        => '50',
        'exonerate_cov_primates'     => '80',
        'exonerate_cov_other'        => '50',


        'create_type' => 'clone',
        'driver' => 'mysql',
        'num_tokens' => 25,
        'user' => 'ensro',

        'pipeline_db' => {
            -dbname => $self->o('pipe_db_name'),
            -host   => $self->o('pipe_db_server'),
            -port   => $self->o('port'),
            -user   => $self->o('user_w'),
            -pass   => $self->o('password'),
            -driver => $self->o('driver'),
        },

        'reference_db' => {
                            -dbname => $self->o('reference_db_name'),
                            -host   => $self->o('reference_db_server'),
                            -port   => $self->o('port'),
                            -user   => $self->o('user_r'),
                          },

        'dna_db' => {
                      -dbname => $self->o('dna_db_name'),
                      -host   => $self->o('dna_db_server'),
                      -port   => $self->o('port'),
                      -user   => $self->o('user_r'),
                    },

        'genblast_output_db' => {
                           -dbname => $self->o('genblast_output_db_name'),
                           -host   => $self->o('genblast_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                         },

        'genewise_output_db' => {
                           -dbname => $self->o('genewise_output_db_name'),
                           -host   => $self->o('genewise_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                        },

        'exonerate_output_db' => {
                           -dbname => $self->o('exonerate_output_db_name'),
                           -host   => $self->o('exonerate_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                         },


        'killlist_db' => {
                           -dbname    => $self->o('killlist_db_name'),
                           -host      => $self->o('killlist_db_server'),
                           -port      => $self->o('port'),
                           -user      => $self->o('user_r'),
                         },
    };
}


sub resource_classes {
    my $self = shift;

  # Note that this builds resource requests off some of the variables at the top. Works off the idea
  # that the references all get put on one server and the pipe db is on another
  my $pipe_db_server = $self->default_options()->{'pipe_db_server'};
  my $dna_db_server = $self->default_options()->{'dna_db_server'};
  my $genblast_output_db_server = $self->default_options()->{'genblast_output_db_server'};
  my $exonerate_output_db_server = $self->default_options()->{'exonerate_output_db_server'};
  my $killlist_db_server = $self->default_options()->{'killlist_db_server'};

  my $pipe_db_server_number;
  my $dna_db_server_number;
  my $genblast_output_db_server_number;
  my $exonerate_output_db_server_number;
  my $killlist_db_server_number;


  my $num_tokens = $self->default_options()->{'num_tokens'};

  unless($pipe_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "pipe_db_server: ".$pipe_db_server;
  }

  $pipe_db_server_number = $1;

    unless($dna_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "dna_db_server: ".$dna_db_server;
  }

  $dna_db_server_number = $1;

  unless($genblast_output_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "genblast_output_db_server: ".$genblast_output_db_server;
  }

  $genblast_output_db_server_number = $1;

    unless($exonerate_output_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "exonerate_output_db_server: ".$exonerate_output_db_server;
  }

  $exonerate_output_db_server_number = $1;

    unless($killlist_db_server=~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "killlist_db_server: ".$killlist_db_server;
  }

  $killlist_db_server_number = $1;

  unless($num_tokens) {
    die "num_tokens is uninitialised or zero. num_tokens needs to be present in default_options and not zero\n".
        "num_tokens: ".$num_tokens;
  }

    return {
      'default' => { LSF => '-q normal -M900 -R"select[mem>900] rusage[mem=900,'.
                            'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'genblast' => { LSF => '-q normal -M1900 -R"select[mem>1900] rusage[mem=1900,'.
                             'myens_build'.$genblast_output_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'genblast_retry' => { LSF => '-q normal -M4900 -R"select[mem>4900] rusage[mem=4900,'.
                                   'myens_build'.$genblast_output_db_server_number.'tok='.$num_tokens.','.
                                   'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                   'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'exonerate' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,'.
                              'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'exonerate_retry' => { LSF => '-q normal -M4900 -R"select[mem>4900] rusage[mem=4900,'.
                                    'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    }
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

      {
        -logic_name => 'create_genblast_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('genblast_output_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -input_ids => [{}],
      },

            {
        -logic_name => 'create_exonerate_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('exonerate_output_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -input_ids => [{}],
      },

      {
        -logic_name => 'download_uniprot_human_pe12',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
             -parameters => {
                         file_name => 'human_pe12.fasta',
                         taxon_id => '9606',
                         dest_dir => $self->o('output_path'),
                         compress => 1,
                         pe_level => '1,2',
                       },

        -rc_name          => 'default',
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['clean_uniprot_human'],
                      },
      },

      {
        -logic_name => 'download_uniprot_primate',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         file_name => 'primate_pe12.fasta',
                         taxon_id => '9443',
                         exclude_group => 'human',
                         dest_dir => $self->o('output_path'),
                         compress => 1,
                         pe_level => '1,2',
                       },

        -rc_name          => 'default',
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['clean_uniprot_primates'],
                      },

      },

      {
        -logic_name => 'download_uniprot_primate_pe345',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         file_name => 'primate_pe345.fasta',
                         taxon_id => '9443',
                         exclude_group => 'human',
                         dest_dir => $self->o('output_path'),
                         compress => 1,
                         pe_level => '3,4,5',
                       },

        -rc_name          => 'default',
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['clean_uniprot_primates_345'],
                      },

      },

      {
        -logic_name => 'download_uniprot_mammal',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         file_name => 'mammal_pe12.fasta',
                         taxon_id => '40674',
                         exclude_group => 'primates',
                         dest_dir => $self->o('output_path'),
                         compress => 1,
                         pe_level => '1,2',
                       },

        -rc_name          => 'default',
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['clean_uniprot_mammals'],
                      },

      },

      {
        -logic_name => 'download_uniprot_vert',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         file_name => 'vert_pe12.fasta',
                         taxon_id => '7742',
                         exclude_group => 'mammalia',
                         dest_dir => $self->o('output_path'),
                         compress => 1,
                         pe_level => '1,2',
                       },

        -rc_name          => 'default',
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['clean_uniprot_vert'],
                      },
      },

      {
        -logic_name => 'clean_uniprot_human',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanFastaHeaders',
        -parameters => {
                         min_seq_length => '50',
                         skip_broken_headers => 0,
                         header_source => 'uniprot',
                         use_killlist => 1,
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                       },

        -flow_into => {
                        1 => ['make_chunk_human'],
                      },
        -rc_name          => 'default',
      },

      {
        -logic_name => 'clean_uniprot_primates',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanFastaHeaders',
        -parameters => {
                         min_seq_length => '50',
                         skip_broken_headers => 0,
                         header_source => 'uniprot',
                         use_killlist => 1,
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                       },

        -flow_into => {
                        1 => ['make_chunk_primates'],
                      },
        -rc_name          => 'default',
      },

           {
        -logic_name => 'clean_uniprot_primates_345',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanFastaHeaders',
        -parameters => {
                         min_seq_length => '50',
                         skip_broken_headers => 0,
                         header_source => 'uniprot',
                         use_killlist => 1,
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                       },

        -flow_into => {
                        1 => ['make_chunk_primates_345'],
                      },
        -rc_name          => 'default',
      },

     {
        -logic_name => 'clean_uniprot_mammals',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanFastaHeaders',
        -parameters => {
                         min_seq_length => '50',
                         skip_broken_headers => 0,
                         header_source => 'uniprot',
                         use_killlist => 1,
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                       },

        -flow_into => {
                        1 => ['make_chunk_mammals'],
                      },
        -rc_name          => 'default',
      },

     {
        -logic_name => 'clean_uniprot_vert',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanFastaHeaders',
        -parameters => {
                         min_seq_length => '50',
                         skip_broken_headers => 0,
                         header_source => 'uniprot',
                         use_killlist => 1,
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                       },

        -flow_into => {
                        1 => ['make_chunk_vert'],
                      },
        -rc_name          => 'default',
      },

      {
        -logic_name => 'make_chunk_human',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         chunk_output_dir => $self->o('output_path').'/human_chunks',
                         seqs_per_chunk => $self->o('seqs_per_chunk'),
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },

        -flow_into => {
                        1 => ['genblast_human','exonerate_human'],
                      },
        -rc_name          => 'default',

      },

     {
        -logic_name => 'make_chunk_primates',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         chunk_output_dir => $self->o('output_path').'/primate_chunks',
                         seqs_per_chunk => $self->o('seqs_per_chunk'),
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },

        -flow_into => {
                        1 => ['genblast_primates','exonerate_primates'],
                      },
        -rc_name          => 'default',

      },


     {
        -logic_name => 'make_chunk_primates_345',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         chunk_output_dir => $self->o('output_path').'/primate_345_chunks',
                         seqs_per_chunk => $self->o('seqs_per_chunk'),
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },

        -flow_into => {
                        1 => ['genblast_primates_345','exonerate_primates_345'],
                      },
        -rc_name          => 'default',

      },

      {
        -logic_name => 'make_chunk_mammals',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         chunk_output_dir => $self->o('output_path').'/mammal_chunks',
                         seqs_per_chunk => $self->o('seqs_per_chunk'),
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },

        -flow_into => {
                        1 => ['genblast_mammals','exonerate_mammals'],
                      },
        -rc_name          => 'default',
      },

      {
        -logic_name => 'make_chunk_vert',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         chunk_output_dir => $self->o('output_path').'/vert_chunks',
                         seqs_per_chunk => $self->o('seqs_per_chunk'),
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },

             -flow_into => {
                        1 => ['genblast_vert','exonerate_vert'],
                      },
        -rc_name          => 'default',

      },

      {
        -logic_name => 'genblast_human',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_human',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_primates').' -c '.$self->o('genblast_cov_primates').' ',
                         query_seq_dir => $self->o('output_path').'/human_chunks',
                       },
        -rc_name    => 'genblast',
        -wait_for => ['create_genblast_output_db'],
        -flow_into => {
                        -1 => ['rechunk_genblast_human'],
                        -2 => ['rechunk_genblast_human'],
                      },
      },


     {
        -logic_name => 'rechunk_genblast_human',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/human_chunks',
                         chunk_output_dir => $self->o('output_path').'/human_genblast_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['genblast_human_retry'],
                      },
      },

      {
        -logic_name => 'genblast_human_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_human',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_primates').' -c '.$self->o('genblast_cov_primates').' ',
                         query_seq_dir => $self->o('output_path').'/human_genblast_fail_chunks',
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'genblast_primates',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_primates',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_primates').' -c '.$self->o('genblast_cov_primates').' ',
                         query_seq_dir => $self->o('output_path').'/primate_chunks',
                       },
        -rc_name    => 'genblast',
        -wait_for => ['create_genblast_output_db'],
        -flow_into => {
                        -1 => ['rechunk_genblast_primates'],
                        -2 => ['rechunk_genblast_primates'],
                      },

      },

      {
        -logic_name => 'rechunk_genblast_primates',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/primate_chunks',
                         chunk_output_dir => $self->o('output_path').'/primate_genblast_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['genblast_primates_retry'],
                      },
      },

      {
        -logic_name => 'genblast_primates_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_primates',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_primates').' -c '.$self->o('genblast_cov_primates').' ',
                         query_seq_dir => $self->o('output_path').'/primates_genblast_fail_chunks',
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'genblast_primates_345',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_primates_345',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_primates').' -c '.$self->o('genblast_cov_primates').' ',
                         query_seq_dir => $self->o('output_path').'/primate_345_chunks',
                       },
        -rc_name    => 'genblast',
        -wait_for => ['create_genblast_output_db'],
        -flow_into => {
                        -1 => ['rechunk_genblast_primates_345'],
                        -2 => ['rechunk_genblast_primates_345'],
                      },
      },

      {
        -logic_name => 'rechunk_genblast_primates_345',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/primate_345_chunks',
                         chunk_output_dir => $self->o('output_path').'/primate_345_genblast_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['genblast_primates_345_retry'],
                      },
      },

      {
        -logic_name => 'genblast_primates_345_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_primates_345',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_primates').' -c '.$self->o('genblast_cov_primates').' ',
                         query_seq_dir => $self->o('output_path').'/primates_genblast_fail_chunks',
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'genblast_mammals',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_mammals',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_other').' -c '.$self->o('genblast_cov_other').' ',
                         query_seq_dir => $self->o('output_path').'/mammal_chunks',
                       },
        -rc_name    => 'genblast',
        -wait_for => ['create_genblast_output_db'],
        -flow_into => {
                        -1 => ['rechunk_genblast_mammals'],
                        -2 => ['rechunk_genblast_mammals'],
                      },

      },

      {
        -logic_name => 'rechunk_genblast_mammals',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/mammal_chunks',
                         chunk_output_dir => $self->o('output_path').'/mammal_genblast_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['genblast_mammals_retry'],
                      },
      },

      {
        -logic_name => 'genblast_mammals_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_mammals',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_other').' -c '.$self->o('genblast_cov_other').' ',
                         query_seq_dir => $self->o('output_path').'/mammal_genblast_fail_chunks',
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'genblast_vert',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_vert',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_other').' -c '.$self->o('genblast_cov_other').' ',
                         query_seq_dir => $self->o('output_path').'/vert_chunks',
                       },
        -rc_name    => 'genblast',
        -wait_for => ['create_genblast_output_db'],
        -flow_into => {
                        -1 => ['rechunk_genblast_vert'],
                        -2 => ['rechunk_genblast_vert'],
                      },
      },

      {
        -logic_name => 'rechunk_genblast_vert',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/vert_chunks',
                         chunk_output_dir => $self->o('output_path').'/vert_genblast_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['genblast_vert_retry'],
                      },
      },

      {
        -logic_name => 'genblast_vert_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast_vert',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval_other').' -c '.$self->o('genblast_cov_other').' ',
                         query_seq_dir => $self->o('output_path').'/vert_genblast_fail_chunks',
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'exonerate_human',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name          => 'exonerate',
        -wait_for => ['create_exonerate_output_db'],
        -parameters => {
                         logic_name => 'exonerate_human',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_human'),
                       },
        -flow_into => {
                        -1 => ['rechunk_exonerate_human'],
                        -2 => ['rechunk_exonerate_human'],
                      },

      },


      {
        -logic_name => 'rechunk_exonerate_human',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/human_chunks',
                         chunk_output_dir => $self->o('output_path').'/human_exonerate_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['exonerate_human_retry'],
                      },
      },

      {
        -logic_name => 'exonerate_human_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -parameters => {
                         logic_name => 'exonerate_human',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_human_retry','exonerate_human'),
                       },
        -rc_name          => 'exonerate_retry',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'exonerate_primates',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name          => 'exonerate',
        -wait_for => ['create_exonerate_output_db'],
        -parameters => {
                         logic_name => 'exonerate_primates',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_primates'),
                       },
        -flow_into => {
                        -1 => ['rechunk_exonerate_primates'],
                        -2 => ['rechunk_exonerate_primates'],
                      },
      },

     {
        -logic_name => 'rechunk_exonerate_primates',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/primate_chunks',
                         chunk_output_dir => $self->o('output_path').'/primate_exonerate_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['exonerate_primates_retry'],
                      },
      },

     {
        -logic_name => 'exonerate_primates_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -parameters => {
                         logic_name => 'exonerate_primates',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_primates_retry','exonerate_primates'),
                       },
        -rc_name          => 'exonerate_retry',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'exonerate_primates_345',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name          => 'exonerate',
        -wait_for => ['create_exonerate_output_db'],
        -parameters => {
                         logic_name => 'exonerate_primates_345',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_primates_345'),
                       },
        -flow_into => {
                        -1 => ['rechunk_exonerate_primates_345'],
                        -2 => ['rechunk_exonerate_primates_345'],
                      },

      },

     {
        -logic_name => 'rechunk_exonerate_primates_345',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/primate_345_chunks',
                         chunk_output_dir => $self->o('output_path').'/primate_345_exonerate_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['exonerate_primates_345_retry'],
                      },
      },

     {
        -logic_name => 'exonerate_primates_345_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -parameters => {
                         logic_name => 'exonerate_primates_345',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_primates_345_retry','exonerate_primates_345'),
                       },
        -rc_name          => 'exonerate_retry',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'exonerate_mammals',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name          => 'exonerate',
        -wait_for => ['create_exonerate_output_db'],
        -parameters => {
                         logic_name => 'exonerate_mammals',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_mammals'),
                       },
        -flow_into => {
                        -1 => ['rechunk_exonerate_mammals'],
                        -2 => ['rechunk_exonerate_mammals'],
                      },

      },

     {
        -logic_name => 'rechunk_exonerate_mammals',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/mammal_chunks',
                         chunk_output_dir => $self->o('output_path').'/mammal_exonerate_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['exonerate_mammals_retry'],
                      },
      },

   {
        -logic_name => 'exonerate_mammals_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -parameters => {
                         logic_name => 'exonerate_mammals',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_mammals_retry','exonerate_mammals'),
                       },
        -rc_name          => 'exonerate_retry',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'exonerate_vert',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name          => 'exonerate',
        -wait_for => ['create_exonerate_output_db'],
        -parameters => {
                         logic_name => 'exonerate_vert',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_vert'),
                       },
        -flow_into => {
                        -1 => ['rechunk_exonerate_vert'],
                        -2 => ['rechunk_exonerate_vert'],
                      },

      },

     {
        -logic_name => 'rechunk_exonerate_vert',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -can_be_empty  => 1,
        -parameters => {
                         rechunk => 1,
                         rechunk_dir_path => $self->o('output_path').'/vert_chunks',
                         chunk_output_dir => $self->o('output_path').'/vert_exonerate_fail_chunks',
                         seqs_per_chunk => 1,
                         chunk => 1,
                         fastasplit_random_path => $self->o('fastasplit_random_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        1 => ['exonerate_vert_retry'],
                      },
      },

      {
        -logic_name => 'exonerate_vert_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -parameters => {
                         logic_name => 'exonerate_vert',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate_vert_retry','exonerate_vert'),
                       },
        -rc_name          => 'exonerate_retry',
        -can_be_empty  => 1,
      },

    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    };
}

# override the default method, to force an automatic loading of the registry in all workers
#sub beekeeper_extra_cmdline_options {
#    my $self = shift;
#    return "-reg_conf ".$self->o("registry");
#}


sub get_config_settings {

   # Shift in the group name (a hash that has a collection of logic name hashes and a default hash)
   # Shift in the logic name of the specific analysis
   my $self = shift;
   my $config_group = shift;
   my $config_logic_name = shift;

   # And additional hash keys will be stored in here
   my @additional_configs = @_;

   # Return a ref to the master hash for the group using the group name
   my $config_group_hash = $self->master_config_settings($config_group);
   unless(defined($config_group_hash)) {
     die "You have asked for a group name in master_config_settings that doesn't exist. Group name:\n".$config_group;
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

   $final_hash = $self->add_keys($config_logic_name_hash,$final_hash);

   # Add keys from any additional hashes passed in, keys that are already present will not be overriden
   foreach my $additional_hash (@additional_configs) {
     my $config_additional_hash = $config_group_hash->{$additional_hash};
     $final_hash = $self->add_keys($config_additional_hash,$final_hash);
   }

   # Default is always loaded and has the lowest key value priority
   my $config_default_hash = $config_group_hash->{'Default'};
   $final_hash = $self->add_keys($config_default_hash,$final_hash);

   return($final_hash);
}

sub add_keys {
  my ($self,$hash_to_add,$final_hash) = @_;

  foreach my $key (keys(%$hash_to_add)) {
    unless(exists($final_hash->{$key})) {
      $final_hash->{$key} = $hash_to_add->{$key};
    }
  }

  return($final_hash);
}

sub master_config_settings {

  my ($self,$config_group) = @_;
  my $master_config_settings = {
  exonerate_protein => {
    Default => {
                 IIDREGEXP           => '(\d+):(\d+)',
                 OPTIONS             => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1',
                 COVERAGE_BY_ALIGNED => 0,
                 QUERYTYPE           => 'protein',
                 GENOMICSEQS         => $self->o('genome_file'),
                 OUTDB               => $self->o('exonerate_output_db'),
                 REFDB               => $self->o('reference_db'),
                 PROGRAM             => $self->o('exonerate_path'),
                 SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
               },

    exonerate_human => {
                         QUERYSEQS                     => $self->o('output_path').'/human_chunks',
                         FILTER                        => {
                           OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                           PARAMETERS                  => {
                             -coverage                 => $self->o('exonerate_cov_primates'),
                             -percent_id               => $self->o('exonerate_pid_primates'),
                             -best_in_genome           => 1,
                             -reject_processed_pseudos => 1,
                           },
                         },
                       },

    exonerate_human_retry => {
                               QUERYSEQS                     => $self->o('output_path').'/human_exonerate_fail_chunks',
                             },

    exonerate_primates => {
                         QUERYSEQS                     => $self->o('output_path').'/primate_chunks',
                         FILTER                        => {
                           OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                           PARAMETERS                  => {
                             -coverage                 => $self->o('exonerate_cov_primates'),
                             -percent_id               => $self->o('exonerate_pid_primates'),
                             -best_in_genome           => 1,
                             -reject_processed_pseudos => 1,
                           },
                         },
                       },

    exonerate_primates_retry => {
                               QUERYSEQS                     => $self->o('output_path').'/primate_exonerate_fail_chunks',
                             },

    exonerate_primates_345 => {
                         QUERYSEQS                     => $self->o('output_path').'/primate_345_chunks',
                         FILTER                        => {
                           OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                           PARAMETERS                  => {
                             -coverage                 => $self->o('exonerate_cov_primates'),
                             -percent_id               => $self->o('exonerate_pid_primates'),
                             -best_in_genome           => 1,
                             -reject_processed_pseudos => 1,
                           },
                         },
                       },

    exonerate_primates_345_retry => {
                                     QUERYSEQS                     => $self->o('output_path').'/primate_345_exonerate_fail_chunks',
                                   },

    exonerate_mammals => {
                         QUERYSEQS                     => $self->o('output_path').'/mammal_chunks',
                         FILTER                        => {
                           OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                           PARAMETERS                  => {
                             -coverage                 => $self->o('exonerate_cov_other'),
                             -percent_id               => $self->o('exonerate_pid_other'),
                             -best_in_genome           => 1,
                             -reject_processed_pseudos => 1,
                           },
                         },
                       },

    exonerate_mammals_retry => {
                                 QUERYSEQS                     => $self->o('output_path').'/mammal_exonerate_fail_chunks',
                               },

    exonerate_vert => {
                         QUERYSEQS                     => $self->o('output_path').'/vert_chunks',
                         FILTER                        => {
                           OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                           PARAMETERS                  => {
                             -coverage                 => $self->o('exonerate_cov_other'),
                             -percent_id               => $self->o('exonerate_pid_other'),
                             -best_in_genome           => 1,
                             -reject_processed_pseudos => 1,
                           },
                         },
                       },


   exonerate_vert_retry => {
                             QUERYSEQS                     => $self->o('output_path').'/vert_exonerate_fail_chunks',
                           },

    killlist_protein => {
                         KILLLISTDB          => $self->o('killlist_db'),
                         USE_KILL_LIST       => 1,
                         KILL_TYPE           => 'protein',
                         KILL_LIST_FILTER    => {
                                                  -only_mol_type        => 'protein',
                                                  -user_id              => undef,
                                                  -from_source_species  => undef,
                                                  -before_date          => undef,
                                                  -having_status        => undef,
                                                  -reasons              => [],
                                                  -for_analyses         => [],
                                                  -for_species          => [],
                                                  -for_external_db_ids  => [],
                                                },
                       },
    },
  };

  return($master_config_settings->{$config_group});

}

1;
