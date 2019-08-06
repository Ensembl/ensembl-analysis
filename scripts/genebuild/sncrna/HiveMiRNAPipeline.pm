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

=cut

package HiveSncRNA;


use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
        %{$self->SUPER::default_options()},

        ######################################################
        #
        # Variable settings- You change these!!!
        #
        ######################################################
        ########################
        # Misc setup info
        ########################
        'dbowner'                                => '' || $ENV{EHIVE_USER} || $ENV{USER},
        'pipeline_name'                          => '',
        'user_r'                                 => 'ensro', # read only db user
        'user'                                   => 'ensadmin', # write db user
        'password'                               => '', # password for write db user
        'pipe_db_server'                         => '', # host for pipe db
        'databases_server'                       => '', # host for general output dbs
        'dna_db_server'                          => '', # host for dna db
        'pipe_db_port'                           => '', # port for pipeline host
        'databases_port'                         => '', # port for general output db host
        'dna_db_port'                            => '', # port for dna db host

        'release_number'                         => '' || $self->o('ensembl_release'),
        'species_name'                           => '', # e.g. mus_musculus
        'production_name'                        => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
        'clade'                                  => '', # eg. rodents
        'uniprot_set'                            => 'mammals_basic',
        'output_path'                            => '',


        ########################
        ## Small ncRNAs params
        #########################
        'base_blast_db_path'                     => $ENV{BLASTDB_DIR},
        'mirBase_fasta'                          => 'all_mirnas.fa', # What miRBase file to use, eg. /hps/nobackup2/production/ensembl/genebuild/blastdb/ncrna/mirbase_22/all_mirnas.{embl,fa,fa.nhr,fa.min,fa.nsq}
        'mature_mirnas'                          => 'mature.fa',
        'rfc_scaler'                             => 'filter_dafs_rfc_scaler_human.pkl',
        'rfc_model'                              => 'filter_dafs_rfc_model_human.pkl',

        'full_repbase_logic_name'                => 'repeatmask_repbase_' . $self->o('clade'),
        'blastn_exe_path'                        => catfile($self->o('binary_base'), 'blastn'),
        
        ########################
        # Pipe and ref db info
        ########################
        'pipe_db_name'                           => $self->o('dbowner') . '_' . $self->o('pipeline_name')  . '_pipe_' . $self->o('release_number'),
        'dna_dbname'                             => '',
        'dna_db_name'                            => $self->o('dna_dbname'),
        'ncrna_db_server'                        => $self->o('databases_server'),
        'ncrna_db_port'                          => $self->o('databases_port'),
        'ncrna_db_name'                          => $self->o('dbowner') . '_' . $self->o('species_name') . "_" .  $self->o('pipeline_name') . '_ncrna_' . $self->o('release_number'),

        'ensembl_release'                        => $ENV{ENSEMBL_RELEASE},
        'production_db_server'                   => 'mysql-ens-meta-prod-1',
        'production_db_port'                     => '4483',
        # BLAST db paths
        ########################
        'base_blast_db_path'                     => $ENV{BLASTDB_DIR},
        'ncrna_blast_path'                       => catfile($self->o('base_blast_db_path'), 'ncrna', 'ncrna_2016_05'),
        'mirna_blast_path'                       => catfile($self->o('base_blast_db_path'), 'ncrna', 'mirbase_22'),

        ######################################################
        #
        # Mostly constant settings
        #
        ######################################################

        genome_dumps                             => catdir($self->o('output_path'), 'genome_dumps'),
        # This one is used by most analyses that run against a genome flatfile like exonerate, genblast etc. Has slice name style headers. Is softmasked
        faidx_genome_file             => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),
        ncrna_dir                                => catdir($self->o('output_path'), 'ncrna'),

        ensembl_analysis_script                  => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),

        sequence_dump_script                     => catfile($self->o('ensembl_analysis_script'), 'sequence_dump.pl'),
        sncrna_analysis_script                   => catdir($self->o('ensembl_analysis_script'), 'genebuild', 'sncrna'),

        num_tokens                               => 10,
        mysql_dump_options                       => '--max_allowed_packet=400MB',

        ########################
        # Executable paths
        ########################
        blast_type                               => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'

        trnascan_path                            => catfile($self->o('binary_base'), 'tRNAscan-SE'),
        cmsearch_exe_path                        => catfile($self->o('software_base_path'), 'bin', 'cmsearch'), # #'opt', 'infernal10', 'bin', 'cmsearch'),

        samtools_path                            => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary

        bedtools                                 => catfile($self->o('binary_base'), 'bedtools'),
        bedGraphToBigWig                         => catfile($self->o('binary_base'), 'bedGraphToBigWig'),

        ########################
        # db info
        ########################
        'dna_db'                           => {
            -dbname => $self->o('dna_db_name'),
            -host   => $self->o('dna_db_server'),
            -port   => $self->o('dna_db_port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('hive_driver'),
        },


        'ncrna_db'                               => {
            -dbname => $self->o('ncrna_db_name'),
            -host   => $self->o('ncrna_db_server'),
            -port   => $self->o('ncrna_db_port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('hive_driver'),
        },


    };
}

sub pipeline_create_commands {
    my ($self) = @_;

    return [
        @{$self->SUPER::pipeline_create_commands},
        'mkdir -p ' . $self->o('genome_dumps'),
        'mkdir -p ' . $self->o('ncrna_dir'),
    ];
}


sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{$self->SUPER::pipeline_wide_parameters},
        wide_ensembl_release => $self->o('ensembl_release'),
        genome_file          => $self->o('faidx_genome_file'),
    }
}

## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;
    return [

        {
            -logic_name => 'create_ncrna_db',
            -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
            -parameters => {
                source_db   => $self->o('dna_db'),
                target_db   => $self->o('ncrna_db'),
                create_type => 'clone',
            },
            -input_ids => [{}],
            -rc_name    => 'default',
            -flow_into => { '1' => 'patch_schema' },
        },

        {
            -logic_name => 'patch_schema',
            -module   => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => 'perl ' . $self->o('enscode_root_dir')  . '/ensembl/misc-scripts/schema_patcher.pl --database=' .
                $self->o('ncrna_db_name') . ' --host ' .
                $self->o('ncrna_db_server') . ' --port ' .
                $self->o('ncrna_db_port') . ' --user ' .
                $self->o('user') . ' --pass ' . 
                $self->o('password') . 
                ' --release #wide_ensembl_release# -fix --oldest 94 --nointeractive --type core -verbose',
            },
            -flow_into  => {
                '1->A' => [ 'dump_genome', 'dump_repeats', 'fetch_rfam_accessions' ],
                'A->1' => [ 'create_small_rna_slice_ids' ],
            },

        },

        {
            -logic_name => 'dump_repeats',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'perl ' . catfile($self->o('sncrna_analysis_script'), 'repeats_dump.pl') . ' ' .
                    $self->o('dna_db_name') . " " .
                    $self->o('dna_db_server') . " " .
                    $self->o('dna_db_port') . " " .
                    $self->o('user_r') . " " .
                    $self->o('ncrna_dir') . ' blastmirna',
            },
            -rc_name    => 'filter',
        },


        {
            -logic_name => 'dump_genome',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'perl ' . $self->o('sequence_dump_script') . ' -dbhost ' . $self->o('dna_db_server')
                    . ' -dbname  ' . $self->o('dna_db_name') . ' -dbport ' . $self->o('dna_db_port')
                    . ' -dbuser ' . $self->o('user_r')
                    . ' -coord_system_name toplevel -mask -mask_repeat ' . $self->o('full_repbase_logic_name')
                    . ' -output_dir ' . $self->o('genome_dumps')
                    . ' -softmask -onefile -header rnaseq -filename ' . $self->o('faidx_genome_file'),
            },
            -rc_name    => 'filter',
        },


        {
            -logic_name => 'create_small_rna_slice_ids',
            -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
            -parameters => {
                coord_system_name => 'toplevel',
                iid_type          => 'slice',
                slice_size        => 1000000,
                top_level         => 1,
                target_db         => $self->o('dna_db'),
                batch_slice_ids   => 1,
                batch_target_size => 2000000,
            },
            -flow_into  => {
                '2->A' => [ 'mirna_blast', 'run_cmsearch' ],
                'A->1' => [ 'filter_ncrnas' ],
            },
            -rc_name    => 'default',
        },


        {
            -logic_name    => 'mirna_blast',
            -module        => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastmiRNA',
            -parameters    => {
                repeat_logic_names => [ 'dust' ],
                repeat_db          => $self->o('dna_db'),
                output_db          => $self->o('ncrna_db'),
                dna_db             => $self->o('dna_db'),
                logic_name         => 'blastmirna',
                module             => 'HiveBlastmiRNA',
                blast_db_path      => catfile($self->o('mirna_blast_path'), $self->o('mirBase_fasta')),
                blast_exe_path     => $self->o('blastn_exe_path'),
                commandline_params => ' -num_threads 3 ',
                %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic', 'BlastmiRBase', { BLAST_PARAMS => { type => $self->o('blast_type') } })},
                timer              => '2h',
            },

            -flow_into     => {
                '-1' => [ 'rebatch_mirna' ],
                '-2' => [ 'rebatch_mirna' ],
            },
            -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
            -rc_name       => 'default',
        },


        {
            -logic_name   => 'rebatch_mirna',
            -module       => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
            -parameters   => {
                target_db         => $self->o('dna_db'),
                iid_type          => 'rebatch_and_resize_slices',
                slice_size        => 100000,
                batch_target_size => 100000,
            },
            -flow_into    => {
                '2' => [ 'mirna_blast_retry' ],
            },
            -rc_name      => 'default',
            -can_be_empty => 1,
        },


        {
            -logic_name    => 'mirna_blast_retry',
            -module        => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastmiRNA',
            -parameters    => {
                repeat_logic_names => [ 'dust' ],
                repeat_db          => $self->o('dna_db'),
                output_db          => $self->o('ncrna_db'),
                dna_db             => $self->o('dna_db'),
                logic_name         => 'blastmirna',
                module             => 'HiveBlastmiRNA',
                blast_db_path      => $self->o('mirna_blast_path') . '/' . $self->o('mirBase_fasta'),
                blast_exe_path     => $self->o('blastn_exe_path'),
                commandline_params => ' -num_threads 3 ',
                %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic', 'BlastmiRBase', { BLAST_PARAMS => { type => $self->o('blast_type') } })},
                timer              => '1h',
            },

            -flow_into     => {
                '-1' => [ 'failed_mirna_blast_job' ],
                '-2' => [ 'failed_mirna_blast_job' ],
            },
            -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
            -rc_name       => '5GB',
        },


        {
            -logic_name   => 'failed_mirna_blast_job',
            -module       => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -parameters   => {
            },
            -rc_name      => 'default',
            -can_be_empty => 1,
        },


        {
            -logic_name => 'filter_ncrnas',
            -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFilterncRNAs',
            -parameters => {
                output_db  => $self->o('ncrna_db'),
                dna_db     => $self->o('dna_db'),
                logic_name => 'filter_ncrnas',
                module     => 'HiveFilterncRNAs',
            },
            -rc_name    => 'filter',
            -flow_into  => {
                '2->A' => [ 'run_mirna' ],
                'A->1' => [ 'dump_features' ],

            },
        },


        {
            -logic_name    => 'run_mirna',
            -module        => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivemiRNA',
            -parameters    => {
                output_db         => $self->o('ncrna_db'),
                dna_db            => $self->o('dna_db'),
                logic_name        => 'ncrna',
                module            => 'HivemiRNA',
                blast_db_dir_path => catfile($self->o('mirna_blast_path'), 'all_mirnas.embl'),
                output_dir        => $self->o('ncrna_dir'),
            },
            -batch_size    => 20,
            -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
            -rc_name       => 'filter',

        },


        {
            -logic_name => 'dump_features',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'perl ' . catfile($self->o('sncrna_analysis_script'), 'dump_prefilter_features.pl') . ' ' .
                    $self->o('ncrna_db_name') . " " .
                    $self->o('ncrna_db_server') . " " .
                    $self->o('ncrna_db_port') . " " .
                    $self->o('user_r') . " " .
                    $self->o('ncrna_dir') . ' blastmirna',
            },
            -rc_name    => 'filter',
            -flow_into  => { '1' => 'filter_mirnas' },
        },


        {
            -logic_name => 'filter_mirnas',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'sh ' . catfile($self->o('sncrna_analysis_script'), 'FilterMiRNAs.sh')
                    . ' -d ' . catfile($self->o('ncrna_dir'), 'blastmirna_dafs.bed')
                    . ' -r ' . catfile($self->o('ncrna_dir'), 'repeats.bed')
                    . ' -g ' . $self->o('faidx_genome_file')
                    . ' -w ' . $self->o('ncrna_dir')
                    . ' -m ' . catfile($self->o('mirna_blast_path'), 'rfc_filters', $self->o('rfc_model'))
                    . ' -s ' . catfile($self->o('mirna_blast_path'), 'rfc_filters', $self->o('rfc_scaler'))
                    . ' -c ' . $self->o('ncrna_db_server') . ":" . $self->o('ncrna_db_port') . ":" . $self->o('ncrna_db_name') . ":"
                    . $self->o('user') . ":" . $self->o('password'),
            },
            -rc_name    => 'filter',
            -flow_into  => { '1->A' => 'dump_precursors', 'A->1' => 'ncrna_sanity_checks' },

        },

        {
            -logic_name => 'dump_precursors',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'perl ' . catfile($self->o('sncrna_analysis_script'), 'dump_precursors.pl')
                        . ' -dbname ' . $self->o('ncrna_db_name')
                        . ' -dbhost ' . $self->o('ncrna_db_server')
                        . ' -dbport ' . $self->o('ncrna_db_port')
                        . ' -dbuser ' . $self->o('user_r')
                        . ' -dnadbhost ' . $self->o('dna_db_server')
                        . ' -dnadbport ' . $self->o('dna_db_port')
                        . ' -dnadbname ' . $self->o('dna_db_name')
                        . ' -dnadbuser ' . $self->o('user_r')
                        . ' -output_dir ' . $self->o('ncrna_dir'),
            },
            -rc_name    => 'filter',
            -flow_into  => { '1' => 'map_mature_products' },

        },

        {
            -logic_name => 'map_mature_products',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'perl ' . catfile($self->o('sncrna_analysis_script'), 'map_mature_products.pl')
                        . ' -maturefa ' . catfile($self->o('mirna_blast_path'),$self->o('mature_mirnas'))
                        . ' -precursors ' . $self->o('ncrna_dir') . '/identified_mirnas.fa'
                        . ' -output_dir ' . $self->o('ncrna_dir'),

            },
            -rc_name    => '10GB',
            -flow_into  => { '1' => 'filter_mappings' },

        },

        {
            -logic_name => 'filter_mappings',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'python ' . catfile($self->o('sncrna_analysis_script'), 'filter_mappings.py')
                        . ' -s ' . catfile($self->o('ncrna_dir'), 'mature_mirnas.sam')
                        . ' -w ' . $self->o('ncrna_dir'),
            },
            -rc_name    => 'filter',
            -flow_into  => { '1' => 'populate_rnaproducts' },
        },
        {
            -logic_name => 'populate_rnaproducts',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'perl ' . catfile($self->o('sncrna_analysis_script'), 'populate_rnaproducts_table.pl')
                        . ' -dbname ' . $self->o('ncrna_db_name')
                        . ' -dbhost ' . $self->o('ncrna_db_server')
                        . ' -dbport ' . $self->o('ncrna_db_port')
                        . ' -dbuser ' . $self->o('user')
                        . ' -password ' . $self->o('password')
                        . ' -rnaproducts ' . catfile($self->o('ncrna_dir'), 'mirna_coords.txt')
            },
            -rc_name    => 'default',
        },

        {
            -logic_name => 'ncrna_sanity_checks',
            -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
            -parameters => {
                target_db                  => $self->o('ncrna_db'),
                sanity_check_type          => 'gene_db_checks',
                min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                    'gene_db_checks')->{$self->o('uniprot_set')}->{'ncrna'},
            },

            -rc_name    => '5GB',

        },
    ];
}


sub resource_classes {
    my $self = shift;

    return {
        'default'          => { LSF => $self->lsf_resource_builder('production-rh74', 1000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ]) },
        '5GB'              => { LSF => $self->lsf_resource_builder('production-rh74', 5000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ]) },
        '10GB'             => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ]) },
        '15GB'             => { LSF => $self->lsf_resource_builder('production-rh74', 15000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ]) },
        'filter'           => { LSF => $self->lsf_resource_builder('production-rh74', 3900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ]) },
        '10GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [ $self->default_options->{'pipe_db_server'} ], undef, ($self->default_options->{'use_threads'} + 1)) },

    }
}

sub hive_capacity_classes {
    my $self = shift;

    return {
        'hc_low'    => 200,
        'hc_medium' => 500,
        'hc_high'   => 1000,
    };
}

1;
