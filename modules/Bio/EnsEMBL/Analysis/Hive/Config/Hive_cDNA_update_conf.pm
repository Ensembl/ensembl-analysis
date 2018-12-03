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


package Hive_cDNA_update_conf;

use strict;
use warnings;

use File::Spec::Functions qw(catfile catdir);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;           # Allow this particular config to use conditional dataflow and INPUT_PLUS

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

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
    'release_number'  => '', # What release are you doing this for?
    'species'         => '', # either mus_musculus or homo_sapiens
    'strategy'        => 'update', # set this to update or complete. Update will just align the new cDNAs, complete will realign all of them
    'password'        => '',

    'base_output_dir' => '',

    'refseq_version'  => 90, # look at the refseq home page to get the latest version, this will go into the cDNA database
    'ena_version'     => 136,  # look at the EBI ENA home page - or news

    # Database connection info:
    'pipe_db_server'  => 'mysql-ens-genebuild-prod-7',
    'pipe_db_port'    => 4533,

    'output_db_server' => 'mysql-ens-genebuild-prod-3',
    'output_db_port'   => 4529,

    # details of the last cdna db (eg. on livemirror)
    'core_db_name'     => $self->o('species').'_core_'.$self->o('release_number').'_'.$self->o('coord_system_version'),
    'old_cdna_db_name' => $self->o('species').'_cdna_'.$self->o('release_number').'_'.$self->o('coord_system_version'), # This works because production copies the old DB and patch the schema on staging

    'dna_db_server' => 'mysql-ens-genebuild-prod-2',
    'dna_db_port'   => 4528,
    'dna_db_user'   => 'ensadmin',
    'dna_db_pass'   => $self->o('password'), 


##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

    'pipeline_name' => $self->o('species').'_cdna_update_'.$self->o('release_number'),
    'coord_system_version' => 38,

    'dna_db_name'    => $self->o('dbowner').'_'.$self->o('core_db_name'),
    'output_db_name'   => $self->o('dbowner').'_'.$self->o('species').'_cdna_'.$self->o('release_number').'_'.$self->o('coord_system_version'),

    'hive_capacity' => 100,
    'default_queue'              => 'production-rh7',

    'output_path' => catdir($self->o('base_output_dir'), $self->o('species'), $self->o('release_number')), # output directory you want to place downloaded files and log files
    'genome_file' => catfile($self->o('output_path'), 'genome', $self->o('species').'_softmasked_toplevel.fa'), #The softmasked genome file for either homo_sapiens or mus_musculus
    'repeat_masking_logic_names' => ['dust', 'repeatmask_repbase_'.$self->o('species')], # the repeatmask logic name(s) of the analyses used in the genebuild

    'exonerate_batch_size'       => 100,
    'exonerate_time_limit'       => '2h',

    'target_logic_name'          => 'cdna_update',
    'target_biotype'             => 'cdna_update',
    'cdna_file_name'             => 'cdna_update',

    'cdna_batch_size' => 1,

    'user'                       => 'ensadmin',
    'user_r'                     => 'ensro',
    'password_r'                 => undef,

    'refseq_ftp'                 => 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian',
    'ena_ftp'                    => 'ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence',

    'cdna_table_name'            => 'cdna_sequences',

    'default_mem'                => '900',
    'exonerate_mem'              => '3900',
    'exonerate_retry_mem'        => '5900',
    'exonerate_high_mem'         => '9900',
    'optimise_mem'               => '11900',

    'exonerate_pid'              => '97',
    'exonerate_cov'              => '90',
    'exonerate_version' => '0.9.0', # This is unlikely to change any time soon

    'exonerate_path'             => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),

    'optimize_script'            => catfile($self->o('enscode_root_dir').'', 'ensembl-analysis', 'scripts', 'genebuild', 'load_external_db_ids_and_optimize_af.pl'),

    'analysis_desc_script'       => catfile($self->o('enscode_root_dir').'', 'ensembl-production', 'scripts', 'production_database', 'populate_analysis_description.pl'),
    'populate_production_script' => catfile($self->o('enscode_root_dir').'', 'ensembl-production', 'scripts', 'production_database', 'populate_production_db_tables.pl'),

    'meta_level_script'          => catfile($self->o('enscode_root_dir').'', 'ensembl', 'misc-scripts', 'meta_levels.pl'),
    'meta_coord_script'          => catfile($self->o('enscode_root_dir').'', 'ensembl', 'misc-scripts', 'meta_coord', 'update_meta_coord.pl'),

    'killlist_db_server' => 'mysql-ens-genebuild-prod-6.ebi.ac.uk',
    'killlist_db_port' => 4532,

    'production_db_name' => 'ensembl_production_'.$self->o('release_number'),
    'production_db_server' => 'mysql-ens-sta-1',
    'production_db_port' => 4519,

    'old_cdna_db_server' => $self->o('production_db_server'),
    'old_cdna_db_port' => $self->o('production_db_port'),

    'core_db_server' => $self->o('production_db_server'),
    'core_db_port' => $self->o('production_db_port'),

    'databases_to_delete' => ['target_db'], # Do not put dna_db here as it can be really dnagerous
    'production_db' => {
      -dbname => $self->o('production_db_name'),
      -host => $self->o('production_db_server'),
      -port => $self->o('production_db_port'),
      -user => $self->o('user_r'),
      -pass => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'target_db' => {
      -dbname => $self->o('output_db_name'),
      -host => $self->o('output_db_server'),
      -port => $self->o('output_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host => $self->o('killlist_db_server'),
      -port => $self->o('killlist_db_port'),
      -user => $self->o('user_r'),
      -pass => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'old_cdna_db' => {
      -dbname => $self->o('old_cdna_db_name'),
      -host => $self->o('old_cdna_db_server'),
      -port => $self->o('old_cdna_db_port'),
      -user => $self->o('user_r'),
      -pass => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'core_db' => {
      -dbname => $self->o('core_db_name'),
      -host => $self->o('core_db_server'),
      -port => $self->o('core_db_port'),
      -user => $self->o('user_r'),
      -pass => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },
  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [
  # inheriting database and hive tables' creation
    @{$self->SUPER::pipeline_create_commands},

    $self->db_cmd('CREATE TABLE '.$self->o('cdna_table_name').' ('.
      'accession varchar(50) NOT NULL,'.
      'seq text NOT NULL,'.
      'biotype varchar(50) NOT NULL,'.
      'PRIMARY KEY (accession))'),
  ];
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  my %taxon_id = (
    homo_sapiens => 9606,
    mus_musculus => 10090,
  );
  my %ena_species = (
    homo_sapiens => 'hum',
    mus_musculus => 'mus',
  );
  my %refseq_binomial = (
    homo_sapiens => 'Homo sapiens',
    mus_musculus => 'Mus musculus',
  );

  my @analyses = (
    {
      -logic_name => 'create_output_dir',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mkdir -p '.$self->o('output_path')
      },
      -max_retry_count => 0,
      -input_ids => [{
        cdna_file => catfile('#wide_output_dir#', $self->o('cdna_file_name')),
      }],
      -flow_into => {
        1 => [ 'create_dna_db'],
      }
    },
    {
      -logic_name => 'create_dna_db',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('core_db'),
        target_db => $self->o('dna_db'),
        create_type => 'dna_db',
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['create_output_db'],
      }
    },
    {
      -logic_name => 'create_output_db',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('target_db'),
        create_type => 'clone',
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['output_db_cleanup'],
      }
    },
    {
      -logic_name => 'output_db_cleanup',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('target_db'),
        sql => [
          'DELETE FROM analysis WHERE logic_name != "'.$self->o('target_logic_name').'"',
          'ALTER TABLE analysis AUTO_INCREMENT = 1',
          'DELETE FROM meta WHERE species_id IS NOT NULL AND meta_key NOT IN ('.
            '"assembly.accession", '.
            '"assembly.coverage_depth", '.
            '"assembly.date", '.
            '"assembly.default", '.
            '"assembly.long_name", '.
            '"assembly.mapping", '.
            '"assembly.name", '.
            '"assembly.ucsc_alias", '.
            '"genebuild.id", '.
            '"liftover.mapping", '.
            '"lrg", '.
            '"species.alias", '.
            '"species.classification", '.
            '"species.common_name", '.
            '"species.display_name", '.
            '"species.division", '.
            '"species.production_name", '.
            '"species.scientific_name", '.
            '"species.stable_id_prefix", '.
            '"species.taxonomy_id", '.
            '"species.url"'.
          ')',
        ],
      },
      -flow_into => {
        1 => [ 'populate_production_tables' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'populate_production_tables',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #wide_production_script#'.
               ' -dp #wide_output_dir# --dropbaks'.
               ' -d '.$self->o('target_db','-dbname').
               ' -h '.$self->o('target_db','-host').
               ' -u '.$self->o('target_db','-user').
               ' -p '.$self->o('target_db','-pass').
               ' -P '.$self->o('target_db','-port').
               ' -md '.$self->o('production_db','-dbname').
               ' -mh '.$self->o('production_db','-host').
               ' -mu '.$self->o('production_db','-user').
               ' -mP '.$self->o('production_db','-port').
               ' -t external_db -t attrib_type -t misc_set -t unmapped_reason'
      },
      -max_retry_count => 0,
      -flow_into => {
        '1->A' => ['dump_softmasked_toplevel', 'create_ena_jobs', 'create_refseq_jobs'],
        'A->1' => ['generate_jobs'],
      },
    },
    {
      -logic_name => 'dump_softmasked_toplevel',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
      -parameters => {
                       'coord_system_name'    => 'toplevel',
                       'target_db'            => $self->o('dna_db'),
                       'output_path'          => catfile($self->o('output_path'), 'genome'),
                       'enscode_root_dir'     => $self->o('enscode_root_dir'),
                       'species_name'         => $self->o('species'),
                       'repeat_logic_names'   => $self->o('repeat_masking_logic_names'),
                     },
      -rc_name    => 'exonerate',
    },
    {
      -logic_name => 'create_ena_jobs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFTPURLs',
      -parameters => {
        file_list => ['release/std/rel_htc_'.$ena_species{$self->o('species')}.'_*', 'release/std/rel_std_'.$ena_species{$self->o('species')}.'_*', 'update/std/cum_htc_'.$ena_species{$self->o('species')}.'_*', 'update/std/cum_std_'.$ena_species{$self->o('species')}.'_*'],
        base_url => $self->o('ena_ftp'),
      },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into =>  {
        2 => ['download_ena_files'],
      }
    },
    {
      -logic_name => 'download_ena_files',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData',
      -parameters => {
        output_dir => '#wide_output_dir#',
        download_method => 'ftp',
      },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into =>  {
        2 => ['prepare_ena_file'],
      }
    },
    {
      -logic_name => 'prepare_ena_file',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePreparecDNAs',
      -parameters => {
        format => 'embl',
        data_type => 'ena',
      },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into =>  {
        2 => ['load_cdnas'],
      }
    },
    {
      -logic_name => 'create_refseq_jobs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFTPURLs',
      -parameters => {
        file_list => ['vertebrate_mammalian*rna.fna.gz'],
        base_url => $self->o('refseq_ftp'),
      },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into =>  {
        2 => ['download_refseq_files'],
      }
    },
    {
      -logic_name => 'download_refseq_files',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData',
      -parameters => {
        output_dir => '#wide_output_dir#',
        download_method => 'ftp',
      },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into =>  {
        2 => ['prepare_refseq_file'],
      }
    },
    {
      -logic_name => 'prepare_refseq_file',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePreparecDNAs',
      -parameters => {
        description_filter => $refseq_binomial{$self->o('species')},
        data_type => 'refseq',
      },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into =>  {
        2 => ['load_cdnas'],
      }
    },
    {
      -logic_name => 'load_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadUpdatecDNAs',
      -parameters => {
        dna_db => $self->o('dna_db'),
        old_cdna_db => $self->o('old_cdna_db'),
        new_cdna_db => $self->o('target_db'),
        strategy => $self->o('strategy'),
        killlist_db => $self->o('killlist_db'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
    },
    {
      -logic_name => 'generate_jobs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        batch_size => $self->o('cdna_batch_size'),
        sequence_table_name => $self->o('cdna_table_name'),
        target_db => $self->o('dna_db'),
        iid_type => 'sequence_accession',
      },
      -rc_name => 'default',
      -flow_into => {
        '2->A' => ['exonerate'],
        'A->1' => [ 'find_many_hits' ],
      },
    },
    {
      -logic_name => 'exonerate',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('target_db'),
        logic_name => $self->o('target_logic_name'),
        timer => $self->o('exonerate_time_limit'),
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cdnaupdate')},
        exonerate_pid => $self->o('exonerate_pid'),
        exonerate_cov => $self->o('exonerate_cov'),
        GENOMICSEQS => '#wide_genome_file#',
        PROGRAM => '#wide_exonerate#',
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -flow_into => {
        -1 => ['exonerate_himem'],
        -2 => ['exonerate_second_run'],
        2 => ['exonerate_second_run'],
      },
      -rc_name => 'exonerate',
      -batch_size => $self->o('exonerate_batch_size'),
    },
    {
      -logic_name => 'exonerate_himem',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('target_db'),
        logic_name => $self->o('target_logic_name'),
        timer => $self->o('exonerate_time_limit'),
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cdnaupdate')},
        exonerate_pid => $self->o('exonerate_pid'),
        exonerate_cov => $self->o('exonerate_cov'),
        GENOMICSEQS => '#wide_genome_file#',
        PROGRAM => '#wide_exonerate#',
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -flow_into => {
        2 => ['exonerate_second_run'],
        -1 => ['exonerate_second_run'],
        -2 => ['exonerate_second_run'],
      },
      -rc_name => 'exonerate_himem',
      -batch_size => $self->o('exonerate_batch_size'),
      -can_be_empty => 1,
    },
    {
      -logic_name => 'exonerate_second_run',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('target_db'),
        logic_name => $self->o('target_logic_name'),
        timer => $self->o('exonerate_time_limit'),
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cdnaupdate_loose')},
        exonerate_pid => $self->o('exonerate_pid'),
        exonerate_cov => $self->o('exonerate_cov'),
        GENOMICSEQS => '#wide_genome_file#',
        PROGRAM => '#wide_exonerate#',
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -rc_name => 'exonerate_2',
      -can_be_empty => 1,
      -batch_size => $self->o('exonerate_batch_size'),
      -flow_into => {
        -1 => ['exonerate_second_run_himem'],
        -2 => ['store_unmapped'],
        2 => ['store_unmapped'],
      },
    },
    {
      -logic_name => 'exonerate_second_run_himem',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('target_db'),
        logic_name => $self->o('target_logic_name'),
        timer => $self->o('exonerate_time_limit'),
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cdnaupdate_loose')},
        exonerate_pid => $self->o('exonerate_pid'),
        exonerate_cov => $self->o('exonerate_cov'),
        GENOMICSEQS => '#wide_genome_file#',
        PROGRAM => '#wide_exonerate#',
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -rc_name => 'exonerate_himem',
      -batch_size => $self->o('exonerate_batch_size'),
      -can_be_empty => 1,
      -flow_into => {
        -1 => ['store_unmapped'],
        -2 => ['store_unmapped'],
        2 => ['store_unmapped'],
      },
    },
    {
      -logic_name => 'store_unmapped',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStoreUnmappedcDNAs',
      -parameters => {
        target_db => $self->o('target_db'),
        logic_name => $self->o('target_logic_name'),
      },
      -max_retry_count => 0,
      -rc_name => 'default',
    },
    {
      -logic_name => 'find_many_hits',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivecDNAManyHits',
      -parameters => {
        target_db => $self->o('target_db'),
        old_db => $self->o('old_cdna_db'),
      },
      -rc_name => 'exonerate',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['database_compare'],
        2 => ['check_many_hits'],
      },
    },
    {
      -logic_name => 'check_many_hits',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name => 'default',
      -max_retry_count => 0,
    },
    {
      -logic_name => 'database_compare',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveComparecDNAdbs',
      -parameters => {
        old_cdna_db => $self->o('old_cdna_db'),
        new_cdna_db => $self->o('target_db'),
        output_file => '#wide_output_dir#/comparison.out',
      },
      -flow_into => {
        1 => [ 'comparison_report'],
      },
    },
    {
      -logic_name => 'comparison_report',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
      -parameters => {
        email => $self->o('email_address'),
        subject => 'AUTOMATED REPORT: cDNA update database comparison',
        text => 'Please find below the counts for each toplevel seq_region for the current and the previous cDNA updates:',
        file => '#wide_output_dir#/comparison.out',,
      },
      -flow_into => {
        1 => [ 'null_align_features' ],
      },
    },
    {
      -logic_name => 'null_align_features',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('target_db'),
        sql => [
          "UPDATE dna_align_feature SET external_db_id=NULL",
          "UPDATE protein_align_feature SET external_db_id=NULL",
        ],
      },
      -flow_into => {
        1 => [ 'load_xdbids' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'load_xdbids',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #wide_optimize_script#'.
               ' -output_path #wide_output_dir#/optimise_daf_paf.1'.
               ' -dbname '.$self->o('target_db','-dbname').
               ' -dbhost '.$self->o('target_db','-host').
               ' -dbport '.$self->o('target_db','-port').
               ' -dbuser '.$self->o('target_db','-user').
               ' -dbpass '.$self->o('target_db','-pass').
               ' -prod_dbname '.$self->o('production_db','-dbname').
               ' -prod_dbhost '.$self->o('production_db','-host').
               ' -prod_dbport '.$self->o('production_db','-port').
               ' -prod_dbuser '.$self->o('production_db','-user').
               ' -verbose -clean -nopaf -no_backup'
      },
      -rc_name => 'optimise',
      -max_retry_count => 0,
      -flow_into => {
        1 => [ 'populate_analysis_description' ],
      },
    },
    {
      -logic_name => 'populate_analysis_description',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #wide_analysis_desc#'.
               ' -dp #wide_output_dir# --dropbak'.
               ' -d '.$self->o('target_db','-dbname').
               ' -h '.$self->o('target_db','-host').
               ' -u '.$self->o('target_db','-user').
               ' -p '.$self->o('target_db','-pass').
               ' -P '.$self->o('target_db','-port').
               ' -md '.$self->o('production_db','-dbname').
               ' -mh '.$self->o('production_db','-host').
               ' -mu '.$self->o('production_db','-user').
               ($self->o('production_db','-pass') ? ' -mp '.$self->o('production_db','-pass') : '').
               ' -mP '.$self->o('production_db','-port').
               ' -s #wide_species# -t cdna'
      },
      -flow_into => {
        1 => [ 'handover_preparation' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'handover_preparation',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('target_db'),
        sql => [
          'UPDATE gene SET biotype = "'.$self->o('target_biotype').'"',
          'UPDATE transcript SET biotype = "'.$self->o('target_biotype').'"',
          'UPDATE gene g, transcript t SET g.canonical_transcript_id = t.transcript_id WHERE g.gene_id = t.gene_id',
          'UPDATE analysis SET db = "RefSeq,ENA", db_version = "'.$self->o('refseq_version').','.$self->o('ena_version').'", program_version = "'.$self->o('exonerate_version').'" WHERE logic_name = "'.$self->o('target_logic_name').'"',
        ],
      },
      -flow_into => {
        1 => [ 'meta_levels' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'meta_levels',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('meta_level_script').
               ' --dbpattern '.$self->o('target_db','-dbname').
               ' --host '.$self->o('target_db','-host').
               ' --user '.$self->o('target_db','-user').
               ' --pass '.$self->o('target_db','-pass').
               ' --port '.$self->o('target_db','-port')
      },
      -flow_into => {
        1 => [ 'update_meta_coord' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'update_meta_coord',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('meta_coord_script').
               ' --dbpattern '.$self->o('target_db','-dbname').
               ' --host '.$self->o('target_db','-host').
               ' --user '.$self->o('target_db','-user').
               ' --pass '.$self->o('target_db','-pass').
               ' --port '.$self->o('target_db','-port')
      },
      -max_retry_count => 0,
    },
  );

  foreach my $analysis (@analyses) {
    $analysis->{'-hive_capacity'} = $self->o('hive_capacity') unless (exists $analysis->{'-hive_capacity'});
  }

  return \@analyses;
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    strategy => $self->o('strategy'),
    wide_output_dir => $self->o('output_path'),
    wide_genome_file => $self->o('genome_file'),
    wide_optimize_script => $self->o('optimize_script'),
    wide_exonerate => $self->o('exonerate_path'),
    wide_production_script => $self->o('populate_production_script'),
    wide_analysis_desc => $self->o('analysis_desc_script'),
    wide_species => $self->o('species'),
  };
}

sub resource_classes {
  my $self = shift;

  my $default_queue = $self->default_options()->{'default_queue'};

  my $default_mem = $self->default_options()->{'default_mem'};
  my $exonerate_mem = $self->default_options()->{'exonerate_mem'};
  my $exonerate_retry_mem = $self->default_options()->{'exonerate_retry_mem'};
  my $exonerate_high_mem = $self->default_options()->{'exonerate_high_mem'};
  my $optimise_mem = $self->default_options()->{'optimise_mem'};

  return {
    'default' => { LSF => '-q '.$default_queue.' -M'.$default_mem.' -R"select[mem>'.$default_mem.'] '.
                          'rusage[mem='.$default_mem.']"'},

    'exonerate' => { LSF => '-q '.$default_queue.' -M'.$exonerate_mem.' -R"select[mem>'.$exonerate_mem.'] '.
                            'rusage[mem='.$exonerate_mem.']"'},

    'exonerate_himem' => { LSF => '-q '.$default_queue.' -M'.$exonerate_high_mem.' -R"select[mem>'.$exonerate_high_mem.'] '.
                           'rusage[mem='.$exonerate_high_mem.']"'},

    'exonerate_2' => { LSF => '-q '.$default_queue.' -M'.$exonerate_retry_mem.' -R"select[mem>'.$exonerate_retry_mem.'] '.
                              'rusage[mem='.$exonerate_retry_mem.']"'},

    'optimise' => { LSF => '-q '.$default_queue.' -M'.$optimise_mem.' -R"select[mem>'.$optimise_mem.'] '.
                           'rusage[mem='.$optimise_mem.']"'},

  }
}


1;
