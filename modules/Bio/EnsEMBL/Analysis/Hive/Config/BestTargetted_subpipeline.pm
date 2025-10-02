
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::BestTargetted_subpipeline;


use strict;
use warnings;
use File::Spec::Functions;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');
use Bio::EnsEMBL::Analysis::Tools::SoftwareConfigLoad qw(get_software_path); #Software path config module



sub default_options {
  my ($self) = @_;
  ## Build software path based on new software type
  my $software_type = $ENV{SOFTWARE_TYPE};
  my $exonerate_path = get_software_path($software_type, 'exonerate');
  my $genewise_path = get_software_path($software_type, 'genewise');
  my $indicate_path = get_software_path($software_type, 'indicate');
  my $pmatch_path = get_software_path($software_type, 'pmatch');
  my $exonerate_annotation = get_software_path($software_type, 'exonerate');

  return {

    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

    ######################################################
    #
    # Variable settings- You change these!!!
    #
    ######################################################
    ########################
    # Misc setup info
    ########################
    software_type             => $software_type,
    'dbowner' => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name' => '' || $self->o('production_name') . '_' . $self->o('ensembl_release'),
    'user_r'           => '',    # read only db user
    'user'             => '',    # write db user
    'password'         => '',    # password for write db user
    'pipe_db_host'     => '',    # host for pipe db
    'dna_db_host'      => '',    # host for dna db
    'databases_host'   => '',    # host for general output dbs

    'pipe_db_port'   => '',      # port for pipeline host
    'dna_db_port'    => '',      # port for dna db host
    'databases_port' => '',      # port for general output db host

    'release_number' => '' || $self->o('ensembl_release'),
    'species_name'    => '',     # e.g. mus_musculus
    'production_name' => '',     # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    dbname_accession          => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1

    'taxon_id'    => '',         # should be in the assembly report file
    use_genome_flatfile => 1,
    repeat_logic_names => [],

    'output_path'   => '',                                               # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    targetted_path  => catdir( $self->o('output_path'), 'targetted' ),
    cdna_file       => catfile( $self->o('targetted_path'), 'cdnas' ),
    annotation_file => $self->o('cdna_file') . '.annotation',

    'uniprot_table_name' => 'uniprot_sequences',

    exonerate_logic_name => 'exonerate',
    ncbi_query           => '((txid' . $self->o('taxon_id') . '[Organism:noexp]+AND+biomol_mrna[PROP]))  NOT "tsa"[Properties] NOT EST[keyword]',

    cdna_table_name                             => 'cdna_sequences',
    target_exonerate_calculate_coverage_and_pid => 0,
    cdna_selection_pid                          => '97', # Cut-off for percent id for selecting the cDNAs
    cdna_selection_cov                          => '90', # Cut-off for coverage for selecting the cDNAs
    cdna2genome_region_padding                  => 2000,
    exonerate_max_intron                        => 200000,
    best_targetted_min_coverage                 => 50,                 # This is to avoid having models based on fragment alignment and low identity
    best_targetted_min_identity                 => 50,                 # This is to avoid having models based on fragment alignment and low identity

    ########################
    # Pipe and ref db info
    ########################

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_core_' . $self->o('release_number'),

    'cdna_db_host'   => $self->o('databases_host'),
    'cdna_db_port'   => $self->o('databases_port'),

    'cdna2genome_db_host'   => $self->o('databases_host'),
    'cdna2genome_db_port'   => $self->o('databases_port'),

    'genewise_db_host'   => $self->o('databases_host'),
    'genewise_db_port'   => $self->o('databases_port'),

    'best_targeted_db_host'   => $self->o('databases_host'),
    'best_targeted_db_port'   => $self->o('databases_port'),

    'killlist_db_host'   => $self->o('databases_host'),
    'killlist_db_port'   => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release' => $ENV{ENSEMBL_RELEASE},    # this is the current release version on staging to be able to get the correct database

    
    ########################
    # Executable paths
    ########################
    genewise_path                    => $genewise_path,
    exonerate_path                   => $exonerate_path,
    indicate_path                    => $indicate_path,
    pmatch_path                      => $pmatch_path,
    exonerate_annotation             => $exonerate_annotation,
    
    ######################################################
    #
    # Mostly constant settings
    #
    ######################################################

    genome_dumps => catdir( $self->o('output_path'), 'genome_dumps' ),
    # This one is used by most analyses that run against a genome flatfile like exonerate, genblast etc. Has slice name style headers. Is softmasked
    softmasked_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_softmasked_toplevel.fa' ),
    faidx_genome_file => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),

    # This is used for "messaging" other sub pipeline
    transcript_selection_url => undef,

    ensembl_analysis_script => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    prepare_cdnas_script    => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'prepare_cdnas.pl'),


    ########################
    # db info
    ########################

    'cdna_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_cdna_' . $self->o('release_number'),
      -host   => $self->o('cdna_db_host'),
      -port   => $self->o('cdna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'cdna2genome_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_cdna2genome_' . $self->o('release_number'),
      -host   => $self->o('cdna2genome_db_host'),
      -port   => $self->o('cdna2genome_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genewise_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_genewise_' . $self->o('release_number'),
      -host   => $self->o('genewise_db_host'),
      -port   => $self->o('genewise_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'best_targeted_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_bt_' . $self->o('release_number'),
      -host   => $self->o('best_targeted_db_host'),
      -port   => $self->o('best_targeted_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host   => $self->o('killlist_db_host'),
      -port   => $self->o('killlist_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    databases_to_delete => ['cdna_db', 'cdna2genome_db', 'genewise_db', 'best_targeted_db'],

  };
}

sub pipeline_create_commands {
  my ($self) = @_;

  return [
    @{$self->SUPER::pipeline_create_commands},
    $self->hive_data_table('protein', $self->o('uniprot_table_name')),
    $self->hive_data_table('refseq', $self->o('cdna_table_name')),

  ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    use_genome_flatfile => $self->o('use_genome_flatfile'),
    genome_file => $self->o('faidx_genome_file'),
    wide_repeat_logic_names => $self->o('repeat_logic_names'),
  }
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [

    ######################################################################################
    #
    # cDNA alignment
    #
    ######################################################################################

    {
      -logic_name => 'create_cdna_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('cdna_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -input_ids  => [{}],
      -flow_into => {
        '1->A' => [ 'create_genewise_db', 'download_mRNA' ],
        'A->1' => ['create_besttargetted_db'],
      },
    },

    {
      -logic_name => 'create_genewise_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('genewise_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        1      => ['download_selenocysteines'],
        '1->A' => [ 'download_uniprot_self', 'download_refseq_self' ],
        'A->1' => ['load_self'],
      },
    },

    {
      -logic_name => 'download_refseq_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadNCBImRNA',
      -parameters => {
        query              => '#species_name#[Organism] AND RefSeq[Filter]',
        species_name       => $self->o('species_name'),
        output_file        => catfile( $self->o('targetted_path'), 'ncbi_self.fa' ),
        ncbidb             => 'protein',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['?accu_name=iid&accu_address=[]&accu_input_variable=output_file'],
      },
    },

    {
      -logic_name => 'download_uniprot_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
      -parameters => {
        taxon_id             => $self->o('taxon_id'),
        multi_query_download => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', 'self_isoforms_12' ),
        output_path          => $self->o('targetted_path'),
        _branch_to_flow_to   => 1,
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['?accu_name=iid&accu_address=[]&accu_input_variable=iid'],
      },
    },

    {
      -logic_name => 'load_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
      -parameters => {
        output_file => catfile( $self->o('targetted_path'), 'proteome.fa' ),
        skip_Xs     => 5,
        killlist_db => $self->o('killlist_db'),
        killlist_type       => 'protein',
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['indicate_proteome'],
        2 => ['targetted_exonerate'],
      },
    },

    {
      -logic_name => 'targetted_exonerate',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -parameters => {
        iid_type            => 'db_seq',
        sequence_table_name => $self->o('uniprot_table_name'),
        dna_db              => $self->o('dna_db'),
        target_db           => $self->o('genewise_db'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic', 'exonerate_protein' ) },
        genome_file                => $self->o('softmasked_genome_file'),
        exonerate_path             => $self->o('exonerate_path'),
        repeat_libraries           => '#wide_repeat_logic_names#',
        calculate_coverage_and_pid => $self->o('target_exonerate_calculate_coverage_and_pid'),
        use_genome_flatfile => 0,
      },
      -rc_name => '3GB',
      -hive_capacity => $self->o('hc_normal'),
      -flow_into => {
        -1 => ['targetted_exonerate_retry'],
      },
    },

    {
      -logic_name => 'targetted_exonerate_retry',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -parameters => {
        iid_type => 'db_seq',
        sequence_table_name => $self->o('uniprot_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('genewise_db'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_protein')},
        genome_file      => $self->o('softmasked_genome_file'),
        exonerate_path   => $self->o('exonerate_path'),
        repeat_libraries => '#wide_repeat_logic_names#',
        calculate_coverage_and_pid => $self->o('target_exonerate_calculate_coverage_and_pid'),
        use_genome_flatfile => 0,
      },
      -rc_name => '10GB',
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'indicate_proteome',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd            => '#indicate_path# -d #indicate_dir# -f #proteome# -i #proteome_index# -p singleWordParser',
        indicate_path  => $self->o('indicate_path'),
        proteome       => 'proteome.fa',
        indicate_dir   => $self->o('targetted_path'),
        proteome_index => catdir( $self->o('targetted_path'), 'proteome_index' ),
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['generate_pmatch_jobs'],
      },
    },

    {
      -logic_name => 'generate_pmatch_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db         => $self->o('dna_db'),
        coord_system_name => 'toplevel',
        iid_type          => 'slice',
        top_level         => 1,
      },
      -rc_name   => '2GB',
      -flow_into => {
        '2->A' => ['pmatch'],
        'A->1' => ['bestpmatch'],
      },
    },

    {
      -logic_name => 'pmatch',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePmatch',
      -parameters => {
        target_db         => $self->o('genewise_db'),
        dna_db            => $self->o('dna_db'),
        PROTEIN_FILE      => catfile( $self->o('targetted_path'), 'proteome.fa' ),
        MIN_COVERAGE      => 25,
        BINARY_LOCATION   => $self->o('pmatch_path'),
        REPEAT_MASKING    => [],
        MAX_INTRON_LENGTH => 50000,
        OPTIONS           => '-T 20',                                                # set threshold to 14 for more sensitive search
      },
      -rc_name => '3GB',
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'bestpmatch',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestPmatch',
      -parameters => {
        source_db         => $self->o('genewise_db'),
        target_db         => $self->o('genewise_db'),
        PMATCH_LOGIC_NAME => ['pmatch'],
        MIN_COVERAGE      => 50,
      },
      -rc_name   => '2GB',
      -flow_into => {
        1 => ['generate_targetted_jobs'],
      },
    },

    {
      -logic_name => 'generate_targetted_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type          => 'feature_region',
        feature_type      => 'protein_align_feature',
        target_db         => $self->o('genewise_db'),
        logic_name        => ['bestpmatch'],
        coord_system_name => 'toplevel',
      },
      -rc_name   => 'default',
      -flow_into => {
        2 => [ 'targetted_genewise_gtag', 'targetted_exo' ],
      },
    },

    {
      -logic_name => 'targetted_genewise_gtag',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise',
      -parameters => {
        source_db         => $self->o('genewise_db'),
        target_db         => $self->o('genewise_db'),
        killlist_db       => $self->o('killlist_db'),
        dna_db            => $self->o('dna_db'),
        gtag              => 0,                                                             # 0 is for gtag, 1 is for non canonical
        biotype           => 'gw_gtag',
        max_intron_length => 200000,
        disconnect_jobs   => 1,
        seqfetcher_index  => [ catfile( $self->o('targetted_path'), 'proteome_index' ) ],
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_genewise' ) },
      },
      -rc_name => '3GB',
      -hive_capacity => $self->o('hc_normal'),
      -flow_into => {
        MEMLIMIT => ['targetted_genewise_gtag_6GB'],
      },
    },

    {
      -logic_name => 'targetted_genewise_gtag_6GB',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise',
      -parameters => {
        source_db         => $self->o('genewise_db'),
        target_db         => $self->o('genewise_db'),
        killlist_db       => $self->o('killlist_db'),
        dna_db            => $self->o('dna_db'),
        gtag              => 0,                                                             # 0 is for gtag, 1 is for non canonical
        biotype           => 'gw_gtag',
        max_intron_length => 200000,
        seqfetcher_index  => [ catfile( $self->o('targetted_path'), 'proteome_index' ) ],
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_genewise' ) },
      },
      -rc_name => '6GB',
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'targetted_exo',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerateForGenewise',
      -parameters => {
        source_db         => $self->o('genewise_db'),
        target_db         => $self->o('genewise_db'),
        killlist_db       => $self->o('killlist_db'),
        dna_db            => $self->o('dna_db'),
        biotype           => 'gw_exo',
        seqfetcher_index  => [ catfile( $self->o('targetted_path'), 'proteome_index' ) ],
        max_intron_length => 700000,
        program_file      => $self->o('exonerate_path'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_exonerate' ) },
      },
      -rc_name => '3GB',
      -hive_capacity => $self->o('hc_normal'),
      -flow_into => {
        MEMLIMIT => ['targetted_exo_6GB'],
      },
    },

    {
      -logic_name => 'targetted_exo_6GB',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerateForGenewise',
      -parameters => {
        source_db         => $self->o('genewise_db'),
        target_db         => $self->o('genewise_db'),
        killlist_db       => $self->o('killlist_db'),
        dna_db            => $self->o('dna_db'),
        biotype           => 'gw_exo',
        seqfetcher_index  => [ catfile( $self->o('targetted_path'), 'proteome_index' ) ],
        max_intron_length => 700000,
        program_file      => $self->o('exonerate_path'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_exonerate' ) },
      },
      -rc_name => '6GB',
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'download_mRNA',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadNCBImRNA',
      -parameters => {
        output_file => $self->o('cdna_file'),
        filetype    => 'gb',
        query       => $self->o('ncbi_query'),
      },
      -rc_name   => 'default',
      -flow_into => {
        2 => ['prepare_cdna'],
      },
    },

    {
      -logic_name => 'prepare_cdna',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('prepare_cdnas_script') .
          ' -killdbnames ' . $self->o( 'killlist_db', '-dbname' ) .
          ' -killdbhost ' . $self->o( 'killlist_db', '-host' ) .
          ' -killdbuser ' . $self->o( 'killlist_db', '-user' ) .
          ' -killdbport ' . $self->o( 'killlist_db', '-port' ) .
          ( $self->o( 'killlist_db', '-pass' ) ? ' -killdbpass ' . $self->o( 'killlist_db', '-pass' ) : '' ) .
          ' -infile #sequence_file#' .
          ' -outfile #sequence_file#.clipped' .
          ( $self->o('taxon_id') ? ' -tax_id ' . $self->o('taxon_id') : '' ) .
          ' -nomole',
        sequence_file => $self->o('cdna_file'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['load_cdna_file'],
      },
    },

    {
      -logic_name => 'load_cdna_file',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadmRNAs',
      -parameters => {
        sequence_table_name => $self->o('cdna_table_name'),
        filetype            => 'fasta',
        sequence_file       => $self->o('cdna_file') . '.clipped',
      },
      -flow_into => {
        '2->A' => ['exonerate'],
        'A->1' => ['prepare_cdna2genome'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'exonerate',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '3GB',
      -parameters => {
        iid_type            => 'db_seq',
        sequence_table_name => $self->o('cdna_table_name'),
        dna_db              => $self->o('dna_db'),
        target_db           => $self->o('cdna_db'),
        logic_name          => $self->o('exonerate_logic_name'),
        genome_file         => $self->o('softmasked_genome_file'),
        exonerate_path      => $self->o('exonerate_path'),
        repeat_libraries    => '#wide_repeat_logic_names#',
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic', 'cdna_est2genome' ) },
        exonerate_cdna_pid         => 50,
        exonerate_cdna_cov         => 50,
        calculate_coverage_and_pid => 0,
      },
      -batch_size => 100,
      -hive_capacity => $self->o('hc_normal'),
      -flow_into  => {
        -1 => ['exonerate_retry'],
      },
      -batch_size           => 100,
      -failed_job_tolerance => 5,
    },

    {
      -logic_name => 'exonerate_retry',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '6GB',
      -parameters => {
        iid_type            => 'db_seq',
        sequence_table_name => $self->o('cdna_table_name'),
        dna_db              => $self->o('dna_db'),
        target_db           => $self->o('cdna_db'),
        logic_name          => $self->o('exonerate_logic_name'),
        genome_file         => $self->o('softmasked_genome_file'),
        exonerate_path      => $self->o('exonerate_path'),
        repeat_libraries    => '#wide_repeat_logic_names#',
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic', 'cdna_est2genome' ) },
        exonerate_cdna_pid         => 50,
        exonerate_cdna_cov         => 50,
        calculate_coverage_and_pid => 0,
      },
      -batch_size           => 100,
      -hive_capacity => $self->o('hc_normal'),
      -failed_job_tolerance => 100,
      -can_be_empty         => 1,
    },

    {
      -logic_name => 'prepare_cdna2genome',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('prepare_cdnas_script') .
          ' -killdbnames ' . $self->o( 'killlist_db', '-dbname' ) .
          ' -killdbhost ' . $self->o( 'killlist_db', '-host' ) .
          ' -killdbuser ' . $self->o( 'killlist_db', '-user' ) .
          ' -killdbport ' . $self->o( 'killlist_db', '-port' ) .
          ( $self->o( 'killlist_db', '-pass' ) ? ' -killdbpass ' . $self->o( 'killlist_db', '-pass' ) : '' ) .
          ' -infile #sequence_file#' .
          ' -outfile #sequence_file#.cds' .
          ' -annotation ' . $self->o('annotation_file') .
          ( $self->o('taxon_id') ? ' -tax_id ' . $self->o('taxon_id') : '' ) .
          ' -nomole',
        sequence_file => $self->o('cdna_file'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_cdna2genome_db'],
      },
    },

    {
      -logic_name => 'create_cdna2genome_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('cdna2genome_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['create_cdna_toplevel_slices'],
      },
    },

    {
      -logic_name => 'create_cdna_toplevel_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('cdna_db'),
        iid_type              => 'slice',
        coord_system_name     => 'toplevel',
        include_non_reference => 0,
        top_level             => 1,
        feature_constraint    => 1,
        feature_type          => 'gene',
      },
      -flow_into => {
        '2' => ['apply_threshold'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'apply_threshold',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectGeneOnFilter',
      -parameters => {
        dna_db     => $self->o('dna_db'),
        source_db  => $self->o('cdna_db'),
        logic_name => 'cdna_alignment',
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic', 'cdna_est2genome' ) },
        exonerate_cdna_pid => $self->o('cdna_selection_pid'),
        exonerate_cdna_cov => $self->o('cdna_selection_cov'),
      },
      -rc_name           => 'default',
      -analysis_capacity => 5,
      -batch_size        => 10,
      -flow_into         => {
        '1' => ['create_cdna2genome_slices'],
      },
    },

    {
      -logic_name => 'create_cdna2genome_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('cdna_db'),
        iid_type              => 'feature_region',
        feature_type          => 'gene',
        logic_name            => ['cdna_alignment'],
        coord_system_name     => 'toplevel',
        include_non_reference => 0,
        top_level             => 1,
        use_annotation        => 1,
        # These options will create only slices that have a gene on the slice in one of the feature dbs
        annotation_file => $self->o('annotation_file'),
        region_padding  => $self->o('cdna2genome_region_padding'),
      },
      -flow_into => {
        '2->A' => ['cdna2genome'],
        'A->1' => ['internal_stop'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'cdna2genome',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion',
      -rc_name    => '3GB',
      -parameters => {
        iid_type            => 'db_seq',
        dna_db              => $self->o('dna_db'),
        sequence_table_name => $self->o('cdna_table_name'),
        source_db           => $self->o('cdna_db'),
        target_db           => $self->o('cdna2genome_db'),
        logic_name          => 'cdna2genome',
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic', 'cdna2genome' ) },
        calculate_coverage_and_pid => 1,
        exonerate_path             => $self->o('exonerate_annotation'),
        annotation_file            => $self->o('annotation_file'),
        SOFT_MASKED_REPEATS        => '#wide_repeat_logic_names#',
      },
      -batch_size => 10,
      -hive_capacity => $self->o('hc_normal'),
      -flow_into  => {
        '-1' => ['cdna2genome_himem'],
      },
    },

    {
      -logic_name => 'cdna2genome_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion',
      -rc_name    => '6GB',
      -parameters => {
        iid_type            => 'db_seq',
        dna_db              => $self->o('dna_db'),
        sequence_table_name => $self->o('cdna_table_name'),
        source_db           => $self->o('cdna_db'),
        target_db           => $self->o('cdna2genome_db'),
        logic_name          => 'cdna2genome',
        module              => 'HiveExonerate2Genes',
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic', 'cdna2genome' ) },
        calculate_coverage_and_pid => 1,
        exonerate_path             => $self->o('exonerate_annotation'),
        annotation_file            => $self->o('annotation_file'),
        repeat_libraries           => '#wide_repeat_logic_names#',
      },
      -batch_size => 10,
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'internal_stop',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveInternalStopFix',
      -parameters => {
        dna_db             => $self->o('dna_db'),
        source_db          => $self->o('cdna2genome_db'),
        edited_biotype     => 'edited',
        stop_codon_biotype => 'stop_codon',
        logic_name         => 'cdna2genome',
        biotype            => undef,
        source             => undef,
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'download_selenocysteines',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
      -parameters => {
        taxon_id => $self->o('taxon_id'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', 'selenocysteine' ) },
        output_path => $self->o('targetted_path'),
      },
      -rc_name   => 'default',
      -flow_into => {
        2 => ['load_selenocysteine'],
      },
    },

    {
      -logic_name => 'load_selenocysteine',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
      -parameters => {
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name   => 'default',
      -flow_into => {
        2 => ['process_selenocysteine'],
      },
    },

    {
      -logic_name => 'process_selenocysteine',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelenocysteineFinder',
      -parameters => {
        target_db           => $self->o('genewise_db'),
        dna_db              => $self->o('dna_db'),
        genome              => $self->o('softmasked_genome_file'),
        biotype             => 'seleno_self',
        exonerate           => $self->o('exonerate_path'),
        genewise            => $self->o('genewise_path'),
        iid_type            => 'db_seq',
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'create_besttargetted_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('best_targeted_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['generate_besttargetted_index'],
      },
    },

    {
      -logic_name => 'generate_besttargetted_index',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateBestTargettedIndex',
      -parameters => {
        source_db        => $self->o('cdna2genome_db'),
        seqfetcher_index => [ catfile( $self->o('targetted_path'), 'proteome_index' ) ],
        fasta_filename   => catfile( $self->o('targetted_path'), 'best_targetted.fa' ),
        email            => $self->o('email_address'),
        genbank_file     => $self->o('cdna_file'),
        protein_files    => [ catfile( $self->o('targetted_path'), 'proteome.fa' ), catfile( $self->o('targetted_path'), $self->o('taxon_id') . '_seleno.fa' ) ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['indicate_BT'],
      },
    },

    {
      -logic_name => 'indicate_BT',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd           => '#indicate_path# -d #indicate_dir# -f #proteome# -i #indicate_dir#/best_targetted_index -p singleWordParser -M BTMultiParser',
        indicate_path => $self->o('indicate_path'),
        proteome      => 'best_targetted.fa',
        indicate_dir  => $self->o('targetted_path'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['generate_besttargetted_jobs'],
      },
    },

    {
      -logic_name => 'generate_besttargetted_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db          => $self->o('genewise_db'),
        coord_system_name  => 'toplevel',
        iid_type           => 'slice',
        feature_constraint => 1,
        feature_type       => 'gene',
        top_level          => 1,
        feature_dbs        => [ $self->o('genewise_db'), $self->o('cdna2genome_db') ],
      },
      -rc_name   => 'default',
      -flow_into => {
        '2->A' => ['best_targetted'],
        'A->1' => ['best_targetted_healthchecks'],
      },
    },

    {
      -logic_name => 'best_targetted',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestTargetted',
      -parameters => {
        target_db => $self->o('best_targeted_db'),
        dna_db    => $self->o('dna_db'),
        source_db => { protein_db => $self->o('genewise_db'), cdna2genome_db => $self->o('cdna2genome_db') },
        SEQFETCHER_DIR => [ catfile( $self->o('targetted_path'), 'proteome_index' ),
          catfile( $self->o('targetted_path'), 'best_targetted_index' ) ],
        SEQFETCHER_OBJECT   => 'Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher',
        INPUT_DATA_FROM_DBS => {
          protein_db => [ 'seleno_self', 'gw_gtag', 'gw_nogtag', 'gw_exo', 'targetted_exonerate' ],
          cdna2genome_db => [ 'cdna2genome', 'edited' ],
        },
        BIOTYPES => [ 'seleno_self', 'cdna2genome', 'edited', 'gw_gtag', 'gw_nogtag', 'gw_exo', 'targetted_exonerate' ],    # sorted list, preferred is first
        protein_min_coverage => $self->o('best_targetted_min_coverage'),
        protein_min_identity => $self->o('best_targetted_min_identity'),
      },
      -rc_name => '3GB',
    },

    {
      -logic_name => 'best_targetted_healthchecks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
      -parameters => {
        input_db => $self->o('best_targeted_db'),
        species  => $self->o('species_name'),
        group    => 'protein_cdna',
      },
      -max_retry_count => 0,
      -rc_name         => 'default',
      -flow_into  => {
        '1'  => ['notification_pipeline_is_done'],
      },
    },

    {
      -logic_name => 'notification_pipeline_is_done',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::MessagePipeline',
      -parameters => {
        messages   => [
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'create_toplevel_slices',
          param => 'feature_dbs',
          data => $self->o('best_targeted_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'split_slices_on_intergenic',
          param => 'input_gene_dbs',
          data => $self->o('best_targeted_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'layer_annotation',
          param => 'SOURCEDB_REFS',
          data => $self->o('best_targeted_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition',
          param => 'donor_dbs',
          data => $self->o('cdna_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition_10GB',
          param => 'donor_dbs',
          data => $self->o('cdna_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition_30GB',
          param => 'donor_dbs',
          data => $self->o('cdna_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'utr_memory_failover',
          param => 'donor_dbs',
          data => $self->o('cdna_db'),
          update => 1,
        }],
        tweak_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'tweak_pipeline.pl'),
      },
      -rc_name    => 'default',
    },

  ];
}


sub resource_classes {
  my $self = shift;

  return {
    #inherit other stuff from the base class
     %{ $self->SUPER::resource_classes() },
    }
}

1;
