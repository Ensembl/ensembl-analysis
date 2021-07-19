=head1 LICENSE

Copyright [2021] EMBL-European Bioinformatics Institute

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

package RepeatsOnly__conf;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

### Not been fully tested but worked for Leps
### Creates core DBs with no gene set but containing repreat features.

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },
    'num_threads' => 20,

    ########### CHANGE THESE 3 #############
    'dbowner' => '', #FILL OUT
    'base_output_dir' => '', #FILL OUT
    'registry_file'             => '/hps/nobackup2/production/ensembl/fergal/production/test_runs/non_verts/production_test/Databases.pm', #CHANGE
    ########################################

    'release_number'            => '001' || $self->o('ensembl_release'),
    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'             => 'gbiab_prod_test' || $self->o('production_name').$self->o('production_name_modifier').'_'.$self->o('ensembl_release'),
    'user_r'                    => 'ensro', # read only db user
    'user'                      => 'ensadmin', # write db user
    'password'                  => 'ensembl', # password for write db user
    'server_set'                => '', # What server set to user, e.g. set1
    'pipe_db_server'            => $ENV{GBS7}, # host for pipe db
    'databases_server'          => $ENV{GBS5}, # host for general output dbs
    'dna_db_server'             => $ENV{GBS6}, # host for dna db
    'pipe_db_port'              => $ENV{GBP7}, # port for pipeline host
    'databases_port'            => $ENV{GBP5}, # port for general output db host
    'dna_db_port'               => $ENV{GBP6}, # port for dna db host
    'registry_db_server'        => $ENV{GBS1}, # host for registry db
    'registry_db_port'          => $ENV{GBP1}, # port for registry db
    'registry_db_name'          => 'gb_assembly_registry', # name for registry db
    'repbase_logic_name'        => '', # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    'repbase_library'           => '', # repbase library name, this is the actual repeat repbase library to use, e.g. "Mus musculus"
    'species_name'              => '', # e.g. mus_musculus
    'production_name'           => 'test_scale_gbiab' || $self->o('species_name'), # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'                  => '', # should be in the assembly report file
    'species_taxon_id'          => '' || $self->o('taxon_id'), # Species level id, could be different to taxon_id if we have a subspecies, used to get species level RNA-seq CSV data
    'genus_taxon_id'            => '' || $self->o('taxon_id'), # Genus level taxon id, used to get a genus level csv file in case there is not enough species level transcriptomic data

    'output_path'               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'             => '', # Name (as it appears in the assembly report file)
    'assembly_accession'        => '', # Versioned GCA assembly accession, e.g. GCA_001857705.1

    'use_genome_flatfile'       => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    'species_url'               => '' || $self->o('production_name').$self->o('production_name_modifier'), # sets species.url meta key
    'species_division'          => '', # sets species.division meta key


    'load_toplevel_only'           => '1', # This will not load the assembly info and will instead take any chromosomes, unplaced and unlocalised scaffolds directly in the DNA table
    'repeatmodeler_library'        => '', # This should be the path to a custom repeat library, leave blank if none exists
    'use_repeatmodeler_to_mask'    => '0', # Setting this will include the repeatmodeler library in the masking process

    'softmask_logic_names' => [],

########################
# Pipe and ref db info
########################
    'provider_name'                => 'Ensembl',
    'provider_url'                 => 'www.ensembl.org',

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_pipe_'.$self->o('release_number'),
    'dna_db_name'                   => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_core_'.$self->o('release_number'),

    'core_db_name'     => $self->o('dna_db_name'),
    'core_db_server'  => $self->o('dna_db_server'),
    'core_db_port'    => $self->o('dna_db_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    'production_db_server'         => 'mysql-ens-meta-prod-1',
    'production_db_port'           => '4483',

    databases_to_delete => ['reference_db', 'cdna_db', 'genblast_db', 'genewise_db', 'projection_db', 'selected_projection_db', 'layering_db', 'utr_db', 'genebuilder_db', 'pseudogene_db', 'ncrna_db', 'final_geneset_db', 'refseq_db', 'cdna2genome_db', 'rnaseq_blast_db', 'rnaseq_refine_db', 'rnaseq_rough_db', 'lincrna_db', 'otherfeatures_db', 'rnaseq_db'],#, 'projection_realign_db'

######################################################
#
# Mostly constant settings
#
######################################################

    genome_dumps                  => catdir($self->o('output_path'), 'genome_dumps'),
    # This one is used by most analyses that run against a genome flatfile like exonerate, genblast etc. Has slice name style headers. Is softmasked
    softmasked_genome_file        => catfile($self->o('genome_dumps'), $self->o('species_name').'_softmasked_toplevel.fa'),
    # This one is used in replacement of the dna table in the core db, so where analyses override slice->seq. Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file             => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),
    # This one is a cross between the two above, it has the seq_region name header but is softmasked. It is used by things that would both want to skip using the dna table and also want to avoid the repeat_feature table, e.g. bam2introns
    faidx_softmasked_genome_file  => catfile($self->o('genome_dumps'), $self->o('species_name').'_softmasked_toplevel.fa.reheader'),
    'primary_assembly_dir_name' => 'Primary_Assembly',
    'refseq_cdna_calculate_coverage_and_pid' => '0',
    'contigs_source'            => 'ena',
    full_repbase_logic_name => "repeatmask_repbase_".$self->o('repbase_logic_name'),
    'min_toplevel_slice_length'   => 250,
    'repeatmodeler_logic_name'    => 'repeatmask_repeatmodeler',
    'homology_models_path'        => catdir($self->o('output_path'),'homology_models'),

    ensembl_analysis_script           => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    remove_duplicates_script_path     => catfile($self->o('ensembl_analysis_script'), 'find_and_remove_duplicates.pl'),
    flag_potential_pseudogenes_script => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'flag_potential_pseudogenes.pl'),
    load_optimise_script              => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl'),
    prepare_cdnas_script              => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'prepare_cdnas.pl'),
    load_fasta_script_path            => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'load_fasta_to_db_table.pl'),
    loading_report_script             => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'report_genome_prep_stats.pl'),
    refseq_synonyms_script_path       => catfile($self->o('ensembl_analysis_script'), 'refseq', 'load_refseq_synonyms.pl'),
    refseq_import_script_path         => catfile($self->o('ensembl_analysis_script'), 'refseq', 'parse_ncbi_gff3.pl'),
    sequence_dump_script              => catfile($self->o('ensembl_analysis_script'), 'sequence_dump.pl'),
    sncrna_analysis_script             => catdir($self->o('ensembl_analysis_script'), 'genebuild', 'sncrna'),

    ensembl_misc_script        => catdir($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts'),
    repeat_types_script        => catfile($self->o('ensembl_misc_script'), 'repeats', 'repeat-types.pl'),
    meta_coord_script          => catfile($self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl'),
    meta_levels_script         => catfile($self->o('ensembl_misc_script'), 'meta_levels.pl'),
    assembly_name_script       => catfile($self->o('ensembl_analysis_script'), 'update_assembly_name.pl'),

########################
# Extra db settings
########################

    'num_tokens' => 10,
    mysql_dump_options => '--max_allowed_packet=1000MB',

########################
# Executable paths
########################
    'dust_path' => catfile($self->o('binary_base'), 'dustmasker'),
    'trf_path' => catfile($self->o('binary_base'), 'trf'),
    'eponine_java_path' => catfile($self->o('binary_base'), 'java'),
    'eponine_jar_path' => catfile($self->o('software_base_path'), 'opt', 'eponine', 'libexec', 'eponine-scan.jar'),
    'cpg_path' => catfile($self->o('binary_base'), 'cpg_lh'),
    'trnascan_path' => catfile($self->o('binary_base'), 'tRNAscan-SE'),
    'repeatmasker_path' => catfile($self->o('binary_base'), 'RepeatMasker'),

########################
# Misc setup info
########################
    'repeatmasker_engine'       => 'crossmatch',
    'masking_timer_long'        => '5h',
    'masking_timer_short'       => '2h',


########################################################
# URLs for retrieving the INSDC contigs and RefSeq files
########################################################
    'ncbi_base_ftp'           => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
    'insdc_base_ftp'          => $self->o('ncbi_base_ftp').'/#expr(substr(#assembly_accession#, 0, 3))expr#/#expr(substr(#assembly_accession#, 4, 3))expr#/#expr(substr(#assembly_accession#, 7, 3))expr#/#expr(substr(#assembly_accession#, 10, 3))expr#/#assembly_accession#_#assembly_name#',
    'assembly_ftp_path'       => $self->o('insdc_base_ftp'),
    'refseq_base_ftp'         => $self->o('ncbi_base_ftp').'/#expr(substr(#assembly_refseq_accession#, 0, 3))expr#/#expr(substr(#assembly_refseq_accession#, 4, 3))expr#/#expr(substr(#assembly_refseq_accession#, 7, 3))expr#/#expr(substr(#assembly_refseq_accession#, 10, 3))expr#/#assembly_refseq_accession#_#assembly_name#',
    'refseq_import_ftp_path'  => $self->o('refseq_base_ftp').'/#assembly_refseq_accession#_#assembly_name#_genomic.gff.gz',
    'refseq_mrna_ftp_path'    => $self->o('refseq_base_ftp').'/#assembly_refseq_accession#_#assembly_name#_rna.fna.gz',
    'refseq_report_ftp_path' => $self->o('refseq_base_ftp').'/#assembly_refseq_accession#_#assembly_name#_assembly_report.txt',

##################################
# Memory settings for the analyses
##################################
    'default_mem'          => '900',
    'refseq_mem'           => '9900',

########################
# db info
########################
    'core_db' => {
      -dbname => $self->o('dna_db_name'),
      -host   => $self->o('dna_db_server'),
      -port   => $self->o('dna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'production_db' => {
      -host   => $self->o('production_db_server'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ensembl_production',
      -driver => $self->o('hive_driver'),
    },

    'taxonomy_db' => {
      -host   => $self->o('production_db_server'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ncbi_taxonomy',
      -driver => $self->o('hive_driver'),
    },

    'registry_db' => {
      -host   => $self->o('registry_db_server'),
      -port   => $self->o('registry_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => $self->o('registry_db_name'),
      -driver => $self->o('hive_driver'),
    },
  };
}

sub pipeline_create_commands {
    my ($self) = @_;

    my $tables;
    my %small_columns = (
        paired => 1,
        read_length => 1,
        is_13plus => 1,
        is_mate_1 => 1,
        );

    # We need to store the values of the csv file to easily process it. It will be used at different stages
    foreach my $key (@{$self->default_options->{'file_columns'}}) {
        if (exists $small_columns{$key}) {
            $tables .= $key.' SMALLINT UNSIGNED NOT NULL,'
        }
        elsif ($key eq 'DS') {
            $tables .= $key.' VARCHAR(255) NOT NULL,'
        }
        else {
            $tables .= $key.' VARCHAR(50) NOT NULL,'
        }
    }
    $tables .= ' KEY(SM), KEY(ID)';

    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},

      $self->hive_data_table('protein', $self->o('uniprot_table_name')),

      $self->hive_data_table('refseq', $self->o('cdna_table_name')),

      $self->db_cmd('CREATE TABLE '.$self->o('realign_table_name').' ('.
                    'accession varchar(50) NOT NULL,'.
                    'seq text NOT NULL,'.
                    'PRIMARY KEY (accession))'),

      $self->db_cmd('CREATE TABLE '.$self->o('summary_csv_table')." ($tables)"),

      $self->db_cmd('CREATE TABLE '.$self->o('read_length_table').' ('.
                    'fastq varchar(50) NOT NULL,'.
                    'read_length int(50) NOT NULL,'.
                    'PRIMARY KEY (fastq))'),
    ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    wide_ensembl_release => $self->o('ensembl_release'),
  }
}

=head2 create_header_line

 Arg [1]    : Arrayref String, it will contains the values of 'file_columns'
 Example    : create_header_line($self->o('file_columns');
 Description: It will create a RG line using only the keys present in your csv file
 Returntype : String representing the RG line in a BAM file
 Exceptions : None

=cut

sub create_header_line {
    my ($items) = shift;

    my @read_tags = qw(ID SM DS CN DT FO KS LB PG PI PL PM PU);
    my $read_line = '@RG';
    foreach my $rt (@read_tags) {
        $read_line .= "\t$rt:#$rt#" if (grep($rt eq $_, @$items));
    }
    return $read_line."\n";
}

## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    my %commandline_params = (
      'ncbi' => '-num_threads 3 -window_size 40',
      'wu' => '-cpus 3 -hitdist 40',
      'legacy_ncbi' => '-a 3 -A 40',
      );
    my $header_line = create_header_line($self->default_options->{'file_columns'});

    return [

###############################################################################
#
# ASSEMBLY LOADING ANALYSES
#
###############################################################################
# 1) Process GCA - works out settings, flows them down the pipeline -> this should be seeded by another analysis later
# 2) Standard create core, populate tables, download data etc
# 3) Either run gbiab or setup gbiab
# 4) Finalise steps


      {
        # Creates a reference db for each species
        -logic_name => 'process_gca',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessGCA',
        -parameters => {
                         'num_threads'      => $self->o('num_threads'),
                         'dbowner'          => $self->o('dbowner'),
                         'core_db'          => $self->o('core_db'),
                         'ensembl_release'  => $self->o('ensembl_release'),
                         'base_output_dir'  => $self->o('base_output_dir'),
                         'registry_db'      => $self->o('registry_db'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'registry_file'      => $self->o('registry_file'),
                       },
        -rc_name    => 'default',

        -flow_into  => {
                          1 => ['create_core_db'],
                       },
        -analysis_capacity => 1,
        -input_ids  => [],
      },

      {
        # Creates a reference db for each species
        -logic_name => 'create_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'target_db'        => '#core_db#',
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'create_type'      => 'core_only',
                       },
        -rc_name    => 'default',

        -flow_into  => {
                         1 => ['populate_production_tables'],
                       },
      },

      {
        # Load production tables into each reference
        -logic_name => 'populate_production_tables',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
        -parameters => {
                         'target_db'        => '#core_db#',
                         'output_path'      => '#output_path#',
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'production_db'    => $self->o('production_db'),
                       },
        -rc_name    => 'default',

        -flow_into  => {
                         1 => ['process_assembly_info'],
                       },
      },

      {
        # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
        -logic_name => 'process_assembly_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessAssemblyReport',
        -parameters => {
          full_ftp_path => $self->o('assembly_ftp_path'),
          output_path   => '#output_path#',
          target_db     => '#core_db#',
        },
        -rc_name    => '8GB',
      	-max_retry_count => 0,
        -flow_into  => {
          1 => ['load_meta_info'],
        },
      },

      {
        # Load some meta info and seq_region_synonyms
        -logic_name => 'load_meta_info',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => '#core_db#',
          sql => [
            'DELETE FROM meta WHERE meta_key="species.display_name"',
            'INSERT IGNORE INTO meta (species_id, meta_key, meta_value) VALUES '.
              '(1, "annotation.provider_name", "Ensembl"),'.
              '(1, "annotation.provider_url", "www.ensembl.org"),'.
              '(1, "assembly.coverage_depth", "high"),'.
              '(1, "assembly.provider_name", NULL),'.
              '(1, "assembly.provider_url", NULL),'.
              '(1, "assembly.ucsc_alias", NULL),'.
              '(1, "species.stable_id_prefix", "#stable_id_prefix#"),'.
              '(1, "species.url", "#species_url#"),'.
              '(1, "species.display_name", "#species_display_name#"),'.
              '(1, "species.division", "#species_division#"),'.
              '(1, "species.strain", "#species_strain#"),'.
              '(1, "species.strain_group", "#species_strain_group#"),'.
              '(1, "species.production_name", "#production_name#"),'.
              '(1, "strain.type", "#strain_type#"),'.
              '(1, "repeat.analysis", "repeatdetector"),'.
              '(1, "repeat.analysis", "dust"),'.
              '(1, "repeat.analysis", "trf"),'.
              '(1, "genebuild.initial_release_date", NULL),'.
              '(1, "genebuild.projection_source_db", NULL),'.
              '(1, "genebuild.id", '.$self->o('genebuilder_id').'),'.
              '(1, "genebuild.method", "full_genebuild")'
          ],
        },
        -max_retry_count => 0,
        -rc_name    => 'default',
        -flow_into  => {
          1 => ['load_taxonomy_info'],
        },
      },

      {
        -logic_name => 'load_taxonomy_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo',
        -parameters => {
                         'target_db'        => '#core_db#',
                         'taxonomy_db'      => $self->o('taxonomy_db'),
                       },
        -rc_name    => 'default',

        -flow_into  => {
                          1 => ['dump_toplevel_file'],#['load_windowmasker_repeats'],# 'fan_refseq_import'],
                       },
      },

      {
        -logic_name => 'dump_toplevel_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
        -parameters => {
                         'coord_system_name'    => 'toplevel',
                         'target_db'            => '#core_db#',
                         'output_path'          => '#output_path#',
                         'enscode_root_dir'     => $self->o('enscode_root_dir'),
                         'species_name'         => '#species_name#',
                         'repeat_logic_names'   => $self->o('softmask_logic_names'), # This is emtpy as we just use masking present in downloaded file
                       },
        -flow_into => {
          1 => ['reheader_toplevel_file'],
        },
        -rc_name    => '3GB',
      },

      {
        -logic_name => 'reheader_toplevel_file',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'convert_genome_dump.pl').
                                      ' -conversion_type slice_name_to_seq_region_name'.
                                      ' -input_file #toplevel_genome_file#'.
                                      ' -output_file #reheadered_toplevel_genome_file#',
                       },
        -rc_name => 'default',
        -flow_into => {
          1 => ['run_gbiab'],
        },
      },

      {
        -logic_name => 'run_gbiab',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'python3.7 '.catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'gbiab', 'gbiab.py').' #gbiab_commandline#',
                       },
        -rc_name => 'gbiab',
        -max_retry_count => 0,
        -flow_into => {
          1 => ['set_meta_coords'],
        },
      },

      {
        -logic_name => 'set_meta_coords',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_coord_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('core_db','-host').
                                ' -port '.$self->o('core_db','-port').
                                ' -dbpattern '.'#core_dbname#'
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_meta_levels'],
                      },
      },

      {
        -logic_name => 'set_meta_levels',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_levels_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('core_db','-host').
                                ' -port '.$self->o('core_db','-port').
                                ' -dbname '.'#core_dbname#'
                       },
        -rc_name => 'default',
        -flow_into => {
          1 => ['null_columns'],
        },
      },

      {
        -logic_name => 'null_columns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => '#core_db#',
          sql => [
            'UPDATE gene SET stable_id = NULL',
            'UPDATE transcript SET stable_id = NULL',
            'UPDATE translation SET stable_id = NULL',
            'UPDATE exon SET stable_id = NULL',
            'UPDATE protein_align_feature set external_db_id = NULL',
            'UPDATE dna_align_feature set external_db_id = NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['final_meta_updates'],
                      },
      },

      {
        -logic_name => 'final_meta_updates',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => '#core_db#',
          sql => [
            'INSERT INTO meta (species_id, meta_key, meta_value) VALUES '.
              '(1, "genebuild.last_geneset_update", (SELECT CONCAT((EXTRACT(YEAR FROM now())),"-",(LPAD(EXTRACT(MONTH FROM now()),2,"0")))))'
          ],
        },
        -rc_name    => 'default',
        -flow_into  => {
                          1 => ['populate_analysis_descriptions'],
                       },
      },

      {
        -logic_name => 'populate_analysis_descriptions',
        -module     => 'Bio::EnsEMBL::Production::Pipeline::ProductionDBSync::PopulateAnalysisDescription',
        -parameters => {
                        species => '#production_name#',
                        group => 'core',
                       },
        -rc_name    => 'default_registry',
      },
    ];
}


sub resource_classes {
  my $self = shift;

  return {
    'default_registry' => { LSF => [$self->lsf_resource_builder('production-rh74', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{''}]},
    'gbiab' => { LSF => $self->lsf_resource_builder('production-rh74', 20000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, $self->default_options->{'num_threads'})},
    '3GB' => { LSF => $self->lsf_resource_builder('production-rh74', 3000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '8GB' => { LSF => $self->lsf_resource_builder('production-rh74', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
      'default' => { LSF => $self->lsf_resource_builder('production-rh74', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
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

sub check_file_in_ensembl {
  my ($self, $file_path) = @_;
  push @{$self->{'_ensembl_file_paths'}}, $file_path;
  return $self->o('enscode_root_dir').'/'.$file_path;
}
1;
