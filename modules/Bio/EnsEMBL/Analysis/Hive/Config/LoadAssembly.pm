=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::LoadAssembly;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf; # Needed for the WHEN ELSE structure
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
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
    dbowner                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    pipeline_name             => '' || $self->o('production_name').'_'.$self->o('ensembl_release'),
    user_r                    => '', # read only db user
    user                      => '', # write db user
    password                  => '', # password for write db user
    pipe_db_server            => '', # host for pipe db
    dna_db_server             => '', # host for dna db
    pipe_db_port              => '', # port for pipeline host
    dna_db_port               => '', # port for dna db host
    repbase_logic_name        => '', # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    release_number            => '' || $self->o('ensembl_release'),
    species_name              => '', # e.g. mus_musculus
    production_name           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    taxon_id                  => '', # should be in the assembly report file
    output_path               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    wgs_id                    => '', # Can be found in assembly report file on ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
    assembly_name             => '', # Name (as it appears in the assembly report file)
    assembly_accession        => '', # Versioned GCA assembly accession, e.g. GCA_001857705.1
    stable_id_prefix          => '', # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    species_url               => '', # sets species.url meta key
    species_division          => 'EnsemblVertebrates', # sets species.division meta key

    # Keys for custom loading, only set/modify if that's what you're doing
    load_toplevel_only        => '1', # This will not load the assembly info and will instead take any chromosomes, unplaced and unlocalised scaffolds directly in the DNA table
    custom_toplevel_file_path => '', # Only set this if you are loading a custom toplevel, requires load_toplevel_only to also be set to 2
    use_repeatmodeler_to_mask => '0', # Setting this will include the repeatmodeler library in the masking process

    red_logic_name                   => 'repeatdetector', # logic name for the Red repeat finding analysis
    replace_repbase_with_red_to_mask => '0', # Setting this will replace 'full_repbase_logic_name' with 'red_logic_name' repeat features in the masking process

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is
    assembly_provider_name   => '',
    assembly_provider_url    => '',
    annotation_provider_name => 'Ensembl',
    annotation_provider_url  => 'www.ensembl.org',


########################
# Pipe and ref db info
########################
    pipe_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_load_assembly_pipe_'.$self->o('release_number'),
    dna_db_name  => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),

    reference_db_name   => $self->o('dna_db_name'),
    reference_db_server => $self->o('dna_db_server'),
    reference_db_port   => $self->o('dna_db_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    ensembl_release      => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    production_db_server => 'mysql-ens-meta-prod-1',
    production_db_port   => '4483',
    production_db_name   => 'ensembl_production',

    taxonomy_db_server => $self->o('production_db_server'),
    taxonomy_db_port   => $self->o('production_db_port'),
    taxonomy_db_name   => 'ncbi_taxonomy',

    projection_source_db_name => '', # This is generally a pre-existing db, like the current human/mouse core for example

    databases_to_delete => ['reference_db'],#, 'projection_realign_db'


######################################################
#
# Mostly constant settings
#
######################################################
    genome_dumps => catdir($self->o('output_path'), 'genome_dumps'),
    # This one is used in replacement of the dna table in the core db, so where analyses override slice->seq. Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),

    primary_assembly_dir_name => 'Primary_Assembly',
    contigs_source            => 'ena',

    full_repbase_logic_name  => "repeatmask_repbase_".$self->o('repbase_logic_name'),
    repeatmodeler_logic_name => 'repeatmask_repeatmodeler',

    ensembl_analysis_script => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    sequence_dump_script    => catfile($self->o('ensembl_analysis_script'), 'sequence_dump.pl'),


########################
# Extra db settings
########################
    num_tokens => 10,


########################
# Executable paths
########################
    samtools_path => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
########################################################
# URLs for retrieving the INSDC contigs and RefSeq files
########################################################
    ncbi_base_ftp           => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
    insdc_base_ftp          => $self->o('ncbi_base_ftp').'/#expr(substr(#assembly_accession#, 0, 3))expr#/#expr(substr(#assembly_accession#, 4, 3))expr#/#expr(substr(#assembly_accession#, 7, 3))expr#/#expr(substr(#assembly_accession#, 10, 3))expr#/#assembly_accession#_#assembly_name#',
    assembly_ftp_path       => $self->o('insdc_base_ftp'),


########################
# db info
########################
    reference_db => {
      -dbname => $self->o('reference_db_name'),
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    production_db => {
      -host   => $self->o('production_db_server'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => $self->o('production_db_name'),
      -driver => $self->o('hive_driver'),
    },

    taxonomy_db => {
      -host   => $self->o('taxonomy_db_server'),
      -port   => $self->o('taxonomy_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => $self->o('taxonomy_db_name'),
      -driver => $self->o('hive_driver'),
    },
  };
}

sub pipeline_create_commands {
  my ($self) = @_;

  return [
  # inheriting database and hive tables' creation
    @{$self->SUPER::pipeline_create_commands},
    'mkdir -p '.$self->o('genome_dumps'),
  ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    load_toplevel_only => $self->o('load_toplevel_only'),
  }
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'create_core_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        'target_db'        => $self->o('reference_db'),
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'create_type'      => 'core_only',
      },
      -rc_name    => 'default',
      -input_ids  => [
        {
          assembly_name => $self->o('assembly_name'),
          assembly_accession => $self->o('assembly_accession'),
        },
      ],
      -flow_into  => {
        1 => ['populate_production_tables'],
      },
    },

    {
      -logic_name => 'populate_production_tables',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
      -parameters => {
        'target_db'        => $self->o('reference_db'),
        'output_path'      => $self->o('output_path'),
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'production_db'    => $self->o('production_db'),
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => WHEN ('#load_toplevel_only# == 1' => ['process_assembly_info'],
              '#load_toplevel_only# == 2' => ['custom_load_toplevel'],
            ELSE ['download_assembly_info']),
      },
    },


    ####
    # Loading custom assembly where the user provide a FASTA file, probably a repeat library
    ####
    {
      -logic_name => 'custom_load_toplevel',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'assembly_loading', 'load_seq_region.pl').
          ' -dbhost '.$self->o('reference_db_server').
          ' -dbuser '.$self->o('user').
          ' -dbpass '.$self->o('password').
          ' -dbport '.$self->o('reference_db_port').
          ' -dbname '.$self->o('reference_db_name').
          ' -coord_system_version '.$self->o('assembly_name').
          ' -default_version'.
          ' -coord_system_name primary_assembly'.
          ' -rank 1'.
          ' -fasta_file '. $self->o('custom_toplevel_file_path').
          ' -sequence_level'.
          ' -noverbose',
      },
      -rc_name => '4GB',
      -flow_into => {
        1 => ['custom_set_toplevel'],
      },
    },

    {
      -logic_name => 'custom_set_toplevel',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'assembly_loading', 'set_toplevel.pl').
          ' -dbhost '.$self->o('reference_db_server').
          ' -dbuser '.$self->o('user').
          ' -dbpass '.$self->o('password').
          ' -dbport '.$self->o('reference_db_port').
          ' -dbname '.$self->o('reference_db_name'),
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => ['custom_add_meta_keys'],
      },
    },

    {
      -logic_name => 'custom_add_meta_keys',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql => [
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"assembly.default","'.$self->o('assembly_name').'")',
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"assembly.name","'.$self->o('assembly_name').'")',
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"species.taxonomy_id","'.$self->o('taxon_id').'")',
        ],
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['load_meta_info'],
      },
    },


    ####
    # Loading assembly with only the toplevel sequences loaded, it fetches the data from INSDC databases
    ####
    # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
    {
      -logic_name => 'process_assembly_info',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessAssemblyReport',
      -parameters => {
        full_ftp_path => $self->o('assembly_ftp_path'),
        output_path   => $self->o('output_path'),
        target_db     => $self->o('reference_db'),
      },
      -rc_name    => '8GB',
      -max_retry_count => 0,
      -flow_into  => {
        1 => ['load_meta_info'],
      },
    },

    # Load some meta info and seq_region_synonyms
    {
      -logic_name => 'load_meta_info',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql => [
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES '.
          '(1, "species.stable_id_prefix", "'.$self->o('stable_id_prefix').'"),'.
          '(1, "species.url", "'.ucfirst($self->o('species_url')).'"),'.
          '(1, "species.division", "'.$self->o('species_division').'"),'.
          '(1, "genebuild.initial_release_date", NULL),'.
          '(1, "assembly.coverage_depth", "high"),'.
          '(1, "genebuild.id", '.$self->o('genebuilder_id').'),'.
          '(1, "genebuild.method", "full_genebuild"),'.
          '(1, "genebuild.projection_source_db", "'.$self->o('projection_source_db_name').'"),'.
          '(1, "assembly.provider_name", "'.$self->o('assembly_provider_name').'"),'.
          '(1, "assembly.provider_url", "'.$self->o('assembly_provider_url').'"),'.
          '(1, "annotation.provider_name", "'.$self->o('annotation_provider_name').'"),'.
          '(1, "annotation.provider_url", "'.$self->o('annotation_provider_url').'"),'.
          '(1, "species.production_name", "'.$self->o('production_name').'"),'.
          ($self->o('replace_repbase_with_red_to_mask') ? '(1, "repeat.analysis", "'.$self->o('red_logic_name').'"),' :
            '(1, "repeat.analysis", "'.$self->o('full_repbase_logic_name').'"),').
          ($self->o('use_repeatmodeler_to_mask') ? '(1, "repeat.analysis", "'.$self->o('repeatmodeler_logic_name').'"),': '').
          '(1, "repeat.analysis", "dust"),'.
          '(1, "repeat.analysis", "trf"),'.
          '(1, "species.strain_group", "'.$self->o('species_name').'")',
        ],
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => ['load_taxonomy_info'],
      },
    },


    ####
    # Loading assembly with the full assembly representation, it fetches data from INSDC databases
    ####
    # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
    {
      -logic_name => 'download_assembly_info',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadNCBIFtpFiles',
      -parameters => {
        'full_ftp_path'             => $self->o('assembly_ftp_path'),
        'output_path'               => $self->o('output_path'),
        'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => ['find_contig_accessions'],
      },
    },

    # Get the prefixes for all contigs from the AGP files
    {
      -logic_name => 'find_contig_accessions',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveFindContigAccessions',
      -parameters => {
        'contigs_source'            => $self->o('contigs_source'),
        'wgs_id'                    => $self->o('wgs_id'),
        'output_path'               => $self->o('output_path'),
        'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['download_contigs'],
      },
    },

    # Download contig from NCBI
    {
      -logic_name => 'download_contigs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadContigs',
      -parameters => {
        'contigs_source' => $self->o('contigs_source'),
        'wgs_id' => $self->o('wgs_id'),
        'output_path' => $self->o('output_path'),
        'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['load_contigs'],
      },
    },


    # Load the contigs into each reference db
    {
      -logic_name => 'load_contigs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadSeqRegions',
      -parameters => {
        'coord_system_version' => $self->o('assembly_name'),
        'target_db' => $self->o('reference_db'),
        'output_path' => $self->o('output_path'),
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
      },
      -rc_name => '4GB',
      -flow_into => {
        1 => ['load_assembly_info'],
      },
    },

    # Load the AGP files
    {
      -logic_name => 'load_assembly_info',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadAssembly',
      -parameters => {
        'target_db' => $self->o('reference_db'),
        'output_path' => $self->o('output_path'),
        'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['set_toplevel'],
      },
    },

    # Set the toplevel
    {
      -logic_name => 'set_toplevel',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetAndCheckToplevel',
      -parameters => {
        'target_db' => $self->o('reference_db'),
        'output_path' => $self->o('output_path'),
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
      },
      -rc_name => '2GB',
      -flow_into => {
        1 => ['load_meta_info_full'],
      },
    },

    # Load some meta info and seq_region_synonyms
    {
      -logic_name => 'load_meta_info_full',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetMetaAndSeqRegionSynonym',
      -parameters => {
        'target_db' => $self->o('reference_db'),
        'output_path' => $self->o('output_path'),
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
        'meta_key_list' => {
          'assembly.accession' => $self->o('assembly_accession'),
          'assembly.coverage_depth' => 'high',
          'assembly.default' => $self->o('assembly_name'),
          'assembly.name' => $self->o('assembly_name'),
          'assembly.web_accession_source' => 'NCBI',
          'assembly.web_accession_type' => 'GenBank Assembly ID',
          'genebuild.id' => $self->o('genebuilder_id'),
          'genebuild.method' => 'full_genebuild',
          'genebuild.projection_source_db' => $self->o('projection_source_db_name'),
          'assembly.provider_name' => $self->o('assembly_provider_name'),
          'assembly.provider_url' => $self->o('assembly_provider_url'),
          'annotation.provider_name' => $self->o('annotation_provider_name'),
          'annotation.provider_url' => $self->o('annotation_provider_url'),
          'repeat.analysis' => [$self->o('full_repbase_logic_name'), 'dust', 'trf'],
          'species.production_name' => $self->o('production_name'),
          'species.taxonomy_id' => $self->o('taxon_id'),
        }
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['load_taxonomy_info'],
      },
    },


    {
      -logic_name => 'load_taxonomy_info',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo',
      -parameters => {
        'target_db'        => $self->o('reference_db'),
        'taxonomy_db'      => $self->o('taxonomy_db'),
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => ['load_mitochondrion'],
      },
    },

    {
      -logic_name => 'load_mitochondrion',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion',
      -parameters => {
        'target_db'        => $self->o('reference_db'),
        'output_path'      => $self->o('output_path'),
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'species_name'     => $self->o('species_name'),
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => ['create_faidx_genome_file'],
      },
    },

    {
      -logic_name => 'create_faidx_genome_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => '1GB',
      -parameters => {
        cmd => 'if [ ! -s "'.$self->o('faidx_genome_file').'" ]; then perl '.$self->o('sequence_dump_script').' -dbhost '.$self->o('dna_db_server').' -dbuser '.$self->o('dna_db_user').' -dbport '.$self->o('dna_db_port').' -dbname '.$self->o('dna_db_name').' -coord_system_name '.$self->o('assembly_name').' -toplevel -onefile -header rnaseq -filename '.$self->o('faidx_genome_file').';fi',
      },
      -flow_into => {
        1 => [ 'create_faidx'],
      },
    },

    {
      -logic_name => 'create_faidx',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => '5GB',
      -parameters => {
        cmd => 'if [ ! -e "'.$self->o('faidx_genome_file').'.fai" ]; then '.$self->o('samtools_path').' faidx '.$self->o('faidx_genome_file').';fi',
      },
    },
  ];
}


sub resource_classes {
  my $self = shift;

  return {
    '1GB' => { LSF => $self->lsf_resource_builder('production', 1000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '2GB' => { LSF => $self->lsf_resource_builder('production', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB' => { LSF => $self->lsf_resource_builder('production', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB' => { LSF => $self->lsf_resource_builder('production', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '8GB' => { LSF => $self->lsf_resource_builder('production', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default' => { LSF => $self->lsf_resource_builder('production', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
  }
}


1;
