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

package Bio::EnsEMBL::Analysis::Hive::Config::Genome_annotation_conf;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
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
    dbowner                          => '',
    pipeline_name                    => '' || $self->o('production_name').'_'.$self->o('release_number'),
    user_r                           => '', # read only db user
    user                             => '', # write db user
    password                         => '', # password for write db user
    server_set                       => '', # What server set to user, e.g. set1
    pipe_db_host                     => '', # host for pipe db
    databases_host                 => '', # host for general output dbs
    dna_db_host                      => '', # host for dna db
    pipe_db_port                     => '', # port for pipeline host
    databases_port                   => '', # port for general output db host
    dna_db_port                      => '', # port for dna db host
                                    
    registry_host                    => '', # host for registry db
    registry_port                    => '', # port for registry db
    registry_db                      => '', # name for registry db
                                    
    release_number                   => '' || $ENV{ENSEMBL_RELEASE},
    species_name                     => '', # e.g. mus_musculus
    production_name                  => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    annotation_source                => 'ensembl', # Used to define the structure of the ftp, default value ensembl 
    dbname_accession                 => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1
    taxon_id                         => '', # should be in the assembly report file
    species_taxon_id                 => '' || $self->o('taxon_id'), # Species level id, could be different to taxon_id if we have a subspecies, used to get species level RNA-seq CSV data
    genus_taxon_id                   => '' || $self->o('taxon_id'), # Genus level taxon id, used to get a genus level csv file in case there is not enough species level transcriptomic data
    uniprot_set                      => '', # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    sanity_set                       => '', #sanity checks
    ig_tr_fasta_file                 => '', # file containing ig and tr proteins to be used during the IGTR subpipeline. This would come from the clade settings defined in "create_annotation_configs.pl" (ie 'fish_ig_tr.fa')
    output_path                      => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    wgs_id                           => '', # Can be found in assembly report file on https://ftp.ncbi.nlm.nih.gov/genomes/genbank/
    assembly_name                    => '', # Name (as it appears in the assembly report file)
    assembly_accession               => '', # Versioned GCA assembly accession, e.g. GCA_001857705.1
    assembly_refseq_accession        => '', # Versioned GCF accession, e.g. GCF_001857705.1
    registry_file                    => '', # Path to databse registry for LastaZ and Production sync
    use_genome_flatfile              => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    repeatmasker_slice_size          => '1000000',# This is the default value for creating repeatmasker slice sizes
    batch_target_size                => '500000',# This is the default value for batching repeatmasker slice jobs
    species_url                      => '', # sets species.url meta key
    species_division                 => 'EnsemblVertebrates', # sets species.division meta key
    is_non_vert                      => '0', # Setting this will indicate that the assembly corresponds to a non-vertebrate species.
    protein_blast_db_file            => '', # use PE12 for non-vertebrates. Note there must also be a PE12_index file available in the same directory.
    protein_entry_loc_file           => 'entry_loc',

    repbase_logic_name               => '', # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    repbase_library                  => '', # repbase library name, this is the actual repeat repbase library to use, e.g. "Mus musculus"
    repeatmodeler_library            => '', # This should be the path to a custom repeat library, leave blank if none exists

    stable_id_start                  => '0', # When mapping is not required this is usually set to 0
    stable_id_prefix                 => '', # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    mapping_required                 => '0', # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    mapping_db                       => 'undef', # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now

    paired_end_only                  => '1', # Will only use paired-end rnaseq data if 1
    rnaseq_study_accession           => '', # A study accession for a transcriptomic dataset, if provided, only this data will be used
    long_read_study_accession        => '', # A study accession for a transcriptomic dataset, if provided, only this data will be used
    rnaseq_summary_file              => '', # Set this if you have a pre-existing cvs file with the expected columns
    rnaseq_summary_file_genus        => '', # Set this if you have a pre-existing genus level cvs file with the expected columns
    long_read_summary_file           => '', # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    long_read_summary_file_genus     => '', # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    long_read_fastq_dir              => '',
	
    skip_repeatmodeler               => '', # Skip using our repeatmodeler library for the species with repeatmasker, will still run standard repeatmasker
    skip_post_repeat_analyses        => '0', # Will skip everything after the repreats (rm, dust, trf) in the genome prep phase if 1, i.e. skips cpg, eponine, genscan, genscan blasts etc.
    skip_projection                  => '0', # Will skip projection process if 1
    skip_lastz                       => '0', # Will skip lastz if 1 (if skip_projection is enabled this is irrelevant)
    skip_rnaseq                      => '0', # Will skip rnaseq analyses if 1
    skip_long_read                   => '0', # Will skip long read analyses if 1
    skip_ncrna                       => '0', # Will skip ncrna process if 1
    skip_cleaning                    => '0', # Will skip the cleaning phase, will keep more genes/transcripts but some lower quality models may be kept

    # Keys for custom loading, only set/modify if that's what you're doing
    skip_genscan_blasts              => '1',
    load_toplevel_only               => '1', # This will not load the assembly info and will instead take any chromosomes, unplaced and unlocalised scaffolds directly in the DNA table
    custom_toplevel_file_path        => '', # Only set this if you are loading a custom toplevel, requires load_toplevel_only to also be set to 2

    email_address					 => '' || $ENV{USER}.'@ebi.ac.uk',

###############################
# Sub pipeline configurations
###############################
    hive_load_assembly_config => 'Bio::EnsEMBL::Analysis::Hive::Config::LoadAssembly',
    hive_refseq_import_config => 'Bio::EnsEMBL::Analysis::Hive::Config::Refseq_import_subpipeline',
    hive_repeat_masking_config => 'Bio::EnsEMBL::Analysis::Hive::Config::RepeatMasking',
    hive_lastz_config => 'Bio::EnsEMBL::Analysis::Hive::Config::LastZ',
    hive_projection_config => 'Bio::EnsEMBL::Analysis::Hive::Config::Projection_subpipeline_conf',
    hive_homology_config => 'Bio::EnsEMBL::Analysis::Hive::Config::GenblastHomology',
    hive_best_targeted_config => 'Bio::EnsEMBL::Analysis::Hive::Config::BestTargetted_subpipeline',
    hive_igtr_config => 'Bio::EnsEMBL::Analysis::Hive::Config::IGTR_subpipeline',
    hive_short_ncrna_config => 'Bio::EnsEMBL::Analysis::Hive::Config::ShortncRNA',
    hive_rnaseq_config => 'Bio::EnsEMBL::Analysis::Hive::Config::StarScallopRnaseq',
    hive_long_read_config => 'Bio::EnsEMBL::Analysis::Hive::Config::long_read',
    hive_homology_rnaseq_config => 'Bio::EnsEMBL::Analysis::Hive::Config::HomologyRnaseq',
    hive_transcript_finalisation_config => 'Bio::EnsEMBL::Analysis::Hive::Config::TranscriptSelection',
    hive_core_db_finalisation_config => 'Bio::EnsEMBL::Analysis::Hive::Config::FinaliseCoreDB',
    hive_otherfeatures_db_config => 'Bio::EnsEMBL::Analysis::Hive::Config::OtherFeatureDb',
    hive_rnaseq_db_config => 'Bio::EnsEMBL::Analysis::Hive::Config::production_rnaseq_db',


########################
# Pipe and ref db info
########################
    pipe_db_name                  => '',
    dna_db_name                   => $self->o('dbowner').'_'.$self->o('dbname_accession').'_core_'.$self->o('release_number'),

    reference_db_name            => $self->o('dna_db_name'),
    reference_db_host            => $self->o('dna_db_host'),
    reference_db_port            => $self->o('dna_db_port'),

    refseq_db_host => $self->o('databases_host'),
    refseq_db_port => $self->o('databases_port'),

    projection_source_db_name    => '', # This is generally a pre-existing db, like the current human/mouse core for example
    projection_source_db_host    => '', 
    projection_source_db_port    => '',
    projection_source_production_name => '',

    compara_db_name => 'gb_ensembl_compara_112',
    compara_db_host => 'mysql-ens-genebuild-prod-1',
    compara_db_port => 4527,

    projection_db_host    => $self->o('databases_host'),
    projection_db_port    => $self->o('databases_port'),

    genblast_db_host             => $self->o('databases_host'),
    genblast_db_port             => $self->o('databases_port'),

    genblast_nr_db_host             => $self->o('databases_host'),
    genblast_nr_db_port             => $self->o('databases_port'),

    cdna_db_host                 => $self->o('databases_host'),
    cdna_db_port                 => $self->o('databases_port'),

    best_targeted_db_host                 => $self->o('databases_host'),
    best_targeted_db_port                 => $self->o('databases_port'),

    ig_tr_db_host                => $self->o('databases_host'),
    ig_tr_db_port                => $self->o('databases_port'),

    ncrna_db_host                => $self->o('databases_host'),
    ncrna_db_port                => $self->o('databases_port'),

    rnaseq_db_host               => $self->o('databases_host'),
    rnaseq_db_port               => $self->o('databases_port'),

    rnaseq_refine_db_host         => $self->o('databases_host'),
    rnaseq_refine_db_port         => $self->o('databases_port'),

    rnaseq_blast_db_host         => $self->o('databases_host'),
    rnaseq_blast_db_port         => $self->o('databases_port'),

    long_read_initial_db_host    => $self->o('databases_host'),
    long_read_initial_db_port    => $self->o('databases_port'),

    long_read_final_db_host      => $self->o('databases_host'),
    long_read_final_db_port      => $self->o('databases_port'),

    rnaseq_for_layer_db_host     => $self->o('databases_host'),
    rnaseq_for_layer_db_port     => $self->o('databases_port'),

    rnaseq_for_layer_nr_db_host     => $self->o('databases_host'),
    rnaseq_for_layer_nr_db_port     => $self->o('databases_port'),

    pcp_db_host                     => $self->o('databases_host'),
    pcp_db_port                     => $self->o('databases_port'),

    genblast_rnaseq_support_db_host    => $self->o('databases_host'),
    genblast_rnaseq_support_db_port    => $self->o('databases_port'),

    # Layering is one of the most intesnive steps, so separating it off the main output server helps
    # Have also set module to use flatfile seq retrieval, so even if it's on the same server as the
    # core, the core should not be accessed
    layering_db_host             => $self->o('dna_db_host'),
    layering_db_port             => $self->o('dna_db_port'),

    utr_db_host                  => $self->o('databases_host'),
    utr_db_port                  => $self->o('databases_port'),

    genebuilder_db_host         => $self->o('databases_host'),
    genebuilder_db_port          => $self->o('databases_port'),

    pseudogene_db_host           => $self->o('databases_host'),
    pseudogene_db_port           => $self->o('databases_port'),

    final_geneset_db_host        => $self->o('databases_host'),
    final_geneset_db_port        => $self->o('databases_port'),

    killlist_db_host             => $self->o('databases_host'),
    killlist_db_port             => $self->o('databases_port'),

    otherfeatures_db_host        => $self->o('databases_host'),
    otherfeatures_db_port        => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    production_db_host           => 'mysql-ens-meta-prod-1',
    production_db_port           => '4483',

    taxonomy_db_host        => $self->o('production_db_host'),
    taxonomy_db_port        => $self->o('production_db_port'),

    databases_to_delete => [],


######################################################
#
# Mostly constant settings
#
######################################################
    assembly_provider_name   => '',
    assembly_provider_url    => '',
    annotation_provider_name => 'Ensembl',
    annotation_provider_url  => 'www.ensembl.org',

    full_repbase_logic_name  => "repeatmask_repbase_".$self->o('repbase_logic_name'),
    red_logic_name           => 'repeatdetector', # logic name for the Red repeat finding analysis
    repeatmodeler_logic_name => 'repeatmask_repeatmodeler',
    first_choice_repeat => $self->o('full_repbase_logic_name'),
    second_choice_repeat => $self->o('repeatmodeler_logic_name'),
    third_choice_repeat => $self->o('red_logic_name'),

    ensembl_analysis_script => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    loading_report_script   => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'report_genome_prep_stats.pl'),
    ensembl_misc_script     => catdir($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts'),
    repeat_types_script     => catfile($self->o('ensembl_misc_script'), 'repeats', 'repeat-types.pl'),
    registry_status_update_script     => catfile($self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'info_from_registry', 'update_assembly_registry.py'),

    hive_beekeeper_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'beekeeper.pl'),

    rnaseq_dir    => '',
    long_read_dir => '',
    gst_dir       => '',
	
    blast_type => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    cdna_threshold => 10, # The lowest number of genes using the cdna_db in the otherfeatures_db

########################
# Extra db settings
########################
    mysql_dump_options => '--max_allowed_packet=1000MB --quick',


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
########################
# db info
########################
    reference_db => {
      -dbname => $self->o('reference_db_name'),
      -host   => $self->o('reference_db_host'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    refseq_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_refseq_'.$self->o('release_number'),
      -host   => $self->o('refseq_db_host'),
      -port   => $self->o('refseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    compara_db => {
      -dbname => $self->o('compara_db_name'),
      -host   => $self->o('compara_db_host'),
      -port   => $self->o('compara_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    projection_source_db => {
      -dbname => $self->o('projection_source_db_name'),
      -host   => $self->o('projection_source_db_host'),
      -port   => $self->o('projection_source_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    selected_projection_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_sel_proj_'.$self->o('release_number'),
      -host   => $self->o('projection_db_host'),
      -port   => $self->o('projection_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    genblast_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_genblast_'.$self->o('release_number'),
      -host   => $self->o('genblast_db_host'),
      -port   => $self->o('genblast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    genblast_nr_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_genblast_nr_'.$self->o('release_number'),
      -host   => $self->o('genblast_nr_db_host'),
      -port   => $self->o('genblast_nr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    cdna_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_cdna_'.$self->o('release_number'),
      -host   => $self->o('cdna_db_host'),
      -port   => $self->o('cdna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    best_targeted_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_bt_'.$self->o('release_number'),
      -host   => $self->o('best_targeted_db_host'),
      -port   => $self->o('best_targeted_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    ig_tr_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_igtr_'.$self->o('release_number'),
      -host   => $self->o('ig_tr_db_host'),
      -port   => $self->o('ig_tr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    ncrna_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_ncrna_'.$self->o('release_number'),
      -host   => $self->o('ncrna_db_host'),
      -port   => $self->o('ncrna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    rnaseq_refine_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_refine_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_refine_db_host'),
      -port   => $self->o('rnaseq_refine_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    rnaseq_blast_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_scallop_blast_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_blast_db_host'),
      -port   => $self->o('rnaseq_blast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    rnaseq_for_layer_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_rnaseq_layer_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_db_host'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    rnaseq_for_layer_nr_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_rnalayer_nr_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_nr_db_host'),
      -port   => $self->o('rnaseq_for_layer_nr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    pcp_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_pcp_'.$self->o('release_number'),
      -host   => $self->o('pcp_db_host'),
      -port   => $self->o('pcp_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    pcp_nr_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_pcp_nr_'.$self->o('release_number'),
      -host   => $self->o('pcp_db_host'),
      -port   => $self->o('pcp_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    long_read_final_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_lrfinal_'.$self->o('release_number'),
      -host => $self->o('long_read_final_db_host'),
      -port => $self->o('long_read_final_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    genblast_rnaseq_support_nr_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_gb_rnaseq_nr_'.$self->o('release_number'),
      -host   => $self->o('genblast_rnaseq_support_db_host'),
      -port   => $self->o('genblast_rnaseq_support_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    final_geneset_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_final_'.$self->o('release_number'),
      -host   => $self->o('final_geneset_db_host'),
      -port   => $self->o('final_geneset_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    rnaseq_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_rnaseq_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_db_host'),
      -port   => $self->o('rnaseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    otherfeatures_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_otherfeatures_'.$self->o('release_number'),
      -host   => $self->o('otherfeatures_db_host'),
      -port   => $self->o('otherfeatures_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    killlist_db => {
      -dbname => $self->o('killlist_db_name'),
      -host   => $self->o('killlist_db_host'),
      -port   => $self->o('killlist_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    production_db => {
      -host   => $self->o('production_db_host'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ensembl_production',
      -driver => $self->o('hive_driver'),
    },

    taxonomy_db => {
      -host   => $self->o('production_db_host'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ncbi_taxonomy',
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_create_commands {
    my ($self) = @_;

    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      'mkdir -p '.$self->o('rnaseq_dir'),
      'mkdir -p '.$self->o('long_read_fastq_dir'),
      'mkdir -p '.$self->o('gst_dir'),
    ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    skip_projection => $self->o('skip_projection'),
    skip_rnaseq => $self->o('skip_rnaseq'),
    skip_ncrna => $self->o('skip_ncrna'),
    skip_long_read => $self->o('skip_long_read'),
    skip_lastz => $self->o('skip_lastz'),
    skip_repeatmodeler => $self->o('skip_repeatmodeler'),
    skip_post_repeat_analyses => $self->o('skip_post_repeat_analyses'),
  }
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  my $db_name_template = "%s_%s_%s_%s";
  my ($load_assembly_pipe_db, $load_assembly_pipe_url, $load_assembly_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'load_assembly_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($repeat_masking_pipe_db, $repeat_masking_pipe_url, $repeat_masking_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'repeat_masking_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($refseq_import_pipe_db, $refseq_import_pipe_url, $refseq_import_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'refseq_import_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($lastz_pipe_db, $lastz_pipe_url, $lastz_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'lastz_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($projection_pipe_db, $projection_pipe_url, $projection_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'projection_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($homology_pipe_db, $homology_pipe_url, $homology_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'homology_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($best_targeted_pipe_db, $best_targeted_pipe_url, $best_targeted_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'best_targeted_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($igtr_pipe_db, $igtr_pipe_url, $igtr_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'igtr_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($short_ncrna_pipe_db, $short_ncrna_pipe_url, $short_ncrna_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'short_ncrna_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($rnaseq_pipe_db, $rnaseq_pipe_url, $rnaseq_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'rnaseq_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($long_read_pipe_db, $long_read_pipe_url, $long_read_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'long_read_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($homology_rnaseq_pipe_db, $homology_rnaseq_pipe_url, $homology_rnaseq_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'homology_rnaseq_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($transcript_selection_pipe_db, $transcript_selection_pipe_url, $transcript_finalisation_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'set_finalisation_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($core_db_finalisation_pipe_db, $core_db_finalisation_pipe_url, $core_db_finalisation_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'db_finalisation_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($otherfeatures_db_pipe_db, $otherfeatures_db_pipe_url, $otherfeatures_db_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'otherfeatures_db_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );
  my ($rnaseq_db_pipe_db, $rnaseq_db_pipe_url, $rnaseq_db_guihive) = $self->get_meta_db_information(
      undef,
      sprintf($db_name_template, $self->o('dbowner'), $self->o('dbname_accession'), 'rnaseq_db_pipe', $self->o('release_number')),
      $self->o('pipe_db_host'),
      $self->o('pipe_db_port'),
      $self->o('pipe_db_user'),
      $self->o('pipe_db_pass'),
    );


  return [
    {
      -logic_name => 'download_rnaseq_csv_from_transcriptomic_registry',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '1GB',
      -parameters => {
        cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes','transcriptomic_data','select_transcriptomic_data.py' ) . ' --taxon_id  ' . $self->o('species_taxon_id') .' --file_name ' . $self->o('rnaseq_summary_file') . ' --csv_for_main --limit 200',
      },
      -flow_into => {
        1 => ['check_rnaseq_summary_file'],
      },
      -input_ids  => [
        {
          assembly_name => $self->o('assembly_name'),
          assembly_accession => $self->o('assembly_accession'),
          assembly_refseq_accession => $self->o('assembly_refseq_accession'),
          old_core_dbname => '',
        },
      ],
    },
    {
      -logic_name => 'check_rnaseq_summary_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -s "'.$self->o('rnaseq_summary_file').'" ] ; then exit 0; else  exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -flow_into => {
        1 => ['create_load_assembly_pipeline_job'],
        2 => ['download_rnaseq_csv'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'download_rnaseq_csv',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => '1GB',
      -parameters => {
        cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes','transcriptomic_data','get_transcriptomic_data.py' ) . ' -t  ' . $self->o('species_taxon_id') .' -f ' . $self->o('rnaseq_summary_file') . ' --read_type short --tree -l 250' ,
      },
      -flow_into => {
        1 => ['create_load_assembly_pipeline_job'],
      },
    },

    {
      -logic_name => 'create_load_assembly_pipeline_job',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => ['meta_pipeline_db', 'ehive_url', 'pipeline_name', 'external_link', 'old_core_dbname'],
        inputlist => [[$load_assembly_pipe_db, $load_assembly_pipe_url, 'load_assembly_'.$self->o('production_name'), $load_assembly_guihive, '#old_core_dbname#']],
        guihive_host => $self->o('guihive_host'),
        guihive_port => $self->o('guihive_port'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        '2->A' => ['initialise_load_assembly'],
        'A->1' => ['create_repeat_masking_pipeline_jobs']
      }
    },

    {
      -logic_name => 'initialise_load_assembly',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit',
      -parameters => {
        hive_config => $self->o('hive_load_assembly_config'),
        databases => ['reference_db', 'dna_db'],
        reference_db => $self->o('reference_db'),
        dna_db => $self->o('dna_db'),
        enscode_root_dir => $self->o('enscode_root_dir'),
        extra_parameters => {
          output_path => $self->o('output_path'),
          user_r => $self->o('user_r'),
          dna_db_host => $self->o('dna_db_host'),
          dna_db_port => $self->o('dna_db_port'),
          release_number => $self->o('release_number'),
          species_name => $self->o('species_name'),
          production_name => $self->o('production_name'),
          dbname_accession => $self->o('dbname_accession'),
          taxon_id => $self->o('taxon_id'),
          wgs_id => $self->o('wgs_id'),
          assembly_name => $self->o('assembly_name'),
          assembly_accession => $self->o('assembly_accession'),
          stable_id_prefix => $self->o('stable_id_prefix'),
          species_url => $self->o('species_url'),
          load_toplevel_only => $self->o('load_toplevel_only'),
          custom_toplevel_file_path => $self->o('custom_toplevel_file_path'),
          old_core_dbname => '#old_core_dbname#',
        },
      },
      -rc_name      => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['run_load_assembly'],
      },
    },

    {
      -logic_name => 'run_load_assembly',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun',
      -parameters => {
        beekeeper_script => $self->o('hive_beekeeper_script'),
      },
      -rc_name      => 'default',
      -max_retry_count => 1,
    },

    {
      -logic_name => 'create_repeat_masking_pipeline_jobs',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => ['meta_pipeline_db', 'ehive_url', 'pipeline_name', 'external_link'],
        inputlist => [[$repeat_masking_pipe_db, $repeat_masking_pipe_url, 'repeat_masking_'.$self->o('production_name'), $repeat_masking_guihive]],
        guihive_host => $self->o('guihive_host'),
        guihive_port => $self->o('guihive_port'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        '2->A' => ['initialise_repeat_masking'],
        'A->1' => ['create_transcript_selection_pipeline_job'],
      }
    },

    {
      -logic_name => 'initialise_repeat_masking',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit',
      -parameters => {
        hive_config => $self->o('hive_repeat_masking_config'),
        databases => ['reference_db', 'dna_db'],
        reference_db => $self->o('reference_db'),
        dna_db => $self->o('dna_db'),
        enscode_root_dir => $self->o('enscode_root_dir'),
        extra_parameters => {
          output_path => $self->o('output_path'),
          user_r => $self->o('user_r'),
          dna_db_host => $self->o('dna_db_host'),
          dna_db_port => $self->o('dna_db_port'),
          release_number => $self->o('release_number'),
          production_name => $self->o('production_name'),
          dbname_accession => $self->o('dbname_accession'),
          repbase_logic_name => $self->o('repbase_logic_name'),
          repbase_library => $self->o('repbase_library'),
          species_name => $self->o('species_name'),
          use_genome_flatfile => $self->o('use_genome_flatfile'),
          skip_repeatmodeler => $self->o('skip_repeatmodeler'),
          red_logic_name => $self->o('red_logic_name'),
          repeatmodeler_library => $self->o('repeatmodeler_library'),
          skip_post_repeat_analyses => $self->o('skip_post_repeat_analyses'),
          batch_target_size => $self->o('batch_target_size'),
          repeatmasker_slice_size => $self->o('repeatmasker_slice_size'),
        },
      },
      -rc_name      => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['run_repeat_masking'],
      },
    },

    {
      -logic_name => 'run_repeat_masking',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun',
      -parameters => {
        beekeeper_script => $self->o('hive_beekeeper_script'),
      },
      -rc_name      => 'default',
      -max_retry_count => 1,
    },

    {
      -logic_name => 'create_transcript_selection_pipeline_job',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => ['meta_pipeline_db', 'ehive_url', 'pipeline_name', 'external_link'],
        inputlist => [[$transcript_selection_pipe_db, $transcript_selection_pipe_url, 'transcript_finalisation_'.$self->o('production_name'), $transcript_finalisation_guihive]],
        guihive_host => $self->o('guihive_host'),
        guihive_port => $self->o('guihive_port'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        '2->A' => ['create_rnaseq_pipeline_job'],
        'A->1' => ['create_rnaseq_db_pipeline_job']
      }
    },

    {
      -logic_name => 'create_rnaseq_pipeline_job',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => ['meta_pipeline_db', 'ehive_url', 'pipeline_name', 'external_link'],
        inputlist => [[$rnaseq_pipe_db, $rnaseq_pipe_url, 'rnaseq_'.$self->o('production_name'), $rnaseq_guihive]],
        guihive_host => $self->o('guihive_host'),
        guihive_port => $self->o('guihive_port'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        2 => ['initialise_rnaseq'],
      }
    },

    {
      -logic_name => 'initialise_rnaseq',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit',
      -parameters => {
        hive_config => $self->o('hive_rnaseq_config'),
        databases => ['rnaseq_refine_db', 'rnaseq_for_layer_nr_db', 'rnaseq_for_layer_db', 'dna_db', 'pcp_db', 'pcp_nr_db'],
        rnaseq_refine_db => $self->o('rnaseq_refine_db'),
        rnaseq_for_layer_db => $self->o('rnaseq_for_layer_db'),
        rnaseq_for_layer_nr_db => $self->o('rnaseq_for_layer_nr_db'),
        dna_db => $self->o('dna_db'),
	      pcp_db => $self->o('pcp_db'),
	      pcp_nr_db => $self->o('pcp_nr_db'),
        enscode_root_dir => $self->o('enscode_root_dir'),
        extra_parameters => {
          output_path => $self->o('output_path'),
          user_r => $self->o('user_r'),
          dna_db_host => $self->o('dna_db_host'),
          dna_db_port => $self->o('dna_db_port'),
          databases_host => $self->o('databases_host'),
          databases_port => $self->o('databases_port'),
          release_number => $self->o('release_number'),
          assembly_name => $self->o('assembly_name'),
          production_name => $self->o('production_name'),
          dbname_accession => $self->o('dbname_accession'),
          species_name => $self->o('species_name'),
          taxon_id => $self->o('taxon_id'),
          uniprot_set => $self->o('uniprot_set'),
          sanity_set => $self->o('sanity_set'),
          use_genome_flatfile => $self->o('use_genome_flatfile'),
          main_pipeline_url => $self->pipeline_url,
          transcript_selection_url => $transcript_selection_pipe_url,
          homology_rnaseq_url => $homology_rnaseq_pipe_url,
          is_non_vert => $self->o('is_non_vert'),
          protein_blast_db_file => $self->o('protein_blast_db_file'),
          protein_entry_loc_file => $self->o('protein_entry_loc_file'),
        },
      },
      -rc_name      => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['run_rnaseq'],
      },
    },

    {
      -logic_name => 'run_rnaseq',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun',
      -parameters => {
        beekeeper_script => $self->o('hive_beekeeper_script'),
      },
      -rc_name      => 'default',
      -max_retry_count => 1,
    },
    {
      -logic_name => 'create_rnaseq_db_pipeline_job',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => ['meta_pipeline_db', 'ehive_url', 'pipeline_name', 'external_link'],
        inputlist => [[$rnaseq_db_pipe_db, $rnaseq_db_pipe_url, 'rnaseq_db_'.$self->o('production_name'), $rnaseq_db_guihive]],
        guihive_host => $self->o('guihive_host'),
        guihive_port => $self->o('guihive_port'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        2 => ['initialise_rnaseq_db'],
      }
    },

    {
      -logic_name => 'initialise_rnaseq_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit',
      -parameters => {
        hive_config => $self->o('hive_rnaseq_db_config'),
        databases => ['rnaseq_db', 'rnaseq_refine_db', 'rnaseq_blast_db', 'dna_db'],
        rnaseq_db => $self->o('rnaseq_db'),
        rnaseq_refine_db => $self->o('rnaseq_refine_db'),
        rnaseq_blast_db => $self->o('rnaseq_blast_db'),
        dna_db => $self->o('dna_db'),
		ehive_root_dir => $ENV{EHIVE_ROOT_DIR},
		enscode_root_dir => $self->o('enscode_root_dir'),
        extra_parameters => {
          output_path => $self->o('output_path'),
          user_r => $self->o('user_r'),
          dna_db_host => $self->o('dna_db_host'),
          dna_db_port => $self->o('dna_db_port'),
          databases_host => $self->o('databases_host'),
          databases_port => $self->o('databases_port'),
          release_number => $self->o('release_number'),
          production_name => $self->o('production_name'),
          dbname_accession => $self->o('dbname_accession'),
          species_name => $self->o('species_name'),
          registry_host => $self->o('registry_host'),
          registry_port => $self->o('registry_port'),
          registry_db => $self->o('registry_db'),
          assembly_name => $self->o('assembly_name'),
          assembly_accession => $self->o('assembly_accession'),
          annotation_source => $self->o('annotation_source'),
          sanity_set => $self->o('sanity_set'),
        },
      },
      -rc_name      => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['run_rnaseq_db'],
      },
    },

    {
      -logic_name => 'run_rnaseq_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun',
      -parameters => {
        beekeeper_script => $self->o('hive_beekeeper_script'),
      },
      -rc_name      => 'default',
      -max_retry_count => 1,
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

