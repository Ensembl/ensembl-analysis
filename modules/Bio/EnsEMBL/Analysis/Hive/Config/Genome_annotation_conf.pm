=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
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
    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'             => '' || $self->o('production_name').$self->o('production_name_modifier').'_'.$self->o('ensembl_release'),
    'user_r'                    => '', # read only db user
    'user'                      => '', # write db user
    'password'                  => '', # password for write db user
    'server_set'                => '', # What server set to user, e.g. set1
    'pipe_db_server'            => '', # host for pipe db
    'databases_server'          => '', # host for general output dbs
    'dna_db_server'             => '', # host for dna db
    'pipe_db_port'              => '', # port for pipeline host
    'databases_port'            => '', # port for general output db host
    'dna_db_port'               => '', # port for dna db host
    'registry_host'             => '', # host for registry db
    'registry_port'             => '', # port for registry db
    'registry_db'               => '', # name for registry db
    'repbase_logic_name'        => '', # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    'repbase_library'           => '', # repbase library name, this is the actual repeat repbase library to use, e.g. "Mus musculus"
    'rnaseq_summary_file'       => '' || catfile($self->o('rnaseq_dir'), $self->o('species_name').'.csv'), # Set this if you have a pre-existing cvs file with the expected columns
    'rnaseq_summary_file_genus' => '' || catfile($self->o('rnaseq_dir'), $self->o('species_name').'_gen.csv'), # Set this if you have a pre-existing genus level cvs file with the expected columns
    'long_read_summary_file'    => '' || catfile($self->o('long_read_dir'), $self->o('species_name').'_long_read.csv'), # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_summary_file_genus' => '' || catfile($self->o('long_read_dir'), $self->o('species_name').'_long_read_gen.csv'), # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_fastq_dir'       => '' || catdir($self->o('long_read_dir'),'input'),
    'release_number'            => '' || $self->o('ensembl_release'),
    'species_name'              => '', # e.g. mus_musculus
    'production_name'           => '' || $self->o('species_name'), # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'                  => '', # should be in the assembly report file
    'species_taxon_id'          => '' || $self->o('taxon_id'), # Species level id, could be different to taxon_id if we have a subspecies, used to get species level RNA-seq CSV data
    'genus_taxon_id'            => '' || $self->o('taxon_id'), # Genus level taxon id, used to get a genus level csv file in case there is not enough species level transcriptomic data
    'uniprot_set'               => '', # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'wgs_id'                    => '', # Can be found in assembly report file on ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
    'assembly_name'             => '', # Name (as it appears in the assembly report file)
    'assembly_accession'        => '', # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'assembly_refseq_accession' => '', # Versioned GCF accession, e.g. GCF_001857705.1
    'stable_id_prefix'          => '', # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    'use_genome_flatfile'       => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    'species_url'               => '' || $self->o('production_name').$self->o('production_name_modifier'), # sets species.url meta key
    'species_division'          => 'EnsemblVertebrates', # sets species.division meta key
    'stable_id_start'           => '0', # When mapping is not required this is usually set to 0
    'skip_repeatmodeler'        => '0', # Skip using our repeatmodeler library for the species with repeatmasker, will still run standard repeatmasker
    'skip_post_repeat_analyses' => '0', # Will everything after the repreats (rm, dust, trf) in the genome prep phase if 1, i.e. skips cpg, eponine, genscan, genscan blasts etc.
    'skip_projection'           => '0', # Will skip projection process if 1
    'skip_lastz'                => '0', # Will skip lastz if 1 (if skip_projection is enabled this is irrelevant)
    'skip_rnaseq'               => '0', # Will skip rnaseq analyses if 1
    'skip_long_read'            => '1', # Will skip long read analyses if 1
    'skip_ncrna'                => '0', # Will skip ncrna process if 1
    'skip_cleaning'             => '0', # Will skip the cleaning phase, will keep more genes/transcripts but some lower quality models may be kept
    'mapping_required'          => '0', # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'                => '', # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now
    'uniprot_version'           => 'uniprot_2018_07', # What UniProt data dir to use for various analyses
    'vertrna_version'           => '136', # The version of VertRNA to use, should correspond to a numbered dir in VertRNA dir

    'ig_tr_fasta_file'          => 'human_ig_tr.fa', # What IMGT fasta file to use. File should contain protein segments with appropriate headers
    'mt_accession'              => undef, # This should be set to undef unless you know what you are doing. If you specify an accession, then you need to add the parameters to the load_mitochondrion analysis
    'production_name_modifier'  => '', # Do not set unless working with non-reference strains, breeds etc. Must include _ in modifier, e.g. _hni for medaka strain HNI

    # Keys for custom loading, only set/modify if that's what you're doing
    'skip_genscan_blasts'          => '1',
    'load_toplevel_only'           => '1', # This will not load the assembly info and will instead take any chromosomes, unplaced and unlocalised scaffolds directly in the DNA table
    'custom_toplevel_file_path'    => '', # Only set this if you are loading a custom toplevel, requires load_toplevel_only to also be set to 2
    'repeatmodeler_library'        => '', # This should be the path to a custom repeat library, leave blank if none exists
    'use_repeatmodeler_to_mask'    => '0', # Setting this will include the repeatmodeler library in the masking process
    'protein_blast_db'             => '' || catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata'), # Blast database for comparing the final models to.
    'protein_blast_index'          => '' || catdir($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata_index'), # Indicate Index for the blast database.
    'protein_entry_loc'            => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc'), # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing


########################
## Small ncRNAs params
#########################
    'mirBase_fasta'             => 'all_mirnas.fa', # What mirBase file to use. It is currently best to use on with the most appropriate set for your species
    'rfc_scaler'                => 'filter_dafs_rfc_scaler_human.pkl',
    'rfc_model'                 => 'filter_dafs_rfc_model_human.pkl',

    # Clade-based filtering on rfam accessions
    # Rfam db details should stay constant but check periodically
    'rfam_user' => 'rfamro',
    'rfam_dbname' => 'Rfam',
    'rfam_host' => 'mysql-rfam-public.ebi.ac.uk',
    'rfam_port' => 4497,

    'rfam_path' => catfile($self->o('base_blast_db_path'), 'ncrna', 'Rfam_14.0'),
    'rfam_seeds' => $self->o('rfam_path') . "/Rfam.seed",
    'rfam_cm' => $self->o('rfam_path') . "/Rfam.cm",
    'filtered_rfam_cm' => $self->o('output_path') . '/Rfam.cm',
    'clade' => $self->o('repbase_logic_name'),


########################
# Pipe and ref db info
########################

    'projection_source_db_name'    => 'homo_sapiens_core_96_38', # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_server'  => 'mysql-ens-mirror-1',
    'projection_source_db_port'    => '4240',
    'projection_source_production_name' => 'homo_sapiens',

    'compara_db_name'     => 'leanne_ensembl_compara_95',
    'compara_db_server'  => 'mysql-ens-genebuild-prod-5',
    'compara_db_port'    => 4531,

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is
    'projection_lastz_db_name'     => $self->o('pipe_db_name'),
    'projection_lastz_db_server'   => $self->o('pipe_db_server'),
    'projection_lastz_db_port'     => $self->o('pipe_db_port'),

    'provider_name'                => 'Ensembl',
    'provider_url'                 => 'www.ensembl.org',

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_pipe_'.$self->o('release_number'),
    'dna_db_name'                   => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_core_'.$self->o('release_number'),

    'reference_db_name'            => $self->o('dna_db_name'),
    'reference_db_server'          => $self->o('dna_db_server'),
    'reference_db_port'            => $self->o('dna_db_port'),

    'cdna_db_server'               => $self->o('databases_server'),
    'cdna_db_port'                 => $self->o('databases_port'),

    cdna2genome_db_server          => $self->o('databases_server'),
    cdna2genome_db_port            => $self->o('databases_port'),

    'genblast_db_server'           => $self->o('databases_server'),
    'genblast_db_port'             => $self->o('databases_port'),

    'genblast_rnaseq_support_db_server'  => $self->o('databases_server'),
    'genblast_rnaseq_support_db_port'    => $self->o('databases_port'),

    'ig_tr_db_server'              => $self->o('databases_server'),
    'ig_tr_db_port'                => $self->o('databases_port'),

    'genewise_db_server'           => $self->o('databases_server'),
    'genewise_db_port'             => $self->o('databases_port'),

    'projection_db_server'  => $self->o('databases_server'),
    'projection_db_port'    => $self->o('databases_port'),

    'projection_realign_db_server' => $self->o('databases_server'),
    'projection_realign_db_port'   => $self->o('databases_port'),

    'projection_lincrna_db_server' => $self->o('databases_server'),
    'projection_lincrna_db_port'   => $self->o('databases_port'),

    'projection_pseudogene_db_server' => $self->o('databases_server'),
    'projection_pseudogene_db_port'   => $self->o('databases_port'),

    'long_read_initial_db_server'  => $self->o('databases_server'),
    'long_read_initial_db_port'    => $self->o('databases_port'),

    'long_read_blast_db_server'    => $self->o('databases_server'),
    'long_read_blast_db_port'      => $self->o('databases_port'),

    'long_read_collapse_db_server' => $self->o('databases_server'),
    'long_read_collapse_db_port'   => $self->o('databases_port'),

    'long_read_final_db_server'    => $self->o('databases_server'),
    'long_read_final_db_port'      => $self->o('databases_port'),

    'rnaseq_for_layer_db_server'   => $self->o('databases_server'),
    'rnaseq_for_layer_db_port'     => $self->o('databases_port'),

    'rnaseq_db_server'             => $self->o('databases_server'),
    'rnaseq_db_port'               => $self->o('databases_port'),

    'rnaseq_rough_db_server'       => $self->o('databases_server'),
    'rnaseq_rough_db_port'         => $self->o('databases_port'),

    'rnaseq_refine_db_server'       => $self->o('databases_server'),
    'rnaseq_refine_db_port'         => $self->o('databases_port'),

    'rnaseq_blast_db_server'       => $self->o('databases_server'),
    'rnaseq_blast_db_port'         => $self->o('databases_port'),

    'lincrna_db_server'            => $self->o('databases_server'),
    'lincrna_db_port'              => $self->o('databases_port'),

    # Layering is one of the most intesnive steps, so separating it off the main output server helps
    # Have also set module to use flatfile seq retrieval, so even if it's on the same server as the
    # core, the core should not be accessed
    'layering_db_server'           => $self->o('dna_db_server'),
    'layering_db_port'             => $self->o('dna_db_port'),

    'utr_db_server'                => $self->o('databases_server'),
    'utr_db_port'                  => $self->o('databases_port'),

    'genebuilder_db_server'        => $self->o('databases_server'),
    'genebuilder_db_port'          => $self->o('databases_port'),

    'pseudogene_db_server'         => $self->o('databases_server'),
    'pseudogene_db_port'           => $self->o('databases_port'),

    'ncrna_db_server'              => $self->o('databases_server'),
    'ncrna_db_port'                => $self->o('databases_port'),
    ncrna_db_name                  => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_ncrna_'.$self->o('release_number'),

    'final_geneset_db_server'      => $self->o('databases_server'),
    'final_geneset_db_port'        => $self->o('databases_port'),

    'refseq_db_server'             => $self->o('databases_server'),
    'refseq_db_port'               => $self->o('databases_port'),

    'killlist_db_server'           => $self->o('databases_server'),
    'killlist_db_port'             => $self->o('databases_port'),

    'otherfeatures_db_server'      => $self->o('databases_server'),
    'otherfeatures_db_port'        => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    'production_db_server'         => 'mysql-ens-meta-prod-1',
    'production_db_port'           => '4483',


    databases_to_delete => ['reference_db', 'cdna_db', 'genblast_db', 'genewise_db', 'projection_db', 'selected_projection_db', 'layering_db', 'utr_db', 'genebuilder_db', 'pseudogene_db', 'ncrna_db', 'final_geneset_db', 'refseq_db', 'cdna2genome_db', 'rnaseq_blast_db', 'rnaseq_refine_db', 'rnaseq_rough_db', 'lincrna_db', 'otherfeatures_db', 'rnaseq_db'],#, 'projection_realign_db'

########################
# BLAST db paths
########################
    'base_blast_db_path'        => $ENV{BLASTDB_DIR},
    'vertrna_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'vertrna', $self->o('vertrna_version'), 'embl_vertrna-1'),
    'unigene_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'unigene', 'unigene'),
    'ncrna_blast_path'          => catfile($self->o('base_blast_db_path'), 'ncrna', 'ncrna_2016_05'),
    'mirna_blast_path'          => catfile($self->o('base_blast_db_path'), 'ncrna', 'mirbase_22'),
    'ig_tr_blast_path'          => catfile($self->o('base_blast_db_path'), 'ig_tr_genes'),

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

    'layering_input_gene_dbs' => [
                                   $self->o('genblast_nr_db'),
                                   $self->o('genblast_rnaseq_support_nr_db'),
                                   $self->o('rnaseq_for_layer_nr_db'),
                                   $self->o('selected_projection_db'),
                                   $self->o('ig_tr_db'),
                                   $self->o('best_targeted_db'),
                                   $self->o('long_read_final_db'),
                                 ],


    utr_donor_dbs => [
      $self->o('cdna_db'),
      $self->o('rnaseq_for_layer_db'),
      $self->o('long_read_final_db'),
    ],

    utr_acceptor_dbs => [
      $self->o('layering_db'),
    ],

    'utr_biotype_priorities'  => {
                                   'rnaseq' => 2,
                                   'cdna' => 1,
                                 },

    'cleaning_blessed_biotypes' => {
                                     'pseudogene' => 1,
                                     'processed_pseudogene' => 1,
                                     'IG_C_gene' => 1,
                                     'IG_V_gene' => 1,
                                     'TR_C_gene' => 1,
                                     'TR_D_gene' => 1,
                                     'TR_V_gene' => 1,
                                     'lncRNA'    => 1,
                                   },

    'min_toplevel_slice_length'   => 250,

    'repeatmodeler_logic_name'    => 'repeatmask_repeatmodeler',
    'homology_models_path'        => catdir($self->o('output_path'),'homology_models'),

    ncrna_dir => catdir($self->o('output_path'), 'ncrna'),
    targetted_path => catdir($self->o('output_path'), 'targetted'),
    cdna_file      => catfile($self->o('targetted_path'), 'cdnas'),
    annotation_file => $self->o('cdna_file').'.annotation',

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
    frameshift_attrib_script   => catfile($self->o('ensembl_misc_script'), 'frameshift_transcript_attribs.pl'),
    select_canonical_script    => catfile($self->o('ensembl_misc_script'),'canonical_transcripts', 'select_canonical_transcripts.pl'),
    assembly_name_script       => catfile($self->o('ensembl_analysis_script'), 'update_assembly_name.pl'),

    rnaseq_daf_introns_file => catfile($self->o('output_dir'), 'rnaseq_daf_introns.dat'),

    # Genes biotypes to ignore from the final db when copying to core
    copy_biotypes_to_ignore => {
                                 'low_coverage' => 1,
                                 'CRISPR' => 1,
                               },
########################
# Extra db settings
########################

    'num_tokens' => 10,
    mysql_dump_options => '--max_allowed_packet=1000MB',

########################
# Executable paths
########################
    'minimap2_genome_index'  => $self->o('faidx_genome_file').'.mmi',
    'minimap2_path'          => '/hps/nobackup2/production/ensembl/fergal/coding/long_read_aligners/new_mm2/minimap2/minimap2',
    'paftools_path'          => '/hps/nobackup2/production/ensembl/fergal/coding/long_read_aligners/new_mm2/minimap2/misc/paftools.js',
    'minimap2_batch_size'    => '5000',

    'blast_type' => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'dust_path' => catfile($self->o('binary_base'), 'dustmasker'),
    'trf_path' => catfile($self->o('binary_base'), 'trf'),
    'eponine_java_path' => catfile($self->o('binary_base'), 'java'),
    'eponine_jar_path' => catfile($self->o('software_base_path'), 'opt', 'eponine', 'libexec', 'eponine-scan.jar'),
    'cpg_path' => catfile($self->o('binary_base'), 'cpg_lh'),
    'trnascan_path' => catfile($self->o('binary_base'), 'tRNAscan-SE'),
    'repeatmasker_path' => catfile($self->o('binary_base'), 'RepeatMasker'),
    'genscan_path' => catfile($self->o('binary_base'), 'genscan'),
    'genscan_matrix_path' => catfile($self->o('software_base_path'), 'share', 'HumanIso.smat'),
    'uniprot_blast_exe_path' => catfile($self->o('binary_base'), 'blastp'),
    'blastn_exe_path' => catfile($self->o('binary_base'), 'blastn'),
    'vertrna_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    'unigene_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    genewise_path => catfile($self->o('binary_base'), 'genewise'),
    'exonerate_path'         => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
    'cmsearch_exe_path'    => catfile($self->o('software_base_path'), 'bin', 'cmsearch'), # #'opt', 'infernal10', 'bin', 'cmsearch'),
    indicate_path  => catfile($self->o('binary_base'), 'indicate'),
    pmatch_path  => catfile($self->o('binary_base'), 'pmatch'),
    exonerate_annotation => catfile($self->o('binary_base'), 'exonerate'),
    samtools_path => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary
    picard_lib_jar => catfile($self->o('software_base_path'), 'Cellar', 'picard-tools', '2.6.0', 'libexec', 'picard.jar'), #You need to specify the full path to the picard library
    bwa_path => catfile($self->o('software_base_path'), 'opt', 'bwa-051mt', 'bin', 'bwa'), #You may need to specify the full path to the bwa binary
    refine_ccode_exe => catfile($self->o('binary_base'), 'RefineSolexaGenes'), #You may need to specify the full path to the RefineSolexaGenes binary
    interproscan_exe => catfile($self->o('binary_base'), 'interproscan.sh'),
    bedtools => catfile($self->o('binary_base'), 'bedtools'),
    bedGraphToBigWig => catfile($self->o('binary_base'), 'bedGraphToBigWig'),
    'cesar_path' => catdir($self->o('software_base_path'),'opt','cesar','bin'),

    'uniprot_genblast_batch_size' => 15,
    'uniprot_table_name'          => 'uniprot_sequences',

    'genblast_path'     => catfile($self->o('binary_base'), 'genblast'),
    'genblast_eval'     => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
    'genblast_cov'      => '0.5',
    'genblast_pid'      => '30',
    'genblast_max_rank' => '5',
    'genblast_flag_small_introns' => 1,
    'genblast_flag_subpar_models' => 0,

    'ig_tr_table_name'    => 'ig_tr_sequences',
    'ig_tr_genblast_cov'  => '0.8',
    'ig_tr_genblast_pid'  => '70',
    'ig_tr_genblast_eval' => '1',
    'ig_tr_genblast_max_rank' => '5',
    'ig_tr_batch_size'    => 10,

    'exonerate_cdna_pid' => '95', # Cut-off for percent id
    'exonerate_cdna_cov' => '50', # Cut-off for coverage

    'cdna_selection_pid' => '97', # Cut-off for percent id for selecting the cDNAs
    'cdna_selection_cov' => '90', # Cut-off for coverage for selecting the cDNAs

# Best targetted stuff
    exonerate_logic_name => 'exonerate',
    ncbi_query => '((txid'.$self->o('taxon_id').'[Organism:noexp]+AND+biomol_mrna[PROP]))  NOT "tsa"[Properties]', 

    cdna_table_name    => 'cdna_sequences',
    target_exonerate_calculate_coverage_and_pid => 0,
    exonerate_protein_pid => 95,
    exonerate_protein_cov => 50,
    cdna2genome_region_padding => 2000,
    exonerate_max_intron => 200000,

    best_targetted_min_coverage => 50, # This is to avoid having models based on fragment alignment and low identity
    best_targetted_min_identity => 50, # This is to avoid having models based on fragment alignment and low identity


# RNA-seq pipeline stuff
    # You have the choice between:
    #  * using a csv file you already created
    #  * using a study_accession like PRJEB19386
    #  * using the taxon_id of your species
    # 'rnaseq_summary_file' should always be set. If 'taxon_id' or 'study_accession' are not undef
    # they will be used to retrieve the information from ENA and to create the csv file. In this case,
    # 'file_columns' and 'summary_file_delimiter' should not be changed unless you know what you are doing
    'study_accession'        => '',

    'max_reads_per_split' => 2500000, # This sets the number of reads to split the fastq files on
    'max_total_reads'     => 200000000, # This is the total number of reads to allow from a single, unsplit file

    'summary_file_delimiter' => '\t', # Use this option to change the delimiter for your summary data file
    'summary_csv_table' => 'csv_data',
    'read_length_table' => 'read_length',
    'rnaseq_data_provider' => 'ENA', #It will be set during the pipeline or it will use this value

    'rnaseq_dir'    => catdir($self->o('output_path'), 'rnaseq'),
    'input_dir'     => catdir($self->o('rnaseq_dir'),'input'),
    'output_dir'    => catdir($self->o('rnaseq_dir'),'output'),
    'merge_dir'     => catdir($self->o('rnaseq_dir'),'merge'),
    'sam_dir'       => catdir($self->o('rnaseq_dir'),'sams'),
    'header_file'   => catfile($self->o('output_dir'), '#'.$self->o('read_id_tag').'#_header.h'),

    'rnaseq_ftp_base' => 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/',

    'long_read_dir'       => catdir($self->o('output_path'),'long_read'),
    'long_read_fastq_dir' => catdir($self->o('long_read_dir'),'input'),
    'use_ucsc_naming' => 0,

    # If your reads are unpaired you may want to run on slices to avoid
    # making overlong rough models.  If you want to do this, specify a
    # slice length here otherwise it will default to whole chromosomes.
    slice_length => 10000000,

    # Regular expression to allow FastQ files to be correctly paired,
    # for example: file_1.fastq and file_2.fastq could be paired using
    # the expression "\S+_(\d)\.\S+".  Need to identify the read number
    # in brackets; the name the read number (1, 2) and the
    # extension.
    pairing_regex => '\S+_(\d)\.\S+',
    paired => 1,

    # Do you want to make models for the each individual sample as well
    # as for the pooled samples (1/0)?
    single_tissue => 1,

    # What Read group tag would you like to group your samples
    # by? Default = ID
    read_group_tag => 'SM',
    read_id_tag => 'ID',

    use_threads => 3,
    rnaseq_merge_threads => 12,
    rnaseq_merge_type => 'samtools',
    read_min_paired => 50,
    read_min_mapped => 50,
    other_isoforms => 'other', # If you don't want isoforms, set this to undef
    maxintron => 200000,

    # Please assign some or all columns from the summary file to the
    # some or all of the following categories.  Multiple values can be
    # separted with commas. ID, SM, DS, CN, is_paired, filename, read_length, is_13plus,
    # is_mate_1 are required. If pairing_regex can work for you, set is_mate_1 to -1.
    # You can use any other tag specified in the SAM specification:
    # http://samtools.github.io/hts-specs/SAMv1.pdf

    ####################################################################
    # This is just an example based on the file snippet shown below.  It
    # will vary depending on how your data looks.
    ####################################################################
    file_columns      => ['SM', 'ID', 'is_paired', 'filename', 'is_mate_1', 'read_length', 'is_13plus', 'CN', 'PL', 'DS'],
    long_read_columns => ['sample','filename'],

# lincRNA pipeline stuff
    'lncrna_dir' => catdir($self->o('output_path'), 'lincrna'),
    registry_file => catfile($self->o('lncrna_dir'), 'registry.pm'),
    'file_translations' => catfile($self->o('lncrna_dir'), 'hive_dump_translations.fasta'),
    'file_for_length' => catfile($self->o('lncrna_dir'), 'check_lincRNA_length.out'),  # list of genes that are smaller than 200bp, if any
    'file_for_biotypes' => catfile($self->o('lncrna_dir'), 'check_lincRNA_need_to_update_biotype_antisense.out'), # mysql queries that will apply or not in your dataset (check update_database) and will update biotypes
    'file_for_introns_support' => catfile($self->o('lncrna_dir'), 'check_lincRNA_Introns_supporting_evidence.out'), # for debug
    biotype_output => 'rnaseq',
    lincrna_protein_coding_set => [
      'rnaseq_merged_1',
      'rnaseq_merged_2',
      'rnaseq_merged_3',
      'rnaseq_merged_4',
      'rnaseq_merged_5',
      'rnaseq_tissue_1',
      'rnaseq_tissue_2',
      'rnaseq_tissue_3',
      'rnaseq_tissue_4',
      'rnaseq_tissue_5',
    ],

########################
# SPLIT PROTEOME File
########################
    'max_seqs_per_file' => 20,
    'max_seq_length_per_file' => 20000, # Maximum sequence length in a file
    'max_files_per_directory' => 1000, # Maximum number of files in a directory
    'max_dirs_per_directory'  => $self->o('max_files_per_directory'),

########################
# FINAL Checks parameters - Update biotypes to lincRNA, antisense, sense, problem ...
########################

     update_database => 'yes', # Do you want to apply the suggested biotypes? yes or no.

########################
# Interproscan
########################
    required_externalDb => '',
    interproscan_lookup_applications => [
      'PfamA',
    ],
    required_externalDb => [],
    pathway_sources => [],
    required_analysis => [
      {
        'logic_name'    => 'pfam',
        'db'            => 'Pfam',
        'db_version'    => '31.0',
        'ipscan_name'   => 'Pfam',
        'ipscan_xml'    => 'PFAM',
        'ipscan_lookup' => 1,
      },
    ],




# Max internal stops for projected transcripts
    'projection_pid'                        => '50',
    'projection_cov'                        => '50',
    'projection_max_internal_stops'         => '1',
    'projection_calculate_coverage_and_pid' => '1',

    'projection_lincrna_percent_id'         => 90,
    'projection_lincrna_coverage'           => 90,
    'projection_pseudogene_percent_id'      => 60,
    'projection_pseudogene_coverage'        => 75,
    'projection_ig_tr_percent_id'           => 70,
    'projection_ig_tr_coverage'             => 90,
    'projection_exonerate_padding'          => 5000,

    'realign_table_name'                    => 'projection_source_sequences',
    'max_projection_structural_issues'      => 1,

## Add in genewise path and put in matching code
    'genewise_pid'                        => '50',
    'genewise_cov'                        => '50',
    'genewise_region_padding'             => '50000',
    'genewise_calculate_coverage_and_pid' => '1',

########################
# Misc setup info
########################
    'repeatmasker_engine'       => 'crossmatch',
    'masking_timer_long'        => '5h',
    'masking_timer_short'       => '2h',

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
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
    'genblast_mem'         => '1900',
    'genblast_retry_mem'   => '4900',
    'genewise_mem'         => '3900',
    'genewise_retry_mem'   => '5900',
    'refseq_mem'           => '9900',
    'projection_mem'       => '1900',
    'layer_annotation_mem' => '3900',
    'genebuilder_mem'      => '1900',



########################
# LastZ
########################

    'compara_master'             => 'compara_master',
    'compara_conf_file'             => '',
    'compara_innodb_schema'         => 1,
    'compara_genome_db_update_path' => catfile($self->o('enscode_root_dir'),'/ensembl-compara/scripts/pipeline/update_genome.pl'),
    'compara_mlss_script_path'      => catfile($self->o('enscode_root_dir'),'/ensembl-compara/scripts/pipeline/create_mlss.pl'),
    'compara_mlss_reg_conf_path'    => catfile($self->o('enscode_root_dir'),'/ensembl-compara/scripts/pipeline/production_reg_ensembl_conf.pl'),
    'compara_populate_new_database_exe' => catfile($self->o('enscode_root_dir'),'ensembl-compara/scripts/pipeline/populate_new_database.pl'),
    'compara_only_cellular_component' => undef,
    'compara_dump_dir'              => catdir($self->o('output_path'),'lastz'),

    'mlss_id_list' => undef,
    'compara_collection' => '',

    'compara_ref_species'       => $self->o('projection_source_production_name'),
    'compara_non_ref_species'   => $self->o('production_name'),
    'only_cellular_component'   => undef,   # Do we load *all* the dnafrags or only the ones from a specific cellular-component ?
    'mix_cellular_components'   => 0,       # Do we try to allow the nuclear genome vs MT, etc ?
    'dump_min_nib_size'         => 11500000,
    'dump_min_chunk_size'       => 1000000,
    'dump_min_chunkset_size'    => 1000000,
    'quick' => 1,
    'default_chunks' => {
      'reference'   => {
        'homo_sapiens' => {
          'chunk_size' => 30000000,
          'overlap'    => 0,
          'include_non_reference' => -1, #1  => include non_reference regions (eg human assembly patches)
                                         #0  => do not include non_reference regions
                                         #-1 => auto-detect (only include non_reference regions if the non-reference species is high-coverage
                                         #ie has chromosomes since these analyses are the only ones we keep up-to-date with the patches-pipeline)
          'masking_options' => '{default_soft_masking => 1}',
           # if you have a specific selection of repeat elements for the masking
           #'masking_options_file' => $self->check_file_in_ensembl('ensembl-compara/scripts/pipeline/human36.spec'),
        },
        #non human example
        'default' => {
          'chunk_size'      => 10000000,
          'overlap'         => 0,
          'masking_options' => '{default_soft_masking => 1}'
        },
      },
      'non_reference' => {
        'chunk_size'      => 10100000,
        'group_set_size'  => 10100000,
        'overlap'         => 100000,
        'masking_options' => '{default_soft_masking => 1}'
      },
    },

    'compara_window_size' => 10000,
    'filter_duplicates_rc_name' => '2GB_lastz',
    'filter_duplicates_himem_rc_name' => '8GB_lastz',

   #
    #Default pair_aligner
    #
    'pair_aligner_method_link' => [1001, 'LASTZ_RAW'],
    'pair_aligner_logic_name' => 'LastZ',
    'pair_aligner_module' => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::LastZ',

    'pair_aligner_options' => {
       default => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac', # ensembl genomes settings
       7742    => 'T=1 K=3000 L=3000 H=2200 O=400 E=30 --ambiguous=iupac', # vertebrates - i.e. ensembl-specific
       9526    => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 Q=' . $self->check_file_in_ensembl('ensembl-compara/scripts/pipeline/primate.matrix').' --ambiguous=iupac', # primates
       33554   => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 --ambiguous=iupac', # carnivora
       3913    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
       4070    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
    },

    #
    #Default chain
    #
    'chain_input_method_link' => [1001, 'LASTZ_RAW'],
    'chain_output_method_link' => [1002, 'LASTZ_CHAIN'],

    #linear_gap=>medium for more closely related species, 'loose' for more distant
    'linear_gap' => 'medium',

    'chain_parameters' => {'max_gap'=>'50','linear_gap'=> $self->o('linear_gap'), 'faToNib' => $self->o('faToNib_exe'), 'lavToAxt'=> $self->o('lavToAxt_exe'), 'axtChain'=>$self->o('axtChain_exe'), 'max_blocks_for_chaining' => 100000},

    #
    #Default patch_alignments
    #
    'patch_alignments' => 0,  #set to 1 to align the patches of a species to many other species

    #
    #Default net
    #
    'net_input_method_link' => [1002, 'LASTZ_CHAIN'],
    'net_output_method_link' => [16, 'LASTZ_NET'],
    'net_ref_species' => $self->o('compara_ref_species'),  #default to ref_species
    'net_parameters' => {'max_gap'=>'50', 'chainNet'=>$self->o('chainNet_exe')},
    'bidirectional' => 0,

    #
    #Default healthcheck
    #
    'previous_db' => 'compara_prev',
    'prev_release' => 0,   # 0 is the default and it means "take current release number and subtract 1"
    'max_percent_diff' => 20,
    'max_percent_diff_patches' => 99.99,
    'do_pairwise_gabs' => 1,
    'do_compare_to_previous_db' => 0,

    'compara_bed_dir' => $self->o('compara_dump_dir').'/bed_dir',
    'compara_feature_dir' => $self->o('compara_dump_dir').'/feature_dumps',

    #
    #Default pairaligner config
    #
    'skip_pairaligner_stats' => 1, #skip this module if set to 1

    'pair_aligner_method_link' => [1001, 'LASTZ_RAW'],
    'pair_aligner_logic_name' => 'LastZ',
    'pair_aligner_module' => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::LastZ',
    'chain_input_method_link' => [1001, 'LASTZ_RAW'],
    'chain_output_method_link' => [1002, 'LASTZ_CHAIN'],
    'linear_gap' => 'medium',
    'net_input_method_link' => [1002, 'LASTZ_CHAIN'],
    'net_output_method_link' => [16, 'LASTZ_NET'],

    # Capacities
    'pair_aligner_analysis_capacity' => 700,
    'pair_aligner_batch_size' => 40,
    'chain_hive_capacity' => 200,
    'chain_batch_size' => 10,
    'net_hive_capacity' => 300,
    'net_batch_size' => 10,
    'filter_duplicates_hive_capacity' => 200,
    'filter_duplicates_batch_size' => 10,

    # LastZ is used to align the genomes
    'pair_aligner_exe'  => $self->o('lastz_exe'),
    'cellar_dir'                        => '/nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/Cellar/',
    'lastz_exe'                         => catfile($self->o('cellar_dir'),'lastz/1.04.00/bin/lastz'),
    'axtChain_exe'                      => catfile($self->o('cellar_dir'),'kent/v335_1/bin/axtChain'),
    'chainNet_exe'                      => catfile($self->o('cellar_dir'),'kent/v335_1/bin/chainNet'),
    'faToNib_exe'                       => catfile($self->o('cellar_dir'),'kent/v335_1/bin/faToNib'),
    'lavToAxt_exe'                      => catfile($self->o('cellar_dir'),'kent/v335_1/bin/lavToAxt'),
    'compare_beds_exe'                  => catfile($self->o('enscode_root_dir'),'ensembl-compara/scripts/pipeline/compare_beds.pl'),
    'create_pair_aligner_page_exe'      => catfile($self->o('enscode_root_dir'),'ensembl-compara/scripts/report/create_pair_aligner_page.pl'),
    'dump_features_exe'                 => catfile($self->o('enscode_root_dir'),'ensembl-compara/scripts/dumps/DumpMultiAlign.pl'),


########################
# db info
########################
    'reference_db' => {
      -dbname => $self->o('reference_db_name'),
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


   'compara_db' => {
      -dbname => $self->o('compara_db_name'),
      -host   => $self->o('compara_db_server'),
      -port   => $self->o('compara_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'cdna_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_cdna_'.$self->o('release_number'),
      -host   => $self->o('cdna_db_server'),
      -port   => $self->o('cdna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'genblast_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_genblast_'.$self->o('release_number'),
      -host   => $self->o('genblast_db_server'),
      -port   => $self->o('genblast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'genblast_nr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_genblast_nr_'.$self->o('release_number'),
      -host   => $self->o('genblast_db_server'),
      -port   => $self->o('genblast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'genblast_rnaseq_support_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_genblast_rnaseq_'.$self->o('release_number'),
      -host   => $self->o('genblast_rnaseq_support_db_server'),
      -port   => $self->o('genblast_rnaseq_support_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'genblast_rnaseq_support_nr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_genblast_rnaseq_nr_'.$self->o('release_number'),
      -host   => $self->o('genblast_rnaseq_support_db_server'),
      -port   => $self->o('genblast_rnaseq_support_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'ig_tr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_igtr_'.$self->o('release_number'),
      -host   => $self->o('ig_tr_db_server'),
      -port   => $self->o('ig_tr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    cdna2genome_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_cdna2genome_'.$self->o('release_number'),
      -host   => $self->o('cdna2genome_db_server'),
      -port   => $self->o('cdna2genome_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genewise_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_genewise_'.$self->o('release_number'),
      -host   => $self->o('genewise_db_server'),
      -port   => $self->o('genewise_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'best_targeted_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_bt_'.$self->o('release_number'),
      -host   => $self->o('genewise_db_server'),
      -port   => $self->o('genewise_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    long_read_initial_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_lrinitial_'.$self->o('release_number'),
      -host   => $self->o('long_read_initial_db_server'),
      -port   => $self->o('long_read_initial_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    long_read_collapse_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_lrcollapse_'.$self->o('release_number'),
      -host => $self->o('long_read_collapse_db_server'),
      -port => $self->o('long_read_collapse_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    long_read_blast_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_lrblast_'.$self->o('release_number'),
      -host => $self->o('long_read_blast_db_server'),
      -port => $self->o('long_read_blast_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    long_read_final_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_lrfinal_'.$self->o('release_number'),
      -host => $self->o('long_read_final_db_server'),
      -port => $self->o('long_read_final_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_proj_'.$self->o('release_number'),
      -host   => $self->o('projection_db_server'),
      -port   => $self->o('projection_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'selected_projection_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_sel_proj_'.$self->o('release_number'),
      -host   => $self->o('projection_db_server'),
      -port   => $self->o('projection_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_realign_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_realign_'.$self->o('release_number'),
      -host   => $self->o('projection_realign_db_server'),
      -port   => $self->o('projection_realign_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_lincrna_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_proj_linc_'.$self->o('release_number'),
      -host   => $self->o('projection_lincrna_db_server'),
      -port   => $self->o('projection_lincrna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_pseudogene_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_proj_pseudo_'.$self->o('release_number'),
      -host   => $self->o('projection_pseudogene_db_server'),
      -port   => $self->o('projection_pseudogene_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_source_db' => {
      -dbname => $self->o('projection_source_db_name'),
      -host   => $self->o('projection_source_db_server'),
      -port   => $self->o('projection_source_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'projection_lastz_db' => {
      -dbname => $self->o('projection_lastz_db_name'),
      -host   => $self->o('projection_lastz_db_server'),
      -port   => $self->o('projection_lastz_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_for_layer_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_rnaseq_layer_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_db_server'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'rnaseq_for_layer_nr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_rnaseq_layer_nr_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_db_server'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_rnaseq_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_db_server'),
      -port   => $self->o('rnaseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_blast_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_rnaseq_blast_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_blast_db_server'),
      -port   => $self->o('rnaseq_blast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_refine_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_refine_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_refine_db_server'),
      -port   => $self->o('rnaseq_refine_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_rough_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_rough_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_rough_db_server'),
      -port   => $self->o('rnaseq_rough_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    lincrna_db => {
      -host   => $self->o('lincrna_db_server'),
      -port   => $self->o('lincrna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_lincrna_'.$self->o('release_number'),
      -driver => $self->o('hive_driver'),
    },

    'layering_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_layer_'.$self->o('release_number'),
      -host   => $self->o('layering_db_server'),
      -port   => $self->o('layering_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'utr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_utr_'.$self->o('release_number'),
      -host   => $self->o('utr_db_server'),
      -port   => $self->o('utr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genebuilder_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_gbuild_'.$self->o('release_number'),
      -host   => $self->o('genebuilder_db_server'),
      -port   => $self->o('genebuilder_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'pseudogene_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_pseudo_'.$self->o('release_number'),
      -host   => $self->o('pseudogene_db_server'),
      -port   => $self->o('pseudogene_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'ncrna_db' => {
      -dbname => $self->o('ncrna_db_name'),
      -host   => $self->o('ncrna_db_server'),
      -port   => $self->o('ncrna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'final_geneset_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_final_'.$self->o('release_number'),
      -host   => $self->o('final_geneset_db_server'),
      -port   => $self->o('final_geneset_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'refseq_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_refseq_'.$self->o('release_number'),
      -host   => $self->o('refseq_db_server'),
      -port   => $self->o('refseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host   => $self->o('killlist_db_server'),
      -port   => $self->o('killlist_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
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

    'otherfeatures_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_otherfeatures_'.$self->o('release_number'),
      -host   => $self->o('otherfeatures_db_server'),
      -port   => $self->o('otherfeatures_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
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

################
# LastZ
################

    my $second_pass     = exists $self->{'_is_second_pass'};
    $self->{'_is_second_pass'} = $second_pass;
    return $self->SUPER::pipeline_create_commands if $self->can('no_compara_schema');
    my $pipeline_url    = $self->pipeline_url();
    my $parsed_url      = $second_pass && Bio::EnsEMBL::Hive::Utils::URL::parse( $pipeline_url );
    my $driver          = $second_pass ? $parsed_url->{'driver'} : '';

################
# /LastZ
################

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

      'mkdir -p '.$self->o('rnaseq_dir'),
      'mkdir -p '.$self->o('long_read_fastq_dir'),
      'mkdir -p '.$self->o('genome_dumps'),

# Commenting out lincRNA pfam pipeline commands until we put that bit back in
#"cat <<EOF > ".$self->o('registry_file')."
#{
#package reg;

#Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#-host => '".$self->o('lincrna_db', '-host')."',
#-port => ".$self->o('lincrna_db', '-port').",
#-user => '".$self->o('lincrna_db', '-user')."',
#-pass => '".$self->o('lincrna_db', '-pass')."',
#-dbname => '".$self->o('lincrna_db', '-dbname')."',
#-species => '".$self->o('species_name')."',
#-WAIT_TIMEOUT => undef,
#-NO_CACHE => undef,
#-VERBOSE => '1',
#);

#Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#-host => '".$self->o('production_db', '-host')."',
#-port => ".$self->o('production_db', '-port').",
#-user => '".$self->o('production_db', '-user')."',
#-dbname => '".$self->o('production_db', '-dbname')."',
#-species => 'multi',
#-group => 'production'
#);

#1;
#}
#EOF",

#################
# LastZ
#################

     'mkdir -p '.$self->o('compara_dump_dir'),
     'mkdir -p '.$self->o('compara_bed_dir'),
      # Compara 'release' tables will be turned from MyISAM into InnoDB on the fly by default:
      ($self->o('compara_innodb_schema') ? "sed 's/ENGINE=MyISAM/ENGINE=InnoDB/g' " : 'cat ')
      . $self->check_file_in_ensembl('ensembl-compara/sql/table.sql').' | '.$self->db_cmd(),

      # Compara 'pipeline' tables are already InnoDB, but can be turned to MyISAM if needed:
      ($self->o('compara_innodb_schema') ? 'cat ' : "sed 's/ENGINE=InnoDB/ENGINE=MyISAM/g' ")
      . $self->check_file_in_ensembl('ensembl-compara/sql/pipeline-tables.sql').' | '.$self->db_cmd(),

      # MySQL specific procedures
      $driver eq 'mysql' ? ($self->db_cmd().' < '.$self->check_file_in_ensembl('ensembl-compara/sql/procedures.'.$driver)) : (),

#################
# /LastZ
#################


    ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    skip_post_repeat_analyses => $self->o('skip_post_repeat_analyses'),
    skip_projection => $self->o('skip_projection'),
    skip_rnaseq => $self->o('skip_rnaseq'),
    skip_ncrna => $self->o('skip_ncrna'),
    skip_long_read => $self->o('skip_long_read'),
    skip_lastz => $self->o('skip_lastz'),
    skip_repeatmodeler => $self->o('skip_repeatmodeler'),
    load_toplevel_only => $self->o('load_toplevel_only'),
    wide_repeat_logic_names => $self->o('use_repeatmodeler_to_mask') ? [$self->o('full_repbase_logic_name'),$self->o('repeatmodeler_logic_name'),'dust'] :
                                                                                       [$self->o('full_repbase_logic_name'),'dust'],
    wide_ensembl_release => $self->o('ensembl_release'),
    use_genome_flatfile  => $self->o('use_genome_flatfile'),
    genome_file          => $self->o('faidx_genome_file'),
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

    my %genblast_params = (
      wu    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
      ncbi  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
      wu_genome    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
      ncbi_genome  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
      wu_projection    => '-P wublast -gff -e #blast_eval# -c #blast_cov# -n 100 -x 5 ',
      ncbi_projection  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
      );
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

     {
        -logic_name => 'download_rnaseq_csv',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
        -rc_name => '1GB',
        -parameters => {
          study_accession => $self->o('study_accession'),
          taxon_id => $self->o('species_taxon_id'),
          inputfile => $self->o('rnaseq_summary_file'),
        },

        -flow_into => {
           1 => ['download_genus_rnaseq_csv'],
         },

        -input_ids  => [
			{
            assembly_name => $self->o('assembly_name'),
            assembly_accession => $self->o('assembly_accession'),
            assembly_refseq_accession => $self->o('assembly_refseq_accession'),
          },
        ]
      },


     {
        -logic_name => 'download_genus_rnaseq_csv',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
        -rc_name => '1GB',
        -parameters => {
          study_accession => $self->o('study_accession'),
          taxon_id => $self->o('genus_taxon_id'),
          inputfile => $self->o('rnaseq_summary_file_genus'),
        },

        -flow_into => {
           1 => ['download_long_read_csv'],
         },
      },


     {
        -logic_name => 'download_long_read_csv',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
        -rc_name => '1GB',
        -parameters => {
          study_accession => $self->o('study_accession'),
          taxon_id => $self->o('species_taxon_id'),
          inputfile => $self->o('long_read_summary_file'),
          read_type => 'isoseq',
        },

        -flow_into => {
           1 => ['download_genus_long_read_csv'],
         },
      },


     {
        -logic_name => 'download_genus_long_read_csv',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
        -rc_name => '1GB',
        -parameters => {
          study_accession => $self->o('study_accession'),
          taxon_id => $self->o('genus_taxon_id'),
          inputfile => $self->o('long_read_summary_file_genus'),
          read_type => 'isoseq',
        },

        -flow_into => {
           1 => ['create_core_db'],
         },
      },


      {
        # Creates a reference db for each species
        -logic_name => 'create_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
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
      {
        # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
        -logic_name => 'process_assembly_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessAssemblyReport',
        -parameters => {
          full_ftp_path => $self->o('assembly_ftp_path'),
          output_path   => $self->o('output_path'),
          target_db     => $self->o('reference_db'),
        },
        -rc_name    => '8GB',
        -flow_into  => {
          1 => ['load_meta_info'],
        },
      },

      {
        # Load some meta info and seq_region_synonyms
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
              '(1, "provider.name", "'.$self->o('provider_name').'"),'.
              '(1, "provider.url", "'.$self->o('provider_url').'"),'.
              '(1, "species.production_name", "'.$self->o('production_name').$self->o('production_name_modifier').'"),'.
              '(1, "repeat.analysis", "'.$self->o('full_repbase_logic_name').'"),'.
              ($self->o('use_repeatmodeler_to_mask') ? '(1, "repeat.analysis", "'.$self->o('repeatmodeler_logic_name').'"),': '').
              '(1, "repeat.analysis", "dust"),'.
              '(1, "repeat.analysis", "trf")',
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
      {
#       Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
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

      {
# Get the prefixes for all contigs from the AGP files
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

      {
# Download contig from NCBI
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


      {
# Load the contigs into each reference db
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


      {
# Load the AGP files
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


      {
# Set the toplevel
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


      {
# Load some meta info and seq_region_synonyms
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
            'provider.name' => $self->o('provider_name'),
            'provider.url' => $self->o('provider_url'),
            'repeat.analysis' => [$self->o('full_repbase_logic_name'), 'dust', 'trf'],
            'species.production_name' => $self->o('production_name').$self->o('production_name_modifier'),
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
                          1 => ['load_mitochondrion', 'fan_refseq_import'],
                       },
      },

###############################################################################
#
# REFSEQ ANNOTATION
#
###############################################################################
      {
        -logic_name => 'fan_refseq_import',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ -n "#assembly_refseq_accession#" ]; then exit 0; else exit 42;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['load_refseq_synonyms'],
                       },
      },

      {
        -logic_name => 'load_refseq_synonyms',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadRefSeqSynonyms',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
			 'output_dir'       => $self->o('output_path'),
			 'url'              => $self->o('refseq_report_ftp_path'),
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['download_refseq_gff'],
                       },
      },


      {
        -logic_name => 'download_refseq_gff',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData',
        -parameters => {
          output_dir => catdir($self->o('output_path'), 'refseq_import'),
          url => $self->o('refseq_import_ftp_path'),
          download_method => 'ftp',
        },
        -rc_name => 'default',
        -flow_into  => {
                         1 => ['create_refseq_db'],
                       },
      },


      {
        -logic_name => 'create_refseq_db',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('refseq_db'),
                         create_type => 'clone',
                       },
       -rc_name => 'default',
       -flow_into  => {
                          1 => ['load_refseq_gff'],
                      },
      },

      {
        -logic_name => 'load_refseq_gff',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('refseq_import_script_path').
                                      ' -dnahost '.$self->o('dna_db','-host').
                                      ' -dnadbname '.$self->o('dna_db','-dbname').
                                      ' -dnaport '.$self->o('dna_db','-port').
                                      ' -dnauser '.$self->o('dna_db','-user').
                                      ' -user '.$self->o('user').
                                      ' -pass '.$self->o('password').
                                      ' -host '.$self->o('refseq_db','-host').
                                      ' -port '.$self->o('refseq_db','-port').
                                      ' -dbname '.$self->o('refseq_db','-dbname').
                                      ' -write'.
                                      ' -file '.catfile($self->o('output_path'), 'refseq_import', '#assembly_refseq_accession#_#assembly_name#_genomic.gff'),
                       },
        -rc_name => 'refseq_import',
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

        -flow_into  => {
          '1->A' => ['create_10mb_slice_ids'],
          'A->1' => ['genome_prep_sanity_checks'],
        },

      },


###############################################################################
#
# REPEATMASKER ANALYSES
#
###############################################################################
      {
        # Create 10mb toplevel slices, these will be split further for repeatmasker
        -logic_name => 'create_10mb_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 10000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                         batch_slice_ids => 1,
                         batch_target_size => 10000000,
                       },
        -rc_name    => '2GB',
        -flow_into => {
                         '2'    => ['semaphore_10mb_slices'],
                      },
      },


      {
         # Wait for repeatmasker to complete all the sub slices for a 10mb slice before passing to dust
        -logic_name => 'semaphore_10mb_slices',
         -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
         -parameters => {},
         -flow_into => {
                        '1->A' => ['create_repeatmasker_slices'],
                        'A->1' => ['run_dust'],
                       },
      },


      {
        -logic_name => 'create_repeatmasker_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 1000000,
                         batch_target_size => 500000,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '2' => ['run_repeatmasker','fan_repeatmodeler'],
                      },
      },


      {
        -logic_name => 'run_repeatmasker',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => $self->o('masking_timer_long'),
                         target_db => $self->o('reference_db'),
                         logic_name => $self->o('full_repbase_logic_name'),
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -species "'.$self->o('repbase_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                         use_genome_flatfile => 1,
                         genome_file => $self->o('faidx_genome_file'),
                       },
        -rc_name    => 'repeatmasker',
        -flow_into => {
                        '-1' => ['rebatch_repeatmasker'],
                        '-2' => ['rebatch_repeatmasker'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'rebatch_repeatmasker',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 10000,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '2' => ['run_repeatmasker_small_batch'],
                      },
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'run_repeatmasker_small_batch',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => $self->o('masking_timer_short'),
                         target_db => $self->o('reference_db'),
                         logic_name => $self->o('full_repbase_logic_name'),
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -species "'.$self->o('repbase_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                         use_genome_flatfile => 1,
                         genome_file => $self->o('faidx_genome_file'),
                       },
        -rc_name    => 'repeatmasker_rebatch',
        -flow_into => {
                         -1 => ['failed_repeatmasker_batches'],
                         -2 => ['failed_repeatmasker_batches'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'failed_repeatmasker_batches',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'fan_repeatmodeler',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ -n "'.$self->o('repeatmodeler_library').'" ] && [ "#skip_repeatmodeler#" == "0" ]; then exit 0; else exit 42;fi',
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => { '1' => ['run_repeatmasker_repeatmodeler'] },
      },


      {
        -logic_name => 'run_repeatmasker_repeatmodeler',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => $self->o('masking_timer_long'),
                         target_db => $self->o('reference_db'),
                         logic_name => 'repeatmask_repeatmodeler',
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -lib "'.$self->o('repeatmodeler_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                         use_genome_flatfile => 1,
                         genome_file => $self->o('faidx_genome_file'),
                       },
        -rc_name    => 'repeatmasker',
        -flow_into => {
                        '-1' => ['rebatch_repeatmasker_repeatmodeler'],
                        '-2' => ['rebatch_repeatmasker_repeatmodeler'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'rebatch_repeatmasker_repeatmodeler',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 10000,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '2' => ['run_repeatmasker_repeatmodeler_small_batch'],
                      },
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'run_repeatmasker_repeatmodeler_small_batch',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => $self->o('masking_timer_short'),
                         target_db => $self->o('reference_db'),
                         logic_name => $self->o('repeatmodeler_logic_name'),
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -lib "'.$self->o('repeatmodeler_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                         use_genome_flatfile => 1,
                         genome_file => $self->o('faidx_genome_file'),
                       },
        -rc_name    => 'repeatmasker_rebatch',
        -flow_into => {
                         -1 => ['failed_repeatmasker_repeatmodeler_batches'],
                         -2 => ['failed_repeatmasker_repeatmodeler_batches'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'failed_repeatmasker_repeatmodeler_batches',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },



      {
        # Set the toplevel
        -logic_name => 'dump_softmasked_toplevel',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
        -parameters => {
                         'coord_system_name'    => 'toplevel',
                         'target_db'            => $self->o('reference_db'),
                         'output_path'          => $self->o('genome_dumps'),
                         'enscode_root_dir'     => $self->o('enscode_root_dir'),
                         'species_name'         => $self->o('species_name'),
                         'repeat_logic_names'   => '#wide_repeat_logic_names#',
                       },
        -input_ids => [{}],
        -wait_for => ['run_dust'],
        -flow_into => {
          1 => ['format_softmasked_toplevel'],
        },
        -rc_name    => '3GB',
      },


      {
        # This should probably be a proper module
        -logic_name => 'format_softmasked_toplevel',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'    => 'if [ "'.$self->o('blast_type').'" = "ncbi" ]; then convert2blastmask -in '.$self->o('softmasked_genome_file').' -parse_seqids -masking_algorithm repeatmasker -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out '.$self->o('softmasked_genome_file').'.asnb;makeblastdb -in '.$self->o('softmasked_genome_file').' -dbtype nucl -parse_seqids -mask_data '.$self->o('softmasked_genome_file').'.asnb -title "'.$self->o('species_name').'"; else xdformat -n '.$self->o('softmasked_genome_file').';fi',
                       },
        -rc_name    => '5GB',
        -flow_into => {
          1 => ['create_reheadered_softmasked_file'],
        },
      },


      {
        -logic_name => 'create_reheadered_softmasked_file',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'convert_genome_dump.pl').
                                      ' -conversion_type slice_name_to_seq_region_name'.
                                      ' -input_file '.$self->o('softmasked_genome_file').
                                      ' -output_file '.$self->o('faidx_softmasked_genome_file'),
                       },
        -rc_name => 'default',
        -flow_into => {
          1 => ['create_softmasked_faidx'],
        },
      },


      {
        -logic_name => 'create_softmasked_faidx',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	     -rc_name => '5GB',
	     -parameters => {
          cmd => 'if [ ! -e "'.$self->o('faidx_softmasked_genome_file').'.fai" ]; then '.$self->o('samtools_path').' faidx '.$self->o('faidx_softmasked_genome_file').';fi',
			     },

      },

###############################################################################
#
# SIMPLE FEATURE AND OTHER REPEAT ANALYSES
#
###############################################################################

      {
        # Run dust
        -logic_name => 'run_dust',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDust',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'dust',
                         module => 'HiveDust',
                         dust_path => $self->o('dust_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['run_trf'],
                         -1 => ['run_trf'],
		         -2 => ['run_trf'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        # Run TRF
        -logic_name => 'run_trf',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRF',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'trf',
                         module => 'HiveTRF',
                         trf_path => $self->o('trf_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['fan_post_repeat_analyses'],
                        -1 => ['fan_post_repeat_analyses'],
                        -2 => ['fan_post_repeat_analyses'],
                      },
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -batch_size => 20,
      },


      {
        # This will skip downstream analyses like cpg, eponine, genscan etc. if the flag is set
        -logic_name => 'fan_post_repeat_analyses',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ "#skip_post_repeat_analyses#" -ne "0" ]; then exit 42; else exit 0;fi',
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => { '1' => ['run_eponine'] },
      },



      {
        # Run eponine
        -logic_name => 'run_eponine',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveEponine',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'eponine',
                         module => 'HiveEponine',
                         eponine_path => $self->o('eponine_java_path'),
                         commandline_params => '-epojar => '.$self->o('eponine_jar_path').', -threshold => 0.999',
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['run_cpg'],
                         -1 => ['run_cpg'],
                         -2 => ['run_cpg'],
                      },
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -batch_size => 20,
      },


      {
        # Run CPG
        -logic_name => 'run_cpg',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveCPG',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'cpg',
                         module => 'HiveCPG',
                         cpg_path => $self->o('cpg_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                        1 => ['run_trnascan'],
                        -1 => ['run_trnascan'],
                        -2 => ['run_trnascan'],
                      },
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -batch_size => 20,
      },


      {
        # Run tRNAscan
        -logic_name => 'run_trnascan',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRNAScan',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'trnascan',
                         module => 'HiveTRNAScan',
                         trnascan_path => $self->o('trnascan_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['create_genscan_slices'],
                         -1 => ['create_genscan_slices'],
                         -2 => ['create_genscan_slices'],

                      },
         -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
         -batch_size => 20,
      },

###############################################################################
#
# GENSCAN ANALYSIS
#
##############################################################################
      {
         # Genscan has issues with large slices, so instead rebatch the 10mb slices
         # into 1mb slices with a target batch size of 10mb
        -logic_name => 'create_genscan_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 1000000,
                         batch_target_size => 10000000,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '2' => ['run_genscan'],
                      },
      },


      {
        # Run genscan, uses 1mb slices from repeatmasker. Flows into create_prediction_transcript_ids which
        # then takes these 1mb slices and converts them into individual prediction transcript input ids based
        # on the dbID of each feature generate by this analysis
        -logic_name => 'run_genscan',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'genscan',
                         module => 'HiveGenscan',
                         genscan_path => $self->o('genscan_path'),
                         genscan_matrix_path => $self->o('genscan_matrix_path'),
                         repeat_masking_logic_names => [$self->o('full_repbase_logic_name')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        # No need to semaphore the jobs with issues as the blast analyses work off prediction transcript
                        # ids from slices that genscan succeeds on. So passing small slices in and ignore failed slices is fine
                        1 => ['fan_genscan_blasts'],
                        -1 => ['decrease_genscan_slice_size'],
                        -2 => ['decrease_genscan_slice_size'],
                        -3 => ['decrease_genscan_slice_size'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        -logic_name => 'decrease_genscan_slice_size',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 100000,
                       },
        -flow_into => {
                        '2' => ['run_genscan_short_slice'],
                      },
        -rc_name    => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'run_genscan_short_slice',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'genscan',
                         module => 'HiveGenscan',
                         genscan_path => $self->o('genscan_path'),
                         genscan_matrix_path => $self->o('genscan_matrix_path'),
                         repeat_masking_logic_names => [$self->o('full_repbase_logic_name')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        1 => ['fan_genscan_blasts'],
                        -1 => ['failed_genscan_slices'],
                        -2 => ['failed_genscan_slices'],
                        -3 => ['failed_genscan_slices'],
                      },
        -rc_name    => 'genscan_short',
        -can_be_empty  => 1,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'failed_genscan_slices',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        # This will skip downstream analyses like cpg, eponine, genscan etc. if the flag is set
        -logic_name => 'fan_genscan_blasts',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ "'.$self->o('skip_genscan_blasts').'" -ne "0" ]; then exit 42; else exit 0;fi',
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => { '1' => ['create_prediction_transcript_ids'] },
      },


      {
        # Create input ids for prediction transcripts. Takes a slice as an input id and converts it
        # to a set of input ids that are individual dbIDs for the prediction transcripts. This avoids empty slices
        # being submitted as jobs and also means one feature corresponds to one job. Each species flows into this
        # independantly with 1mb slices. Also has the advantage that downstream analyses can start working as soon
        # as a single slice is finished
        -logic_name => 'create_prediction_transcript_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         feature_type => 'prediction_transcript',
                         iid_type => 'slice_to_feature_ids',
                         prediction_transcript_logic_names => ['genscan'],
                         batch_size => 50,
                       },
        -flow_into => {
                        2 => ['run_uniprot_blast'],
                      },
        -rc_name    => 'default',
      },

##############################################################################
#
# BLAST ANALYSES
#
##############################################################################

      {
        # BLAST individual prediction transcripts against uniprot. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_uniprot_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan',
        -parameters => {
                         timer => '30m',
                         sequence_type => 'peptide',
                         prediction_transcript_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         repeat_masking_logic_names => [$self->o('full_repbase_logic_name')], # not sure if this is used
                         blast_db_path => $self->o('protein_blast_db'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscan',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic',
                                                 'BlastGenscanPep',
                                                 {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -flow_into => {
                        1 => ['run_vertrna_blast'],
                        -1 => ['split_blast_jobs'],
                        -2 => ['split_blast_jobs'],
                      },
        -rc_name    => 'blast',
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        # Only do the split on the uniprot ones as they're the only ones that seem to have any issues
        -logic_name => 'split_blast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk',
                         batch_size => 1,
                       },
        -rc_name      => 'default',
        -can_be_empty  => 1,
        -flow_into => {
                        2 => ['run_uniprot_blast_retry'],
                      },
      },


      {
        # BLAST individual prediction transcripts against uniprot. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_uniprot_blast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan',
        -parameters => {
                         timer => '45m',
                         sequence_type => 'peptide',
                         prediction_transcript_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         repeat_masking_logic_names => [$self->o('full_repbase_logic_name')], # not sure if this is used
                         blast_db_path => $self->o('protein_blast_db'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscan',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic',
                                                 'BlastGenscanPep',
                                                 {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -can_be_empty  => 1,
        -flow_into => {
                        1 => ['run_vertrna_blast'],
                        -1 => ['run_vertrna_blast'],
                        -2 => ['run_vertrna_blast'],
                      },
        -rc_name    => 'blast_retry',
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        # BLAST individual prediction transcripts against vertRNA. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_vertrna_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan',
        -parameters => {
                         timer => '10m',
                         sequence_type => 'dna',
                         prediction_transcript_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         repeat_masking_logic_names => [$self->o('full_repbase_logic_name')], # not sure if this is used
                         blast_db_path => $self->o('vertrna_blast_db_path'),
                         blast_exe_path => $self->o('vertrna_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         iid_type => 'feature_id',
                         logic_name => 'vertrna',
                         module => 'HiveBlastGenscan',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic',
                                                 'BlastGenscanVertRNA',
                                                 {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -flow_into => {
                        1 => ['run_unigene_blast'],
                        -1 => ['run_unigene_blast'],
                        -2 => ['run_unigene_blast'],
                      },
        -rc_name    => 'blast',
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -batch_size => 20,
      },


      {
        # BLAST individual prediction transcripts against unigene. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_unigene_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan',
        -parameters => {
                         timer => '10m',
                         sequence_type => 'dna',
                         prediction_transcript_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         repeat_masking_logic_names => [$self->o('full_repbase_logic_name')], # not sure if this is used
                         blast_db_path => $self->o('unigene_blast_db_path'),
                         blast_exe_path => $self->o('unigene_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         iid_type => 'feature_id',
                         logic_name => 'unigene',
                         module => 'HiveBlastGenscan',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic',
                                                 'BlastGenscanUnigene',
                                                 {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                       },
        -flow_into => {
                        -1 => ['failed_blast_job'],
                        -2 => ['failed_blast_job'],
                      },
        -rc_name    => 'blast',
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        -logic_name => 'failed_blast_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

      {
        -logic_name => 'genome_prep_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         sanity_check_type => 'genome_preparation_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'genome_preparation_checks')->{$self->o('uniprot_set')},
                       },

        -flow_into =>  {
                         1 => ['set_repeat_types'],
                       },
        -rc_name    => '15GB',
     },


     {
        -logic_name => 'set_repeat_types',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('repeat_types_script').
                                ' -user '.$self->o('reference_db', '-user').
                                ' -pass '.$self->o('reference_db', '-pass').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbpattern '.$self->o('reference_db','-dbname')
                       },
         -rc_name => 'default',
         -flow_into => { 1 => ['backup_core_db'] },
     },


     {
        -logic_name => 'backup_core_db',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
        -parameters => {
                         src_db_conn => $self->o('dna_db'),
                         output_file => catfile($self->o('output_path'), 'core_bak.sql.gz'),
                         dump_options => $self->o('mysql_dump_options'),
                       },
        -rc_name    => '2GB',
        -flow_into => { 1 => ['assembly_loading_report','create_selected_projection_db'] },
      },


      {
        -logic_name => 'assembly_loading_report',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('loading_report_script').
                                ' -user '.$self->o('user_r').
                                ' -host '.$self->o('dna_db','-host').
                                ' -port '.$self->o('dna_db','-port').
                                ' -dbname '.$self->o('dna_db','-dbname').
                                ' -report_type assembly_loading'.
                                ' > '.catfile($self->o('output_path'), 'loading_report.txt'),
                       },
         -rc_name => 'default',
         -flow_into => { 1 => ['email_loading_report'] },
      },

      {
        -logic_name => 'email_loading_report',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
        -parameters => {
                         email => $self->o('email_address'),
                         subject => 'AUTOMATED REPORT: assembly loading and feature annotation for '.$self->o('dna_db','-dbname').' completed',
                         text => 'Assembly loading and feature annotation have completed for '.$self->o('dna_db','-dbname').". Basic stats can be found below",
                         file => catfile($self->o('output_path'), 'loading_report.txt'),
                       },
        -rc_name => 'default',
        -flow_into => {
          1 =>['repeat_masker_coverage'],
        },
      },


      {
        -logic_name => 'repeat_masker_coverage',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatCoverage',
        -parameters => {
          source_db => $self->o('reference_db'),
          repeat_logic_names => '#wide_repeat_logic_names#',
          coord_system_version => '#assembly_name#',
          email => $self->o('email_address'),
        },
        -rc_name => '10GB',
        -flow_into => {
          1 =>['repeatmodeler_coverage'],
        },
      },


      {
        -logic_name => 'repeatmodeler_coverage',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatCoverage',
        -parameters => {
          source_db => $self->o('reference_db'),
          repeat_logic_names => [$self->o('repeatmodeler_logic_name')],
          coord_system_version => '#assembly_name#',
          email => $self->o('email_address'),
        },
        -rc_name => '4GB',
      },



########################################################################
#
# Projection analyses
#
########################################################################

      {
        -logic_name => 'create_selected_projection_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db   => $self->o('dna_db'),
          target_db   => $self->o('selected_projection_db'),
          create_type => 'clone',
        },
        -rc_name    => 'default',
        -flow_into  => {
          '1->A' => ['fan_projection'],
          'A->1' => ['create_genblast_output_db'],
        },
      },



     {
        -logic_name => 'fan_projection',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ "#skip_projection#" -ne "0" ]; then exit 42; else exit 0;fi',
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => {
           '1->A' => ['fan_lastz'],
           'A->1' => ['create_projection_db'],
        },
     },



      {
        -logic_name => 'create_projection_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('projection_db'),
                         create_type => 'clone',
                         #script_path => $self->o('clone_db_script_path'),
                         #user_r => $self->o('user_r'),
                         #user_w => $self->o('user_w'),
                         #pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
          '1->A' => ['cesar_create_projection_input_ids','wga_create_projection_input_ids'],
          'A->1' => ['classify_projected_genes'],
        },
      },

      {
        -logic_name => 'wga_create_projection_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          target_db           => $self->o('projection_source_db'),
          iid_type            => 'feature_id',
          feature_type        => 'transcript',
          feature_restriction => 'projection',
          biotypes            => {
            'protein_coding' => 1,
          },
          batch_size          => 100,
        },
        -rc_name    => '4GB',
        -flow_into => {
          2 => ['wga_project_transcripts'],
        },
      },


      {
        -logic_name => 'wga_project_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
        -parameters => {
          logic_name => 'project_transcripts',
          module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
          source_dna_db         => $self->default_options()->{'projection_source_db'},
          target_dna_db         => $self->o('dna_db'),
          source_transcript_db  => $self->default_options()->{'projection_source_db'},
          target_transcript_db   => $self->o('projection_db'),
          compara_db            => $self->o('projection_lastz_db'),
          method_link_type => 'LASTZ_NET',
          max_exon_readthrough_dist => 15,
          TRANSCRIPT_FILTER => {
            OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
            PARAMETERS => {
              -coverage => $self->o('projection_cov'),
              -percent_id => $self->o('projection_pid'),
            },
          },
          iid_type => 'feature_id',
          feature_type => 'transcript',
          calculate_coverage_and_pid => $self->o('projection_calculate_coverage_and_pid'),
          max_internal_stops => $self->o('projection_max_internal_stops'),
          timer => '30m',
        },
        -flow_into => {
          -3 => ['wga_failed_projection'],
        },
        -rc_name    => 'project_transcripts',
        -batch_size => 100,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'wga_failed_projection',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {},
        -rc_name    => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

      {
        -logic_name => 'cesar_create_projection_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('projection_source_db'),
                         iid_type => 'feature_id',
                         batch_size => 20,
                         feature_type => 'gene',
                         feature_id_include_non_reference => 0,
                         #feature_restriction => 'projection',
                         feature_restriction => 'biotype',
                         biotypes => {'protein_coding' => 1},
                       },

        -flow_into => {
                        2 => ['cesar'],
                      },
         -rc_name    => '5GB',
      },

      {
        -logic_name => 'cesar',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar',
        -parameters => {
                         'output_path' => $self->o('output_path')."/cesar_projection/",
                         'source_dna_db' => $self->default_options()->{'projection_source_db'},
                         'target_dna_db' => $self->o('dna_db'),
                         'source_db' => $self->o('projection_source_db'),
                         'target_db' => $self->o('projection_db'),
                         'compara_db' => $self->o('projection_lastz_db'),
                         'method_link_type' => 'LASTZ_NET',
                         'cesar_path' => $self->o('cesar_path'),
                         'cesar_mem' => '3', # mem in GB to be used by cesar (parameter --max-memory)
                         #TRANSCRIPT_FILTER => {
                         #  OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                         #  PARAMETERS => {
                         #    -coverage => 50,
                         #    -percent_id => 50,
                         #  },
                         #},
                         #'canonical' => 1,
                         'canonical_or_longest' => 1,
                         'stops2introns' => 1,
                       },
        -rc_name    => '3GB',
        -analysis_capacity => 50,
        -max_retry_count => 1,
        -failed_job_tolerance => 5,
        -flow_into => {
                        15 => ['cesar_15'],
                        30 => ['cesar_30'],
                        -1 => ['cesar_30'], # some jobs with batches still fail with -1, try cesar_30 in these cases
                      },
      },

      {
        -logic_name => 'cesar_15',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar',
        -parameters => {
                         'output_path' => $self->o('output_path')."/cesar_projection/",
                         'source_dna_db' => $self->default_options()->{'projection_source_db'},
                         'target_dna_db' => $self->o('dna_db'),
                         'source_db' => $self->o('projection_source_db'),
                         'target_db' => $self->o('projection_db'),
                         'compara_db' => $self->o('projection_lastz_db'),
                         'method_link_type' => 'LASTZ_NET',
                         'cesar_path' => $self->o('cesar_path'),
                         'cesar_mem' => '15', # mem in GB to be used by cesar (parameter --max-memory)
                         #TRANSCRIPT_FILTER => {
                         #  OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                         #  PARAMETERS => {
                         #    -coverage => 50,
                         #    -percent_id => 50,
                         #  },
                         #},
                         #'canonical' => 1,
                         'canonical_or_longest' => 1,
                         'stops2introns' => 1,
                       },
        -rc_name    => '15GB',
        -analysis_capacity => 50,
        -max_retry_count => 1,
        -can_be_empty  => 1,
        -failed_job_tolerance => 5,
        -flow_into => {
                        -1 => ['cesar_30'],
                      },
      },

      {
        -logic_name => 'cesar_30',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar',
        -parameters => {
                         'output_path' => $self->o('output_path')."/cesar_projection/",
                         'source_dna_db' => $self->default_options()->{'projection_source_db'},
                         'target_dna_db' => $self->o('dna_db'),
                         'source_db' => $self->o('projection_source_db'),
                         'target_db' => $self->o('projection_db'),
                         'compara_db' => $self->o('projection_lastz_db'),
                         'method_link_type' => 'LASTZ_NET',
                         'cesar_path' => $self->o('cesar_path'),
                         'cesar_mem' => '30', # mem in GB to be used by cesar (parameter --max-memory)
                         #TRANSCRIPT_FILTER => {
                         #  OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                         #  PARAMETERS => {
                         #    -coverage => 50,
                         #    -percent_id => 50,
                         #  },
                         #},
                         #'canonical' => 1,
                         'canonical_or_longest' => 1,
                         'stops2introns' => 1,
                       },
        -rc_name    => '30GB',
        -analysis_capacity => 50,
        -max_retry_count => 1,
        -can_be_empty  => 1,
        -failed_job_tolerance => 10,
        -flow_into => {
                        -1 => ['cesar_failed_projection'],
                      },
      },

      {
        -logic_name => 'cesar_failed_projection',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -failed_job_tolerance => 100,
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'classify_projected_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
          skip_analysis => $self->o('skip_projection'),
          classification_type => 'standard',
          update_gene_biotype => 1,
          target_db => $self->o('projection_db'),
        },
        -rc_name    => '2GB',
        -flow_into => {
          1 => ['flag_problematic_projections'],
#          When the realign part is fix, the line above needs to be deleted and the lines below uncommented
#          All analysis below create_projection_realign_db need to be uncommented too
#          '1->A' => ['fix_unaligned_protein_hit_names'],
#          'A->1' => ['create_projection_realign_db'],
        },
      },

      {
        -logic_name => 'flag_problematic_projections',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('flag_potential_pseudogenes_script').
                                ' -host '.$self->o('projection_db','-host').
                                ' -port '.$self->o('projection_db','-port').
                                ' -user_w '.$self->o('projection_db','-user').
                                ' -pass '.$self->o('projection_db','-pass').
                                ' -dbname '.$self->o('projection_db','-dbname').
                                ' -dna_host '.$self->o('dna_db','-host').
                                ' -dna_port '.$self->o('dna_db','-port').
                                ' -user_r '.$self->o('dna_db','-user').
                                ' -dna_dbname '.$self->o('dna_db','-dbname'),
                       },
        -rc_name => '2GB',
        -flow_into  => {
          1 => ['fix_projection_db_issues'],
        },
      },

      {
        # This will fix issues when proteins that were too long for MUSCLE alignment didn't have proper ENST accessions
        # Will also update the transcript table to match the gene table once the pseudo/canon genes are tagged in the
        # previous analysis
        -logic_name => 'fix_projection_db_issues',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn    => $self->o('projection_db'),
          sql => [
            'UPDATE gene JOIN transcript USING(gene_id) SET transcript.biotype=gene.biotype',
            'UPDATE protein_align_feature JOIN transcript_supporting_feature ON feature_id = protein_align_feature_id'.
              ' JOIN transcript USING(transcript_id) SET hit_name = stable_id',
            'UPDATE protein_align_feature JOIN supporting_feature ON feature_id = protein_align_feature_id'.
              ' JOIN exon_transcript USING(exon_id) JOIN transcript USING(transcript_id) SET hit_name = stable_id',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
          1 => ['projection_sanity_checks'],
        },
      },

      {
        -logic_name => 'projection_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
          target_db => $self->o('projection_db'),
          sanity_check_type => 'gene_db_checks',
          min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
            'gene_db_checks')->{$self->o('uniprot_set')}->{'projection_coding'},
        },
        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['select_projected_genes'],
                      },
      },

      {
        -logic_name => 'select_projected_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectProjectedGenes',
        -parameters => {
          wga_db   => $self->o('projection_db'),
          cesar_db   => $self->o('projection_db'),
          output_db => $self->o('selected_projection_db'),
          dna_db => $self->o('dna_db'),
        },
        -rc_name    => '4GB',
#        -flow_into => {
#                        1 => ['cluster_ig_tr_genes'],
#                      },
      },

#      When the realign part is fixed, uncomment the whole block
#      {
#        -logic_name => 'create_projection_realign_db',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
#        -parameters => {
#          source_db => $self->o('dna_db'),
#          target_db => $self->o('projection_realign_db'),
#          create_type => 'clone',
#        },
#        -rc_name    => 'default',
#        -flow_into => {
#          '1->A' => ['create_ids_for_evaluate_projection'],
#          'A->1' => ['update_realign_biotypes'],
#        },
#      },
#
#
#      {
#        -logic_name => 'create_ids_for_evaluate_projection',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#          target_db => $self->o('selected_projection_db'),
#          iid_type => 'feature_id',
#          feature_type => 'transcript',
#          batch_size => 500,
#        },
#        -rc_name    => 'default',
#        -flow_into => {
#          2 => ['evaluate_coding_transcripts'],
#        },
#      },
#
#
#      {
#        -logic_name => 'evaluate_coding_transcripts',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionRealignment',
#        -parameters => {
#          dna_db => $self->o('dna_db'),
#          projection_db => $self->o('selected_projection_db'),
#          projection_source_db => $self->o('projection_source_db'),
#          projection_realign_db => $self->o('projection_realign_db'),
#          protein_table_name => $self->o('realign_table_name'),
#          max_feature_issues => $self->o('max_projection_structural_issues'),
#        },
#        -rc_name    => 'default',
#        -flow_into  => {
#          2 => ['realign_projection'],
#        },
#      },
#
#
#      {
#        -logic_name => 'realign_projection',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
#        -parameters => {
#          iid_type => 'projection_transcript_id',
#          sequence_table_name => 'projection_source_sequences',
#          projection_padding => 0,
#          dna_db => $self->o('dna_db'),
#          target_db => $self->o('projection_realign_db'),
#          projection_db => $self->o('selected_projection_db'),
#          logic_name => 'genblast',
#          module => 'HiveGenblast',
#          genblast_path => $self->o('genblast_path'),
#          commandline_params => $genblast_params{$self->o('blast_type').'_projection'},
#          query_seq_dir => $self->o('homology_models_path'),
#          max_rank => 1,
#          genblast_pid => $self->o('genblast_pid'),
#          blast_eval => $self->o('genblast_eval'),
#          blast_cov => $self->o('genblast_cov'),
#          timer => '30m',
#        },
#        -rc_name    => 'default',
#        -flow_into  => {
#          -1 => ['failed_realignments'],
#          -2 => ['failed_realignments'],
#          -3 => ['failed_realignments'],
#        },
#      },
#
#
#      {
#        -logic_name => 'failed_realignments',
#        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
#        -parameters => {},
#        -rc_name    => 'default',
#        -can_be_empty  => 1,
#      },
#
#
#      {
#        -logic_name => 'update_realign_biotypes',
#        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
#        -parameters => {
#          db_conn    => $self->o('projection_realign_db'),
#          sql => [
#            'UPDATE gene SET biotype = "realign"',
#            'UPDATE transcript SET biotype = "realign"',
#          ],
#        },
#        -rc_name    => 'default',
#        -flow_into => {
#          1 => ['classify_realigned_models'],
#        },
#      },
#
#
#      {
#        -logic_name => 'classify_realigned_models',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
#        -parameters => {
#          classification_type => 'standard',
#          update_gene_biotype => 1,
#          target_db => $self->o('projection_realign_db'),
#        },
#        -rc_name    => 'default',
#        -flow_into => {
#          1 => ['realign_sanity_checks'],
#        },
#      },
#
#
#      {
#        -logic_name => 'realign_sanity_checks',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
#        -parameters => {
#                         target_db => $self->o('projection_realign_db'),
#                         sanity_check_type => 'gene_db_checks',
#                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
#                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'realign'},
#                       },
#
#        -rc_name    => '4GB',
#        -flow_into => {
#                       1 => ['create_projection_lincrna_db'],
#                      },
#      },
#
#
#      {
#        -logic_name => 'create_projection_lincrna_db',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
#        -parameters => {
#                         source_db => $self->o('dna_db'),
#                         target_db => $self->o('projection_lincrna_db'),
#                         create_type => 'clone',
#                       },
#        -flow_into => {
#                       1 => ['create_projection_pseudogene_db'],
#                      },
#      },
#
#
#      {
#        -logic_name => 'create_projection_pseudogene_db',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
#        -parameters => {
#                         source_db => $self->o('dna_db'),
#                         target_db => $self->o('projection_pseudogene_db'),
#                         create_type => 'clone',
#                       },
#
#        -flow_into => {
#                        1 => ['generate_projection_lincrna_ids','generate_projection_pseudogene_ids','generate_projection_ig_tr_ids'],
#                      },
#      },
#
#
#      {
#        -logic_name => 'generate_projection_lincrna_ids',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#                         target_db => $self->o('projection_source_db'),
#                         iid_type => 'feature_id',
#                         feature_type => 'transcript',
#                         feature_restriction => 'biotype',
#                         biotypes => {
#                                       'lincRNA' => 1,
#                                     },
#                         batch_size => 100,
#        },
#        -flow_into => {
#                       '2->A' => ['run_project_lincrnas'],
#                       'A->1' => ['project_lincrna_sanity_checks'],
#                      },
#        -rc_name    => 'default',
#      },
#
#
#      {
#        -logic_name => 'run_project_lincrnas',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate',
#        -parameters => {
#                         'logic_name'                 => 'project_lincrna',
#                         'source_dna_db'              => $self->o('projection_source_db'),
#                         'target_dna_db'              => $self->o('dna_db'),
#                         'source_db'                  => $self->o('projection_source_db'),
#                         'target_db'                  => $self->o('projection_lincrna_db'),
#                         'compara_db'                 => $self->o('projection_lastz_db'),
#                         'method_link_type'           => 'LASTZ_NET',
#                         'exon_region_padding'        => $self->o('projection_exonerate_padding'),
#                         'exonerate_path'             => $self->o('exonerate_path'),
#                         'exonerate_coverage'         => $self->o('projection_lincrna_coverage'),
#                         'exonerate_percent_id'       => $self->o('projection_lincrna_percent_id'),
#                         'calculate_coverage_and_pid' => 0,
#                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_projection_dna')},
#			},
#        -rc_name    => 'default',
#        -hive_capacity => 900,
#      },
#
#
#      {
#        -logic_name => 'project_lincrna_sanity_checks',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
#        -parameters => {
#                         target_db => $self->o('projection_lincrna_db'),
#                         sanity_check_type => 'gene_db_checks',
#                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
#                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'projection_lincrna'},
#                       },
#        -rc_name    => '4GB',
#      },
#
#
#      {
#        -logic_name => 'generate_projection_pseudogene_ids',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#                         target_db => $self->o('projection_source_db'),
#                         iid_type => 'feature_id',
#                         feature_type => 'transcript',
#                         feature_restriction => 'biotype',
#                         biotypes => {
#                                       'pseudogene' => 1,
#                                       'processed_pseudogene' => 1,
#                                       'unprocessed_pseudogene' => 1,
#                                     },
#                         batch_size => 100,
#        },
#        -flow_into => {
#                       '2->A' => ['run_project_pseudogenes'],
#                       'A->1' => ['project_pseudogene_sanity_checks'],
#        },
#        -rc_name    => 'default',
#      },
#
#
#      {
#        -logic_name => 'run_project_pseudogenes',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate',
#        -parameters => {
#                         'logic_name'                 => 'project_pseudogene',
#                         'source_dna_db'              => $self->o('projection_source_db'),
#                         'target_dna_db'              => $self->o('dna_db'),
#                         'source_db'                  => $self->o('projection_source_db'),
#                         'target_db'                  => $self->o('projection_pseudogene_db'),
#                         'compara_db'                 => $self->o('projection_lastz_db'),
#                         'method_link_type'           => 'LASTZ_NET',
#                         'exon_region_padding'        => $self->o('projection_exonerate_padding'),
#                         'exonerate_path'             => $self->o('exonerate_path'),
#                         'exonerate_coverage'         => $self->o('projection_pseudogene_coverage'),
#                         'exonerate_percent_id'       => $self->o('projection_pseudogene_percent_id'),
#                         'calculate_coverage_and_pid' => 0,
#                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_projection_pseudogene')},
#			},
#        -rc_name    => 'default',
#        -hive_capacity => 900,
#      },
#
#
#      {
#        -logic_name => 'project_pseudogene_sanity_checks',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
#        -parameters => {
#                         target_db => $self->o('projection_pseudogene_db'),
#                         sanity_check_type => 'gene_db_checks',
#                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
#                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'projection_pseudogene'},
#                       },
#        -rc_name    => '4GB',
#      },
#
#
#      {
#        -logic_name => 'generate_projection_ig_tr_ids',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#                         target_db => $self->o('projection_source_db'),
#                         iid_type => 'feature_id',
#                         feature_type => 'transcript',
#                         feature_restriction => 'biotype',
#                         biotypes => {
#                                       'IG_V_gene' => 1,
#                                       'IG_D_gene' => 1,
#                                       'IG_J_gene' => 1,
#                                       'IG_C_gene' => 1,
#                                       'TR_V_gene' => 1,
#                                       'TR_J_gene' => 1,
#                                       'TR_D_gene' => 1,
#                                       'TR_C_gene' => 1,
#                                     },
#                         batch_size => 100,
#        },
#        -flow_into => {
#                       '2->A' => ['run_project_ig_tr'],
#                       'A->1' => ['project_ig_tr_sanity_checks'],
#        },
#        -rc_name    => 'default',
#      },
#
#
#      {
#        -logic_name => 'run_project_ig_tr',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate',
#        -parameters => {
#                         'logic_name'                 => 'project_ig_tr',
#                         'source_dna_db'              => $self->o('projection_source_db'),
#                         'target_dna_db'              => $self->o('dna_db'),
#                         'source_db'                  => $self->o('projection_source_db'),
#                         'target_db'                  => $self->o('ig_tr_db'),
#                         'compara_db'                 => $self->o('projection_lastz_db'),
#                         'method_link_type'           => 'LASTZ_NET',
#                         'exon_region_padding'        => $self->o('projection_exonerate_padding'),
#                         'exonerate_path'             => $self->o('exonerate_path'),
#                         'exonerate_coverage'         => $self->o('projection_ig_tr_coverage'),
#                         'exonerate_percent_id'       => $self->o('projection_ig_tr_percent_id'),
#                         'calculate_coverage_and_pid' => 1,
#                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_projection_ig_tr_protein')},
#			},
#        -rc_name    => 'default',
#        -hive_capacity => 900,
#      },
#
#
#      {
#        -logic_name => 'project_ig_tr_sanity_checks',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
#        -parameters => {
#                         target_db => $self->o('ig_tr_db'),
#                         sanity_check_type => 'gene_db_checks',
#                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
#                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'projection_ig_tr'},
#                       },
#        -rc_name    => '2GB',
#      },



######################################################################################
#
# LastZ pipeline
#
######################################################################################

	    {
        -logic_name => 'fan_lastz',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_lastz#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => '2GB_lastz',
        -flow_into  => {
          '1' => ['setup_lastz'],
        },
      },




	    {
        -logic_name => 'setup_lastz',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::LastZSetup',
        -parameters => {
                         compara_genome_db_update_path => $self->o('compara_genome_db_update_path'),
                         compara_mlss_script_path => $self->o('compara_mlss_script_path'),
                         compara_mlss_reg_conf_path => $self->o('compara_mlss_reg_conf_path'),
                         compara_db => $self->o('compara_db'),
                         projection_source_db => $self->o('projection_source_db'),
                         target_db => $self->o('dna_db'),
                         pipeline_db => $self->o('pipeline_db'),
                         output_path => $self->o('output_path'),
                         compara_db_url => 'mysql://'.$self->o('user').':'.$self->o('password').'@'.$self->o('compara_db_server').':'.$self->o('compara_db_port').'/'.$self->o('compara_db_name'),
                       },

        -rc_name => '2GB_lastz',
             -flow_into      => {
                             1 => {'get_species_list' => {'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#'}},
                           },
       },


	    {
         -logic_name    => 'get_species_list',
         -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ParsePairAlignerConf',
         -parameters    => {
                            'master_db' => $self->o('compara_master'),
                            'get_species_list' => 1,
                            'core_dbs' => undef,
                          },

         -rc_name => '2GB_lastz',
         -flow_into      => {
                              1 => {'populate_new_database' => {'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#'}},
                            },
       },


	    {
         -logic_name => 'populate_new_database',
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::PopulateNewDatabase',
         -parameters    => {
                             'program'        => $self->o('compara_populate_new_database_exe'),
                             'mlss_id_list'   => undef,
                             'collection'     => $self->o('compara_collection'),
                             'master_db'      => $self->o('compara_master'),
                             'only_cellular_component' => $self->o('only_cellular_component'),
                           },
         -flow_into => {
           1 => {'add_method_link_species_link_tag' => {'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#'}},
         },
         -rc_name => '2GB_lastz',
       },


	    {
        -logic_name => 'add_method_link_species_link_tag',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('pipeline_db'),
          sql => [
            'INSERT INTO method_link_species_set_tag VALUES (#mlss_id#,"reference_species","'.$self->o('compara_ref_species').'")',
          ],
        },
        -rc_name    => '2GB_lastz',
        -flow_into => {
          1 => {'parse_pair_aligner_conf' => {'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#'}}
        },
      },


	    {
         -logic_name    => 'parse_pair_aligner_conf',
         -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ParsePairAlignerConf',
         -parameters    => {
                             'conf_file' => $self->o('compara_conf_file'),
                             'ref_species' => $self->o('compara_ref_species'),
                             'non_ref_species' => $self->o('compara_non_ref_species'),
                             'dump_dir' => $self->o('compara_dump_dir'),
                             'default_chunks' => $self->o('default_chunks'),
                             'default_pair_aligner' => $self->o('pair_aligner_method_link'),
                             'default_parameters' => $self->o('pair_aligner_options'),
                             'default_chain_output' => $self->o('chain_output_method_link'),
                             'default_net_output' => $self->o('net_output_method_link'),
                             'default_chain_input' => $self->o('chain_input_method_link'),
                             'default_net_input' => $self->o('net_input_method_link'),
                             'net_ref_species' => $self->o('net_ref_species'),
                             'mlss_id_list' => $self->o('mlss_id_list'),
                             'collection' => $self->o('compara_collection'),
                             'master_db' => $self->o('compara_master'),
                             'do_pairwise_gabs' => $self->o('do_pairwise_gabs'), #healthcheck options
                             'do_compare_to_previous_db' => $self->o('do_compare_to_previous_db'), #healthcheck options
                             'bidirectional' => $self->o('bidirectional'),
                           },
         -flow_into => {
                         1 => [ 'create_pair_aligner_jobs'],
                         2 => [ 'chunk_and_group_dna' ],
                         3 => [ 'create_filter_duplicates_jobs' ],
                         4 => [ 'no_chunk_and_group_dna' ],
                         5 => [ 'create_alignment_chains_jobs' ],
                         6 => [ 'create_alignment_nets_jobs' ],
                        10 => [ 'create_filter_duplicates_net_jobs' ],
                         7 => [ 'pairaligner_stats' ],
                         8 => [ 'healthcheck' ],
                       },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => 'chunk_and_group_dna',
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ChunkAndGroupDna',
         -parameters => {
                          'only_cellular_component' => $self->o('only_cellular_component'),
                          'mix_cellular_components' => $self->o('mix_cellular_components'),
                        },
         -flow_into => {
           2 => [ 'store_sequence' ],
         },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => 'store_sequence',
         -hive_capacity => 100,
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::StoreSequence',
         -parameters => {
           'dump_min_chunkset_size' => $self->o('dump_min_chunkset_size'),
           'dump_min_chunk_size' => $self->o('dump_min_chunk_size'),
         },
         -flow_into => {
           -1 => [ 'store_sequence_again' ],
         },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => 'store_sequence_again',
         -hive_capacity => 50,
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::StoreSequence',
         -parameters => {
           'dump_min_chunkset_size' => $self->o('dump_min_chunkset_size'),
           'dump_min_chunk_size' => $self->o('dump_min_chunk_size'),
          },
          -can_be_empty  => 1,
          -rc_name => '4GB_lastz',
       },


	    {
         -logic_name => 'create_pair_aligner_jobs',  #factory
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreatePairAlignerJobs',
         -parameters => {
                          'mix_cellular_components' => $self->o('mix_cellular_components'),
                        },
         -hive_capacity => 10,
         -wait_for => [ 'store_sequence', 'store_sequence_again', 'chunk_and_group_dna'  ],
         -flow_into => {
                         1 => [ 'check_no_partial_gabs' ],
                         2 => [ $self->o('pair_aligner_logic_name')  ],
                       },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => $self->o('pair_aligner_logic_name'),
         -module     => $self->o('pair_aligner_module'),
         -analysis_capacity => $self->o('pair_aligner_analysis_capacity'),
         -batch_size => $self->o('pair_aligner_batch_size'),
         -parameters => {
                          'pair_aligner_exe' => $self->o('pair_aligner_exe'),
                        },
         -wait_for  => [ 'create_pair_aligner_jobs'  ],
         -flow_into => {
                         -1 => [ $self->o('pair_aligner_logic_name') . '_himem1' ],  # MEMLIMIT
                       },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => $self->o('pair_aligner_logic_name') . "_himem1",
         -module     => $self->o('pair_aligner_module'),
         -analysis_capacity => $self->o('pair_aligner_analysis_capacity'),
         -parameters => {
                          'pair_aligner_exe' => $self->o('pair_aligner_exe'),
                        },
         -wait_for   => [ 'create_pair_aligner_jobs'  ],
         -batch_size => $self->o('pair_aligner_batch_size'),
         -can_be_empty  => 1,
         -rc_name => '8GB_lastz',
       },


	    {
         -logic_name => 'check_no_partial_gabs',
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SqlHealthChecks',
         -parameters => {
           'mode'          => 'gab_inconsistencies',
         },
         -wait_for =>  [ $self->o('pair_aligner_logic_name'), $self->o('pair_aligner_logic_name') . "_himem1" ],
         -flow_into => {
           1 => [ 'update_max_alignment_length_before_FD' ],
         },
       },


	    {
         -logic_name => 'update_max_alignment_length_before_FD',
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
         -parameters => {
           'quick' => $self->o('quick'),
         },
         -flow_into => {
           1 => [ 'update_max_alignment_length_after_FD' ],
         },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => 'create_filter_duplicates_jobs', #factory
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateFilterDuplicatesJobs',
         -parameters => { },
         -wait_for =>  [ 'update_max_alignment_length_before_FD' ],
         -flow_into => {
           2 => { 'filter_duplicates' => INPUT_PLUS() },
         },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name   => 'filter_duplicates',
         -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
         -parameters    => {
           'window_size' => $self->o('compara_window_size')
         },
         -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
         -batch_size    => $self->o('filter_duplicates_batch_size'),
         -can_be_empty  => 1,
         -flow_into => {
                         -1 => [ 'filter_duplicates_himem' ], # MEMLIMIT
         },
         -rc_name => $self->o('filter_duplicates_rc_name'),
       },


	    {
        -logic_name   => 'filter_duplicates_himem',
        -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
        -parameters    => {
                            'window_size' => $self->o('compara_window_size')
         },
         -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
         -batch_size    => $self->o('filter_duplicates_batch_size'),
         -can_be_empty  => 1,
         -rc_name => $self->o('filter_duplicates_himem_rc_name'),
      },


	    {
        -logic_name => 'update_max_alignment_length_after_FD',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
        -parameters => {
          'quick' => $self->o('quick'),
        },
        -wait_for =>  [ 'create_filter_duplicates_jobs', 'filter_duplicates', 'filter_duplicates_himem' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'no_chunk_and_group_dna',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ChunkAndGroupDna',
        -parameters => {
          'only_cellular_component' => $self->o('only_cellular_component'),
          'mix_cellular_components' => $self->o('mix_cellular_components'),
        },
        -flow_into => {
          2 => [ 'dump_large_nib_for_chains' ],
        },
        -wait_for  => ['update_max_alignment_length_after_FD' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'dump_large_nib_for_chains',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::DumpDnaCollection',
        -parameters => {
          'faToNib_exe' => $self->o('faToNib_exe'),
          'dump_min_nib_size' => $self->o('dump_min_nib_size'),
          'overwrite'=>1,
        },
        -hive_capacity => 10,
        -flow_into => {
          -1 => [ 'dump_large_nib_for_chains_himem' ],  # MEMLIMIT
        },
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'dump_large_nib_for_chains_himem',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::DumpDnaCollection',
        -parameters => {
                         'faToNib_exe' => $self->o('faToNib_exe'),
                         'dump_min_nib_size' => $self->o('dump_min_nib_size'),
                         'overwrite'=>1,
        },
        -hive_capacity => 10,
        -can_be_empty  => 1,
        -rc_name => '8GB_lastz',
      },



	    {
        -logic_name => 'create_alignment_chains_jobs',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateAlignmentChainsJobs',
        -parameters => { },
        -flow_into => {
                        1 => [ 'remove_inconsistencies_after_chain' ],
                        2 => [ 'alignment_chains' ],
        },
        -wait_for => [ 'no_chunk_and_group_dna', 'dump_large_nib_for_chains', 'dump_large_nib_for_chains_himem' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'alignment_chains',
        -hive_capacity => $self->o('chain_hive_capacity'),
        -batch_size => $self->o('chain_batch_size'),
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
        -parameters => $self->o('chain_parameters'),
        -max_retry_count => 10,
        -flow_into => {
          -1 => [ 'alignment_chains_himem' ],  # MEMLIMIT
        },
        -wait_for   => [ 'create_alignment_chains_jobs' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'alignment_chains_himem',
        -hive_capacity => $self->o('chain_hive_capacity'),
        -batch_size => 1,
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
        -parameters => $self->o('chain_parameters'),
        -can_be_empty  => 1,
        -max_retry_count => 10,
        -rc_name => '8GB_lastz',
        -flow_into => {
          -1 => [ 'alignment_chains_super_himem' ],  # MEMLIMIT
         },
         -wait_for => ['alignment_chains'],
         -can_be_empty  => 1,
      },


	    {
        -logic_name => 'alignment_chains_super_himem',
        -hive_capacity => $self->o('chain_hive_capacity'),
        -batch_size => 1,
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
        -parameters => $self->o('chain_parameters'),
        -can_be_empty  => 1,
        -max_retry_count => 10,
        -rc_name => '15GB_lastz',
        -wait_for => ['alignment_chains'],
         -can_be_empty  => 1,
      },

	    {
        -logic_name => 'remove_inconsistencies_after_chain',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::RemoveAlignmentDataInconsistencies',
        -flow_into => {
          1 => [ 'update_max_alignment_length_after_chain' ],
        },
        -wait_for =>  [ 'alignment_chains', 'alignment_chains_himem', 'alignment_chains_super_himem' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'update_max_alignment_length_after_chain',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
        -parameters => {
          'quick' => $self->o('quick'),
        },
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'create_alignment_nets_jobs',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateAlignmentNetsJobs',
        -parameters => { },
        -flow_into => {
          1 => [ 'remove_inconsistencies_after_net' ],
          2 => [ 'alignment_nets' ],
        },
        -wait_for => [ 'update_max_alignment_length_after_chain', 'remove_inconsistencies_after_chain' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'alignment_nets',
        -hive_capacity => $self->o('net_hive_capacity'),
        -batch_size => $self->o('net_batch_size'),
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
        -parameters => $self->o('net_parameters'),
        -flow_into => {
          -1 => [ 'alignment_nets_himem' ],  # MEMLIMIT
         },
         -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'alignment_nets_himem',
        -hive_capacity => $self->o('net_hive_capacity'),
        -batch_size => $self->o('net_batch_size'),
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
        -parameters => $self->o('net_parameters'),
        -can_be_empty => 1,
        -flow_into => {
          -1 => [ 'alignment_nets_hugemem' ],  # MEMLIMIT
        },
        -rc_name => '4GB_lastz',
      },


	    {
        -logic_name     => 'alignment_nets_hugemem',
        -hive_capacity  => $self->o('net_hive_capacity'),
        -batch_size     => $self->o('net_batch_size'),
        -module         => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
        -parameters     => $self->o('net_parameters'),
        -can_be_empty   => 1,
        -rc_name        => '8GB_lastz',
      },


	    {
        -logic_name => 'remove_inconsistencies_after_net',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::RemoveAlignmentDataInconsistencies',
        -flow_into => {
          1 => [ 'update_max_alignment_length_after_net' ],
        },
        -wait_for =>  [ 'alignment_nets', 'alignment_nets_himem', 'alignment_nets_hugemem', 'create_alignment_nets_jobs' ],    # Needed because of bi-directional netting: 2 jobs in create_alignment_nets_jobs can result in 1 job here
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'create_filter_duplicates_net_jobs', #factory
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateFilterDuplicatesJobs',
        -parameters => { },
        -wait_for =>  [ 'remove_inconsistencies_after_net' ],
        -flow_into => {
          2 => { 'filter_duplicates_net' => INPUT_PLUS() },
        },
        -can_be_empty  => 1,
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name   => 'filter_duplicates_net',
        -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
        -parameters    => {
          'window_size' => $self->o('compara_window_size'),
          'filter_duplicates_net' => 1,
        },
        -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
        -batch_size    => $self->o('filter_duplicates_batch_size'),
        -flow_into => {
          -1 => [ 'filter_duplicates_net_himem' ], # MEMLIMIT
        },
        -can_be_empty  => 1,
        -rc_name => $self->o('filter_duplicates_rc_name'),
      },


	    {
        -logic_name   => 'filter_duplicates_net_himem',
        -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
        -parameters    => {
          'window_size' => $self->o('compara_window_size'),
          'filter_duplicates_net' => 1,
        },
       -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
       -batch_size    => $self->o('filter_duplicates_batch_size'),
       -can_be_empty  => 1,
       -rc_name => $self->o('filter_duplicates_himem_rc_name'),
      },


	    {
        -logic_name => 'update_max_alignment_length_after_net',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
        -rc_name => '2GB_lastz',
        -wait_for =>  [ 'create_filter_duplicates_net_jobs', 'filter_duplicates_net', 'filter_duplicates_net_himem' ],
        -flow_into => [ 'set_internal_ids_collection' ],
      },

	    {
        -logic_name => 'set_internal_ids_collection',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SetInternalIdsCollection',
        -parameters => {
          'skip' => $self->o('patch_alignments'),
         },
         -flow_into => {
           2 => [ 'set_internal_ids_slow' ],
          },
         -analysis_capacity => 1,
      },


	    {
        -logic_name => 'set_internal_ids_slow',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SetInternalIdsSlow',
        -analysis_capacity => 1,
        -rc_name => '8GB_lastz',
      },

	    {
        -logic_name => 'healthcheck',
        -module => 'Bio::EnsEMBL::Compara::RunnableDB::HealthCheck',
        -parameters => {
          'previous_db' => $self->o('previous_db'),
          'ensembl_release' => $self->o('ensembl_release'),
          'prev_release' => $self->o('prev_release'),
          'max_percent_diff' => $self->o('patch_alignments') ? $self->o('max_percent_diff_patches') : $self->o('max_percent_diff'),
        },
        -wait_for => [ 'set_internal_ids_collection' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'pairaligner_stats',
        -module => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::PairAlignerStats',
        -parameters => {
          'skip' => $self->o('skip_pairaligner_stats'),
          'dump_features' => $self->o('dump_features_exe'),
          'compare_beds' => $self->o('compare_beds_exe'),
          'create_pair_aligner_page' => $self->o('create_pair_aligner_page_exe'),
          'bed_dir' => $self->o('compara_bed_dir'),
          'ensembl_release' => $self->o('ensembl_release'),
          'output_dir' => $self->o('compara_feature_dir'),
        },
        -wait_for =>  [ 'healthcheck' ],
        -rc_name => '2GB_lastz',
      },



######################################################################################
#
# Protein models (genblast and genewise)
#
######################################################################################

      {
        -logic_name => 'create_genblast_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -wait_for => ['create_softmasked_faidx'],
        -flow_into => {
                        '1->A' => ['download_uniprot_files'],
                        'A->1' => ['classify_genblast_models'],
                      },
      },


      {

        -logic_name => 'download_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         multi_query_download => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', $self->o('uniprot_set')),
                         taxon_id => $self->o('taxon_id'),
                         output_path => $self->o('homology_models_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        '2->A' => ['process_uniprot_files'],
                        'A->1' => ['generate_genblast_jobs'],
                      },
      },

      {
        -logic_name => 'process_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
        -parameters => {
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                         sequence_table_name => $self->o('uniprot_table_name'),
                      },
        -rc_name => 'default',
      },



      {
        -logic_name => 'generate_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'sequence_accession',
                         batch_size => $self->o('uniprot_genblast_batch_size'),
                         sequence_table_name => $self->o('uniprot_table_name'),
                       },
        -rc_name      => '3GB',
        -flow_into => {
                        2 => ['genblast'],
                        1 => ['create_seleno_homology_jobs'],
                      },
      },

      {
        -logic_name => 'create_seleno_homology_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputquery => 'SELECT accession FROM '.$self->o('uniprot_table_name').' WHERE source_db = "seleno"',
          column_names => ['iid'],
        },
        -rc_name          => 'default',
        -flow_into => {
          2 => ['process_homology_selenocysteine'],
        },
      },

      {
        -logic_name => 'process_homology_selenocysteine',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelenocysteineFinder',
        -parameters => {
          target_db => $self->o('genblast_db'),
          dna_db => $self->o('dna_db'),
          genome => $self->o('softmasked_genome_file'),
          exonerate => $self->o('exonerate_path'),
          genewise => $self->o('genewise_path'),
          iid_type => 'db_seq',
          sequence_table_name => $self->o('uniprot_table_name'),
          biotype => 'seleno_other',
          missmatch_allowed => 10,
        },
        -rc_name          => '3GB',
      },



      {
        -logic_name => 'genblast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('softmasked_genome_file'),
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '2h',
                         blast_eval => $self->o('genblast_eval'),
                         blast_cov  => $self->o('genblast_cov'),
                         flag_small_introns => $self->o('genblast_flag_small_introns'),
                         flag_subpar_models => $self->o('genblast_flag_subpar_models'),
                       },
        -rc_name    => 'genblast',
        -flow_into => {
                        -1 => ['split_genblast_jobs'],
                        -2 => ['split_genblast_jobs'],
                        -3 => ['split_genblast_jobs'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },

      {
        -logic_name => 'split_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk',
                         batch_size => 1,
                       },
        -rc_name      => 'default',
        -can_be_empty  => 1,
        -flow_into => {
                        2 => ['genblast_retry'],
                      },
      },

      {
        -logic_name => 'genblast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('softmasked_genome_file'),
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '1h',
                         blast_eval => $self->o('genblast_eval'),
                         blast_cov  => $self->o('genblast_cov'),
                         flag_small_introns => $self->o('genblast_flag_small_introns'),
                         flag_subpar_models => $self->o('genblast_flag_subpar_models'),
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
        -flow_into => {
                        -1 => ['failed_genblast_proteins'],
                        -2 => ['failed_genblast_proteins'],
                        -3 => ['failed_genblast_proteins'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },

      {
        -logic_name => 'failed_genblast_proteins',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },


      {
        -logic_name => 'classify_genblast_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                         target_db => $self->o('genblast_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['genblast_sanity_checks'],
                      },
      },


      {
        -logic_name => 'genblast_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('genblast_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'genblast'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        '1->A' => ['create_cdna_db','create_genblast_nr_db'],
                        'A->1' => ['create_ig_tr_db'],
                      },

      },


      {
        -logic_name => 'create_genblast_nr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('genblast_db'),
                         target_db => $self->o('genblast_nr_db'),
                         create_type => 'copy',
                         force_drop => 1,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['create_genblast_nr_slices'],
                      },
      },


      {
        -logic_name => 'create_genblast_nr_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 20000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                         batch_slice_ids => 1,
                         batch_target_size => 20000000,
                       },
        -rc_name    => '2GB',
        -flow_into => {
                         '2'    => ['remove_redundant_genblast_genes'],
                      },
      },


    {
      -logic_name => 'remove_redundant_genblast_genes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db => $self->o('genblast_nr_db'),
        target_type => 'biotype_priority',
        layers => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic', $self->o('uniprot_set'), undef, 'ARRAY'),
      },
      -rc_name          => '5GB',
    },

######################################################################################
#
# cDNA alignment
#
######################################################################################

    {
      -logic_name => 'create_cdna_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('cdna_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        '1->A' => ['create_genewise_db', 'download_mRNA'],
        'A->1' => ['create_besttargetted_db'],
      },
    },

    {
      -logic_name => 'create_genewise_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('genewise_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['download_selenocysteines'],
        '1->A' => ['download_uniprot_self', 'download_refseq_self'],
        'A->1' => ['load_self'],
      },
    },
    {

      -logic_name => 'download_refseq_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadNCBImRNA',
      -parameters => {
        query => '#species_name#[Organism] AND RefSeq[Filter]',
        species_name => $self->o('species_name'),
        output_file => catfile($self->o('targetted_path'), 'ncbi_self.fa'),
        ncbidb => 'protein',
        _branch_to_flow_to => 1,
      },
      -rc_name          => 'default',
      -flow_into => {
        1 => [':////accu?iid=[]'],
      },
    },
    {

      -logic_name => 'download_uniprot_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
      -parameters => {
        taxon_id => $self->o('taxon_id'),
        multi_query_download => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', 'self_isoforms_12'),
        output_path => $self->o('targetted_path'),
        _branch_to_flow_to => 1,
      },
      -rc_name          => 'default',
      -flow_into => {
        1 => [':////accu?iid=[]'],
      },
    },

    {

      -logic_name => 'load_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
      -parameters => {
        output_file => catfile($self->o('targetted_path'), 'proteome.fa'),
        skip_Xs => 5,
        killlist_db => $self->o('killlist_db'),
        killlist_type => 'protein',
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name          => 'default',
      -flow_into => {
        1 => ['indicate_proteome'],
        2 => ['targetted_exonerate'],
      },
    },
    {

      -logic_name => 'targetted_exonerate',
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
      },
      -rc_name          => '3GB',
    },


    {
      -logic_name => 'indicate_proteome',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => '#indicate_path# -d #indicate_dir# -f #proteome# -i #proteome_index# -p singleWordParser',
        indicate_path => $self->o('indicate_path'),
        proteome => 'proteome.fa',
        indicate_dir => $self->o('targetted_path'),
        proteome_index => catdir($self->o('targetted_path'), 'proteome_index'),
      },
      -rc_name => 'default',
      -flow_into => {
        '1' => ['generate_pmatch_jobs'],
      },
    },


    {
      -logic_name => 'generate_pmatch_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db        => $self->o('dna_db'),
        coord_system_name => 'toplevel',
        iid_type => 'slice',
        top_level => 1,
      },
      -rc_name      => '2GB',
      -flow_into => {
        '2->A' => ['pmatch'],
        'A->1' => ['bestpmatch'],
      },
    },
    {

      -logic_name => 'pmatch',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePmatch',
      -parameters => {
        target_db => $self->o('genewise_db'),
        dna_db => $self->o('dna_db'),
        PROTEIN_FILE => catfile($self->o('targetted_path'), 'proteome.fa'),
        MIN_COVERAGE => 25,
        BINARY_LOCATION => $self->o('pmatch_path'),
        REPEAT_MASKING => [],
        MAX_INTRON_LENGTH => 50000,
        OPTIONS => '-T 20', # set threshold to 14 for more sensitive search
      },
      -rc_name          => '2GB',
    },
    {

      -logic_name => 'bestpmatch',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestPmatch',
      -parameters => {
        source_db => $self->o('genewise_db'),
        target_db => $self->o('genewise_db'),
        PMATCH_LOGIC_NAME => ['pmatch'],
        MIN_COVERAGE => 50,
      },
      -rc_name          => '2GB',
      -flow_into => {
        1 => ['generate_targetted_jobs'],
      },
    },
    {
      -logic_name => 'generate_targetted_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type => 'feature_region',
        feature_type => 'protein_align_feature',
        target_db        => $self->o('genewise_db'),
        logic_name => ['bestpmatch'],
        coord_system_name => 'toplevel',
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['targetted_genewise_gtag', 'targetted_genewise_nogtag', 'targetted_exo'],
      },
    },
    {

      -logic_name => 'targetted_genewise_gtag',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise',
      -parameters => {
        source_db => $self->o('genewise_db'),
        target_db => $self->o('genewise_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('dna_db'),
        gtag => 0, # 0 is for gtag, 1 is for non canonical
        biotype => 'gw_gtag',
        max_intron_length => 200000,
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_genewise')},
      },
      -rc_name          => '3GB',
    },
    {

      -logic_name => 'targetted_genewise_nogtag',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise',
      -parameters => {
        source_db => $self->o('genewise_db'),
        target_db => $self->o('genewise_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('dna_db'),
        gtag => 1, # 0 is for gtag, 1 is for non canonical
        biotype => 'gw_nogtag',
        max_intron_length => 200000,
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_genewise')},
      },
      -rc_name          => '3GB',
    },
    {

      -logic_name => 'targetted_exo',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerateForGenewise',
      -parameters => {
        source_db => $self->o('genewise_db'),
        target_db => $self->o('genewise_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('dna_db'),
        biotype => 'gw_exo',
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        max_intron_length => 700000,
        program_file => $self->o('exonerate_path'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_exonerate')},
      },
      -rc_name          => '3GB',
    },
    {
      -logic_name => 'download_mRNA',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadNCBImRNA',
      -parameters => {
        output_file => $self->o('cdna_file'),
        filetype => 'gb',
        query => $self->o('ncbi_query'),
      },
      -rc_name    => 'default',
      -flow_into => {
        2 => ['prepare_cdna'],
      },
    },
    {
      -logic_name => 'prepare_cdna',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('prepare_cdnas_script').
          ' -killdbnames '.$self->o('killlist_db','-dbname').
          ' -killdbhost '.$self->o('killlist_db','-host').
          ' -killdbuser '.$self->o('killlist_db','-user').
          ' -killdbport '.$self->o('killlist_db','-port').
          ($self->o('killlist_db','-pass') ? ' -killdbpass '.$self->o('killlist_db','-pass') : '').
          ' -infile #sequence_file#'.
          ' -outfile #sequence_file#.clipped'.
          ($self->o('taxon_id') ? ' -tax_id '.$self->o('taxon_id') : '').
          ' -nomole',
        sequence_file => $self->o('cdna_file'),
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['load_cdna_file'],
      },
    },
    {
      -logic_name => 'load_cdna_file',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadmRNAs',
      -parameters => {
        sequence_table_name => $self->o('cdna_table_name'),
        filetype => 'fasta',
        sequence_file => $self->o('cdna_file').'.clipped',
      },
      -flow_into => {
        '2->A' => ['exonerate'],
        'A->1' => ['prepare_cdna2genome'],
      },
      -rc_name    => 'default',
    },

    {
      -logic_name => 'exonerate',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '3GB',
      -parameters => {
        iid_type => 'db_seq',
        sequence_table_name => $self->o('cdna_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('cdna_db'),
        logic_name => $self->o('exonerate_logic_name'),
        genome_file      => $self->o('softmasked_genome_file'),
        exonerate_path   => $self->o('exonerate_path'),
        repeat_libraries => '#wide_repeat_logic_names#',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna_est2genome')},
        exonerate_cdna_pid => 50,
        exonerate_cdna_cov => 50,
        calculate_coverage_and_pid => 0,
      },
      -batch_size => 100,
      -flow_into => {
        -1 => ['exonerate_retry'],
      },
      -batch_size => 100,
      -failed_job_tolerance => 5,
    },

    {
      -logic_name => 'exonerate_retry',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '6GB',
      -parameters => {
        iid_type => 'db_seq',
        sequence_table_name => $self->o('cdna_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('cdna_db'),
        logic_name => $self->o('exonerate_logic_name'),
        genome_file      => $self->o('softmasked_genome_file'),
        exonerate_path   => $self->o('exonerate_path'),
        repeat_libraries => '#wide_repeat_logic_names#',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna_est2genome')},
        exonerate_cdna_pid => 50,
        exonerate_cdna_cov => 50,
        calculate_coverage_and_pid => 0,
      },
      -batch_size => 100,
      -failed_job_tolerance => 100,
      -can_be_empty => 1,
    },
    {
      -logic_name => 'prepare_cdna2genome',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('prepare_cdnas_script').
          ' -killdbnames '.$self->o('killlist_db','-dbname').
          ' -killdbhost '.$self->o('killlist_db','-host').
          ' -killdbuser '.$self->o('killlist_db','-user').
          ' -killdbport '.$self->o('killlist_db','-port').
          ($self->o('killlist_db','-pass') ? ' -killdbpass '.$self->o('killlist_db','-pass') : '').
          ' -infile #sequence_file#'.
          ' -outfile #sequence_file#.cds'.
          ' -annotation '.$self->o('annotation_file').
          ($self->o('taxon_id') ? ' -tax_id '.$self->o('taxon_id') : '').
          ' -nomole',
        sequence_file => $self->o('cdna_file'),
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['create_cdna2genome_db'],
      },
    },
    {
      -logic_name => 'create_cdna2genome_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('cdna2genome_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        '1' => ['create_cdna_toplevel_slices'],
      },
    },
    {
      -logic_name => 'create_cdna_toplevel_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db => $self->o('cdna_db'),
        iid_type => 'slice',
        coord_system_name => 'toplevel',
        include_non_reference => 0,
        top_level => 1,
        feature_constraint => 1,
        feature_type => 'gene',
      },
      -flow_into => {
        '2' => ['apply_threshold'],
      },
      -rc_name    => 'default',
    },


    {
      -logic_name => 'apply_threshold',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectGeneOnFilter',
      -parameters => {
        dna_db => $self->o('dna_db'),
        source_db => $self->o('cdna_db'),
        logic_name => 'cdna_alignment',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna_est2genome')},
        exonerate_cdna_pid => $self->o('cdna_selection_pid'),
        exonerate_cdna_cov => $self->o('cdna_selection_cov'),
      },
      -rc_name => 'default',
      -analysis_capacity => 5,
      -batch_size => 10,
      -flow_into => {
        '1' => ['create_cdna2genome_slices'],
      },
    },


    {
      -logic_name => 'create_cdna2genome_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db => $self->o('cdna_db'),
        iid_type => 'feature_region',
        feature_type => 'gene',
        logic_name => ['cdna_alignment'],
        coord_system_name => 'toplevel',
        include_non_reference => 0,
        top_level => 1,
        use_annotation => 1,
# These options will create only slices that have a gene on the slice in one of the feature dbs
        annotation_file => $self->o('annotation_file'),
        region_padding => $self->o('cdna2genome_region_padding'),
      },
      -flow_into => {
        '2->A' => ['cdna2genome'],
        'A->1' => ['internal_stop'],
      },

      -rc_name    => 'default',
    },


    {
      -logic_name => 'cdna2genome',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion',
      -rc_name    => '3GB',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        sequence_table_name => $self->o('cdna_table_name'),
        source_db => $self->o('cdna_db'),
        target_db => $self->o('cdna2genome_db'),
        logic_name => 'cdna2genome',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna2genome')},
        calculate_coverage_and_pid => 1,
        exonerate_path => $self->o('exonerate_annotation'),
        annotation_file => $self->o('annotation_file'),
        SOFT_MASKED_REPEATS => '#wide_repeat_logic_names#',
      },
      -batch_size => 10,
      -flow_into => {
        '-1' => ['cdna2genome_himem'],
      },
    },

    {
      -logic_name => 'cdna2genome_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion',
      -rc_name    => '6GB',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        sequence_table_name => $self->o('cdna_table_name'),
        source_db => $self->o('cdna_db'),
        target_db => $self->o('cdna2genome_db'),
        logic_name => 'cdna2genome',
        module     => 'HiveExonerate2Genes',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna2genome')},
        calculate_coverage_and_pid => 1,
        exonerate_path => $self->o('exonerate_annotation'),
        annotation_file => $self->o('annotation_file'),
        repeat_libraries => '#wide_repeat_logic_names#',
      },
      -batch_size => 10,
    },
    {
      -logic_name => 'internal_stop',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveInternalStopFix',
      -parameters => {
        dna_db => $self->o('dna_db'),
        source_db => $self->o('cdna2genome_db'),
        edited_biotype => 'edited',
        stop_codon_biotype => 'stop_codon',
        logic_name => 'cdna2genome',
        biotype => undef,
        source => undef,
      },
      -rc_name    => 'default',
    },
    {

      -logic_name => 'download_selenocysteines',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
      -parameters => {
        taxon_id => $self->o('taxon_id'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', 'selenocysteine')},
        output_path => $self->o('targetted_path'),
      },
      -rc_name          => 'default',
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
      -rc_name          => 'default',
      -flow_into => {
        2 => ['process_selenocysteine'],
      },
    },
    {

      -logic_name => 'process_selenocysteine',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelenocysteineFinder',
      -parameters => {
        target_db => $self->o('genewise_db'),
        dna_db => $self->o('dna_db'),
        genome => $self->o('softmasked_genome_file'),
        biotype => 'seleno_self',
        exonerate => $self->o('exonerate_path'),
        genewise => $self->o('genewise_path'),
        iid_type => 'db_seq',
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name          => 'default',
    },
    {
      -logic_name => 'create_besttargetted_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('best_targeted_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['generate_besttargetted_index'],
      },
    },
    {
      -logic_name => 'generate_besttargetted_index',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateBestTargettedIndex',
      -parameters => {
        source_db        => $self->o('cdna2genome_db'),
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        fasta_filename => catfile($self->o('targetted_path'), 'best_targetted.fa'),
        email => $self->o('email_address'),
        genbank_file => $self->o('cdna_file'),
        protein_files => [catfile($self->o('targetted_path'), 'proteome.fa')],
      },
      -rc_name      => 'default',
      -flow_into => {
        1 => ['indicate_BT'],
      },
    },
    {
      -logic_name => 'indicate_BT',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => '#indicate_path# -d #indicate_dir# -f #proteome# -i #indicate_dir#/best_targetted_index -p singleWordParser -M BTMultiParser',
        indicate_path => $self->o('indicate_path'),
        proteome => 'best_targetted.fa',
        indicate_dir => $self->o('targetted_path'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['generate_besttargetted_jobs'],
      },
    },
    {
      -logic_name => 'generate_besttargetted_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db        => $self->o('genewise_db'),
        coord_system_name => 'toplevel',
        iid_type => 'slice',
        feature_constraint => 1,
        feature_type => 'gene',
        top_level => 1,
        feature_dbs => [$self->o('genewise_db'), $self->o('cdna2genome_db')],
      },
      -rc_name      => 'default',
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
        dna_db => $self->o('dna_db'),
        source_db => { protein_db => $self->o('genewise_db'), cdna2genome_db => $self->o('cdna2genome_db')},
        SEQFETCHER_DIR => [catfile($self->o('targetted_path'), 'proteome_index'),
        catfile($self->o('targetted_path'), 'best_targetted_index')],
        SEQFETCHER_OBJECT => 'Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher',
        INPUT_DATA_FROM_DBS  => {
          protein_db => ['seleno_self', 'gw_gtag', 'gw_nogtag', 'gw_exo', 'targetted_exonerate'],
          cdna2genome_db => ['cdna2genome', 'edited'],
        },
        BIOTYPES => ['seleno_self', 'cdna2genome', 'edited', 'gw_gtag', 'gw_nogtag', 'gw_exo', 'targetted_exonerate'], # sorted list, preferred is first
        protein_min_coverage => $self->o('best_targetted_min_coverage'),
        protein_min_identity => $self->o('best_targetted_min_identity'),
      },
      -rc_name          => '3GB',
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
      -rc_name    => 'default',
    },


######################################################################################
#
# IG and TR genes
#
######################################################################################

      {
        -logic_name => 'create_ig_tr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                     source_db => $self->o('dna_db'),
                     target_db => $self->o('ig_tr_db'),
                     create_type => 'clone',
                   },
         -flow_into => {
                         1 => ['load_ig_tr_seqs'],
                       },
      },


      {
        -logic_name => 'load_ig_tr_seqs',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('load_fasta_script_path')
                                ." -dbhost ".$self->o('pipeline_db','-host')
                                ." -dbuser ".$self->o('pipeline_db','-user')
                                ." -dbpass ".$self->o('pipeline_db','-pass')
                                ." -dbname ".$self->o('pipeline_db','-dbname')
                                ." -dbport ".$self->o('pipeline_db','-port')
                                ." -fasta_file ".$self->o('ig_tr_blast_path')."/".$self->o('ig_tr_fasta_file')
                                ." -sequence_table_name ".$self->o('ig_tr_table_name')
                                ." -create_table 1"
                                ." -force_uniq_hitnames 1",
                       },

         -rc_name => 'default',
         -flow_into => {
                         1 => ['generate_ig_tr_jobs'],
                       },
      },


      {
        -logic_name => 'generate_ig_tr_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'sequence_accession',
                         batch_size => $self->o('ig_tr_batch_size'),
                         sequence_table_name => $self->o('ig_tr_table_name'),
                       },
        -rc_name      => 'default',
        -flow_into => {
                        '2->A' => ['ig_tr_genblast'],
                        'A->1' => ['update_ig_tr_hitnames'],
                      },
      },


      {
        -logic_name => 'ig_tr_genblast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('ig_tr_db'),
                         logic_name => 'ig_tr_gene',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('softmasked_genome_file'),
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('ig_tr_table_name'),
                         max_rank => $self->o('ig_tr_genblast_max_rank'),
                         genblast_pid => $self->o('ig_tr_genblast_pid'),
                         timer => '2h',
                         blast_eval => $self->o('ig_tr_genblast_eval'),
                         blast_cov  => $self->o('ig_tr_genblast_cov'),
                       },
        -rc_name    => 'genblast',
        -flow_into => {
                        -1 => ['split_ig_tr_genblast_jobs'],
                        -2 => ['split_ig_tr_genblast_jobs'],
                        -3 => ['split_ig_tr_genblast_jobs'],
                      },
      },


      {
        -logic_name => 'split_ig_tr_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk',
                         batch_size => 1,
                       },
        -rc_name      => 'default',
        -can_be_empty  => 1,
        -flow_into => {
                        2 => ['ig_tr_genblast_retry'],
                      },
      },


      {
        -logic_name => 'ig_tr_genblast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('ig_tr_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('softmasked_genome_file'),
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('ig_tr_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '1h',
                         blast_eval => $self->o('ig_tr_genblast_eval'),
                         blast_cov  => $self->o('ig_tr_genblast_cov'),
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
        -flow_into => {
                        -1 => ['failed_ig_tr_genblast_proteins'],
                        -2 => ['failed_ig_tr_genblast_proteins'],
                        -3 => ['failed_ig_tr_genblast_proteins'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'failed_ig_tr_genblast_proteins',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },


      {
        -logic_name => 'update_ig_tr_hitnames',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          sql => 'UPDATE protein_align_feature set hit_name = concat("IMGT_", hit_name) where hit_name not like "%ENS%";',
          db_conn => $self->o('ig_tr_db'),
        },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['ig_tr_sanity_checks'],
                      },
      },


      {
        -logic_name => 'ig_tr_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('ig_tr_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'ig_tr'},
                       },

        -rc_name    => '2GB',

        -flow_into => {
                        '1->A' => ['fan_ncrna'],
                        'A->1' => ['create_rnaseq_for_layer_db'],
                      },
      },


############################################################################
#
# ncRNA pipeline
#
############################################################################
      {
        -logic_name => 'fan_ncrna',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_ncrna#" -ne "0" ]; then exit 42; else mkdir -p '.$self->o('ncrna_dir').'; exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['create_ncrna_db'],
        },
      },


      {
        -logic_name => 'create_ncrna_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('ncrna_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                          '1->A' => ['dump_genome', 'dump_repeats', 'fetch_rfam_accessions'],
                          'A->1' => ['create_small_rna_slice_ids'],
                      },

      },

      {
        -logic_name => 'dump_repeats',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                    cmd => 'perl '.catfile($self->o('sncrna_analysis_script'), 'repeats_dump.pl').' '.
                            $self->o('dna_db_name'). " " .
                            $self->o('dna_db_server') . " " .
                            $self->o('dna_db_port') . " " .
                            $self->o('user_r') . " " .
                            $self->o('ncrna_dir').' blastmirna',
                      },
       -rc_name => '6GB',
      },


      {
        -logic_name => 'dump_genome',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                          cmd => 'perl '.$self->o('sequence_dump_script').' -dbhost '.$self->o('dna_db_server')
                                  .' -dbname  '.$self->o('dna_db_name').' -dbport '.$self->o('dna_db_port')
                                  .' -dbuser '.$self->o('user_r')
                                  .' -coord_system_name toplevel -mask -mask_repeat '.$self->o('full_repbase_logic_name')
                                  .' -output_dir '.$self->o('genome_dumps')
                                  .' -softmask -onefile -header rnaseq -filename '.$self->o('faidx_genome_file'),
                      },
         -rc_name   => 'filter',
      },

      {
        -logic_name => 'fetch_rfam_accessions',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
            cmd => "perl " . $self->o('sncrna_analysis_script') . "/fetch_rfam_accessions.pl " .
                    " -h " . $self->o('rfam_host') . 
                    " -u " . $self->o('rfam_user') .
                    " -p " . $self->o('rfam_port') .
                    " -d " . $self->o('rfam_dbname') .
                    " -c " . $self->o('clade') .
                    " -o " . $self->o('output_path'),
        },
        -rc_name => 'default',
        -flow_into => { '1' => 'extract_rfam_cm'},
      },

      {
        -logic_name => 'extract_rfam_cm',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
            cmd => "perl " . $self->o('sncrna_analysis_script') . "/filter_cm.pl " .
                             $self->o('rfam_cm') . ' ' .
                             $self->o('output_path') . '/accessions.txt ' .
                             $self->o('output_path'),
        },
        -rc_name => 'filter',
      },

      {
        -logic_name => 'create_small_rna_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 1000000,
                         top_level => 1,
                         target_db => $self->o('dna_db'),
                         batch_slice_ids => 1,
                         batch_target_size => 2000000,
                      },
        -flow_into => {
                       '2->A' => ['mirna_blast','run_cmsearch'],
                       'A->1' => ['filter_ncrnas'],
                      },
        -rc_name    => 'default',
      },

      {
        -logic_name => 'run_cmsearch',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCMSearch',
        -parameters => {
                    output_db => $self->o('ncrna_db'),
                    dna_db => $self->o('dna_db'),
                    logic_name => 'ncrna',
                    module     => 'HiveCMSearch',
                    cmsearch_exe_path => $self->o('cmsearch_exe_path'),
                    rfam_seeds => $self->o('rfam_seeds'),
                    rfam_cm => $self->o('filtered_rfam_cm'),
                    blast_db_dir_path => $self->o('ncrna_blast_path').'/',
                    output_dir => $self->o('output_path'),
        },
        -hive_capacity => 900,
        -rc_name    => '5GB',
        -flow_into => { '-1' => 'run_cmsearch_highmem', '-2' => 'run_cmsearch_highmem'},
      },

      {
        -logic_name => 'run_cmsearch_highmem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCMSearch',
        -parameters => {
                    output_db => $self->o('ncrna_db'),
                    dna_db => $self->o('dna_db'),
                    logic_name => 'ncrna',
                    module     => 'HiveCMSearch',
                    cmsearch_exe_path => $self->o('cmsearch_exe_path'),
                    rfam_seeds => $self->o('rfam_seeds'),
                    rfam_cm => $self->o('filtered_rfam_cm'),
                    blast_db_dir_path => $self->o('ncrna_blast_path').'/',
                    output_dir => $self->o('output_path'),
        },
        -hive_capacity => 900,
        -rc_name    => '10GB',
      },


      {
        -logic_name => 'mirna_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastmiRNA',
        -parameters => {
                         repeat_logic_names => ['dust'],
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'blastmirna',
                         module     => 'HiveBlastmiRNA',
                         blast_db_path => catfile($self->o('mirna_blast_path'), $self->o('mirBase_fasta')),
                         blast_exe_path => $self->o('blastn_exe_path'),
                         commandline_params => ' -num_threads 3 ',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastmiRBase', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '2h',
                       },

        -flow_into => {
                        '-1' => ['rebatch_mirna'],
                        '-2' => ['rebatch_mirna'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'blast',
      },


      {
        -logic_name => 'rebatch_mirna',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 100000,
                       },
       -flow_into => {
                       '2' => ['mirna_blast_retry'],
                     },
       -rc_name    => 'default',
       -can_be_empty  => 1,
      },


      {
        -logic_name => 'mirna_blast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastmiRNA',
        -parameters => {
                         repeat_logic_names => ['dust'],
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'blastmirna',
                         module     => 'HiveBlastmiRNA',
                         blast_db_path => $self->o('mirna_blast_path') . '/' . $self->o('mirBase_fasta') ,
                         blast_exe_path => $self->o('blastn_exe_path'),
                         commandline_params => ' -num_threads 3 ',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastmiRBase', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '1h',
                       },

        -flow_into => {
                        '-1' => ['failed_mirna_blast_job'],
                        '-2' => ['failed_mirna_blast_job'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'blast_retry',
      },


      {
        -logic_name => 'failed_mirna_blast_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'filter_ncrnas',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFilterncRNAs',
        -parameters => {
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'filter_ncrnas',
                         module     => 'HiveFilterncRNAs',
                       },
        -rc_name    => 'filter',
        -flow_into => {
                        '2->A' => ['run_mirna'],
                        'A->1' => ['dump_features'],
                        #'3' => ['run_infernal'],
                      },
      },


      {
        -logic_name => 'run_mirna',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivemiRNA',
        -parameters => {
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'ncrna',
                         module     => 'HivemiRNA',
                         blast_db_dir_path => catfile($self->o('mirna_blast_path'), 'all_mirnas.embl'),
                         output_dir => $self->o('ncrna_dir'),
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'filter',
        #-flow_into  => { '1->A' => ['dump_repeats', 'dump_features', 'dump_genome'], 'A->1' => ['filter_mirnas']},
      },


      {
        -logic_name => 'dump_features',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.catfile($self->o('sncrna_analysis_script'), 'dump_prefilter_features.pl').' '.
                                $self->o('ncrna_db_name'). " " .
                                $self->o('ncrna_db_server') . " " .
                                $self->o('ncrna_db_port') . " " .
                                $self->o('user_r') . " " .
                                $self->o('ncrna_dir').' blastmirna',
                        },
         -rc_name   => 'filter',
         -flow_into => {'1' => 'filter_mirnas'},
      },


      {
        -logic_name => 'filter_mirnas',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                          cmd => 'sh '.catfile($self->o('sncrna_analysis_script'), 'FilterMiRNAs.sh')
                                  .' -d '.catfile($self->o('ncrna_dir'), 'blastmirna_dafs.bed')
                                  .' -r '.catfile($self->o('ncrna_dir'), 'repeats.bed')
                                  .' -g '.$self->o('faidx_genome_file')
                                  .' -w '.$self->o('ncrna_dir')
                                  .' -m '.catfile($self->o('mirna_blast_path'), 'rfc_filters', $self->o('rfc_model'))
                                  .' -s '.catfile($self->o('mirna_blast_path'), 'rfc_filters', $self->o('rfc_scaler'))
                                  .' -c '.$self->o('ncrna_db_server').":".$self->o('ncrna_db_port').":".$self->o('ncrna_db_name').":"
                                  .$self->o('user').":".$self->o('password'),
                        },
        -rc_name   => 'filter',
        -flow_into => {'1' => 'ncrna_sanity_checks' },

      },


      {
        -logic_name => 'ncrna_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('ncrna_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'ncrna'},
                       },

        -rc_name    => '4GB',
#        -flow_into => {
#          1 => ['transfer_ncrnas'],
#        },
      },


############################################################################
#
# RNA-seq analyses
#
############################################################################
      {
        -logic_name => 'create_rnaseq_for_layer_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('dna_db'),
          target_db => $self->o('rnaseq_for_layer_db'),
          create_type => 'clone',
        },
        -rc_name    => 'default',
        -flow_into => {
          '1->A' => ['fan_rnaseq_for_layer_db'],
          'A->1' => ['create_long_read_final_db'],
        },
      },


      {
        -logic_name => 'fan_rnaseq_for_layer_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_rnaseq#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['checking_file_path'],
        },
      },


      {
        -logic_name => 'checking_file_path',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '1GB',
        -parameters => {
          cmd => 'EXIT_CODE=0; for F in '.join (' ',
              $self->o('bwa_path'),
              $self->o('samtools_path'),
              $self->o('exonerate_path'),
              $self->o('sequence_dump_script'),
              $self->o('uniprot_blast_exe_path')
              ).'; do which "$F"; if [ "$?" == 1 ]; then EXIT_CODE=1;fi; done; '
            .'if [ $EXIT_CODE -eq 1 ];then exit $EXIT_CODE;fi; '
            .'for D in '.join(' ',
              $self->o('output_dir'),
              $self->o('input_dir'),
              $self->o('merge_dir'),
              $self->o('sam_dir')
              ).'; do mkdir -p "$D"; done; '
            .'which lfs > /dev/null; if [ $? -eq 0 ]; then for D in '.join(' ',
              $self->o('output_dir'),
              $self->o('input_dir'),
              $self->o('merge_dir')
              ).'; do lfs getdirstripe -q $D > /dev/null; if [ $? -eq 0 ]; then lfs setstripe -c -1 $D;fi;done;fi',
        },
        -flow_into => {
          '1->A' => ['create_fastq_download_jobs','index_rnaseq_genome_file'],
          'A->1' => ['create_rough_db'],
        },
      },


      {
        -logic_name => 'create_fastq_download_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputfile => $self->o('rnaseq_summary_file'),
          column_names => $self->o('file_columns'),
          delimiter => '\t',
        },
        -flow_into => {
          2 => {'download_RNASeq_fastqs' => {'iid' => '#filename#'}},
        },
      },


      {
        -logic_name => 'download_RNASeq_fastqs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs',
        -parameters =>{
          ftp_base_url => $self->o('rnaseq_ftp_base'),
          input_dir => $self->o('input_dir'),
        },
        -flow_into => {
          1 => ['get_read_lengths'],
        },
        -analysis_capacity => 50,
      },


      {
        -logic_name => 'get_read_lengths',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCalculateReadLength',
        -parameters =>{
          input_dir => $self->o('input_dir'),
          read_length_table => $self->o('read_length_table'),
        },
        -flow_into => {
          1 => ['split_fastq_files'],
        },
      },


      {
        -logic_name => 'split_fastq_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SplitFastQFiles',
        -parameters => {
          'max_reads_per_split' => $self->o('max_reads_per_split'),
          'max_total_reads'     => $self->o('max_total_reads'),
          'rnaseq_summary_file' => $self->o('rnaseq_summary_file'),
          'fastq_dir'           => $self->o('input_dir'),
        },
      },


     {
        -logic_name => 'index_rnaseq_genome_file',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '5GB',
        -parameters => {
          cmd => 'if [ ! -e "'.$self->o('faidx_genome_file').'.ann" ]; then '.$self->o('bwa_path').' index -a bwtsw '.$self->o('faidx_genome_file').';fi',
        },
      },


      {
        -logic_name => 'create_rough_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('dna_db'),
          target_db => $self->o('rnaseq_rough_db'),
          create_type => 'clone',
        },
        -rc_name => '1GB',
        -flow_into => {
          1 => ['backup_original_csv'],
        },
      },


      {
        -logic_name => 'backup_original_csv',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cp '.$self->o('rnaseq_summary_file').' '.$self->o('rnaseq_summary_file').'_orig_bak',
        },
        -flow_into => {
          1 => ['create_updated_csv'],
        },
      },


      {
        -logic_name => 'create_updated_csv',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cat '.catfile($self->o('input_dir'), '*_new.csv').' > '.$self->o('rnaseq_summary_file'),
        },
        -flow_into => {
          1 => ['parse_summary_file'],
        },
      },

      {
        -logic_name => 'parse_summary_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable',
        -rc_name => '1GB',
        -parameters => {
          column_names => $self->o('file_columns'),
          sample_column => $self->o('read_group_tag'),
          inputfile => $self->o('rnaseq_summary_file'),
          delimiter => $self->o('summary_file_delimiter'),
          csvfile_table => $self->o('summary_csv_table'),
          pairing_regex => $self->o('pairing_regex'),
          read_length_table => $self->o('read_length_table'),
        },
        -flow_into => {
          '2->A' => [ 'create_tissue_jobs'],
          'A->1' => [ 'merged_bam_file' ],
        },
      },
      {
        -logic_name => 'create_tissue_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputquery => join(' ', 'SELECT', $self->o('read_group_tag'), ',', $self->o('read_id_tag'), ', is_paired, CN', 'FROM', $self->o('summary_csv_table'), 'WHERE', $self->o('read_group_tag'), '= "#sample_name#"'),
          column_names => [$self->o('read_group_tag'), $self->o('read_id_tag'), 'is_paired', 'rnaseq_data_provider'],
        },
        -rc_name    => '1GB',
        -flow_into => {
          '2->A' => ['create_bwa_jobs'],
          'A->1' => ['merged_tissue_file'],
        },
      },
      {
        -logic_name => 'create_bwa_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateBWAJobs',
        -parameters => {
          sample_column => $self->o('read_group_tag'),
          sample_id_column => $self->o('read_id_tag'),
          csvfile_table => $self->o('summary_csv_table'),
          column_names => $self->o('file_columns'),
          use_threading => $self->o('use_threads'),
        },
        -rc_name    => '1GB',
        -flow_into => {
          '2->A' => ['bwa', 'create_header_files'],
          'A->1' => ['bwa2bam'],
        },
      },
      {
        -logic_name => 'create_header_files',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '1GB',
        -parameters => {
          cmd => 'if [ ! -e "'.$self->o('header_file').'" ]; then printf "'.$header_line.'" > '.$self->o('header_file').'; fi',
        },
      },
      {
        -logic_name => 'bwa',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -parameters => {
          disconnect_jobs => 1,
          short_read_aligner => $self->o('bwa_path'),
          input_dir => $self->o('input_dir'),
          genome_file => $self->o('faidx_genome_file'),
          output_dir => $self->o('output_dir'),
        },
        -flow_into => {
          1 => [ ':////accu?fastq=[]' ],
          -1 => [ 'bwa_20GB' ],
          -2 => [ 'bwa_20GB' ],
        },
        -rc_name    => '10GB_multithread',
      },
      {
        -logic_name => 'bwa_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -can_be_empty => 1,
        -parameters => {
          disconnect_jobs => 1,
          short_read_aligner => $self->o('bwa_path'),
          input_dir => $self->o('input_dir'),
          genome_file => $self->o('faidx_genome_file'),
          output_dir => $self->o('output_dir'),
        },
        -flow_into => {
          1 => [ ':////accu?fastq=[]' ],
        },
        -rc_name    => '20GB_multithread',
      },
      {
        -logic_name => 'bwa2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -parameters => {
          sampe_options => '-A -a '.$self->o('maxintron'),
          samse_options => '',
          min_paired => $self->o('read_min_paired'),
          min_mapped => $self->o('read_min_mapped'),
          header_file => $self->o('header_file'),
          bam_prefix => $self->o('read_id_tag'),
          email => $self->o('email_address'),
          disconnect_jobs => 1,
          short_read_aligner => $self->o('bwa_path'),
          input_dir => $self->o('input_dir'),
          genome_file => $self->o('faidx_genome_file'),
          output_dir => $self->o('output_dir'),
          samtools => $self->o('samtools_path'),
        },
        -flow_into => {
          1 => [ ':////accu?filename=[]' ],
          -1 => {'bwa2bam_20GB' => { fastq => '#fastq#', is_paired => '#is_paired#', $self->o('read_id_tag') => '#'.$self->o('read_id_tag').'#'}},
          -2 => {'bwa2bam_20GB' => { fastq => '#fastq#', is_paired => '#is_paired#', $self->o('read_id_tag') => '#'.$self->o('read_id_tag').'#'}},
        },
        -rc_name    => '10GB',
      },
      {
        -logic_name => 'bwa2bam_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -can_be_empty => 1,
        -parameters => {
          sampe_options => '-A -a '.$self->o('maxintron'),
          samse_options => '',
          min_paired => $self->o('read_min_paired'),
          min_mapped => $self->o('read_min_mapped'),
          header_file => $self->o('header_file'),
          bam_prefix => $self->o('read_id_tag'),
          email => $self->o('email_address'),
          disconnect_jobs => 1,
          short_read_aligner => $self->o('bwa_path'),
          input_dir => $self->o('input_dir'),
          genome_file => $self->o('faidx_genome_file'),
          output_dir => $self->o('output_dir'),
          samtools => $self->o('samtools_path'),
        },
        -flow_into => {
          1 => [ ':////accu?filename=[]' ],
        },
        -rc_name    => '20GB',
      },

      {
        -logic_name => 'merged_tissue_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
        -parameters => {
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BamMergeStatic', $self->o('rnaseq_merge_type'))},
          # target_db is the database where we will write the files in the data_file table
          # You can use store_datafile => 0, if you don't want to store the output file
          target_db => $self->o('rnaseq_rough_db'),
          assembly_name => $self->o('assembly_name'),
          rnaseq_data_provider => $self->o('rnaseq_data_provider'),
          disconnect_jobs => 1,
          species => $self->o('species_name'),
          output_dir => $self->o('merge_dir'),
          input_dir => $self->o('output_dir'),
          samtools => $self->o('samtools_path'),
          picard_lib_jar => $self->o('picard_lib_jar'),
          use_threads => $self->o('rnaseq_merge_threads'),
        },
        -rc_name    => '3GB_merged_multithread',
        -flow_into => {
          1 => ['create_analyses_type_job', '?accu_name=filename&accu_address=[]&accu_input_variable=alignment_bam_file'],
        },
      },
      {
        -logic_name => 'create_analyses_type_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -rc_name    => '1GB',
        -parameters => {
          inputlist => ['gene', 'daf', 'ise'],
          column_names => ['type'],
          species => $self->o('species_name'),
        },
        -flow_into => {
          2 => {'create_rnaseq_tissue_analyses' => {analyses => [{'-logic_name' => '#species#_#sample_name#_rnaseq_#type#'}]}},
        },
      },
      {
        -logic_name => 'create_rnaseq_tissue_analyses',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
        -rc_name    => '1GB',
        -parameters => {
          source_type => 'list',
          target_db => $self->o('rnaseq_rough_db'),
        },
      },
      {
        -logic_name => 'merged_bam_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
        -parameters => {
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BamMergeStatic', $self->o('rnaseq_merge_type'))},
          # target_db is the database where we will write the files in the data_file table
          # You can use store_datafile => 0, if you don't want to store the output file
          target_db => $self->o('rnaseq_rough_db'),
          assembly_name => $self->o('assembly_name'),
          rnaseq_data_provider => $self->o('rnaseq_data_provider'),
          disconnect_jobs => 1,
          alignment_bam_file => catfile($self->o('merge_dir'), '#assembly_name#.#rnaseq_data_provider#.merged.1.bam'),
          species => $self->o('species_name'),
          output_dir => $self->o('merge_dir'),
          input_dir => $self->o('merge_dir'),
          samtools => $self->o('samtools_path'),
          picard_lib_jar => $self->o('picard_lib_jar'),
          use_threads => $self->o('rnaseq_merge_threads'),
        },
        -rc_name    => '5GB_merged_multithread',
        -flow_into => {
          1 => ['fan_merge_analyses'],
          2 => ['create_header_intron'],
        },
      },

     {
        -logic_name => 'fan_merge_analyses',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
	  cmd => "if [[ \$(cut -d\$'\\t' -f1 ".$self->o('rnaseq_summary_file')." | sort | uniq | wc -l) == 1 ]]; then exit 42; else exit 0;fi",
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => {
	  1 => ['create_merge_analyses_type_job'],
        },
      },

      {
        -logic_name => 'create_merge_analyses_type_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -rc_name    => '1GB',
        -parameters => {
          inputlist => ['gene', 'daf', 'ise'],
          column_names => ['type'],
          species => $self->o('species_name'),
        },
        -flow_into => {
          2 => {'create_rnaseq_merge_analyses' => {analyses => [{'-logic_name' => '#species#_merged_rnaseq_#type#'}]}},
        },
      },
      {
        -logic_name => 'create_rnaseq_merge_analyses',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
        -rc_name    => '1GB',
        -parameters => {
          source_type => 'list',
          target_db => $self->o('rnaseq_rough_db'),
        },
      },
      {
        -logic_name => 'create_header_intron',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name    => '1GB',
        -parameters => {
          cmd => $self->o('samtools_path').' view -H #filename# | grep -v @SQ | grep -v @HD > '.catfile($self->o('rnaseq_dir'),'merged_header.h'),
        },
        -flow_into => {
          '1->A' => [ 'create_toplevel_input_ids'],
          'A->1' => ['sam2bam'],
        },
      },


      {
        -logic_name => 'create_toplevel_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'slice',
          batch_slice_ids => 1,
          batch_target_size => 1000000,
          target_db => $self->o('rnaseq_rough_db'),
        },
        -flow_into => {
          2 => {'create_overlapping_slices' => {'iid' => '#iid#', alignment_bam_file => '#filename#'}},
        },
      },

      {
        -logic_name => 'create_overlapping_slices', # Hopefully this can be removed when the new Bam2Genes module is ready
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'split_slice',
          slice_size => 5000000,
          slice_overlaps => 2500000,
          target_db => $self->o('rnaseq_rough_db'),
        },
        -flow_into => {
          '2->A' => {'split_on_low_coverage' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          'A->1' => ['check_and_delete_broken_duplicated'],
        },
      },

      {
        -logic_name => 'split_on_low_coverage',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::CalculateLowCoverageSlices',
        -parameters => {
          'disconnect_jobs' => 1,
          'dna_db'   => $self->o('dna_db'),
        },

        -rc_name   => '10GB',
        -flow_into => {
          2 => ['rough_transcripts'],
          -1 => ['split_on_low_coverage_20GB'],
        },

      },

      {
        -logic_name => 'split_on_low_coverage_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::CalculateLowCoverageSlices',
        -parameters => {
          'disconnect_jobs' => 1,
          'dna_db'   => $self->o('dna_db'),
        },

        -rc_name          => '20GB',
        -flow_into => {
          2 => ['rough_transcripts'],
        },

      },


      {
        -logic_name => 'rough_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -parameters => {
          logic_name => 'rough_transcripts',
          target_db    => $self->o('rnaseq_rough_db'),
          dna_db    => $self->o('dna_db'),
          min_length => 300,
          min_exons  =>   1,
          max_intron_length => $self->o('maxintron'),
          min_single_exon_length => 1000,
          min_span   =>   1.5,
          paired => $self->o('paired'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
        },
        -rc_name    => '2GB',
        -flow_into => {
          -1 => {'rough_transcripts_10GB' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -2 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -3 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
        },
      },


      {
        -logic_name => 'rough_transcripts_10GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -can_be_empty => 1,
        -parameters => {
          logic_name => 'rough_transcripts',
          target_db    => $self->o('rnaseq_rough_db'),
          dna_db    => $self->o('dna_db'),
          min_length => 300,
          min_exons  =>   1,
          max_intron_length => $self->o('maxintron'),
          min_single_exon_length => 1000,
          min_span   =>   1.5,
          paired => $self->o('paired'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
        },
        -rc_name    => '10GB',
        -flow_into => {
          -1 => {'rough_transcripts_30GB' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -2 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -3 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
        },
      },


      {
        -logic_name => 'rough_transcripts_30GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -can_be_empty => 1,
        -parameters => {
          logic_name => 'rough_transcripts',
          target_db    => $self->o('rnaseq_rough_db'),
          dna_db    => $self->o('dna_db'),
          min_length => 300,
          min_exons  =>   1,
          max_intron_length => $self->o('maxintron'),
          min_single_exon_length => 1000,
          min_span   =>   1.5,
          paired => $self->o('paired'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
        },
        -rc_name    => '30GB',
        -flow_into => {
          -1 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -2 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -3 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
        },
      },


      {
        -logic_name => 'failed_rough_transcript_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'check_and_delete_broken_duplicated',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects',
        -parameters => {
          target_db => $self->o('rnaseq_rough_db'),
          check_support => 0,
        },
        -rc_name    => '5GB',
        -analysis_capacity => 1, # Because there is slice overlap, having parallel jobs can cause problems
        -flow_into => {
          1 => ['create_bam2introns_input_ids'],
        },
      },

      {
        -logic_name => 'create_bam2introns_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          iid_type => 'slice_to_feature_ids',
          target_db => $self->o('rnaseq_rough_db'),
          feature_type => 'gene',
          logic_name => ['rough_transcripts'],
          use_stable_ids => 1,
          create_stable_ids => 1,
          stable_id_prefix => 'RNASEQ',
        },
        -rc_name    => '1GB',
        -batch_size => 100,
        -flow_into => {
          2 => {'bam2introns' => {iid => '#iid#', bam_file => '#alignment_bam_file#'}},
        },
      },

      {
        -logic_name => 'bam2introns',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -parameters => {
          program_file => $self->o('exonerate_path'),
          source_db => $self->o('rnaseq_rough_db'),
          dna_db => $self->o('dna_db'),
          missmatch => 6,
          word_length => 10,
          saturate_threshold => 10000,
          mask => 1,
          percent_id => 97,
          coverage => 90,
          fullseq   => 1,
          max_transcript => 1000000,
          batch_size => 10000,
          maxintron => $self->o('maxintron'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
          output_dir => $self->o('sam_dir'),
          use_genome_flatfile => 1,
          flatfile_masked => 1,
          genome_file => $self->o('faidx_softmasked_genome_file'),
          timer => '15m',
        },
        -rc_name    => '5GB',
        -batch_size => 200,
        -analysis_capacity => 500,
        -flow_into => {
          1 => [':////accu?filename=[]'],
          2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -1 => {'bam2introns_20GB' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -2 => ['failed_bam2introns_jobs'],
          -3 => ['failed_bam2introns_jobs'],
        },
      },


      {
        -logic_name => 'bam2introns_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -can_be_empty => 1,
        -parameters => {
          program_file => $self->o('exonerate_path'),
          source_db => $self->o('rnaseq_rough_db'),
          dna_db => $self->o('dna_db'),
          missmatch => 6,
          word_length => 10,
          saturate_threshold => 10000,
          mask => 1,
          percent_id => 97,
          coverage => 90,
          fullseq   => 1,
          max_transcript => 1000000,
          batch_size => 10000,
          maxintron => $self->o('maxintron'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
          output_dir => $self->o('sam_dir'),
          use_genome_flatfile => 1,
          flatfile_masked => 1,
          genome_file => $self->o('faidx_softmasked_genome_file'),
          timer => '30m',
        },
        -rc_name    => '20GB',
        -flow_into => {
          1 => [':////accu?filename=[]'],
          2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -1 => {'bam2introns_30GB' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -2 => ['failed_bam2introns_jobs'],
          -3 => ['failed_bam2introns_jobs'],
        },
      },

      {
        -logic_name => 'bam2introns_30GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -can_be_empty => 1,
        -parameters => {
          program_file => $self->o('exonerate_path'),
          source_db => $self->o('rnaseq_rough_db'),
          dna_db => $self->o('dna_db'),
          missmatch => 6,
          word_length => 10,
          saturate_threshold => 10000,
          mask => 1,
          percent_id => 97,
          coverage => 90,
          fullseq   => 1,
          max_transcript => 1000000,
          batch_size => 10000,
          maxintron => $self->o('maxintron'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
          output_dir => $self->o('sam_dir'),
          use_genome_flatfile => 1,
          flatfile_masked => 1,
          genome_file => $self->o('faidx_softmasked_genome_file'),
          timer => '30m',
        },
        -rc_name    => '30GB',
        -flow_into => {
          1 => [':////accu?filename=[]'],
          2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -1 => ['failed_bam2introns_jobs'],
          -2 => ['failed_bam2introns_jobs'],
          -3 => ['failed_bam2introns_jobs'],
        },
      },

      {
        -logic_name => 'failed_bam2introns_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'sam2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSam2Bam',
        -parameters => {
          headerfile => catfile($self->o('rnaseq_dir'), 'merged_header.h'),
          disconnect_jobs => 1,
          samtools => $self->o('samtools_path'),
          intron_bam_file => catfile($self->o('output_dir'), 'introns'),
          genome_file => $self->o('faidx_genome_file'),
          use_threading => $self->o('use_threads'),
          sam_dir => $self->o('sam_dir'),
        },
        -rc_name    => '5GB_multithread',
        -flow_into => ['create_refine_db'],
      },


     {
        -logic_name => 'create_refine_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('rnaseq_rough_db'),
          target_db => $self->o('rnaseq_refine_db'),
          create_type => 'clone',
          extra_data_tables => ['data_file'],
        },
        -rc_name => '1GB',
        -flow_into => ['create_blast_db'],
      },

      {
        -logic_name => 'create_blast_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('rnaseq_refine_db'),
          target_db => $self->o('rnaseq_blast_db'),
          create_type => 'clone',
          extra_data_tables => ['data_file'],
        },
        -rc_name => '1GB',
        -flow_into => ['create_ccode_config'],
      },
      {
        -logic_name => 'create_ccode_config',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRefineConfig',
        -parameters => {
          species => $self->o('species_name'),
          output_dir => $self->o('output_dir'),
          intron_bam_file => catfile($self->o('output_dir'), 'introns'),
          single_tissue => $self->o('single_tissue'),
          sample_column => $self->o('read_group_tag'),
          sample_id_column => $self->o('read_id_tag'),
          csvfile_table => $self->o('summary_csv_table'),
          source_db => $self->o('rnaseq_rough_db'),
          dna_db => $self->o('dna_db'),
          target_db => $self->o('rnaseq_refine_db'),
          # write the intron features into the OUTPUT_DB along with the models
          write_introns => 1,
          # maximum number of times to loop when building all possible paths through the transcript
          max_recursions => 10000000000000,
          # analysis logic_name for the dna_align_features to fetch from the INTRON_DB
          # If left blank all features will be fetched
          logicname => [],
          # logic name of the gene models to fetch
          model_ln  => '',
          # penalty for removing a retined intron
          retained_intron_penalty => 2,
          #Remove introns that overlap X introns
          filter_on_overlap => 0,
          # minimum size for an intron
          min_intron_size  => 30,
          max_intron_size  => $self->o('maxintron'),
          # biotype to give to single exon models if left blank single exons are ignored
          # minimum single exon size (bp)
          min_single_exon => 1000,
          # minimum percentage of single exon length that is coding
          single_exon_cds => 66,
          # Intron with most support determines the splice sites for an internal exon
          # lower scoring introns with different splice sites are rejected
          strict_internal_splice_sites => 1,
          # In some species alternate splice sites for end exons seem to be common
          strict_internal_end_exon_splice_sites => 1,
          # biotypes to give gene models if left blank these models will not get written to the output database
          # best score - model with most supporting intron features
          # all other possible models
          # max number of other models to make - blank = all
          other_num      => '10',
          # max number of other models to process - blank = all
          max_num      => '1000',
          other_isoforms => $self->o('other_isoforms'),
          # biotype to label bad models ( otherwise they are not written )
          # do you want to trim UTR
          trim_utr => 1,
          # config for trimming UTR
          max_3prime_exons => 2,
          max_3prime_length => 5000,
          max_5prime_exons => 3,
          max_5prime_length => 1000,
          # % of average intron score that a UTR intron must have
          reject_intron_cutoff => 5,
        },

        -rc_name          => '1GB',
        -flow_into => {
          '2->A' => ['create_ccode_input_ids'],
          'A->1' => ['copy_rnaseq_blast_db'],
        },
      },
      {
        -logic_name => 'create_ccode_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'slice',
          coord_system_name => 'toplevel',
          slice => 1,
          include_non_reference => 0,
          top_level => 1,
          feature_constraint => 1,
          feature_type => 'gene',
          target_db => $self->o('rnaseq_rough_db'),
        },
        -flow_into => {
          2 => {'refine_genes' => {iid => '#iid#', logic_name => '#logic_name#', config_file => '#config_file#'}},
        },
      },

      {
        -logic_name => 'refine_genes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => $self->o('refine_ccode_exe').($self->o('use_ucsc_naming') ? ' -u ' : ' ').($self->o('use_threads') ? ' -t '.$self->o('use_threads').' ' : ' ').'-c #config_file# -i #iid# -l #logic_name# -v 0',
          return_codes_2_branches => {
            42 => 2,
          },
        },
        -rc_name => '2GB_multithread',
        -flow_into => {
          1 => ['create_gene_id_input_ids'],
          -1 => {'refine_genes_20GB' => {iid => '#iid#', config_file => '#config_file#', logic_name => '#logic_name#'}},
        },
      },
      {
        -logic_name => 'refine_genes_20GB',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => $self->o('refine_ccode_exe').($self->o('use_ucsc_naming') ? ' -u ' : ' ').($self->o('use_threads') ? ' -t '.$self->o('use_threads').' ' : ' ').'-c #config_file# -i #iid# -l #logic_name# -v 0',
          return_codes_2_branches => {
            42 => 2,
          },
        },
        -rc_name => '20GB_multithread',
        -flow_into => {
          1 => ['create_gene_id_input_ids'],
        },
      },


      {
        -logic_name => 'create_gene_id_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'feature_id',
          coord_system_name => 'toplevel',
          target_db => $self->o('rnaseq_refine_db'),
          feature_logic_names => ['#logic_name#'],
          feature_type => 'gene',
          batch_size => 50,
        },
        -flow_into => {
          2 => {'blast_rnaseq' => {iid => '#iid#', logic_name => '#logic_name#'}},
        },
      },


      {
        -logic_name => 'blast_rnaseq',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
        -parameters => {
          source_db => $self->o('rnaseq_refine_db'),
          target_db => $self->o('rnaseq_blast_db'),
          dna_db => $self->o('dna_db'),
          iid_type => 'object_id',
          # path to index to fetch the sequence of the blast hit to calculate % coverage
          indicate_index => $self->o('protein_blast_index'),
          uniprot_index => [$self->o('protein_blast_db')],
          blast_program => $self->o('uniprot_blast_exe_path'),
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
          commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40',
        },
        -rc_name => 'blast',
      },


      {
        -logic_name => 'copy_rnaseq_blast_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('rnaseq_blast_db'),
                         target_db => $self->o('rnaseq_for_layer_db'),
                         create_type => 'copy',
                         force_drop => 1,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['update_rnaseq_for_layer_biotypes'],
                      },
      },


      {
        -logic_name => 'update_rnaseq_for_layer_biotypes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('rnaseq_for_layer_db'),
          sql => [
            'UPDATE gene SET biotype = "rnaseq_merged" WHERE biotype IN ("best","single","other_merged")',
            'UPDATE gene SET biotype = "rnaseq_tissue" WHERE biotype != "rnaseq_merged"',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['remove_rnaseq_for_layer_daf_features'],
                      },
      },


      {
        -logic_name => 'remove_rnaseq_for_layer_daf_features',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('rnaseq_for_layer_db'),
          sql => [
            'TRUNCATE dna_align_feature',
            'DELETE transcript_supporting_feature FROM transcript_supporting_feature WHERE feature_type = "dna_align_feature"',
            'DELETE supporting_feature FROM supporting_feature WHERE feature_type = "dna_align_feature"',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['classify_rnaseq_for_layer_models'],
                      },
      },


      {
        -logic_name => 'classify_rnaseq_for_layer_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         update_gene_biotype => 1,
                         classification_type => 'standard',
                         target_db => $self->o('rnaseq_for_layer_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['rnaseq_for_layer_sanity_checks'],
                      },

      },


      {
        -logic_name => 'rnaseq_for_layer_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('rnaseq_for_layer_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'rnaseq_blast'},
                       },

        -flow_into => {
          1 => ['create_rnaseq_layer_nr_db'],
        },

#        -flow_into => {
#          1 => ['create_lincrna_db'],
#        },

        -rc_name    => '4GB',
      },


      {
        -logic_name => 'create_rnaseq_layer_nr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('rnaseq_for_layer_db'),
                         target_db => $self->o('rnaseq_for_layer_nr_db'),
                         create_type => 'copy',
                         force_drop => 1,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['create_rnaseq_layer_nr_slices'],
                      },
      },


      {
        -logic_name => 'create_rnaseq_layer_nr_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 20000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                         batch_slice_ids => 1,
                         batch_target_size => 20000000,
                       },
        -rc_name    => '2GB',
        -flow_into => {
                         '2'    => ['remove_redundant_rnaseq_layer_genes'],
                      },
      },


    {
      -logic_name => 'remove_redundant_rnaseq_layer_genes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db => $self->o('rnaseq_for_layer_nr_db'),
        target_type => 'generic',
      },
      -rc_name          => '5GB',
    },

#      {
#        -logic_name => 'create_lincrna_db',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
#        -parameters => {
#                           source_db => $self->o('rnaseq_for_layer_db'),
#                           target_db => $self->o('lincrna_db'),
#                           create_type => 'clone',
#                       },
#        -rc_name    => 'default',
#        -max_retry_count => 0,
#        -flow_into => {
#          '1' => ['create_pfam_analysis'],
#        },
#      },

#      {
#        -logic_name => 'create_pfam_analysis',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
#        -rc_name    => 'default',
#        -parameters => {
#          source_type => 'list',
#          target_db => $self->o('rnaseq_rough_db'),
#          analyses => $self->o('required_analysis'),
#        },
#        -flow_into => {
#          '1' => ['create_lincrna_toplevel_slices'],
#        },
#      },
#      {
#        -logic_name => 'create_lincrna_toplevel_slices',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#                           target_db => $self->o('dna_db'),
#                           coord_system_name => 'toplevel',
#                           iid_type => 'slice',
#                           include_non_reference => 0,
#                           top_level => 1,
#                           slice_size => 2000000,  # this is for the size of the slice
#                         },
#        -flow_into => {
#          '2->A' => ['Hive_LincRNARemoveDuplicateGenes'],
#          'A->1' => ['delete_duplicate_lincrna_genes'],

#        },
#        -rc_name    => 'default',
#      },

#      {
#        -logic_name => 'Hive_LincRNARemoveDuplicateGenes',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveDuplicateGenes',
#        -max_retry_count => 1,
#        -hive_capacity => 200,
#        -batch_size    => 100,
#        -parameters => {
#          RNA_DB => {
#            'cdna_db' => ['fetch_all_biotypes'],  # fetch_all_biotypes for all biotypes
#          },
#          output_db => $self->o('lincrna_db'),
#          dna_db => $self->o('dna_db'),
#          cdna_db => $self->o('rnaseq_refine_db'),
#          biotype_output => $self->o('biotype_output'),
#        },
#        -rc_name    => '5GB',
#        -flow_into => {
#          1 => ['Hive_LincRNAFinder'],
#        }
#      },

##############################################################################
# LincRNA ANALYSES
##############################################################################

 #     {
 #       -logic_name => 'Hive_LincRNAFinder',
 #       -hive_capacity => 200,
 #       -batch_size    => 100,
 #       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAFinder',
 #       -max_retry_count => 1,
 #       -parameters => {
 #         NEW_SET_1_CDNA => {
 #           'output_db'  => ['fetch_all_biotypes'],
 #         },
 #         NEW_SET_2_PROT  => {
 #           'protein_coding_db' => $self->o('lincrna_protein_coding_set'),
 #         },
 #         FIND_SINGLE_EXON_LINCRNA_CANDIDATES => 1000, # I don't want single exon candidates!
 #         CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND => 1,
 #         MAXIMUM_TRANSLATION_LENGTH_RATIO => 99,
 #         MAX_TRANSLATIONS_PER_GENE => 20,
 #         OUTPUT_DB => 'lincrna_db',
 #         OUTPUT_BIOTYPE => 'lincRNA_finder_2round',
 #         WRITE_DEBUG_OUTPUT => 0,     # Set this to "0" to turn off debugging OR to "1000" to set it on.
 #         DEBUG_OUTPUT_DB    => 'output_db',    # where debug output (if any) will be written to
 #         protein_coding_db => $self->o('rnaseq_for_layer_db'),
 #         output_db => $self->o('lincrna_db'),
 #         dna_db => $self->o('dna_db'),
 #       },
#        -rc_name    => '5GB',
#        -flow_into => {
#                        1  => ['HiveDumpTranslations'],
#                       -1 => ['Hive_LincRNAFinder_himem'],
#                      },
#      },

#      {
#        -logic_name => 'Hive_LincRNAFinder_himem',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAFinder',
#        -batch_size => 5,
#        -parameters => {
#          NEW_SET_1_CDNA => {
#            'output_db'  => ['fetch_all_biotypes'],
#          },
#          NEW_SET_2_PROT  => {
#            'protein_coding_db' => $self->o('lincrna_protein_coding_set'),
#          },
#          FIND_SINGLE_EXON_LINCRNA_CANDIDATES => 1000, # I don't want single exon candidates!
#            CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND => 1,
#          MAXIMUM_TRANSLATION_LENGTH_RATIO => 99,
#          MAX_TRANSLATIONS_PER_GENE => 20,
#          OUTPUT_DB => 'lincrna_db',
#          OUTPUT_BIOTYPE => 'lincRNA_finder_2round',
#          WRITE_DEBUG_OUTPUT => 0,     # Set this to "0" to turn off debugging OR to "1000" to set it on.
#          DEBUG_OUTPUT_DB    => 'output_db',    # where debug output (if any) will be written to
#          protein_coding_db => $self->o('rnaseq_for_layer_db'),
#          output_db => $self->o('lincrna_db'),
#          dna_db => $self->o('dna_db'),
#        },
#        -rc_name    => '20GB',
#        -can_be_empty  => 1,
#        -flow_into => {
#                        1  => ['HiveDumpTranslations'],
#                      },
#      },

#      {
#        -logic_name => 'HiveDumpTranslations',
#        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpTranslations',
#        -batch_size    => 100,
#        -parameters => {
#                          dump_translations_script => catfile($self->o('ensembl_analysis_script'), 'protein', 'dump_translations.pl'),
#                          dna_db => $self->o('dna_db'),
#                          source_db => $self->o('lincrna_db'),
#                          file => $self->o('file_translations'),
#                       },
#        -rc_name    => '2GB',
#        -flow_into => {
#                       '2->A' => ['SplitDumpFiles'],
#                       'A->1' => ['Hive_LincRNAEvaluator'],
#                      },
#       },
#       {
#         -logic_name => 'SplitDumpFiles',
#         -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplitFasta',
#         -parameters => {
#                          fasta_file              => $self->o('file_translations'),
#                          out_dir                 => $self->o('lncrna_dir'),
#                          max_seqs_per_file       => $self->o('max_seqs_per_file'),
#                          max_seq_length_per_file => $self->o('max_seq_length_per_file'),
#                          max_files_per_directory => $self->o('max_files_per_directory'),
#                          max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
#                         },
#         -rc_name    => '2GB',
#         -flow_into     => {
#            2 => ['RunI5Lookup'],
#         }
#      },

       ### Here begins the running InterproScan
#      {
#        -logic_name    => 'RunI5Lookup',
#        -module        => 'Bio::EnsEMBL::Production::Pipeline::ProteinFeatures::InterProScan',
#        -hive_capacity => 250,
#        -batch_size    => 10,
#        -parameters    => {
#                             input_file                => '#split_file#',
#                             run_mode                  => 'lookup',
#                             interproscan_exe          => $self->o('interproscan_exe'),
#                             interproscan_applications => $self->o('interproscan_lookup_applications'),
#                             run_interproscan          => 1,
#                             escape_branch             => -1,
#                           },
#        -rc_name       => '6GB_registry',
#        -flow_into     => ['StoreProteinFeatures'],
#      },

#      {
#        -logic_name    => 'StoreProteinFeatures',
#        -max_retry_count => 1,
#        -module        => 'Bio::EnsEMBL::Production::Pipeline::ProteinFeatures::StoreProteinFeatures',
#        -parameters    => {
#                             species             => $self->o('species_name'),
#                             required_externalDb => $self->o('required_externalDb'),
#                            analyses        => $self->o('required_analysis'),
#                            pathway_sources => $self->o('pathway_sources'),
#                          },
#        -hive_capacity => 50,
#        -rc_name => '6GB_registry',
#      },


#      {
#        -logic_name => 'Hive_LincRNAEvaluator',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAEvaluator',
#        -max_retry_count => 0,
#        -hive_capacity => 100,
#        -batch_size    => 150,
#        -parameters => {
#          LINCRNA_DB => {
#            output_db => ['lincRNA_finder_2round'],
#          },
#          VALIDATION_DBS => {
#            protein_coding_db => $self->o('lincrna_protein_coding_set'),
#          },
#          EXCLUDE_SINGLE_EXON_LINCRNAS => 1000,
#          MAX_FRAMESHIFT_INTRON_LEN => 9,
#          EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS => 0,
#          MAX_TRANSCRIPTS_PER_CLUSTER => 3,
#          FINAL_OUTPUT_BIOTYPE => "lincRNA_pass_Eval_no_pfam",
#          FINAL_OUTPUT_DB      => 'lincrna_db',
#          MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB => 10000, # no validation db check for now. if you say yes here, you have to change the following parameters about validation
#          PROC_TRANS_HAVANA_LOGIC_NAME_STRING => 'havana',
#          OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME => 'ensembl_havana_gene',
#          UPDATE_SOURCE_DB => 'lincrna_db',
#          WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS => 1,
#          MARK_EXISTING_LINCRNA_IN_VALIDATION_DB => 1,
#          WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS => 1,
#          WRITE_REJECTED_NCRNAS => 1,
#          protein_coding_db => $self->o('rnaseq_for_layer_db'),
#          output_db => $self->o('lincrna_db'),
#          dna_db => $self->o('dna_db'),
#          cdna_db => $self->o('rnaseq_refine_db'),
#        },
#        -rc_name    => '5GB',

#       },


#       {
#         -logic_name => 'delete_duplicate_lincrna_genes',
#         -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
#         -parameters => {
#                          cmd => "perl ".$self->o('remove_duplicates_script_path')
#                                 ." -dbhost ".$self->o('lincrna_db','-host')
#                                 ." -dbuser ".$self->o('user')
#                                 ." -dbpass ".$self->o('password')
#                                 ." -dbname ".$self->o('lincrna_db','-dbname')
#                                 ." -dbport ".$self->o('lincrna_db','-port')
#                                 ." -dnadbhost ".$self->o('dna_db','-host')
#                                 ." -dnadbuser ".$self->o('user_r')
#                                 ." -dnadbname ".$self->o('dna_db','-dbname')
#                                 ." -dnadbport ".$self->o('dna_db','-port'),
#                      },
#          -max_retry_count => 0,
#          -rc_name => 'default',
#          -flow_into     => {
#            1 => ['Hive_LincRNAAftCheck_pi'],
#          },

#        },


#      {
#        -logic_name => 'Hive_LincRNAAftCheck_pi',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincAfterChecks',
#        -parameters => {
#                           Final_BIOTYPE_TO_CHECK  => 'lincRNA_pass_Eval_no_pfam',
#                           output_db => $self->o('lincrna_db'),
#                           dna_db => $self->o('dna_db'),
#                           file_l => $self->o('file_for_length'),
#                           file_is => $self->o('file_for_introns_support'),
#                           file_b => $self->o('file_for_biotypes'),
#                           assembly_name => $self->o('assembly_name'),
#                           update_database => $self->o('update_database'),
#                        },
#        -rc_name    => '8GB',
#      },


########################################################################
#
# Minimap2 long read analyses
#
########################################################################


      {
        -logic_name => 'create_long_read_final_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db   => $self->o('dna_db'),
          target_db   => $self->o('long_read_final_db'),
          create_type => 'clone',
        },
        -rc_name    => 'default',
        -flow_into  => {
          '1->A' => ['fan_long_read'],
          'A->1' => ['cluster_ig_tr_genes'],
        },
      },


      {
        -logic_name => 'fan_long_read',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_long_read#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
	},
        -flow_into  => {
          1 => ['create_long_read_dir'],
	},
        -rc_name => 'default',
      },


      {
        -logic_name => 'create_long_read_dir',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => 'default',
        -parameters => {
          cmd => 'if [ ! -e "'.$self->o('long_read_fastq_dir').'" ]; then mkdir -p '.$self->o('long_read_fastq_dir').';fi',
        },
        -flow_into => {
          '1' => ['create_long_read_initial_db'],
        },
        -rc_name => 'default',
      },


      {
        -logic_name => 'create_long_read_initial_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
	     -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('long_read_initial_db'),
                         create_type => 'clone',
                       },
        -rc_name    => '1GB',
        -flow_into => {
                        '1' => ['create_minimap2_index'],
                      },
      },


      {
        -logic_name => 'create_minimap2_index',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ ! -e "'.$self->o('minimap2_genome_index').'" ]; then '.$self->o('minimap2_path').
                 ' -d '.$self->o('minimap2_genome_index').' '.$self->o('faidx_genome_file').';fi',
        },
        -flow_into  => {
          1 => ['check_index_not_empty'],
	},
        -rc_name => '20GB',
      },


      {
        -logic_name => 'check_index_not_empty',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ -s "'.$self->o('minimap2_genome_index').'" ]; then exit 0; else exit 42;fi',
                         return_codes_2_branches => {'42' => 2},
	},
        -flow_into  => {
         '1->A' => ['create_lr_fastq_download_jobs'],
         'A->1' => ['create_collapse_db'],
        },
        -rc_name => 'default',
      },


      {
        -logic_name => 'create_lr_fastq_download_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputfile    => $self->o('long_read_summary_file'),
          column_names => $self->o('long_read_columns'),
          delimiter => '\t',
        },
        -flow_into => {
          2 => {'download_long_read_fastq' => {'iid' => '#filename#'}},
        },
      },


      {
        -logic_name => 'download_long_read_fastq',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs',
        -parameters =>{
          ftp_base_url => $self->o('rnaseq_ftp_base'),
          input_dir => $self->o('long_read_fastq_dir'),
          samtools_path => $self->o('samtools_path'),
          decompress => 1,
          create_faidx => 1,
        },
        -rc_name => '1GB',
        -analysis_capacity => 50,
        -flow_into => {
          1 => {'generate_minimap2_jobs' => {'fastq_file' => $self->o('long_read_fastq_dir').'/'.'#fastq_file#'}},
        },
      },


      {
        -logic_name => 'generate_minimap2_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'fastq_range',
                         batch_size => $self->o('minimap2_batch_size'),
                       },
        -rc_name      => '2GB',
        -flow_into => {
                        2 => {'minimap2' => {'input_file' => '#fastq_file#','iid' => '#iid#'}},
                      },
      },


      {
        -logic_name => 'minimap2',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2',
        -parameters => {
                         genome_file => $self->o('faidx_genome_file'),
                         long_read_summary_file => $self->o('long_read_summary_file'),
                         long_read_summary_file_genus => $self->o('long_read_summary_file_genus'),
                         minimap2_genome_index => $self->o('minimap2_genome_index'),
                         minimap2_path => $self->o('minimap2_path'),
                         paftools_path => $self->o('paftools_path'),
                         target_db => $self->o('long_read_initial_db'),
                         logic_name => 'minimap2',
                         module     => 'Minimap2',
                       },
        -rc_name => '15GB',
        -flow_into => {
                        -1 => {'minimap2_himem' => {'input_file' => '#input_file#','iid' => '#iid#'}},
                      },
      },


      {
        -logic_name => 'minimap2_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2',
        -parameters => {
                         genome_file => $self->o('faidx_genome_file'),
                         long_read_summary_file => $self->o('long_read_summary_file'),
                         long_read_summary_file_genus => $self->o('long_read_summary_file_genus'),
                         minimap2_genome_index => $self->o('minimap2_genome_index'),
                         minimap2_path => $self->o('minimap2_path'),
                         paftools_path => $self->o('paftools_path'),
                         target_db => $self->o('long_read_initial_db'),
                         logic_name => 'minimap2',
                         module     => 'Minimap2',
                       },
        -rc_name => '25GB',
      },


      {
        -logic_name => 'create_collapse_db',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('dna_db'),
          target_db => $self->o('long_read_collapse_db'),
          create_type => 'clone',
         },
         -rc_name => 'default',
         -max_retry_count => 0,
         -flow_into => {
           1 => ['generate_collapse_jobs'],
         }
      },


#      {
#        -logic_name => 'create_long_read_blast_db',
#        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
#        -parameters => {
#          source_db => $self->o('dna_db'),
#          target_db => $self->o('long_read_blast_db'),
#          create_type => 'clone',
#        },
#        -rc_name => 'default',
#        -max_retry_count => 0,
#        -flow_into => {
#          1 => ['generate_collapse_jobs'],
#        }
#      },


      {
        -logic_name => 'generate_collapse_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          target_db        => $self->o('long_read_collapse_db'),
          feature_dbs => [$self->o('long_read_initial_db')],
          coord_system_name => 'toplevel',
          iid_type => 'stranded_slice',
          feature_constraint => 1,
          feature_type => 'gene',
          top_level => 1,
        },
        -rc_name      => 'default',
        -max_retry_count => 1,
        -flow_into => {
          '2->A' => ['split_lr_slices_on_intergenic'],
	  'A->1' => ['classify_long_read_models'],
        },
      },


      {
        -logic_name => 'split_lr_slices_on_intergenic',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions',
        -parameters => {
	  dna_db => $self->o('dna_db'),
          input_gene_dbs => [$self->o('long_read_initial_db')],
          iid_type => 'slice',
          use_strand => 1,
        },
        -batch_size => 100,
        -rc_name    => '5GB',
        -flow_into => {
          2 => {'collapse_transcripts' => {'slice_strand' => '#slice_strand#','iid' => '#iid#'}},
        },
      },


      {
        -logic_name => 'collapse_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer',
        -parameters => {
	  target_db        => $self->o('long_read_collapse_db'),
          dna_db        => $self->o('dna_db'),
          source_dbs        => [$self->o('long_read_initial_db')],
          biotypes => ["isoseq","cdna"],
          reduce_large_clusters => 1,
        },
        -rc_name      => '5GB',
        -flow_into => {
           1 => ['blast_long_read'],
          -1 => {'collapse_transcripts_20GB' => {'slice_strand' => '#slice_strand#','iid' => '#iid#'}},
        },
        -batch_size => 100,
        -analysis_capacity => 1000,
      },


      {
        -logic_name => 'collapse_transcripts_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer',
        -parameters => {
          target_db        => $self->o('long_read_collapse_db'),
          dna_db        => $self->o('dna_db'),
          source_dbs        => [$self->o('long_read_initial_db')],
          biotypes => ["isoseq","cdna"],
          reduce_large_clusters => 1,
        },
        -rc_name      => '20GB',
        -flow_into => {
          1 => {'blast_long_read' => {'slice_strand' => '#slice_strand#','iid' => '#iid#'}},
         -1 => {'failed_collapse' => {'slice_strand' => '#slice_strand#','iid' => '#iid#'}},
        },
        -batch_size => 10,
        -analysis_capacity => 1000,
      },


      {
        -logic_name => 'failed_collapse',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer',
        -parameters => {
          target_db         => $self->o('long_read_collapse_db'),
          dna_db            => $self->o('dna_db'),
          source_dbs        => [$self->o('long_read_initial_db')],
          biotypes => ["isoseq","cdna"],
          copy_only => 1,
        },
        -rc_name      => '10GB',
        -flow_into => {
          1 => {'blast_long_read' => {'slice_strand' => '#slice_strand#','iid' => '#iid#'}},
        },
      },


      {
        -logic_name => 'blast_long_read',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
        -parameters => {
          input_db => $self->o('long_read_collapse_db'),
          output_db => $self->o('long_read_final_db'),
          source_db => $self->o('long_read_collapse_db'),
          target_db => $self->o('long_read_final_db'),
          dna_db => $self->o('dna_db'),
          indicate_index => $self->o('protein_blast_index'),
          uniprot_index => [$self->o('protein_blast_db')],
          blast_program => $self->o('uniprot_blast_exe_path'),
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
          commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40 -seg no',
        },
        -rc_name => 'blast',
        -flow_into => {
          -1 => ['blast_long_read_10G'],
#          1 => ['intron_check'],
        },
      },


      {
        -logic_name => 'blast_long_read_10G',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
        -parameters => {
          input_db => $self->o('long_read_collapse_db'),
          output_db => $self->o('long_read_final_db'),
          dna_db => $self->o('dna_db'),
          source_db => $self->o('long_read_collapse_db'),
          target_db => $self->o('long_read_final_db'),
          indicate_index => $self->o('protein_blast_index'),
          uniprot_index => [$self->o('protein_blast_db')],
          blast_program => $self->o('uniprot_blast_exe_path'),
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
          commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40 -seg no',
        },
        -rc_name => 'blast10GB',
#        -flow_into => {
#          1 => ['intron_check'],
#        },
      },


#      {
#        -logic_name => 'intron_check',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck',
#        -parameters => {
#          source_db => $self->o('long_read_blast_db'),
#          target_db => $self->o('long_read_final_db'),
#          intron_db => $self->o('rnaseq_refine_db'),
#          dna_db => $self->o('dna_db'),
#        },
#        -rc_name    => '2GB',
#        -flow_into => {
#          1 => ['intron_check_10GB'],
#        },
#      },


#      {
#        -logic_name => 'intron_check_10GB',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck',
#        -parameters => {
#          source_db => $self->o('long_read_blast_db'),
#          target_db => $self->o('long_read_final_db'),
#          intron_db => $self->o('rnaseq_refine_db'),
#          dna_db => $self->o('dna_db'),
#        },
#        -rc_name    => '10GB',
#      },


      {
        -logic_name => 'classify_long_read_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
          classification_type => 'long_read',
          update_gene_biotype => 1,
          target_db => $self->o('long_read_final_db'),
        },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'cluster_ig_tr_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCollapseIGTR',
        -parameters => {
                         target_db => $self->o('ig_tr_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'ig_tr_gene',
                         logic_names_to_cluster => ['ig_tr_gene','ig_tr_gene_not_best'],
                       },
        -rc_name    => 'genblast',
        -flow_into => {
                        1 => ['update_ig_tr_biotypes'],
                      },
      },


      {
        -logic_name => 'update_ig_tr_biotypes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('ig_tr_db'),
          sql => [
            'UPDATE transcript JOIN analysis USING(analysis_id) SET biotype = CONCAT(biotype, "_pre_collapse")'.
              ' WHERE logic_name != "ig_tr_collapse"',
            'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype = transcript.biotype',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
          '1->A' => ['fan_genblast_rnaseq_support'],
          'A->1' => ['create_layering_output_db'],
        },
      },


      {
        -logic_name => 'fan_genblast_rnaseq_support',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_rnaseq#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
         '1->A' => ['create_genblast_rnaseq_support_db'],
         'A->1' => ['create_genblast_rnaseq_nr_db'],
        },
      },


     {
        -logic_name => 'create_genblast_rnaseq_support_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_rnaseq_support_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_genblast_rnaseq_slice_ids'],
                      },
     },


     {
        -logic_name => 'create_genblast_rnaseq_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db  => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => 0,
                       },
        -rc_name => 'default',
        -flow_into => {
                        '2' => ['genblast_rnaseq_support'],
                      },
     },


    {
        -logic_name => 'genblast_rnaseq_support',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         source_db => $self->o('genblast_db'),
                         intron_db => $self->o('rnaseq_refine_db'),
                         target_db => $self->o('genblast_db'),
                         logic_name => 'genblast_rnaseq_support',
                         classify_by_count => 1,
                         update_genes => 0,
                         module => 'HiveHomologyRNASeqIntronsCheck',
                       },
        -rc_name => '4GB',
        -flow_into => {
                        '-1' => ['genblast_rnaseq_support_himem'],
                      },

     },


     {
        -logic_name => 'genblast_rnaseq_support_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         source_db => $self->o('genblast_db'),
                         intron_db => $self->o('rnaseq_refine_db'),
                         target_db => $self->o('genblast_db'),
                         logic_name => 'genblast_rnaseq_support',
                         classify_by_count => 1,
                         update_genes => 0,
                         module => 'HiveHomologyRNASeqIntronsCheck',
                       },
        -rc_name => '8GB',
      },


      {
        -logic_name => 'create_genblast_rnaseq_nr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('genblast_db'),
                         target_db => $self->o('genblast_rnaseq_support_nr_db'),
                         create_type => 'copy',
                         force_drop => 1,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['create_genblast_rnaseq_nr_slices'],
                      },
      },


      {
        -logic_name => 'create_genblast_rnaseq_nr_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 20000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                         batch_slice_ids => 1,
                         batch_target_size => 20000000,
                       },
        -rc_name    => '2GB',
        -flow_into => {
                         '2'    => ['remove_redundant_genblast_rnaseq_genes'],
                      },
     },


    {
      -logic_name => 'remove_redundant_genblast_rnaseq_genes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db => $self->o('genblast_rnaseq_support_nr_db'),
        target_type => 'biotype_priority',
        layers => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic', $self->o('uniprot_set'), undef, 'ARRAY'),
      },
      -rc_name          => '5GB',
    },


############################################################################
#
# Finalisation analyses
#
############################################################################
      {
        -logic_name => 'create_layering_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('layering_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_utr_db'],
                      },
      },


      {
        -logic_name => 'create_utr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('utr_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_genebuilder_db'],
                      },
      },


      {
        -logic_name => 'create_genebuilder_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('genebuilder_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                       '1->A' => ['create_toplevel_slices'],
                       'A->1' => ['layer_annotation_sanity_checks'],
                      },
      },


      {
        -logic_name => 'create_toplevel_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'slice',
                         coord_system_name => 'toplevel',
                         include_non_reference => 0,
                         top_level => 1,
                         # These options will create only slices that have a gene on the slice in one of the feature dbs
                         feature_constraint => 1,
                         feature_type => 'gene',
                         feature_dbs => [$self->o('genblast_db'),$self->o('selected_projection_db'),$self->o('rnaseq_for_layer_db')],
                      },
        -flow_into => {
                       '2'    => ['split_slices_on_intergenic'],
                      },

        -rc_name    => 'default',
      },


      {
        -logic_name => 'split_slices_on_intergenic',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         input_gene_dbs => [@{$self->default_options->{'layering_input_gene_dbs'}},$self->default_options->{'cdna_db'}],
                         iid_type => 'slice',
                       },
        -batch_size => 100,
        -hive_capacity => $self->hive_capacity_classes->{'hc_medium'},
        -rc_name    => '5GB',
        -flow_into => {
                       '2' => ['layer_annotation'],
                      },
      },


      {
        -logic_name => 'layer_annotation',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLayerAnnotation',
        -parameters => {
                         use_genome_flatfile => 1,
                         genome_file => $self->o('faidx_genome_file'),
                         logic_name => 'layer_annotation',
                         module     => 'HiveLayerAnnotation',
                         TARGETDB_REF => $self->o('layering_db'),
                         SOURCEDB_REFS => $self->o('layering_input_gene_dbs'),

                         # Filtering is using done at the exon-overlap level
                         # When no FILTER exists in this file, this is the default behaviour
                         # If you would like to filter in a different way, please specify filter
                         #FILTER => 'Bio::EnsEMBL::Analysis::Tools::GenomeOverlapFilter',
                         #FILTER => 'Bio::EnsEMBL::Analysis::Tools::AllExonOverlapFilter',
                         FILTER => 'Bio::EnsEMBL::Analysis::Tools::CodingExonOverlapFilter',
                         # ordered list of annotation layers. Genes from lower layers
                         # are only retained if they do not "interfere" with genes from
                         # higher layers. Genes in "Discard" layers are when assessing
                         # interference, but are not written to the final database
                         LAYERS => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic', $self->o('uniprot_set'), undef, 'ARRAY'),
                       },
        -rc_name    => '4GB',
        -flow_into  => {
                         '1->A' => ['run_utr_addition'],
                         'A->1' => ['genebuilder'],
                       },
        -batch_size => 50,
      },


      {
        -logic_name => 'run_utr_addition',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
        -parameters => {
                         logic_name => 'utr_addition',
                         dna_db => $self->o('dna_db'),
                         donor_dbs => $self->o('utr_donor_dbs'),
                         acceptor_dbs => $self->o('utr_acceptor_dbs'),
                         utr_biotype_priorities => $self->o('utr_biotype_priorities'),
                         target_db => $self->o('utr_db'),
                         iid_type => 'slice',
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => '5GB',
        -flow_into => {
                        -1 => ['run_utr_addition_10GB'],
                      },

     },


     {
        -logic_name => 'run_utr_addition_10GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
        -parameters => {
                         logic_name => 'utr_addition',
                         dna_db => $self->o('dna_db'),
                         donor_dbs => $self->o('utr_donor_dbs'),
                         acceptor_dbs => $self->o('utr_acceptor_dbs'),
                         utr_biotype_priorities => $self->o('utr_biotype_priorities'),
                         target_db => $self->o('utr_db'),
                         iid_type => 'slice',
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => '10GB',
        -flow_into => {
                        -1 => ['run_utr_addition_30GB'],
                      },
     },


     {
        -logic_name => 'run_utr_addition_30GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
        -parameters => {
                         logic_name => 'utr_addition',
                         dna_db => $self->o('dna_db'),
                         donor_dbs => $self->o('utr_donor_dbs'),
                         acceptor_dbs => $self->o('utr_acceptor_dbs'),
                         utr_biotype_priorities => $self->o('utr_biotype_priorities'),
                         target_db => $self->o('utr_db'),
                         iid_type => 'slice',
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => '30GB',
        -flow_into => {
                        -1 => ['utr_memory_failover'],
                      },
     },


     {
        -logic_name => 'utr_memory_failover',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
        -parameters => {
                         logic_name => 'utr_addition',
                         dna_db => $self->o('dna_db'),
                         donor_dbs => $self->o('utr_donor_dbs'),
                         acceptor_dbs => $self->o('utr_acceptor_dbs'),
                         utr_biotype_priorities => $self->o('utr_biotype_priorities'),
                         target_db => $self->o('utr_db'),
                         iid_type => 'slice',
                         copy_only => 1,
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => '10GB',
     },


     {
        -logic_name => 'genebuilder',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGeneBuilder',
        -parameters => {
                         source_db => $self->o('utr_db'),
                         target_db => $self->o('genebuilder_db'),
                         dna_db     => $self->o('dna_db'),
                         logic_name => 'ensembl',
                         module     => 'HiveGeneBuilder',
                         INPUT_GENES => {
                           source_db => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic',
                                                               $self->o('uniprot_set'), undef, 'ARRAY'),
                         },
                         OUTPUT_BIOTYPE => 'protein_coding',
                         MAX_TRANSCRIPTS_PER_CLUSTER => 10,
                         MIN_SHORT_INTRON_LEN => 7, #introns shorter than this seem
                         #to be real frame shifts and shoudn't be ignored
                         MAX_SHORT_INTRON_LEN => 15,
                         BLESSED_BIOTYPES => {
                                              'ccds_gene' => 1,
                                              'IG_C_gene' => 1,
                                              'IG_J_gene' => 1,
                                              'IG_V_gene' => 1,
                                              'IG_D_gene' => 1,
                                              'TR_C_gene' => 1,
                                              'TR_J_gene' => 1,
                                              'TR_V_gene' => 1,
                                              'TR_D_gene' => 1,
                                             },
                         #the biotypes of the best genes always to be kept
                         MAX_EXON_LENGTH => 20000,
                         #if the coding_only flag is set to 1, the transcript clustering into genes is done over coding exons only
                         # the current standard way is to cluster only on coding exons
                         CODING_ONLY => 1,
                       },
        -rc_name    => '4GB',
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'layer_annotation_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('layering_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'layer'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['genebuilder_sanity_checks'],
                      },
      },


      {
        -logic_name => 'genebuilder_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('genebuilder_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'genebuilder'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['restore_ig_tr_biotypes'],
                      },
      },


      {
        -logic_name => 'restore_ig_tr_biotypes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('genebuilder_db'),
          sql => [
            'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype = transcript.biotype'.
              ' WHERE transcript.biotype LIKE "IG\_%" OR transcript.biotype LIKE "TR\_%"',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_pseudogene_db'],
                      },
      },


      {
        -logic_name => 'create_pseudogene_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('pseudogene_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['update_lncrna_biotypes'],
                      },
      },


      {
        -logic_name => 'update_lncrna_biotypes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('genebuilder_db'),
          sql => [
            'UPDATE transcript SET biotype="pre_lncRNA" WHERE biotype IN ("rnaseq_merged","rnaseq_tissue","cdna")',
            'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype="pre_lncRNA" WHERE transcript.biotype="pre_lncRNA"',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['pseudogenes'],
                      },
      },


      {
        -logic_name => 'pseudogenes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePseudogenes',
        -parameters => {
                         single_multi_file => 1,
                         output_path => $self->o('output_path').'/pseudogenes/',
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'pseudogenes',
                         module     => 'HivePseudogenes',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic','pseudogenes')},
                       },

	     -rc_name    => '30GB',
	     -flow_into => {
			    1 => ['format_blast_db'],
                      },
      },


      {
        -logic_name => 'format_blast_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "'.$self->o('blast_type').
                                '" = "ncbi" ];then makeblastdb -dbtype nucl -in '.
			        $self->o('output_path').'/pseudogenes/all_multi_exon_genes.fasta;'.
	                        ' else xdformat -n '.$self->o('output_path').'/pseudogenes/all_multi_exon_genes.fasta;fi'
		       },
         -rc_name => 'default',
         -flow_into => {
                         1 => ['spliced_elsewhere'],
                       },
      },


      {
        -logic_name => 'spliced_elsewhere',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplicedElsewhere',
        -parameters => {
                         multi_exon_db_path => $self->o('output_path').'/pseudogenes/',
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'spliced_elsewhere',
                         module     => 'HiveSplicedElsewhere',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic','pseudogenes')},
                       },
        -rc_name          => '5GB',
        -flow_into => {
                        1 => ['create_final_geneset_db'],
                      },
      },


      {
        -logic_name => 'create_final_geneset_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('pseudogene_db'),
                         target_db => $self->o('final_geneset_db'),
                         create_type => 'copy',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['filter_lncrnas'],
                      },
      },


      {
        -logic_name => 'filter_lncrnas',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFilterlncRNAs',
        -parameters => {
                         input_gene_db => $self->o('final_geneset_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'filter_lncrnas',
                         module     => 'HiveFilterlncRNAs',
                       },
        -rc_name          => '3GB',
        -flow_into => {
                        1 => ['change_biotype_for_weak_cds'],
                      },
      },


      {
        -logic_name => 'change_biotype_for_weak_cds',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('final_geneset_db'),
          sql => [
            'UPDATE transcript JOIN transcript_supporting_feature USING(transcript_id) JOIN protein_align_feature ON feature_id=protein_align_feature_id SET biotype="low_coverage" WHERE feature_type="protein_align_feature" AND hcoverage < 50 AND biotype not like "pseudo%"',
            'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype="low_coverage" WHERE transcript.biotype="low_coverage"',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['update_rnaseq_ise_logic_names'],
                      },
      },


      {
        -logic_name => 'update_rnaseq_ise_logic_names',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('final_geneset_db'),
          sql => [
            'UPDATE analysis SET logic_name = REPLACE(logic_name, "_rnaseq_gene", "_rnaseq_ise")',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['run_cleaner'],
                      },
      },


      {
        -logic_name => 'run_cleaner',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanGeneset',
        -parameters => {
                         skip_analysis => $self->o('skip_cleaning'),
                         input_db => $self->o('final_geneset_db'),
                         dna_db => $self->o('dna_db'),
                         output_path => $self->o('output_path').'/clean_genes/',
                         blessed_biotypes => $self->o('cleaning_blessed_biotypes'),
                         flagged_redundancy_coverage_threshold => 95,
                         general_redundancy_coverage_threshold => 95,
                       },
        -rc_name    => '4GB',
        -flow_into => {
                        '1' => ['delete_flagged_transcripts'],
                      },
      },


      {
       -logic_name => 'delete_flagged_transcripts',
       -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts',
       -parameters => {
                        skip_analysis => $self->o('skip_cleaning'),
                        dbhost => $self->o('final_geneset_db','-host'),
                        dbname => $self->o('final_geneset_db','-dbname'),
                        dbuser => $self->o('user'),
                        dbpass => $self->o('password'),
                        dbport => $self->o('final_geneset_db','-port'),
                        transcript_ids_file => catfile($self->o('output_path'), 'clean_genes', 'transcript_ids_to_remove.txt'),
                        delete_transcripts_path => catdir($self->o('ensembl_analysis_script'), 'genebuild/'),
                        delete_genes_path => catdir($self->o('ensembl_analysis_script'), 'genebuild/'),
                        delete_transcripts_script_name => '/delete_transcripts.pl',
                        delete_genes_script_name => '/delete_genes.pl',
                        output_path => catdir($self->o('output_path'), 'clean_genes'),
                        output_file_name => 'delete_transcripts.out',
                      },
        -max_retry_count => 0,
        -flow_into => {
          1 => ['transfer_ncrnas'],
        },
     },



     {
       -logic_name => 'transfer_ncrnas',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
                        cmd => 'perl '.catfile($self->o('ensembl_analysis_script'), 'genebuild', 'copy_genes.pl').
                                     ' -sourcehost '.$self->o('ncrna_db','-host').
                                     ' -sourceuser '.$self->o('user_r').
                                     ' -sourceport '.$self->o('ncrna_db','-port').
                                     ' -sourcedbname '.$self->o('ncrna_db','-dbname').
                                     ' -dnauser '.$self->o('user_r').
                                     ' -dnahost '.$self->o('dna_db','-host').
                                     ' -dnaport '.$self->o('dna_db','-port').
                                     ' -dnadbname '.$self->o('dna_db','-dbname').
                                     ' -targetuser '.$self->o('user').
                                     ' -targetpass '.$self->o('password').
                                     ' -targethost '.$self->o('final_geneset_db','-host').
                                     ' -targetport '.$self->o('final_geneset_db','-port').
                                     ' -targetdbname '.$self->o('final_geneset_db','-dbname').
                                     ' -all'
                      },
        -rc_name => 'default',
        -flow_into => {
          '1' => ['delete_duplicate_genes'],
        },
     },



     {
       -logic_name => 'delete_duplicate_genes',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
                        cmd => "perl ".$self->o('remove_duplicates_script_path')
                               ." -dbhost ".$self->o('final_geneset_db','-host')
                               ." -dbuser ".$self->o('user')
                               ." -dbpass ".$self->o('password')
                               ." -dbname ".$self->o('final_geneset_db','-dbname')
                               ." -dbport ".$self->o('final_geneset_db','-port')
                               ." -dnadbhost ".$self->o('dna_db','-host')
                               ." -dnadbuser ".$self->o('user_r')
                               ." -dnadbname ".$self->o('dna_db','-dbname')
                               ." -dnadbport ".$self->o('dna_db','-port'),
                     },
        -max_retry_count => 0,
        -rc_name => '2GB',
        -flow_into => {
                        '1' => ['final_db_sanity_checks'],
                      },
     },



############################################################################
#
# Transfer projected lincRNA and pseudogenes
#
############################################################################
#     {
#       -logic_name => 'dummy_projected_lincrna_pseudos',
#       -module => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
#       -parameters => {},
#       -rc_name => 'default',
#       -flow_into  => {
#         '1->A' => ['fan_transfer_projected_genes'],
#         'A->1' => ['delete_duplicate_genes'],
#       },
#     },
#
#
#     {
#       -logic_name => 'fan_transfer_projected_genes',
#       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
#       -parameters => {
#         cmd => 'if [ "#skip_projection#" -ne -1 ]; then exit 42; else exit 0;fi', # To enable this part, replace -1 with 0
#         return_codes_2_branches => {'42' => 2},
#       },
#       -rc_name => 'default',
#       -flow_into  => {
#         1 => ['create_projected_lincrna_ids_to_copy'],
#       },
#     },
#
#
#     {
#       -logic_name => 'create_projected_lincrna_ids_to_copy',
#       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#       -parameters => {
#                        target_db    => $self->o('projection_lincrna_db'),
#                        iid_type     => 'feature_id',
#                        feature_type => 'gene',
#                        batch_size   => 500,
#                     },
#       -flow_into => {
#         '2->A' => ['transfer_projected_lincrnas'],
#         'A->1' => ['create_projected_pseudogene_ids_to_copy'],
#       },
#       -rc_name    => 'default',
#     },
#
#
#     {
#       -logic_name => 'transfer_projected_lincrnas',
#       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
#       -parameters => {
#         copy_genes_directly => 1,
#         source_db => $self->o('projection_lincrna_db'),
#         dna_db => $self->o('dna_db'),
#         target_db => $self->o('final_geneset_db'),
#         filter_on_overlap => 1,
#         filter_on_strand  => 0,
#         overlap_filter_type => 'genomic_overlap',
#         filter_against_biotypes => {
#           protein_coding => 1,
#           IG_C_gene => 1,
#           IG_V_gene => 1,
#           TR_C_gene => 1,
#           TR_J_gene => 1,
#           TR_V_gene => 1,
#           processed_pseudogene => 1,
#           pseudogene => 1,
#         },
#       },
#       -rc_name    => 'default',
#     },
#
#
#     {
#       -logic_name => 'create_projected_pseudogene_ids_to_copy',
#       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#       -parameters => {
#                        target_db    => $self->o('projection_pseudogene_db'),
#                        iid_type     => 'feature_id',
#                        feature_type => 'gene',
#                        batch_size   => 500,
#                      },
#       -flow_into => {
#                       '2' => ['transfer_projected_pseudogenes'],
#                     },
#       -rc_name    => 'default',
#     },
#
#
#      {
#        -logic_name => 'transfer_projected_pseudogenes',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
#        -parameters => {
#                         copy_genes_directly => 1,
#                         source_db => $self->o('projection_pseudogene_db'),
#                         dna_db => $self->o('dna_db'),
#                         target_db => $self->o('final_geneset_db'),
#                         filter_on_overlap => 1,
#                         filter_on_strand  => 1,
#                         overlap_filter_type => 'genomic_overlap',
#                         filter_against_biotypes => {
#                                                      protein_coding => 1,
#                                                      IG_C_gene => 1,
#                                                      IG_V_gene => 1,
#                                                      TR_C_gene => 1,
#                                                      TR_J_gene => 1,
#                                                      TR_V_gene => 1,
#                                                      processed_pseudogene => 1,
#                                                      pseudogene => 1,
#                                                      lincRNA => 1,
#                                                    },
#                       },
#        -rc_name    => 'default',
#      },


      {
        -logic_name => 'final_db_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('final_geneset_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'final'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        '1->A' => ['create_gene_ids_to_copy'],
                        'A->1' => ['update_biotypes_and_analyses'],
                      },
      },


      {
        -logic_name => 'create_gene_ids_to_copy',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db    => $self->o('final_geneset_db'),
                         iid_type     => 'feature_id',
                         feature_type => 'gene',
                         batch_size   => 100,
                      },
        -flow_into => {
                       '2' => ['copy_genes_to_core'],
                      },

        -rc_name    => 'default',
      },


      {
        -logic_name => 'copy_genes_to_core',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
        -parameters => {
                         copy_genes_directly => 1,
                         source_db => $self->o('final_geneset_db'),
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         copy_biotypes_to_ignore => $self->o('copy_biotypes_to_ignore'),
                       },
        -rc_name    => '1GB',
      },


      {
        -logic_name => 'update_biotypes_and_analyses',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'UPDATE gene SET biotype = "protein_coding" WHERE biotype = "ensembl"',
            'UPDATE gene SET biotype = "vaultRNA" WHERE biotype = "Vault_RNA"',
            'UPDATE analysis set logic_name="cdna2genome" where logic_name="best_targetted"',
            'UPDATE gene SET analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = "ensembl")'.
              ' WHERE analysis_id IN'.
              ' (SELECT analysis_id FROM analysis'.
              ' WHERE logic_name IN ("spliced_elsewhere","pseudogenes","genblast","genblast_not_best","project_transcripts","project_pseudogene","project_lincrna","cesar"))',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.analysis_id = gene.analysis_id',
            'INSERT IGNORE into analysis (created,logic_name,db) VALUES (NOW(),"other_protein","uniprot")',
            'INSERT IGNORE into analysis (created,logic_name,db) VALUES (NOW(),"projected_transcript","'.$self->o('projection_source_db_name').'")',
            'UPDATE protein_align_feature SET analysis_id ='.
              '(SELECT analysis_id FROM analysis WHERE logic_name = "projected_transcript") WHERE analysis_id IN '.
              '(SELECT analysis_id FROM analysis WHERE logic_name IN ("project_transcripts","cesar"))',
            'UPDATE protein_align_feature SET analysis_id ='.
              '(SELECT analysis_id FROM analysis WHERE logic_name = "other_protein") WHERE analysis_id NOT IN'.
              '(SELECT analysis_id FROM analysis WHERE logic_name IN ("uniprot","projected_transcript"))',
            'UPDATE dna_align_feature SET analysis_id ='.
              '(SELECT analysis_id FROM analysis WHERE logic_name = "projected_transcript") WHERE analysis_id IN'.
              '(SELECT analysis_id FROM analysis WHERE logic_name IN ("project_lincrna","project_pseudogene"))',
            'UPDATE repeat_feature SET repeat_start = 1 WHERE repeat_start < 1',
            'UPDATE repeat_feature SET repeat_end = 1 WHERE repeat_end < 1',
            'UPDATE gene SET analysis_id=(select analysis_id from analysis where logic_name="ensembl") WHERE analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="filter_lncrnas")',
            'UPDATE transcript SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ensembl") WHERE analysis_id=(SELECT analysis_id from analysis WHERE logic_name="filter_lncrnas")',
            'DELETE FROM analysis WHERE logic_name="filter_lncrnas"',
            'UPDATE gene SET display_xref_id=NULL',
            'UPDATE transcript SET display_xref_id=NULL',
          ],
        },
        -rc_name    => 'default',
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
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbpattern '.$self->o('reference_db','-dbname')
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
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbname '.$self->o('reference_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_frameshift_introns'] },
      },

      {
        -logic_name => 'set_frameshift_introns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('frameshift_attrib_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbpattern '.$self->o('reference_db','-dbname')
                       },
        -rc_name => '10GB',
        -flow_into => { 1 => ['set_canonical_transcripts'] },
      },

      {
        -logic_name => 'set_canonical_transcripts',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('select_canonical_script').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbhost '.$self->o('reference_db','-host').
                                ' -dbport '.$self->o('reference_db','-port').
                                ' -dbname '.$self->o('reference_db','-dbname').
                                ' -coord toplevel -write'
                       },
        -rc_name => '10GB',
        -flow_into => { 1 => ['null_columns'] },
      },


      {
        -logic_name => 'null_columns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'UPDATE transcript set stable_id = NULL',
            'UPDATE translation set stable_id = NULL',
            'UPDATE exon set stable_id = NULL',
            'UPDATE protein_align_feature set external_db_id = NULL',
            'UPDATE dna_align_feature set external_db_id = NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['run_stable_ids'],
                      },
      },


      {
        -logic_name => 'run_stable_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
        -parameters => {
                         enscode_root_dir => $self->o('enscode_root_dir'),
                         mapping_required => $self->o('mapping_required'),
                         target_db => $self->o('reference_db'),
                         mapping_db => $self->o('mapping_db'),
                         id_start => $self->o('stable_id_prefix').$self->o('stable_id_start'),
                         output_path => $self->o('output_path'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['backup_core_db_pre_optimise'],
                      },
      },


      {
        # Creates a reference db for each species
        -logic_name => 'backup_core_db_pre_optimise',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
        -parameters => {
                         src_db_conn => $self->o('dna_db'),
                         output_file => catfile($self->o('output_path'), 'core_post_stable_idsbak.sql.gz'),
                         dump_options => $self->o('mysql_dump_options'),
                       },
        -rc_name    => 'default',
        -flow_into => { 1 => ['load_external_db_ids_and_optimise_af_tables'] },
      },


      {
        -logic_name => 'load_external_db_ids_and_optimise_af_tables',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('load_optimise_script').
                                ' -output_path '.catdir($self->o('output_path'), 'optimise').
                                ' -uniprot_filename '.$self->o('protein_entry_loc').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbport '.$self->o('reference_db','-port').
                                ' -dbhost '.$self->o('reference_db','-host').
                                ' -dbname '.$self->o('reference_db','-dbname').
                                ' -prod_dbuser '.$self->o('user_r').
                                ' -prod_dbhost '.$self->o('production_db','-host').
                                ' -prod_dbname '.$self->o('production_db','-dbname').
                                ' -prod_dbport '.$self->o('production_db','-port').
                                ' -ise -core'
                       },
        -max_retry_count => 0,
        -rc_name => '8GB',
        -flow_into => {
                        1 => ['clean_unused_analyses'],
                      },
      },


      {
        -logic_name => 'clean_unused_analyses',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'DELETE FROM analysis WHERE logic_name IN'.
              ' ("spliced_elsewhere","pseudogenes","genblast","genblast_not_best","project_pseudogene",'.
              ' "project_lincrna","project_transcripts","ig_tr_collapse")',
            'DELETE FROM ad USING analysis_description ad LEFT JOIN analysis a ON ad.analysis_id = a.analysis_id'.
              ' WHERE a.analysis_id IS NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['drop_backup_tables_job'],
                      },
      },


      {
        -logic_name => 'drop_backup_tables_job',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputquery => 'SHOW TABLES LIKE "%bak%"',
          column_names => ['table'],
          db_conn =>$self->o('reference_db'),
        },
        -rc_name    => 'default',
        -flow_into => {
          '2->A' => ['drop_backup_tables'],
          'A->1' => ['final_meta_updates'],
        },
      },


      {
        -logic_name => 'drop_backup_tables',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          sql => 'DROP TABLE #table#',
          db_conn =>$self->o('reference_db'),
        },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'final_meta_updates',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'INSERT INTO meta (species_id, meta_key, meta_value) VALUES '.
              '(1, "genebuild.last_geneset_update", (SELECT CONCAT((EXTRACT(YEAR FROM now())),"-",(LPAD(EXTRACT(MONTH FROM now()),2,"0")))))'
          ],
        },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['final_cleaning'],
                       },
      },


      {
        -logic_name => 'final_cleaning',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'TRUNCATE associated_xref',
            'TRUNCATE dependent_xref',
            'TRUNCATE identity_xref',
            'TRUNCATE object_xref',
            'TRUNCATE ontology_xref',
            'TRUNCATE xref',
            'DELETE exon FROM exon LEFT JOIN exon_transcript ON exon.exon_id = exon_transcript.exon_id WHERE exon_transcript.exon_id IS NULL',
            'DELETE supporting_feature FROM supporting_feature LEFT JOIN exon ON supporting_feature.exon_id = exon.exon_id WHERE exon.exon_id IS NULL',
            'DELETE supporting_feature FROM supporting_feature LEFT JOIN dna_align_feature ON feature_id = dna_align_feature_id WHERE feature_type="dna_align_feature" AND dna_align_feature_id IS NULL',
            'DELETE supporting_feature FROM supporting_feature LEFT JOIN protein_align_feature ON feature_id = protein_align_feature_id WHERE feature_type="protein_align_feature" AND protein_align_feature_id IS NULL',
            'DELETE transcript_supporting_feature FROM transcript_supporting_feature LEFT JOIN dna_align_feature ON feature_id = dna_align_feature_id WHERE feature_type="dna_align_feature" AND dna_align_feature_id IS NULL',
            'DELETE transcript_supporting_feature FROM transcript_supporting_feature LEFT JOIN protein_align_feature ON feature_id = protein_align_feature_id WHERE feature_type="protein_align_feature" AND protein_align_feature_id IS NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['core_gene_set_sanity_checks'],
                       },
      },


      {
        -logic_name => 'core_gene_set_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('reference_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'core'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['core_healthchecks'],
                      },
      },


      {
        -logic_name => 'core_healthchecks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
        -parameters => {
                         input_db         => $self->o('reference_db'),
                         species          => $self->o('species_name'),
                         group            => 'core_handover',
                         enscode_root_dir => $self->o('enscode_root_dir'),
                       },
        -rc_name    => 'default',
        -max_retry_count => 0,
        -flow_into => {
                        '1->A' => ['fan_otherfeatures_db', 'fan_rnaseq_db'],
                        'A->1' => ['update_assembly_name'],},
       },


      {
        -logic_name => 'fan_otherfeatures_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ -n "#assembly_refseq_accession#" ]; then exit 0; else exit 42;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['create_otherfeatures_db'],
                       },
      },


      {
        -logic_name => 'create_otherfeatures_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('cdna_db'),
                         target_db => $self->o('otherfeatures_db'),
                         create_type => 'copy',
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['update_cdna_analyses'],
                       },
      },



      {
        -logic_name => 'update_cdna_analyses',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('otherfeatures_db'),
          sql => [
            'UPDATE gene, analysis SET gene.analysis_id = analysis.analysis_id WHERE analysis.logic_name = "cdna_alignment"',
            'UPDATE transcript join gene using(gene_id) set transcript.analysis_id=gene.analysis_id',
            'UPDATE gene set biotype="cdna"',
            'UPDATE transcript set biotype="cdna"',
            'UPDATE dna_align_feature, analysis SET dna_align_feature.analysis_id = analysis.analysis_id WHERE analysis.logic_name = "cdna_alignment"',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['create_refseq_import_ids_to_copy'],
                        'A->1' => ['update_otherfeatures_db'],
                      },
      },



      {
        -logic_name => 'create_refseq_import_ids_to_copy',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db    => $self->o('refseq_db'),
                         iid_type     => 'feature_id',
                         feature_type => 'gene',
                         batch_size   => 500,
                      },
        -flow_into => {
                       '2' => ['copy_refseq_genes_to_otherfeatures'],
                      },

        -rc_name    => 'default',
      },


      {
        -logic_name => 'copy_refseq_genes_to_otherfeatures',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
        -parameters => {
                         copy_genes_directly => 1,
                         source_db => $self->o('refseq_db'),
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('otherfeatures_db'),
                       },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'update_otherfeatures_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('otherfeatures_db'),
          sql => [
            'DELETE analysis_description FROM analysis_description join analysis using(analysis_id)'.
              ' WHERE logic_name NOT IN ("refseq_import","cdna_alignment")',
            'DELETE FROM analysis WHERE logic_name NOT IN ("refseq_import","cdna_alignment")',
            'DELETE FROM meta WHERE meta_key LIKE "%.level"',
            'DELETE FROM meta WHERE meta_key LIKE "provider.%"',
            'DELETE FROM meta WHERE meta_key LIKE "assembly.web_accession%"',
            'DELETE FROM meta WHERE meta_key LIKE "removed_evidence_flag.%"',
            'DELETE FROM meta WHERE meta_key LIKE "marker.%"',
            'DELETE FROM meta WHERE meta_key = "repeat.analysis"',
            'DELETE FROM meta WHERE meta_key IN'.
              ' ("genebuild.last_geneset_update","genebuild.method","genebuild.projection_source_db","genebuild.start_date")',
            'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"genebuild.last_otherfeatures_update",NOW())',
            'UPDATE transcript JOIN transcript_supporting_feature USING(transcript_id)'.
              ' JOIN dna_align_feature ON feature_id = dna_align_feature_id SET stable_id = hit_name',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['set_otherfeatures_meta_coords'],
                      },
      },


      {
        -logic_name => 'set_otherfeatures_meta_coords',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_coord_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('otherfeatures_db','-host').
                                ' -port '.$self->o('otherfeatures_db','-port').
                                ' -dbpattern '.$self->o('otherfeatures_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_otherfeatures_meta_levels'],
                      },
      },


      {
        -logic_name => 'set_otherfeatures_meta_levels',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_levels_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('otherfeatures_db','-host').
                                ' -port '.$self->o('otherfeatures_db','-port').
                                ' -dbname '.$self->o('otherfeatures_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_otherfeatures_frameshift_introns'] },
      },


      {
        -logic_name => 'set_otherfeatures_frameshift_introns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('frameshift_attrib_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('otherfeatures_db','-host').
                                ' -port '.$self->o('otherfeatures_db','-port').
                                ' -dbpattern '.$self->o('otherfeatures_db','-dbname')
                       },
        -rc_name => '4GB',
        -flow_into => { 1 => ['set_otherfeatures_canonical_transcripts'] },
      },


      {
        -logic_name => 'set_otherfeatures_canonical_transcripts',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('select_canonical_script').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbhost '.$self->o('otherfeatures_db','-host').
                                ' -dbport '.$self->o('otherfeatures_db','-port').
                                ' -dbname '.$self->o('otherfeatures_db','-dbname').
                                ' -dnadbuser '.$self->o('user_r').
                                ' -dnadbhost '.$self->o('dna_db','-host').
                                ' -dnadbport '.$self->o('dna_db','-port').
                                ' -dnadbname '.$self->o('dna_db','-dbname').
                                ' -coord toplevel -write'
                       },
        -rc_name => '2GB',
        -flow_into => { 1 => ['populate_production_tables_otherfeatures'] },
      },


      {
        -logic_name => 'populate_production_tables_otherfeatures',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
        -parameters => {
                         'target_db'        => $self->o('otherfeatures_db'),
                         'output_path'      => $self->o('output_path'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'production_db'    => $self->o('production_db'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['null_otherfeatures_columns'],
                       },
      },


      {
        -logic_name => 'null_otherfeatures_columns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('otherfeatures_db'),
          sql => [
            'UPDATE dna_align_feature SET external_db_id = NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['load_external_db_ids_and_optimise_otherfeatures'],
                      },
      },


      {
        -logic_name => 'load_external_db_ids_and_optimise_otherfeatures',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('load_optimise_script').
                                ' -output_path '.catdir($self->o('output_path'), 'optimise_otherfeatures').
                                ' -uniprot_filename '.$self->o('protein_entry_loc').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbport '.$self->o('otherfeatures_db','-port').
                                ' -dbhost '.$self->o('otherfeatures_db','-host').
                                ' -dbname '.$self->o('otherfeatures_db','-dbname').
                                ' -prod_dbuser '.$self->o('user_r').
                                ' -prod_dbhost '.$self->o('production_db','-host').
                                ' -prod_dbname '.$self->o('production_db','-dbname').
                                ' -prod_dbport '.$self->o('production_db','-port').
                                ' -verbose'
                       },
        -max_retry_count => 0,
        -rc_name => '4GB',
        -flow_into => {
                        1 => ['otherfeatures_sanity_checks'],
                      },
      },


      {
        -logic_name => 'otherfeatures_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('otherfeatures_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{'otherfeatures'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['otherfeatures_healthchecks'],
                      },
      },


      {
        -logic_name => 'otherfeatures_healthchecks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
        -parameters => {
          input_db         => $self->o('otherfeatures_db'),
          species          => $self->o('species_name'),
          group            => 'otherfeatures_handover',
        },
        -max_retry_count => 0,

        -rc_name    => '4GB',
      },


      {
        -logic_name => 'fan_rnaseq_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_rnaseq#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['create_rnaseq_db'],
                       },

      },


      {
        -logic_name => 'create_rnaseq_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('rnaseq_blast_db'),
                         target_db => $self->o('rnaseq_db'),
                         create_type => 'copy',
                       },
        -rc_name    => 'default',

        -flow_into => {
                        '1' => ['prepare_rnaseq_meta_data'],
                      },
      },


      {
        -logic_name => 'prepare_rnaseq_meta_data',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
                         db_conn => $self->o('rnaseq_db'),
                         sql => [
                          'TRUNCATE dna_align_feature',
                          'DELETE FROM transcript_supporting_feature WHERE feature_type = "dna_align_feature"',
                          'DELETE FROM supporting_feature WHERE feature_type = "dna_align_feature"',
                          'DELETE FROM analysis WHERE logic_name NOT LIKE "%rnaseq%"',
                          'INSERT INTO analysis (logic_name, module, created, db_version) VALUES ("other_protein", "HiveBlastRNAseq", NOW(), "#uniprot_version#")',
                          'UPDATE protein_align_feature paf, analysis a SET paf.analysis_id = a.analysis_id WHERE a.logic_name = "other_protein"',
                          'DELETE FROM meta WHERE meta_key LIKE "provider\.%"',
                          'DELETE FROM meta WHERE meta_key LIKE "assembly.web_accession%"',
                          'DELETE FROM meta WHERE meta_key LIKE "removed_evidence_flag\.%"',
                          'DELETE FROM meta WHERE meta_key LIKE "marker\.%"',
                          'DELETE FROM meta WHERE meta_key IN ("genebuild.method","genebuild.projection_source_db","genebuild.start_date","repeat.analysis")',
                          'UPDATE gene JOIN transcript USING(gene_id) SET canonical_transcript_id = transcript_id',
                          'UPDATE transcript JOIN translation USING(transcript_id) SET canonical_translation_id = translation_id',
                          'UPDATE intron_supporting_evidence ise, analysis a1, analysis a2 SET ise.analysis_id = a2.analysis_id WHERE ise.analysis_id = a1.analysis_id AND a2.logic_name = REPLACE(a1.logic_name, "rnaseq_gene", "rnaseq_ise")',
                          'UPDATE gene SET source = "ensembl", biotype = "protein_coding", stable_id = NULL',
                          'UPDATE transcript SET source = "ensembl", biotype = "protein_coding", stable_id = NULL',
                          'UPDATE translation SET stable_id = NULL',
                          'UPDATE exon SET stable_id = NULL',
                         ],
                         uniprot_version => $self->o('uniprot_version'),
                       },
        -rc_name    => 'default',

        -flow_into => {
                        '1' => ['generate_rnaseq_stable_ids'],
                      },
      },


      {
        -logic_name => 'generate_rnaseq_stable_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
        -parameters => {
                         enscode_root_dir => $self->o('enscode_root_dir'),
                         mapping_required => 0,
                         target_db => $self->o('rnaseq_db'),
                         id_start => 'RNASEQ',
                         output_path => $self->o('output_path'),
                         _stable_id_file => 'rnaseq_stable_ids.sql',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['populate_production_tables_rnaseq'],
                      },
      },


      {
        -logic_name => 'populate_production_tables_rnaseq',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
        -parameters => {
                         'target_db'        => $self->o('rnaseq_db'),
                         'output_path'      => $self->o('output_path'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'production_db'    => $self->o('production_db'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['dump_daf_introns'],
                       },
      },


      {
        -logic_name => 'dump_daf_introns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn => $self->o('rnaseq_refine_db'),
                         input_query => 'SELECT daf.* FROM dna_align_feature daf, analysis a WHERE daf.analysis_id = a.analysis_id AND a.logic_name != "rough_transcripts"',
                         command_out => q(sort -nk2 -nk3 -nk4 | sed 's/NULL/\\N/g;s/^[0-9]\+/\\N/' > #daf_file#),
                         daf_file => $self->o('rnaseq_daf_introns_file'),
                         prepend => ['-NB', '-q'],
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['load_daf_introns'],
                      },
      },


      {
        -logic_name => 'load_daf_introns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn => $self->o('rnaseq_db'),
                         input_query => 'LOAD DATA LOCAL INFILE "#daf_file#" INTO TABLE dna_align_feature',
                         daf_file => $self->o('rnaseq_daf_introns_file'),
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_rnaseq_meta_coords'],
                      },
      },


      {
        -logic_name => 'set_rnaseq_meta_coords',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_coord_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('rnaseq_db','-host').
                                ' -port '.$self->o('rnaseq_db','-port').
                                ' -dbpattern '.$self->o('rnaseq_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_rnaseq_meta_levels'],
                      },
      },


      {
        -logic_name => 'set_rnaseq_meta_levels',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_levels_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('rnaseq_db','-host').
                                ' -port '.$self->o('rnaseq_db','-port').
                                ' -dbname '.$self->o('rnaseq_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['optimise_rnaseq'] },
      },


      {
        -logic_name => 'optimise_rnaseq',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('load_optimise_script').
                                ' -output_path '.catfile($self->o('rnaseq_dir'), 'optimise_rnaseq').
                                ' -uniprot_filename '.$self->o('protein_entry_loc').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbport '.$self->o('rnaseq_db','-port').
                                ' -dbhost '.$self->o('rnaseq_db','-host').
                                ' -dbname '.$self->o('rnaseq_db','-dbname').
                                ' -prod_dbuser '.$self->o('user_r').
                                ' -prod_dbhost '.$self->o('production_db','-host').
                                ' -prod_dbname '.$self->o('production_db','-dbname').
                                ' -prod_dbport '.$self->o('production_db','-port').
                                ' -nodaf -ise'
                       },
        -max_retry_count => 0,
        -rc_name => '8GB',
        -flow_into => {
                        1 => ['rnaseq_gene_sanity_checks'],
                      },
      },


      {
        -logic_name => 'rnaseq_gene_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('rnaseq_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{'rnaseq_final'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['create_bam_file_job'],
                      },
      },
      {
        -logic_name => 'create_bam_file_job',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputcmd => 'cd #working_dir#; ls *.bam',
          column_names => ['bam_file'],
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => 'default',
        -flow_into  => {
          '2->A' => ['create_chromosome_file'],
          'A->1' => ['concat_md5_sum'],
        },
      },

      {
        -logic_name => 'create_chromosome_file',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#samtools# view -H '.catfile('#working_dir#', '#bam_file#').q( | grep \@SQ |cut -f2,3 | sed 's/[SL]N://g' > ).catfile('#working_dir#', '#bam_file#.txt'),
          samtools => $self->o('samtools_path'),
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => '3GB',
        -flow_into  => {
          1 => ['bam2bedgraph'],
        },
      },

      {
        -logic_name => 'bam2bedgraph',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#bedtools# genomecov -ibam '.catfile('#working_dir#', '#bam_file#').' -bg -split > '.catfile('#working_dir#', '#bam_file#.bg.unsorted').' ; LC_COLLATE=C sort -k1,1 -k2,2n '.catfile('#working_dir#', '#bam_file#.bg.unsorted').' > '.catfile('#working_dir#', '#bam_file#.bg'),
          bedtools => $self->o('bedtools'),
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => '3GB',
        -flow_into  => {
          1 => ['bedgrap2bigwig'],
          -1 => ['bam2bedgraph_himem'],
        },
      },

      {
        -logic_name => 'bam2bedgraph_himem',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#bedtools# genomecov -ibam '.catfile('#working_dir#', '#bam_file#').' -bg -split > '.catfile('#working_dir#', '#bam_file#.bg.unsorted').' ; LC_COLLATE=C sort -k1,1 -k2,2n '.catfile('#working_dir#', '#bam_file#.bg.unsorted').' > '.catfile('#working_dir#', '#bam_file#.bg'),
          bedtools => $self->o('bedtools'),
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => '40GB',
        -flow_into  => {
          1 => ['bedgrap2bigwig_himem'],
        },
      },

      {
        -logic_name => 'bedgrap2bigwig',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#bedGraphToBigWig# '.catfile('#working_dir#', '#bam_file#.bg').' '.catfile('#working_dir#', '#bam_file#.txt').' '.catfile('#working_dir#', '#bam_file#.bw'),
          bedGraphToBigWig => $self->o('bedGraphToBigWig'),
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => '3GB',
        -flow_into  => {
          1 => ['clean_bg_files'],
          -1 => ['bedgrap2bigwig_himem'],
        },
      },

      {
        -logic_name => 'bedgrap2bigwig_himem',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#bedGraphToBigWig# '.catfile('#working_dir#', '#bam_file#.bg').' '.catfile('#working_dir#', '#bam_file#.txt').' '.catfile('#working_dir#', '#bam_file#.bw'),
          bedGraphToBigWig => $self->o('bedGraphToBigWig'),
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => '8GB',
        -flow_into  => {
          1 => ['clean_bg_files'],
        },
      },

      {
        -logic_name => 'clean_bg_files',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'rm  '.catfile('#working_dir#', '#bam_file#.bg.unsorted').' '.catfile('#working_dir#', '#bam_file#.bg').' '.catfile('#working_dir#', '#bam_file#.txt'),
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['md5_sum'],
        },
      },

      {
        -logic_name => 'md5_sum',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cd #working_dir#;md5sum #bam_file#.bw #bam_file# #bam_file#.bai > #bam_file#.md5',
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => 'default',
      },

      {
        -logic_name => 'concat_md5_sum',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cat '.catfile('#working_dir#', '*md5').' > '.catfile('#working_dir#', 'md5sum.txt.1'),
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['clean_concat_md5_sum'],
        },
      },

      {
        -logic_name => 'clean_concat_md5_sum',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'rm '.catfile('#working_dir#', '*md5'),
          working_dir => $self->o('merge_dir'),
        },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['create_readme'],
        },
      },

      {
        -logic_name => 'create_readme',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cd #working_dir#;FILES=($(ls *.bam));printf "#free_text#" | sed "s/NUM/$((${#FILES[*]}-1))/g;s/ \([a-z]\)\([a-z]\+_\)/ \U\1\E\2/;s/_/ /g" > README.1; IFS=$\'\n\';echo "${FILES[*]}" >> README.1',
          working_dir => $self->o('merge_dir'),
          species_name  => $self->o('species_name'),
          free_text => 'Note\n------\n\n'.
                       'The RNASeq data for #species_name# consists of NUM individual samples and one merged set containing all NUM samples.\n\n'.
                       'All files have an index file (.bai) and a BigWig file (.bw) which contains the coverage information.\n\n'.
                       'Use the md5sum.txt file to check the integrity of the downloaded files.\n\n'.
                       'Files\n-----\n',
        },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['rnaseq_healthchecks'],
        },
      },


      {
        -logic_name => 'rnaseq_healthchecks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
        -parameters => {
          input_db => $self->o('rnaseq_db'),
          species  => $self->o('species_name'),
          group    => 'rnaseq_handover',
        },
        -max_retry_count => 0,

        -rc_name    => '4GB',
      },


      {
        -logic_name => 'update_assembly_name',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => { 
                         cmd => 'perl '.$self->o('assembly_name_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbname '.$self->o('reference_db','-dbname').
                                ' -driver '.$self->o('hive_driver').
                                ' -assembly_accession '.$self->o('assembly_accession').
                                ' -assembly_name '.$self->o('assembly_name').
                                ' -registry_host '.$self->o('registry_host').
                                ' -registry_port '.$self->o('registry_port').
                                ' -registry_db '.$self->o('registry_db'),
                       },
        -rc_name => 'default',
       },

    ];
}


sub resource_classes {
  my $self = shift;

  return {
    '1GB' => { LSF => $self->lsf_resource_builder('production-rh74', 1000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '2GB_lastz' => { LSF => [$self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '2GB' => { LSF => $self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '3GB' => { LSF => $self->lsf_resource_builder('production-rh74', 3000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB_lastz' => { LSF => [$self->lsf_resource_builder('production-rh74', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]}, 
    '4GB' => { LSF => $self->lsf_resource_builder('production-rh74', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB' => { LSF => $self->lsf_resource_builder('production-rh74', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '6GB' => { LSF => $self->lsf_resource_builder('production-rh74', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '6GB_registry' => { LSF => [$self->lsf_resource_builder('production-rh74', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{registry_file}]},
    '7GB' => { LSF => $self->lsf_resource_builder('production-rh74', 7000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '8GB_lastz' => { LSF => [$self->lsf_resource_builder('production-rh74', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '8GB' => { LSF => $self->lsf_resource_builder('production-rh74', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '9GB' => { LSF => $self->lsf_resource_builder('production-rh74', 9000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '10GB' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '15GB_lastz' => { LSF => [$self->lsf_resource_builder('production-rh74', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '15GB' => { LSF => $self->lsf_resource_builder('production-rh74', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '20GB' => { LSF => $self->lsf_resource_builder('production-rh74', 20000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '25GB' => { LSF => $self->lsf_resource_builder('production-rh74', 25000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '30GB' => { LSF => $self->lsf_resource_builder('production-rh74', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '35GB' => { LSF => $self->lsf_resource_builder('production-rh74', 35000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '40GB' => { LSF => $self->lsf_resource_builder('production-rh74', 40000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '50GB' => { LSF => $self->lsf_resource_builder('production-rh74', 50000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '75GB' => { LSF => $self->lsf_resource_builder('production-rh74', 75000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '80GB' => { LSF => $self->lsf_resource_builder('production-rh74', 80000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '100GB' => { LSF => $self->lsf_resource_builder('production-rh74', 100000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default' => { LSF => $self->lsf_resource_builder('production-rh74', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker' => { LSF => $self->lsf_resource_builder('production-rh74', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker_rebatch' => { LSF => $self->lsf_resource_builder('production-rh74', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'simple_features' => { LSF => $self->lsf_resource_builder('production-rh74', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan' => { LSF => $self->lsf_resource_builder('production-rh74', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan_short' => { LSF => $self->lsf_resource_builder('production-rh74', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'blast' => { LSF => $self->lsf_resource_builder('production-rh74', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast10GB' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast_retry' => { LSF => $self->lsf_resource_builder('production-rh74', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'genblast' => { LSF => $self->lsf_resource_builder('production-rh74', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genblast_retry' => { LSF => $self->lsf_resource_builder('production-rh74', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'project_transcripts' => { LSF => $self->lsf_resource_builder('production-rh74', 7200, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'projection_db_server'}, $self->default_options->{'projection_lastz_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'refseq_import' => { LSF => $self->lsf_resource_builder('production-rh74', 9900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'refseq_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'layer_annotation' => { LSF => $self->lsf_resource_builder('production-rh74', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genebuilder' => { LSF => $self->lsf_resource_builder('production-rh74', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'transcript_finalisation' => { LSF => $self->lsf_resource_builder('production-rh74', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'filter' => { LSF => $self->lsf_resource_builder('production-rh74', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '2GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}], undef, $self->default_options->{'use_threads'})},
    '3GB_merged_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 3000, [$self->default_options->{'pipe_db_server'}], undef, $self->default_options->{'rnaseq_merge_threads'})},
    '5GB_merged_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 5000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'rnaseq_merge_threads'}))},
    '5GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 5000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
    '10GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
    '20GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 20000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
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
