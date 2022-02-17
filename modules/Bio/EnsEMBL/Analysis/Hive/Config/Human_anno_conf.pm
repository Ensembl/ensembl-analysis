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

package Human_anno_conf;

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
    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'user'                      => '', # write db user
    'password'                  => '', # password for write db user
    'base_output_dir'           => '', # Where to write the files to
    'registry_file'             => '', # This needs to be a standard registry with the production/meta data/taxonomy db adaptors in there
    'pipeline_name'             => '' || $self->o('production_name').$self->o('production_name_modifier').'_'.$self->o('ensembl_release'),
    'release_number'            => '' || $self->o('ensembl_release'),
    'xy_scanner_path'           => '/hps/software/users/ensembl/repositories/fergal/ensembl-analysis/scripts/genebuild/xy_scanner.py',
    'x_marker_fasta_path'       => '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/hprc/x_markers.fa',
    'y_marker_fasta_path'       => '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/hprc/y_markers.fa',
    'ref_db_server'             => 'mysql-ens-mirror-1', # host for dna db
    'ref_db_port'               => '4240',
    'ref_db_name'               => 'homo_sapiens_core_104_38',
    'user_r'                    => 'ensro', # read only db user
    'num_threads' => 20,
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
    'output_path'               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'             => '', # Name (as it appears in the assembly report file)
    'assembly_accession'        => '', # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'assembly_refseq_accession' => '', # Versioned GCF accession, e.g. GCF_001857705.1
    'stable_id_prefix'          => '', # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    'use_genome_flatfile'       => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    'species_url'               => '' || $self->o('production_name').$self->o('production_name_modifier'), # sets species.url meta key
    'species_division'          => '', # sets species.division meta key
    'stable_id_start'           => '', # When mapping is not required this is usually set to 0
    'mapping_required'          => '0', # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'                => '', # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now
    'uniprot_version'           => 'uniprot_2018_07', # What UniProt data dir to use for various analyses
    'vertrna_version'           => '136', # The version of VertRNA to use, should correspond to a numbered dir in VertRNA dir
    'production_name_modifier'  => '', # Do not set unless working with non-reference strains, breeds etc. Must include _ in modifier, e.g. _hni for medaka strain HNI
    'compara_registry_file'     => '',

    # Keys for custom loading, only set/modify if that's what you're doing
    'repeatmodeler_library'        => '', # This should be the path to a custom repeat library, leave blank if none exists
    'use_repeatmodeler_to_mask'    => '0', # Setting this will include the repeatmodeler library in the masking process
    'protein_blast_db'             => '' || catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata'), # Blast database for comparing the final models to.
    'protein_blast_index'          => '' || catdir($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata_index'), # Indicate Index for the blast database.
    'protein_entry_loc'            => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc'), # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

    'softmask_logic_names' => [],

########################
# Pipe and ref db info
########################

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is
    'projection_lastz_db_server'   => $self->o('pipe_db_server'),

    'provider_name'                => 'Ensembl',
    'provider_url'                 => 'www.ensembl.org',

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_pipe_'.$self->o('release_number'),
    'dna_db_name'                   => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_core_'.$self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    'production_db_server'         => 'mysql-ens-meta-prod-1',
    'production_db_port'           => '4483',

    databases_to_delete => ['reference_db'],

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
    load_optimise_script              => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl'),

    ensembl_misc_script        => catdir($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts'),
    repeat_types_script        => catfile($self->o('ensembl_misc_script'), 'repeats', 'repeat-types.pl'),
    meta_coord_script          => catfile($self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl'),
    meta_levels_script         => catfile($self->o('ensembl_misc_script'), 'meta_levels.pl'),
    frameshift_attrib_script   => catfile($self->o('ensembl_misc_script'), 'frameshift_transcript_attribs.pl'),
    select_canonical_script    => catfile($self->o('ensembl_misc_script'),'canonical_transcripts', 'select_canonical_transcripts.pl'),
    assembly_name_script       => catfile($self->o('ensembl_analysis_script'), 'update_assembly_name.pl'),

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
    'minimap2_path'          => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/minimap2',
    'paftools_path'          => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/paftools.js',
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
    samtools_path => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary

    'uniprot_genblast_batch_size' => 15,
    'uniprot_table_name'          => 'uniprot_sequences',

    'genblast_path'     => catfile($self->o('binary_base'), 'genblast'),
    'genblast_eval'     => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
    'genblast_cov'      => '0.5',
    'genblast_pid'      => '30',
    'genblast_max_rank' => '5',
    'genblast_flag_small_introns' => 1,
    'genblast_flag_subpar_models' => 0,

# Best targetted stuff
    cdna_table_name    => 'cdna_sequences',

# RNA-seq pipeline stuff
    # You have the choice between:
    #  * using a csv file you already created
    #  * using a study_accession like PRJEB19386
    #  * using the taxon_id of your species
    # 'rnaseq_summary_file' should always be set. If 'taxon_id' or 'study_accession' are not undef
    # they will be used to retrieve the information from ENA and to create the csv file. In this case,
    # 'file_columns' and 'summary_file_delimiter' should not be changed unless you know what you are doing
    'summary_csv_table' => 'csv_data',
    'read_length_table' => 'read_length',

    use_threads => 3,
    rnaseq_merge_threads => 12,

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

# lincRNA pipeline stuff
    'lncrna_dir' => catdir($self->o('output_path'), 'lincrna'),
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
    'realign_table_name'                    => 'projection_source_sequences',

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

    'reference_db' => {
      -dbname => $self->o('ref_db_name'),
      -host   => $self->o('ref_db_server'),
      -port   => $self->o('ref_db_port'),
      -user   => $self->o('user_r'),
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
        -analysis_capacity => 20,
      },


      {
        # Load some meta info and seq_region_synonyms
        -logic_name => 'load_meta_info',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => '#core_db#',
          sql => [
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
        -analysis_capacity => 20,
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
          1 => ['create_faidx'],
        },
     },


     {
        -logic_name => 'create_faidx',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '5GB',
        -parameters => {
          cmd => 'if [ ! -e "'.'#reheadered_toplevel_genome_file#'.'.fai" ]; then '.$self->o('samtools_path').' faidx '.'#reheadered_toplevel_genome_file#'.';fi',
        },

       -flow_into  => {
          1 => ['create_minimap2_index'],
        },

      },


      {
        -logic_name => 'create_minimap2_index',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ ! -e "'.'#reheadered_toplevel_genome_file#.mmi'.'" ]; then '.$self->o('minimap2_path').
                 ' -d '.'#reheadered_toplevel_genome_file#.mmi'.' '.'#reheadered_toplevel_genome_file#'.';fi',
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
                         cmd => 'if [ -s "'.'#reheadered_toplevel_genome_file#.mmi'.'" ]; then exit 0; else exit 42;fi',
                         return_codes_2_branches => {'42' => 2},
        },
        -flow_into  => {
          1 => ['xy_scanner'],
        },
        -rc_name => 'default',
      },


      {
        -logic_name => 'xy_scanner',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::XYScanner',
        -parameters => {
                         xy_scanner_path     => $self->o('xy_scanner_path'),
                         target_genome_index => '#reheadered_toplevel_genome_file#'.'.mmi',
                         x_marker_fasta_path => $self->o('x_marker_fasta_path'),
                         y_marker_fasta_path => $self->o('y_marker_fasta_path'),
                         output_dir          => '#output_path#',
                       },
        -rc_name    => '30GB',
        -flow_into  => {
          '1->A' => ['create_remap_jobs'],
          'A->1' => ['create_paralogue_jobs'],
        },
      },


      {
        -logic_name => 'create_remap_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          target_db           => $self->o('reference_db'),
          iid_type            => 'feature_id',
          feature_type        => 'gene',
          batch_size          => 100,
        },
        -rc_name    => '4GB',
        -flow_into => {
          2 => {'minimap2genomic' => {'xy_scanner' => '#xy_scanner#', 'core_db' => '#core_db#','genome_index' => '#reheadered_toplevel_genome_file#'.'.mmi','iid' => '#iid#'}},
        },
      },


      {
        -logic_name => 'minimap2genomic',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2Remap',
        -parameters => {
                         genome_index   => '#genome_index#',
                         source_dna_db => $self->o('reference_db'),
                         source_gene_db => $self->o('reference_db'),
                         target_dna_db => '#core_db#',
                         target_gene_db => '#core_db#',
                         paftools_path => $self->o('paftools_path'),
                         minimap2_path => $self->o('minimap2_path'),
                       },
        -rc_name    => '25GB',
      },


      {
        -logic_name => 'create_paralogue_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          target_db           => '#core_db#',
          iid_type            => 'feature_id',
          feature_type        => 'gene',
          batch_size          => 1000,
        },
        -rc_name    => '4GB',
        -flow_into => {
          '2->A' => {'find_paralogues' => {'core_db' => '#core_db#','genome_file' => '#reheadered_toplevel_genome_file#','genome_index' => '#reheadered_toplevel_genome_file#'.'.mmi','iid' => '#iid#'}},
          'A->1' => ['collapse_paralogues'],
        },
      },


      {
        -logic_name => 'find_paralogues',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::FindRecentParalogues',
        -parameters => {
                         genome_index    => '#genome_index#',
                         genome_file     => '#reheadered_toplevel_genome_file#',
                         target_dna_db   => '#core_db#',
                         target_db  => '#core_db#',
                         paftools_path => $self->o('paftools_path'),
                         minimap2_path => $self->o('minimap2_path'),
                       },
        -rc_name    => '25GB',
      },


      {
        -logic_name => 'collapse_paralogues',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::CollapseParalogues',
        -parameters => {
                         genome_file     => '#reheadered_toplevel_genome_file#',
                         target_dna_db   => '#core_db#',
                         target_db  => '#core_db#',
                       },
        -flow_into => {
                        1 => ['update_biotypes_and_analyses'],
                      },

        -rc_name    => '25GB',
        -flow_into  => {
          1 => ['finalise_geneset'],
        },
      },


      {
        -logic_name => 'finalise_geneset',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::FinaliseRemappedGeneset',
        -parameters => {
                         genome_file     => '#reheadered_toplevel_genome_file#',
                         source_gene_db => $self->o('reference_db'),
                         target_dna_db   => '#core_db#',
                         target_db  => '#core_db#',
                       },
        -flow_into => {
                        1 => ['update_biotypes_and_analyses'],
                      },

        -rc_name    => '25GB',
        -flow_into  => {
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
        -flow_into => { 1 => ['set_frameshift_introns'] },
      },


      {
        -logic_name => 'set_frameshift_introns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('frameshift_attrib_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('core_db','-host').
                                ' -port '.$self->o('core_db','-port').
                                ' -dbpattern '.'#core_dbname#'
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
                                ' -dbhost '.$self->o('core_db','-host').
                                ' -dbport '.$self->o('core_db','-port').
                                ' -dbname '.'#core_dbname#'.
                                ' -coord toplevel -write'
                       },
        -rc_name => '10GB',
        -flow_into => { 1 => ['null_columns'] },
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
                        1 => ['run_stable_ids'],
                      },
      },


      {
        -logic_name => 'run_stable_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
        -parameters => {
                         enscode_root_dir => $self->o('enscode_root_dir'),
                         mapping_required => 0,
                         target_db => '#core_db#',
                         id_start => '#stable_id_prefix#'.'#stable_id_start#',
                         output_path => '#output_path#',
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
                         1 => ['final_cleaning'],
                       },
      },


      {
        -logic_name => 'final_cleaning',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => '#core_db#',
          sql => [
            'DELETE exon FROM exon LEFT JOIN exon_transcript ON exon.exon_id = exon_transcript.exon_id WHERE exon_transcript.exon_id IS NULL',
            'DELETE supporting_feature FROM supporting_feature LEFT JOIN exon ON supporting_feature.exon_id = exon.exon_id WHERE exon.exon_id IS NULL',
            'DELETE supporting_feature FROM supporting_feature LEFT JOIN dna_align_feature ON feature_id = dna_align_feature_id WHERE feature_type="dna_align_feature" AND dna_align_feature_id IS NULL',
            'DELETE supporting_feature FROM supporting_feature LEFT JOIN protein_align_feature ON feature_id = protein_align_feature_id WHERE feature_type="protein_align_feature" AND protein_align_feature_id IS NULL',
            'DELETE transcript_supporting_feature FROM transcript_supporting_feature LEFT JOIN dna_align_feature ON feature_id = dna_align_feature_id WHERE feature_type="dna_align_feature" AND dna_align_feature_id IS NULL',
            'DELETE transcript_supporting_feature FROM transcript_supporting_feature LEFT JOIN protein_align_feature ON feature_id = protein_align_feature_id WHERE feature_type="protein_align_feature" AND protein_align_feature_id IS NULL',
            'UPDATE analysis SET logic_name="ensembl" WHERE logic_name="minimap2remap"',
            'UPDATE gene SET analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = "ensembl")'.
            ' WHERE analysis_id IN'.
            ' (SELECT analysis_id FROM analysis WHERE logic_name IN ("find_paralogues","collapse_paralogues"))',
            'DELETE FROM analysis WHERE logic_name IN ("find_paralogues","collapse_paralogues")',
            'DELETE FROM ad USING analysis_description ad LEFT JOIN analysis a ON ad.analysis_id = a.analysis_id WHERE a.analysis_id IS NULL',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.analysis_id = gene.analysis_id',
            'UPDATE repeat_feature SET repeat_start = 1 WHERE repeat_start < 1',
            'UPDATE repeat_feature SET repeat_end = 1 WHERE repeat_end < 1',
            'UPDATE repeat_feature JOIN seq_region USING(seq_region_id) SET repeat_end = length WHERE repeat_end > length',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                       1 => ['add_placeholder_sample_location'],
        },

      },


      {
        -logic_name => 'add_placeholder_sample_location',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddPlaceholderLocation',
        -parameters => {
                        input_db => '#core_db#',
                       },
        -rc_name    => 'default',
        -flow_into => {
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
    'default_registry' => { LSF => [$self->lsf_resource_builder('production', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{'registry_file'}]},
    'gbiab' => { LSF => $self->lsf_resource_builder('production', 50000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, $self->default_options->{'num_threads'})},
    '1GB' => { LSF => $self->lsf_resource_builder('production', 1000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '2GB_lastz' => { LSF => [$self->lsf_resource_builder('production', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '2GB' => { LSF => $self->lsf_resource_builder('production', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '3GB' => { LSF => $self->lsf_resource_builder('production', 3000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB_lastz' => { LSF => [$self->lsf_resource_builder('production', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '4GB' => { LSF => $self->lsf_resource_builder('production', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB' => { LSF => $self->lsf_resource_builder('production', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '6GB' => { LSF => $self->lsf_resource_builder('production', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '6GB_registry' => { LSF => [$self->lsf_resource_builder('production', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{registry_file}]},
    '7GB' => { LSF => $self->lsf_resource_builder('production', 7000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '8GB_lastz' => { LSF => [$self->lsf_resource_builder('production', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '8GB' => { LSF => $self->lsf_resource_builder('production', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '9GB' => { LSF => $self->lsf_resource_builder('production', 9000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '10GB' => { LSF => $self->lsf_resource_builder('production', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '15GB_lastz' => { LSF => [$self->lsf_resource_builder('production', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '15GB' => { LSF => $self->lsf_resource_builder('production', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '20GB' => { LSF => $self->lsf_resource_builder('production', 20000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '25GB' => { LSF => $self->lsf_resource_builder('production', 25000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '30GB' => { LSF => $self->lsf_resource_builder('production', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '35GB' => { LSF => $self->lsf_resource_builder('production', 35000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '40GB' => { LSF => $self->lsf_resource_builder('production', 40000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '50GB' => { LSF => $self->lsf_resource_builder('production', 50000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '75GB' => { LSF => $self->lsf_resource_builder('production', 75000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '80GB' => { LSF => $self->lsf_resource_builder('production', 80000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '100GB' => { LSF => $self->lsf_resource_builder('production', 100000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default' => { LSF => $self->lsf_resource_builder('production', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker' => { LSF => $self->lsf_resource_builder('production', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker_rebatch' => { LSF => $self->lsf_resource_builder('production', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'simple_features' => { LSF => $self->lsf_resource_builder('production', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan' => { LSF => $self->lsf_resource_builder('production', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan_short' => { LSF => $self->lsf_resource_builder('production', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'blast' => { LSF => $self->lsf_resource_builder('production', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast10GB' => { LSF => $self->lsf_resource_builder('production', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast_retry' => { LSF => $self->lsf_resource_builder('production', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'genblast' => { LSF => $self->lsf_resource_builder('production', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genblast_retry' => { LSF => $self->lsf_resource_builder('production', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'project_transcripts' => { LSF => $self->lsf_resource_builder('production', 7200, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'projection_db_server'}, $self->default_options->{'projection_lastz_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'refseq_import' => { LSF => $self->lsf_resource_builder('production', 9900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'refseq_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'layer_annotation' => { LSF => $self->lsf_resource_builder('production', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genebuilder' => { LSF => $self->lsf_resource_builder('production', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'transcript_finalisation' => { LSF => $self->lsf_resource_builder('production', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'filter' => { LSF => $self->lsf_resource_builder('production', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '2GB_multithread' => { LSF => $self->lsf_resource_builder('production', 2000, [$self->default_options->{'pipe_db_server'}], undef, $self->default_options->{'use_threads'})},
    '3GB_merged_multithread' => { LSF => $self->lsf_resource_builder('production', 3000, [$self->default_options->{'pipe_db_server'}], undef, $self->default_options->{'rnaseq_merge_threads'})},
    '5GB_merged_multithread' => { LSF => $self->lsf_resource_builder('production', 5000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'rnaseq_merge_threads'}))},
    '5GB_multithread' => { LSF => $self->lsf_resource_builder('production', 5000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
    '10GB_multithread' => { LSF => $self->lsf_resource_builder('production', 10000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
    '20GB_multithread' => { LSF => $self->lsf_resource_builder('production', 20000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
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
