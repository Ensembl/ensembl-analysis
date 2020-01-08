=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::Config::BestTargetted_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::BestTargetted_conf;

use strict;
use warnings;


use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;

  return {
    %{$self->SUPER::default_options()},
    species_name => '',
    production_name => '',
    taxon_id => 0,

    output_path => '',

    repeatmasker_library => '',

    dna_db_name   => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),
    dna_db_server => '',
    dna_db_port   => 3306,

    datas_db_server => '',
    datas_db_port   => 3306,

###
#  The values below may need to be changed
###
    genewise_db_server => $self->o('datas_db_server'),
    genewise_db_port   => $self->o('datas_db_port'),

    cdna_db_server => $self->o('datas_db_server'),
    cdna_db_port   => $self->o('datas_db_port'),

    cdna2genome_db_server => $self->o('datas_db_server'),
    cdna2genome_db_port   => $self->o('datas_db_port'),

    indicate_path  => catfile($self->o('binary_base'), 'indicate'),
    pmatch_path  => catfile($self->o('binary_base'), 'pmatch'),
    exonerate_path => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
    exonerate_annotation => catfile($self->o('binary_base'), 'exonerate'),
    prepare_cdnas_script => catfile($self->o('enscode_root_dir'), 'ensembl-analysis','scripts', 'genebuild', 'prepare_cdnas.pl'),

    repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library'), 'dust'],

    targetted_path => $self->o('output_path'),
    genome_file    => catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name').'_softmasked_toplevel.fa'),
    cdna_file      => catfile($self->o('output_path'), 'cdnas'),
    annotation_file => $self->o('cdna_file').'.annotation',
###
#  All data below should not need to be changed
###
    killlist_db_server => 'mysql-ens-genebuild-prod-6',
    killlist_db_port   => '4532',

    exonerate_logic_name => 'exonerate',
    ncbi_query => 'txid'.$self->o('taxon_id').'[Organism:noexp]+AND+biomol_mrna[PROP]',
    uniprot_table_name => 'uniprot_sequences',
    cdna_table_name    => 'cdna_sequences',
    refseq_cdna_calculate_coverage_and_pid => 0,
    exonerate_protein_pid => 95,
    exonerate_protein_cov => 50,
    region_padding => 2000,

    databases_to_delete => ['genewise_db', 'cdna_db', 'cdna2genome_db'],

    genewise_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_genewise_'.$self->o('release_number'),
      -host   => $self->o('genewise_db_server'),
      -port   => $self->o('genewise_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    killlist_db => {
      -dbname => $self->o('killlist_db_name'),
      -host   => $self->o('killlist_db_server'),
      -port   => $self->o('killlist_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    cdna_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_cdna_'.$self->o('release_number'),
      -host   => $self->o('cdna_db_server'),
      -port   => $self->o('cdna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    cdna2genome_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_cdna2genome_'.$self->o('release_number'),
      -host   => $self->o('cdna2genome_db_server'),
      -port   => $self->o('cdna2genome_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

  }
}


sub pipeline_analyses {
  my ($self) = @_;

  return [
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
        '1->A' => ['create_genewise_db', 'download_mRNA', 'download_selenocysteines'],
        'A->1' => ['generate_besttargetted_index'],
      },
      -input_ids => [{}],
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
        fan_branch_code => 1,
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
        '1' => ['indicate_proteome', 'generate_exonerate_jobs'],
      },
    },
    {
      -logic_name => 'generate_exonerate_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type => 'sequence_accession',
        sequence_table_name => $self->o('uniprot_table_name'),
        batch_size => 1,
        constraint => 'biotype IN ("self_pe12_sp", "self_pe12_tr", "self_isoforms_12_sp", "self_isoforms_12_tr", "ncbi_self_refseq")',
      },
      -rc_name      => 'default',
      -flow_into => {
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
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_sub')},
        GENOMICSEQS         => $self->o('genome_file'),
        PROGRAM             => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
        calculate_coverage_and_pid => $self->o('refseq_cdna_calculate_coverage_and_pid'),
        exonerate_cdna_pid => $self->o('exonerate_protein_pid'),
        exonerate_cdna_cov => $self->o('exonerate_protein_cov'),
      },
      -rc_name          => 'exonerate',
    },

    {
      -logic_name => 'indicate_proteome',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => '#indicate_path# -d #indicate_dir# -f #proteome# -o proteome_index -p singleWordParser',
        indicate_path => $self->o('indicate_path'),
        proteome => 'proteome.fa',
        indicate_dir => $self->o('targetted_path'),
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
      -rc_name      => 'default',
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
      -rc_name          => 'default',
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
      -rc_name          => 'default',
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
        logic_name => 'bestpmatch',
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
        target_db => $self->o('genewise_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('dna_db'),
        gtag => 0, # 0 is for gtag, 1 is for non canonical
        biotype => 'gw_gtag',
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_genewise')},
      },
      -rc_name          => 'exonerate',
    },
    {

      -logic_name => 'targetted_genewise_nogtag',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise',
      -parameters => {
        target_db => $self->o('genewise_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('dna_db'),
        gtag => 1, # 0 is for gtag, 1 is for non canonical
        biotype => 'gw_nogtag',
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_genewise')},
      },
      -rc_name          => 'exonerate',
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
        program_file => $self->o('exonerate_path'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_exonerate')},
      },
      -rc_name          => 'exonerate',
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
        table_name => $self->o('cdna_table_name'),
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
      -rc_name    => 'exonerate',
      -parameters => {
        iid_type => 'db_seq',
        query_table_name => $self->o('cdna_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('cdna_db'),
        logic_name => $self->o('exonerate_logic_name'),
        GENOMICSEQS         => $self->o('genome_file'),
        PROGRAM             => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_bestn_sub')},
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
      -rc_name    => 'exonerate_6G',
      -parameters => {
        iid_type => 'db_seq',
        query_table_name => $self->o('cdna_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('cdna_db'),
        logic_name => $self->o('exonerate_logic_name'),
        GENOMICSEQS         => $self->o('genome_file'),
        PROGRAM             => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_bestn_sub')},
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
        logic_name => $self->o('species_name').'_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna_selection')},
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
        logic_name => [$self->o('exonerate_logic_name')],
        coord_system_name => 'toplevel',
        include_non_reference => 0,
        top_level => 1,
        use_annotation => 1,
# These options will create only slices that have a gene on the slice in one of the feature dbs
        annotation_file => $self->o('annotation_file'),
        region_padding => $self->o('region_padding'),
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
      -rc_name    => 'exonerate',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        query_table_name => $self->o('cdna_table_name'),
        source_db => $self->o('cdna_db'),
        target_db => $self->o('cdna2genome_db'),
        logic_name => 'cdna2genome',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna2genome')},
        calculate_coverage_and_pid => 1,
        program => $self->o('exonerate_annotation'),
        annotation_file => $self->o('annotation_file'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -batch_size => 10,
      -flow_into => {
        '-1' => ['cdna2genome_himem'],
      },
    },

    {
      -logic_name => 'cdna2genome_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion',
      -rc_name    => 'exonerate_6G',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        query_table_name => $self->o('cdna_table_name'),
        source_db => $self->o('cdna_db'),
        target_db => $self->o('cdna2genome_db'),
        logic_name => 'cdna2genome',
        module     => 'HiveExonerate2Genes',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna2genome')},
        calculate_coverage_and_pid => 1,
        program => $self->o('exonerate_annotation'),
        annotation_file => $self->o('annotation_file'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
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
        multi_query_download => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', 'selenocysteine'),
        dest_dir => $self->o('targetted_path'),
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
        source_db => $self->o('genewise_db'),
        target_db => $self->o('genewise_db'),
        dna_db => $self->o('dna_db'),
        genome => $self->o('genome_file'),
        biotype => 'seleno_self',
      },
      -rc_name          => 'default',
    },
    {
      -logic_name => 'generate_besttargetted_index',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateBestTargettedIndex',
      -parameters => {
        source_db        => $self->o('cdna2genome_db'),
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        fasta_file => catfile($self->o('targetted_path'), 'best_targetted.fa'),
        email => $self->o('email_address'),
        genbank_file => $self->o('cdna_file'),
        files => [catfile($self->o('targetted_path'), 'proteome.fa')],
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
        2 => ['best_targetted'],
      },
    },
    {

      -logic_name => 'best_targetted',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestTargetted',
      -parameters => {
        target_db => $self->o('genewise_db'),
        dna_db => $self->o('dna_db'),
        source_db => { protein_db => $self->o('genewise_db'), cdna2genome_db => $self->o('cdna2genome_db')},
        SEQFETCHER_DIR => [catfile($self->o('targetted_path'), 'proteome_index'),
        catfile($self->o('targetted_path'), 'best_targetted_index')],
        INPUT_DATA_FROM_DBS  => {
          protein_db => ['seleno_self', 'gw_gtag', 'gw_nogtag', 'gw_exo'],
          cdna2genome_db => ['cdna2genome', 'edited'],
        },
        BIOTYPES => ['seleno_self', 'cdna2genome', 'edited', 'gw_gtag', 'gw_nogtag', 'gw_exo'], # sorted list, preferred is first
      },
      -rc_name          => 'exonerate',
    },
  ]
}


sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder('production-rh7', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'exonerate' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'exonerate_6G' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
  }
}

1;
