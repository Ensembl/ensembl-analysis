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

package Bio::EnsEMBL::Analysis::Hive::Config::ShortncRNA;

use strict;
use warnings;
use File::Spec::Functions;

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
    dbowner                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    pipeline_name             => '' || $self->o('production_name').'_'.$self->o('ensembl_release'),
    user_r                    => '', # read only db user
    user                      => '', # write db user
    password                  => '', # password for write db user
    server_set                => '', # What server set to user, e.g. set1
    pipe_db_server            => '', # host for pipe db
    databases_server          => '', # host for general output dbs
    dna_db_server             => '', # host for dna db
    pipe_db_port              => '', # port for pipeline host
    databases_port            => '', # port for general output db host
    dna_db_port               => '', # port for dna db host
    repbase_logic_name        => '', # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    release_number            => '' || $self->o('ensembl_release'),
    species_name              => '', # e.g. mus_musculus
    production_name           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    taxon_id                  => '', # should be in the assembly report file
    uniprot_set               => '', # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    output_path               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses


########################
## Small ncRNAs params
#########################
    mirBase_fasta             => 'all_mirnas.fa', # What mirBase file to use. It is currently best to use on with the most appropriate set for your species
    rfc_scaler                => 'filter_dafs_rfc_scaler_human.pkl',
    rfc_model                 => 'filter_dafs_rfc_model_human.pkl',

    # Clade-based filtering on rfam accessions
    # Rfam db details should stay constant but check periodically
    rfam_user => 'rfamro',
    rfam_dbname => 'Rfam',
    rfam_host => 'mysql-rfam-public.ebi.ac.uk',
    rfam_port => 4497,

    rfam_path => catfile($self->o('base_blast_db_path'), 'ncrna', 'Rfam_14.0'),
    rfam_seeds => catfile($self->o('rfam_path'), 'Rfam.seed'),
    rfam_cm => catfile($self->o('rfam_path'), 'Rfam.cm'),
    filtered_rfam_cm => catfile($self->o('ncrna_dir'), 'Rfam.cm'),
    clade => $self->o('repbase_logic_name'),


########################
# Pipe and ref db info
########################
    pipe_db_name                  => $self->o('dbowner').'_'.$self->o('production_name').'_pipe_'.$self->o('release_number'),
    dna_db_name                   => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),

    ncrna_db_server              => $self->o('databases_server'),
    ncrna_db_port                => $self->o('databases_port'),
    ncrna_db_name                  => $self->o('dbowner').'_'.$self->o('production_name').'_ncrna_'.$self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    ensembl_release              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database


    databases_to_delete => ['ncrna_db'],

########################
# BLAST db paths
########################
    base_blast_db_path        => $ENV{BLASTDB_DIR},
    ncrna_blast_path          => catdir($self->o('base_blast_db_path'), 'ncrna', 'ncrna_2016_05'),
    mirna_blast_path          => catdir($self->o('base_blast_db_path'), 'ncrna', 'mirbase_22'),

######################################################
#
# Mostly constant settings
#
######################################################

    genome_dumps                  => catdir($self->o('output_path'), 'genome_dumps'),
    # This one is used in replacement of the dna table in the core db, so where analyses override slice->seq. Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file             => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),
    use_genome_flatfile           => '1',# This will read sequence where possible from a dumped flatfile instead of the core db

    ncrna_dir => catdir($self->o('output_path'), 'ncrna'),

    ensembl_analysis_script           => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    sncrna_analysis_script             => catdir($self->o('ensembl_analysis_script'), 'genebuild', 'sncrna'),


########################
# Extra db settings
########################
    num_tokens => 10,

########################
# Executable paths
########################
    blast_type => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    blastn_exe_path => catfile($self->o('binary_base'), 'blastn'),
    cmsearch_exe_path    => catfile($self->o('binary_base'), 'cmsearch'), # #'opt', 'infernal10', 'bin', 'cmsearch'),

########################
# Misc setup info
########################
    repeatmasker_engine       => 'crossmatch',
    masking_timer_long        => '5h',
    masking_timer_short       => '2h',

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

########################
# db info
########################
    ncrna_db => {
      -dbname => $self->o('ncrna_db_name'),
      -host   => $self->o('ncrna_db_server'),
      -port   => $self->o('ncrna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },
  };
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    use_genome_flatfile => $self->o('use_genome_flatfile'),
    genome_file          => $self->o('faidx_genome_file'),
  }
}


sub pipeline_create_commands {
  my ($self) = @_;

  return [

    # inheriting database and hive tables' creation
    @{ $self->SUPER::pipeline_create_commands },

    'mkdir -p ' . $self->o('ncrna_dir'),
  ];
}


sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'create_ncrna_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('ncrna_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -input_ids  => [{}],
      -flow_into => {
        '1->A' => ['dump_repeats', 'fetch_rfam_accessions'],
        'A->1' => ['create_small_rna_slice_ids'],
      },
    },

    {
      -logic_name => 'dump_repeats',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.catfile($self->o('sncrna_analysis_script'), 'repeats_dump.pl').' '.
          $self->o('dna_db_name'). ' ' .
          $self->o('dna_db_server') . ' ' .
          $self->o('dna_db_port') . ' ' .
          $self->o('user_r') . ' ' .
          $self->o('ncrna_dir').' blastmirna',
      },
      -rc_name => '6GB',
    },

    {
      -logic_name => 'fetch_rfam_accessions',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . catfile($self->o('sncrna_analysis_script'), 'fetch_rfam_accessions.pl') .
          ' -h ' . $self->o('rfam_host') .
          ' -u ' . $self->o('rfam_user') .
          ' -p ' . $self->o('rfam_port') .
          ' -d ' . $self->o('rfam_dbname') .
          ' -c ' . $self->o('clade') .
          ' -o ' . $self->o('ncrna_dir'),
      },
      -rc_name => 'default',
      -flow_into => { '1' => 'extract_rfam_cm'},
    },

    {
      -logic_name => 'extract_rfam_cm',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.catfile($self->o('sncrna_analysis_script'), 'filter_cm.pl').
          ' '.$self->o('rfam_cm').
          ' '.catfile($self->o('ncrna_dir'), 'accessions.txt').
          ' '.$self->o('filtered_rfam_cm'),
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
        output_dir => $self->o('ncrna_dir'),
      },
      -hive_capacity => 900,
      -rc_name    => '5GB',
      -flow_into => {
        '-1' => 'run_cmsearch_highmem',
        '-2' => 'run_cmsearch_highmem'},
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
        output_dir => $self->o('ncrna_dir'),
      },
      -hive_capacity => 900,
      -rc_name    => '10GB',
    },

    {
      -logic_name => 'mirna_blast',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastmiRNA',
      -parameters => {
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
        output_db => $self->o('ncrna_db'),
        dna_db => $self->o('dna_db'),
        logic_name => 'blastmirna',
        module     => 'HiveBlastmiRNA',
        blast_db_path => catfile($self->o('mirna_blast_path'), $self->o('mirBase_fasta')),
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
      -parameters => {},
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
        'A->1' => ['concat_rnafold_result'],
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
    },

    {
      -logic_name => 'concat_rnafold_result',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cat #data_file_pattern# > #output_file#',
        data_file_pattern => catfile($self->o('ncrna_dir'), 'rna_fold_*.part'),
        output_file => catfile($self->o('ncrna_dir'), 'rna_fold_results.txt'),
      },
      -rc_name   => 'filter',
      -flow_into => {
        1 => 'fan_dump_features',
      },
    },

    {
      -logic_name => 'fan_dump_features',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => WHEN(
            '-e "'.catfile($self->o('ncrna_dir'),'rna_fold_results.txt').'"' => ['dump_features']
            ),
      },
    },

    {
      -logic_name => 'dump_features',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.catfile($self->o('sncrna_analysis_script'), 'dump_prefilter_features.pl')
          .' '.$self->o('ncrna_db_name')
          .' '.$self->o('ncrna_db_server')
          .' '.$self->o('ncrna_db_port')
          .' '.$self->o('user_r')
          .' '.$self->o('ncrna_dir')
          .' blastmirna',
      },
      -rc_name   => 'filter',
      -flow_into => {
        1 => 'dump_annotated_dafs',
      },
    },

    {
      -logic_name => 'dump_annotated_dafs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => "annotateBed -i #bedfile_dafs# -files #repeats_file# | bedtools nuc -fi #genome_file# -bed stdin -s | cut -f1-11,13 | grep -v '#' > #output_file# ",
        output_file => catfile($self->o('ncrna_dir'), 'annotated_dafs.tsv'),
        bedfile_dafs => catfile($self->o('ncrna_dir'), 'blastmirna_dafs.bed'),
        repeats_file => catfile($self->o('ncrna_dir'), 'repeats.bed'),
        genome_file => $self->o('faidx_genome_file'),
      },
      -rc_name   => 'filter',
      -flow_into => {
        1 => 'filter_mirnas'
      },
    },

    {
      -logic_name => 'filter_mirnas',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'PYENV_VERSION="#pyenv_virtualenv#" python '.catfile($self->o('sncrna_analysis_script'), 'FilterDafs.py')
          .' '.catfile($self->o('mirna_blast_path'), 'rfc_filters', $self->o('rfc_model'))
          .' '.catfile($self->o('mirna_blast_path'), 'rfc_filters', $self->o('rfc_scaler'))
          .' '.$self->o('ncrna_dir')
          .' '.catfile($self->o('ncrna_dir'), 'annotated_dafs.tsv')
          .' '.catfile($self->o('ncrna_dir'), 'rna_fold_results.txt')
          .' '.catfile($self->o('ncrna_dir'), 'identified_mirnas.bed')
          .' '.catfile($self->o('ncrna_dir'), 'mirnas_to_delete.txt'),
        pyenv_virtualenv => 'genebuild-mirna',
      },
      -rc_name   => 'filter',
      -flow_into => {
        1 => 'delete_flagged_mirnas',
      },
    },

    {
      -logic_name => 'delete_flagged_mirnas',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.catfile($self->o('ensembl_analysis_script'), 'genebuild', 'delete_genes.pl')
          .' -dbhost #expr(#target_db#->{-host})expr#'
          .' -dbport #expr(#target_db#->{-port})expr#'
          .' -dbuser #expr(#target_db#->{-user})expr#'
          .' -dbpass #expr(#target_db#->{-pass})expr#'
          .' -dbname #expr(#target_db#->{-dbname})expr#'
          .' -idfile '.catfile($self->o('ncrna_dir'), 'mirnas_to_delete.txt'),
        target_db => $self->o('ncrna_db'),
      },
      -rc_name   => 'filter',
      -flow_into => {
        1 => 'ncrna_sanity_checks',
      },
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
    },
  ];
}


sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder('production', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB' => { LSF => $self->lsf_resource_builder('production', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB' => { LSF => $self->lsf_resource_builder('production', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '6GB' => { LSF => $self->lsf_resource_builder('production', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '10GB' => { LSF => $self->lsf_resource_builder('production', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'blast' => { LSF => $self->lsf_resource_builder('production', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast_retry' => { LSF => $self->lsf_resource_builder('production', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'filter' => { LSF => $self->lsf_resource_builder('production', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
  }
}

1;
