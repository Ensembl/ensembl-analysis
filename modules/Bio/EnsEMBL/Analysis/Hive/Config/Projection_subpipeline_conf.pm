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

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::Projection_subpipeline_conf;

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
    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'             => '' || $self->o('production_name').'_'.$self->o('ensembl_release'),
    'production_name'           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'release_number'            => '' || $self->o('ensembl_release'),
    'user_r'                    => '', # read only db user
    'user'                      => '', # write db user
    'password'                  => '', # password for write db user
    'pipe_db_server'            => '', # host for pipe db
    'dna_db_server'             => '', # host for dna db
    'databases_server'          => '', # host for general output dbs
    'pipe_db_port'              => '', # port for pipeline host
    'dna_db_port'               => '', # port for dna db host
    'databases_port'            => '', # port for general output db host
    'uniprot_set'               => '', # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'skip_projection'           => '0', # Will skip projection process if 1

########################
# Pipe and ref db info
########################

    'projection_source_db_name'         => '', # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_server'       => 'mysql-ens-mirror-1',
    'projection_source_db_port'         => '4240',

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is
    'pipe_db_name'  => $self->o('dbowner').'_'.$self->o('production_name').'_pipe_'.$self->o('release_number'),

    'projection_lastz_db_name'     => $self->o('dbowner').'_'.$self->o('production_name').'_lastz_pipe_'.$self->o('release_number'),
    'projection_lastz_db_server'   => $self->o('pipe_db_server'),
    'projection_lastz_db_port'     => $self->o('pipe_db_port'),

    'dna_db_name'   => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),
    faidx_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_toplevel.fa' ),
    faidx_softmasked_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_softmasked_toplevel.fa.reheader' ),

    'projection_db_server'  => $self->o('databases_server'),
    'projection_db_port'    => $self->o('databases_port'),

    'selected_projection_db_server'  => $self->o('databases_server'),
    'selected_projection_db_port'    => $self->o('databases_port'),

    databases_to_delete => ['projection_db', 'selected_projection_db',],

######################################################
#
# Mostly constant settings
#
######################################################

    # This is used for "messaging" other sub pipeline
    transcript_selection_url => undef,

    ensembl_analysis_script           => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    flag_potential_pseudogenes_script => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'flag_potential_pseudogenes.pl'),

########################
# Extra db settings
########################

    'num_tokens' => 10,

########################
# Executable paths
########################

    'cesar_path' => catdir($self->o('software_base_path'),'opt','cesar','bin'),

# Max internal stops for projected transcripts
    'projection_pid'                        => '50',
    'projection_cov'                        => '50',
    'projection_max_internal_stops'         => '1',
    'projection_calculate_coverage_and_pid' => '1',

########################
# db info
########################

    'projection_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_proj_'.$self->o('release_number'),
      -host   => $self->o('projection_db_server'),
      -port   => $self->o('projection_db_port'),
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

    'selected_projection_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_sel_proj_'.$self->o('release_number'),
      -host   => $self->o('selected_projection_db_server'),
      -port   => $self->o('selected_projection_db_port'),
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
    genome_file => $self->o('faidx_genome_file'),
  }
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;
  return [
    {
      -logic_name => 'create_selected_projection_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('selected_projection_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -input_ids  => [{}],
      -flow_into  => {
        '1->A' => ['create_projection_db'],
        'A->1' => ['classify_projected_genes'],
      },
    },

    {
      -logic_name => 'create_projection_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('projection_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => ['cesar_create_projection_input_ids','wga_create_projection_input_ids'],
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
        logic_name                  => 'project_transcripts',
        module                      => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
        source_dna_db               => $self->default_options()->{'projection_source_db'},
        target_dna_db               => $self->o('dna_db'),
        source_transcript_db        => $self->default_options()->{'projection_source_db'},
        target_transcript_db        => $self->o('projection_db'),
        compara_db                  => $self->o('projection_lastz_db'),
        method_link_type            => 'LASTZ_NET',
        max_exon_readthrough_dist   => 15,
        TRANSCRIPT_FILTER => {
          OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
          PARAMETERS => {
            -coverage   => $self->o('projection_cov'),
            -percent_id => $self->o('projection_pid'),
          },
        },
        iid_type                    => 'feature_id',
        feature_type                => 'transcript',
        calculate_coverage_and_pid  => $self->o('projection_calculate_coverage_and_pid'),
        max_internal_stops          => $self->o('projection_max_internal_stops'),
        timer                       => '30m',
      },
      -flow_into  => {
        -1 => ['wga_project_transcripts_himem'],
        -3 => ['wga_failed_projection'],
      },
      -rc_name        => 'project_transcripts',
      -batch_size     => 100,
      -hive_capacity  => $self->hive_capacity_classes->{'hc_high'},
    },

    {
      -logic_name => 'wga_project_transcripts_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
      -parameters => {
        logic_name => 'project_transcripts',
        module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
        source_dna_db             => $self->default_options()->{'projection_source_db'},
        target_dna_db             => $self->o('dna_db'),
        source_transcript_db      => $self->default_options()->{'projection_source_db'},
        target_transcript_db      => $self->o('projection_db'),
        compara_db                => $self->o('projection_lastz_db'),
        method_link_type          => 'LASTZ_NET',
        max_exon_readthrough_dist => 15,
        TRANSCRIPT_FILTER => {
          OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
          PARAMETERS => {
            -coverage   => $self->o('projection_cov'),
            -percent_id => $self->o('projection_pid'),
          },
        },
        iid_type                    => 'feature_id',
        feature_type                => 'transcript',
        calculate_coverage_and_pid  => $self->o('projection_calculate_coverage_and_pid'),
        max_internal_stops          => $self->o('projection_max_internal_stops'),
        timer                       => '30m',
      },
      -flow_into => {
        -3 => ['wga_failed_projection'],
      },
      -rc_name        => 'project_transcripts_himem',
      -batch_size     => 100,
      -hive_capacity  => $self->hive_capacity_classes->{'hc_high'},
      -can_be_empty   => 1,
    },

    {
      -logic_name           => 'wga_failed_projection',
      -module               => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters           => {},
      -rc_name              => 'default',
      -can_be_empty         => 1,
      -failed_job_tolerance => 100,
    },

    {
      -logic_name => 'cesar_create_projection_input_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db                        => $self->o('projection_source_db'),
        iid_type                         => 'feature_id',
        batch_size                       => 20,
        feature_type                     => 'gene',
        feature_id_include_non_reference => 0,
        feature_restriction              => 'biotype',
        biotypes                         => {'protein_coding' => 1},
      },
      -flow_into  => {
        2 => ['cesar'],
      },
      -rc_name    => '5GB',
    },

    {
      -logic_name => 'cesar',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar',
      -parameters => {
        'output_path'           => $self->o('output_path')."/cesar_projection/",
        'source_dna_db'         => $self->default_options()->{'projection_source_db'},
        'target_dna_db'         => $self->o('dna_db'),
        'source_db'             => $self->o('projection_source_db'),
        'target_db'             => $self->o('projection_db'),
        'compara_db'            => $self->o('projection_lastz_db'),
        'method_link_type'      => 'LASTZ_NET',
        'cesar_path'            => $self->o('cesar_path'),
        'cesar_mem'             => '3',
        'canonical_or_longest'  => 1,
        'stops2introns'         => 1,
      },
      -rc_name              => '3GB',
      -analysis_capacity    => 50,
      -max_retry_count      => 1,
      -failed_job_tolerance => 5,
      -flow_into => {
        15 => ['cesar_15'],
        30 => ['cesar_30'],
        -1 => ['cesar_30'],
      },
    },

    {
      -logic_name => 'cesar_15',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar',
      -parameters => {
        'output_path'           => $self->o('output_path')."/cesar_projection/",
        'source_dna_db'         => $self->default_options()->{'projection_source_db'},
        'target_dna_db'         => $self->o('dna_db'),
        'source_db'             => $self->o('projection_source_db'),
        'target_db'             => $self->o('projection_db'),
        'compara_db'            => $self->o('projection_lastz_db'),
        'method_link_type'      => 'LASTZ_NET',
        'cesar_path'            => $self->o('cesar_path'),
        'cesar_mem'             => '15',
        'canonical_or_longest'  => 1,
        'stops2introns'         => 1,
      },
      -rc_name              => '15GB',
      -analysis_capacity    => 50,
      -max_retry_count      => 1,
      -can_be_empty         => 1,
      -failed_job_tolerance => 5,
      -flow_into => {
        -1 => ['cesar_30'],
      },
    },

    {
      -logic_name => 'cesar_30',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar',
      -parameters => {
        'output_path'           => $self->o('output_path')."/cesar_projection/",
        'source_dna_db'         => $self->default_options()->{'projection_source_db'},
        'target_dna_db'         => $self->o('dna_db'),
        'source_db'             => $self->o('projection_source_db'),
        'target_db'             => $self->o('projection_db'),
        'compara_db'            => $self->o('projection_lastz_db'),
        'method_link_type'      => 'LASTZ_NET',
        'cesar_path'            => $self->o('cesar_path'),
        'cesar_mem'             => '30',
        'canonical_or_longest'  => 1,
        'stops2introns'         => 1,
      },
      -rc_name              => '30GB',
      -analysis_capacity    => 50,
      -max_retry_count      => 1,
      -can_be_empty         => 1,
      -failed_job_tolerance => 10,
      -flow_into => {
        -1 => ['cesar_failed_projection'],
      },
    },

    {
      -logic_name           => 'cesar_failed_projection',
      -module               => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters           => {},
      -rc_name              => 'default',
      -failed_job_tolerance => 100,
      -can_be_empty         => 1,
    },

    {
      -logic_name => 'classify_projected_genes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
      -parameters => {
        skip_analysis       => $self->o('skip_projection'),
        classification_type => 'standard',
        update_gene_biotype => 1,
        target_db           => $self->o('projection_db'),
      },
      -rc_name    => '2GB',
      -flow_into  => {
        1 => ['flag_problematic_projections'],
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
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn   => $self->o('projection_db'),
        sql => [
          'UPDATE gene JOIN transcript USING(gene_id) SET transcript.biotype=gene.biotype',
          'UPDATE protein_align_feature JOIN transcript_supporting_feature ON feature_id = protein_align_feature_id'.
            ' JOIN transcript USING(transcript_id) SET hit_name = stable_id',
          'UPDATE protein_align_feature JOIN supporting_feature ON feature_id = protein_align_feature_id'.
            ' JOIN exon_transcript USING(exon_id) JOIN transcript USING(transcript_id) SET hit_name = stable_id',
        ],
      },
      -rc_name    => 'default',
      -flow_into  => {
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
          data => $self->o('selected_projection_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'split_slices_on_intergenic',
          param => 'input_gene_dbs',
          data => $self->o('selected_projection_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'layer_annotation',
          param => 'SOURCEDB_REFS',
          data => $self->o('selected_projection_db'),
          update => 1,
        }],
        tweak_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'tweak_pipeline.pl'),
      },
      -rc_name    => 'default',
    },
  ];
}

sub hive_capacity_classes {
  my $self = shift;
  return {
    'hc_high' => 1000,
  };
}

sub resource_classes {
  my $self = shift;
  return {
    'default' => { LSF => $self->lsf_resource_builder('production', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '2GB' => { LSF => $self->lsf_resource_builder('production', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '3GB' => { LSF => $self->lsf_resource_builder('production', 3000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB' => { LSF => $self->lsf_resource_builder('production', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB' => { LSF => $self->lsf_resource_builder('production', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '15GB' => { LSF => $self->lsf_resource_builder('production', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '30GB' => { LSF => $self->lsf_resource_builder('production', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'project_transcripts' => { LSF => $self->lsf_resource_builder('production', 3000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'projection_db_server'}, $self->default_options->{'projection_lastz_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'project_transcripts_himem' => { LSF => $self->lsf_resource_builder('production', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'projection_db_server'}, $self->default_options->{'projection_lastz_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
  }
}
1;
