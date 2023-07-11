
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::IGTR_subpipeline;

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
    'user_r'                    => '', # read only db user
    'user'                      => '', # write db user
    'password'                  => '', # password for write db user
    'pipe_db_host'              => '', # host for pipe db
    'dna_db_host'               => '', # host for dna db
    'pipe_db_port'              => '', # port for pipeline host
    'dna_db_port'               => '', # port for dna db host

    'release_number'            => '' || $self->o('ensembl_release'),
    'species_name'              => '', # e.g. mus_musculus
    'production_name'           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    dbname_accession          => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1

    # need for IGTR part:
    'ig_tr_fasta_file'          => 'human_ig_tr.fa', # What IMGT fasta file to use. File should contain protein segments with appropriate headers
    'uniprot_set'               => '', # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'sanity_set'                => '', #sanity checks

########################
    # Pipe and ref db info
########################

# The following might not be known in advance, since the come from other pipelines
# These values can be replaced in the analysis_base table if they're not known yet
# If they are not needed (i.e. no projection or rnaseq) then leave them as is

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('dbname_accession').'_pipe_'.$self->o('release_number'),
    'dna_db_name'                   => $self->o('dbowner').'_'.$self->o('dbname_accession').'_core_'.$self->o('release_number'),

    'ig_tr_db_host'                => $self->o('databases_host'),
    'ig_tr_db_port'                => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database

########################
    # BLAST db paths
########################
    'base_blast_db_path'        => $ENV{BLASTDB_DIR},
    'ig_tr_blast_path'          => catfile($self->o('base_blast_db_path'), 'ig_tr_genes'),

######################################################
    #
    # Mostly constant settings
    #
######################################################

    genome_dumps                  => catdir($self->o('output_path'), 'genome_dumps'),

    # This is used for "messaging" other sub pipeline
    transcript_selection_url => undef,

# This one is used by most analyses that run against a genome flatfile like exonerate, genblast etc. Has slice name style headers. Is softmasked
    genome_file        => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),
    softmasked_genome_file        => catfile($self->o('genome_dumps'), $self->o('species_name').'_softmasked_toplevel.fa'),
    ensembl_analysis_script           => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    load_fasta_script_path            => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'load_fasta_to_db_table.pl'),

########################
    # Executable paths
########################
    'blast_type' => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'

    'genblast_path'     => catfile($self->o('binary_base'), 'genblast'),
    'genblast_eval'     => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
    'genblast_pid'      => '30',
    'genblast_max_rank' => '5',

    'ig_tr_table_name'    => 'ig_tr_sequences',
    'ig_tr_genblast_cov'  => '0.8',
    'ig_tr_genblast_pid'  => '70',
    'ig_tr_genblast_eval' => '1',
    'ig_tr_genblast_max_rank' => '5',
    'ig_tr_batch_size'    => 10,

########################
    # db info
########################

    'ig_tr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_igtr_'.$self->o('release_number'),
      -host   => $self->o('ig_tr_db_host'),
      -port   => $self->o('ig_tr_db_port'),
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
    genome_file => $self->o('genome_file'),
  }
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  my %genblast_params = (
    wu_genome    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
    ncbi_genome  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
    );
  my %commandline_params = (
    'ncbi' => '-num_threads 3 -window_size 40',
    'wu' => '-cpus 3 -hitdist 40',
    'legacy_ncbi' => '-a 3 -A 40',
    );

  return [

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
      -input_ids  => [{}],
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
      -rc_name    => '4GB',
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
      -rc_name          => '5GB',
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
          'gene_db_checks')->{$self->o('sanity_set')}->{'ig_tr'},
        },

      -rc_name    => '2GB',

      -flow_into => {
        '1' => ['cluster_ig_tr_genes'],
        },
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
      -rc_name    => '4GB',
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
          data => $self->o('ig_tr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'split_slices_on_intergenic',
          param => 'input_gene_dbs',
          data => $self->o('ig_tr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'layer_annotation',
          param => 'SOURCEDB_REFS',
          data => $self->o('ig_tr_db'),
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
