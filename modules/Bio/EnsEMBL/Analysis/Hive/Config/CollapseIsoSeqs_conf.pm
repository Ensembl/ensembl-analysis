=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

CollapseIsoSeqs_conf

=head1 SYNOPSIS


=head1 DESCRIPTION

Pipeline to be run once all IsoSeq tissue samples have been aligned.
It will collapse the transcripts and run blast to check the coding
potential of the models then classify them

=cut

package CollapseIsoSeqs_conf;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw (get_analysis_settings);
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;


sub default_options {
  my ($self) = @_;

  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

##########################################################################
#                                                                        #
# CHANGE STUFF HERE                                                      #
#                                                                        #
##########################################################################

    'pipeline_name' => '',

    'user_r'   => '',
    'user'     => '',
    'password' => '',
    'port'     => '',
    'host'     => '',

    'dna_db_name'    => '',
    'dna_db_server' => '',
    'dna_db_port'   => '',

    'data_dbs_server' => '',
    'data_dbs_port'   => '',

    'intron_db_name'   => '',
    'intron_db_server' => '',
    'intron_db_port'   => 4529,

    'isoseq1_db_name'   => '',
    'isoseq1_db_server' => '',
    'isoseq1_db_port'   => 4528,

    'isoseq2_db_name'   => '',
    'isoseq2_db_server' => '',
    'isoseq2_db_port'   => 4529,

    'isoseq_dbs' => [$self->o('isoseq1_db'), $self->o('isoseq2_db')],

    'uniprotdb'    => '',# Uniprot blast index
    'uniprotindex' => '', #Uniprot OBDA (indicate) index
    'blastp'       => catfile($self->o('binary_base'), 'blastp'), #You may need to specify the full path to the blastp binary
    'blast_type'   => 'ncbi',
    'use_threads'  => 3,

    'collapse_db_name'   => $self->o('dbowner').'_'.$self->o('pipeline_name').'_collapse',
    'collapse_db_server' => $self->o('data_dbs_server'),
    'collapse_db_port'   => $self->o('data_dbs_port'),

    'blast_db_name'   => $self->o('dbowner').'_'.$self->o('pipeline_name').'_blast',
    'blast_db_server' => $self->o('data_dbs_server'),
    'blast_db_port'   => $self->o('data_dbs_port'),

    'check_db_name'   => $self->o('dbowner').'_'.$self->o('pipeline_name').'_check',
    'check_db_server' => $self->o('data_dbs_server'),
    'check_db_port'   => $self->o('data_dbs_port'),

##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

    'create_type' => 'clone',

    'collapse_db' => {
      -dbname => $self->o('collapse_db_name'),
      -host => $self->o('collapse_db_server'),
      -port => $self->o('collapse_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'blast_db' => {
      -dbname => $self->o('blast_db_name'),
      -host => $self->o('blast_db_server'),
      -port => $self->o('blast_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'check_db' => {
      -dbname => $self->o('check_db_name'),
      -host => $self->o('check_db_server'),
      -port => $self->o('check_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'intron_db' => {
      -dbname => $self->o('intron_db_name'),
      -host => $self->o('intron_db_server'),
      -port => $self->o('intron_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'isoseq1_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('isoseq1_db_name').'_best1_exonerate',
      -host   => $self->o('isoseq1_db_server'),
      -port   => $self->o('isoseq1_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'isoseq2_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('isoseq2_db_name').'_best1_exonerate',
      -host   => $self->o('isoseq2_db_server'),
      -port   => $self->o('isoseq2_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },


    databases_to_delete => ['collapse_db', 'blast_db', 'check_db'],
  };
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      # need to make sure the database actually copies as if it doesn't the job does appear to complete according to eHive
      -logic_name => 'create_collapse_db',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('collapse_db'),
        create_type => $self->o('create_type'),
      },
      -rc_name => 'default',
      -input_ids => [{}],
      -max_retry_count => 0,
      -flow_into => {
        1 => ['create_blast_db'],
      }
    },
    {
      # need to make sure the database actually copies as if it doesn't the job does appear to complete according to eHive
      -logic_name => 'create_blast_db',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('blast_db'),
        create_type => $self->o('create_type'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['create_check_db'],
      }
    },
    {
      -logic_name => 'create_check_db',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('check_db'),
        create_type => $self->o('create_type'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['generate_collapse_jobs'],
      }
    },
    {
      -logic_name => 'generate_collapse_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
                       target_db        => $self->o('collapse_db'),
                       feature_dbs => $self->o('isoseq_dbs'),
                       coord_system_name => 'toplevel',
                       iid_type => 'slice',
                       feature_constraint => 1,
                       feature_type => 'gene',
                       top_level => 1,
                     },
      -rc_name      => 'default',
      -max_retry_count => 1,
      -flow_into => {
                      '2->A' => ['split_slices_on_intergenic'],
                      'A->1' => ['classify_isoseq_models'],
                    },
    },
    {
      -logic_name => 'split_slices_on_intergenic',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions',
      -parameters => {
        dna_db => $self->o('dna_db'),
        input_gene_dbs => $self->o('isoseq_dbs'),
        iid_type => 'slice',
      },
      -batch_size => 100,
      -rc_name    => '5GB',
      -flow_into => {
        2 => ['collapse_transcripts'],
      },
    },
    {
      -logic_name => 'collapse_transcripts',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer',
      -parameters => {
                       target_db        => $self->o('collapse_db'),
                       dna_db        => $self->o('dna_db'),
                       source_dbs        => $self->o('isoseq_dbs'),
                       biotypes => ["isoseq"],
                     },
      -rc_name      => '5GB',
      -flow_into => {
                      1 => ['blast_isoseq'],
                      -1 => ['collapse_transcripts_20GB'],
                    },
    },
    {
      -logic_name => 'collapse_transcripts_20GB',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer',
      -parameters => {
                       target_db        => $self->o('collapse_db'),
                       dna_db        => $self->o('dna_db'),
                       source_dbs        => $self->o('isoseq_dbs'),
                       biotypes => ["isoseq"],
                     },
      -rc_name      => '20GB',
      -flow_into => {
                      1 => ['blast_isoseq'],
                    },
    },
    {
      -logic_name => 'blast_isoseq',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
      -parameters => {

          input_db => $self->o('collapse_db'),
          output_db => $self->o('blast_db'),
          dna_db => $self->o('dna_db'),

          # path to index to fetch the sequence of the blast hit to calculate % coverage
          indicate_index => $self->o('uniprotindex'),
          uniprot_index => [$self->o('uniprotdb')],
          blast_program => $self->o('blastp'),
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
          commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40 -seg no',
                    },
      -rc_name => '2GB_blast',
      -flow_into => {
                      -1 => ['blast_isoseq_10G'],
                      1 => ['intron_check'],
                    },
    },
    {
      -logic_name => 'blast_isoseq_10G',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
      -parameters => {

          input_db => $self->o('collapse_db'),
          output_db => $self->o('blast_db'),
          dna_db => $self->o('dna_db'),

          # path to index to fetch the sequence of the blast hit to calculate % coverage
          indicate_index => $self->o('uniprotindex'),
          uniprot_index => [$self->o('uniprotdb')],
          blast_program => $self->o('blastp'),
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
          commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40 -seg no',
                    },
      -rc_name => '10GB_blast',
      -flow_into => {
                      1 => ['intron_check'],
                    },
    },
    {
      -logic_name => 'intron_check',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck',
      -parameters => {
        source_db => $self->o('blast_db'),
        target_db => $self->o('check_db'),
        intron_db => $self->o('intron_db'),
        dna_db => $self->o('dna_db'),
      },
      -rc_name    => 'default',
    },
    {
      -logic_name => 'classify_isoseq_models',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
      -parameters => {
        classification_type => 'standard',
        update_gene_biotype => 1,
        target_db => $self->o('check_db'),
      },
      -rc_name    => 'default',
    },
  ];
}


sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder('production-rh7', 1000) },
    '5GB' => { LSF => $self->lsf_resource_builder('production-rh7', 5000) },
    '20GB' => { LSF => $self->lsf_resource_builder('production-rh7', 20000) },
    '2GB_blast' => { LSF => $self->lsf_resource_builder('production-rh7', 2000, undef, undef, ($self->default_options->{'use_threads'}+1))},
    '10GB_blast' => { LSF => $self->lsf_resource_builder('production-rh7', 10000, undef, undef, ($self->default_options->{'use_threads'}+1))},
  }
}


1;
