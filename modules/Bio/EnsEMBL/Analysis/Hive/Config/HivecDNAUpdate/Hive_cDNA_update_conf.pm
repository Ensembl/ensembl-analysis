=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

# this draft will just do the first part of the pipeline - downlaod the fasta files, remove kill list obj and then clean and clip


package Hive_cDNA_update_conf;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::Utilities qw (get_analysis_settings);
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

my %taxon_id;
$taxon_id {"human"} = 9606;
$taxon_id {"mouse"} = 10090;

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
    
    'recipient_email'            => 'dmurphy@ebi.ac.uk', # email address where reports will be sent

    'species'                    => 'human',

    'ensembl_release'            => '87',

    'gb_user'                    => 'dm15',

    'pipeline_name'              => $self->o('species').'_cdna_update_'.$self->o('ensembl_release'),
    'pipe_db_name'               => $self->o('gb_user').'_'.$self->o('species').'_cdna_hive_'.$self->o('ensembl_release'),
    'pipe_db_server'             => 'genebuild11',

    # the dna db should be the one on livemirror so need to change this
    'dna_db_name'                => 'dm15_human_84_copy',
    'dna_db_server'              => 'genebuild12',

    'killlist_db_server'         => 'genebuild6',

    'exonerate_output_db_server' => 'genebuild13',

    'last_release'               => '86',

    'old_cdna_db_name'           => 'homo_sapiens_cdna_84_38',

    'old_cdna_file_name'         => '/lustre/scratch109/ensembl/dm15/humancdna_86/cdna_update.clipped',

    'output_path'                => '/lustre/scratch109/ensembl/dm15/humancdna_87/',

    'refseq_path'                => '/data/blastdb/Ensembl/RefSeq_2016_01/',
    'refseq_file'                => 'hs.fna',

    'ensembl_repo_root'          => $ENV{ENSCODE},
    'clone_db_script_path'       => $self->o('ensembl_repo_root').'/ensembl-analysis/scripts/clone_database.ksh',

    'genome_file'                => '/data/blastdb/Ensembl/human/GRCh38/genome/softmasked_dusted/toplevel.with_nonref_and_GRCh38_p7.no_duplicate.softmasked_dusted.fa',

    'repeat_masking_logic_names' => ['repeatmask_repbase_human'],

    'gss_file'                   => $self->o('ensembl_repo_root').'/ensembl-personal/genebuilders/cDNA_update/gss_acc.txt',

    'optimize_script'            => $self->o('ensembl_repo_root').'/ensembl-personal/genebuilders/scripts/load_external_db_ids_and_optimize_af.pl',

    'refseq_version'             => '77',

    'exonerate_version'          => 'exonerate-0.9.0',

##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

    'pipeline_name'              => $self->o('species').'_cdna_update_'.$self->o('ensembl_release'),
    'pipe_db_name'               => $self->o('gb_user').'_'.$self->o('species').'_cdna_'.$self->o('ensembl_release'),

    'exonerate_output_db_name'   => $self->o('gb_user').'_'.$self->o('species').'_cdna_exonerate_'.$self->o('ensembl_release'),

    'old_cdna_db_server'         => 'ens-livemirror',
     
    'production_db_name'         => 'ensembl_production',
    'production_db_server'       => 'ens-staging1',

    'exonerate_batch_size'       => '50',

    'fastasplit_random_path'     => '/software/ensembl/bin/fastasplit_random',

    'killlist_db_name'           => 'gb_kill_list',

    'cdna_file_name'             => 'cdna_update',

    'user_r'                     => 'ensro',
    'user_w'                     => 'ensadmin',
    'password'                   => '',
    'port'                       => '3306',

    'cdna_query_dir_name'        => 'cdna_temp',

    'retired_cdnas_file'         => $self->o('output_path').'/cdna_update.retired',
    'genes_to_delete'            => $self->o('output_path').'/genes_to_delete.ids',

    'many_hits_dir'              => 'many_hits',

    'cdna_table_name'            => 'cdna_sequences',
    'cdna_batch_size'            => '1',

    'default_mem'                => '900',
    'exonerate_mem'              => '3900',
    'exonerate_retry_mem'        => '5900',
    'optimise_mem'               => '11900',
    'download_mem'               => '8000',

    'exonerate_path'             => '/software/ensembl/genebuild/usrlocalensemblbin/exonerate-0.9.0',
    'exonerate_pid'              => '97',
    'exonerate_cov'              => '90',

    'polyA_script'               => $self->o('ensembl_repo_root').'/ensembl-pipeline/scripts/EST/new_polyA_clipping.pl',

    'delete_genes_script'        => $self->o('ensembl_repo_root').'/ensembl-analysis/scripts/genebuild/delete_genes.pl',

    'populate_production_script' => $self->o('ensembl_repo_root').'/ensembl-production/scripts/production_database/populate_production_db_tables.pl',

    'findN_script'               => $self->o('ensembl_repo_root').'/ensembl-pipeline/scripts/cDNA_update/find_N.pl',

    'analysis_desc_script'       => $self->o('ensembl_repo_root').'/ensembl-production/scripts/production_database/populate_analysis_description.pl'.

    'driver'                     => 'mysql',
    'num_tokens'                 => 10,

    'create_type'                => 'copy',

    'pipeline_db' => {
      -dbname => $self->o('pipe_db_name'),
      -host => $self->o('pipe_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user_w'),
      -pass => $self->o('password'),
      -driver => $self->o('driver'),
    },

    'production_db' => {
      -dbname => $self->o('production_db_name'),
      -host => $self->o('production_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user_r'),
    },
 
    'exonerate_output_db' => {
      -dbname => $self->o('exonerate_output_db_name'),
      -host => $self->o('exonerate_output_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user_w'),
      -pass => $self->o('password'),
    },

    'dna_db' => {
      -dbname => $self->o('dna_db_name'),
      -host => $self->o('dna_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user_r'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host => $self->o('killlist_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user_r'),
    },

    'old_cdna_db' => {
      -dbname => $self->o('old_cdna_db_name'),
      -host => $self->o('old_cdna_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user_r'),
    },
  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [
  # inheriting database and hive tables' creation
    @{$self->SUPER::pipeline_create_commands},

    $self->db_cmd('CREATE TABLE '.$self->o('cdna_table_name').' ('.
      'accession varchar(50) NOT NULL,'.
      'seq text NOT NULL,'.
      'biotype varchar(50) NOT NULL,'.
      'PRIMARY KEY (accession))'),
  ];
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      # need to make sure the database actually copies as if it doesn't the job does appear to complete according to eHive
      -logic_name => 'create_output_db',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('old_cdna_db'),
        target_db => $self->o('exonerate_output_db'),
        create_type => 'copy',
        script_path => $self->o('clone_db_script_path'),
        pass_w => $self->o('password'),
        user_w => $self->o('user_w'),
        user_r => $self->o('user_r'),
      },
      -rc_name => 'default',
      -input_ids => [{
        cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name'),
      }],
      -max_retry_count => 0,
      -flow_into => {
        1 => ['populate_production'],
      }
    },
    {
      -logic_name => 'populate_production',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('populate_production_script').
               ' -dp '.$self->o('output_path').
               ' -d '.$self->o('exonerate_output_db','-dbname').
               ' -h '.$self->o('exonerate_output_db','-host').
               ' -u '.$self->o('exonerate_output_db','-user').
               ' -p '.$self->o('exonerate_output_db','-pass').
               ' -md '.$self->o('production_db','-dbname').
               ' -mh '.$self->o('production_db','-host').
               ' -mu '.$self->o('production_db','-user').
               ' -t external_db -t attrib_type -t misc_set -t unmapped_reason'
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['download_cdnas'],
      }
    },
    {
      -logic_name => 'download_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadcDNAFiles',
      -parameters => {
        embl_sequences => {
          output_path => $self->o('output_path'),
          output_file => $self->o('cdna_file_name'),
          species => $self->o('species'),
        },
        refseq_sequences => {
          refseq_path => $self->o('refseq_path'),
          refseq_file => $self->o('refseq_file'),
        },
      },  
      -max_retry_count => 0,
      -rc_name => 'download',
      #-input_ids => [{}],
      -flow_into => {
        1 => ['prepare_cdnas'],
      },
    },
    {
      -logic_name => 'prepare_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePreparecDNAs',
      -parameters => {
        prepare_seqs => {
          dest_dir => $self->o('output_path'),
          embl_file => 'embl_' . $taxon_id {$self->o('species')} . '.fa',
          refseq_file => $self->o('refseq_file'),
          gss_file => $self->o('gss_file'),
          killlist_type => 'cdna_update',
          killlist_db => $self->o('killlist_db'),
          polyA_script => $self->o('polyA_script'),
          cdna_file => $self->o('cdna_file_name'),
          species => $self->o('species'),
        },
      },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into => {
        1 => ['compare_cdna_files'],
      },
    },
    {
      # there should probably be a check here to make sure that we get roughly the number of retired sequences we expect
      -logic_name => 'compare_cdna_files',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveComparecDNAfiles',
      -parameters => {
        compare_files => {
          dest_dir => $self->o('output_path'),
          new_cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name').'.clipped',
          old_cdna_file => $self->o('old_cdna_file_name'),
          filtered_file => $self->o('output_path').'/'.$self->o('cdna_file_name').'.filtered',
          retired_file => $self->o('retired_cdnas_file'),
        },
      },
      -max_retry_count => 0,
      -rc_name => 'download',
      #-input_ids => [{}],
      #-wait_for => ['prepare_cdnas'],
      -flow_into => {
        1 => ['list_retired_genes'],
      },
    },
    {
      -logic_name => 'list_retired_genes',
      #-wait_for => ['compare_cdna_files'],
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'for cdna in `cat '.$self->o('retired_cdnas_file').'`; do mysql -h '.$self->o('exonerate_output_db','-host'). ' -u '
               .$self->o('exonerate_output_db','-user').' -p'.$self->o('exonerate_output_db','-pass').' -D '
               .$self->o('exonerate_output_db','-dbname').' -NB -e"DELETE FROM unmapped_object where identifier = \'${cdna}\'"; '
               .'mysql -h '.$self->o('exonerate_output_db','-host'). ' -u '.$self->o('exonerate_output_db','-user').' -p'
               .$self->o('exonerate_output_db','-pass').' -D '.$self->o('exonerate_output_db','-dbname').' -NB  -e"SELECT DISTINCT(g.gene_id)'
               .' FROM gene g left join transcript t on g.gene_id = t.gene_id left join exon_transcript et on t.transcript_id '
               .'= et.transcript_id left join supporting_feature sf on et.exon_id = sf.exon_id left join dna_align_feature daf on ' 
               .'daf.dna_align_feature_id = sf.feature_id where daf.hit_name = \'${cdna}\'" >> '.$self->o('genes_to_delete').'; done'
      },
      #-input_ids => [{}],
      -flow_into => {
        1 => ['delete_retired_genes'],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'delete_retired_genes',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('delete_genes_script').' -h '.$self->o('exonerate_output_db','-host').' -u '.$self->o('exonerate_output_db','-user')   
               .' -p '.$self->o('exonerate_output_db','-pass').' -D '.$self->o('exonerate_output_db','-dbname'). ' -P '.$self->o('exonerate_output_db','-port')
               .' -idfile '.$self->o('genes_to_delete')
      },
      #-input_ids => [{}],
      -max_retry_count => 0,
      #-wait_for => ['list_retired_genes'],
    },
    {
      -logic_name => 'load_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadcDNAs',
      -parameters => {
        cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name').'.filtered',
        species => $self->o('species'),
      },
      -rc_name => 'default',
      -input_ids => [{}],
      -max_retry_count => 0,
      -wait_for => ['delete_retired_genes'],
      -flow_into => {
        1 => ['generate_jobs'],
      },
    },
    {
      -logic_name => 'generate_jobs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        cdna_accession => 1,
        cdna_batch_size => $self->o('cdna_batch_size'),
        cdna_table_name => $self->o('cdna_table_name'),
      },
      -rc_name => 'default',
      #-wait_for => ['load_cdnas'],
      #-input_ids => [{}],
      -flow_into => {
        2 => ['exonerate'],
      },
    },
    {
      -logic_name => 'exonerate',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      #-wait_for => ['populate_production'],
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        logic_name => 'cdna_update',
        module => 'HiveExonerate2Genes_cdna',
        #config_settings => $self->get_config_settings('exonerate_cdna','exonerate'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate')},
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_1.out',
        GENOMICSEQS => $self->o('genome_file'),
        PROGRAM => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -flow_into => {
        2 => [ 'exonerate_second_run' ],
        -1 => ['exonerate_himem'],
      },
      -rc_name => 'exonerate',
      -failed_job_tolerance => 50,
      -batch_size => $self->o('exonerate_batch_size'),
    },
    {
      -logic_name => 'exonerate_himem',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        logic_name => 'cdna_update',
        module => 'HiveExonerate2Genes_cdna',
        #config_settings => $self->get_config_settings('exonerate_cdna','exonerate'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate')},
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_himem.out',
        GENOMICSEQS => $self->o('genome_file'),
        PROGRAM => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
        },
      },
      -flow_into => {
        2 => ['exonerate_second_run'],
      },
      -rc_name => 'exonerate_himem',
      -can_be_empty => 1,
      -failed_job_tolerance => 100,
    },
    {
      # it may be better to carry out this analysis as a series of sql commands like handover_preparation
      # instead of creating a module. I'll leave it like this for now but I'll come back to this
      -logic_name => 'find_missing_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMissingcDNAs',
      -parameters => {
        dest_dir => $self->o('output_path'),
        query_db => $self->o('exonerate_output_db'),
        cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name').'.clipped',
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -wait_for => ['exonerate','exonerate_himem','exonerate_second_run'],
      -input_ids => [{}],
    },
    {
      -logic_name => 'exonerate_second_run',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        logic_name => 'cdna_update',
        module => 'HiveExonerate2Genes_cdna',
        #config_settings => $self->get_config_settings('exonerate_cdna','exonerate_2'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_2')},
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_2.out',
        GENOMICSEQS => $self->o('genome_file'),
        PROGRAM => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -rc_name => 'exonerate_2',
      -can_be_empty => 1,
      -failed_job_tolerance => 100,
    },
    {
      # it may be better to carry out this analysis as a series of sql commands like handover_preparation
      # instead of creating a module. I'll leave it like this for now but I'll come back to this
      -logic_name => 'find_many_hits',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivecDNAManyHits',
      -parameters => {
        dest_dir => $self->o('output_path'),
        query_db => $self->o('exonerate_output_db'),
        file_dir => $self->o('many_hits_dir'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -wait_for => ['exonerate','exonerate_himem','exonerate_second_run'],
      -input_ids => [{}],
      -flow_into => {
        1 => [ 'filter_output' ],
      },
    },
    {
      -logic_name => 'filter_output',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'grep -r "rpp" ' . $self->o('output_path') . '/exonerate_2*out | awk \'{split($0,a,":"); print a[2]}\' >> ' 
               . $self->o('output_path') . '/failed_hits.out;'
               . ' grep -r "only" ' . $self->o('output_path') . '/exonerate_2*out | awk \'{split($0,a,":"); print a[2]}\' >> '
               . $self->o('output_path') . '/failed_hits.out;'
               . ' grep -r "reject" ' . $self->o('output_path') . '/exonerate_2*out | awk \'{split($0,a,":"); print a[2] ": " a[3]}\' >> '
               . $self->o('output_path') . '/failed_hits.out;'
               . ' grep -r "max_coverage" ' . $self->o('output_path') . '/exonerate_2*out | awk \'{split($0,a,":"); print a[1]}\' >> '
               . $self->o('output_path') . '/failed_hits.out'
      },
      #-input_ids => [{}],
      -max_retry_count => 0,
      #-wait_for => ['find_many_hits'],
      -flow_into => {
        1 => [ 'store_unmapped' ],
      },
    },
    {
      -logic_name => 'store_unmapped',
      #-wait_for => ['find_missing_cdnas','filter_output'],
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStoreUnmappedcDNAs',
      -parameters => {
        gss => $self->o('gss_file'),
        seq_file => $self->o('output_path').'/missing_cdnas.fasta',
        user => $self->o('user_w'),
        pass => $self->o('password'),
        host => $self->o('exonerate_output_db_server'),
        port => $self->o('port'),
        dbname => $self->o('exonerate_output_db_name'),
        species => $self->o('species'),
        vertrna => $self->o('output_path').'/embl_'.$taxon_id{$self->o('species')}.'.fa',
        refseq => $self->o('refseq_path').'/'.$self->o('refseq_file'),
        infile => $self->o('output_path').'/failed_hits.out',
        findN_prog => $self->o('findN_script'),
        reasons => $self->o('output_path').'/unmapped_reasons.txt',
        pid => $self->o('exonerate_pid'),
    	cov => $self->o('exonerate_cov'), 
        outdir => $self->o('output_path'),
     },
      #-input_ids => [{}],
      -max_retry_count => 0,
    },
    {
      -logic_name => 'database_compare',
      -wait_for => ['store_unmapped','delete_retired_genes'],
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveComparecDNAdbs',
      -parameters => {
        old_cdna_db => $self->o('old_cdna_db'),
        new_cdna_db => $self->o('exonerate_output_db'),
        output_file => $self->o('output_path').'/comparison.out',
      },
      -input_ids => [{}],
      -failed_job_tolerance => 0,
      -flow_into => {
        1 => [ 'comparison_report','load_xdbids' ],
      },
    },
    {
      -logic_name => 'comparison_report',
      #-wait_for => ['database_compare'],
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
      -parameters => {
        email => $self->o('recipient_email'),
        subject => 'AUTOMATED REPORT: cDNA update database comparison',
        text => 'Please find below the counts for each toplevel seq_region for the current and the previous cDNA updates:',
        file => $self->o('output_path').'/comparison.out',,
      },
      #-input_ids => [{}],
      -failed_job_tolerance => 0,
    },
    {
      -logic_name => 'load_xdbids',
      #-wait_for => ['database_compare'],
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '. $self->o('optimize_script').
               ' -output_path '.$self->o('output_path').'/optimise_daf_paf.1'.
               ' -dbname '.$self->o('exonerate_output_db','-dbname').
               ' -dbhost '.$self->o('exonerate_output_db','-host').
               ' -dbport '.$self->o('exonerate_output_db','-port').
               ' -dbuser '.$self->o('exonerate_output_db','-user').
               ' -dbpass '.$self->o('exonerate_output_db','-pass').
               ' -verbose -clean -no_external_db'
      },
      -rc_name => 'optimise',
      #-input_ids => [{}],
      -max_retry_count => 0,
      -flow_into => {
        1 => [ 'populate_analysis_description' ],
      },
    },
    {
      -logic_name => 'populate_analysis_description',
      #-wait_for => ['load_xdbids'],
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('analysis_desc_script').
               ' -dp '.$self->o('output_path').
               ' -d '.$self->o('exonerate_output_db','-dbname').
               ' -h '.$self->o('exonerate_output_db','-host').
               ' -u '.$self->o('exonerate_output_db','-user').
               ' -p '.$self->o('exonerate_output_db','-pass').
               ' -s homo_sapiens -t cdna'
      },
      #-input_ids => [{}],
      -flow_into => {
        1 => [ 'handover_preparation' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'handover_preparation',
      -wait_for => ['populate_analysis_description'],
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => 'mysql://'.$self->o('exonerate_output_db','-user').':'.
                   $self->o('exonerate_output_db','-pass').'@'.$self->o('exonerate_output_db','-host').
                   ':'.$self->o('exonerate_output_db','-port').'/'.$self->o('exonerate_output_db','-dbname'),
        sql => [
          "DELETE FROM meta WHERE meta_key = 'progress_status'",
          "UPDATE analysis SET db_version = '".$self->o('refseq_version')."', db = 'RefSeq' WHERE logic_name = 'cdna_update'",
          "UPDATE analysis SET program_file= '".$self->o('exonerate_version')."' WHERE logic_name = 'cdna_update'", 
          "DELETE FROM analysis WHERE logic_name != 'cdna_update'",
          "DROP TABLE analysis_description_bak",
          "DROP TABLE attrib_type_bak",
          "DROP TABLE external_db_bak",
          "DROP TABLE misc_set_bak",
          "DROP TABLE unmapped_reason_bak"
        ],
      },
      #-input_ids => [{}],
      -max_retry_count => 0,
    },
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
  };
}

sub resource_classes {
  my $self = shift;
  # Note that this builds resource requests off some of the variables at the top. Works off the idea
  # that the references all get put on one server and the pipe db is on another
  my $pipe_db_server = $self->default_options()->{'pipe_db_server'};
  my $dna_db_server = $self->default_options()->{'dna_db_server'};
  my $exonerate_output_db_server = $self->default_options()->{'exonerate_output_db_server'};
  my $killlist_db_server = $self->default_options()->{'killlist_db_server'};

  my $default_mem = $self->default_options()->{'default_mem'};
  my $exonerate_mem = $self->default_options()->{'exonerate_mem'};
  my $exonerate_retry_mem = $self->default_options()->{'exonerate_retry_mem'};
  my $optimise_mem = $self->default_options()->{'optimise_mem'};
  my $download_mem = $self->default_options()->{'download_mem'};

  my $pipe_db_server_number;
  my $dna_db_server_number;
  my $exonerate_output_db_server_number;
  my $killlist_db_server_number;

  my $num_tokens = $self->default_options()->{'num_tokens'};

  unless($pipe_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "pipe_db_server: ".$pipe_db_server;
  }

  $pipe_db_server_number = $1;

  unless($dna_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "dna_db_server: ".$dna_db_server;
  }

  $dna_db_server_number = $1;

  unless($exonerate_output_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "exonerate_output_db_server: ".$exonerate_output_db_server;
  }

  $exonerate_output_db_server_number = $1;

  unless($killlist_db_server=~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "killlist_db_server: ".$killlist_db_server;
  }

  $killlist_db_server_number = $1;

  unless($num_tokens) {
    die "num_tokens is uninitialised or zero. num_tokens needs to be present in default_options and not zero\n".
        "num_tokens: ".$num_tokens;
  }
  return {
    'default' => { LSF => '-q normal -M'.$default_mem.' -R"select[mem>'.$default_mem.'] '.
                          'rusage[mem='.$default_mem.','.
                          'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'exonerate' => { LSF => '-q normal -M'.$exonerate_mem.' -R"select[mem>'.$exonerate_mem.'] '.
                            'rusage[mem='.$exonerate_mem.','.
                            'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                            'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                            'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'exonerate_himem' => { LSF => '-q normal -M'.$exonerate_retry_mem.
                                  ' -R"select[mem>'.$exonerate_retry_mem.'] '.
                                  'rusage[mem='.$exonerate_retry_mem.','.
                                  'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                                  'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                  'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'exonerate_2' => { LSF => '-q normal -M'.$exonerate_retry_mem.
                              ' -R"select[mem>'.$exonerate_retry_mem.'] '.
                              'rusage[mem='.$exonerate_retry_mem.','.
                              'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'optimise' => { LSF => '-q normal -M'.$optimise_mem.' -R"select[mem>'.$optimise_mem.'] '.
                           'rusage[mem='.$optimise_mem.','.
                           'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.']"' },

    'download' => { LSF => '-q normal -M'.$download_mem.' -R"select[mem>'.$download_mem.'] '.
                           'rusage[mem='.$download_mem.']"' },
  }
}


sub get_config_settings {
  # Shift in the group name (a hash that has a collection of logic name hashes and a default hash)
  # Shift in the logic name of the specific analysis
  my $self = shift;
  my $config_group = shift;
  my $config_logic_name = shift;

  # And additional hash keys will be stored in here
  my @additional_configs = @_;

  # Return a ref to the master hash for the group using the group name
  my $config_group_hash = $self->master_config_settings($config_group);
  unless(defined($config_group_hash)) {
    die "You have asked for a group name in master_config_settings that doesn't exist. Group name:\n".$config_group;
  }
  # Final hash is the hash reference that gets returned. It is important to note that the keys added have
  # priority based on the call to this subroutine, with priority from left to right. Keys assigned to
  # $config_logic_name will have most priority, then keys in any additional hashes, then keys from the
  # default hash. A default hash key will never override a $config_logic_name key
  my $final_hash;

  # Add keys from the logic name hash
  my $config_logic_name_hash = $config_group_hash->{$config_logic_name};
  unless(defined($config_logic_name_hash)) {
    die "You have asked for a logic name hash that doesn't exist in the group you specified.\n".
        "Group name:\n".$config_group."\nLogic name:\n".$config_logic_name;
  }
  $final_hash = $self->add_keys($config_logic_name_hash,$final_hash);

  # Add keys from any additional hashes passed in, keys that are already present will not be overriden
  foreach my $additional_hash (@additional_configs) {
    my $config_additional_hash = $config_group_hash->{$additional_hash};
    $final_hash = $self->add_keys($config_additional_hash,$final_hash);
  }

  # Default is always loaded and has the lowest key value priority
  my $config_default_hash = $config_group_hash->{'Default'};
  $final_hash = $self->add_keys($config_default_hash,$final_hash);

  return($final_hash);
}

sub add_keys {
  my ($self,$hash_to_add,$final_hash) = @_;

  foreach my $key (keys(%$hash_to_add)) {
    unless(exists($final_hash->{$key})) {
      $final_hash->{$key} = $hash_to_add->{$key};
    }
  }
  return($final_hash);
}

sub master_config_settings {
  my ($self,$config_group) = @_;
  my $master_config_settings = {
    exonerate_cdna => {
      Default => {
        IIDREGEXP           => '(\d+):(\d+)',
        OPTIONS             => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE',
        COVERAGE_BY_ALIGNED => 0,
        QUERYTYPE           => 'dna',
        GENOMICSEQS         => $self->o('genome_file'),
        PROGRAM             => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      exonerate => {
        COVERAGE_BY_ALIGNED => 1,
        FILTER => {
          OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
          PARAMETERS => {
            -best_in_genome => 0,
            -coverage => $self->o('exonerate_cov'),
            -percent_id => $self->o('exonerate_pid'),
            -reject_processed_pseudos => 1,
            -verbosity => 1,
          }
        },
        KILL_TYPE => undef,
        USE_KILL_LIST => 0,
        OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 100000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      },
      exonerate_2 => {
        COVERAGE_BY_ALIGNED => 1,
        FILTER => {
          OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
          PARAMETERS => {
            -best_in_genome => 10,
            -coverage => $self->o('exonerate_cov'),
            -percent_id => $self->o('exonerate_pid'),
            -reject_processed_pseudos => 1,
            -verbosity => 1,
          }
        },
        KILL_TYPE => undef,
        OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 100000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      },

      # it looks like the filter option is set, along with killlist, for the second run of the cdna update in the original code
      # I need to therefore look into this - probably best to set in the hash of the analysis in the config settings using killlist_cdna, exonerate etc.
      # so I set the filter => 1 here and then use the appropriate name in the config settings. However as this is only done in the second run I need to
      # think more carefully about how I want to run exonertae - no jobs are failing due to mem. cdna_update_2 analysis seems to be about more than
      # just an increase in mem so I think I need another exonerate run after to do what the rest of cdna_update_2 does 
      killlist_cdna => {
        KILLLISTDB          => $self->o('killlist_db'),
        USE_KILL_LIST       => 1,
        KILL_TYPE           => 'cdna',
        KILL_LIST_FILTER    => {
          -only_mol_type        => 'cdna',
          -user_id              => undef,
          -from_source_species  => undef,
          -before_date          => undef,
          -having_status        => undef,
          -reasons              => [],
          -for_analyses         => [],
          -for_species          => [],
          -for_external_db_ids  => [],
        },
      },
    },
  };
  return($master_config_settings->{$config_group});
}

1;
