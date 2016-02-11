=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

my %taxon_id;
$taxon_id {"human"} = 9606;
$taxon_id {"mouse"} = 10090;


sub default_options {
  my ($self) = @_;

  #my %taxon_id;
  #$taxon_id {'human'} = 9606;
  #$taxon_id {'mouse'} = 10090;

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

    'pipeline_name'              => 'human_cdna_84_rerun_2',
    'pipe_db_name'               => 'dm15_human_cdna_84_rerun_2',
    'pipe_db_server'             => 'genebuild11',

    'dna_db_name'                => 'dm15_human_83_copy',
    'dna_db_server'              => 'genebuild12',

    'killlist_db_server'         => 'genebuild6',

    'exonerate_output_db_name'   => 'dm15_human_cdna_exonerate_84_rerun_2',
    'exonerate_output_db_server' => 'genebuild13',

    'old_cdna_db_name'           => 'homo_sapiens_cdna_83_38',
    'old_cdna_db_server'         => 'ens-livemirror',

    'output_path'                => '/lustre/scratch109/ensembl/dm15/hive_humancdna_84_rerun_2/',

    'refseq_path'                => '/data/blastdb/Ensembl/RefSeq_2016_01/',
    'refseq_file'                => 'hs.fna',

    'exonerate_batch_size'       => '50',

    'clone_db_script_path'       => '/nfs/users/nfs_d/dm15/cvs_checkout_head/ensembl-personal/genebuilders/scripts/clone_database.ksh',

    'genome_file'                => '/data/blastdb/Ensembl/human/GRCh38/genome/softmasked_dusted/toplevel.with_nonref_and_GRCh38_p5.no_duplicate.softmasked_dusted.fa',

    'repeat_masking_logic_names' => ['repeatmask_repbase_human'],

    'gss_file'                   => '/nfs/users/nfs_d/dm15/cvs_checkout_head/ensembl-personal/genebuilders/cDNA_update/gss_acc.txt',

    'uniprot_filename'           => '/data/blastdb/Ensembl/uniprot_2015_12/entry_loc',

##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

    'fastasplit_random_path'     => '/software/ensembl/bin/fastasplit_random',
    'killlist_db_name'           => 'gb_kill_list',

    'cdna_file_name'             => 'cdna_update',

    'human_taxon_id'             => '9606',
    'mouse_taxon_id'             => '10090',

    'user_r'                     => 'ensro',
    'user_w'                     => 'ensadmin',
    'password'                   => 'ensembl',
    'port'                       => '3306',

    'cdna_query_dir_name'        => 'cdna_temp',

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

    'polyA_script'               => 'ENSEMBL_REPO_ROOT/ensembl-pipeline/scripts/EST/new_polyA_clipping.pl',

    'driver'                     => 'mysql',
    'num_tokens'                 => 10,
    'user'                       => 'ensro',

    'create_type'                => 'clone',

    'pipeline_db' => {
      -dbname => $self->o('pipe_db_name'),
      -host => $self->o('pipe_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user_w'),
      -pass => $self->o('password'),
      -driver => $self->o('driver'),
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
        source_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        create_type => $self->o('create_type'),
        script_path => $self->o('clone_db_script_path'),
      },
      -rc_name => 'default',
      -input_ids => [{}],
      -max_retry_count => 0,
    },
    {
      -logic_name => 'populate_production',
      -wait_for => ['create_output_db'],
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl $ENSEMBL_REPO_ROOT/ensembl-production/scripts/production_database/populate_production_db_tables.pl'.
               ' -dp '.$self->o('output_path').
               ' -d '.$self->o('exonerate_output_db','-dbname').
               ' -h '.$self->o('exonerate_output_db','-host').
               ' -u '.$self->o('exonerate_output_db','-user').
               ' -p '.$self->o('exonerate_output_db','-pass').
               ' -md ensembl_production -mh ens-staging1 -mu ensro -t external_db -t attrib_type -t misc_set -t unmapped_reason'
      },
      -input_ids => [{}],
      -max_retry_count => 0,
    },
    {
      -logic_name => 'download_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadcDNAFiles',
      -parameters => {
        embl_sequences => {
          output_path => $self->o('output_path'),
          output_file => $self->o('cdna_file_name'),
          species => $self->o('species'),
          #taxon_id => $self->o('human_taxon_id'),
        },
        refseq_sequences => {
          refseq_path => $self->o('refseq_path'),
          refseq_file => $self->o('refseq_file'),
        },
      },  
      -max_retry_count => 0,
      -rc_name => 'download',
      -input_ids => [{}],
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
          #embl_file => $self->o('embl_file'),
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
        1 => ['load_cdnas'],
      },
    },
    {
      -logic_name => 'load_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadcDNAs',
      -parameters => {
        cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name').'.clipped',
        species => $self->o('species'),
      },
      -rc_name => 'default',
      -input_ids => [{}],
      -max_retry_count => 0,
      -wait_for => ['prepare_cdnas'],
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
      -wait_for => ['load_cdnas'],
      -input_ids => [{}],
      -flow_into => {
        1 => ['exonerate'],
      },
    },
##    {
##      -logic_name => 'list_cdnas',
##      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
##      -parameters => {
##        inputcmd => 'mysql -NB -u'.$self->o('user_r').
##                    ' -h'.$self->o('pipeline_db','-host').
##                    ' -D'.$self->o('pipeline_db','-dbname').
##                    ' -P'.$self->o('pipeline_db','-port').
##                    ' -e"select accession from '.$self->o('cdna_table_name').';"',
##        column_names => ['accession'],
##      },
##      -wait_for => ['load_cdnas'],
##      -input_ids => [{}],
##      -flow_into => { 
##        2 => [ 'exonerate' ],
##      },
##      -rc_name => 'default',
##    },
    {
      -logic_name => 'exonerate',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -wait_for => ['populate_production'],
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        logic_name => 'cdna_update',
        module => 'HiveExonerate2Genes_cdna',
        config_settings => $self->get_config_settings('exonerate_cdna','exonerate'),
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_1.out',
      },
      -flow_into => {
        2 => [ 'exonerate_second_run' ],
        -1 => ['exonerate_himem'],
      },
      -rc_name => 'exonerate',
      -failed_job_tolerance => 50,
      -max_retry_count => 1,
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
        config_settings => $self->get_config_settings('exonerate_cdna','exonerate'),
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_himem.out',
      },
      -flow_into => {
        2 => ['exonerate_second_run'],
      },
      -rc_name => 'exonerate_himem',
      -max_retry_count => 1,
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
        config_settings => $self->get_config_settings('exonerate_cdna','exonerate_2'),
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_2.out',
      },
      -rc_name => 'exonerate_2',
      -can_be_empty => 1,
      -max_retry_count => 1,
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
               . ' grep -r "max_coverage" ' . $self->o('output_path') . '/exonerate_2*out | awk \'{split($0,a,":"); print a[2]}\' >> '
               . $self->o('output_path') . '/failed_hits.out'
      },
      -input_ids => [{}],
      -max_retry_count => 0,
      -wait_for => ['find_many_hits'],
    },
    {
      -logic_name => 'store_unmapped',
      -wait_for => ['find_missing_cdnas','filter_output'],
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
        findN_prog => '$ENSEMBL_REPO_ROOT/ensembl-pipeline/scripts/cDNA_update/find_N.pl',
        reasons => $self->o('output_path').'/unmapped_reasons.txt',
        pid => $self->o('exonerate_pid'),
    	cov => $self->o('exonerate_cov'), 
        outdir => $self->o('output_path'),
     },
      -input_ids => [{}],
      -max_retry_count => 0,
    },

##    {
##      -logic_name => 'store_unmapped',
##      -wait_for => ['filter_output'],
##      -wait_for => ['find_missing_cdnas'],
##      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
##      -parameters => {
##        cmd => 'perl $ENSEMBL_REPO_ROOT/ensembl-pipeline/scripts/cDNA_update/store_unmapped_cdnas.pl'
##               . ' -gss ' . $self->o('gss_file')
##               . ' -seq_file ' . $self->o('output_path') . '/missing_cdnas.fasta'
##               . ' -user ' . $self->o('user_w')
##               . ' -pass ' . $self->o('password')
##               . ' -host ' . $self->o('exonerate_output_db_server')
##               . ' -port ' . $self->o('port')
##               . ' -dbname ' . $self->o('exonerate_output_db_name')
##               . ' -species "' . $self->o('species') . '"'
##               . ' -vertrna ' . $self->o('output_path') . '/embl_' . $taxon_id {$self->o('species')} . '.fa'
               #. ' -vertrna ' . $self->o('output_path') . '/' . $self->o('embl_file')
               #. ' -vertrna ' . $self->o('output_path') . '/embl_' . $taxon_id {"$self->o('species')"} . '.fa'
##               . ' -refseq ' . $self->o('refseq_path') . '/' . $self->o('refseq_file')
##               . ' -infile ' . $self->o('output_path') . '/failed_hits.out'
##               . ' -findN_prog $ENSEMBL_REPO_ROOT/ensembl-pipeline/scripts/cDNA_update/find_N.pl'
##               . ' -reasons_file ' . $self->o('output_path') . '/unmapped_reasons.txt'
##      },
##      -input_ids => [ {} ],
##      -max_retry_count => 0,
##    },
    {
      -logic_name => 'database_compare',
      -wait_for => ['store_unmapped'],
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveComparecDNAdbs',
      -parameters => {
        old_cdna_db => $self->o('old_cdna_db'),
        new_cdna_db => $self->o('exonerate_output_db'),
        output_file => $self->o('output_path').'/comparison.out',
      },
      -input_ids => [{}],
      -failed_job_tolerance => 0,
    },
    {
      -logic_name => 'comparison_report',
      -wait_for => ['database_compare'],
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
      -parameters => {
        email => $self->o('recipient_email'),
        subject => 'AUTOMATED REPORT: cDNA update database comparison',
        text => 'Please find below the counts for each toplevel seq_region for the current and the previous cDNA updates:',
        file => $self->o('output_path').'/comparison.out',,
      },
      -input_ids => [{}],
      -failed_job_tolerance => 0,
    },
##    {
##      -logic_name => 'detailed_compare',
##      -wait_for => ['database_compare'],
##      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
##      -parameters => {
##        cmd => 'perl $ENSGBSCRIPTS/load_external_db_ids_and_optimize_af.pl'.
##      },
##      -can_be_empty  => 1,
##      -rc_name => 'default',
##      -failed_job_tolerance => 0,
##    },
    {
      -logic_name => 'load_xdbids',
      -wait_for => ['database_compare'],
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl $ENSGBSCRIPTS/load_external_db_ids_and_optimize_af.pl'.
               ' -output_path '.$self->o('output_path').'/optimise_daf_paf.1'.
               ' -dbname '.$self->o('exonerate_output_db','-dbname').
               ' -dbhost '.$self->o('exonerate_output_db','-host').
               ' -dbport '.$self->o('exonerate_output_db','-port').
               ' -dbuser '.$self->o('exonerate_output_db','-user').
               ' -dbpass '.$self->o('exonerate_output_db','-pass').
               ' -ensgbscripts $ENSGBSCRIPTS -verbose'.
               ' -uniprot_filename '.$self->o('uniprot_filename')
      },
      -rc_name => 'optimise',
      -input_ids => [{}],
      -max_retry_count => 0,
    },
##    {
##      -logic_name => 'populate_production_tables',
##      -wait_for => ['load_xdbids'],
##      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
##      -parameters => {
##        cmd => 'perl $ENSEMBL_REPO_ROOT/ensembl-production/scripts/production_database/populate_production_db_tables.pl'.
##               ' -dp '.$self->o('output_path').
##               ' -d '.$self->o('exonerate_output_db','-dbname').
##               ' -h '.$self->o('exonerate_output_db','-host').
##               ' -u '.$self->o('exonerate_output_db','-user').
##               ' -p '.$self->o('exonerate_output_db','-pass').
##               ' -md ensembl_production -mh ens-staging1 -mu ensro -t external_db -t attrib_type -t misc_set -t unmapped_reason'
##      },
##      -input_ids => [ {} ],
##      -max_retry_count => 0,
##    },
    {
      -logic_name => 'populate_analysis_description',
      -wait_for => ['load_xdbids'],
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl $ENSEMBL_REPO_ROOT/ensembl-production/scripts/production_database/populate_analysis_description.pl'.
               ' -dp '.$self->o('output_path').
               ' -d '.$self->o('exonerate_output_db','-dbname').
               ' -h '.$self->o('exonerate_output_db','-host').
               ' -u '.$self->o('exonerate_output_db','-user').
               ' -p '.$self->o('exonerate_output_db','-pass').
               ' -s homo_sapiens -t cdna'
      },
      -input_ids => [{}],
      -max_retry_count => 0,
    },
    {
      -logic_name => 'handover_preparation',
      -wait_for => ['populate_analysis_description'],
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => 'mysql://'.$self->o('exonerate_output_db','-user').':'.$self->o('exonerate_output_db','-pass').'@'.$self->o('exonerate_output_db','-host').':'.$self->o('exonerate_output_db','-port').'/'.$self->o('exonerate_output_db','-dbname'),
        sql => [
          "DELETE FROM meta WHERE meta_key = 'progress_status'",
          "UPDATE analysis SET db_version = 'refseq_70', db = 'RefSeq' WHERE logic_name = 'cdna_update'",
          "UPDATE analysis SET program_file='exonerate-0.9.0' WHERE logic_name = 'cdna_update'", 
          "DELETE FROM analysis WHERE logic_name != 'cdna_update'",
          "DROP TABLE analysis_description_bak",
          "DROP TABLE attrib_type_bak",
          "DROP TABLE external_db_bak",
          "DROP TABLE misc_set_bak",
          "DROP TABLE unmapped_reason_bak"
        ],
      },
      -input_ids => [{}],
      -max_retry_count => 0,
    },
##    {
##      -logic_name => 'completion_report',
##      -wait_for => ['handover_preparation'],
##      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
##      -parameters => {
##        email => $self->o('recipient_email'),
##        subject => 'AUTOMATED REPORT: cDNA update has completed',
##        text => 'The cDNA update is now complete and the database is ready for handover.',
##      },
##      -input_ids => [{}],
##      -failed_job_tolerance => 0,
##    },
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

    'exonerate' => { LSF => '-q normal -W6:00 -M'.$exonerate_mem.' -R"select[mem>'.$exonerate_mem.'] '.
                            'rusage[mem='.$exonerate_mem.','.
                            'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                            'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                            'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'exonerate_himem' => { LSF => '-q normal -W6:00 -M'.$exonerate_retry_mem.
                                  ' -R"select[mem>'.$exonerate_retry_mem.'] '.
                                  'rusage[mem='.$exonerate_retry_mem.','.
                                  'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                                  'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                  'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'exonerate_2' => { LSF => '-q normal -W6:00 -M'.$exonerate_retry_mem.
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
        OPTIONS             => '',
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
        #OPTIONS => '--model est2genome --forwardcoordinates FALSE --exhaustive FALSE --bestn 0',
        OPTIONS => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
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
        OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --bestn 10 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
        #OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      #OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --bestn 10 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      },

#      exonerate => {
#        FILTER                        => {
#          OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
#          PARAMETERS                  => {
#            -coverage                 => $self->o('exonerate_cov'),
#            -percent_id               => $self->o('exonerate_pid'),
#            -best_in_genome           => 1,
#            -reject_processed_pseudos => 1,
#          },
#        },
#      },
#      exonerate_second_run => {
#        OPTIONS                       => '--model est2genome --forwardcoordinates FALSE --softmasktarget FALSE --maxintron 400000 --bestn 10 --exhaustive FALSE',
#        FILTER                        => {
#          OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
#          PARAMETERS                  => {
#            -verbosity                => 1,
#            -coverage                 => $self->o('exonerate_cov'),
#            -percent_id               => $self->o('exonerate_pid'),
#            -best_in_genome           => 1,
#            -reject_processed_pseudos => 1,
#          },
#        },
#      },

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
