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

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::Config::EnsemblCore2Otter_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::EnsemblCore2Otter_conf;

use strict;
use warnings;

use File::Spec::Functions qw(catfile catdir);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
# This is to be able to use WHEN ELSE
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');


=head2 default_options

 Arg [1]    : None
 Description: Create default hash for this configuration file
 Returntype : Hashref
 Exceptions : None

=cut

sub default_options {
  my ($self) = @_;

  return {
    %{$self->SUPER::default_options()},

    species => 'clupea_harengus',

    output_dir => catdir('/hps/nobackup2/production/ensembl/thibaut/loutre', $self->o('species'), $self->o('db_prefix')),
    blast_db_path        => catfile($self->o('output_dir'), $self->o('species').'_softmasked_toplevel.fa'),
    password => 'ensembl',
    user => 'ensadmin',
    user_r => 'ensro',
    pipe_db_host => 'mysql-ens-genebuild-prod-7',
    pipe_db_port => 4533,

    db_prefix => '1',
    pipeline_name => 'loutre_ensembl_'.$self->o('species').'_'.$self->o('db_prefix'),

    havana_db_host => 'mysql-ens-genebuild-prod-2',
    havana_db_port => 4528,

    ensembl_db_name => 'leanne_clupea_harengus_core_98',
    ensembl_db_host => 'mysql-ens-genebuild-prod-2',
    ensembl_db_port => 4528,

    rnaseq_gene_db_name => 'leanne_clupea_harengus_rnaseq_layer_nr_98',
    rnaseq_gene_db_host => 'mysql-ens-genebuild-prod-3',
    rnaseq_gene_db_port => 4529,

    refine_db_name => 'leanne_clupea_harengus_refine_98',
    refine_db_host => 'mysql-ens-genebuild-prod-3',
    refine_db_port => 4529,

    rnaseq_db_host => 'mysql-ens-genebuild-prod-4',
    rnaseq_db_port => 4530,

    do_uniprot_run => 1,
    uniprot_set => 'havana_human_blast',
    blast_type => 'ncbi',
    full_repbase_logic_name => 'repeatmask_repbase_'.$self->o('species'),

    meta_pipeline_db_host => 'mysql-ens-genebuild-prod-7',
    meta_pipeline_db_port => 4533,

    blast_db_host => 'mysql-ens-genebuild-prod-3',
    blast_db_port => 4529,

    killlist_db_host => 'mysql-ens-genebuild-prod-6',
    killlist_db_port => 4532,

    hive_init_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'init_pipeline.pl'),
    hive_beekeeper_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'beekeeper.pl'),
    hive_uniprot_config => 'Bio::EnsEMBL::Analysis::Hive::Config::UniProt_BLAST_genome_conf',
    meta_hive_capacity => 100,

    email_suffix => 'ebi.ac.uk',
    havana_group_name => 'Havana',
    havana_email_prefix => 'vega',
    havana_email_suffix => $self->o('email_suffix'),
    genebuild_group_name => 'GeneBuild',
    genebuild_email_prefix => 'ensembl-genebuild',
    genebuild_email_suffix => $self->o('email_suffix'),

    otter_archive_rank => 100,
    otter_archive_coord_system_id => 10010,

    meta_pipeline_db_head => "'-dbname' => '#expr(#pipe_db#->{-dbname})expr#', '-host' => '#expr(#pipe_db#->{-host})expr#.ebi.ac.uk', '-port' => '#expr(#pipe_db#->{-port})expr#', '-user' => '".$self->o('user_r')."'",
    meta_pipeline_db_head_rw => "'-dbname' => '#expr(#pipe_db#->{-dbname})expr#', '-host' => '#expr(#pipe_db#->{-host})expr#.ebi.ac.uk', '-port' => '#expr(#pipe_db#->{-port})expr#', '-user' => '#expr(#pipe_db#->{-user})expr#', '-pass' => '#expr(#pipe_db#->{-pass})expr#'",

    ensembl_dir => catdir($self->o('enscode_root_dir'), 'ensembl'),
    ensembl_scripts => catdir($self->o('ensembl_dir'), 'misc-scripts'),
    otter_schema => catfile($self->o('enscode_root_dir'), 'ensembl-otter', 'sql', 'otter_schema.sql'),

    havana_db_name => $self->o('dbowner').'_'.$self->o('species').'_havana_'.$self->o('db_prefix'),
    havana_db_user => $self->o('user'),
    havana_db_password => $self->o('password'),
    havana_db_driver => $self->o('hive_driver'),

    ensembl_db_user => $self->o('user_r'),
    ensembl_db_password => $self->o('password_r'),
    ensembl_db_driver => $self->o('hive_driver'),

    meta_pipeline_db_name => $self->o('dbowner').'_'.$self->o('species').'_uniprot_'.$self->o('db_prefix').'_pipe',
    meta_pipeline_db_user => $self->o('user'),
    meta_pipeline_db_password => $self->o('password'),
    meta_pipeline_db_driver => $self->o('hive_driver'),

    blast_db_name => $self->o('dbowner').'_'.$self->o('species').'_uniprot_blast_'.$self->o('db_prefix'),
    blast_db_user   => $self->o('user'),
    blast_db_pass   => $self->o('password'),
    blast_db_driver => $self->o('hive_driver'),

    refine_db_name => $self->o('dbowner').'_'.$self->o('species').'_refine_'.$self->o('db_prefix'),
    refine_db_user   => $self->o('user_r'),
    refine_db_pass   => $self->o('password_r'),
    refine_db_driver => $self->o('hive_driver'),

    rnaseq_gene_db_name => $self->o('dbowner').'_'.$self->o('species').'_rnaseq_blast_'.$self->o('db_prefix'),
    rnaseq_gene_db_user   => $self->o('user_r'),
    rnaseq_gene_db_pass   => $self->o('password_r'),
    rnaseq_gene_db_driver => $self->o('hive_driver'),

    rnaseq_db_name => $self->o('dbowner').'_'.$self->o('species').'_rnaseq_'.$self->o('db_prefix'),
    rnaseq_db_user   => $self->o('user'),
    rnaseq_db_pass   => $self->o('password'),
    rnaseq_db_driver => $self->o('hive_driver'),

    killlist_db_name => 'gb_kill_list',
    killlist_db_user   => $self->o('user_r'),
    killlist_db_pass   => $self->o('password_r'),
    killlist_db_driver => $self->o('hive_driver'),

    havana_db => {
                 -dbname => $self->o('havana_db_name'),
                 -host   => $self->o('havana_db_host'),
                 -port   => $self->o('havana_db_port'),
                 -user   => $self->o('havana_db_user'),
                 -pass   => $self->o('havana_db_password'),
                 -driver => $self->o('havana_db_driver'),
                 },

    ensembl_db => {
                 -dbname => $self->o('ensembl_db_name'),
                 -host   => $self->o('ensembl_db_host'),
                 -port   => $self->o('ensembl_db_port'),
                 -user   => $self->o('ensembl_db_user'),
                 -pass   => $self->o('ensembl_db_password'),
                 -driver => $self->o('ensembl_db_driver'),
                 },

    meta_pipeline_db => {
                 -dbname => $self->o('meta_pipeline_db_name'),
                 -host   => $self->o('meta_pipeline_db_host'),
                 -port   => $self->o('meta_pipeline_db_port'),
                 -user   => $self->o('meta_pipeline_db_user'),
                 -pass   => $self->o('meta_pipeline_db_password'),
                 -driver => $self->o('meta_pipeline_db_driver'),
                 },

    blast_db => {
                 -dbname => $self->o('blast_db_name'),
                 -host   => $self->o('blast_db_host'),
                 -port   => $self->o('blast_db_port'),
                 -user   => $self->o('blast_db_user'),
                 -pass   => $self->o('blast_db_pass'),
                 -driver => $self->o('blast_db_driver'),
    },

    refine_db => {
                 -dbname => $self->o('refine_db_name'),
                 -host   => $self->o('refine_db_host'),
                 -port   => $self->o('refine_db_port'),
                 -user   => $self->o('refine_db_user'),
                 -pass   => $self->o('refine_db_pass'),
                 -driver => $self->o('refine_db_driver'),
    },

    rnaseq_gene_db => {
                 -dbname => $self->o('rnaseq_gene_db_name'),
                 -host   => $self->o('rnaseq_gene_db_host'),
                 -port   => $self->o('rnaseq_gene_db_port'),
                 -user   => $self->o('rnaseq_gene_db_user'),
                 -pass   => $self->o('rnaseq_gene_db_pass'),
                 -driver => $self->o('rnaseq_gene_db_driver'),
    },

    rnaseq_db => {
                 -dbname => $self->o('rnaseq_db_name'),
                 -host   => $self->o('rnaseq_db_host'),
                 -port   => $self->o('rnaseq_db_port'),
                 -user   => $self->o('rnaseq_db_user'),
                 -pass   => $self->o('rnaseq_db_pass'),
                 -driver => $self->o('rnaseq_db_driver'),
    },

    killlist_db => {
                 -dbname => $self->o('killlist_db_name'),
                 -host   => $self->o('killlist_db_host'),
                 -port   => $self->o('killlist_db_port'),
                 -user   => $self->o('killlist_db_user'),
                 -pass   => $self->o('killlist_db_pass'),
                 -driver => $self->o('killlist_db_driver'),
    },

    databases_to_delete => ['havana_db', 'rnaseq_db'],
  };
}


sub pipeline_analyses {
  my ($self) = @_;

  my $pipedb;
  my $url;
  my $guiurl;
  my $meta_pipeline_db = $self->o('meta_pipeline_db');
  if ($self->_is_second_pass('do_uniprot_run')) {
    if ($self->o('do_uniprot_run')) {
      ($pipedb, $url, $guiurl) = $self->get_meta_db_information($meta_pipeline_db);
    }
    else {
      $meta_pipeline_db = $self->o('havana_db');
    }
  }
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


  my @analyses = (
    {
      -logic_name => 'create_working_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mkdir -p #output_dir#',
        output_dir => $self->o('output_dir'),
      },
      -input_ids => [{}],
      -flow_into => {
        1 => ['create_havana_db'],
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'create_havana_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('ensembl_db'),
        target_db => $self->o('havana_db'),
        create_type => 'copy',
        _lock_tables => 'false',
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['fan_uniprot_run', 'patch_havana_db'],
      },
    },
    {
      -logic_name => 'fan_uniprot_run',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ #do_uniprot_run# -eq 0 ];then exit 42;fi',
        do_uniprot_run => $self->o('do_uniprot_run'),
        return_codes_2_branches => {'42' => 2},
        blast_db_path => $self->o('blast_db_path'),
      },
      -flow_into  => {
        1 => WHEN ('-e #blast_db_path#' => ['create_uniprot_pipeline_job'],
             ELSE ['dump_softmasked_toplevel']),
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'dump_softmasked_toplevel',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
      -parameters => {
        coord_system_name    => 'toplevel',
        target_db            => $self->o('havana_db'),
        output_path          => $self->o('output_dir'),
        enscode_root_dir     => $self->o('enscode_root_dir'),
        species_name         => $self->o('species'),
        repeat_logic_names   => ['dust', $self->o('full_repbase_logic_name')],
      },
      -flow_into => {
        1 => ['format_softmasked_toplevel'],
      },
      -rc_name    => '4GB',
    },
    {
      -logic_name => 'format_softmasked_toplevel',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        'cmd'    => 'if [ "'.$self->o('blast_type').'" = "ncbi" ]; then convert2blastmask -in '.$self->o('blast_db_path').' -parse_seqids -masking_algorithm repeatmasker -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out '.$self->o('blast_db_path').'.asnb;makeblastdb -in '.$self->o('blast_db_path').' -dbtype nucl -parse_seqids -mask_data '.$self->o('blast_db_path').'.asnb -title "'.$self->o('species').'"; else xdformat -n '.$self->o('blast_db_path').';fi',
      },
      -rc_name    => '4GB',
      -flow_into => {
        1 => ['create_uniprot_pipeline_job'],
      },
    },
    {
      -logic_name => 'create_uniprot_pipeline_job',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => ['meta_pipeline_db', 'ehive_url', 'pipeline_name', 'external_link'],
        inputlist => [[$self->o('meta_pipeline_db'), $url, 'havana_uniprot_blast_'.$self->o('species'), $guiurl]],
        guihive_host => $self->o('guihive_host'),
        guihive_port => $self->o('guihive_port'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        2 => ['init_uniprot_pipeline'],
      }
    },
    {
      -logic_name => 'init_uniprot_pipeline',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit',
      -parameters => {
        hive_config => $self->o('hive_uniprot_config'),
        databases => ['blast_db', 'killlist_db', 'dna_db'],
        blast_db => $self->o('blast_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('havana_db'),
        enscode_root_dir => $self->o('enscode_root_dir'),
        extra_parameters => {
          blast_db_path => $self->o('blast_db_path'),
          output_path => $self->o('output_dir'),
          user_r => $self->o('user_r'),
          meta_hive_capacity => $self->o('meta_hive_capacity'),
          uniprot_set => $self->o('uniprot_set'),
        },
      },
      -rc_name      => 'default',
      -max_retry_count => 1,
      -flow_into => {
        1 => ['run_uniprot_pipeline'],
      },
    },
    {
      -logic_name => 'run_uniprot_pipeline',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun',
      -parameters => {
        beekeeper_script => $self->o('hive_beekeeper_script'),
        ehive_url => $url,
      },
      -rc_name      => 'default',
      -max_retry_count => 1,
    },
    {
      -logic_name => 'patch_havana_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        ensembl_scripts => $self->o('ensembl_scripts'),
        cmd => 'perl '.catfile('#ensembl_scripts#', 'schema_patcher.pl')
          .' --type core --nointeractive'
          .' --host #expr(#db_conn#->{-host})expr#'
          .' --port #expr(#db_conn#->{-port})expr#'
          .' --database #expr(#db_conn#->{-dbname})expr#'
          .' --user #expr(#db_conn#->{-user})expr#'
          .' --pass #expr(#db_conn#->{-pass})expr#',
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['load_otter_schema'],
      },
    },
    {
      -logic_name => 'load_otter_schema',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        input_file => $self->o('otter_schema'),
      },
      -flow_into => {
        1 => ['load_otter_archive_coord_system'],
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'load_otter_archive_coord_system',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        sql => ['INSERT INTO coord_system (coord_system_id, species_id, name, version, rank) VALUES (#coord_system_id#, 1, "chromosome", "OtterArchive", #rank#)'],
        rank => $self->o('otter_archive_rank'),
        coord_system_id => $self->o('otter_archive_coord_system_id'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'load_author_information',
      }
    },
    {
      -logic_name => 'load_author_information',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        havana_group_name => $self->o('havana_group_name'),
        havana_email => $self->o('havana_email_prefix').'@'.$self->o('havana_email_suffix'),
        genebuild_group_name => $self->o('genebuild_group_name'),
        genebuild_email => $self->o('genebuild_email_prefix').'@'.$self->o('genebuild_email_suffix'),
        sql => [
          'INSERT INTO author_group (group_name, group_email) VALUES ("#havana_group_name#", "#havana_email#")',
          'INSERT INTO author_group (group_name, group_email) VALUES ("#genebuild_group_name#", "#genebuild_email#")',
          'INSERT INTO author (author_email, author_name, group_id) VALUES ("ensembl-genebuild@ebi.ac.uk", "genebuild", (SELECT group_id FROM author_group WHERE group_name = "GeneBuild"))',
        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'load_seq_region_attributes',
      }
    },
    {
      -logic_name => 'load_seq_region_attributes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        sql => [
          'INSERT INTO seq_region_attrib (seq_region_id, attrib_type_id, value) SELECT sr.seq_region_id, at1.attrib_type_id, 1 FROM seq_region sr LEFT JOIN seq_region_attrib sra ON sr.seq_region_id = sra.seq_region_id LEFT JOIN attrib_type at ON sra.attrib_type_id = at.attrib_type_id, attrib_type at1 WHERE at.code = "toplevel" AND at1.code = "write_access"',
          'INSERT INTO seq_region_attrib (seq_region_id, attrib_type_id, value) SELECT sr.seq_region_id, at1.attrib_type_id, 0 FROM seq_region sr LEFT JOIN seq_region_attrib sra ON sr.seq_region_id = sra.seq_region_id LEFT JOIN attrib_type at ON sra.attrib_type_id = at.attrib_type_id, attrib_type at1 WHERE at.code = "toplevel" AND at1.code = "hidden"',
          'INSERT INTO seq_region_attrib (seq_region_id, attrib_type_id, value) SELECT sr.seq_region_id, at1.attrib_type_id, sr.name FROM seq_region sr LEFT JOIN seq_region_attrib sra ON sr.seq_region_id = sra.seq_region_id LEFT JOIN attrib_type at ON sra.attrib_type_id = at.attrib_type_id, attrib_type at1 WHERE at.code = "toplevel" AND at1.code = "description"',
          'UPDATE seq_region sr, seq_region_attrib sra, seq_region_synonym srs, external_db ed, attrib_type at SET sra.value = srs.synonym WHERE sr.seq_region_id = sra.seq_region_id AND sr.seq_region_id = srs.seq_region_id AND srs.external_db_id = ed.external_db_id AND sra.attrib_type_id = at.attrib_type_id AND at.code = "description" AND ed.db_name = "INSDC"',
        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'load_feature_attributes',
      }
    },
    {
      -logic_name => 'load_feature_attributes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        sql => [
          'INSERT INTO gene_attrib (gene_id, attrib_type_id, value) SELECT gene_id, attrib_type_id, "PUTATIVE" FROM gene, attrib_type WHERE attrib_type.code = "status"',
          'INSERT INTO transcript_attrib (transcript_id, attrib_type_id, value) SELECT transcript_id, attrib_type_id, "PUTATIVE" FROM transcript, attrib_type WHERE attrib_type.code = "status"',
          'INSERT INTO transcript_attrib SELECT transcript_id, 4, CONCAT(gene.gene_id, ".", gene.version, "-", SUBSTR(transcript.transcript_id, -3)) FROM transcript, gene WHERE transcript.gene_id = gene.gene_id',
        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'load_pipeline_database_information',
      }
    },
    {
      -logic_name => 'load_pipeline_database_information',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        meta_pipeline_db_head => $self->o('meta_pipeline_db_head'),
        meta_pipeline_db_head_rw => $self->o('meta_pipeline_db_head_rw'),
        pipe_db => $meta_pipeline_db,
        sql => [
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, "pipeline_db_head", "#meta_pipeline_db_head#")',
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, "pipeline_db_rw_head", "#meta_pipeline_db_head_rw#")',
        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'set_pool_id_max_id',
      }
    },
    {
      -logic_name => 'set_pool_id_max_id',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        pipe_db => $meta_pipeline_db,
        sql => [
          'SET @MAXID=(SELECT SUBSTRING(MAX(stable_id), 8) FROM gene)',
          'INSERT INTO gene_stable_id_pool VALUES(@MAXID)',
          'SET @MAXID=(SELECT SUBSTRING(MAX(stable_id), 8) FROM transcript)',
          'INSERT INTO transcript_stable_id_pool VALUES(@MAXID)',
          'SET @MAXID=(SELECT SUBSTRING(MAX(stable_id), 8) FROM translation)',
          'INSERT INTO translation_stable_id_pool VALUES(@MAXID)',
          'SET @MAXID=(SELECT SUBSTRING(MAX(stable_id), 8) FROM exon)',
          'INSERT INTO exon_stable_id_pool VALUES(@MAXID)',
        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'update_meta_levels',
      }
    },
    {
      -logic_name => 'update_meta_levels',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        ensembl_scripts => $self->o('ensembl_scripts'),
        cmd => 'perl '.catfile('#ensembl_scripts#', 'meta_levels.pl')
          .' --host #expr(#db_conn#->{-host})expr#'
          .' --port #expr(#db_conn#->{-port})expr#'
          .' --dbname #expr(#db_conn#->{-dbname})expr#'
          .' --user #expr(#db_conn#->{-user})expr#'
          .' --pass #expr(#db_conn#->{-pass})expr#',
      },
      -flow_into => {
        1 => 'update_meta_coords',
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'update_meta_coords',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        ensembl_scripts => $self->o('ensembl_scripts'),
        cmd => 'perl '.catfile('#ensembl_scripts#', 'meta_coord', 'update_meta_coord.pl')
          .' --host #expr(#db_conn#->{-host})expr#'
          .' --port #expr(#db_conn#->{-port})expr#'
          .' --dbname #expr(#db_conn#->{-dbname})expr#'
          .' --user #expr(#db_conn#->{-user})expr#'
          .' --pass #expr(#db_conn#->{-pass})expr#;'
          .' rm #expr(#db_conn#->{-dbname})expr#*.meta_coord.backup',
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'create_rnaseq_db',
      },
    },
    {
      -logic_name => 'create_rnaseq_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('rnaseq_gene_db'),
        target_db => $self->o('rnaseq_db'),
        create_type => 'copy',
        _lock_tables => 'false',
      },
      -rc_name => 'default',
      -flow_into => {
        '1->A' => ['fan_rnaseq_intron'],
        'A->1' => ['update_rnaseq_tables'],
      },
    },
    {
      -logic_name => 'fan_rnaseq_intron',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ "#rnaseq_gene_db#" = "#refine_db#" ]; then exit 42; fi',
        rnaseq_gene_db => $self->o('rnaseq_gene_db_name'),
        refine_db => $self->o('refine_db'),
      },
      -flow_into => {
        1 => ['add_rnaseq_intron'],
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'add_rnaseq_intron',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
      -parameters => {
        db_conn => $self->o('refine_db'),
        output_db => $self->o('rnaseq_db'),
        table_list => [
          'dna_align_feature',
        ],
      },
      -rc_name    => 'default',
    },
    {
      -logic_name => 'update_rnaseq_tables',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_db'),
        sql => [
          'DELETE FROM daf USING dna_align_feature daf, analysis a WHERE daf.analysis_id = a.analysis_id AND a.logic_name = "rough_transcripts"',
          'UPDATE analysis SET logic_name = CONCAT(logic_name, "_plus") WHERE logic_name LIKE "%_daf"',
          'INSERT INTO analysis (logic_name, module) SELECT REPLACE(logic_name, "_plus", "_minus"), module FROM analysis WHERE logic_name LIKE "%_daf_plus"',
          'UPDATE dna_align_feature daf, analysis ap, analysis am SET daf.analysis_id = am.analysis_id WHERE daf.analysis_id = ap.analysis_id AND ap.logic_name LIKE "%_daf_plus" AND am.logic_name = REPLACE(ap.logic_name, "_plus", "_minus")',
          'INSERT INTO analysis (logic_name) VALUES("non_canonical_plus"), ("non_canonical_minus")',
          'UPDATE dna_align_feature daf, analysis ap, analysis am SET daf.analysis_id = am.analysis_id WHERE daf.analysis_id = ap.analysis_id AND ap.logic_name LIKE "%_daf_plus" AND am.logic_name = "non_canonical_plus" AND hit_name LIKE "%:non %"',
          'UPDATE dna_align_feature daf, analysis ap, analysis am SET daf.analysis_id = am.analysis_id WHERE daf.analysis_id = ap.analysis_id AND ap.logic_name LIKE "%_daf_minus" AND am.logic_name = "non_canonical_minus" AND hit_name LIKE "%:non %"',

        ],
      },
      -rc_name => 'default',
    },

  );

  foreach my $analysis (@analyses) {
    $analysis->{'-max_retry_count'} = 0 unless (exists $analysis->{'-max_retry_count'});
  }

  return \@analyses;
}


sub resource_classes {
  my ($self) = @_;

  return {
    'default' => { LSF => '-q production-rh74 -M 1000 -R"select[mem>1000] rusage[mem=1000]"'},
    '4GB' => { LSF => '-q production-rh74 -M 4000 -R"select[mem>4000] rusage[mem=4000]"'},
    '16GB' => { LSF => '-q production-rh74 -M 16000 -R"select[mem>16000] rusage[mem=16000]"'},
  };
}

1;

