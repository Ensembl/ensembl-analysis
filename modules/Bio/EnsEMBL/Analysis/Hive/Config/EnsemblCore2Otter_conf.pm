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

    base_dir => '/hps/nobackup2/production/ensembl/thibaut/loutre',
    output_dir => catdir($self->o('base_dir'), $self->o('species'), $self->o('db_prefix')),
    blast_db_path => catfile($self->o('output_dir'), $self->o('species').'_softmasked_toplevel.fa'),
    password => '',
    user => 'ensadmin',
    user_r => 'ensro',
    pipe_db_host => 'mysql-ens-genebuild-prod-7',
    pipe_db_port => 4533,

    db_prefix => $self->o('assembly_version'),
    pipeline_name => 'loutre_ensembl_'.$self->o('species').'_'.$self->o('db_prefix'),

    assembly_version => 1,
    current_release => 101,
    current_db_host => 'mysql-ens-mirror-1',
    current_db_port => 4240,

    havana_db_host => 'mysql-ens-genebuild-prod-2',
    havana_db_port => 4528,

    ensembl_db_name => $self->o('species').'_core_'.$self->o('current_release').'_'.$self->o('assembly_version'),
    ensembl_db_host => $self->o('current_db_host'),
    ensembl_db_port => $self->o('current_db_port'),

    rnaseq_gene_db_name => $self->o('species').'_rnaseq_'.$self->o('current_release').'_'.$self->o('assembly_version'),
    rnaseq_gene_db_host => $self->o('current_db_host'),
    rnaseq_gene_db_port => $self->o('current_db_port'),

    refine_db_name => $self->o('species').'_rnaseq_'.$self->o('current_release').'_'.$self->o('assembly_version'),
    refine_db_host => $self->o('current_db_host'),
    refine_db_port => $self->o('current_db_port'),

    rnaseq_db_host => $self->o('havana_db_host'),
    rnaseq_db_port => $self->o('havana_db_port'),

    do_uniprot_run => 1,
    uniprot_set => 'havana_human_blast',
    blast_type => 'ncbi',
    protein_entry_loc => '/hps/nobackup2/production/ensembl/genebuild/blastdb/uniprot/uniprot_2019_04/entry_loc',

    meta_pipeline_db_host => $self->o('pipe_db_host'),
    meta_pipeline_db_port => $self->o('pipe_db_port'),

    blast_db_host => $self->o('havana_db_host'),
    blast_db_port => $self->o('havana_db_port'),

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
    gene_number_delimiter => '#',

    meta_pipeline_db_head => "'-dbname' => '#expr(#pipe_db#->{-dbname})expr#', '-host' => '#expr(#pipe_db#->{-host})expr#.ebi.ac.uk', '-port' => '#expr(#pipe_db#->{-port})expr#', '-user' => '".$self->o('user_r')."'",
    meta_pipeline_db_head_rw => "'-dbname' => '#expr(#pipe_db#->{-dbname})expr#', '-host' => '#expr(#pipe_db#->{-host})expr#.ebi.ac.uk', '-port' => '#expr(#pipe_db#->{-port})expr#', '-user' => '#expr(#pipe_db#->{-user})expr#', '-pass' => '#expr(#pipe_db#->{-pass})expr#'",
    meta_core_db_head => "'-dbname' => '#expr(#havana_db#->{-dbname})expr#', '-host' => '#expr(#havana_db#->{-host})expr#.ebi.ac.uk', '-port' => '#expr(#havana_db#->{-port})expr#', '-user' => '".$self->o('user_r')."'",
    meta_intron_db_head => "'-dbname' => '#expr(#rnaseq_db#->{-dbname})expr#', '-host' => '#expr(#rnaseq_db#->{-host})expr#.ebi.ac.uk', '-port' => '#expr(#rnaseq_db#->{-port})expr#', '-user' => '".$self->o('user_r')."'",

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

    refine_db_name => $self->o('species').'_rnaseq_'.$self->o('current_release').'_'.$self->o('assembly_version'), # this would be _refine_ if the database has not been handed over yet
    refine_db_user   => $self->o('user_r'),
    refine_db_pass   => $self->o('password_r'),
    refine_db_driver => $self->o('hive_driver'),

    rnaseq_gene_db_name => $self->o('species').'_rnaseq_'.$self->o('current_release').'_'.$self->o('assembly_version'), # this would be _rnaseq_blast_ if the database has not been handed over yet
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
      $meta_pipeline_db = $self->o('blast_db');
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
             ELSE ['retrieve_repeat_analyses']),
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'retrieve_repeat_analyses',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => ['logic_name'],
        db_conn => $self->o('ensembl_db'),
        inputquery => "SELECT meta_value FROM meta WHERE meta_key = 'repeat.analysis' AND meta_value NOT IN ('trf', 'repeatmask_repeatmodeler')",
        step => 8, # This is a trick to have the multiline result in on array and to have only one input_id. THe value should be higher than the number of expected lines
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        2 => ['dump_softmasked_toplevel'],
      }
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
        repeat_logic_names   => '#_range_list#',
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
          protein_entry_loc => $self->o('protein_entry_loc'),
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
        1 => 'load_feature_status_attributes',
      }
    },
    {
      -logic_name => 'load_feature_status_attributes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        sql => [
          'INSERT INTO gene_attrib SELECT gene_id, attrib_type_id, "PREDICTED" FROM gene, attrib_type WHERE attrib_type.code = "status"',
          'INSERT INTO transcript_attrib SELECT transcript_id, attrib_type_id, "PREDICTED" FROM transcript, attrib_type WHERE attrib_type.code = "status"',
        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'load_evidence_table',
      }
    },
    {
      -logic_name => 'load_evidence_table',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        sql => [
          'INSERT IGNORE INTO evidence SELECT t.transcript_id, xaf.hit_name, "cDNA" FROM transcript t LEFT JOIN transcript_supporting_feature tsf ON t.transcript_id = tsf.transcript_id LEFT JOIN dna_align_feature xaf ON tsf.feature_id = xaf.dna_align_feature_id WHERE tsf.feature_type = "dna_align_feature"',
          'INSERT IGNORE INTO evidence SELECT et.transcript_id, xaf.hit_name, "cDNA" FROM exon_transcript et LEFT JOIN supporting_feature sf ON et.exon_id = sf.exon_id LEFT JOIN dna_align_feature xaf ON sf.feature_id = xaf.dna_align_feature_id WHERE sf.feature_type = "dna_align_feature"',
          'INSERT IGNORE INTO evidence SELECT t.transcript_id, xaf.hit_name, "Protein" FROM transcript t LEFT JOIN transcript_supporting_feature tsf ON t.transcript_id = tsf.transcript_id LEFT JOIN protein_align_feature xaf ON tsf.feature_id = xaf.protein_align_feature_id WHERE tsf.feature_type = "protein_align_feature"',
          'INSERT IGNORE INTO evidence SELECT et.transcript_id, xaf.hit_name, "Protein" FROM exon_transcript et LEFT JOIN supporting_feature sf ON et.exon_id = sf.exon_id LEFT JOIN protein_align_feature xaf ON sf.feature_id = xaf.protein_align_feature_id WHERE sf.feature_type = "protein_align_feature"',
        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'load_feature_gene_name_attributes',
      }
    },
    {
      -logic_name => 'load_feature_gene_name_attributes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        sql => [
          'INSERT INTO gene_attrib SELECT gene_id, attrib_type_id, stable_id FROM gene, attrib_type WHERE attrib_type.code = "name"',
          'UPDATE gene g, gene_attrib ga, xref x, attrib_type at SET ga.value = x.display_label WHERE ga.attrib_type_id = at.attrib_type_id AND ga.gene_id = g.gene_id AND g.display_xref_id = x.xref_id AND at.code = "name"',
        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => 'fan_multiple_gene_names',
        'A->1' => 'load_feature_transcript_name_attributes',
      }
    },
    {
      -logic_name => 'fan_multiple_gene_names',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery => 'SELECT ga.value FROM gene_attrib ga, attrib_type at WHERE ga.attrib_type_id = at.attrib_type_id AND at.code = "name" GROUP BY ga.value HAVING COUNT(*) > 2',
        db_conn => $self->o('havana_db'),
        column_names => ['gnee_name'],
      },
      -rc_name => 'default',
      -flow_into => {
        '2->A' => 'make_gene_names_unique',
        'A->1' => 'load_feature_transcript_name_attributes',
      }
    },
    {
      -logic_name => 'make_gene_names_unique',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        sql => [
          'SET @i:=0',
          'UPDATE gene_attrib ga, attrib_type at SET ga.value = CONCAT(ga.value, "#gene_number_delimiter#", (@i:=@i+1)) WHERE ga.attrib_type_id = at.attrib_type_id AND at.code = "name" AND value = "#gene_name#"',
        ],
        gene_number_delimiter => $self->o('gene_number_delimiter'),
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'load_feature_transcript_name_attributes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('havana_db'),
        sql => [
          'INSERT INTO transcript_attrib SELECT t.transcript_id, at.attrib_type_id, CONCAT(ga.value, "-", SUBSTRING(LPAD(t.transcript_id, 10, 0), -3, 3)) FROM transcript t, gene_attrib ga, attrib_type at WHERE t.gene_id = ga.gene_id AND at.code = "name" AND ga.attrib_type_id = at.attrib_type_id',
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
        meta_core_db_head => $self->o('meta_core_db_head'),
        meta_intron_db_head => $self->o('meta_intron_db_head'),
        pipe_db => $meta_pipeline_db,
        sql => [
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, "pipeline_db_head", "#meta_pipeline_db_head#")',
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, "pipeline_db_rw_head", "#meta_pipeline_db_head_rw#")',
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, "ensembl_core_db_head", "#meta_core_db_head#")',
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, "ensembl_intron_db_head", "#meta_intron_db_head#")',
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
        sql => [
          'INSERT INTO gene_stable_id_pool SELECT SUBSTRING(MAX(stable_id), 8) FROM gene',
          'INSERT INTO transcript_stable_id_pool SELECT SUBSTRING(MAX(stable_id), 8) FROM transcript',
          'INSERT INTO translation_stable_id_pool SELECT SUBSTRING(MAX(stable_id), 8) FROM translation',
          'INSERT INTO exon_stable_id_pool SELECT SUBSTRING(MAX(stable_id), 8) FROM exon',
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
        1 => 'load_otter_misc_set',
      },
      -rc_name => 'default',
    },
    # It needs to be after meta_levels as the script delete all meta key but does not set the misc_feature key
    {
      -logic_name => 'load_otter_misc_set',
      -module     => 'AddVirtualClones',
      -parameters => {
        target_db => $self->o('havana_db'),
        feature_dbs => [$self->o('havana_db')],
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
        refine_db => $self->o('refine_db_name'),
        return_codes_2_branches => {'42' => 2},
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
          'UPDATE dna_align_feature daf, analysis ap, analysis am SET daf.analysis_id = am.analysis_id WHERE daf.analysis_id = ap.analysis_id AND ap.logic_name LIKE "%_daf_plus" AND am.logic_name = REPLACE(ap.logic_name, "_plus", "_minus") AND daf.seq_region_strand = -1',
          'INSERT INTO analysis (logic_name) VALUES("non_canonical_plus"), ("non_canonical_minus")',
          'UPDATE dna_align_feature daf, analysis ap, analysis am SET daf.analysis_id = am.analysis_id WHERE daf.analysis_id = ap.analysis_id AND ap.logic_name LIKE "%_daf_plus" AND am.logic_name = "non_canonical_plus" AND hit_name LIKE "%:non %"',
          'UPDATE dna_align_feature daf, analysis ap, analysis am SET daf.analysis_id = am.analysis_id WHERE daf.analysis_id = ap.analysis_id AND ap.logic_name LIKE "%_daf_minus" AND am.logic_name = "non_canonical_minus" AND hit_name LIKE "%:non %"',

        ],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['fan_rnaseq_config'],
      },
    },
    {
      -logic_name => 'fan_rnaseq_config',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => ['logic_name', 'otter_name', 'zmap_type'],
        db_conn => $self->o('rnaseq_db'),
        inputquery => 'SELECT logic_name, REPLACE(logic_name, "_daf", ""), SUBSTRING_INDEX(logic_name, '_', -1) FROM analysis WHERE logic_name LIKE "%daf%"',
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        '2->A' => ['write_anaysis_config'],
        'A->1' => ['write_swissprot_anaysis_config'],
      }
    },
    {
      -logic_name => 'write_anaysis_config',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'echo "#config_text# > #output_file#',
        output_file => catfile($self->o('output_dir'), 'swissprot.part'),
        species => $self->o('species'),
        config_text => '[#species#.filter.#otter_name#]\n
        analysis=#logic_name#\n
        classification=RNA-seq > Introns\n
        server_script=get_gff/features\n
        metakey=ensembl_rnaseq_intron\n
        description=Introns confirmed by spanning RNASeq reads with number of supporting reads\n
        feature_kind=DnaDnaAlignFeature\n
        zmap_column=rnaseq_introns_#zmap_type#\n
        zmap_style=ensembl_rnaseq_intron\n',
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'write_swissprot_anaysis_config',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'echo "#config_text# > #output_file#',
        output_file => catfile($self->o('output_dir'), 'swissprot.part'),
        species => $self->o('species'),
        config_text => '[#species#.filter.SwissProt]\n
        classification=core\n
        server_script=get_gff/features\n
        analysis=uniprot_sp\n
        description=SwissProt protein tblastn hits\n
        feature_kind=DnaPepAlignFeature\n
        content_type=alignment_feature\n
        resource_bin=core_with_seq\n
        sequence_db=uniprot,uniprot_archive\n
        blixem_data_type=protein-match\n\n',
      },
      -flow_into => {
        1 => ['write_trembl_anaysis_config'],
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'write_trembl_anaysis_config',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'echo "#config_text# > #output_file#',
        output_file => catfile($self->o('output_dir'), 'trembl.part'),
        species => $self->o('species'),
        config_text => '[#species#.filter.TrEMBL]\n
        classification=core\n
        server_script=get_gff/features\n
        analysis=uniprot_tr\n
        description=TrEMBL protein tblastn hits\n
        feature_kind=DnaPepAlignFeature\n
        content_type=alignment_feature\n
        resource_bin=core_with_seq\n
        sequence_db=uniprot,uniprot_archive\n
        blixem_data_type=protein-match\n\n',
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['concat_anaysis_config'],
      },
    },
    {
      -logic_name => 'concat_anaysis_config',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cat "#output_dir#/*.part > #output_file#',
        output_file => catfile($self->o('output_dir'), 'otter_config.ini'),
        output_dir => $self->o('output_dir'),
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

