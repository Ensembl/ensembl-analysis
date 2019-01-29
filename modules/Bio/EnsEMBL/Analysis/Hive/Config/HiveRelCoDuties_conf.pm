=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Hive::Config::HiveRelCoDuties_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::HiveRelCoDuties_conf;

use strict;
use warnings;

use 5.014000;
use File::Spec::Functions qw(catfile catdir);
use MIME::Base64 qw(encode_base64url);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(first_upper_case);
use parent ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

=head2 default_options

 Description: It returns a hashref containing the default options for HiveGeneric_conf
 Returntype : Hashref
 Exceptions : None


=cut

sub default_options {
    my ($self) = @_;

    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },
        use_jira => 1,
        jira_user => '' || $ENV{JIRA_USER},
        jira_password => '' || $ENV{JIRA_PASS},
        ensembl_release => $ENV{ENSEMBL_RELEASE}, # Use it on the commandline: -ensembl_release XX
        #working_dir => '', # Use it on the command line or uncomment: -working_dir /path/to/my/scratch/relco_duties_87
        #user => '', # Use it on the command line or uncomment: -user mysql_rw_user
        #password => '', # Use it on the command line or uncomment: -password mysql_rw_password
        #user_r => '', # Use it on the command line or uncomment: -user_r mysql_ro_user
        pass_r => undef, # (optional) Use it on the command line: -password mysql_ro_password
        human_gencode_version => '27', # Use it on the command line: -human_gencode_version XX
        mouse_gencode_version => 'M15', # Use it on the command line: -password mysql_rw_password
        do_human => 0,
        do_mouse => 1,
        do_pig => 0,
        do_chimpanzee => 1,
        do_rat => 0,
        do_zebrafish => 0,
#################
#        Everything below should not need modification
#################
        pipeline_name => 'relco_duties_'.$self->o('ensembl_release'),
        ensembl_analysis_dir => catdir($self->o('enscode_root_dir'), 'ensembl-analysis'),
        production_dir => catdir($self->o('enscode_root_dir'), 'ensembl-production'),
        human_cs_version => '38',
        mouse_cs_version => '38',
        human_cs_name => 'GRCh'.$self->o('human_cs_version'),
        mouse_cs_name => 'GRCm'.$self->o('mouse_cs_version'),
        human_alias => 'homo_sapiens',
        mouse_alias => 'mus_musculus',
        binary_base => catdir($ENV{LINUXBREW_HOME}, 'bin'),
        samtools => catfile($self->o('binary_base'), 'samtools'), # It might be better to give the absolute path, but I think it's ok
        webdev_nfs => '/nfs/production/panda/ensembl/production/ensemblftp/data_files',

        enscode_root_dir => $ENV{ENSCODE},
        ensembl_analysis_dir => catdir($self->o('enscode_root_dir'), 'ensembl-analysis'),
        production_dir => catdir($self->o('enscode_root_dir'), 'ensembl-production'),
        check_datafiles_script => catfile($self->o('production_dir'), 'scripts', 'datafiles', 'check_datafiles.pl'),
        compare_gencode_refseq_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'compare_gencode_refseq_genes.pl'),
        overlap_gencode_refseq_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'refseq', 'ens_refseq_comparison.pl'),
        label_gencode_basic_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'label_gencode_basic_transcripts.pl'),
        update_gencode_version => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'create_gencode_web_data_patch.pl'),

        query_delete_attributes => ['DELETE FROM transcript_attrib WHERE attrib_type_id = 510'],

        staging1_db => {
            -host   => 'mysql-ens-sta-1',
            -port   => 4519,
            -user   => $self->o('user_r'),
            -pass   => $self->o('pass_r'),
            -driver => $self->o('hive_driver'),
        },
        pre_staging1_db => {
            -host   => 'mysql-ens-vertannot-staging',
            -port   => 4573,
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('hive_driver'),
        },
        human_ensembl_db => {
            -dbname => $self->o('human_alias').'_core_'.$self->o('ensembl_release').'_'.$self->o('human_cs_version'),
            -host   => $self->o('pre_staging1_db', '-host'),
            -port   => $self->o('pre_staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver   => $self->o('pre_staging1_db', '-driver'),
        },
        mouse_ensembl_db => {
            -dbname => $self->o('mouse_alias').'_core_'.$self->o('ensembl_release').'_'.$self->o('mouse_cs_version'),
            -host   => $self->o('pre_staging1_db', '-host'),
            -port   => $self->o('pre_staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('pre_staging1_db', '-driver'),
        },
        human_refseq_db => {
            -dbname => $self->o('human_alias').'_otherfeatures_'.$self->o('ensembl_release').'_'.$self->o('human_cs_version'),
            -host   => $self->o('pre_staging1_db', '-host'),
            -port   => $self->o('pre_staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('pre_staging1_db', '-driver'),
        },
        mouse_refseq_db => {
            -dbname => $self->o('mouse_alias').'_otherfeatures_'.$self->o('ensembl_release').'_'.$self->o('mouse_cs_version'),
            -host   => $self->o('pre_staging1_db', '-host'),
            -port   => $self->o('pre_staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('pre_staging1_db', '-driver'),
        },
        production_db => {
            -dbname => 'ensembl_production_'.$self->o('ensembl_release'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user_r'),
            -pass   => $self->o('pass_r'),
            -driver => $self->o('staging1_db', '-driver'),
        },
    };
}


=head2 pipeline_analyses

 Arg [1]    : None
 Description: Returns a hashref containing the analyses to run
 Returntype : Hashref
 Exceptions : None

=cut

sub pipeline_analyses {
  my ($self) = @_;

  my %jira_tickets;
  my @default_tickets;
  my $user = $self->o('jira_user');
  my $password = $self->o('jira_password');
  if ($self->o('use_jira') eq '1') {
    %jira_tickets = (
      version_ticket_name => 'Update release ENSEMBL_RELEASE in genebuild.sh to '.$self->o('ensembl_release'),
      release_ticket_name  => 'Genebuild Relco release '.$self->o('ensembl_release'),
      species_ticket_name  => 'GENCODE scripts release '.$self->o('ensembl_release'),
      datafile_ticket_name => 'Check data_file scripts release '.$self->o('ensembl_release'),
      base64 => encode_base64url($user.':'.$password),
    );
    @default_tickets = (
      'Update ENSEMBL_RELEASE to '.$self->o('ensembl_release').' in genebuild.sh',
      'Delete databases from previous release from mysql-ens-vertannot-staging',
      'Clean old healthchecks before handover',
      'Human cDNA update release '.$self->o('ensembl_release'),
      'Mouse cDNA update release '.$self->o('ensembl_release'),
    );
    foreach my $species ('human', 'mouse') {
      if ($self->o('do_'.$species)) {
        push(@default_tickets, $species.' CCDS update release '.$self->o('ensembl_release'));
      }
    }
  }
  my @refseq_list;
  foreach my $species ('human', 'mouse') {
    if ($self->o('do_'.$species)) {
      push(@refseq_list, [$species, $self->o($species.'_ensembl_db'), $self->o($species.'_cs_name'), $self->o($species.'_refseq_db'), "refseq_${species}_import"]);
    }
  }
  my @analysis = (
    {
      -logic_name => 'create_release_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'create',
        type => 'Task',
        ticket_name => $jira_tickets{release_ticket_name},
        base64 => $jira_tickets{base64},
        assignee => 'self',
        params => {description => $self->pipeline_url},
      },
      -input_ids => [{working_dir => $self->o('working_dir')}],
      -flow_into => {
          1 => ['create_genebuld_env_jira'],
      },
    },
    {
      -logic_name => 'create_genebuld_env_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'create',
        type => 'Sub-task',
        ticket_name => $jira_tickets{version_ticket_name},
        parent => $jira_tickets{release_ticket_name},
        base64 => $jira_tickets{base64},
        assignee => 'self',
      },
      -flow_into => {
          1 => ['create_jira_tickets'],
      },
    },
    {
      -logic_name => 'create_jira_tickets',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
        inputlist => \@default_tickets,
        column_names => ['ticket_name'],
      },
      -flow_into => {
        '2->A' => ['create_default_jira'],
        'A->1' => ['create_working_directory'],
      },
    },
    {
      -logic_name => 'create_default_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'create',
        type => 'Sub-task',
        base64 => $jira_tickets{base64},
        parent => $jira_tickets{release_ticket_name},
      },
    },

    {
      -logic_name => 'create_working_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'mkdir #working_dir#',
      },
      -flow_into => {
          '1' => ['create_species_jira', 'create_datafile_jira'],
      },
    },

    {
      -logic_name => 'create_species_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'create',
        type => 'Sub-task',
        ticket_name => $jira_tickets{species_ticket_name},
        base64 => $jira_tickets{base64},
        parent => $jira_tickets{release_ticket_name},
        assignee => 'self',
      },
      -flow_into => {
        1 => ['create_species_input'],
      },
    },

    {
      -logic_name => 'create_species_input',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputlist => \@refseq_list,
          column_names => ['species', 'target_db', 'cs_name', 'refseq_db', 'refseq_logicname'],
      },
      -flow_into => {
          '2' => ['dump_db_backup'],
      },
    },

    {
      -logic_name => 'dump_db_backup',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
      -rc_name    => 'default',
      -parameters => {
          src_db_conn => '#target_db#',
          output_file => catfile($self->o('working_dir'), '#species#_#cs_name#.sql'),
      },
      -flow_into => {
          1 => ['comment_species_jira'],
      },
    },

    {
      -logic_name => 'comment_species_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'comment',
        type => 'Sub-task',
        ticket_name => $jira_tickets{species_ticket_name},
        base64 => $jira_tickets{base64},
        comment => 'Back up done for #species# #cs_name#',
      },
      -flow_into => {
          1 => ['create_top_level_input_ids', 'label_gencode_basic', 'remove_comparison_attributes_core'],
      },
    },

    {
      -logic_name => 'create_top_level_input_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -rc_name    => 'default',
      -parameters => {
        iid_type => 'slice',
        coord_system_name => 'chromosome',
        slice => 1,
        include_non_reference => 0,
        top_level => 1,
      },
      -flow_into => {
        '2->A' => {'compare_gencode_refseq' => {'iid' => '#iid#', species => '#species#', target_db => '#target_db#', refseq_db => '#refseq_db#', cs_name => '#cs_name#', refseq_logicname => '#refseq_logicname#'}},
        'A->1' => ['comment_compare_jira'],
      },
    },

    {
      -logic_name => 'compare_gencode_refseq',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => 'default',
      -parameters => {
          cmd => 'perl '.$self->o('compare_gencode_refseq_script').
            ' -ensemblhost #expr(#target_db#->{-host})expr#'.
            ' -ensembldbname #expr(#target_db#->{-dbname})expr#'.
            ' -ensemblport #expr(#target_db#->{-port})expr#'.
            ' -ensembluser #expr(#target_db#->{-user})expr#'.
            ' -ensemblpass #expr(#target_db#->{-pass})expr#'.
            ' -dnahost #expr(#target_db#->{-host})expr#'.
            ' -dnauser '.$self->o('user_r').
            ' -dnadbname #expr(#target_db#->{-dbname})expr#'.
            ' -dnaport #expr(#target_db#->{-port})expr#'.
            ' -prodhost '.$self->o('production_db', '-host').
            ' -prodport '.$self->o('production_db', '-port').
            ' -produser '.$self->o('production_db', '-user').
            ' -proddbname '.$self->o('production_db', '-dbname').
            ' -refseqhost #expr(#refseq_db#->{-host})expr#'.
            ' -refsequser #expr(#refseq_db#->{-user})expr#'.
            ' -refseqport #expr(#refseq_db#->{-port})expr#'.
            ' -refseqpass #expr(#refseq_db#->{-pass})expr#'.
            ' -refseqdbname #expr(#refseq_db#->{-dbname})expr#'.
            ' -coord_system_version #cs_name#'.
            ' -refseq_logicname #refseq_logicname#'.
            ' -chr_num #iid# -store',
      },
    },

    {
      -logic_name => 'comment_compare_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'comment',
        type => 'Sub-task',
        ticket_name => $jira_tickets{species_ticket_name},
        base64 => $jira_tickets{base64},
        comment => 'Gencode compare done for #species# #cs_name#',
      },
    },

    {
      -logic_name => 'label_gencode_basic',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'perl '.$self->o('label_gencode_basic_script').
            ' -h #expr(#target_db#->{-host})expr#'.
            ' -D #expr(#target_db#->{-dbname})expr#'.
            ' -P #expr(#target_db#->{-port})expr#'.
            ' -u #expr(#target_db#->{-user})expr#'.
            ' -p #expr(#target_db#->{-pass})expr#'.
            ' -dnadbhost #expr(#target_db#->{-host})expr#'.
            ' -dnadbuser '.$self->o('user_r').
            ' -dnadbname #expr(#target_db#->{-dbname})expr#'.
            ' -dnadbport #expr(#target_db#->{-port})expr#'.
            ' -prodhost '.$self->o('production_db', '-host').
            ' -prodport '.$self->o('production_db', '-port').
            ' -produser '.$self->o('production_db', '-user').
            ' -proddbname '.$self->o('production_db', '-dbname').
            ' -path #cs_name#'.
            ' -write',
      },
      -flow_into => {
        1 => ['comment_basic_jira'],
      }
    },

    {
      -logic_name => 'comment_basic_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'comment',
        type => 'Sub-task',
        ticket_name => $jira_tickets{species_ticket_name},
        base64 => $jira_tickets{base64},
        comment => 'Gencode basic done for #species# #cs_name#',
      },
    },

    {
      -logic_name => 'remove_comparison_attributes_core',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -rc_name => 'default',
      -parameters => {
          db_conn => '#target_db#',
          sql => $self->o('query_delete_attributes'),
      },
      -flow_into => {
        1 => ['remove_comparison_attributes_otherfeatures'],
      }
    },

    {
      -logic_name => 'remove_comparison_attributes_otherfeatures',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -rc_name => 'default',
      -parameters => {
          db_conn => '#refseq_db#',
          sql => $self->o('query_delete_attributes'),
      },
      -flow_into => {
        1 => ['comparison_top_levels_ids'],
      }
    },
 
    {
      -logic_name => 'comparison_top_levels_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -rc_name    => 'default',
      -parameters => {
        iid_type => 'slice',
        coord_system_name => 'chromosome',
        slice => 1,
        include_non_reference => 1,
        top_level => 1,
      },
      -flow_into => {
        '2->A' => {'comparison_attributes_core_otherfeature' => {'iid' => '#iid#', species => '#species#', target_db => '#target_db#', refseq_db => '#refseq_db#', cs_name => '#cs_name#', refseq_logicname => '#refseq_logicname#'}},
        'A->1' => ['comment_overlap_jira'],
      },
    },

    {
      -logic_name => 'comparison_attributes_core_otherfeature',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'perl '.$self->o('overlap_gencode_refseq_script').
            ' -primary_host #expr(#target_db#->{-host})expr#'.
            ' -primary_dbname #expr(#target_db#->{-dbname})expr#'.
            ' -primary_port #expr(#target_db#->{-port})expr#'.
            ' -primary_pass #expr(#target_db#->{-pass})expr#'.
            ' -primary_user #expr(#target_db#->{-user})expr#'.
            ' -secondary_host #expr(#refseq_db#->{-host})expr#'.
            ' -secondary_dbname #expr(#refseq_db#->{-dbname})expr#'.
            ' -secondary_port #expr(#refseq_db#->{-port})expr#'.
            ' -secondary_user '.$self->o('user_r').
            ' -secondary_logic_name refseq_import'.
            ' -dna_host #expr(#target_db#->{-host})expr#'.
            ' -dna_user '.$self->o('user_r').
            ' -dna_port #expr(#target_db#->{-port})expr#'.
            ' -dna_dbname #expr(#target_db#->{-dbname})expr#'.
            ' -primary_set_name #primary_set_name#'.
            ' -secondary_set_name #secondary_set_name#'.
            ' -iid #iid#'.
            ' -write',
          primary_set_name => 'Ensembl',
          secondary_set_name => 'RefSeq',
      },
      -flow_into => {
        1 => ['comparison_attributes_otherfeature_core'],
      }
    },

    {
      -logic_name => 'comparison_attributes_otherfeature_core',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'perl '.$self->o('overlap_gencode_refseq_script').
            ' -primary_host #expr(#refseq_db#->{-host})expr#'.
            ' -primary_dbname #expr(#refseq_db#->{-dbname})expr#'.
            ' -primary_port #expr(#refseq_db#->{-port})expr#'.
            ' -primary_user #expr(#refseq_db#->{-user})expr#'.
            ' -primary_pass #expr(#refseq_db#->{-pass})expr#'.
            ' -primary_logic_name refseq_import'.
            ' -secondary_host #expr(#target_db#->{-host})expr#'.
            ' -secondary_dbname #expr(#target_db#->{-dbname})expr#'.
            ' -secondary_port #expr(#target_db#->{-port})expr#'.
            ' -secondary_user '.$self->o('user_r').
            ' -dna_host #expr(#target_db#->{-host})expr#'.
            ' -dna_user '.$self->o('user_r').
            ' -dna_dbname #expr(#target_db#->{-dbname})expr#'.
            ' -dna_port #expr(#target_db#->{-port})expr#'.
            ' -primary_set_name #primary_set_name#'.
            ' -secondary_set_name #secondary_set_name#'.
            ' -iid #iid#'.
            ' -write',
          primary_set_name => 'RefSeq',
          secondary_set_name => 'Ensembl',
      },
    },

    {
      -logic_name => 'comment_overlap_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'comment',
        type => 'Sub-task',
        ticket_name => $jira_tickets{species_ticket_name},
        base64 => $jira_tickets{base64},
        comment => 'Ensembl Refseq overlap done for #species# #cs_name#',
      },
    },

    {
      -logic_name => 'create_datafile_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'create',
        type => 'Sub-task',
        parent => $jira_tickets{release_ticket_name},
        ticket_name => $jira_tickets{datafile_ticket_name},
        base64 => $jira_tickets{base64},
        assignee => 'self',
      },
      -flow_into => {
          1 => ['create_datafile_input'],
      },
    },

    {
      -logic_name => 'create_datafile_input',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputlist => [[$self->o('staging1_db', '-host'), $self->o('staging1_db', '-user'), $self->o('staging1_db', '-port')]],
          column_names => ['host', 'user', 'port'],
      },
      -flow_into => {
          '2' => ['check_datafiles'],
      },
    },

# The script is really verbose
    {
      -logic_name => 'check_datafiles',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'perl '.$self->o('check_datafiles_script').' -release '.$self->o('ensembl_release').' -host #host# -user #user# -port #port# -datafile_dir '.$self->o('webdev_nfs').' -group rnaseq --force_samtools_binary --samtools_binary '.$self->o('samtools').' 1>/dev/null',
      },
      -flow_into => {
          1 => ['comment_check_datafile_jira'],
      },
    },

    {
      -logic_name => 'comment_check_datafile_jira',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket',
      -rc_name => 'default',
      -parameters => {
        cmd => 'comment',
        type => 'Sub-task',
        ticket_name => $jira_tickets{datafile_ticket_name},
        base64 => $jira_tickets{base64},
        comment => 'Data files are all OK',
      },
    },
  );

  foreach my $analysis (@analysis) {
    $analysis->{-max_retry_count} = 1 unless (exists $analysis->{-max_retry_count});
  }
  return \@analysis;
}


=head2 resource_classes

 Arg [1]    : None
 Description: Resources needed for the pipeline, it uses the default one at 1GB and one requesting 4GB if needed
 Returntype : Hashref
 Exceptions : None

=cut

sub resource_classes {
    my $self = shift;

    return {
        %{ $self->SUPER::resource_classes() },  # inherit other stuff from the base class
      'default' => { LSF => '-M1000 -R"select[mem>1000] rusage[mem=1000]"'},
      '4GB' => { LSF => '-M4000 -R"select[mem>4000] rusage[mem=4000]"'},
    };
}

1;
