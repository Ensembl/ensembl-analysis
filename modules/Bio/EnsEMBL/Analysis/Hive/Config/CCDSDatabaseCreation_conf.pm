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

=cut

package CCDSDatabaseCreation_conf;

use strict;
use warnings;

use File::Spec::Functions;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use POSIX 'strftime';

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },

# Current date (YYYY_MM and YYYY_MM_DD). To be used for directories, filenames and metadata.
'YYYY_MM' => strftime('%Y',localtime).'_'.strftime('%m',localtime),
'YYYY_MM_DD' => strftime('%Y',localtime).'-'.strftime('%m',localtime).'-'.strftime('%d',localtime),

##################################
# Variable settings. Change these.
##################################

# /path/to/ensembl/genebuild/env/genebuild.sh which defines the variable BIOPERL_LIB
'genebuild_environment_file' => '',

# main working directory which will contain the code and data (must not exist)
'ccds_update_dir'  => '/path/to/scratch/ccds_database_creation_'.$self->o('YYYY_MM'),

# warehouse directory for backup (must not exist)
'warehouse_dir'    => '/path/to/warehouse/ccds_database_creation_'.$self->o('YYYY_MM'),

# CCDS FTP parameters for the CCDS data downloading step
'ccds_ftp_user'    => '',
'ccds_ftp_password'=> '',
'ccds_ftp_url'     => '',
'ccds_ftp_filename'=> 'CCDS.*.tar.gz',

# Directory that contains the CCDS scripts required by the analyses of this pipeline.
'ccds_scripts_dir'   => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'ccds'),

# previous/latest Ensembl release number
'release' => 0, # Example: 90

# Ensembl latest existing human core database
'human_core_db_host' => '',
'human_core_db_port' => '',

# Ensembl latest existing mouse core database
'mouse_core_db_host' => '',
'mouse_core_db_port' => '',

# Ensembl production database details
'production_db_name' => 'ensembl_production_'.$self->o('release_number'),
'production_db_host' => '',
'production_db_port' => '',

# Databases to be created
'pipeline_db_host'=>'',
'pipeline_db_port'=>'',
'ccds_track_db_host'=>'',
'ccds_track_db_port'=>'',
'ccds_human_db_host'=>'',
'ccds_human_db_port'=>'',
'ccds_mouse_db_host'=>'',
'ccds_mouse_db_port'=>'',

# read-only and write users and passwords
'user_r'            => '',
'pass_r'            => '',
'user_w'            => '',
'pass_w'            => '',

########################################################
# Static settings (do not need to change the following).
########################################################

'human_assembly_version' => 38, # Example: 38
'mouse_assembly_version' => 38, # Example: 38

'human_assembly_path' => 'GRCh38',
'mouse_assembly_path' => 'GRCm38',

'human_taxon_id' => '9606',
'mouse_taxon_id' => '10090',

# maximum number of tries for each table in the load_data step
'load_data_max_tries' => 4,

'human_core_db_name' => 'homo_sapiens_core_'.$self->o('release').'_'.$self->o('human_assembly_version'),
'mouse_core_db_name' => 'mus_musculus_core_'.$self->o('release').'_'.$self->o('mouse_assembly_version'),

'data_dir'         => $self->o('ccds_update_dir').'/data',
'human_dir'        => $self->o('ccds_update_dir').'/human',
'mouse_dir'        => $self->o('ccds_update_dir').'/mouse',
'backups_dir'      => $self->o('ccds_update_dir').'/backups',
'code_dir'         => $self->o('ccds_update_dir').'/code',
'docfile'          => $self->o('ccds_update_dir').'/CCDSDatabaseCreation_'.$self->o('YYYY_MM').'.txt',

'driver'            => $self->o('hive_driver'),

 # Databases

'pipeline_db' => {
  -dbname => $self->o('dbowner').'_ccdsdatabasecreation_pipeline_'.$self->o('YYYY_MM'),
  -host   => $self->o('pipeline_db_host'),
  -port   => $self->o('pipeline_db_port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('pass_w'),
  -driver => $self->o('driver'),
},

'ccds_track_db' => {
  -dbname => $self->o('dbowner').'_cdstrack_'.$self->o('YYYY_MM'),
  -host   => $self->o('ccds_track_db_host'),
  -port   => $self->o('ccds_track_db_port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('pass_w'),
},

'ccds_human_db' => {
  -dbname => $self->o('dbowner').'_human_cdsonly_'.$self->o('YYYY_MM'),
  -host   => $self->o('ccds_human_db_host'),
  -port   => $self->o('ccds_human_db_port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('pass_w'),
},

'ccds_mouse_db' => {
  -dbname => $self->o('dbowner').'_mouse_cdsonly_'.$self->o('YYYY_MM'),
  -host   => $self->o('ccds_mouse_db_host'),
  -port   => $self->o('ccds_mouse_db_port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('pass_w'),
},

'human_core_db' => {
  -dbname => $self->o('human_core_db_name'),
  -host   => $self->o('human_core_db_host'),
  -port   => $self->o('human_core_db_port'),
  -user   => $self->o('user_r'),
  -pass   => $self->o('pass_r'),
},

'mouse_core_db' => {
  -dbname => $self->o('mouse_core_db_name'),
  -host   => $self->o('mouse_core_db_host'),
  -port   => $self->o('mouse_core_db_port'),
  -user   => $self->o('user_r'),
  -pass   => $self->o('pass_r'),
},

'production_db' => {
  -dbname => $self->o('production_db_name'),
  -host   => $self->o('production_db_host'),
  -port   => $self->o('production_db_port'),
  -user   => $self->o('user_r'),
  -pass   => $self->o('pass_r'),
},

    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

###############################################################################
#
# SET UP CODE BASE AND ENVIRONMENT VARIABLES
#
# Taking out code and setting up data directories.
#
###############################################################################

      {
        # Make the directories that will contain the code, data and backups.
        -logic_name => 'make_dirs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'mkdir -p '.$self->o('ccds_update_dir').' ;'.
                                    'mkdir -p '.$self->o('data_dir').' ;'.
                                    'mkdir -p '.$self->o('human_dir').' ;'.
                                    'mkdir -p '.$self->o('mouse_dir').' ;'.
                                    'mkdir -p '.$self->o('backups_dir').' ;'.
                                    'mkdir -p '.$self->o('code_dir').' ;'.
                                    'mkdir -p '.$self->o('warehouse_dir')     
                       },
        -flow_into => { 1 => ['download_code'] },
        -meadow_type => 'LOCAL',
        -input_ids => [{}],
      },

      {
        # Download the Ensembl code from GitHub.
        -logic_name => 'download_code',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'git clone --depth=1 https://github.com/Ensembl/ensembl-analysis.git '.$self->o('code_dir').'/ensembl-analysis & '.
                                    'git clone --depth=1 https://github.com/Ensembl/ensembl.git '.$self->o('code_dir').'/ensembl & '.
                                    'git clone --depth=1 https://github.com/Ensembl/ensembl-external.git '.$self->o('code_dir').'/ensembl-external'
                       },
        -flow_into => { 1 => ['record_code_versions'] },
        -meadow_type => 'LOCAL',
      },
      
      {
        # Record the code versions.
        -logic_name => 'record_code_versions',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'for ec in '.$self->o('code_dir').'/ens* ; do '.
                                    '  cd $ec ; '.
                                    '  echo $ec ; '.
                                    '  git log -n 1 --pretty=format:"%H" ; printf "\n"'.
                                    '  cd .. ; '.
                                    'done '.
                                    '> '.$self->o('docfile')
                       },
        -flow_into => { 1 => ['set_genebuild_environment_variables'] },
        -meadow_type => 'LOCAL',
      },

      {
        # Set the environment variables. Note BIOPERL_LIB is set in the genebuild source file.
        -logic_name => 'set_genebuild_environment_variables',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'source '.$self->o('genebuild_environment_file').' ; '.
                                    'export PERL5LIB=$BIOPERL_LIB:'.
                                                     $self->o('code_dir').'/ensembl/modules/:'.
                                                     $self->o('code_dir').'/ensembl-analysis/modules:'.
                                                     $self->o('code_dir').'/ensembl-external/modules'
                       },
        -flow_into => { 1 => ['download_ccds_data'] },
        -meadow_type => 'LOCAL',
      },

###############################################################################
#
# DATA DOWNLOAD AND DATABASE SETUP
#
###############################################################################

      {
        # Get CCDS data via FTP, copy the data file to the warehouse and unzip it in the data dir.
        -logic_name => 'download_ccds_data',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'wget --user='.$self->o('ccds_ftp_user').
                                    '     --password='.$self->o('ccds_ftp_password').
                                    '    '.$self->o('ccds_ftp_url').'/'.$self->o('ccds_ftp_filename').
                                    '     --directory-prefix='.$self->o('data_dir').' && '.
                                    'cp '.$self->o('data_dir').'/'.$self->o('ccds_ftp_filename').' '.$self->o('warehouse_dir').' && '.
                                    'tar -xvzf '.$self->o('data_dir').'/'.$self->o('ccds_ftp_filename').' --directory='.$self->o('data_dir').' && '.
                                    'rm '.$self->o('data_dir').'/'.$self->o('ccds_ftp_filename')
                       },
        -flow_into => { 1 => ['prepare_cdstrack'] },
        -meadow_type => 'LOCAL',
      },

      {
        # Get CCDS data via FTP, copy the data file to the warehouse and unzip it in the data dir.
        -logic_name => 'prepare_cdstrack',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'perl '.$self->o('ccds_scripts_dir').'/prepare_cdstrack.pl -dir '.$self->o('data_dir')
                       },
        -flow_into => { 1 => ['create_ccds_track_db'] },
        -meadow_type => 'LOCAL',
      },

      {
        # Create the CCDS track database.
        -logic_name => 'create_ccds_track_db',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'mysql -NB -u'.$self->o('user_w').
                                             ' -p'.$self->o('pass_w').
                                             ' -h'.$self->o('ccds_track_db','-host').
                                             ' -P'.$self->o('ccds_track_db','-port').
                                             ' -e "CREATE DATABASE '.$self->o('ccds_track_db','-dbname').'"'
                       },
        -flow_into => { 1 => ['load_new_create_tables'] },
        -meadow_type => 'LOCAL',
      },

      {
        # Create the tables in the CCDS track database.
        -logic_name => 'load_new_create_tables',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'mysql -NB -u'.$self->o('user_w').
                                             ' -p'.$self->o('pass_w').
                                             ' -h'.$self->o('ccds_track_db','-host').
                                             ' -P'.$self->o('ccds_track_db','-port').
                                             ' -D'.$self->o('ccds_track_db','-dbname').
                                             ' < '.$self->o('data_dir').'/sql/new_createTables.sql'
                       },
        -flow_into => { 1 => ['load_data'] },
        -meadow_type => 'LOCAL',
      },

      {
        # Load the data into the CCDS track database.
        -logic_name => 'load_data',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'mv '.$self->o('data_dir').'/data/new_Interpretations.txt '.$self->o('data_dir').'/data/Interpretations.txt ;'.
                                    'for file in `ls '.$self->o('data_dir').'/data/*.txt`; do '.
                                    '  mysqlimport -u'.$self->o('user_w').
                                                 ' -p'.$self->o('pass_w').
                                                 ' -h'.$self->o('ccds_track_db','-host').
                                                 ' -P'.$self->o('ccds_track_db','-port').
                                                 ' --local '.
                                                 $self->o('ccds_track_db','-dbname').
                                                 ' $file ;'.
                                    'done'
                       },
        -flow_into => { 1 => ['load_new_create_keys'] },
        -meadow_type => 'LOCAL',
      },

      {
        # Create the keys for the CCDS track database.
        -logic_name => 'load_new_create_keys',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'mysql -NB -u'.$self->o('user_w').
                                             ' -p'.$self->o('pass_w').
                                             ' -h'.$self->o('ccds_track_db','-host').
                                             ' -P'.$self->o('ccds_track_db','-port').
                                             ' -D'.$self->o('ccds_track_db','-dbname').
                                             ' < '.$self->o('data_dir').'/sql/new_createKeys.sql'
                       },
        -flow_into => { 1 => ['remove_public_notes'] },
        -meadow_type => 'LOCAL',
      },

      {
        # Remove items with multiple public notes as it causes things to fall over when using the API.
        -logic_name => 'remove_public_notes',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'for id in `mysql -NB -u'.$self->o('user_w').
                                                        ' -p'.$self->o('pass_w').
                                                        ' -h'.$self->o('ccds_track_db','-host').
                                                        ' -P'.$self->o('ccds_track_db','-port').
                                                        ' -D'.$self->o('ccds_track_db','-dbname').
                                                        ' -e"SELECT i.ccds_uid '.
                                                        '    FROM Interpretations i '.
                                                        '    WHERE i.interpretation_subtype_uid = 17 '.
                                                        '    GROUP BY i.ccds_uid HAVING COUNT(*) > 1`; do '.
                                       'mysql -NB -u'.$self->o('user_w').
                                                ' -p'.$self->o('pass_w').
                                                ' -h'.$self->o('ccds_track_db','-host').
                                                ' -P'.$self->o('ccds_track_db','-port').
                                                ' -D'.$self->o('ccds_track_db','-dbname').
                                                ' -e"'.
                                       'DELETE FROM Interpretations '.
                                       'WHERE interpretation_subtype_uid = 17 AND '.
                                       '      integer_val = 1 AND '.
                                       '      ccds_uid = $id "; '.
                                     'done'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['dump_external_db_table'] },
      },

      {
        # Dump the external_db table from the production database.
        -logic_name => 'dump_external_db_table',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'   => 'mysql -NB -u'.$self->o('user_r').
                                             ' -h'.$self->o('production_db','-host').
                                             ' -P'.$self->o('production_db','-port').
                                             ' -D'.$self->o('production_db','-dbname').
                                             ' -e"SELECT * FROM external_db"'.
                                             ' > '.$self->o('ccds_update_dir').'/external_db.txt'
                       },
        -flow_into => { 1 => ['create_human_ccds_db'] },
        -meadow_type => 'LOCAL',
      },
      
      {
        # Create the human CCDS database from the previous core database.
              -logic_name => 'create_human_ccds_db',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
              -parameters => {
                                create_type => 'copy',
                                source_db => $self->o('human_core_db'),
                                target_db => $self->o('ccds_human_db'),
                                db_dump_file => $self->o('backups_dir').'/human_core_db.tmp',
                             },
             -flow_into => { 1 => ['create_mouse_ccds_db'] },
      },

      {
        # Create the mouse CCDS database from the previous core database.
              -logic_name => 'create_mouse_ccds_db',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
              -parameters => {
                                create_type => 'copy',
                                source_db => $self->o('mouse_core_db'),
                                target_db => $self->o('ccds_mouse_db'),
                                db_dump_file => $self->o('backups_dir').'/mouse_core_db.tmp',
                             },
             -flow_into => { 1 => ['truncate_human_external_db'] },
      },

      {
        # Truncate the human external db table in preparation for loading the table.
        -logic_name => 'truncate_human_external_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_human_db'),
                         input_query => 'TRUNCATE external_db;',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['truncate_mouse_external_db'] },
      },
      
      {
        # Truncate the mouse external db table in preparation for loading the table.
        -logic_name => 'truncate_mouse_external_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_mouse_db'),
                         input_query => 'TRUNCATE external_db;',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['load_external_db_into_human_ccds_db'] },
      },

      {
        # Load the external db table into the human CCDS db.
        -logic_name => 'load_external_db_into_human_ccds_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                          db_conn => $self->o('ccds_human_db'),
                          executable => 'mysqlimport',
                          prepend => ["--local"],
                          append => [$self->o('ccds_update_dir').'/external_db.txt'],
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['load_external_db_into_mouse_ccds_db'] },
      },

      {
        # Load the external db table into the mouse CCDS db.
        -logic_name => 'load_external_db_into_mouse_ccds_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                          db_conn => $self->o('ccds_mouse_db'),
                          executable => 'mysqlimport',
                          prepend => ["--local"],
                          append => [$self->o('ccds_update_dir').'/external_db.txt'],
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_human_external_db_null'] },
      },

      {
        # Set the human external db table db versions to NULL rather than 'NULL'
        -logic_name => 'set_human_external_db_null',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_human_db'),
                         input_query => 'UPDATE external_db SET db_release=\N WHERE db_release like "NULL";'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_mouse_external_db_null'] },
      },
      
      {
        # Set the mouse external db table db versions to NULL rather than 'NULL'
        -logic_name => 'set_mouse_external_db_null',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_mouse_db'),
                         input_query => 'UPDATE external_db SET db_release=\N WHERE db_release like "NULL";'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['load_ccds_transcripts_into_human_ccds_db'] },
      },

      {
        # Load the CCDS transcripts into the human CCDS database.
        -logic_name => 'load_ccds_transcripts_into_human_ccds_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                          'cmd'   => 'ncbi_build_number=`mysql -NB -u'.$self->o('user_r').
                                                                 ' -h'.$self->o('ccds_track_db','-host').
                                                                 ' -P'.$self->o('ccds_track_db','-port').
                                                                 ' -D'.$self->o('ccds_track_db','-dbname').
                                                                 ' -e"SELECT max(gv.ncbi_build_number) '.
                                                                 '    FROM GroupVersions gv,Groups g '.
                                                                 '    WHERE gv.group_uid = g.group_uid AND '.
                                                                 '          gv.version = g.current_version AND '.
                                                                 '          g.tax_id = '.$self->o('human_taxon_id').';"`; '.
                                     'perl '.$self->o('ccds_scripts_dir').'/load_ccds_transcripts_using_api.pl '.
                                     ' -cds_host '.$self->o('ccds_track_db','-host').
                                     ' -cds_port '.$self->o('ccds_track_db','-port').
                                     ' -cds_user '.$self->o('user_r').
                                     ' -cds_dbname '.$self->o('ccds_track_db','-dbname').
                                     ' -host '.$self->o('ccds_human_db','-host').
                                     ' -port '.$self->o('ccds_human_db','-port').
                                     ' -user '.$self->o('ccds_human_db','-user').
                                     ' -pass '.$self->o('ccds_human_db','-pass').
                                     ' -dbname '.$self->o('ccds_human_db','-dbname').
                                     ' -dna_host '.$self->o('human_core_db','-host').
                                     ' -dna_user '.$self->o('human_core_db','-user').
                                     ' -dna_port '.$self->o('human_core_db','-port').
                                     ' -dna_dbname '.$self->o('human_core_db','-dbname').
                                     ' -tax_id '.$self->o('human_taxon_id').
                                     ' -ncbi_build_number $ncbi_build_number'.
                                     ' -analtype ccds_gene '.
                                     ' -path '.$self->o('human_assembly_path')
                       },
        -rc_name => 'default_himem',
        -flow_into => { 1 => ['load_ccds_transcripts_into_mouse_ccds_db'] },
      },

      {
        # Load the CCDS transcripts into the mouse CCDS database.
        -logic_name => 'load_ccds_transcripts_into_mouse_ccds_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                          'cmd'   => 'ncbi_build_number=`mysql -NB -u'.$self->o('user_r').
                                                                 ' -h'.$self->o('ccds_track_db','-host').
                                                                 ' -P'.$self->o('ccds_track_db','-port').
                                                                 ' -D'.$self->o('ccds_track_db','-dbname').
                                                                 ' -e"SELECT max(gv.ncbi_build_number) '.
                                                                 '    FROM GroupVersions gv,Groups g '.
                                                                 '    WHERE gv.group_uid = g.group_uid AND '.
                                                                 '          gv.version = g.current_version AND '.
                                                                 '          g.tax_id = '.$self->o('mouse_taxon_id').';"`; '.
                                     'perl '.$self->o('ccds_scripts_dir').'/load_ccds_transcripts_using_api.pl '.
                                     ' -cds_host '.$self->o('ccds_track_db','-host').
                                     ' -cds_port '.$self->o('ccds_track_db','-port').
                                     ' -cds_user '.$self->o('user_r').
                                     ' -cds_dbname '.$self->o('ccds_track_db','-dbname').
                                     ' -host '.$self->o('ccds_mouse_db','-host').
                                     ' -port '.$self->o('ccds_mouse_db','-port').
                                     ' -user '.$self->o('ccds_mouse_db','-user').
                                     ' -pass '.$self->o('ccds_mouse_db','-pass').
                                     ' -dbname '.$self->o('ccds_mouse_db','-dbname').
                                     ' -dna_host '.$self->o('mouse_core_db','-host').
                                     ' -dna_user '.$self->o('mouse_core_db','-user').
                                     ' -dna_port '.$self->o('mouse_core_db','-port').
                                     ' -dna_dbname '.$self->o('mouse_core_db','-dbname').
                                     ' -tax_id '.$self->o('mouse_taxon_id').
                                     ' -ncbi_build_number $ncbi_build_number'.
                                     ' -analtype ccds_gene '.
                                     ' -path '.$self->o('mouse_assembly_path')
                       },
        -rc_name => 'default_himem',
        -flow_into => { 1 => ['set_human_canonical_transcripts'] },
      },
      
      {
        # Set human canonical transcripts (CCDS genes have only 1 transcript per gene).
        -logic_name => 'set_human_canonical_transcripts',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_human_db'),
                         input_query => 'UPDATE gene g,transcript t,analysis a '.
                                        'SET g.canonical_transcript_id=t.transcript_id '.
                                        'WHERE g.gene_id=t.gene_id AND '.
                                        '      g.analysis_id=a.analysis_id AND '.
                                        '      logic_name="ccds_gene";',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_mouse_canonical_transcripts'] },
      },
      
      {
        # Set mouse canonical transcripts (CCDS genes have only 1 transcript per gene).
        -logic_name => 'set_mouse_canonical_transcripts',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_mouse_db'),
                         input_query => 'UPDATE gene g,transcript t,analysis a '.
                                        'SET g.canonical_transcript_id=t.transcript_id '.
                                        'WHERE g.gene_id=t.gene_id AND '.
                                        '      g.analysis_id=a.analysis_id AND '.
                                        '      logic_name="ccds_gene";',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['update_human_gene_sources'] },
      },
      
      {
        # Update human ccds db gene sources to "ccds".
        -logic_name => 'update_human_gene_sources',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_human_db'),
                         input_query => 'UPDATE gene SET source="ccds";',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['update_mouse_gene_sources'] },
      },
      
      {
        # Update mouse ccds db gene sources to "ccds".
        -logic_name => 'update_mouse_gene_sources',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_mouse_db'),
                         input_query => 'UPDATE gene SET source="ccds";',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['update_human_metadata'] },
      },
           
      {
        # Update human ccds db metadata including release dates and genebuild id.
        -logic_name => 'update_human_metadata',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_human_db'),
                         input_query => 'UPDATE meta SET meta_value='.$self->o('genebuilder_id').
                                        '            WHERE meta_key="genebuild.id";'.
                                        'UPDATE meta SET meta_value="'.$self->o('YYYY_MM_DD').'" '.
                                        '            WHERE meta_key="genebuild.initial_release_date";'.
                                        'UPDATE meta SET meta_value="'.$self->o('YYYY_MM_DD').'" '.
                                        'WHERE meta_key="genebuild.last_geneset_update";'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['update_mouse_metadata'] },
      },
                  
      {
        # Update mouse ccds db metadata including release dates and genebuild id.
        -logic_name => 'update_mouse_metadata',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_mouse_db'),
                         input_query => 'UPDATE meta SET meta_value='.$self->o('genebuilder_id').
                                        '            WHERE meta_key="genebuild.id";'.
                                        'UPDATE meta SET meta_value="'.$self->o('YYYY_MM_DD').'" '.
                                        '            WHERE meta_key="genebuild.initial_release_date";'.
                                        'UPDATE meta SET meta_value="'.$self->o('YYYY_MM_DD').'" '.
                                        'WHERE meta_key="genebuild.last_geneset_update";'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['delete_human_xrefs'] },
      },

      {
        # Delete Ens_Hs_transcript and Ens_Hs_translation xrefs for human
        -logic_name => 'delete_human_xrefs',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_human_db'),
                         input_query => 'DELETE FROM xref '.
                                        'WHERE external_db_id IN '.
                                        '(SELECT external_db_id '.
                                        ' FROM external_db '.
                                        ' WHERE db_name '.
                                        ' IN("Ens_Hs_transcript","Ens_Hs_translation"));'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['delete_mouse_xrefs'] },
      },
     
      {
        # Delete Ens_Mm_transcript and Ens_Mm_translation xrefs for mouse
        -logic_name => 'delete_mouse_xrefs',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     =>  $self->o('ccds_mouse_db'),
                         input_query => 'DELETE FROM xref '.
                                        'WHERE external_db_id IN '.
                                        '(SELECT external_db_id '.
                                        ' FROM external_db '.
                                        ' WHERE db_name '.
                                        ' IN("Ens_Mm_transcript","Ens_Mm_translation"));'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['backup_ccds_track_db'] },
      },

      {
        -logic_name => 'backup_ccds_track_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     => $self->o('ccds_track_db'),
                         executable    => 'mysqldump',
                         output_file   => $self->o('warehouse_dir').'/'.$self->o('ccds_track_db','-dbname').'.sql',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['backup_ccds_human_db'] },
      },
      
      {
        -logic_name => 'backup_ccds_human_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     => $self->o('ccds_human_db'),
                         executable    => 'mysqldump',
                         output_file   => $self->o('warehouse_dir').'/'.$self->o('ccds_human_db','-dbname').'.sql',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['backup_ccds_mouse_db'] },
      },
      
      {
        -logic_name => 'backup_ccds_mouse_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn     => $self->o('ccds_mouse_db'),
                         executable    => 'mysqldump',
                         output_file   => $self->o('warehouse_dir').'/'.$self->o('ccds_mouse_db','-dbname').'.sql',
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['email'] },
      },

      {
        -logic_name => 'email',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
        -parameters => {
                         email => $self->o('email_address'),
                         subject => 'AUTOMATED EMAIL: CCDS databases creation finished',
                         text => 'The following CCDS databases have been created: '.
                                 $self->o('ccds_human_db','-dbname').'@'.
                                 $self->o('ccds_human_db','-host').':'.
                                 $self->o('ccds_human_db','-port').' and '.
                                 $self->o('ccds_mouse_db','-dbname').'@'.
                                 $self->o('ccds_mouse_db','-host').':'.
                                 $self->o('ccds_mouse_db','-port').' . '.
                                 'The corresponding backup files containing the database dumps have been created at: '.
                                 $self->o('warehouse_dir').
                                 ' . Code versions recorded as below:',
                         file => $self->o('docfile'),
                         command => q{cat},
                       },
        -rc_name => 'default',
      },
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    };
}

# override the default method, to force an automatic loading of the registry in all workers
#sub beekeeper_extra_cmdline_options {
#    my $self = shift;
#    return "-reg_conf ".$self->o("registry");
#}

sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => '-M900 -R"select[mem>900]"' },
    'default_himem' => { LSF => '-M5000 -R"select[mem>5000]"' },
  }
}

1;
