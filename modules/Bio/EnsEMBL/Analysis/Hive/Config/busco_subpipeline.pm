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

package Bio::EnsEMBL::Analysis::Hive::Config::busco_subpipeline;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
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
    # need for pipeline set up
    
    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'             => '' || $self->o('production_name').'_busco_'.$self->o('ensembl_release'),
    'user_r'                    => 'ensro', # read only db user
    # 'user'                      => '', # write db user. You may not need that for this pipeline
    # 'password'                  => '', # password for write db user. You may not need that for this pipeline
    'pipe_db_server'            => '', # host for pipe db
    'pipe_db_port'              => '', # port for pipeline host
    'release_number'            => '' || $self->o('ensembl_release'),
    'ensembl_release'           => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    'species_name'              => '', # e.g. mus_musculus
    'production_name'           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names


    # need for busco part:
    ## Name of the output subdirectory.
    'busco_output_name'         => '', # You need to specify this and be careful. Just put the name where your output file is going to be.   
    'load_inputfile_only'       => 1, 
    'run_busco_singularity'     => 1, 
    ## What fasta file to check (if you provide it). File should contain protein segments with appropriate headers. 
    'busco_input_file'          => '/hps/nobackup2/production/ensembl/kbillis/research/busco_QC/busco/fergal_test_butterfly/input/non_verts/butterfly/lep_prots_unique_ids.fa', 
    ## where is your output/working dir?
    'output_path'               => '', 
    'input_dir'                 => catdir($self->o('output_path'), 'input_dir'),  
    'working_dir'               => catdir($self->o('output_path'), 'working_dir'),
    'busco_input_file_std'      => catdir($self->o('input_dir'), 'input_file.fa' ), 
    'pipe_db_name'              => $self->o('dbowner').'_'.$self->o('production_name').'_pipe_'.$self->o('release_number'), 
    ## If you want me to check a db, provide one (run cannonical transcripts first.): 
    'busco_db_name'             => '', 
    'busco_db_server'           => '', 
    'busco_db_port'             => '', 
    'busco_params'              => ' -f ', # -f is for delete and overwrite existing dirs 
    ## in case normal busco. This can be name of dataset or location where dataset stored in the farm. 
    'busco_set'                 => 'lepidoptera_odb10',    # e.g. glires_odb10 , 
    ## in case of singularity: 
    'output_path_singlularity'               => '/hps/nobackup2/singularity/kbillis/busco_test_ehive/', 
    'busco_singlularity_default_config'      => '/hps/nobackup2/production/ensembl/kbillis/research/busco_QC/busco/docker/busco_config_docker.ini',
    'busco_singularity_image'                => '/hps/nobackup2/singularity/kbillis/busco-v5.0.0_cv1.simg',
    'busco_singularity_dataset'              => 'lepidoptera_odb10', 
    # this is the output file. May copy it afterwards to another location 
    'busco_output_file'  => $self->o('output_path_singlularity').$self->o('busco_output_name').'/short_summary.specific.'.$self->o('busco_singularity_dataset').'.txt.',

########################
    # Pipe and busco db info
########################

# The following might not be known in advance, since the come from other pipelines
# These values can be replaced in the analysis_base table if they're not known yet
# If they are not needed (i.e. no projection or rnaseq) then leave them as is

    'busco_input_file_stid'         => 'stable_id_to_dump.txt',
    # 'base_busco_path'        => '/hps/nobackup2/production/ensembl/kbillis/research/busco_QC/busco/bin/busco', 
    'base_busco_path'        => '/homes/kbillis/.local/bin/busco', 
    'default_busco_config'   => '/hps/nobackup2/production/ensembl/kbillis/research/busco_QC/busco/busco_config_bk.ini', 

######################################################
    #
    # Mostly constant settings
    #
######################################################

    ensembl_analysis_script    => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    print_protein_script_path  => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'print_translations.pl'),


########################
    # db info
########################

    'busco_db' => {
      -dbname => $self->o('busco_db_name'),
      -host   => $self->o('busco_db_server'),
      -port   => $self->o('busco_db_port'),
      -user   => $self->o('user_r'),
      # -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
      },

#######################
    # Extra db settings
########################
    num_tokens => 10,

    };
}


sub pipeline_wide_parameters {
  my ($self) = @_;


  return {
  	%{$self->SUPER::pipeline_wide_parameters},
  	load_inputfile_only => $self->o('load_inputfile_only'),
  	run_busco_singularity => $self->o('run_busco_singularity'),
  }
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;


  return [


# Is there an input file? else dump translations not dump all canonical transcript translations



       {
        -logic_name => 'create_input_dir',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mkdir -p '.$self->o('input_dir'),
                       },
        -input_ids => [{}],
        -flow_into     => WHEN ( '#load_inputfile_only# == 1'  => [ 'copy_input_file' ], 
                                 '#load_inputfile_only# == 0' => ['dump_canonical_stable_ids']),
                          # ELSE [ 'dump_canonical_stable_ids' ] ),
        -rc_name => 'default',

      },


      
      {
        -logic_name => 'dump_canonical_stable_ids',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn => $self->o('busco_db'),
                         input_query => 'SELECT transcript.stable_id from gene, transcript '.
                          ' WHERE gene.gene_id = transcript.gene_id '.
                          ' AND gene.canonical_transcript_id = transcript.transcript_id '.
                          ' AND transcript.biotype = "protein_coding" limit 2 ', 
                         command_out => q( grep 'ENS' > #busco_input_file_m#), 
                         busco_input_file_m => catfile($self->o('input_dir'),  $self->o('busco_input_file_stid') ),
                         prepend => ['-NB', '-q'],
                       },
        -rc_name => '2GB',
        -flow_into => {
                        1 => ['print_translations'],
                      },
      },

      {
        -logic_name => 'print_translations',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('print_protein_script_path').
                                ' --user='.$self->o('user_r').
                                # ' -pass='.$self->o('password').
                                ' --host='.$self->o('busco_db','-host').
                                ' --port='.$self->o('busco_db','-port').
                                ' --dbname='.$self->o('busco_db','-dbname').
                                ' --id_file='.catfile($self->o('input_dir'), $self->o('busco_input_file_stid') ).
                                ' --output_file='.$self->o('busco_input_file_std'), 
                       },
        -rc_name => 'default',
        -flow_into     => WHEN ( '#run_busco_singularity# == 1'  => [ 'run_busco_singularity' ], 
                                 '#load_inputfile_only# == 0' => ['run_busco_normal']),
      },

      {
        -logic_name => 'copy_input_file',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'cp '.$self->o('busco_input_file').' '.catfile($self->o('input_dir'),'input_file.fa'), 
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['run_busco'],
                      },
        -flow_into     => WHEN ( '#run_busco_singularity# == 1'  => [ 'run_busco_singularity' ], 
                                 '#load_inputfile_only# == 0' => ['run_busco_normal']),
                      
      },
      
      # run busco in normal mode 
      {
        -logic_name => 'run_busco_normal',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'busco '. 
                         $self->o('busco_params').
                         ' -i '.catfile($self->o('input_dir'),'input_file.fa').
                         ' -m prot '.
                         ' -l '.$self->o('busco_set').
                         ' -o '.$self->o('busco_output_name').
                         ' --config '.$self->o('default_busco_config'),
                       },
        -rc_name => 'busco_rc',
      },

      # run busco in singularity mode 
      {
        -logic_name => 'run_busco_singularity',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
        	             cmd => 'mkdir '.$self->o('output_path_singlularity').';'. 
        	                    'cd '.$self->o('output_path_singlularity').';'.
        	                    'cp '.$self->o('busco_singlularity_default_config').' ./ '.' ;'. 
        	                    'cp '.catfile($self->o('input_dir'),'input_file.fa').' ./ '.' ;'. 
        	                    'singularity exec '.$self->o('busco_singularity_image').' busco -i input_file.fa -m prot -l '.
        	                        $self->o('busco_singularity_dataset').' -o '. $self->o('busco_output_name').' --config busco_config_docker.ini'.' ;'. 
        	                    'cd - '.';',  
        	             },
        -rc_name => 'busco_singularity_rc',
        -flow_into => {
                        1 => ['check_cp_file'],
                      },

      },

      {
        -logic_name => 'check_cp_file',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => { 
          cmd => 'if [ ! -f "'.$self->o('busco_output_file').'" ] ; then die;fi',
        },        
        -rc_name => 'default',
      }


    ];
}



sub resource_classes {
  my $self = shift;

  return {
    '2GB' => { LSF => $self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default' => { LSF => $self->lsf_resource_builder('production-rh74', 900, [$self->default_options->{'pipe_db_server'}], [$self->default_options->{'num_tokens'}])},
    'busco_rc' => { LSF => $self->lsf_resource_builder('production-rh74', 23900, [$self->default_options->{'pipe_db_server'}], [$self->default_options->{'num_tokens'}])},
    'busco_singularity_rc' =>  { LSF => $self->lsf_resource_builder('production-rh74', 35000, [$self->default_options->{'pipe_db_server'}], [$self->default_options->{'num_tokens'}])},
    }
}

sub hive_capacity_classes {
  my $self = shift;

  return {
    'hc_very_low'    => 35,
    'hc_low'    => 200,
    'hc_medium' => 500,
    'hc_high'   => 1000,
    };
}

1;
