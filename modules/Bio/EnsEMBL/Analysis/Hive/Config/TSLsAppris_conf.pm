=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::Config::TSLsAppris_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::TSLsAppris_conf;

use strict;
use warnings;

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
        ensembl_release => 91, # Use it on the commandline: -ensembl_release XX
        pass_r => undef, # (optional) Use it on the command line: -password mysql_ro_password
        user_r => 'ensro',
        human_gencode_version => '27', # Use it on the command line: -human_gencode_version XX
        mouse_gencode_version => 'M16', # Use it on the command line: -password mysql_rw_password
        do_human => 0,
        do_mouse => 1,
        do_pig => 0,
        do_chimpanzee => 1,
        do_rat => 0,
        do_zebrafish => 0,
#################
#        Everything below should not need modification
#################
        pipeline_name => 'tsl_appris_'.$self->o('ensembl_release'),
        enscode_root_dir => '', #path to your perl modules
        ensembl_analysis_dir => catdir($self->o('enscode_root_dir'), 'ensembl-analysis'),
        production_dir => catdir($self->o('enscode_root_dir'), 'ensembl-production'),
        human_cs_version => '38',
        mouse_cs_version => '38',
        human_cs_name => 'GRCh'.$self->o('human_cs_version'),
        mouse_cs_name => 'GRCm'.$self->o('mouse_cs_version'),
        human_alias => 'homo_sapiens',
        mouse_alias => 'mus_musculus',
        rat_cs_name => 'Rnor_6.0',
        rat_cs_version => 6,
        zebrafish_cs_name => 'GRCz10',
        zebrafish_cs_version => 10,
        chimpanzee_cs_name => 'CHIMP2.1.4',
        chimpanzee_cs_version => 3,
        pig_cs_name => 'Sscrofa11.1',
        pig_cs_version => 111,
        rat_alias => 'rattus_norvegicus',
        zebrafish_alias => 'danio_rerio',
        chimpanzee_alias => 'pan_troglotydes',
        pig_alias => 'sus_scrofa',
        tsl_ftp_base => 'http://hgwdev.cse.ucsc.edu/~markd/gencode/tsl-handoff',
        appris_ftp_base => 'http://apprisws.bioinfo.cnio.es/forEnsembl',

        load_tsl_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'import_transcript_support_levels.pl'),
        load_appris_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'import_appris.pl'),

        staging1_db => {
            -host   => 'mysql-ens-sta-1',
            -port   => 4519,
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('hive_driver'),
        },
        human_ensembl_db => {
            -dbname => $self->o('human_alias').'_core_'.$self->o('ensembl_release').'_'.$self->o('human_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver   => $self->o('staging1_db', '-driver'),
        },
        mouse_ensembl_db => {
            -dbname => $self->o('mouse_alias').'_core_'.$self->o('ensembl_release').'_'.$self->o('mouse_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        human_refseq_db => {
            -dbname => $self->o('human_alias').'_otherfeatures_'.$self->o('ensembl_release').'_'.$self->o('human_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        mouse_refseq_db => {
            -dbname => $self->o('mouse_alias').'_otherfeatures_'.$self->o('ensembl_release').'_'.$self->o('mouse_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        production_db => {
            -dbname => 'ensembl_production_'.$self->o('ensembl_release'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user_r'),
            -pass   => $self->o('pass_r'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        zebrafish_ensembl_db => {
            -dbname => $self->o('zebrafish_alias').'_core_'.$self->o('ensembl_release').'_'.$self->o('zebrafish_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        rat_ensembl_db => {
            -dbname => $self->o('rat_alias').'_core_'.$self->o('ensembl_release').'_'.$self->o('rat_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        pig_ensembl_db => {
            -dbname => $self->o('pig_alias').'_core_'.$self->o('ensembl_release').'_'.$self->o('pig_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        chimpanzee_ensembl_db => {
            -dbname => $self->o('chimpanzee_alias').'_core_'.$self->o('ensembl_release').'_'.$self->o('chimpanzee_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
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

  my @tsl_list;
  my @appris_list;
  foreach my $species ('human', 'mouse') {
    if ($self->o('do_'.$species)) {
      push(@tsl_list, [$species, $self->o($species.'_cs_name'), $self->o($species.'_gencode_version'), $self->o($species.'_ensembl_db')]);
    }
  }
  foreach my $species ('human', 'mouse', 'pig', 'zebrafish', 'chimpanzee', 'rat') {
    if ($self->o('do_'.$species)) {
      push(@appris_list, [$species, first_upper_case($self->o($species.'_alias')), $self->o($species.'_cs_name'), $self->o($species.'_ensembl_db')]);
    }
  }
  my @analysis = (
    {
      -logic_name => 'create_working_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'mkdir -p #working_dir#',
      },
      -input_ids => [{working_dir => $self->o('working_dir')}],
      -flow_into => {
          '1' => ['create_tsl_directory', 'create_appris_directory'],
      },
    },

    {
      -logic_name => 'create_tsl_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
        cmd => 'EXIT_CODE=1; if [ -e "'.$self->o('working_dir').'/TSL_'.$self->o('ensembl_release').'" ] ;then EXIT_CODE=4; else mkdir "'.$self->o('working_dir').'/TSL_'.$self->o('ensembl_release').'"; EXIT_CODE=$?; fi; exit $EXIT_CODE',
        return_codes_2_branches => {4 => 2},
      },
      -flow_into => {
        1 => ['create_tsl_species_input'],
      },
    },

    {
      -logic_name => 'create_tsl_species_input',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputlist => \@tsl_list,
          column_names => ['species', 'cs_name', 'gencode_version', 'target_db'],
      },
      -flow_into => {
        2 => ['download_tsl'],
      },
    },

    {
      -logic_name => 'download_tsl',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'cd '.catfile($self->o('working_dir'), 'TSL_'.$self->o('ensembl_release')).'; wget -qq "'.$self->o('tsl_ftp_base').'/gencode.v#gencode_version#.transcriptionSupportLevel.tsv.gz"; gunzip gencode.v#gencode_version#.transcriptionSupportLevel.tsv.gz',
      },
      -flow_into => {
        1 => ['load_tsl'],
      },
    },

    {
      -logic_name => 'load_tsl',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
        cmd => 'perl '.$self->o('load_tsl_script').
          ' -h #expr(#target_db#->{-host})expr#'.
          ' -D #expr(#target_db#->{-dbname})expr#'.
          ' -P #expr(#target_db#->{-port})expr#'.
          ' -u #expr(#target_db#->{-user})expr#'.
          ' -p #expr(#target_db#->{-pass})expr#'.
          ' -path #cs_name#'.
          ' -write -verbose'.
          ' -file '.catfile($self->o('working_dir'), 'TSL_'.$self->o('ensembl_release'), 'gencode.v#gencode_version#.transcriptionSupportLevel.tsv'),
      },
    },

    {
      -logic_name => 'create_appris_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'EXIT_CODE=1; if [ -e "'.catdir($self->o('working_dir'), 'appris_'.$self->o('ensembl_release')).'" ] ;then EXIT_CODE=4; else mkdir "'.catdir($self->o('working_dir'), 'appris_'.$self->o('ensembl_release')).'"; EXIT_CODE=$?; fi; exit $EXIT_CODE',
          return_codes_2_branches => {4 => 2},
      },
      -flow_into => {
          '1' => ['create_appris_species_input'],
      },
    },

    {
      -logic_name => 'create_appris_species_input',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputlist => \@appris_list,
          column_names => ['species', 'species_alias', 'cs_name', 'target_db'],
      },
      -flow_into => {
          '2' => ['download_appris'],
      },
    },

    {
      -logic_name => 'download_appris',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'wget -P'.catfile($self->o('working_dir'), 'appris_'.$self->o('ensembl_release')).' -qq "'.$self->o('appris_ftp_base').'/#species_alias#.#cs_name#.e'.$self->o('ensembl_release').'appris_data.principal.txt"',
          return_codes_2_branches => {'8' => 2},
      },
      -flow_into => {
          '1' => ['load_appris'],
      },
    },

    {
      -logic_name => 'load_appris',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'perl '.$self->o('load_appris_script').
            ' -h #expr(#target_db#->{-host})expr#'.
            ' -D #expr(#target_db#->{-dbname})expr#'.
            ' -P #expr(#target_db#->{-port})expr#'.
            ' -u #expr(#target_db#->{-user})expr#'.
            ' -p #expr(#target_db#->{-pass})expr#'.
            ' -cs_version #cs_name#'.
            ' -write -verbose'.
            ' -infile '.catfile($self->o('working_dir'), 'appris_'.$self->o('ensembl_release'), '#species_alias#.#cs_name#.e'.$self->o('ensembl_release').'appris_data.principal.txt'),
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
