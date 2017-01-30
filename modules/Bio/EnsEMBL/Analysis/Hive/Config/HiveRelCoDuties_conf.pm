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

Bio::EnsEMBL::Analysis::Hive::Config::HiveRelCoDuties_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::HiveRelCoDuties_conf;

use strict;
use warnings;
use File::Spec::Functions qw(catfile catdir);

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
        release => 87, # Use it on the commandline: -release XX
        #working_dir => '', # Use it on the command line or uncomment: -working_dir /path/to/my/scratch/relco_duties_87
        #user => '', # Use it on the command line or uncomment: -user mysql_rw_user
        #password => '', # Use it on the command line or uncomment: -password mysql_rw_password
        #user_r => '', # Use it on the command line or uncomment: -user_r mysql_ro_user
        pass_r => undef, # (optional) Use it on the command line: -password mysql_ro_password
        human_gencode_version => '23', # Use it on the command line: -human_gencode_version XX
        mouse_gencode_version => 'M12', # Use it on the command line: -password mysql_rw_password
#################
#        Everything below should not need modification
#################
        pipeline_name => 'relco_duties_'.$self->o('release'),
        ensembl_root_dir => $ENV{ENSCODE},
        ensembl_analysis_dir => catdir($self->o('ensembl_root_dir'), 'ensembl-analysis'),
        production_dir => catdir($self->o('ensembl_root_dir'), 'ensembl-production'),
        human_cs_version => '38',
        mouse_cs_version => '38',
        human_cs_name => 'GRCh'.$self->o('human_cs_version'),
        mouse_cs_name => 'GRCm'.$self->o('mouse_cs_version'),
        human_alias => 'homo_sapiens',
        mouse_alias => 'mus_musculus',
        rat_cs_name => 'Rnor_6.0',
        zebrafish_cs_name => 'GRCz10',
        chimpanzee_cs_name => 'CHIMP2.1.4',
        pig_cs_name => 'Sscrofa10.2',
        rat_alias => 'rattus_norvegicus',
        zebrafish_alias => 'danio_rerio',
        chimpanzee_alias => 'pan_troglotydes',
        pig_alias => 'sus_scrofa',
        samtools => 'samtools', # It might be better to give the absolute path, but I think it's ok
        webdev_nfs => '/nfs/ensnfs/webdev/staging',
        blastdb_basedir => $ENV{BLASTDB_DIR}, # This can be anywhere, the path has to be given to the peson doing the cDNA update
        refseq_ftp_basedir => 'ftp://ftp.ncbi.nlm.nih.gov/refseq',
        tsl_ftp_base => 'http://hgwdev.cse.ucsc.edu/~markd/gencode/tsl-handoff',
        appris_ftp_base => 'http://apprisws.bioinfo.cnio.es/forEnsembl',
        check_datafiles_script => catfile($self->o('production_dir'), 'scripts', 'datafiles', 'check_datafiles.pl'),
        compare_gencode_refseq_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'compare_gencode_refseq_genes.pl'),
        label_gencode_basic_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'label_gencode_basic_transcripts.pl'),
        load_tsl_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'import_transcript_support_levels.pl'),
        load_appris_script => catfile($self->o('ensembl_analysis_dir'), 'scripts', 'Merge', 'import_appris.pl'),
        staging1_db => {
            -host   => 'ens-staging1',
            -port   => 3306,
            -user   => $self->o('user_r'),
            -pass   => $self->o('pass_r'),
            -driver => $self->o('hive_driver'),
        },
        staging2_db => {
            -host   => 'ens-staging2',
            -port   => 3306,
            -user   => $self->o('user_r'),
            -pass   => $self->o('pass_r'),
            -driver => $self->o('hive_driver'),
        },
        human_ensembl_db => {
            -dbname => $self->o('human_alias').'_core_'.$self->o('release').'_'.$self->o('human_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver   => $self->o('staging1_db', '-driver'),
        },
        mouse_ensembl_db => {
            -dbname => $self->o('mouse_alias').'_core_'.$self->o('release').'_'.$self->o('mouse_cs_version'),
            -host   => $self->o('staging2_db', '-host'),
            -port   => $self->o('staging2_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging2_db', '-driver'),
        },
        human_refseq_db => {
            -dbname => $self->o('human_alias').'_otherfeatures_'.$self->o('release').'_'.$self->o('human_cs_version'),
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user_r'),
            -pass   => $self->o('pass_r'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        mouse_refseq_db => {
            -dbname => $self->o('mouse_alias').'_otherfeatures_'.$self->o('release').'_'.$self->o('mouse_cs_version'),
            -host   => $self->o('staging2_db', '-host'),
            -port   => $self->o('staging2_db', '-port'),
            -user   => $self->o('user_r'),
            -pass   => $self->o('pass_r'),
            -driver => $self->o('staging2_db', '-driver'),
        },
        production_db => {
            -dbname => 'ensembl_production',
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user_r'),
            -pass   => $self->o('pass_r'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        zebrafish_ensembl_db => {
            -dbname => $self->o('zebrafish_alias').'_core_'.$self->o('release').'_10',
            -host   => $self->o('staging1_db', '-host'),
            -port   => $self->o('staging1_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging1_db', '-driver'),
        },
        rat_ensembl_db => {
            -dbname => $self->o('rat_alias').'_core_'.$self->o('release').'_6',
            -host   => $self->o('staging2_db', '-host'),
            -port   => $self->o('staging2_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging2_db', '-driver'),
        },
        pig_ensembl_db => {
            -dbname => $self->o('pig_alias').'_core_'.$self->o('release').'_102',
            -host   => $self->o('staging2_db', '-host'),
            -port   => $self->o('staging2_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging2_db', '-driver'),
        },
        chimpanzee_ensembl_db => {
            -dbname => $self->o('chimpanzee_alias').'_core_'.$self->o('release').'_214',
            -host   => $self->o('staging2_db', '-host'),
            -port   => $self->o('staging2_db', '-port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -driver => $self->o('staging2_db', '-driver'),
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

  my @analysis = (
    {
      -logic_name => 'create_working_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'mkdir #working_dir#',
      },
      -input_ids => [{working_dir => $self->o('working_dir')}],
      -flow_into => {
          '1' => ['create_species_input', 'create_tsl_directory', 'create_appris_directory', 'check_refseq_directory', 'create_datafile_input'],
      },
    },

    {
      -logic_name => 'check_refseq_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'EXIT_CODE=1;REFSEQ_DATE="RefSeq_"`date "+%Y_%m"`; if [ -e "'.$self->o('blastdb_basedir').'/$REFSEQ_DATE" ] ;then echo "'.$self->o('blastdb_basedir').'/$REFSEQ_DATE already exists"; EXIT_CODE=4; else mkdir "'.$self->o('blastdb_basedir').'/$REFSEQ_DATE"; EXIT_CODE=$?;fi; exit $EXIT_CODE',
          return_codes_2_branches => {4 => 2},
      },
      -flow_into => {
          '1' => ['create_refseq_species_input'],
      },
    },

    {
      -logic_name => 'create_refseq_species_input',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputlist => [['human', 'H_sapiens'],
                        ['mouse', 'M_musculus']],
          column_names => ['species', 'species_alias'],
      },
      -flow_into => {
          '2' => ['download_refseq_cdnas'],
      },
    },

    {
      -logic_name => 'download_refseq_cdnas',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'REFSEQ_DATE="RefSeq_"`date "+%Y_%m"`; cd "'.$self->o('blastdb_basedir').'/$REFSEQ_DATE"; wget -qq "'.$self->o('refseq_ftp_basedir').'/#species_alias#/mRNA_Prot/#species#.*.rna.fna.gz"; gunzip #species#*.gz; cat #species#*.rna.fna > #species#.fna;COUNTT=`grep -c \> #species#.fna | cut -f1`;COUNTA=0;for F in #species#*.rna.fna;do COUNTA=$((COUNTA+`grep -c \> $F | cut -f1`));done; if [ $COUNTT -ne $COUNTA ]; then exit 2;fi; LASTREFSEQ=`ls -td /data/blastdb/Ensembl/RefSeq* | head -n 2 | tail -n 1`; if [ $COUNTT -gt "`grep -c \> $LASTREFSEQ | cut -f1`" ]; then rm -f #species#*.rna.fna; else exit 2; fi',
      },
    },

    {
      -logic_name => 'create_species_input',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputlist => [['human', $self->o('human_ensembl_db'), $self->o('human_cs_name'), $self->o('human_refseq_db'), 'refseq_human_import'],
                        ['mouse', $self->o('mouse_ensembl_db'), $self->o('mouse_cs_name'), $self->o('mouse_refseq_db'), 'refseq_mouse_import']],
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
          '1' => ['create_top_level_input_ids', 'label_gencode_basic'],
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
                      2 => {'compare_gencode_refseq' => {'iid' => '#iid#', species => '#species#', target_db => '#target_db#', refseq_db => '#refseq_db#', cs_name => '#cs_name#', refseq_logicname => '#refseq_logicname#'}},
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
            ' -refseqdbname #expr(#refseq_db#->{-dbname})expr#'.
            ' -coord_system_version #cs_name#'.
            ' -refseq_logicname #refseq_logicname#'.
            ' -chr_num #iid# -store',
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
    },

    {
      -logic_name => 'create_tsl_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'EXIT_CODE=1; if [ -e "'.$self->o('working_dir').'/TSL_'.$self->o('release').'" ] ;then EXIT_CODE=4; else mkdir "'.$self->o('working_dir').'/TSL_'.$self->o('release').'"; EXIT_CODE=$?; fi; exit $EXIT_CODE',
          return_codes_2_branches => {4 => 2},
      },
      -flow_into => {
          '1' => ['create_tsl_species_input'],
      },
    },

    {
      -logic_name => 'create_tsl_species_input',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputlist => [['human', $self->o('human_cs_name'), $self->o('human_gencode_version'), $self->o('human_ensembl_db')],
                        ['mouse', $self->o('mouse_cs_name'), $self->o('mouse_gencode_version'), $self->o('mouse_ensembl_db')]],
          column_names => ['species', 'cs_name', 'gencode_version', 'target_db'],
      },
      -flow_into => {
          '2' => ['download_tsl'],
      },
    },

    {
      -logic_name => 'download_tsl',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'cd '.catfile($self->o('working_dir'), 'TSL_'.$self->o('release')).'; wget -qq "'.$self->o('tsl_ftp_base').'/gencode.v#gencode_version#.transcriptionSupportLevel.tsv.gz"; gunzip gencode.v#gencode_version#.transcriptionSupportLevel.tsv.gz',
      },
      -flow_into => {
          '1' => ['load_tsl'],
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
            ' -file '.catfile($self->o('working_dir'), 'TSL_'.$self->o('release'), 'gencode.v#gencode_version#.transcriptionSupportLevel.tsv'),
      },
    },

    {
      -logic_name => 'create_appris_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'EXIT_CODE=1; if [ -e "'.catdir($self->o('working_dir'), 'appris_'.$self->o('release')).'" ] ;then EXIT_CODE=4; else mkdir "'.catdir($self->o('working_dir'), 'appris_'.$self->o('release')).'"; EXIT_CODE=$?; fi; exit $EXIT_CODE',
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
          inputlist => [['human', first_upper_case($self->o('human_alias')), $self->o('human_cs_name'), $self->o('human_ensembl_db')],
                        ['pig', first_upper_case($self->o('pig_alias')), $self->o('pig_cs_name'), $self->o('pig_ensembl_db')],
                        ['rat', first_upper_case($self->o('rat_alias')), $self->o('rat_cs_name'), $self->o('rat_ensembl_db')],
                        ['chimpanzee', first_upper_case($self->o('chimpanzee_alias')), $self->o('chimpanzee_cs_name'), $self->o('chimpanzee_ensembl_db')],
                        ['zebrafish', first_upper_case($self->o('zebrafish_alias')), $self->o('zebrafish_cs_name'), $self->o('zebrafish_ensembl_db')],
                        ['mouse', first_upper_case($self->o('mouse_alias')), $self->o('mouse_cs_name'), $self->o('mouse_ensembl_db')]],
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
          cmd => 'wget -P'.catfile($self->o('working_dir'), 'appris_'.$self->o('release')).' -qq "'.$self->o('appris_ftp_base').'/#species_alias#.#cs_name#.e'.$self->o('release').'appris_data.principal.txt"',
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
            ' -infile '.catfile($self->o('working_dir'), 'appris_'.$self->o('release'), '#species_alias#.#cs_name#.e'.$self->o('release').'appris_data.principal.txt'),
      },
    },

    {
      -logic_name => 'create_datafile_input',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputlist => [[$self->o('staging1_db', '-host'), $self->o('staging1_db', '-user'), $self->o('staging1_db', '-port')],
                        [$self->o('staging2_db', '-host'), $self->o('staging2_db', '-user'), $self->o('staging2_db', '-port')]],
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
          cmd => 'perl '.$self->o('check_datafiles_script').' -release '.$self->o('release').' -host #host# -user #user# -port #port# -datafile_dir '.$self->o('webdev_nfs').' --force_samtools_binary --samtools_binary '.$self->o('samtools').' 1>/dev/null',
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
