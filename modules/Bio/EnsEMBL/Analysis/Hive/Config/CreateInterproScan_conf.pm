=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

=head1 DESCRIPTION

This configuration will create a working Interproscan with Panther data and Interpro2GO

=cut

package CreateInterproScan_conf;

use strict;
use warnings;

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
        pipeline_name => 'hive_interproscan_'.`printf "%s" "\`date +%Y_%m_%d\`"`,
        base_directory  => '/software/ensembl/genebuild/bin/interproscan',
        working_dir => '',
        interproscan_http => 'https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5',
        interproscan_version => '',
        tmhmm_model_path => '/software/ensembl/genebuild/bin/tmhmm/tmhmm-2.0c/lib/TMHMM2.0.model',
        binary_signalp_path => '/software/ensembl/genebuild/bin/signalp-4.1/signalp',
        signalp_perl_library_dir => '/software/ensembl/genebuild/bin/signalp-4.1/lib',
        binary_tmhmm_path => '/software/ensembl/genebuild/bin/tmhmm/tmhmm-2.0c/bin/decodeanhmm',
        temp_dir => '/tmp',
        java_path => '/software/java/bin/java',
        interproscan_conf => '/software/ensembl/genebuild/lib/eg_protein_annotation/InterProScanSeg_conf.pm',
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
      -logic_name => 'check_directories',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -meadow_type => 'LOCAL',
      -parameters => {
          cmd => 'EXIT_CODE=0;COUNT=0;for F in '.$self->o('tmhmm_model_path').' '.$self->o('signalp_perl_library_dir').' '.$self->o('temp_dir').' '.$self->o('interproscan_conf').'; do if [ ! -e "$F" ]; then echo "$F does not exist";EXIT_CODE=1;fi;COUNT=$((COUNT+1));done; if [ $COUNT -ne 4 ];then echo "Some files are not set";EXIT_CODE=1;fi;COUNT=0; for E in '.$self->o('binary_tmhmm_path').' '.$self->o('binary_signalp_path').' '.$self->o('java_path').'; do if [ ! -x "$E" ]; then echo "$E is not executable";EXIT_CODE=1;fi;COUNT=$((COUNT+1));done; if [ $COUNT -ne 3 ];then echo "Some executables are not set";EXIT_CODE=1;fi;COUNT=0; for D in '.$self->o('base_directory').' '.$self->o('working_dir').'; do if [ ! -d "$D" ]; then echo "Creating $D";mkdir "$D";EXIT_CODE=1;fi;COUNT=$((COUNT+1));done; if [ $COUNT -ne 2 ];then echo "Some directories are not set";EXIT_CODE=1;fi;exit $EXIT_CODE',
      },
      -input_ids => [{filename => 'interproscan-'.$self->o('interproscan_version').'-64-bit.tar.gz'}],
      -flow_into => {
          '1' => {'download_interproscan' => {filename => '#filename#', iprscan_dir => $self->o('working_dir').'/interproscan-'.$self->o('interproscan_version')}},
      },
    },
    {
      -logic_name => 'download_interproscan',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'wget -q -P '.$self->o('working_dir').' "'.$self->o('interproscan_http').'/'.$self->o('interproscan_version').'/#filename#" "'.$self->o('interproscan_http').'/'.$self->o('interproscan_version').'/#filename#.md5"',
      },
      -flow_into => {
          '1' => ['check_downloaded_file'],
      },
    },
    {
      -logic_name => 'check_downloaded_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'cd '.$self->o('working_dir').'; md5sum -c #filename#.md5',
      },
      -flow_into => {
          '1' => ['decompress_file'],
      },
    },
    {
      -logic_name => 'decompress_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'cd '.$self->o('working_dir').'; tar -zxvf #filename#',
      },
      -flow_into => {
          '1->A' => ['create_input_panther_file', 'download_interpro2go'],
          'A->1' => ['update_properties'],
      },
    },
    {
      -logic_name => 'create_input_panther_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name => 'default',
      -parameters => {
          inputcmd => 'grep "panther.signature.library.release" #iprscan_dir#/interproscan.properties | awk -F= \'{print "panther-data-"$2".tar.gz"}\'',
          column_names => ['filename'],
      },
      -flow_into => {
          '2' => {'download_panther_file' => {filename => '#filename#', iprscan_dir => '#iprscan_dir#'}},
      },
    },
    {
      -logic_name => 'download_panther_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'wget -q -P #iprscan_dir#/data "'.$self->o('interproscan_http').'/data/#filename#" "'.$self->o('interproscan_http').'/data/#filename#.md5"',
      },
      -flow_into => {
          '1' => ['check_panther_file'],
      },
    },
    {
      -logic_name => 'check_panther_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'cd #iprscan_dir#/data; md5sum -c #filename#.md5',
      },
      -flow_into => {
          '1' => ['decompress_panther_file'],
      },
    },
    {
      -logic_name => 'decompress_panther_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'cd #iprscan_dir#/data; tar -zxvf #filename#',
      },
    },
    {
      -logic_name => 'download_interpro2go',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -meadow_type => 'LOCAL',
      -rc_name => 'default',
      -parameters => {
          cmd => 'wget -q -N -P '.$self->o('base_directory').'/interpro2go "http://geneontology.org/external2go/interpro2go"',
      },
      -flow_into => {
          '1' => ['process_interpro2go'],
      },
    },
    {
      -logic_name => 'process_interpro2go',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -meadow_type => 'LOCAL',
      -parameters => {
          cmd => 'cd '.$self->o('base_directory').'/interpro2go; INTERPRO2GO_DIR=`ls --full-time interpro2go |awk \'{print $6}\' | sed \'s/-/_/g\'`; mkdir ${INTERPRO2GO_DIR};mv interpro2go ${INTERPRO2GO_DIR}/',
      },
    },
    {
      -logic_name => 'update_properties',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'sed -i.bak \'s"\(tmhmm.model.path\)=.*"\1='.$self->o('tmhmm_model_path').'";s"\(binary.tmhmm.path\)=.*"\1='.$self->o('binary_tmhmm_path').'";s"\(binary.signalp.path\)=.*"\1='.$self->o('binary_signalp_path').'";s"\(signalp.perl.library.dir\)=.*"\1='.$self->o('signalp_perl_library_dir').'";s"\(temporary.file.directory\)=\w\+"\1='.$self->o('temp_dir').'";s"\(panther.temporary.file.directory\)=.*"\1=temp";s"\(precalculated.match.lookup.service.url\)=.*"\1="\' #iprscan_dir#/interproscan.properties',
      },
      -flow_into => {
          '1' => ['update_java_path'],
      },
    },
    {
      -logic_name => 'update_java_path',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'sed -i.bak \'s"JAVA=.*"JAVA='.$self->o('java_path').'"\' #iprscan_dir#/interproscan.sh',
      },
      -flow_into => {
          '1' => ['test_interproscan'],
      },
    },
    {
      -logic_name => 'test_interproscan',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => '2GB',
      -parameters => {
          cmd => 'cd #iprscan_dir#; ./interproscan.sh -i test_proteins.fasta -f tsv',
      },
      -flow_into => {
          '1' => ['clean_directory'],
      },
    },
    {
      -logic_name => 'clean_directory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
          cmd => 'cd '.$self->o('working_dir').'; rm interproscan*.gz* #iprscan_dir#/data/panther*gz* #iprscan_dir#/interproscan*bak',
      },
      -flow_into => {
          '1' => ['copy_to_basedir'],
      },
    },
    {
      -logic_name => 'copy_to_basedir',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -meadow_type => 'LOCAL',
      -rc_name => 'default',
      -parameters => {
          cmd => 'cp -r #iprscan_dir# '.$self->o('base_directory'),
      },
      -flow_into => {
          '1' => ['test_installed_interproscan'],
      },
    },
    {
      -logic_name => 'test_installed_interproscan',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => '2GB',
      -parameters => {
          cmd => 'cd '.$self->o('base_directory').'/interproscan-'.$self->o('interproscan_version').'; ./interproscan.sh -i test_proteins.fasta -f tsv -o ${HOME}/interpro_'.$self->o('interproscan_version').'_test.out',
      },
      -flow_into => {
          '1' => ['delete_working_dir'],
      },
    },
    {
      -logic_name => 'delete_working_dir',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -meadow_type => 'LOCAL',
      -rc_name => 'default',
      -parameters => {
          cmd => 'rm -r #iprscan_dir#; rm ${HOME}/interpro_'.$self->o('interproscan_version').'_test.out',
      },
    },
  );
  foreach my $analyses (@analysis) {
      $analyses->{-max_retry_count} = 0 unless (exists $analyses->{-max_retry_count});
  }
  return \@analysis;
}


=head2 resource_classes

 Arg [1]    : None
 Description: Resources needed for the pipeline, it uses the default one and one requesting 2GB when testing interproscan
 Returntype : Hashref
 Exceptions : None

=cut

sub resource_classes {
    my $self = shift;

    return {
        %{ $self->SUPER::resource_classes() },  # inherit other stuff from the base class
      '2GB' => { LSF => '-M2000 -R"select[mem>2000] rusage[mem=2000]"'},
    };
}

1;
