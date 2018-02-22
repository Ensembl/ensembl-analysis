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

Bio::EnsEMBL::Analysis::Config::VertRNA_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::VertRNA_conf;

use strict;
use warnings;

use File::Spec::Functions;

use Bio::EnsEMBL::Hive::Version 2.4;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf; # This is needed for WHEN() ELSE
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');


=head2 default_options

    Description : Interface method that should return a hash of option_name->default_option_value pairs.
                  Please see existing PipeConfig modules for examples.

=cut

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },

        base_ena_ftp  => $ENV{'BASE_ENA_FTP'} || 'ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/std',
        vertrna_file  => $ENV{'VERTRNA_FILE'} || 'embl_vertrna-1',
        embl2fasta_script  => $ENV{'EMBL2FASTA_SCRIPT'} || catfile($ENV{ENSCODE}, 'ensembl-analysis', 'scripts', 'databases', 'EMBL2fasta.pl'),
        vertrna_version  => $ENV{'VERTRNA_VERSION'},
        slack_token   => $ENV{SLACK_GENEBUILD},
        slack_hook_url => 'https://hooks.slack.com/services/'.$self->o('slack_token'),
    };
}


sub resource_classes {
    my $self = shift @_;

    my $resources = $self->SUPER::resource_classes();
    $resources->{'default'}->{LSF} = '-M4000 -R"select[mem>4000] rusage[mem=4000]"';
    return $resources;
}


sub pipeline_analyses {
    my ($self) = @_;

    my @analyses = (
        {
          -logic_name => 'setup_directory',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd     => 'mkdir #vertrna_dir#; which lfs &> /dev/null; if [ $? -eq 0 ]; then lfs getstripe #vertrna_dir# &> /dev/null; if [ $? -eq 0 ];then lfs setstripe -c -1 #vertrna_dir#;fi;fi',
          },
          -input_ids  => [
            { vertrna_dir => $self->o('vertrna_dir'), blast_type => $self->o('blast_type'), vertrna_file => $self->o('vertrna_file')},
          ],
          -flow_into => {
          '1' => ['create_tax_ids'],
          },
        },

        {
          -logic_name => 'create_tax_ids',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
            inputlist    => ['hum', 'vrt', 'mus', 'mam', 'rod', 'syn'],
            column_names => ['tax_id'],
          },
          -flow_into => {
            '2->A' => ['create_type_ids'],
            'A->1' => ['concat_file'],
          },
        },

        {
          -logic_name => 'create_type_ids',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
            inputlist     => ['std', 'htc'],
            column_names => ['type'],
          },
          -analysis_capacity  => 6,
          -flow_into => {
            2 => {'download_files' => {tax_id => '#tax_id#', type => '#type#'}},
          },
        },

        {
          -logic_name => 'download_files',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'wget -q -P '.$self->o('vertrna_dir').' "'.$self->o('base_ena_ftp').'/rel_#type#_#tax_id#_*_r'.$self->o('vertrna_version').'.dat.gz"',
          },
          -analysis_capacity  => 12,
          -flow_into => {
            1 => ['create_file_ids'],
          },
        },

        {
          -logic_name => 'create_file_ids',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
            inputcmd     => 'ls '.$self->o('vertrna_dir').'/rel_#type#_#tax_id#_*_r'.$self->o('vertrna_version').".dat.gz | sed 's/.dat.gz//'",
            column_names => ['input_file'],
          },
          -analysis_capacity  => 12,
          -flow_into => {
            2 => ['gunzip'],
          },
        },

        {
          -logic_name => 'gunzip',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'gunzip #input_file#.dat.gz',
          },
          -analysis_capacity  => 12,
          -flow_into => {
            1 => [ 'embl2fasta' ],
          },
        },

        {
          -logic_name => 'embl2fasta',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'perl '.$self->o('embl2fasta_script').' -vertrna -embl #input_file#.dat -fasta #input_file#',
          },
          -analysis_capacity  => 30,
          -flow_into => {
            1 => ['delete_embl_file'],
          },
        },

        {
          -logic_name => 'delete_embl_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'rm #input_file#.dat',
          },
          -analysis_capacity  => 12,
          -flow_into => {
            1 => ['?accu_name=input_file&accu_input_variable=input_file&accu_address=[]'],
          },
        },

        {
          -logic_name => 'concat_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'cat #expr(join(" ", @{#input_file#}))expr# > #vertrna_dir#/#vertrna_file#; rm #expr(join(" ", @{#input_file#}))expr#',
          },
          -flow_into => {
            1 => ['check_against_previous'],
          }
        },

        {
          -logic_name => 'check_against_previous',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'DIR="#vertrna_dir#";VERSION=`ls $DIR/.. | sort -nr | head -n 2 | tail -n1`;OLDDIR="${DIR}/../$VERSION";if [ -d "$OLDDIR" ]; then OLD=`grep -c \> $OLDDIR/#vertrna_file#`;NEW=`grep -c \> #vertrna_dir#/#vertrna_file#`; if [ $NEW -lt $OLD ];then exit 2;fi;fi',
          },
          -flow_into => WHEN( '#blast_type# eq "ncbi"' => 'ncbi_format', ELSE 'wu_format',),
        },

        {
          -logic_name => 'ncbi_format',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'makeblastdb -dbtype nucl -in #vertrna_dir#/#vertrna_file# -title "EMBL VertRNA"',
          },
          -flow_into => {
            1 => ['finalise_dbs'],
          }
        },

        {
          -logic_name => 'wu_format',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'xdformat -n #vertrna_dir#/#vertrna_file# -t "EMBL VertRNA"',
          },
          -flow_into => {
            1 => ['finalise_dbs'],
          }
        },

        {
          -logic_name => 'finalise_dbs',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'chmod g-w #vertrna_dir#; find #vertrna_dir# -type f -execdir chmod g-w {} \;',
          },
          -flow_into => {
            1 => ['propagate_the_word'],
          }
        },

        {
          -logic_name => 'propagate_the_word',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SlackNotification',
          -parameters => {
            slack_webhook => $self->o('slack_hook_url'),
            text => 'The VertRNA blast database has been updated to '.$self->o('vertrna_version'),
          },
        },
    );

    foreach my $analysis (@analyses) {
        $analysis->{'-max_retry_count'} = 0 unless (exists $analysis->{'-max_retry_count'});
    }
    return \@analyses;
}

1;
