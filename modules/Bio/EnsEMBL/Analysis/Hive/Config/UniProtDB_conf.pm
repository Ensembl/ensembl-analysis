=pod

=head1 NAME

    Bio::EnsEMBL::Analysis::Hive::Config::UniProtDB_conf

=head1 DESCRIPTION

    Configuration to create multiple blast databases which can be shared by Enesmbl and Havana
    It has to be used in conjonction with the create_uniprot_database.sh in ensembl-analysis/scripts

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2018] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    Please subscribe to the Hive mailing list:  http://listserver.ebi.ac.uk/mailman/listinfo/ehive-users  to discuss Hive-related questions or to be notified of our updates
    For any other matter, please subscribe to the Dev mailing list: http://lists.ensembl.org/mailman/listinfo/dev

=cut


package Bio::EnsEMBL::Analysis::Hive::Config::UniProtDB_conf;

use strict;
use warnings;

use 5.014000;
use Bio::EnsEMBL::ApiVersion ();
use File::Spec::Functions;

use  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf; # This is needed for WHEN() ELSE
use parent ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');


=head2 default_options

    Description : Interface method that should return a hash of option_name->default_option_value pairs.
                  Please see existing PipeConfig modules for examples.

=cut

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },

        pipeline_name     => $ENV{pipeline_name},
        base_uniprot_ftp  => $ENV{'BASE_UNIPROT_FTP'} || 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase',
        base_uniprot_url  => $ENV{'BASE_UNIPROT_URL'} || 'http://www.uniprot.org/uniprot/?query=',
        uniprot_vert_file => $ENV{'UNIPROT_VERT_FILE'} || 'uniprot_vertebrata',
        uniprot_file      => $ENV{'UNIPROT_FILE'} || 'uniprot',
        uniprot_dir       => $ENV{'UNIPROT_DIR'} || catdir($ENV{BLASTDB_DIR}, 'uniprot', 'uniprot_'.$self->o('uniprot_version')),
        varsplic_file     => $ENV{'VARSPLIC_FILE'} || 'uniprot_sprot_varsplic',
        fasta_suffix      => $ENV{'FASTA_SUFFIX'} || '.fasta',
        uniprot_version   => $ENV{'UNIPROT_VERSION'},
        uniprot_date      => $ENV{'UNIPROT_DATE'} || localtime,
        indicate          => 'indicate',
        blast_type        => 'ncbi', # Should be either 'ncbi' or 'wu'
        slack_token       => $ENV{SLACK_GENEBUILD},
        slack_hook_url    => 'https://hooks.slack.com/services/'.$self->o('slack_token'),
        embl2fasta_script => $ENV{'EMBL2FASTA_SCRIPT'} || catfile($ENV{ENSCODE}, 'ensembl-analysis', 'scripts', 'databases', 'EMBL2fasta.pl'),
        process_isoforms_script => $ENV{'PROCESS_ISOFORMS_SCRIPT'} || catfile($ENV{ENSCODE}, 'ensembl-analysis', 'scripts', 'databases', 'process_uniprot_isoforms.pl'),
        kill_list_host => '',
        kill_list_port => '',
        kill_list_user => '',
        kill_list_pass => '',
        kill_list_dbname => '',
    };
}


sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        base_uniprot_ftp        => $self->o('base_uniprot_ftp'),
        base_uniprot_url        => $self->o('base_uniprot_url'),
        uniprot_vert_file       => $self->o('uniprot_vert_file'),
        uniprot_file            => $self->o('uniprot_file'),
        embl2fasta_script       => $self->o('embl2fasta_script'),
        uniprot_dir             => $self->o('uniprot_dir'),
        varsplic_file           => $self->o('varsplic_file'),
        fasta_suffix            => $self->o('fasta_suffix'),
        uniprot_version         => $self->o('uniprot_version'),
        uniprot_date            => $self->o('uniprot_date'),
        process_isoforms_script => $self->o('process_isoforms_script'),
        blast_type              => $self->o('blast_type'),
    };
}


sub resource_classes {
    my $self = shift @_;

    my $resources = $self->SUPER::resource_classes();
    $resources->{'default'}->{LSF} = '-M1000 -R"select[mem>1000] rusage[mem=1000]"';
    return $resources;
}


sub pipeline_analyses {
    my ($self) = @_;
    my @analyses = (

        {
          -logic_name => 'setup_directory',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'mkdir #uniprot_dir#; which lfs &> /dev/null; if [ $? -eq 0 ]; then lfs getstripe #uniprot_dir# &> /dev/null; if [ $? -eq 0 ];then lfs setstripe -c -1 #uniprot_dir#;fi;fi',
          },
          -analysis_capacity  => 1,
          -input_ids  => [ { uniprot_dir => $self->o('uniprot_dir')} ],
          -flow_into => {
            '1->A' => ['create_fasta_ids', 'dummy_fan'],
            'A->1' => ['finalise_dbs'],
          },
        },


#################
# Fasta based DB
#################
        {
          -logic_name => 'create_fasta_ids',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
            inputlist     => [
              ['uniprot_trembl.fasta', 'http://www.ebi.ac.uk/uniprot/TrEMBLstats/'],
              ['uniprot_sprot.fasta', 'http://web.expasy.org/docs/relnotes/relstat.html'],
            ],
            column_names => ['input_file', 'check_file'],
          },
          -analysis_capacity  => 1,
          -flow_into => {
            '2->A' => ['download_fasta'],
            'A->1' => ['concat_fasta'],
          },
        },

        {
          -logic_name => 'download_fasta',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'wget -O #uniprot_dir#/#input_file#.gz -q "#base_uniprot_ftp#/complete/#input_file#.gz"',
          },
          -analysis_capacity  => 2,
          -flow_into => {
            1 => [ 'gunzip_fasta_file' ],
          },
        },

        {
          -logic_name => 'gunzip_fasta_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'gunzip -c #uniprot_dir#/#input_file#.gz | perl -ne \'if (/>/) { s/>...([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*PE=([0-9]).*SV=([0-9]{1,}).*/>$1.$4 $3/} else {s/O/K/g} print;\' > #uniprot_dir#/#input_file#',
          },
          -analysis_capacity  => 2,
          -flow_into => {
            1 => [ 'check_fasta_file' ],
          },
        },

        {
          -logic_name => 'check_fasta_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'EXIT_CODE=1; COUNT=`grep -c \> #uniprot_dir#/#input_file#`;if [ "`wget -q -O - #check_file# | perl -ne \'if(/(\d+)\s+sequence\s+entries/) {print $1}\'`" -eq "$COUNT" ]; then rm #uniprot_dir#/#input_file#.gz; echo $COUNT > #uniprot_dir#/#input_file#.count;EXIT_CODE=0; fi; exit $EXIT_CODE',
          },
          -analysis_capacity  => 2,
          -flow_into => {
            1 => [ 'pre_entry_loc', '?accu_name=filename&accu_input_variable=input_file&accu_address=[]' ],
          },
        },

        {
          -logic_name => 'pre_entry_loc',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'EXIT_CODE=1; if [ "`echo #input_file# | sed \'s/sprot//\'`" = "#input_file#" ]; then EXIT_CODE=2; else EXIT_CODE=0; fi; exit $EXIT_CODE;',
            return_codes_2_branches => { '2' => 3},
          },
          -analysis_capacity  => 1,
          -flow_into => {
            1 => [ 'entry_loc' ],
          },
        },

        {
          -logic_name => 'entry_loc',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'EXIT_CODE=1;grep \> #uniprot_dir#/#input_file# | sed \'s/>//\' | awk \'{print $1, "STD"}\' > #uniprot_dir#/entry_loc;EXIT_CODE=$?;WCCOUNT=`wc -l #uniprot_dir#/entry_loc | awk \'{print $1}\'`;COUNT=`grep -c \> #uniprot_dir#/#input_file#`; if [ $WCCOUNT -ne $COUNT ]; then EXIT_CODE=1;fi; exit $EXIT_CODE',
          },
          -analysis_capacity  => 1,
        },

        {
          -logic_name => 'concat_fasta',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'EXIT_CODE=1;cd #uniprot_dir#; cat #expr(join(" ", @{#filename#}))expr# > #uniprot_file#;EXIT_CODE=$?; if [ "`awk \'{SUM+=$1} END {print SUM}\' ./*.count`" -eq "`grep -c \> #uniprot_file#`" ]; then rm #expr(join(" ", @{#filename#}))expr#; rm -f ./*.count; EXIT_CODE=0; fi; exit $EXIT_CODE',
          },
          -analysis_capacity  => 1,
          -flow_into => {
            1 => WHEN( '#blast_type# eq "ncbi"' => 'ncbi_format_blast_uniprot', ELSE 'wu_format_blast_uniprot',),
          },
        },

        {
          -logic_name => 'wu_format_blast_uniprot',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'xdformat -p -t "#uniprot_file#" -v "#uniprot_version#" -d "#uniprot_date#" #uniprot_dir#/#uniprot_file#',
          },
          -analysis_capacity  => 1,
        },

        {
          -logic_name => 'ncbi_format_blast_uniprot',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'makeblastdb -dbtype prot -title "#uniprot_file#" -in #uniprot_dir#/#uniprot_file#',
          },
          -analysis_capacity  => 1,
        },

        {
          -logic_name => 'dummy_fan',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -analysis_capacity  => 1,
          -flow_into => {
            '1->A' => ['create_isoforms_ids', 'create_pe_level_ids'],
            'A->1' => ['check_indicate_index'],
          },
        },


#################
# Isoforms
#################
        {
          -logic_name => 'create_isoforms_ids',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
            inputlist    => [['#varsplic_file#', 'http://web.expasy.org/docs/relnotes/relstat.html']],
            column_names => ['input_file', 'check_file'],
          },
          -analysis_capacity  => 1,
          -flow_into => {
            '2->A' => ['download_isoforms'],
            'A->1' => ['create_taxonomy_ids'],
          },
        },

        {
          -logic_name => 'download_isoforms',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'wget -O #uniprot_dir#/#input_file##fasta_suffix#.gz -q "#base_uniprot_ftp#/complete/#input_file##fasta_suffix#.gz"',
          },
          -analysis_capacity  => 1,
          -flow_into => {
            1 => [ 'gunzip_isofasta_file' ],
          },
        },

        {
          -logic_name => 'gunzip_isofasta_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'gunzip -c #uniprot_dir#/#input_file##fasta_suffix#.gz | perl -ne \'if (!/>/) {s/O/K/g} print;\' > #uniprot_dir#/#input_file##fasta_suffix#',
          },
          -analysis_capacity  => 1,
          -flow_into => {
            1 => [ 'check_isofasta_file' ],
          },
        },

        {
          -logic_name => 'check_isofasta_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'EXIT_CODE=1; if [ "`wget -q -O - #check_file# | perl -ne \'if(/alternative\s+splicing.*:\s+(\d+)/) {print $1}\'`" -eq "`grep -c \> #uniprot_dir#/#input_file##fasta_suffix#`" ]; then rm #uniprot_dir#/#input_file##fasta_suffix#.gz; EXIT_CODE=0; fi; exit $EXIT_CODE',
          },
          -analysis_capacity => 1,
        },


#################
# Taxonomy based files
#################
        {
          -logic_name => 'create_taxonomy_ids',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
            inputlist   => [
              'uniprot_trembl_vertebrates',
              'uniprot_trembl_mammals',
              'uniprot_trembl_rodents',
              'uniprot_trembl_human',
              'uniprot_sprot_human',
              'uniprot_sprot_rodents',
              'uniprot_sprot_mammals',
              'uniprot_sprot_vertebrates',
#              'uniprot_trembl_unclassified', # this is the last one that we don't do for now
            ],
            column_names => ['input_file'],
          },
          -analysis_capacity  => 1,
          -flow_into => {
            '2->A' => ['download_by_taxonomy'],
            'A->1' => ['concat_by_taxonomy'],
          },
        },

        {
          -logic_name => 'download_by_taxonomy',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'wget -O #uniprot_dir#/#input_file#.dat.gz -q "#base_uniprot_ftp#/taxonomic_divisions/#input_file#.dat.gz"',
          },
          -analysis_capacity  => 5,
          -flow_into => {
            1 => [ 'gunzip_taxonomy_file' ],
          },
        },

        {
          -logic_name => 'gunzip_taxonomy_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'gunzip #uniprot_dir#/#input_file#.dat.gz',
          },
          -analysis_capacity  => 5,
          -flow_into => {
            1 => [ 'embl2fasta' ],
          },
        },

        {
          -logic_name => 'embl2fasta',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'perl #embl2fasta_script# -embl #uniprot_dir#/#input_file#.dat -format swiss'.
              ' -kill_list_host '.$self->o('kill_list_host').
              ' -kill_list_port '.$self->o('kill_list_port').
              ' -kill_list_user '.$self->o('kill_list_user').
              ' -kill_list_pass '.$self->o('kill_list_pass').
              ' -kill_list_dbname '.$self->o('kill_list_dbname'),
          },
          -analysis_capacity  => 5,
          -flow_into => {
            1 => { 'check_embl2fasta' => {input_file => '#input_file#', filename => '#input_file##fasta_suffix#'}},
          },
        },

        {
          -logic_name => 'check_embl2fasta',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'EXIT_CODE=1; if [ "`grep -c \> #uniprot_dir#/#filename#`" -eq "`grep -c ^SQ #uniprot_dir#/#input_file#.dat`" ]; then rm #uniprot_dir#/#input_file#.dat; EXIT_CODE=0; fi; exit $EXIT_CODE;',
          },
          -analysis_capacity  => 5,
          -flow_into => {
            1 => [ 'pre_process_isoforms', '?accu_name=filename&accu_input_variable=filename&accu_address=[]'],
          },
        },

        {
          -logic_name => 'pre_process_isoforms',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'EXIT_CODE=1; if [ "`echo #input_file# | sed \'s/sprot//\'`" = "#input_file#" ]; then EXIT_CODE=2; else EXIT_CODE=0; fi; exit $EXIT_CODE;',
            return_codes_2_branches => { '2' => 3},
          },
          -analysis_capacity  => 5,
          -flow_into => {
            1 => [ 'process_isoforms' ],
          },
        },

        {
          -logic_name => 'concat_by_taxonomy',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'cd #uniprot_dir#; cat #expr(join(" ", @{#filename#}))expr# > #uniprot_vert_file#; COUNT=0; for F in #expr(join(" ", @{#filename#}))expr#;do FCOUNT=`grep -c \> $F`;COUNT=$((COUNT+FCOUNT));done; if [ "`grep -c \> #uniprot_vert_file#`" -eq "$COUNT" ]; then rm -f #expr(join(" ", @{#filename#}))expr#; else exit 1;fi',
          },
          -analysis_capacity  => 1,
          -flow_into => {
            1 => ['prepare_format_uniprot_and_isoforms'],
          },
        },

        {
          -logic_name => 'prepare_format_uniprot_and_isoforms',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -analysis_capacity  => 1,
          -parameters => {
            inputlist  => [
              ['uniprot_vertebrate_all_isoforms', '#uniprot_vert_file# *.varsplic'],
              ['uniprot_vertebrate', '#uniprot_vert_file#'],
            ],
            column_names => ['db_name', 'input_file'],
          },
          -flow_into => {
            2 => ['concat_uniprot_and_isoforms'],
          }
        },

        {
          -logic_name => 'concat_uniprot_and_isoforms',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'cd #uniprot_dir#; cat #input_file# > #db_name#',
          },
          -analysis_capacity  => 2,
          -flow_into => {
            1 => ['indicate'],
            '1->A' => WHEN( '#blast_type# eq "ncbi"' => 'ncbi_format_uniprot_and_isoforms', ELSE 'wu_format_uniprot_and_isoforms',),
            'A->1' => ['clean_uniprot_and_isoforms'],
          },
        },

        {
          -logic_name => 'clean_uniprot_and_isoforms',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'cd #uniprot_dir#; for F in #input_file# #varsplic_file#; do rm -f "$F"*;done',
          },
          -analysis_capacity  => 1,
        },

        {
          -logic_name => 'indicate',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd   => $self->o('indicate').' -d #uniprot_dir# -f #db_name# --index #uniprot_dir#/#uniprot_file#_index  --parser singleWordParser',
          },
          -analysis_capacity  => 1,
        },

        {
          -logic_name => 'wu_format_uniprot_and_isoforms',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'cd #uniprot_dir#; XDBNAME="`ls #input_file#`"; xdformat -p -o #uniprot_dir#/#xdb_name# -t "`echo $XDBNAME | sed \'s/ /,/g\'`" -v "#uniprot_version#" -d "#uniprot_date#" #input_file#',
          },
          -analysis_capacity  => 2,
        },

        {
          -logic_name => 'ncbi_format_uniprot_and_isoforms',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'makeblastdb -dbtype prot -title "#db_name#" -in #uniprot_dir#/#db_name#',
          },
          -analysis_capacity  => 2,
        },

        {
          -logic_name => 'process_isoforms',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd   => 'perl #process_isoforms_script# -version -infile #uniprot_dir#/#input_file##fasta_suffix# -isofile #uniprot_dir#/#varsplic_file##fasta_suffix# -outfile #uniprot_dir#/#input_file#.varsplic',
          },
          -analysis_capacity  => 5,
        },


#################
# Based on PE levels
#################
        {
          -logic_name  => 'create_pe_level_ids',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters  => {
            inputlist => [
              ['uniprot_PE12_nonvert', '(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+NOT+taxonomy%3aVertebrata+AND+fragment%3ano&compress=yes&format=fasta'],
              ['uniprot_PE12_vertebrata', '(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+AND+taxonomy%3aVertebrata+AND+fragment:no&compress=yes&format=fasta'],
              ['uniprot_PE12_nonvert_frag', '(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+NOT+taxonomy%3aVertebrata+AND+fragment%3ayes&compress=yes&format=fasta'],
              ['uniprot_PE12_vertebrata_frag', '(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+AND+taxonomy%3aVertebrata+AND+fragment:yes&compress=yes&format=fasta'],
            ],
            column_names => ['input_file', 'uniprot_query'],
          },
          -analysis_capacity  => 1,
          -flow_into => {
            '2->A' => ['download_by_pe_level'],
            'A->1' => ['prepare_format_pe_level'],
          },
        },

        {
          -logic_name => 'download_by_pe_level',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'wget -O #uniprot_dir#/#input_file#_initial.gz -q "#base_uniprot_url##uniprot_query#"',
          },
          -analysis_capacity  => 2,
          -flow_into => {
            1 => [ 'gunzip_taxonomy_pe' ],
          },
        },

        {
          -logic_name => 'gunzip_taxonomy_pe',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd         => 'gunzip #uniprot_dir#/#input_file#_initial.gz',
          },
          -analysis_capacity  => 2,
          -flow_into => ['convert_taxonomy_pe'],
        },

        {
          -logic_name => 'convert_taxonomy_pe',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'sed -r -e \'s/^>...([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*PE=([0-9]).*SV=([0-9]{1,}).*/>\1.\4/\' -e \'/^[^>]/ s/O/K/g\' #uniprot_dir#/#input_file#_initial > #uniprot_dir#/#input_file#',
          },
          -analysis_capacity  => 2,
        },

        {
          -logic_name => 'prepare_format_pe_level',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -analysis_capacity  => 1,
          -parameters => {
            inputlist  => [
              ['PE12', 'uniprot_PE12_vertebrata uniprot_PE12_nonvert'],
              ['PE12_wfrag', 'uniprot_PE12_vertebrata uniprot_PE12_nonvert uniprot_PE12_vertebrata_frag uniprot_PE12_nonvert_frag'],
              ['PE12_vertebrata', 'uniprot_PE12_vertebrata'],
              ['PE12_vertebrata_wfrag', 'uniprot_PE12_vertebrata uniprot_PE12_vertebrata_frag'],
            ],
            column_names => ['db_name', 'input_file'],
          },
          -flow_into => {
            2 => ['concat_pe_level_fasta'],
          }
        },

        {
          -logic_name => 'concat_pe_level_fasta',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd         => 'cd #uniprot_dir#; cat #input_file# > #db_name#',
          },
          -analysis_capacity  => 2,
          -flow_into => {
            '1->A' => WHEN( '#blast_type# eq "ncbi"' => 'ncbi_format_pe_level', ELSE 'wu_format_pe_level',),
            'A->1' => ['index_pe12_vertebrata'],
          }
        },

        {
          -logic_name => 'ncbi_format_pe_level',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'makeblastdb -dbtype prot -title "#db_name#" -in #uniprot_dir#/#db_name#',
          },
          -analysis_capacity  => 2,
        },

        {
          -logic_name => 'wu_format_pe_level',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'cd #uniprot_dir#; xdformat -p -o #uniprot_dir#/#xdb_name# -t `echo "#input_file#" | sed "s/ /,/g"` -v "#uniprot_version#" -d "#uniprot_date#" #input_file#',
          },
          -analysis_capacity  => 2,
        },

        {
          -logic_name => 'index_pe12_vertebrata',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => $self->o('indicate').' -d #uniprot_dir# -f PE12_vertebrata --index #uniprot_dir#/PE12_vertebrata_index  --parser singleWordParser',
          },
          -analysis_capacity => 1,
          -flow_into => {
            1 => ['clean_pe_level'],
          },
        },

        {
          -logic_name => 'clean_pe_level',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'cd #uniprot_dir#; for F in #input_file#; do rm -f "$F"*;done',
          },
          -analysis_capacity  => 1,
        },

        {
          -logic_name => 'check_indicate_index',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd   => q(perl -e 'use strict;use warnings;my %acc;open(FH, "#uniprot_dir#/uniprot_vertebrate") || die("Could not open #uniprot_dir#/uniprot_vertebrate");while(my $line = <FH>) {if ($line =~ />(\S+)/) {$acc{$1} = 1}};close(FH) || die("Could not close #uniprot_dir#/uniprot_vertebrate");open(RH, "#uniprot_dir#/PE12_vertebrata") || die("Could not open #uniprot_dir#/PE12_vertebrata");while (my $enil = <RH>) {if ($enil =~ />(\S+)/) {print STDOUT "$1\n" unless (exists $acc{$1});}};close(RH) || die ("Could not close #uniprot_dir#/PE12_vertebrata")'),
          },
          -analysis_capacity  => 1,
        },


        {
          -logic_name => 'finalise_dbs',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'cd #uniprot_dir#; find ./ -type d -execdir chmod g+w {} \;',
          },
          -analysis_capacity  => 1,
          -flow_into => {
            1 => ['propagate_the_word'],
          },
        },

        {
          -logic_name => 'propagate_the_word',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SlackNotification',
          -parameters => {
            slack_webhook => $self->o('slack_hook_url'),
            text => 'The UniProtKB blast database has been updated to '.$self->o('uniprot_version'),
          },
        },
    );

    foreach my $analysis (@analyses) {
        $analysis->{'-max_retry_count'} = 0 unless (exists $analysis->{'-max_retry_count'});
        $analysis->{'-meadow_type'} = 'LSF' unless (exists $analysis->{'-meadow_type'});
    }
    return \@analyses;
}

1;
