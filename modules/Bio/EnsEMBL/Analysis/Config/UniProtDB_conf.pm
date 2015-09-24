=pod

=head1 NAME

    Bio::EnsEMBL::Hive::PipeConfig::UniProtDB_conf

=head1 DESCRIPTION

    Configuration to create multiple blast databases which can be shared by Enesmbl and Havana
    It has to be used in conjonction with the create_uniprot_database.sh in ensembl-analysis/scripts

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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


package UniProtDB_conf;

use strict;
use warnings;

use Bio::EnsEMBL::ApiVersion ();

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');


=head2 default_options

    Description : Interface method that should return a hash of option_name->default_option_value pairs.
                  Please see existing PipeConfig modules for examples.

=cut

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },

        'dbowner' => $self->o('ENV', 'USER'),
        'pipeline_name' => $self->o('ENV', 'pipeline_name'),
        'base_uniprot_ftp'  => $ENV{'BASE_UNIPROT_FTP'} || 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase',
        'base_uniprot_url'  => $ENV{'BASE_UNIPROT_URL'} || 'http://www.uniprot.org/uniprot/?query=',
        'uniprot_vert_file'  => $ENV{'UNIPROT_VERT_FILE'} || $self->o('uniprot_vert_file'),
        'uniprot_file'  => $ENV{'UNIPROT_FILE'} || 'uniprot',
        'uniprot_dir'  => $ENV{'UNIPROT_DIR'},
        'varsplic_file'  => $ENV{'VARSPLIC_FILE'} || 'uniprot_sprot_varsplic',
        'embl2fasta_script'  => $ENV{'EMBL2FASTA_SCRIPT'},
        'process_isoforms_script'  => $ENV{'PROCESS_ISOFORMS_SCRIPT'},
        'fasta_suffix'  => $ENV{'FASTA_SUFFIX'} || '.fasta',
        'uniprot_version'  => $ENV{'UNIPROT_VERSION'},
        'uniprot_date'  => $ENV{'UNIPROT_DATE'} || localtime,
    };
}


sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        'base_uniprot_ftp'  => $self->o('base_uniprot_ftp'),
        'base_uniprot_url'  => $self->o('base_uniprot_url'),
        'uniprot_vert_file'  => $self->o('uniprot_vert_file'),
        'uniprot_file'  => $self->o('uniprot_file'),
        'embl2fasta_script'  => $self->o('embl2fasta_script'),
        'uniprot_dir'  => $self->o('uniprot_dir'),
        'varsplic_file'  => $self->o('varsplic_file'),
        'fasta_suffix'  => $self->o('fasta_suffix'),
        'uniprot_version'  => $self->o('uniprot_version'),
        'uniprot_date'  => $self->o('uniprot_date'),
        'process_isoforms_script'  => $self->o('process_isoforms_script'),
    };
}


sub resource_classes {
    my $self = shift @_;

    my $resources = $self->SUPER::resource_classes();
    $resources->{'1GB_mem'}->{LSF} = '-M1000 -R"select[mem>1000] rusage[mem=1000]"';
    return $resources;
}


sub pipeline_analyses {
    my ($self) = @_;
    my @analyses = (

        {   -logic_name => 'setup_directory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
            -parameters => {
                'cmd'     => 'mkdir #uniprot_dir#',
            },

            -analysis_capacity  => 1,

            -input_ids  => [ { dummy => 'dummy' } ],

        },
        {   -logic_name => 'download_by_taxonomy',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'     => 'wget -O #uniprot_dir#/#input_file#.dat.gz -q "#base_uniprot_ftp#/taxonomic_divisions/#input_file#.dat.gz"',
            },

            -analysis_capacity  => 5,
            -wait_for => ['setup_directory'],

            -input_ids  => [
            {
                input_file => 'uniprot_sprot_human'
            },
            {
                input_file => 'uniprot_sprot_mammals'
            },
            {
                input_file => 'uniprot_sprot_rodents'
            },
            {
                input_file => 'uniprot_sprot_vertebrates'
            },
            {
                input_file => 'uniprot_trembl_human'
            },
            {
                input_file => 'uniprot_trembl_mammals'
            },
            {
                input_file => 'uniprot_trembl_rodents'
            },
            {
                input_file => 'uniprot_trembl_vertebrates'
            },
# The unclassified group will not be used but I'm leaving it here just in case
#            {
#                input_file => 'uniprot_trembl_unclassified'
#            },
          ],

            -flow_into => {
                1 => [ 'gunzip_taxonomy_file' ],
            },
        },

        {   -logic_name => 'gunzip_taxonomy_file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'gunzip #uniprot_dir#/#input_file#.dat.gz',
            },
            -analysis_capacity  => 5,
            -flow_into => {
                1 => [ 'embl2fasta' ],
            },
        },

        {   -logic_name => 'embl2fasta',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'   => 'perl #embl2fasta_script# #uniprot_dir#/#input_file#.dat',
            },
            -analysis_capacity  => 5,
            -flow_into => {
                1 => [ 'check_embl2fasta' ],
            },
        },

        {   -logic_name => 'check_embl2fasta',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'   => 'EXIT_CODE=1; if [ "`grep -c \> #uniprot_dir#/#input_file##fasta_suffix#`" -eq "`grep -c ^SQ #uniprot_dir#/#input_file#.dat`" ]; then rm #uniprot_dir#/#input_file#.dat; EXIT_CODE=0; fi; exit $EXIT_CODE;',
            },
            -analysis_capacity  => 5,
            -flow_into => {
                1 => [ 'concat_by_taxonomy', 'pre_process_isoforms' ],
            },
        },

        {   -logic_name => 'pre_process_isoforms',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type => 'LOCAL',
            -parameters => {
                'cmd'   => 'EXIT_CODE=1; if [ "`echo #input_file# | sed \'s/sprot//\'`" = "#input_file#" ]; then EXIT_CODE=2; else EXIT_CODE=0; fi; exit $EXIT_CODE;',
                'return_codes_2_branches' => { '2' => 3},
            },
            -analysis_capacity  => 5,
            -flow_into => {
                1 => [ 'process_isoforms' ],
            },
        },

        {   -logic_name => 'concat_by_taxonomy',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'   => 'cat #uniprot_dir#/#input_file##fasta_suffix# > #uniprot_dir#/#uniprot_vert_file#',
            },
            -analysis_capacity  => 1,
        },

        {   -logic_name => 'indicate',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'   => '/software/ensembl/bin/indicate -d #uniprot_dir# -f #uniprot_vert_file# --index #uniprot_dir#/#uniprot_file#_index  --parser singleWordParser',
            },
            -rc_name => '1GB_mem',
            -analysis_capacity  => 1,
            -wait_for => ['gunzip_taxonomy_file', 'embl2fasta', 'check_embl2fasta', 'concat_by_taxonomy'],
            -input_ids => [ { dummy => 'dummy' } ],
        },

        {   -logic_name => 'xdformat_uniprot_and_isoforms',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'cd #uniprot_dir#; XDBNAME="`ls #input_file#`"; xdformat -p -o #uniprot_dir#/#xdb_name# -t "`echo $XDBNAME | sed \'s/ /,/g\'`" -v "#uniprot_version#" -d "#uniprot_date#" #input_file#',
            },
            -rc_name => '1GB_mem',
            -input_ids  => [
            {
                xdb_name => 'uniprot_vertebrate_all_isoforms',
                input_file => '#uniprot_vert_file# *.varsplic',
            },
            ],
            -wait_for => [ 'gunzip_taxonomy_file', 'embl2fasta', 'check_embl2fasta', 'process_isoforms', 'concat_by_taxonomy'],
            -analysis_capacity  => 5,
            -flow_into => {
                1 => ['check_xdformat_uniprot_and_isoforms'],
            },
        },

        {   -logic_name => 'check_xdformat_uniprot_and_isoforms',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'EXIT_CODE=1; cd #uniprot_dir#; if [ "`grep -c \> #input_file# | cut -d \':\' -f2 | awk \'{SUM += $1} END {print SUM}\'`" -eq "`xdformat -p -i #uniprot_dir#/#xdb_name# 2>&1 | grep letters | awk \'{print $5}\' | sed \'s/,//g\'`" ]; then EXIT_CODE=0;fi; exit $EXIT_CODE',
            },
            -analysis_capacity  => 1,
        },

        {   -logic_name => 'download_isoforms',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'     => 'wget -O #uniprot_dir#/#input_file##fasta_suffix#.gz -q "#base_uniprot_ftp#/complete/#input_file##fasta_suffix#.gz"',
            },

            -analysis_capacity  => 1,
            -wait_for => ['setup_directory'],

            -input_ids  => [
            {
                input_file => '#varsplic_file#',
                check_file => 'http://web.expasy.org/docs/relnotes/relstat.html',
            },
          ],

            -flow_into => {
                1 => [ 'gunzip_isofasta_file' ],
            },
        },

        {   -logic_name => 'gunzip_isofasta_file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'gunzip -c #uniprot_dir#/#input_file##fasta_suffix#.gz | perl -ne \'if (!/>/) {s/O/K/g} print;\' > #uniprot_dir#/#input_file##fasta_suffix#',
            },
            -analysis_capacity  => 1,
            -flow_into => {
                1 => [ 'check_isofasta_file' ],
            },
        },

        {   -logic_name => 'check_isofasta_file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'EXIT_CODE=1; if [ "`wget -q -O - #check_file# | perl -ne \'if(/alternative\s+splicing.*:\s+(\d+)/) {print $1}\'`" -eq "`grep -c \> #uniprot_dir#/#input_file##fasta_suffix#`" ]; then rm #uniprot_dir#/#input_file##fasta_suffix#.gz; EXIT_CODE=0; fi; exit $EXIT_CODE',
            },
            -analysis_capacity  => 1,
        },

        {   -logic_name => 'process_isoforms',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'   => 'perl #process_isoforms_script# -version -infile #uniprot_dir#/#input_file##fasta_suffix# -isofile #uniprot_dir#/#varsplic_file##fasta_suffix# -outfile #uniprot_dir#/#input_file#.varsplic',
            },
            -wait_for => [ 'check_isofasta_file'],
            -analysis_capacity  => 5,
        },

        {   -logic_name => 'download_fasta',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'     => 'wget -O #uniprot_dir#/#input_file##fasta_suffix#.gz -q "#base_uniprot_ftp#/complete/#input_file##fasta_suffix#.gz"',
            },

            -analysis_capacity  => 2,
            -wait_for => ['setup_directory'],

            -input_ids  => [
            {
                input_file => 'uniprot_sprot',
                check_file => 'http://web.expasy.org/docs/relnotes/relstat.html',
            },
            {
                input_file => 'uniprot_trembl',
                check_file => 'http://www.ebi.ac.uk/uniprot/TrEMBLstats/',
            },
          ],

            -flow_into => {
                1 => [ 'gunzip_fasta_file' ],
            },
        },

        {   -logic_name => 'gunzip_fasta_file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'gunzip -c #uniprot_dir#/#input_file##fasta_suffix#.gz | perl -ne \'if (/>/) { s/>...([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*PE=([0-9]).*SV=([0-9]).*/>\$1.\$4 \$3/} else {s/O/K/g} print;\' > #uniprot_dir#/#input_file##fasta_suffix#',
            },
            -analysis_capacity  => 2,
            -flow_into => {
                1 => [ 'check_fasta_file' ],
            },
        },

        {   -logic_name => 'check_fasta_file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'EXIT_CODE=1; if [ "`wget -q -O - #check_file# | perl -ne \'if(/(\d+)\s+sequence\s+entries/) {print $1}\'`" -eq "`grep -c \> #uniprot_dir#/#input_file##fasta_suffix#`" ]; then rm #uniprot_dir#/#input_file##fasta_suffix#.gz;EXIT_CODE=0; fi; exit $EXIT_CODE',
            },
            -analysis_capacity  => 2,
            -flow_into => {
                1 => [ 'concat_fasta', 'pre_entry_loc' ],
            },
        },

        {   -logic_name => 'pre_entry_loc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type => 'LOCAL',
            -parameters => {
                'cmd'   => 'EXIT_CODE=1; if [ "`echo #input_file# | sed \'s/sprot//\'`" = "#input_file#" ]; then EXIT_CODE=2; else EXIT_CODE=0; fi; exit $EXIT_CODE;',
                'return_codes_2_branches' => { '2' => 3},
            },
            -analysis_capacity  => 1,
            -flow_into => {
                1 => [ 'entry_loc' ],
            },
        },

        {   -logic_name => 'entry_loc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'grep \> #uniprot_dir#/#input_file##fasta_suffix# | awk "{print \$1, "STD"}" > #uniprot_dir#/entry_loc',
            },
            -analysis_capacity  => 1,
            -wait_for => ['download_fasta', 'gunzip_fasta_file', 'check_fasta_file'],
        },

        {   -logic_name => 'concat_fasta',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'   => 'grep -H -c \> #uniprot_dir#/#input_file##fasta_suffix# | sed \'s/:/ /\' >> #uniprot_dir#/#uniprot_file#.a; cat #uniprot_dir#/#input_file##fasta_suffix# >> #uniprot_dir#/#uniprot_file#',
            },
            -analysis_capacity  => 1,
        },

        {   -logic_name => 'check_concat_fasta',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'   => 'EXIT_CODE=1; if [ "`awk \'{SUM+=$2} END {print SUM}\' #uniprot_dir#/#uniprot_file#.a`" -eq "`grep -c \> #uniprot_dir#/#uniprot_file#`" ]; then ECOUNT=`grep -c \> #uniprot_dir#/entry_loc`; for C in `awk \'{print $2}\' #uniprot_dir#/#uniprot_file#.a`; do if [ "$ECOUNT" -eq "$C" ]; then awk \'{print $1}\' #uniprot_dir#/#uniprot_file#.a | xargs -I {} rm {}; rm #uniprot_dir#/#uniprot_file#.a; EXIT_CODE=0; fi; done; fi; exit $EXIT_CODE',
            },
            -analysis_capacity  => 1,
            -wait_for => ['concat_fasta', 'entry_loc'],
            -input_ids => [ { dummy => 'dummy' } ],
        },

        {   -logic_name => 'xdformat_blast_uniprot',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'   => 'xdformat -p -t "#uniprot_file#" -v "#uniprot_version#" -d "#uniprot_date#" #uniprot_dir#/#uniprot_file#',
            },
            -rc_name => '1GB_mem',
            -analysis_capacity  => 1,
            -wait_for => ['check_concat_fasta'],
            -input_ids => [ { dummy => 'dummy' } ],
        },

        {   -logic_name => 'check_xdformat_uniprot',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'EXIT_CODE=1; cd #uniprot_dir#; if [ "`grep -c \> #uniprot_file#`" -eq "`xdformat -p -i #uniprot_dir#/#uniprot_file# 2>&1 | grep letters | awk \'{print $5}\' | sed \'s/,//g\'`" ]; then EXIT_CODE=0;fi; exit $EXIT_CODE',
            },
            -analysis_capacity  => 1,
            -wait_for => ['xdformat_blast_uniprot'],
            -input_ids => [ { dummy => 'dummy' } ],
        },

        {   -logic_name => 'download_by_pe_level',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'     => 'wget -O #uniprot_dir#/#input_file#_initial.gz -q "#base_uniprot_url##uniprot_query#"',
            },

            -analysis_capacity  => 2,
            -wait_for => ['setup_directory'],

            -input_ids  => [
            {
                input_file => 'uniprot_PE12_vertebrata',
                uniprot_query => '(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+AND+taxonomy%3aCraniata+AND+fragment:no&compress=yes&format=fasta',
            },
            {
                input_file => 'uniprot_PE12_nonvert',
                uniprot_query => '(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+NOT+taxonomy%3aCraniata+AND+fragment%3ano&compress=yes&format=fasta',
            },
            {
                input_file => 'uniprot_PE12_vertebrata_frag',
                uniprot_query => '(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+AND+taxonomy%3aCraniata+AND+fragment:yes&compress=yes&format=fasta',
            },
            {
                input_file => 'uniprot_PE12_nonvert_frag',
                uniprot_query => '(existence%3a%22evidence+at+protein+level%22+OR+existence%3a%22evidence+at+transcript+level%22)+NOT+taxonomy%3aCraniata+AND+fragment%3ayes&compress=yes&format=fasta',
            },
          ],

            -flow_into => {
                1 => [ 'gunzip_taxonomy_pe' ],
            },
        },

        {   -logic_name => 'gunzip_taxonomy_pe',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'gunzip #uniprot_dir#/#input_file#_initial.gz',
            },
            -analysis_capacity  => 2,
            -flow_into => ['convert_taxonomy_pe'],
        },

        {   -logic_name => 'convert_taxonomy_pe',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'sed -r -e \'s/^>...([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*PE=([0-9]).*SV=([0-9]).*/>\1.\4/\' -e \'/^[^>]/ s/O/K/g\' #uniprot_dir#/#input_file#_initial > #uniprot_dir#/#input_file#',
            },
            -analysis_capacity  => 2,
        },

        {   -logic_name => 'xdformat_pe_level',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'cd #uniprot_dir#; xdformat -p -o #uniprot_dir#/#xdb_name# -t `echo "#input_file#" | sed "s/ /,/g"` -v "#uniprot_version#" -d "#uniprot_date#" #input_file#',
            },
            -rc_name => '1GB_mem',
            -wait_for => ['download_by_pe_level', 'gunzip_taxonomy_pe', 'convert_taxonomy_pe'],
            -flow_into => ['check_xdformat_pe_level'],
            -analysis_capacity  => 2,
            -input_ids  => [
            {
                xdb_name => 'PE12',
                input_file => 'uniprot_PE12_vertebrata uniprot_PE12_nonvert',
            },
            {
                xdb_name => 'PE12_wfrag',
                input_file => 'uniprot_PE12_vertebrata uniprot_PE12_nonvert uniprot_PE12_vertebrata_frag uniprot_PE12_nonvert_frag',
            },
            {
                xdb_name => 'PE12_vertebrata',
                input_file => 'uniprot_PE12_vertebrata_frag',
            },
            {
                xdb_name => 'PE12_vertebrata_wfrag',
                input_file => 'uniprot_PE12_vertebrata uniprot_PE12_vertebrata_frag',
            },
            ],
        },

        {   -logic_name => 'check_xdformat_pe_level',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd'         => 'EXIT_CODE=1; cd #uniprot_dir#; CFILE=`grep -c \> #input_file# | cut -d \':\' -f2 | awk \'{SUM += $1} END {print SUM}\'`;if [ "$CFILE" -eq "`xdformat -p -i #uniprot_dir#/#xdb_name# 2>&1 | grep letters | awk \'{print $5}\' | sed \'s/,//g\'`" ]; then EXIT_CODE=0;fi; exit $EXIT_CODE',
            },
            -analysis_capacity  => 2,
        },

        {   -logic_name => 'finalise_dbs',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type => 'LOCAL',
            -parameters => {
                'cmd'         => 'cd #uniprot_dir#; find ./ -type d -execdir chmod g+w {} \;',
            },
            -analysis_capacity  => 1,
            -input_ids => [ { dummy => 'dummy' } ],
            -wait_for => ['check_xdformat_uniprot_and_isoforms', 'check_xdformat_uniprot', 'check_xdformat_pe_level' ],
        },

    );

    foreach my $analysis (@analyses) {
        $analysis->{'-max_retry_count'} = 1 unless (exists $analysis->{'-max_retry_count'});
        $analysis->{'-meadow_type'} = 'LSF' unless (exists $analysis->{'-meadow_type'});
    }
    return \@analyses;
}

1;
