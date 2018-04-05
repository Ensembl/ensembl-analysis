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

=cut

package HiveRNASeq_conf;

use strict;
use warnings;

use File::Spec::Functions qw(catfile catdir file_name_is_absolute);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw (get_analysis_settings);
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },

##########################################################################
#                                                                        #
# CHANGE STUFF HERE                                                      #
#                                                                        #
##########################################################################
        'pipeline_name' => '',

        'user_r'   => '',
        'user'     => '',
        'password' => '',
        'port'     => '',
        'species'  => '', # It MUST be the scientific name so the analyses are created correctly
        'assembly_name' => '',
        'email' => '', # Add your email so you can be notified when a bam file is removed

        'pipe_db_name'   => $self->o('dbowner').'_'.$self->o('pipeline_name').'_hive',
        'dna_db_name'    => '',
        'blast_db_name'  => $self->o('dbowner').'_'.$self->o('pipeline_name').'_'.$self->o('species').'_blast',
        'refine_db_name' => $self->o('dbowner').'_'.$self->o('pipeline_name').'_'.$self->o('species').'_refine',
        'rough_db_name'  => $self->o('dbowner').'_'.$self->o('pipeline_name').'_'.$self->o('species').'_rough',

        'pipe_db_server'   => '',
        'dna_db_server'    => '',
        'dna_db_port'   => $self->o('port'),
        'data_db_server'  => '', # Server for the blast, refine and rough DBs, or you can set them down below
        'data_db_port'  => $self->o('port'), # Port for the blast, refine and rough DBs, or you can set them down below
        'blast_db_server'  => $self->o('data_db_server'),
        'refine_db_server' => $self->o('data_db_server'),
        'rough_db_server'  => $self->o('data_db_server'),

        'input_dir'  => '', #You need to specify the full path to the directory where your input files like your fastq files are
        'output_dir' => '', #You need to specify the full path to the directory where your output files will be written
        'merge_dir'  => '', #You can specify a full path to the directory where your BAM files will be, default is <output_dir>
        'sam_dir'    => '', #You can specify a full path to the directory where your SAM files will be, default is <output_dir>/SAM
        'rnaseq_ftp_base' => 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/',

        'genome_file'     => 'genome/genome.fa',
        'use_ucsc_naming' => 0,

        'sequence_dump_script' => catfile($self->o('enscode_root_dir'),'ensembl-analysis', 'scripts', 'sequence_dump.pl'),
        'create_type' => 'clone',

# You have the choice between:
#  * using a csv file you already created
#  * using a study_accession like PRJEB19386
#  * using the taxon_id of your species
# 'rnaseq_summary_file' should always be set. If 'taxon_id' or 'study_accession' are not undef
# they will be used to retrieve the information from ENA and to create the csv file. In this case,
# 'file_columns' and 'summary_file_delimiter' should not be changed unless you know what you are doing
        'rnaseq_summary_file' => catfile($self->o('output_dir'), $self->o('pipeline_name').'.csv'), # You need to specify the full path to your csv summary file
        'taxon_id' => undef,
        'study_accession' => undef,
        # Use this option to change the delimiter for your summary data
        # file.
        summary_file_delimiter => '\t',
        summary_csv_table => 'csv_data',
        rnaseq_data_provider => 'ENA', #It will be set during the pipeline or it will use this value


        'samtools' => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary
        'picard_lib_jar' => catfile($self->o('software_base_path'), 'Cellar', 'picard-tools', '2.6.0', 'libexec', 'picard.jar'), #You need to specify the full path to the picard library
        'short_read_aligner' => catfile($self->o('software_base_path'), 'opt', 'bwa-051mt', 'bin', 'bwa'), #You may need to specify the full path to the bwa binary
        'refine_ccode_exe' => catfile($self->o('binary_base'), 'RefineSolexaGenes'), #You may need to specify the full path to the RefineSolexaGenes binary

        # Blast database for comparing the final models to.
        uniprotdb => '',

        # Indicate Index for the blast database.
        uniprotindex => '',

        blastp => catfile($self->o('binary_base'), 'blastp'), #You may need to specify the full path to the blastp binary
        # blast used, it can be either ncbi or wu, it is overriding the -type value from BLAST_PARAMS
        blast_type => 'ncbi',

        splicing_aligner => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'), #You may need to specify the full path to the exonerate binary version 0.9.0

        # If your reads are unpaired you may want to run on slices to avoid
        # making overlong rough models.  If you want to do this, specify a
        # slice length here otherwise it will default to whole chromosomes.
        slice_length => 10000000,

        # Regular expression to allow FastQ files to be correctly paired,
        # for example: file_1.fastq and file_2.fastq could be paired using
        # the expression "\S+_(\d)\.\S+".  Need to identify the read number
        # in brackets; the name the read number (1, 2) and the
        # extension.
        pairing_regex => '\S+_(\d)\.\S+',
        paired => 1,

        # Do you want to make models for the each individual sample as well
        # as for the pooled samples (1/0)?
        single_tissue => 1,

        # What Read group tag would you like to group your samples
        # by? Default = ID
        read_group_tag => 'SM',
        read_id_tag => 'ID',

        use_threads => 3,
        read_min_paired => 50,
        read_min_mapped => 50,
        other_isoforms => 'other', # If you don't want isoforms, set this to undef
        normal_queue => 'production-rh7', # If LSF/other system has submission queues with multiple run time allowed, you might want to change this
        long_queue => 'production-rh7', # If LSF/other system has submission queues with multiple run time allowed, you might want to change this

        # Please assign some or all columns from the summary file to the
        # some or all of the following categories.  Multiple values can be
        # separted with commas. ID, SM, DS, CN, is_paired, filename, read_length, is_13plus,
        # is_mate_1 are required. If pairing_regex can work for you, set is_mate_1 to -1.
        # You can use any other tag specified in the SAM specification:
        # http://samtools.github.io/hts-specs/SAMv1.pdf

        ####################################################################
        # This is just an example based on the file snippet shown below.  It
        # will vary depending on how your data looks.
        ####################################################################

        file_columns => ['SM', 'ID', 'is_paired', 'filename', 'is_mate_1', 'read_length', 'is_13plus', 'CN', 'PL', 'DS'],


##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

        maxintron => 200000,

        blast_db_port => $self->o('data_db_port'),
        blast_db_user => $self->o('user'),
        blast_db_password => $self->o('password'),
        blast_db_driver => $self->o('hive_driver'),
        refine_db_port => $self->o('data_db_port'),
        refine_db_user => $self->o('user'),
        refine_db_password => $self->o('password'),
        refine_db_driver => $self->o('hive_driver'),
        rough_db_port => $self->o('data_db_port'),
        rough_db_user => $self->o('user'),
        rough_db_password => $self->o('password'),
        rough_db_driver => $self->o('hive_driver'),

        databases_to_delete => ['blast_db', 'refine_db', 'rough_db'],

        'blast_db' => {
                           -dbname => $self->o('blast_db_name'),
                           -host   => $self->o('blast_db_server'),
                           -port   => $self->o('blast_db_port'),
                           -user   => $self->o('blast_db_user'),
                           -pass   => $self->o('blast_db_password'),
                           -driver => $self->o('blast_db_driver'),
                         },

        'refine_db' => {
                           -dbname => $self->o('refine_db_name'),
                           -host   => $self->o('refine_db_server'),
                           -port   => $self->o('refine_db_port'),
                           -user   => $self->o('refine_db_user'),
                           -pass   => $self->o('refine_db_password'),
                           -driver => $self->o('refine_db_driver'),
                        },

        'rough_db' => {
                           -dbname => $self->o('rough_db_name'),
                           -host   => $self->o('rough_db_server'),
                           -port   => $self->o('rough_db_port'),
                           -user   => $self->o('rough_db_user'),
                           -pass   => $self->o('rough_db_password'),
                           -driver => $self->o('rough_db_driver'),
                         },
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    my $output_sam_dir = $self->o('sam_dir') ? $self->o('sam_dir') : catdir($self->o('output_dir'), 'SAM');
    my $merge_dir = $self->o('merge_dir') ? $self->o('merge_dir') : catdir($self->o('output_dir'), 'merge_out');
    my $genome_file = file_name_is_absolute($self->o('genome_file')) ? $self->o('genome_file') : catdir($self->o('input_dir'), $self->o('genome_file'));
    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
                         wide_genome_file => $genome_file,
                         wide_input_dir => $self->o('input_dir'),
                         wide_output_dir => $self->o('output_dir'),
                         wide_merge_dir => $merge_dir,
                         wide_short_read_aligner => $self->o('short_read_aligner'),
                         wide_samtools => $self->o('samtools'),
                         wide_output_sam_dir => $output_sam_dir,
                         wide_species => $self->o('species'),
                         wide_use_ucsc_naming => $self->o('use_ucsc_naming'),
                         wide_intron_bam_file => catfile($self->o('output_dir'), 'introns'),
    };
}


=head2 pipeline_create_commands

 Arg [1]    : None
 Description: Add a table named with 'summary_csv_table' to store information about the reads
              The columns are defined by 'file_columns'
 Returntype : Arrayref String, commands to create/delete databases and/or tables
 Exceptions : None

=cut

sub pipeline_create_commands {
    my ($self) = @_;
    my $tables;
    my %small_columns = (
        paired => 1,
        read_length => 1,
        is_13plus => 1,
        is_mate_1 => 1,
        );
    # We need to store the values of the csv file to easily process it. It will be used at different stages
    foreach my $key (@{$self->default_options->{'file_columns'}}) {
        if (exists $small_columns{$key}) {
            $tables .= $key.' SMALLINT UNSIGNED NOT NULL,'
        }
        elsif ($key eq 'DS') {
            $tables .= $key.' VARCHAR(255) NOT NULL,'
        }
        else {
            $tables .= $key.' VARCHAR(50) NOT NULL,'
        }
    }
    $tables .= ' KEY(SM), KEY(ID)';
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      $self->db_cmd('CREATE TABLE '.$self->o('summary_csv_table')." ($tables)"),
    ];
}

=head2 create_header_line

 Arg [1]    : Arrayref String, it will contains the values of 'file_columns'
 Example    : create_header_line($self->o('file_columns');
 Description: It will create a RG line using only the keys present in your csv file
 Returntype : String representing the RG line in a BAM file
 Exceptions : None


=cut

sub create_header_line {
    my ($items) = shift;

    my @read_tags = qw(ID SM DS CN DT FO KS LB PG PI PL PM PU);
    my $read_line = '@RG';
    foreach my $rt (@read_tags) {
        $read_line .= "\t$rt:#$rt#" if (grep($rt eq $_, @$items));
    }
    return $read_line."\n";
}

## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    my $header_line = create_header_line($self->default_options->{'file_columns'});
    my @analysis = (
 {
      -logic_name => 'checking_file_path',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '1GB',
        -parameters => {
            cmd => 'EXIT_CODE=0; for F in #wide_short_read_aligner# #wide_samtools# '.join (' ', $self->o('splicing_aligner'), $self->o('sequence_dump_script'), $self->o('blastp')).'; do which "$F"; if [ "$?" == 1 ]; then EXIT_CODE=1;fi; done; for D in #wide_output_dir# #wide_input_dir# #wide_merge_dir# #wide_output_sam_dir# `dirname #wide_genome_file#`; do mkdir -p "$D"; done; exit $EXIT_CODE',
        },
        -input_ids => [{
          alignment_bam_file => catfile('#wide_merge_dir#', '#assembly_name#.#rnaseq_data_provider#.merged.1.bam'),
          assembly_name => $self->o('assembly_name'),
          inputfile => $self->o('rnaseq_summary_file'),
          }],
        -flow_into => {
            1 => ['downloading_csv'],
        },
  },
 {
      -logic_name => 'downloading_csv',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
        -rc_name => '1GB',
        -parameters => {
          study_accession => $self->o('study_accession'),
          taxon_id => $self->o('taxon_id'),
        },
        -flow_into => {
            '1->A' => ['create_rnaseq_genome_file', 'create_fastq_download_jobs'],
            'A->1' => ['create_rough_db'],
        },
  },
 {
      -logic_name => 'create_rnaseq_genome_file',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '1GB',
        -parameters => {
            cmd => 'if [ ! -s "#wide_genome_file#" ]; then perl '.$self->o('sequence_dump_script').' -dbhost '.$self->o('dna_db_server').' -dbuser '.$self->o('dna_db_user').' -dbport '.$self->o('dna_db_port').' -dbname '.$self->o('dna_db_name').' -coord_system_name '.$self->o('assembly_name').' -toplevel -onefile -header rnaseq -filename #wide_genome_file#;fi',
        },
        -flow_into => {
            1 => [ 'index_genome_file'],
        },
  },
 {
      -logic_name => 'index_genome_file',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '5GB',
        -parameters => {
            cmd => 'if [ ! -e "#wide_genome_file#.ann" ]; then #wide_short_read_aligner# index -a bwtsw #wide_genome_file#;fi',
        },
  },

  {
   -logic_name => 'create_fastq_download_jobs',
   -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFastqDownloadJobs',
   -parameters => {
     inputfile => $self->o('rnaseq_summary_file'),
   },
   -flow_into => {
     2 => {'download_RNASeq_fastqs' => {'iid' => '#iid#'}},
   },
  },

  {
    -logic_name => 'download_RNASeq_fastqs',
    -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs',
    -parameters =>{
      ftp_base_url => $self->o('rnaseq_ftp_base'),
      input_dir => $self->o('input_dir'),
    },
  },

  {
    -logic_name => 'create_rough_db',
    -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
    -parameters => {
                     source_db => $self->o('dna_db'),
                     target_db => $self->o('rough_db'),
                     create_type => $self->o('create_type'),
                   },
    -rc_name => '1GB',
    -flow_into => {
      1 => [ 'parse_summary_file'],
    },
  },

 {
      -logic_name => 'parse_summary_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable',
        -rc_name => '1GB',
        -parameters => {
            column_names => $self->o('file_columns'),
            sample_column => $self->o('read_group_tag'),
            inputfile => $self->o('rnaseq_summary_file'),
            delimiter => $self->o('summary_file_delimiter'),
            csvfile_table => $self->o('summary_csv_table'),
            pairing_regex => $self->o('pairing_regex'),
        },
        -flow_into => {
            '2->A' => [ 'create_tissue_jobs'],
            'A->1' => [ 'merged_bam_file' ],
        },
  },
            {
        -logic_name => 'create_tissue_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
            inputquery => join(' ', 'SELECT', $self->o('read_group_tag'), ',', $self->o('read_id_tag'), ', is_paired, CN', 'FROM', $self->o('summary_csv_table'), 'WHERE', $self->o('read_group_tag'), '= "#sample_name#"'),
            column_names => [$self->o('read_group_tag'), $self->o('read_id_tag'), 'is_paired', 'rnaseq_data_provider'],
                       },
        -rc_name    => '1GB',
        -flow_into => {
            '2->A' => ['create_bwa_jobs'],
            'A->1' => ['merged_tissue_file'],
            },
      },
            {
        -logic_name => 'create_bwa_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateBWAJobs',
        -parameters => {
            sample_column => $self->o('read_group_tag'),
            sample_id_column => $self->o('read_id_tag'),
            csvfile_table => $self->o('summary_csv_table'),
            column_names => $self->o('file_columns'),
            use_threading => $self->o('use_threads'),
                       },
        -rc_name    => '1GB',
        -flow_into => {
            '2->A' => ['bwa', 'create_header_files'],
            'A->1' => ['bwa2bam'],
            },
      },
  {
      -logic_name => 'create_header_files',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '1GB',
        -parameters => {
            cmd => 'if [ ! -e "#wide_output_dir#/#'.$self->o('read_id_tag').'#_header.h" ]; then printf "'.$header_line.'" > #wide_output_dir#/#'.$self->o('read_id_tag').'#_header.h; fi',
        },
  },
            {
        -logic_name => 'bwa',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -parameters => {
          disconnect_jobs => 1,
        },
        -flow_into => {
            1 => [ ':////accu?fastq=[]' ],
            -1 => [ 'bwa_failed' ],
            -2 => [ 'bwa_failed' ],
            },
        -rc_name    => '5GB_multithread',
      },
            {
        -logic_name => 'bwa_failed',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -can_be_empty => 1,
        -parameters => {
          disconnect_jobs => 1,
        },
        -flow_into => {
            1 => [ ':////accu?fastq=[]' ],
            },
        -rc_name    => '10GB_multithread',
      },
            {
        -logic_name => 'bwa2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -parameters => {
                         sampe_options => '-A -a '.$self->o('maxintron'),
                         samse_options => '',
                         min_paired => $self->o('read_min_paired'),
                         min_mapped => $self->o('read_min_mapped'),
                         header_file => '#wide_output_dir#/#'.$self->o('read_id_tag').'#_header.h',
                         bam_prefix => $self->o('read_id_tag'),
                         email => $self->o('email_address'),
                         disconnect_jobs => 1,
                       },
        -flow_into => {
            1 => [ ':////accu?filename=[]' ],
            -1 => {'bwa2bam_10GB' => { fastq => '#fastq#', is_paired => '#is_paired#', $self->o('read_id_tag') => '#'.$self->o('read_id_tag').'#'}},
            -2 => {'bwa2bam_10GB' => { fastq => '#fastq#', is_paired => '#is_paired#', $self->o('read_id_tag') => '#'.$self->o('read_id_tag').'#'}},
            },
        -rc_name    => '5GB',
      },
            {
        -logic_name => 'bwa2bam_10GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -can_be_empty => 1,
        -parameters => {
                         sampe_options => '-A -a '.$self->o('maxintron'),
                         samse_options => '',
                         min_paired => $self->o('read_min_paired'),
                         min_mapped => $self->o('read_min_mapped'),
                         header_file => '#wide_output_dir#/#'.$self->o('read_id_tag').'#_header.h',
                         bam_prefix => $self->o('read_id_tag'),
                         email => $self->o('email_address'),
                         disconnect_jobs => 1,
                       },
        -flow_into => {
            1 => [ ':////accu?filename=[]' ],
            },
        -rc_name    => '10GB',
      },

            {
        -logic_name => 'merged_tissue_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
        -parameters => {
                         java       => 'java',
                         java_options  => '-Xmx2g',
                         # If 0, do not use multithreading, faster but can use more memory.
                         # If > 0, tells how many cpu to use for samtools or just to use multiple cpus for picard
                         use_threading => $self->o('use_threads'),

                         # Path to MergeSamFiles.jar
                         picard_lib    => $self->o('picard_lib_jar'),
                         # Use this default options for Picard: 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'
                         # You will need to change the options if you want to use samtools for merging
                         options       => 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT',
                         # target_db is the database where we will write the files in the data_file table
                         # You can use store_datafile => 0, if you don't want to store the output file
                         target_db => $self->o('rough_db'),
                         assembly_name => $self->o('assembly_name'),
                         rnaseq_data_provider => $self->o('rnaseq_data_provider'),
                         disconnect_jobs => 1,
                       },
        -rc_name    => '3GB_multithread',
        -flow_into => {
            1 => [ ':////accu?filename=[]' ],
            },
      },
            {
        -logic_name => 'merged_bam_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
        -parameters => {
                         java       => 'java',
                         java_options  => '-Xmx2g',
                         # If 0, do not use multithreading, faster but can use more memory.
                         # If > 0, tells how many cpu to use for samtools or just to use multiple cpus for picard
                         use_threading => $self->o('use_threads'),

                         # Path to MergeSamFiles.jar
                         picard_lib    => $self->o('picard_lib_jar'),
                         # Use this default options for Picard: 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'
                         # You will need to change the options if you want to use samtools for merging
                         options       => 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT',
                         # target_db is the database where we will write the files in the data_file table
                         # You can use store_datafile => 0, if you don't want to store the output file
                         target_db => $self->o('rough_db'),
                         assembly_name => $self->o('assembly_name'),
                         rnaseq_data_provider => $self->o('rnaseq_data_provider'),
                         disconnect_jobs => 1,
                       },
        -rc_name    => '3GB_multithread',
        -flow_into => {
                        1 => ['create_analyses_type_job'],
                      },
      },
            {
              -logic_name => 'create_analyses_type_job',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -rc_name    => '1GB',
              -parameters => {
                inputlist => ['gene', 'daf', 'ise'],
                column_names => ['type'],
              },
              -flow_into => {
                '2->A' => [ 'create_analyses_job'],
                'A->1' => ['create_header_intron'],
              },
            },
            {
              -logic_name => 'create_analyses_job',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -rc_name    => '1GB',
              -parameters => {
                inputquery => q/SELECT DISTINCT(CONCAT('#species#_', LOWER(SM), '_rnaseq_#type#')) FROM csv_data/,
                column_names => ['logic_name'],
                species => $self->o('species'),
              },
              -flow_into => {
                2 => {'create_analyses' => {analyses => [{'-logic_name' => '#logic_name#'}]}},
              },
            },
            {
              -logic_name => 'create_analyses',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
              -rc_name    => '1GB',
              -parameters => {
                source_type => 'list',
                target_db => $self->o('rough_db'),
              },
            },
            {
        -logic_name => 'create_header_intron',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name    => '1GB',
        -parameters => {
                         cmd => '#wide_samtools# view -H #filename# | grep -v @SQ | grep -v @HD > #wide_output_dir#/merged_header.h',
                       },
        -flow_into => {
            '1->A' => [ 'create_top_level_input_ids'],
            'A->1' => ['sam2bam'],
            },
      },
            {
        -logic_name => 'create_top_level_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
                         iid_type => 'slice',
                         coord_system_name => 'toplevel',
                         slice => 1,
                         include_non_reference => 0,
                         top_level => 1,
                         target_db => $self->o('rough_db'),
                       },
        -flow_into => {
                        2 => {'dispatch_toplevel' => {'iid' => '#iid#', alignment_bam_file => '#filename#'}},
                      },
      },
            {
        -logic_name => 'dispatch_toplevel',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name    => '1GB',
        -batch_size => 1000,
        -parameters => {
                         cmd => 'EXIT_CODE=1; if [ "`echo #iid# | cut -d \':\' -f1`" = "chromosome" ]; then EXIT_CODE=2; else EXIT_CODE=0;fi; exit $EXIT_CODE',
                         return_codes_2_branches => {'2' => 2},
                       },
        -flow_into => {
                        1 => ['rough_transcripts'],
                        2 => ['rough_transcripts_5GB'],
                      },
      },
            {
        -logic_name => 'rough_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -parameters => {
                         logic_name => 'rough_transcripts',
                         output_db    => $self->o('rough_db'),
                         dna_db    => $self->o('dna_db'),
                         min_length => 300,
                         min_exons  =>   1,
                         max_intron_length => $self->o('maxintron'),
                         min_single_exon_length => 1000,
                         min_span   =>   1.5,
                         paired => $self->o('paired'),
                         pairing_regex => $self->o('pairing_regex'),
                       },
        -rc_name    => '2GB_rough',
        -flow_into => {
                        1 => ['create_bam2introns_input_ids'],
                        -1 => {'rough_transcripts_5GB' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
                        -2 => {'rough_transcripts_5GB' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
                      },
      },
            {
        -logic_name => 'rough_transcripts_5GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -can_be_empty => 1,
        -parameters => {
                         logic_name => 'rough_transcripts',
                         output_db    => $self->o('rough_db'),
                         dna_db    => $self->o('dna_db'),
                         min_length => 300,
                         min_exons  =>   1,
                         max_intron_length => $self->o('maxintron'),
                         min_single_exon_length => 1000,
                         min_span   =>   1.5,
                         paired => $self->o('paired'),
                         pairing_regex => $self->o('pairing_regex'),
                       },
        -rc_name    => '5GB_rough',
        -flow_into => {
                        1 => ['create_bam2introns_input_ids'],
                        -1 => {'rough_transcripts_15GB' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
                        -2 => {'rough_transcripts_15GB' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
                      },
      },
            {
        -logic_name => 'rough_transcripts_15GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -can_be_empty => 1,
        -parameters => {
                         logic_name => 'rough_transcripts',
                         output_db    => $self->o('rough_db'),
                         dna_db    => $self->o('dna_db'),
                         min_length => 300,
                         min_exons  =>   1,
                         max_intron_length => $self->o('maxintron'),
                         min_single_exon_length => 1000,
                         min_span   =>   1.5,
                         paired => $self->o('paired'),
                         pairing_regex => $self->o('pairing_regex'),
                       },
        -rc_name    => '15GB_rough',
        -flow_into => {
                        1 => ['create_bam2introns_input_ids'],
                      },
      },
            {
        -logic_name => 'create_bam2introns_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'slice_to_feature_ids',
                         target_db => $self->o('rough_db'),
                         feature_type => 'gene',
                         logic_name => ['rough_transcripts'],
                         use_stable_ids => 1,
                         create_stable_ids => 1,
                         stable_id_prefix => 'RNASEQ',
                       },
        -rc_name    => '1GB_rough',
        -batch_size => 100,
        -flow_into => {
                        2 => {'bam2introns' => {iid => '#iid#', bam_file => '#alignment_bam_file#'}},
                      },
      },
            {
        -logic_name => 'bam2introns',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -parameters => {
                         program_file => $self->o('splicing_aligner'),
                         input_db => $self->o('rough_db'),
                         dna_db => $self->o('dna_db'),
                         missmatch => 6,
                         word_length => 10,
                         saturate_threshold => 10000,
                         mask => 1,
                         percent_id => 97,
                         coverage => 90,
                         fullseq   => 1,
                         max_transcript => 1000000,
                         batch_size => 10000,
                         maxintron => $self->o('maxintron'),
                       },
        -rc_name    => '2GB_introns',
        -flow_into => {
                        1 => [':////accu?filename=[]'],
                        2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
                        -1 => {'bam2introns_5GB' => {iid => '#iid#', bam_file => '#bam_file#'}},
                        -2 => {'bam2introns_5GB' => {iid => '#iid#', bam_file => '#bam_file#'}},
                      },
      },
            {
        -logic_name => 'bam2introns_5GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -can_be_empty => 1,
        -parameters => {
                         program_file => $self->o('splicing_aligner'),
                         input_db => $self->o('rough_db'),
                         dna_db => $self->o('dna_db'),
                         missmatch => 6,
                         word_length => 10,
                         saturate_threshold => 10000,
                         mask => 1,
                         percent_id => 97,
                         coverage => 90,
                         fullseq   => 1,
                         max_transcript => 1000000,
                         batch_size => 10000,
                         maxintron => $self->o('maxintron'),
                       },
        -rc_name    => '5GB_introns',
        -flow_into => {
                        1 => [':////accu?filename=[]'],
                        2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
                        -1 => {'bam2introns_10GB' => {iid => '#iid#', bam_file => '#bam_file#'}},
                        -2 => {'bam2introns_10GB' => {iid => '#iid#', bam_file => '#bam_file#'}},
                      },
      },
            {
        -logic_name => 'bam2introns_10GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -can_be_empty => 1,
        -parameters => {
                         program_file => $self->o('splicing_aligner'),
                         input_db => $self->o('rough_db'),
                         dna_db => $self->o('dna_db'),
                         missmatch => 6,
                         word_length => 10,
                         saturate_threshold => 10000,
                         mask => 1,
                         percent_id => 97,
                         coverage => 90,
                         fullseq   => 1,
                         max_transcript => 1000000,
                         batch_size => 10000,
                         maxintron => $self->o('maxintron'),
                       },
        -rc_name    => '10GB_introns',
        -flow_into => {
                        1 => [':////accu?filename=[]'],
                        2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
                      },
      },
            {
        -logic_name => 'sam2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSam2Bam',
        -parameters => {
                         regex => '.sam',
                         headerfile => '#wide_output_dir#/merged_header.h',
                         disconnect_jobs => 1,
                       },
        -rc_name    => '2GB',
        -flow_into => ['check_and_delete_broken_duplicated'],
      },
      {
        -logic_name => 'check_and_delete_broken_duplicated',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects',
        -parameters => {
                         target_db => $self->o('rough_db'),
                         check_support => 0,
                       },
        -rc_name    => '4GB',
        -flow_into => ['create_refine_db'],
      },
      {
        -logic_name => 'create_refine_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('rough_db'),
                         target_db => $self->o('refine_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name => '1GB',
        -flow_into => ['create_blast_db'],
      },

      {
        -logic_name => 'create_blast_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('refine_db'),
                         target_db => $self->o('blast_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name => '1GB',
        -flow_into => ['create_ccode_config'],
      },
      {
              -logic_name => 'create_ccode_config',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRefineConfig',
              -parameters => {
                     single_tissue => $self->o('single_tissue'),
                     sample_column => $self->o('read_group_tag'),
                     sample_id_column => $self->o('read_id_tag'),
                     csvfile_table => $self->o('summary_csv_table'),
                     input_db => $self->o('rough_db'),
                     dna_db => $self->o('dna_db'),
                     output_db => $self->o('refine_db'),
                     # write the intron features into the OUTPUT_DB along with the models
                     write_introns => 1,
                     # maximum number of times to loop when building all possible paths through the transcript
                     max_recursions => 10000000000000,
                     # analysis logic_name for the dna_align_features to fetch from the INTRON_DB
                     # If left blank all features will be fetched
                     logicname => [],
                     # logic name of the gene models to fetch
                     model_ln  => '',
                     # penalty for removing a retined intron
                     retained_intron_penalty => 2,
                     #Remove introns that overlap X introns
                     filter_on_overlap => 0,
                     # minimum size for an intron
                     min_intron_size  => 30,
                     max_intron_size  => $self->o('maxintron'),
                     # biotype to give to single exon models if left blank single exons are ignored
                     # minimum single exon size (bp)
                     min_single_exon => 1000,
                     # minimum percentage of single exon length that is coding
                     single_exon_cds => 66,
                     # Intron with most support determines the splice sites for an internal exon
                     # lower scoring introns with different splice sites are rejected
                     strict_internal_splice_sites => 1,
                     # In some species alternate splice sites for end exons seem to be common
                     strict_internal_end_exon_splice_sites => 1,
                     # biotypes to give gene models if left blank these models will not get written to the output database
                     # best score - model with most supporting intron features
                     # all other possible models
                     # max number of other models to make - blank = all
                     other_num      => '10',
                     # max number of other models to process - blank = all
                     max_num      => '1000',
                     other_isoforms => $self->o('other_isoforms'),
                     # biotype to label bad models ( otherwise they are not written )
                     # do you want to trim UTR
                     trim_utr => 1,
                     # config for trimming UTR
                     max_3prime_exons => 2,
                     max_3prime_length => 5000,
                     max_5prime_exons => 3,
                     max_5prime_length => 1000,
                     # % of average intron score that a UTR intron must have
                     reject_intron_cutoff => 5,
                             },

              -rc_name          => '1GB',
              -flow_into => {
                              2 => ['create_ccode_input_ids'],
                            },
            },
        {
                -logic_name => 'create_ccode_input_ids',
                -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
                -rc_name    => '1GB',
                -parameters => {
                                 iid_type => 'slice',
                                 coord_system_name => 'toplevel',
                                 slice => 1,
                                 include_non_reference => 0,
                                 top_level => 1,
                                 feature_constraint => 1,
                                 feature_type => 'gene',
                                 target_db => $self->o('rough_db'),
                               },
                -flow_into => {
                                2 => {'refine_genes' => {iid => '#iid#', logic_name => '#logic_name#', config_file => '#config_file#'}},
                              },
              },

        {
              -logic_name => 'refine_genes',
                -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
                -parameters => {
                    cmd => $self->o('refine_ccode_exe').($self->o('use_ucsc_naming') ? ' -u ' : ' ').($self->o('use_threads') ? ' -t '.$self->o('use_threads').' ' : ' ').'-c #config_file# -i #iid# -l #logic_name# -v 0',
                    return_codes_2_branches => {
                      42 => 2,
                    },
                },
                -rc_name => '2GB_refine',
                -flow_into => {
                    1 => ['blast_rnaseq'],
                    -1 => {'refine_genes_20GB' => {iid => '#iid#', config_file => '#config_file#', logic_name => '#logic_name#'}},
                },
          },
        {
              -logic_name => 'refine_genes_20GB',
                -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
                -parameters => {
                    cmd => $self->o('refine_ccode_exe').($self->o('use_ucsc_naming') ? ' -u ' : ' ').($self->o('use_threads') ? ' -t '.$self->o('use_threads').' ' : ' ').'-c #config_file# -i #iid# -l #logic_name# -v 0',
                    return_codes_2_branches => {
                      42 => 2,
                    },
                },
                -rc_name => '20GB_refine',
                -flow_into => {
                    1 => ['blast_rnaseq'],
                },
          },

      {
        -logic_name => 'blast_rnaseq',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
        -parameters => {

            input_db => $self->o('refine_db'),
            output_db => $self->o('blast_db'),
            dna_db => $self->o('dna_db'),

            # path to index to fetch the sequence of the blast hit to calculate % coverage
            indicate_index => $self->o('uniprotindex'),
            uniprot_index => [$self->o('uniprotdb')],
            blast_program => $self->o('blastp'),
            %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
            commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40',
                      },
        -rc_name => '2GB_blast',
      },
    );
    foreach my $analyses (@analysis) {
        $analyses->{-max_retry_count} = 1 unless (exists $analyses->{-max_retry_count});
    }
    return \@analysis;
}

sub resource_classes {
    my $self = shift;

    return {
        %{ $self->SUPER::resource_classes() },  # inherit other stuff from the base class
      '1GB' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 1000, [$self->default_options->{'pipe_db_server'}])},
      '1GB_rough' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 1000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}])},
      '2GB_rough' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}])},
      '5GB_rough' => { LSF => $self->lsf_resource_builder($self->default_options->{long_queue}, 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}])},
      '15GB_rough' => { LSF => $self->lsf_resource_builder($self->default_options->{long_queue}, 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}])},
      '2GB_blast' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'refine_db_server'}, $self->default_options->{'blast_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
      '2GB' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 2000, [$self->default_options->{'pipe_db_server'}])},
      '4GB' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}])},
      '5GB' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 5000, [$self->default_options->{'pipe_db_server'}])},
      '2GB_introns' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}, $self->default_options->{'dna_db_server'}])},
      '2GB_refine' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_db_server'}])},
      '5GB_introns' => { LSF => $self->lsf_resource_builder($self->default_options->{long_queue}, 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}, $self->default_options->{'dna_db_server'}])},
      '10GB_introns' => { LSF => $self->lsf_resource_builder($self->default_options->{long_queue}, 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}, $self->default_options->{'dna_db_server'}])},
      '3GB_multithread' => { LSF => $self->lsf_resource_builder($self->default_options->{long_queue}, 3000, [$self->default_options->{'pipe_db_server'}], undef, $self->default_options->{'use_threads'})},
      '5GB_multithread' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 5000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
      '10GB_multithread' => { LSF => $self->lsf_resource_builder($self->default_options->{long_queue}, 10000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
      '20GB_multithread' => { LSF => $self->lsf_resource_builder($self->default_options->{long_queue}, 20000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
      '5GB' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 5000, [$self->default_options->{'pipe_db_server'}])},
      '10GB' => { LSF => $self->lsf_resource_builder($self->default_options->{long_queue}, 10000, [$self->default_options->{'pipe_db_server'}])},
      '5GB_refine' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_db_server'}])},
      '20GB_refine' => { LSF => $self->lsf_resource_builder($self->default_options->{normal_queue}, 20000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_db_server'}])},
    };
}

1;
