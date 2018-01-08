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

=cut

package HiveRNASeq_conf;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

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
        'pipeline_name'              => '',

        'user_r'                     => '',
        'user'                     => '',
        'password'                   => '',
        'port'                       => '',
        species    => '',

        'pipe_dbname'                => $ENV{USER}.'_'.$self->o('pipeline_name').'_hive',
        'dna_dbname'                 => '',
        'blast_output_dbname'     => $ENV{USER}.'_hive_'.$self->o('species').'_blast',
        'refine_output_dbname'     => $ENV{USER}.'_hive_'.$self->o('species').'_refine',
        'rough_output_dbname'    => $ENV{USER}.'_hive_'.$self->o('species').'_rough',

        'pipe_db_server'             => '',
        'dna_db_server'              => '',
        'blast_output_db_server'  => '',
        'refine_output_db_server'  => '',
        'rough_output_db_server' => '',
        'genome_file'                => 'genome/genome.fa',
        'use_ucsc_naming' => 0,

        'clone_db_script_path'       => $ENV{ENSCODE}.'/ensembl-analysis/scripts/clone_database.ksh',
        'create_type'       => 'clone',
        'rnaseq_summary_file'         => '',


        'samtools' => 'samtools',
        'picard_lib_jar' => 'picard.jar',
        'short_read_aligner'    => 'bwa',
        'input_dir'    => '',
         output_dir => '',
        'merge_dir' => '',
        'sam_dir' => '',
        'sequence_dump_script' => $ENV{ENSCODE}.'/ensembl-analysis/scripts/sequence_dump.pl',
        # Use this option to change the delimiter for your summary data
        # file.
        summary_file_delimiter => '\t',
        summary_csv_table => 'csv_data',
        assembly_name => '',

        # Blast database for comparing the final models to.
        uniprotdb => '',

        # Index for the blast database.
        uniprotindex => '',

        blastp => 'blastall',
        # blast used, it can be either ncbi or wu, it is overriding the -type value from BLAST_PARAMS
        blast_type => 'ncbi',

        splicing_aligner => 'exonerate-0.9.0',

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

        'blast_output_db' => {
                           -dbname => $self->o('blast_output_dbname'),
                           -host   => $self->o('blast_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user'),
                           -pass   => $self->o('password'),
                         },

        'refine_output_db' => {
                           -dbname => $self->o('refine_output_dbname'),
                           -host   => $self->o('refine_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user'),
                           -pass   => $self->o('password'),
                        },

        'rough_output_db' => {
                           -dbname => $self->o('rough_output_dbname'),
                           -host   => $self->o('rough_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user'),
                           -pass   => $self->o('password'),
                         },
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    my $output_sam_dir = $self->o('sam_dir') ? $self->o('sam_dir') :$self->o('output_dir').'/SAM';
    my $merge_dir = $self->o('merge_dir') ? $self->o('merge_dir') : $self->o('output_dir').'/merge_out';
    my $genome_file = $self->o('genome_file') =~ /^\// ? $self->o('genome_file') : $self->o('input_dir').'/'.$self->o('genome_file');
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
                         wide_intron_bam_file => $self->o('output_dir').'/introns',
    };
}


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
      $self->db_cmd('CREATE TABLE '.$self->default_options->{'summary_csv_table'}." ($tables)"),
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
        -meadow_type => 'LOCAL',
        -parameters => {
            cmd => 'EXIT_CODE=0; for F in #wide_short_read_aligner# #wide_samtools# '.join (' ', $self->o('splicing_aligner'), $self->o('clone_db_script_path'), $self->o('sequence_dump_script'), $self->o('blastp')).'; do which "$F"; if [ "$?" == 1 ]; then EXIT_CODE=1;fi; done; for D in #wide_output_dir# #wide_input_dir# #wide_merge_dir# #wide_output_sam_dir# `dirname #wide_genome_file#`; do mkdir -p "$D"; done; exit $EXIT_CODE',
        },
        -input_ids => [{}],
        -flow_into => {
            '1->A' => ['create_rnaseq_genome_file'],
            'A->1' => ['parse_summary_file'],
        },
  },
 {
      -logic_name => 'create_rnaseq_genome_file',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '1GB',
        -parameters => {
            cmd => 'if [ ! -e "#wide_genome_file#" ]; then perl '.$self->o('sequence_dump_script').' -dbhost '.$self->o('dna_db_server').' -dbuser '.$self->o('user_r').' -dbport '.$self->o('port').' -dbname '.$self->o('dna_dbname').' -coord_system_name '.$self->o('assembly_name').' -toplevel -onefile -header rnaseq -filename #wide_genome_file#;fi',
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
      -logic_name => 'parse_summary_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable',
        -meadow_type => 'LOCAL',
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
            inputquery => join(' ', 'SELECT', $self->o('read_group_tag'), ',', $self->o('read_id_tag'), ', is_paired', 'FROM', $self->o('summary_csv_table'), 'WHERE', $self->o('read_group_tag'), '= "#sample_name#"'),
            column_names => [$self->o('read_group_tag'), $self->o('read_id_tag'), 'is_paired'],
                       },
        -meadow_type    => 'LOCAL',
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
        -meadow_type    => 'LOCAL',
        -flow_into => {
            '2->A' => ['bwa', 'create_header_files'],
            'A->1' => ['bwa2bam'],
            },
      },
  {
      -logic_name => 'create_header_files',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -meadow_type => 'LOCAL',
        -parameters => {
            cmd => 'if [ ! -e "#wide_output_dir#/#'.$self->o('read_id_tag').'#_header.h" ]; then printf "'.$header_line.'" > #wide_output_dir#/#'.$self->o('read_id_tag').'#_header.h; fi',
        },
  },
            {
        -logic_name => 'bwa',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -flow_into => {
            1 => [ ':////accu?filename=[]' ],
            -1 => [ 'bwa_failed' ],
            -2 => [ 'bwa_failed' ],
            },
        -rc_name    => '5GB_multithread',
      },
            {
        -logic_name => 'bwa_failed',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -flow_into => {
            1 => [ ':////accu?filename=[]' ],
            },
        -rc_name    => '10GB_multithread',
      },
            {
        -logic_name => 'bwa2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -parameters => {
                         sampe_options => '-A -a 200000',
                         samse_options => '',
                         min_paired => $self->o('read_min_paired'),
                         min_mapped => $self->o('read_min_mapped'),
                         header_file => '#wide_output_dir#/#'.$self->o('read_id_tag').'#_header.h',
                         bam_prefix => $self->o('read_id_tag'),
                       },
        -flow_into => {
            1 => [ ':////accu?filename=[]' ],
            -1 => {'bwa2bam_10GB' => { filename => '#filename#', is_paired => '#is_paired#', $self->o('read_id_tag') => '#'.$self->o('read_id_tag').'#', $self->o('read_group_tag') => '#'.$self->o('read_group_tag').'#',}},
            -2 => {'bwa2bam_10GB' => { filename => '#filename#', is_paired => '#is_paired#', $self->o('read_id_tag') => '#'.$self->o('read_id_tag').'#', $self->o('read_group_tag') => '#'.$self->o('read_group_tag').'#',}},
            },
        -rc_name    => '5GB',
      },
            {
        -logic_name => 'bwa2bam_10GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -parameters => {
                         sampe_options => '-A -a 200000',
                         samse_options => '',
                         min_paired => $self->o('read_min_paired'),
                         min_mapped => $self->o('read_min_mapped'),
                         header_file => '#wide_output_dir#/#'.$self->o('read_id_tag').'#_header.h',
                         bam_prefix => $self->o('read_id_tag'),
                       },
        -flow_into => {
            1 => [ ':////accu?filename=[]' ],
            },
        -rc_name    => '10GB',
      },
      {
        -logic_name => 'create_rough_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('rough_output_db'),
                         create_type => $self->o('create_type'),
                         script_path => $self->o('clone_db_script_path'),
                       },
        -meadow => 'LOCAL',
        -input_ids => [{}],
      },

      {
        -logic_name => 'create_refine_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('refine_output_db'),
                         create_type => $self->o('create_type'),
                         script_path => $self->o('clone_db_script_path'),
                       },
        -meadow => 'LOCAL',
        -input_ids => [{}],
      },

      {
        -logic_name => 'create_blast_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('blast_output_db'),
                         create_type => $self->o('create_type'),
                         script_path => $self->o('clone_db_script_path'),
                       },
        -meadow => 'LOCAL',
        -input_ids => [{}],
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
                       },
        -rc_name    => '3GB_multithread',
        -flow_into => {
                        1 => ['create_header_intron', 'clean_sai_files'],
                      },
      },
            {
        -logic_name => 'clean_sai_files',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -meadow_type => 'LOCAL',
        -parameters => {
                         cmd => 'rm #wide_output_dir#/*.sai',
                       },
      },
            {
        -logic_name => 'create_header_intron',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -meadow_type => 'LOCAL',
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
                         target_db => $self->o('rough_output_db'),
                       },
        -wait_for => ['create_rough_output_db'],
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
                         output_db    => $self->o('rough_output_db'),
                         dna_db    => $self->o('dna_db'),
                         alignment_bam_file => '#wide_merge_dir#/merged.bam',
                         min_length => 300,
                         min_exons  =>   1,
                         max_intron_length => 200000,
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
        -parameters => {
                         logic_name => 'rough_transcripts',
                         output_db    => $self->o('rough_output_db'),
                         dna_db    => $self->o('dna_db'),
                         alignment_bam_file => '#wide_merge_dir#/merged.bam',
                         min_length => 300,
                         min_exons  =>   1,
                         max_intron_length => 200000,
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
        -parameters => {
                         logic_name => 'rough_transcripts',
                         output_db    => $self->o('rough_output_db'),
                         dna_db    => $self->o('dna_db'),
                         alignment_bam_file => '#wide_merge_dir#/merged.bam',
                         min_length => 300,
                         min_exons  =>   1,
                         max_intron_length => 200000,
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
                         target_db => $self->o('rough_output_db'),
                         feature_type => 'gene',
                         logic_name => ['rough_transcripts'],
                         use_stable_ids => 1,
                         create_stable_ids => 1,
                         stable_id_prefix => 'RNASEQ',
                       },
        -rc_name    => '1GB_rough',
        -batch_size => 100,
        -flow_into => {
                        2 => ['bam2introns'],
                      },
      },
            {
        -logic_name => 'bam2introns',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -parameters => {
                         program_file => $self->o('splicing_aligner'),
                         input_db => $self->o('rough_output_db'),
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
                         bam_file => '#wide_merge_dir#/merged.bam',
                       },
        -rc_name    => '2GB_introns',
        -flow_into => {
                        1 => [':////accu?filename=[]'],
                        2 => ['bam2introns'],
                        -1 => ['bam2introns_5GB'],
                        -2 => ['bam2introns_5GB'],
                      },
      },
            {
        -logic_name => 'bam2introns_5GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -parameters => {
                         program_file => $self->o('splicing_aligner'),
                         input_db => $self->o('rough_output_db'),
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
                         bam_file => '#wide_merge_dir#/merged.bam',
                       },
        -rc_name    => '5GB_introns',
        -flow_into => {
                        1 => [':////accu?filename=[]'],
                        2 => ['bam2introns'],
                        -1 => ['bam2introns_10GB'],
                        -2 => ['bam2introns_10GB'],
                      },
      },
            {
        -logic_name => 'bam2introns_10GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -parameters => {
                         program_file => $self->o('splicing_aligner'),
                         input_db => $self->o('rough_output_db'),
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
                         bam_file => '#wide_merge_dir#/merged.bam',
                       },
        -rc_name    => '10GB_introns',
        -flow_into => {
                        1 => [':////accu?filename=[]'],
                        2 => ['bam2introns'],
                      },
      },
            {
        -logic_name => 'sam2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSam2Bam',
        -parameters => {
                         regex => '.sam',
                         headerfile => '#wide_output_dir#/merged_header.h',
                       },
        -rc_name    => '2GB',
        -flow_into => ['create_top_level_input_ids_again'],
      },
            {
        -logic_name => 'create_top_level_input_ids_again',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
                         iid_type => 'slice',
                         coord_system_name => 'toplevel',
                         slice => 1,
                         include_non_reference => 0,
                         top_level => 1,
                         target_db => $self->o('rough_output_db'),
                       },
        -flow_into => {
                        2 => ['create_refine_genes_jobs'],
                      },
      },
            {
        -logic_name => 'create_refine_genes_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateRefineGenesJobs',
        -parameters => {
                         single_tissue => $self->o('single_tissue'),
                         sample_column => $self->o('read_group_tag'),
                         sample_id_column => $self->o('read_id_tag'),
                         csvfile_table => $self->o('summary_csv_table'),
                       },
        -rc_name    => '1GB',
        -batch_size => 100,
        -flow_into => {
                        2 => ['refine_genes'],
                      },
      },

      {
        -logic_name => 'refine_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRefineSolexaGenes',
        -parameters => {
               input_db => $self->o('rough_output_db'),
               dna_db => $self->o('dna_db'),
               output_db => $self->o('refine_output_db'),
               # write the intron features into the OUTPUT_DB along with the models
               write_introns => 1,
               # maximum number of times to loop when building all possible paths through the transcript
               max_recursions => 100000,
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
               max_intron_size  => 200000,
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

        -rc_name          => '2GB_refine',
        -wait_for => ['create_refine_output_db'],
        -flow_into => {
                        1 => ['blast_rnaseq'],
                        -1 => ['refine_genes_5GB'],
                        -2 => ['refine_genes_5GB'],
                      },
      },
      {
        -logic_name => 'refine_genes_5GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRefineSolexaGenes',
        -parameters => {
               input_db => $self->o('rough_output_db'),
               dna_db => $self->o('dna_db'),
               output_db => $self->o('refine_output_db'),
               # write the intron features into the OUTPUT_DB along with the models
               write_introns => 1,
               # maximum number of times to loop when building all possible paths through the transcript
               max_recursions => 100000,
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
               max_intron_size  => 200000,
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

        -rc_name          => '5GB_refine',
        -flow_into => {
                        1 => ['blast_rnaseq'],
                        -1 => ['refine_genes_15GB'],
                        -2 => ['refine_genes_15GB_base'],
                      },
      },
      {
        -logic_name => 'refine_genes_15GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRefineSolexaGenes',
        -parameters => {
               input_db => $self->o('rough_output_db'),
               dna_db => $self->o('dna_db'),
               output_db => $self->o('refine_output_db'),
               # write the intron features into the OUTPUT_DB along with the models
               write_introns => 1,
               # maximum number of times to loop when building all possible paths through the transcript
               max_recursions => 100000,
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
               max_intron_size  => 200000,
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

        -rc_name          => '15GB_refine',
        -flow_into => {
                        1 => ['blast_rnaseq'],
                        -1 => ['refine_genes_30GB_base'],
                        -2 => ['refine_genes_30GB_base'],
                      },
      },
      {
        -logic_name => 'refine_genes_15GB_base',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRefineSolexaGenes',
        -parameters => {
               input_db => $self->o('rough_output_db'),
               dna_db => $self->o('dna_db'),
               output_db => $self->o('refine_output_db'),
               # write the intron features into the OUTPUT_DB along with the models
               write_introns => 1,
               # maximum number of times to loop when building all possible paths through the transcript
               max_recursions => 100000,
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
               max_intron_size  => 200000,
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

        -rc_name          => '15GB_refine_base',
        -flow_into => {
                        1 => ['blast_rnaseq'],
                        -1 => ['refine_genes_30GB_base'],
                        -2 => ['refine_genes_30GB_base'],
                      },
      },
      {
        -logic_name => 'refine_genes_30GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRefineSolexaGenes',
        -parameters => {
               input_db => $self->o('rough_output_db'),
               dna_db => $self->o('dna_db'),
               output_db => $self->o('refine_output_db'),
               # write the intron features into the OUTPUT_DB along with the models
               write_introns => 1,
               # maximum number of times to loop when building all possible paths through the transcript
               max_recursions => 100000,
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
               max_intron_size  => 200000,
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

        -rc_name          => '30GB_refine',
        -flow_into => {
                        1 => ['blast_rnaseq'],
                        -2 => ['refine_genes_30GB_base'],
                      },
      },
      {
        -logic_name => 'refine_genes_30GB_base',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRefineSolexaGenes',
        -parameters => {
               input_db => $self->o('rough_output_db'),
               dna_db => $self->o('dna_db'),
               output_db => $self->o('refine_output_db'),
               # write the intron features into the OUTPUT_DB along with the models
               write_introns => 1,
               # maximum number of times to loop when building all possible paths through the transcript
               max_recursions => 100000,
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
               max_intron_size  => 200000,
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

        -rc_name          => '30GB_refine_base',
        -flow_into => {
                        1 => ['blast_rnaseq'],
                      },
      },

      {
        -logic_name => 'blast_rnaseq',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
        -parameters => {

            input_db => $self->o('refine_output_db'),
            output_db => $self->o('blast_output_db'),
            dna_db => $self->o('dna_db'),

            # path to index to fetch the sequence of the blast hit to calculate % coverage
            indicate_index => $self->o('uniprotindex'),
            uniprot_index => [$self->o('uniprotdb')],
            blast_program => $self->o('blastp'),
            type => $self->o('blast_type'),
            config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanPep'),
                      },
        -rc_name => '2GB_blast',
        -wait_for => ['create_blast_output_db'],
      },
    );
    foreach my $analyses (@analysis) {
        $analyses->{-meadow_type} = 'LSF' unless (exists $analyses->{-meadow_type});
        $analyses->{-max_retry_count} = 1 unless (exists $analyses->{-max_retry_count});
    }
    return \@analysis;
}

sub resource_classes {
    my $self = shift;

    return {
        %{ $self->SUPER::resource_classes() },  # inherit other stuff from the base class
      '1GB' => { LSF => $self->lsf_resource_builder('normal', 1000, [$self->default_options->{'pipe_db_server'}])},
      '1GB_rough' => { LSF => $self->lsf_resource_builder('normal', 1000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}])},
      '2GB_rough' => { LSF => $self->lsf_resource_builder('normal', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}])},
      '5GB_rough' => { LSF => $self->lsf_resource_builder('long', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}])},
      '15GB_rough' => { LSF => $self->lsf_resource_builder('long', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}])},
      '2GB_blast' => { LSF => $self->lsf_resource_builder('normal', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'refine_output_db_server'}, $self->default_options->{'blast_output_db_server'}])},
      '2GB' => { LSF => $self->lsf_resource_builder('normal', 2000, [$self->default_options->{'pipe_db_server'}])},
      '4GB' => { LSF => $self->lsf_resource_builder('normal', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}])},
      '5GB' => { LSF => $self->lsf_resource_builder('normal', 5000, [$self->default_options->{'pipe_db_server'}])},
      '2GB_introns' => { LSF => $self->lsf_resource_builder('normal', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}])},
      '2GB_refine' => { LSF => $self->lsf_resource_builder('normal', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_output_db_server'}])},
      '5GB_introns' => { LSF => $self->lsf_resource_builder('long', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}])},
      '10GB_introns' => { LSF => $self->lsf_resource_builder('long', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}])},
      '3GB_multithread' => { LSF => $self->lsf_resource_builder('long', 3000, [$self->default_options->{'pipe_db_server'}], undef, 3)},
      '5GB_multithread' => { LSF => $self->lsf_resource_builder('normal', 5000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
      '10GB_multithread' => { LSF => $self->lsf_resource_builder('long', 10000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
      '20GB_multithread' => { LSF => $self->lsf_resource_builder('long', 20000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
      '5GB' => { LSF => $self->lsf_resource_builder('normal', 5000, [$self->default_options->{'pipe_db_server'}])},
      '10GB' => { LSF => $self->lsf_resource_builder('long', 10000, [$self->default_options->{'pipe_db_server'}])},
      '5GB_refine' => { LSF => $self->lsf_resource_builder('long', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_output_db_server'}])},
      '15GB_refine' => { LSF => $self->lsf_resource_builder('long', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_output_db_server'}])},
      '15GB_refine_base' => { LSF => $self->lsf_resource_builder('basement', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_output_db_server'}])},
      '30GB_refine' => { LSF => $self->lsf_resource_builder('long', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_output_db_server'}])},
      '30GB_refine_base' => { LSF => $self->lsf_resource_builder('basement', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'rough_output_db_server'}, $self->default_options->{'dna_db_server'}, $self->default_options->{'refine_output_db_server'}])},
    };
}

1;
