=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package Hive_primate_basic_conf;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

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
        'user_w'                     => '',
        'password'                   => '',
        'port'                       => '',

        'pipe_dbname'                => '',
        'reference_dbname'           => '',
        'dna_dbname'                 => '',
        'blast_output_dbname'     => '',
        'refine_output_dbname'     => '',
        'rough_output_dbname'    => '',

        'pipe_db_server'             => '',
        'reference_db_server'        => '',
        'dna_db_server'              => '',
        'blast_output_db_server'  => '',
        'refine_output_db_server'  => '',
        'rough_output_db_server' => '',
        'killlist_db_server'         => '',
        'output_path'                => '',
        'genome_file'                => '',

        'clone_db_script_path'       => '',
        'repeat_masking_logic_names' => [],
        'rnaseq_summary_file'         => '',


# ALL
        'samtools_path' => '/samtools',
# BWA
        'genomefile'    => '/path/to/genome.fa',
        'short_read_aligner'    => '/path/to/bwa',
        'input_dir'    => '/path/to/fastq',
        'merge_dir' => "/path/to/my/merge/dir",
        # Use this option to change the delimiter for your summary data
        # file.
        delimiter => "\t",

        # Path to the directory containing your
        # Bio::Ensembl::Analysis::Config.
        analysisconfig_dir => "/path/to/my/Bio/EnsEMBL/Analysis/Config",

        # Path to the directory containing your
        # Bio::Ensembl::Pipeline::Config.
        pipelineconfig_dir => "/path/to/my/Bio/EnsEMBL/Pipeline/Config",

        # Path to base directory to merge the BAM files

        # Blast database for comparing the final models to.
        uniprotdb => '/path/to/my/UniprotDB',

        # Index for the blast database.
        uniprotindex => '/path/to/my/Uniprot/index',

        # blast used, it can be either ncbi or wu
        blastp => 'ncbi',

        splicing_aligner => '/path/to/exonerate-0.9.0',
        # Global read length.
        read_length => 0,

        # This is used by bwa2bam.  Set it to 1 if all the reads are paired.
        # Setting it to 0 means it will treat all reads as unpaired and so make long rough models.
        # Note that the PAIRED column value needs setting independantly of this parameter.
        all_paired => 1,

        # If your reads are unpaired you may want to run on slices to avoid
        # making overlong rough models.  If you want to do this, specify a
        # slice length here otherwise it will default to whole chromosomes.
        slice_length => 10000000,

        # Regular expression to allow FastQ files to be correctly paired,
        # for example: file_1.fastq and file_2.fastq could be paired using
        # the expression "(\S+)_r(\d)\.(\S+)".  Need to identify all 3 parts
        # in brackets; the name the read number (1, 2, a, b etc.) and the
        # extension.
        pairing_regex => '(\S+)_r(\d)\.(\S+)',

        # Do you want to make models for the each individual sample as well
        # as for the pooled samples (1/0)?
        single_tissue => 1,

        # What Read group tag would you like to group your samples
        # by? Default = ID
        read_group_tag => 'SM',

        use_threads => 3,

        # Configure the pipeline to use Gsnap rather than BWA and Exonerate.
        use_gsnap => 0,

        # Path to Gsnap binary.
        gsnap_path => "/path/to/gsnap",

        # Please assign some or all columns from the summary file to the
        # some or all of the following categories.  Multiple values can be
        # separted with commas. only ID, SM and FILE are required.  If you
        # have paired reads you must include the PAIRED category or BWA will
        # not be able to pair the reads correctly.

        ####################################################################
        # This is just an example based on the file snippet shown below.  It
        # will vary depending on how your data looks.
        ####################################################################

        # Unique read group identifier (needs to be unique used as the logic
        # name for the analysis). *required
        ID => "1",

        # Sample (use pool name where a pool is being sequenced). *required
        SM => "3",

        # Library.
        LB => "2",

        # Description.
        DS => "4,5,6,7,8",

        # Platform unit (e.g. lane for Illumina or slide for SOLiD); should
        # be a full, unambiguous identifier.
        PU => "",

        # Name of sequencing center producing the read.
        CN => "9",

        # Date the run was produced (ISO 8601 date or date/time).
        DT => "",

        # Platform/technology used to produce the read.
        PL => "",

        # Path from INPUT_DIR to the filename of the RNAseq FastQ file,
        # comma separated values will be replaced with '/'.
        FILE => "16",

        # Global read length can be overridden if it is different for
        # different lanes.
        LENGTH => "",

        # The read pairing, 1 for paired, 0 for unpaired (independent of
        # ALL_PAIRED parameter above) *required for BWA to pair reads
        # correctly.
        PAIRED => "11",

##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

        'pipeline_db' => {
            -dbname => $self->o('pipe_dbname'),
            -host   => $self->o('pipe_db_server'),
            -port   => $self->o('port'),
            -user   => $self->o('user_w'),
            -pass   => $self->o('password'),
            -driver => $self->o('driver'),
        },

        'reference_db' => {
                            -dbname => $self->o('reference_dbname'),
                            -host   => $self->o('reference_db_server'),
                            -port   => $self->o('port'),
                            -user   => $self->o('user_r'),
                          },

        'dna_db' => {
                      -dbname => $self->o('dna_dbname'),
                      -host   => $self->o('dna_db_server'),
                      -port   => $self->o('port'),
                      -user   => $self->o('user_r'),
                    },

        'blast_output_db' => {
                           -dbname => $self->o('blast_output_dbname'),
                           -host   => $self->o('blast_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                         },

        'refine_output_db' => {
                           -dbname => $self->o('refine_output_dbname'),
                           -host   => $self->o('refine_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                        },

        'rough_output_db' => {
                           -dbname => $self->o('rough_output_dbname'),
                           -host   => $self->o('rough_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                         },


        'killlist_db' => {
                           -dbname    => $self->o('killlist_dbname'),
                           -host      => $self->o('killlist_db_server'),
                           -port      => $self->o('port'),
                           -user      => $self->o('user_r'),
                         },
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
                         genome_file => $self->o('genome_file'),
                         input_dir => $self->o('input_dir'),
                         output_dir => $self->o('output_dir'),
                         short_read_aligner => $self->o('short_read_aligner'),
                         samtools => $self->o('samtools'),
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

    my @analysis = (
  {
      -logic_name => 'parse_summary_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -meadow => 'LOCAL',
  },

      {
        -logic_name => 'create_rough_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
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
                         source_db => $self->o('reference_db'),
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
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('blast_output_db'),
                         create_type => $self->o('create_type'),
                         script_path => $self->o('clone_db_script_path'),
                       },
        -meadow => 'LOCAL',
        -input_ids => [{}],
      },

            {
        -logic_name => 'create_bwa_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
                       },
        -meadow    => 'LOCAL',
        -inputlist => [{input_id => 'file.astq',
                        short_read_aligner_options => '-t 3 -i 50 -n 100',
                       },],
        -flow_into => [ {
            '2->A' => 'bwa',
            'A->1' => 'bwa2bam'
            },],
      },
            {
        -logic_name => 'bwa',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -parameters => {
                       },
        -rc_name    => '5GB_multithread',
        -input_ids => [{input_id => 'file.astq',
                        short_read_aligner_options => '-t 3 -i 50 -n 100',
                       },],
      },
            {
        -logic_name => 'bwa2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -parameters => {
                         sampe_options => '-A -a 200000',
                         samse_options => '',
                       },
        -rc_name    => '5GB_multithread',
        -input_ids => [{input_id => 'file.astq',
                        fastq_pair => '',
                        header_file => '/path/to/header.txt',
                        paired => 1,
                        min_paired => 50,
                        min_mapped => 50,},],
      },
            {
        -logic_name => 'tissue_bam_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
        -parameters => {
                         java       => 'java',
                         java_options  => '-Xmx2g',
                         # If 0, do not use multithreading, faster but can use more memory.
                         # If > 0, tells how many cpu to use for samtools or just to use multiple cpus for picard
                         use_threading => 2,

                         # Path to MergeSamFiles.jar
                         picard_lib    => 'picard.jar',
                         # Use this default options for Picard: 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'
                         # You will need to change the options if you want to use samtools for merging
                         options       => 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT',
                       },
        -rc_name    => '3GB_multithread',
        -input_ids => [{
                        input_id => 'sample name',
                        input_files => [],
                        output_files => '',
                       },
                      ],
      },
            {
        -logic_name => 'merged_bam_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
        -parameters => {
                         java       => 'java',
                         java_options  => '-Xmx2g',
                         # If 0, do not use multithreading, faster but can use more memory.
                         # If > 0, tells how many cpu to use for samtools or just to use multiple cpus for picard
                         use_threading => 2,

                         # Path to MergeSamFiles.jar
                         picard_lib    => 'picard.jar',
                         # Use this default options for Picard: 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'
                         # You will need to change the options if you want to use samtools for merging
                         options       => 'MAX_RECORDS_IN_RAM=20000000 CREATE_INDEX=true SORT_ORDER=coordinate ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT',
                       },
        -rc_name    => '3GB_multithread',
        -input_ids => [{
                        input_id => 'sample name',
                        input_files => [],
                        output_files => '',
                       },
                      ],
      },
            {
        -logic_name => 'create_top_level_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -parameters => {
                         sampe_options => '-A -a 200000',
                         samse_options => '',
                         short_read_aligner => $self->o('short_read_aligner'),
                         short_read_aligner_options => $self->o('short_read_aligner_options'),
                         genome_file => $self->o('genome_file'), #Wide
                         input_dir => $self->o('input_dir'), #Wide
                         output_dir => $self->o('output_dir'), #Wide
                       },
        -rc_name    => '5GB_multithread',
        -input_ids => [{input_id => 'file.astq',
                        fastq_pair => '',
                        header_file => '/path/to/header.txt',
                        paired => 1,
                        min_paired => 50,
                        min_mapped => 50,},],
      },
            {
        -logic_name => 'rough_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -parameters => {
                         output_db    => $self->o('rough_output_db', '-dbname'),
                         alignment_bam_file => '/path/to/my/sorted/indexed/bam_file.bam',
                         min_length => 300,
                         min_exons  =>   1,
                         max_intron_length => 200000,
                         min_single_exon_length => 1000,
                         min_span   =>   1.5,
                       },
        -rc_name    => '2GB',
        -input_ids => [{input_id => 'file.astq',
                        fastq_pair => '',
                        header_file => '/path/to/header.txt',
                        paired => 1,
                        pairing_regex => '_\d',
                        min_paired => 50,
                        min_mapped => 50,},],
      },
            {
        -logic_name => 'create_bam2introns_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -parameters => {
                         sampe_options => '-A -a 200000',
                         samse_options => '',
                         short_read_aligner => $self->o('short_read_aligner'),
                         short_read_aligner_options => $self->o('short_read_aligner_options'),
                         genome_file => $self->o('genome_file'), #Wide
                         input_dir => $self->o('input_dir'), #Wide
                         output_dir => $self->o('output_dir'), #Wide
                       },
        -rc_name    => '1GB',
        -input_ids => [{input_id => 'file.astq',
                        fastq_pair => '',
                        header_file => '/path/to/header.txt',
                        paired => 1,
                        min_paired => 50,
                        min_mapped => 50,},],
      },
            {
        -logic_name => 'bam2introns',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -parameters => {
                         sampe_options => '-A -a 200000',
                         samse_options => '',
                         short_read_aligner => $self->o('short_read_aligner'),
                         short_read_aligner_options => $self->o('short_read_aligner_options'),
                         genome_file => $self->o('genome_file'), #Wide
                         input_dir => $self->o('input_dir'), #Wide
                         output_dir => $self->o('output_dir'), #Wide
                         out_sam_dir => $self->o('output_sam_dir'),
                         missmatch => 6,
                         bam_file  => '/path/to/my/sorted/indexed/bam_file.bam',
                         word_length => 10,
                         saturate_threshold => 10000,
                         mask => 1,
                         percent_id => 97,
                         coverage => 90,
                         fullseq   => 1,
                         max_transcript => 1000000,
                         batch_size => 10000,
                       },
        -rc_name    => '5GB_multithread',
        -input_ids => [{input_id => 'file.astq',
                        fastq_pair => '',
                        header_file => '/path/to/header.txt',
                        paired => 1,
                        min_paired => 50,
                        min_mapped => 50,},],
      },
            {
        -logic_name => 'sam2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSam2Bam',
        -parameters => {
                         sam_dir => '/path/to/directory',
                         bamfile => '/path/to/my/SAM/file/introns.sam',
                         regex => '.sam',
                         headerfile => '/path/to/my/header/file/headers.txt',
                       },
        -rc_name    => '5GB_multithread',
        -input_ids => [{input_id => 'file.astq',
                        fastq_pair => '',
                        header_file => '/path/to/header.txt',
                        paired => 1,
                        min_paired => 50,
                        min_mapped => 50,},],
      },

      {
        -logic_name => 'refine_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRefineSolexaGenes',
        -parameters => {
               intron_bam_files => [
                    {
                     file => "/path/to/introns/bamfile",
                     mixed_bam => "0",
                     depth => "0",
                     groupname => [],
                    },
               ],
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
               single_exon_model => 'single_exon',
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
               best_score => 'best',
               # all other possible models
               other_isoforms => '',
               # max number of other models to make - blank = all
               other_num      => '10',
               # max number of other models to process - blank = all
               max_num      => '1000',
               # biotype to label bad models ( otherwise they are not written )
               bad_models     => '',
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

        -rc_name          => '2GB',
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['process_uniprot_files'],
                      },
      },

      {
        -logic_name => 'blast_rnaseq',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
        -parameters => {

			# If left blank all refined genes will be fetched
			logicname => '',

			# path to index to fetch the sequence of the blast hit to calculate % coverage
			index => '/path/to/indexed/sequences/from/the/blastdb/index',
                      },
        -rc_name => '1GB',
        -flow_into => {
                        1 => ['load_uniprot_seqs'],
                      },
      },


      {
        -logic_name => 'generate_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         uniprot_accession => 1,
                         uniprot_batch_size => $self->o('uniprot_genblast_batch_size'),
                         uniprot_table_name => $self->o('uniprot_table_name'),
                       },
        -rc_name      => 'default',
        -wait_for     => ['load_uniprot_seqs'],
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['genblast'],
                      },
      },


    );
    return \@analysis;
}

# override the default method, to force an automatic loading of the registry in all workers
#sub beekeeper_extra_cmdline_options {
#    my $self = shift;
#    return "-reg_conf ".$self->o("registry");
#}

sub resource_classes {
    my $self = shift;

  # Note that this builds resource requests off some of the variables at the top. Works off the idea
  # that the references all get put on one server and the pipe db is on another
  my $pipe_db_server = $self->default_options()->{'pipe_db_server'};
  my $dna_db_server = $self->default_options()->{'dna_db_server'};
  my $genblast_output_db_server = $self->default_options()->{'genblast_output_db_server'};
  my $exonerate_output_db_server = $self->default_options()->{'exonerate_output_db_server'};
  my $killlist_db_server = $self->default_options()->{'killlist_db_server'};

  my $default_mem = $self->default_options()->{'default_mem'};
  my $genblast_mem = $self->default_options()->{'genblast_mem'};
  my $genblast_retry_mem = $self->default_options()->{'genblast_retry_mem'};
  my $exonerate_mem = $self->default_options()->{'exonerate_mem'};
  my $exonerate_retry_mem = $self->default_options()->{'exonerate_retry_mem'};

  my $pipe_db_server_number;
  my $dna_db_server_number;
  my $genblast_output_db_server_number;
  my $exonerate_output_db_server_number;
  my $killlist_db_server_number;


  my $num_tokens = $self->default_options()->{'num_tokens'};

  unless($pipe_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "pipe_db_server: ".$pipe_db_server;
  }

  $pipe_db_server_number = $1;

    unless($dna_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "dna_db_server: ".$dna_db_server;
  }

  $dna_db_server_number = $1;

  unless($genblast_output_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "genblast_output_db_server: ".$genblast_output_db_server;
  }

  $genblast_output_db_server_number = $1;

    unless($exonerate_output_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "exonerate_output_db_server: ".$exonerate_output_db_server;
  }

  $exonerate_output_db_server_number = $1;

    unless($killlist_db_server=~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "killlist_db_server: ".$killlist_db_server;
  }

  $killlist_db_server_number = $1;

  unless($num_tokens) {
    die "num_tokens is uninitialised or zero. num_tokens needs to be present in default_options and not zero\n".
        "num_tokens: ".$num_tokens;
  }

    return {
      'default' => { LSF => '-q normal -M'.$default_mem.' -R"select[mem>'.$default_mem.'] '.
                            'rusage[mem='.$default_mem.','.
                            'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'genblast' => { LSF => '-q normal -M'.$genblast_mem.' -R"select[mem>'.$genblast_mem.'] '.
                             'rusage[mem='.$genblast_mem.','.
                             'myens_build'.$genblast_output_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'genblast_retry' => { LSF => '-q normal -M'.$genblast_retry_mem.' -R"select[mem>'.$genblast_retry_mem.'] '.
                                   'rusage[mem='.$genblast_retry_mem.','.
                                   'myens_build'.$genblast_output_db_server_number.'tok='.$num_tokens.','.
                                   'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                   'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'exonerate' => { LSF => '-q normal -M'.$exonerate_mem.' -R"select[mem>'.$exonerate_mem.'] '.
                              'rusage[mem='.$exonerate_mem.','.
                              'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'exonerate_retry' => { LSF => '-q normal -M'.$exonerate_retry_mem.' -R"select[mem>'.$exonerate_retry_mem.'] '.
                                    'rusage[mem='.$exonerate_retry_mem.','.
                                    'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    }
}

sub get_config_settings {

   # Shift in the group name (a hash that has a collection of logic name hashes and a default hash)
   # Shift in the logic name of the specific analysis
   my $self = shift;
   my $config_group = shift;
   my $config_logic_name = shift;

   # And additional hash keys will be stored in here
   my @additional_configs = @_;

   # Return a ref to the master hash for the group using the group name
   my $config_group_hash = $self->master_config_settings($config_group);
   unless(defined($config_group_hash)) {
     die "You have asked for a group name in master_config_settings that doesn't exist. Group name:\n".$config_group;
   }
   # Final hash is the hash reference that gets returned. It is important to note that the keys added have
   # priority based on the call to this subroutine, with priority from left to right. Keys assigned to
   # $config_logic_name will have most priority, then keys in any additional hashes, then keys from the
   # default hash. A default hash key will never override a $config_logic_name key
   my $final_hash;

   # Add keys from the logic name hash
   my $config_logic_name_hash = $config_group_hash->{$config_logic_name};
   unless(defined($config_logic_name_hash)) {
     die "You have asked for a logic name hash that doesn't exist in the group you specified.\n".
         "Group name:\n".$config_group."\nLogic name:\n".$config_logic_name;
   }

   $final_hash = $self->add_keys($config_logic_name_hash,$final_hash);

   # Add keys from any additional hashes passed in, keys that are already present will not be overriden
   foreach my $additional_hash (@additional_configs) {
     my $config_additional_hash = $config_group_hash->{$additional_hash};
     $final_hash = $self->add_keys($config_additional_hash,$final_hash);
   }

   # Default is always loaded and has the lowest key value priority
   my $config_default_hash = $config_group_hash->{'Default'};
   $final_hash = $self->add_keys($config_default_hash,$final_hash);

   return($final_hash);
}

sub add_keys {
  my ($self,$hash_to_add,$final_hash) = @_;

  foreach my $key (keys(%$hash_to_add)) {
    unless(exists($final_hash->{$key})) {
      $final_hash->{$key} = $hash_to_add->{$key};
    }
  }

  return($final_hash);
}

sub master_config_settings {

  my ($self,$config_group) = @_;
  my $master_config_settings = {
  exonerate_protein => {
    Default => {
                 IIDREGEXP           => '(\d+):(\d+)',
                 OPTIONS             => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1',
                 COVERAGE_BY_ALIGNED => 0,
                 QUERYTYPE           => 'protein',
                 GENOMICSEQS         => $self->o('genome_file'),
                 PROGRAM             => $self->o('exonerate_path'),
                 SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
               },

    exonerate => {
                   FILTER                        => {
                     OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                     PARAMETERS                  => {
                       -coverage                 => $self->o('exonerate_cov'),
                       -percent_id               => $self->o('exonerate_pid'),
                       -best_in_genome           => 1,
                       -reject_processed_pseudos => 1,
                     },
                   },
                 },

    killlist_protein => {
                         KILLLISTDB          => $self->o('killlist_db'),
                         USE_KILL_LIST       => 1,
                         KILL_TYPE           => 'protein',
                         KILL_LIST_FILTER    => {
                                                  -only_mol_type        => 'protein',
                                                  -user_id              => undef,
                                                  -from_source_species  => undef,
                                                  -before_date          => undef,
                                                  -having_status        => undef,
                                                  -reasons              => [],
                                                  -for_analyses         => [],
                                                  -for_species          => [],
                                                  -for_external_db_ids  => [],
                                                },
                       },
    },
  };

  return($master_config_settings->{$config_group});

}

1;
