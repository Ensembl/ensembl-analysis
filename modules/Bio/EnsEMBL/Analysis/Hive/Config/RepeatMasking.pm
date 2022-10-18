
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::RepeatMasking;

use strict;
use warnings;
use File::Spec::Functions;

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
    # Misc setup info
########################
    'dbowner' => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name' => '' || $self->o('production_name') . '_' . $self->o('ensembl_release'),
    dbname_accession          => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1
    'user_r'                           => '',                                  # read only db user
    'user'                             => '',                                  # write db user
    'password'                         => '',                                  # password for write db user
    'server_set'                       => '',                                  # What server set to user, e.g. set1
    'pipe_db_host'                     => '',                                  # host for pipe db
    'dna_db_host'                      => '',                                  # host for dna db
    'pipe_db_port'                     => '',                                  # port for pipeline host
    'dna_db_port'                      => '',                                  # port for dna db host
    'repbase_logic_name'               => '',                                  # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    'repbase_library'                  => '',                                  # repbase library name, this is the actual repeat repbase library to use, e.g. "Mus musculus"
    'release_number'                   => '' || $self->o('ensembl_release'),
    'species_name'                     => '',                                  # e.g. mus_musculus
    'output_path'                      => '',                                  # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'use_genome_flatfile'              => '1',                                 # This will read sequence where possible from a dumped flatfile instead of the core db
    'skip_repeatmodeler'               => '0',                                 # Skip using our repeatmodeler library for the species with repeatmasker, will still run standard repeatmasker
    'replace_repbase_with_red_to_mask' => '0',                                 # Setting this will replace 'full_repbase_logic_name' with 'red_logic_name' repeat features in the masking process
    'red_logic_name'                   => 'repeatdetector',                    # logic name for the Red repeat finding analysis
                                                                               # Keys for custom loading, only set/modify if that's what you're doing
    'repeatmodeler_library'            => '',                                  # This should be the path to a custom repeat library, leave blank if none exists
    'use_repeatmodeler_to_mask'        => '0',                                 # Setting this will include the repeatmodeler library in the masking process
    'skip_post_repeat_analyses'        => '0',
    first_choice_repeat => $self->o('full_repbase_logic_name'),
    second_choice_repeat => $self->o('repeatmodeler_logic_name'),
    third_choice_repeat => $self->o('red_logic_name'),
	
########################
    # Pipe and ref db info
########################

# The following might not be known in advance, since the come from other pipelines
# These values can be replaced in the analysis_base table if they're not known yet
# If they are not needed (i.e. no projection or rnaseq) then leave them as is

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_core_' . $self->o('release_number'),

    'reference_db_name'   => $self->o('dna_db_name'),
    'reference_db_host'   => $self->o('dna_db_host'),
    'reference_db_port'   => $self->o('dna_db_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_host'   => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',

######################################################
    #
    # Mostly constant settings
    #
######################################################

    genome_dumps => catdir( $self->o('output_path'), 'genome_dumps' ),

# This one is used by most analyses that run against a genome flatfile like exonerate, genblast etc. Has slice name style headers. Is softmasked
    softmasked_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_softmasked_toplevel.fa' ),

# This one is used in replacement of the dna table in the core db, so where analyses override slice->seq. Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_toplevel.fa' ),

# This one is a cross between the two above, it has the seq_region name header but is softmasked. It is used by things that would both want to skip using the dna table and also want to avoid the repeat_feature table, e.g. bam2introns
    faidx_softmasked_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_softmasked_toplevel.fa.reheader' ),

    full_repbase_logic_name => "repeatmask_repbase_" . $self->o('repbase_logic_name'),

    'min_toplevel_slice_length' => 250,

    'repeatmodeler_logic_name' => 'repeatmask_repeatmodeler',

    red_msk => catfile($self->o('genome_dumps'), $self->o('species_name').'_red_msk/'),
    red_rpt => catfile($self->o('genome_dumps'), $self->o('species_name').'_red_rpt/'),

    ensembl_analysis_script => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    sequence_dump_script => catfile( $self->o('ensembl_analysis_script'), 'sequence_dump.pl' ),

########################
    # Executable paths
########################

    'blast_type'        => 'ncbi',                                               # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'dust_path'         => catfile( $self->o('binary_base'), 'dustmasker' ),
    'trf_path'          => catfile( $self->o('binary_base'), 'trf' ),
    'repeatmasker_path' => catfile( $self->o('binary_base'), 'RepeatMasker' ),
    red_path            => catfile($self->o('binary_base'), 'Red'),
    samtools_path       => catfile( $self->o('binary_base'), 'samtools' ),       #You may need to specify the full path to the samtools binary
    'eponine_java_path' => catfile($self->o('binary_base'), 'java'),
    'eponine_jar_path' => catfile($self->o('linuxbrew_home_path'), 'opt', 'eponine', 'libexec', 'eponine-scan.jar'),
    'cpg_path' => catfile($self->o('binary_base'), 'cpg_lh'),
    'trnascan_path' => catfile($self->o('binary_base'), 'tRNAscan-SE'),
	
########################
    # Misc setup info
########################
    'repeatmasker_engine' => 'crossmatch',
    'masking_timer_long'  => '5h',
    'masking_timer_short' => '2h',

########################
    # db info
########################
    'reference_db' => {
      -dbname => $self->o('reference_db_name'),
      -host   => $self->o('reference_db_host'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    skip_repeatmodeler      => $self->o('skip_repeatmodeler'),
    use_genome_flatfile     => $self->o('use_genome_flatfile'),
    genome_file             => $self->o('faidx_genome_file'),
    skip_post_repeat_analyses => $self->o('skip_post_repeat_analyses'),	
    repeatmasker_slice_size   => $self->o('repeatmasker_slice_size'),
    batch_target_size       => $self->o('batch_target_size'),
  }
}

sub pipeline_create_commands {
  my ($self) = @_;

  return [

    # inheriting database and hive tables' creation
    @{ $self->SUPER::pipeline_create_commands },

  ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [
###############################################################################
    #
    # REPEATMASKER ANALYSES
    #
###############################################################################
    {
      # Create 10mb toplevel slices, these will be split further for repeatmasker
      -logic_name => 'create_10mb_slice_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('dna_db'),
        coord_system_name     => 'toplevel',
        iid_type              => 'slice',
        slice_size            => 10000000,
        include_non_reference => 0,
        top_level             => 1,
        min_slice_length      => $self->o('min_toplevel_slice_length'),
        batch_slice_ids       => 1,
        batch_target_size     => 10000000,
      },
      -rc_name   => '2GB',
      -input_ids  => [{}],
      -flow_into => {
        '2->A' => ['semaphore_10mb_slices'],
        'A->1' => ['insert_fixed_repeat_analysis_meta_key_jobs'],
        1 => ['repeatdetector'],
      },
    },

    {
# Wait for repeatmasker to complete all the sub slices for a 10mb slice before passing to dust
      -logic_name => 'semaphore_10mb_slices',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {},
      -flow_into  => {
        '1->A' => ['create_repeatmasker_slices'],
        'A->1' => ['run_dust'],
      },
    },

    {
      -logic_name => 'create_repeatmasker_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db         => $self->o('dna_db'),
        iid_type          => 'rebatch_and_resize_slices',
        slice_size        => $self->o('repeatmasker_slice_size'),
        batch_target_size => $self->o('batch_target_size'),
      },
      -rc_name   => 'default',
      -flow_into => {
        '2' => [ 'run_repeatmasker', 'fan_repeatmodeler' ],
      },
    },

    {
      -logic_name => 'run_repeatmasker',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
      -parameters => {
        timer_batch         => $self->o('masking_timer_long'),
        target_db           => $self->o('reference_db'),
        logic_name          => $self->o('full_repbase_logic_name'),
        module              => 'HiveRepeatMasker',
        repeatmasker_path   => $self->o('repeatmasker_path'),
        commandline_params  => '-nolow -species "' . $self->o('repbase_library') . '" -engine "' . $self->o('repeatmasker_engine') . '"',
        use_genome_flatfile => $self->o('use_genome_flatfile'),
        genome_file         => $self->o('faidx_genome_file'),
        disconnect_jobs     => 1,
      },
      -rc_name   => 'repeatmasker',
      -flow_into => {
        '-1' => ['rebatch_repeatmasker'],
        '-2' => ['rebatch_repeatmasker'],
      },
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'rebatch_repeatmasker',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db         => $self->o('dna_db'),
        iid_type          => 'rebatch_and_resize_slices',
        slice_size        => 100000,
        batch_target_size => 10000,
      },
      -rc_name   => 'default',
      -flow_into => {
        '2' => ['run_repeatmasker_small_batch'],
      },
      -can_be_empty => 1,
    },

    {
      -logic_name => 'run_repeatmasker_small_batch',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
      -parameters => {
        timer_batch         => $self->o('masking_timer_short'),
        target_db           => $self->o('reference_db'),
        logic_name          => $self->o('full_repbase_logic_name'),
        module              => 'HiveRepeatMasker',
        repeatmasker_path   => $self->o('repeatmasker_path'),
        commandline_params  => '-nolow -species "' . $self->o('repbase_library') . '" -engine "' . $self->o('repeatmasker_engine') . '"',
        use_genome_flatfile => $self->o('use_genome_flatfile'),
        genome_file         => $self->o('faidx_genome_file'),
        disconnect_jobs     => 1,
      },
      -rc_name   => 'repeatmasker_rebatch',
      -flow_into => {
        -1 => ['failed_repeatmasker_batches'],
        -2 => ['failed_repeatmasker_batches'],
      },
      -hive_capacity => $self->o('hc_normal'),
      -can_be_empty  => 1,
    },

    {
      -logic_name => 'failed_repeatmasker_batches',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
      },
      -rc_name      => 'default',
      -can_be_empty => 1,
    },

    {
      -logic_name => 'fan_repeatmodeler',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ -n "' . $self->o('repeatmodeler_library') . '" ] && [ "#skip_repeatmodeler#" == "0" ]; then exit 0; else exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -rc_name => 'default',
      -flow_into => { '1' => ['run_repeatmasker_repeatmodeler'] },
    },

    {
      -logic_name => 'run_repeatmasker_repeatmodeler',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
      -parameters => {
        timer_batch         => $self->o('masking_timer_long'),
        target_db           => $self->o('reference_db'),
        logic_name          => 'repeatmask_repeatmodeler',
        module              => 'HiveRepeatMasker',
        repeatmasker_path   => $self->o('repeatmasker_path'),
        commandline_params  => '-nolow -lib "' . $self->o('repeatmodeler_library') . '" -engine "' . $self->o('repeatmasker_engine') . '"',
        use_genome_flatfile => $self->o('use_genome_flatfile'),
        genome_file         => $self->o('faidx_genome_file'),
      },
      -rc_name   => 'repeatmasker',
      -flow_into => {
        '-1' => ['rebatch_repeatmasker_repeatmodeler'],
        '-2' => ['rebatch_repeatmasker_repeatmodeler'],
      },
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'rebatch_repeatmasker_repeatmodeler',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db         => $self->o('dna_db'),
        iid_type          => 'rebatch_and_resize_slices',
        slice_size        => 100000,
        batch_target_size => 10000,
      },
      -rc_name   => 'default',
      -flow_into => {
        '2' => ['run_repeatmasker_repeatmodeler_small_batch'],
      },
      -can_be_empty => 1,
    },

    {
      -logic_name => 'run_repeatmasker_repeatmodeler_small_batch',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
      -parameters => {
        timer_batch         => $self->o('masking_timer_short'),
        target_db           => $self->o('reference_db'),
        logic_name          => $self->o('repeatmodeler_logic_name'),
        module              => 'HiveRepeatMasker',
        repeatmasker_path   => $self->o('repeatmasker_path'),
        commandline_params  => '-nolow -lib "' . $self->o('repeatmodeler_library') . '" -engine "' . $self->o('repeatmasker_engine') . '"',
        use_genome_flatfile => $self->o('use_genome_flatfile'),
        genome_file         => $self->o('faidx_genome_file'),
      },
      -rc_name   => 'repeatmasker_rebatch',
      -flow_into => {
        -1 => ['failed_repeatmasker_repeatmodeler_batches'],
        -2 => ['failed_repeatmasker_repeatmodeler_batches'],
      },
      -hive_capacity => $self->o('hc_normal'),
      -can_be_empty  => 1,
    },

    {
      -logic_name => 'failed_repeatmasker_repeatmodeler_batches',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
      },
      -rc_name      => 'default',
      -can_be_empty => 1,
    },

    {
      -logic_name => 'insert_fixed_repeat_analysis_meta_key_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name    => 'default',
      -parameters => {
        inputlist => ['dust', 'trf'],
        column_names => ['repeat_logic_name'],
      },
      -flow_into => {
        '2->A' => ['check_fixed_repeat_analysis_meta_key'],
        '1->A' => ['check_first_choice_repeat_present'],
        'A->1' => ['dump_softmasked_toplevel'],
      },
    },

    {
      -logic_name => 'check_fixed_repeat_analysis_meta_key',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => 'default',
      -parameters => {
        cmd => q(if [[ `mysql -h #reference_db_host# -P #reference_db_port# -u #user_r# #reference_db_name# -NB -e "SELECT COUNT(*) FROM repeat_feature JOIN analysis USING(analysis_id) WHERE logic_name = '#repeat_logic_name#'"` -eq 0 ]]; then exit 42; fi),
        reference_db_host => $self->o('reference_db_host'),
        reference_db_port => $self->o('reference_db_port'),
        reference_db_name => $self->o('reference_db_name'),
        user_r => $self->o('user_r'),
        return_codes_2_branches => {42 => 2},
      },
      -flow_into => {
        1 => ['insert_repeat_analysis_meta_key'],
      },
    },

    {
      -logic_name => 'check_first_choice_repeat_present',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => 'default',
      -parameters => {
        cmd => q(if [[ `mysql -h #reference_db_host# -P #reference_db_port# -u #user_r# #reference_db_name# -NB -e "SELECT COUNT(*) FROM repeat_feature JOIN analysis USING(analysis_id) WHERE logic_name = '#repeat_logic_name#'"` -eq 0 ]]; then exit 42; fi),
        repeat_logic_name => $self->o('first_choice_repeat'),
        reference_db_host => $self->o('reference_db_host'),
        reference_db_port => $self->o('reference_db_port'),
        reference_db_name => $self->o('reference_db_name'),
        user_r => $self->o('user_r'),
        return_codes_2_branches => {42 => 2},
      },
      -flow_into => {
        1 => {'insert_repeat_analysis_meta_key' => {repeat_logic_name => '#repeat_logic_name#'}},
        2 => ['check_second_choice_repeat_present'],
      },
    },

    {
      -logic_name => 'check_second_choice_repeat_present',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => 'default',
      -parameters => {
        cmd => q(if [[ `mysql -h #reference_db_host# -P #reference_db_port# -u #user_r# #reference_db_name# -NB -e "SELECT COUNT(*) FROM repeat_feature JOIN analysis USING(analysis_id) WHERE logic_name = '#repeat_logic_name#'"` -eq 0 ]]; then exit 42; fi),
        repeat_logic_name => $self->o('second_choice_repeat'),
        reference_db_host => $self->o('reference_db_host'),
        reference_db_port => $self->o('reference_db_port'),
        reference_db_name => $self->o('reference_db_name'),
        user_r => $self->o('user_r'),
        return_codes_2_branches => {42 => 2},
      },
      -flow_into => {
        1 => {'insert_repeat_analysis_meta_key' => {repeat_logic_name => '#repeat_logic_name#'}},
        2 => ['check_third_choice_repeat_present'],
      },
    },

    {
      -logic_name => 'check_third_choice_repeat_present',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => 'default',
      -parameters => {
        cmd => q(if [[ `mysql -h #reference_db_host# -P #reference_db_port# -u #user_r# #reference_db_name# -NB -e "SELECT COUNT(*) FROM repeat_feature JOIN analysis USING(analysis_id) WHERE logic_name = '#repeat_logic_name#'"` -eq 0 ]]; then exit 42; fi),
        repeat_logic_name => $self->o('third_choice_repeat'),
        reference_db_host => $self->o('reference_db_host'),
        reference_db_port => $self->o('reference_db_port'),
        reference_db_name => $self->o('reference_db_name'),
        user_r => $self->o('user_r'),
        return_codes_2_branches => {42 => 2},
      },
      -flow_into => {
        1 => {'insert_repeat_analysis_meta_key' => {repeat_logic_name => '#repeat_logic_name#'}},
      },
    },

    {
      -logic_name => 'insert_repeat_analysis_meta_key',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql => [
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"repeat.analysis","#repeat_logic_name#")',
        ],
      },
      -rc_name    => 'default',
    },

    {
      # Set the toplevel
      -logic_name => 'dump_softmasked_toplevel',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
      -parameters => {
        'coord_system_name'  => 'toplevel',
        'target_db'          => $self->o('reference_db'),
        'output_path'        => $self->o('genome_dumps'),
        'enscode_root_dir'   => $self->o('enscode_root_dir'),
        'species_name'       => $self->o('species_name'),
        mask                 => 1,
      },
      -flow_into => {
        1 => ['format_softmasked_toplevel'],
      },
      -rc_name => '3GB',
    },

    {
      # This should probably be a proper module
      -logic_name => 'format_softmasked_toplevel',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        'cmd' => 'if [ "' . $self->o('blast_type') . '" = "ncbi" ]; then convert2blastmask -in ' . $self->o('softmasked_genome_file') . ' -parse_seqids -masking_algorithm repeatmasker -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out ' . $self->o('softmasked_genome_file') . '.asnb;makeblastdb -in ' . $self->o('softmasked_genome_file') . ' -dbtype nucl  -max_file_sz "10GB"  -parse_seqids -mask_data ' . $self->o('softmasked_genome_file') . '.asnb -title "' . $self->o('species_name') . '"; else xdformat -n ' . $self->o('softmasked_genome_file') . ';fi',
      },
      -rc_name   => '5GB',
      -flow_into => {
        1 => ['create_reheadered_softmasked_file'],
      },
    },

    {
      -logic_name => 'create_reheadered_softmasked_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'convert_genome_dump.pl' ) .
          ' -conversion_type slice_name_to_seq_region_name' .
          ' -input_file ' . $self->o('softmasked_genome_file') .
          ' -output_file ' . $self->o('faidx_softmasked_genome_file'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_softmasked_faidx'],
      },
    },

    {
      -logic_name => 'create_softmasked_faidx',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '5GB',
      -parameters => {
        cmd => 'if [ ! -e "' . $self->o('faidx_softmasked_genome_file') . '.fai" ]; then ' . $self->o('samtools_path') . ' faidx ' . $self->o('faidx_softmasked_genome_file') . ';fi',
      },
    },

###############################################################################
    #
    # SIMPLE FEATURE AND OTHER REPEAT ANALYSES
    #
###############################################################################
    {
      # Run dust
      -logic_name => 'run_dust',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDust',
      -parameters => {
        target_db  => $self->o('reference_db'),
        logic_name => 'dust',
        module     => 'HiveDust',
        dust_path  => $self->o('dust_path'),
      },
      -rc_name   => 'simple_features',
      -flow_into => {
        1  => ['run_trf'],
        -1 => ['run_trf'],
        -2 => ['run_trf'],
      },
      -hive_capacity => $self->o('hc_normal'),
      -batch_size    => 20,
    },

    {
      # Run TRF
      -logic_name => 'run_trf',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRF',
      -parameters => {
        target_db  => $self->o('reference_db'),
        logic_name => 'trf',
        module     => 'HiveTRF',
        trf_path   => $self->o('trf_path'),
      },
      -rc_name       => 'simple_features',
      -flow_into => {
	  1 => ['fan_post_repeat_analyses'],
	  -1 => ['fan_post_repeat_analyses'],
	  -2 => ['fan_post_repeat_analyses'],
      },
      -hive_capacity => $self->o('hc_normal'),
      -batch_size    => 20,
    },

    {
      # This will skip downstream analyses like cpg, eponine, genscan etc. if the flag is set
      -logic_name => 'fan_post_repeat_analyses',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'if [ "#skip_post_repeat_analyses#" -ne "0" ]; then exit 42; else exit 0;fi',
	  return_codes_2_branches => {'42' => 2},
      },
      -rc_name    => 'default',
      -flow_into  => { '1' => ['run_eponine'] },
    },

    {
      # Run eponine
      -logic_name => 'run_eponine',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveEponine',
      -parameters => {
	  target_db => $self->o('reference_db'),
	  logic_name => 'eponine',
	  module => 'HiveEponine',
	  eponine_path => $self->o('eponine_java_path'),
	  commandline_params => '-epojar => '.$self->o('eponine_jar_path').', -threshold => 0.999',
      },
      -rc_name    => 'simple_features',
      -flow_into => {
	  1 => ['run_cpg'],
	  -1 => ['run_cpg'],
	  -2 => ['run_cpg'],
      },
      -hive_capacity => $self->o('hc_normal'),
      -batch_size => 20,
    },

    {
      # Run CPG
      -logic_name => 'run_cpg',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveCPG',
      -parameters => {
	  target_db => $self->o('reference_db'),
	  logic_name => 'cpg',
	  module => 'HiveCPG',
	  cpg_path => $self->o('cpg_path'),
      },
      -rc_name    => 'simple_features',
      -flow_into => {
	  1 => ['run_trnascan'],
	  -1 => ['run_trnascan'],
	  -2 => ['run_trnascan'],
      },
      -hive_capacity => $self->o('hc_normal'),
      -batch_size => 20,
    },

    {
      # Run tRNAscan
      -logic_name => 'run_trnascan',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRNAScan',
      -parameters => {
	  target_db => $self->o('reference_db'),
	  logic_name => 'trnascan',
	  module => 'HiveTRNAScan',
	  trnascan_path => $self->o('trnascan_path'),
      },
      -rc_name    => 'simple_features',
      -hive_capacity => $self->o('hc_normal'),
      -batch_size => 20,
    },
      
# Run Red (REpeat Detector)
    {
      -logic_name => 'repeatdetector',
      -module     => 'Repeatmask_Red',
      -language   => 'python3',
      -parameters => {
        logic_name     => $self->o('red_logic_name'),
        red_path       => $self->o('red_path'),
        genome_file    => $self->o('faidx_genome_file'),
        target_db_url  => $self->o('hive_driver').'://'.$self->o('user').':'.$self->o('password').'@'.$self->o('dna_db_host').':'.$self->o('dna_db_port').'/'.$self->o('dna_db_name'),
        msk            => $self->o('red_msk'),
        rpt            => $self->o('red_rpt'),
        red_meta_key   => $self->o('replace_repbase_with_red_to_mask'),
      },
      -rc_name   => '15GB',
      -flow_into => {
        -1  => ['repeatdetector_50GB'],
      },
    },

    {
      -logic_name => 'repeatdetector_50GB',
      -module     => 'Repeatmask_Red',
      -language   => 'python3',
      -parameters => {
        logic_name     => $self->o('red_logic_name'),
        red_path       => $self->o('red_path'),
        genome_file    => $self->o('faidx_genome_file'),
        target_db_url  => $self->o('hive_driver').'://'.$self->o('user').':'.$self->o('password').'@'.$self->o('dna_db_host').':'.$self->o('dna_db_port').'/'.$self->o('dna_db_name'),
        msk            => $self->o('red_msk'),
        rpt            => $self->o('red_rpt'),
        red_meta_key   => $self->o('replace_repbase_with_red_to_mask'),
      },
      -rc_name => '50GB',
    },
  ];
}

sub resource_classes {
  my $self = shift;
  return {
    '2GB'                  => { LSF => $self->lsf_resource_builder( 'production', 2000 ) },
    '3GB'                  => { LSF => $self->lsf_resource_builder( 'production', 3000 ) },
    '5GB'                  => { LSF => $self->lsf_resource_builder( 'production', 5000 ) },
    '15GB'                 => { LSF => $self->lsf_resource_builder( 'production', 15000 ) },
    '50GB'                 => { LSF => $self->lsf_resource_builder( 'production', 50000 ) },
    'default'              => { LSF => $self->lsf_resource_builder( 'production', 900 ) },
    'repeatmasker'         => { LSF => $self->lsf_resource_builder( 'production', 2900 ) },
    'repeatmasker_rebatch' => { LSF => $self->lsf_resource_builder( 'production', 5900 ) },
    'simple_features'      => { LSF => $self->lsf_resource_builder( 'production', 2900 ) },
    }
}

1;
