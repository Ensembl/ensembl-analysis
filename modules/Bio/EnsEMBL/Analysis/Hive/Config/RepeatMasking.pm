
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
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
    'user_r'                           => '',                                  # read only db user
    'user'                             => '',                                  # write db user
    'password'                         => '',                                  # password for write db user
    'server_set'                       => '',                                  # What server set to user, e.g. set1
    'pipe_db_server'                   => '',                                  # host for pipe db
    'dna_db_server'                    => '',                                  # host for dna db
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

########################
    # Pipe and ref db info
########################

# The following might not be known in advance, since the come from other pipelines
# These values can be replaced in the analysis_base table if they're not known yet
# If they are not needed (i.e. no projection or rnaseq) then leave them as is

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . '_core_' . $self->o('release_number'),

    'reference_db_name'   => $self->o('dna_db_name'),
    'reference_db_server' => $self->o('dna_db_server'),
    'reference_db_port'   => $self->o('dna_db_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_server' => 'mysql-ens-meta-prod-1',
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

    ensembl_analysis_script => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    sequence_dump_script => catfile( $self->o('ensembl_analysis_script'), 'sequence_dump.pl' ),

########################
    # Extra db settings
########################

    'num_tokens' => 10,

########################
    # Executable paths
########################

    'blast_type'        => 'ncbi',                                               # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'dust_path'         => catfile( $self->o('binary_base'), 'dustmasker' ),
    'trf_path'          => catfile( $self->o('binary_base'), 'trf' ),
    'repeatmasker_path' => catfile( $self->o('binary_base'), 'RepeatMasker' ),
    samtools_path       => catfile( $self->o('binary_base'), 'samtools' ),       #You may need to specify the full path to the samtools binary

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
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  # set the logic names for repeat masking
  my $wide_repeat_logic_names;
  if ( $self->o('use_repeatmodeler_to_mask') ) {
    $wide_repeat_logic_names = [ $self->o('full_repbase_logic_name'), $self->o('repeatmodeler_logic_name'), 'dust' ];
  } elsif ( $self->o('replace_repbase_with_red_to_mask') ) {
    $wide_repeat_logic_names = [ $self->o('red_logic_name'), 'dust' ];
  } else {
    $wide_repeat_logic_names = [ $self->o('full_repbase_logic_name'), 'dust' ];
  }

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    skip_repeatmodeler      => $self->o('skip_repeatmodeler'),
    wide_repeat_logic_names => $wide_repeat_logic_names,
    use_genome_flatfile     => $self->o('use_genome_flatfile'),
    genome_file             => $self->o('faidx_genome_file'),
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
      -flow_into => {
        '2' => ['semaphore_10mb_slices'],
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
        slice_size        => 1000000,
        batch_target_size => 500000,
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
        use_genome_flatfile => 1,
        genome_file         => $self->o('faidx_genome_file'),
      },
      -rc_name   => 'repeatmasker',
      -flow_into => {
        '-1' => ['rebatch_repeatmasker'],
        '-2' => ['rebatch_repeatmasker'],
      },
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
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
        use_genome_flatfile => 1,
        genome_file         => $self->o('faidx_genome_file'),
      },
      -rc_name   => 'repeatmasker_rebatch',
      -flow_into => {
        -1 => ['failed_repeatmasker_batches'],
        -2 => ['failed_repeatmasker_batches'],
      },
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
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
        use_genome_flatfile => 1,
        genome_file         => $self->o('faidx_genome_file'),
      },
      -rc_name   => 'repeatmasker',
      -flow_into => {
        '-1' => ['rebatch_repeatmasker_repeatmodeler'],
        '-2' => ['rebatch_repeatmasker_repeatmodeler'],
      },
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
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
        use_genome_flatfile => 1,
        genome_file         => $self->o('faidx_genome_file'),
      },
      -rc_name   => 'repeatmasker_rebatch',
      -flow_into => {
        -1 => ['failed_repeatmasker_repeatmodeler_batches'],
        -2 => ['failed_repeatmasker_repeatmodeler_batches'],
      },
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
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
      # Set the toplevel
      -logic_name => 'dump_softmasked_toplevel',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
      -parameters => {
        'coord_system_name'  => 'toplevel',
        'target_db'          => $self->o('reference_db'),
        'output_path'        => $self->o('genome_dumps'),
        'enscode_root_dir'   => $self->o('enscode_root_dir'),
        'species_name'       => $self->o('species_name'),
        'repeat_logic_names' => '#wide_repeat_logic_names#',
      },
      -input_ids => [ {} ],
      -wait_for  => ['run_dust'],
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
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
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
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      -batch_size    => 20,
    },

  ];
}

sub resource_classes {
  my $self = shift;
  return {
    '2GB'                  => { LSF => $self->lsf_resource_builder( 'production-rh74', 2000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '3GB'                  => { LSF => $self->lsf_resource_builder( 'production-rh74', 3000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '5GB'                  => { LSF => $self->lsf_resource_builder( 'production-rh74', 5000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'default'              => { LSF => $self->lsf_resource_builder( 'production-rh74', 900,  [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'repeatmasker'         => { LSF => $self->lsf_resource_builder( 'production-rh74', 2900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'repeatmasker_rebatch' => { LSF => $self->lsf_resource_builder( 'production-rh74', 5900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'simple_features'      => { LSF => $self->lsf_resource_builder( 'production-rh74', 2900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    }
}

sub hive_capacity_classes {
  my $self = shift;
  return {
    'hc_very_low' => 35,
    'hc_low'      => 200,
    'hc_medium'   => 500,
    'hc_high'     => 1000,
  };
}

1;
