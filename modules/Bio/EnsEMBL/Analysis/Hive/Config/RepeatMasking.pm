
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

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::RepeatMasking;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::Analysis::Tools::SoftwareConfigLoad qw(get_software_path); #Software path config module
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;

  ## Build software path based on new software type variable
    my $software_type = $ENV{SOFTWARE_TYPE};
    my $dust_path =  get_software_path($software_type, 'dust');
    my $trf_path =  get_software_path($software_type, 'trf');
    my $repeatmasker_path = get_software_path($software_type, 'repeatmasker');
    my $red_path =  get_software_path($software_type, 'red');
    my $samtools_path =  get_software_path($software_type, 'samtools');
    my $eponine_java_path =  get_software_path($software_type, 'java');
    my $eponine_jar_path =  get_software_path($software_type, 'eponine_jar');
    my $cpg_path =  get_software_path($software_type, 'cpg');
    my $trnascan_path =  get_software_path($software_type, 'trnascan_path');
    my $blast =  get_software_path($software_type, 'blast');


  return {

    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

######################################################
    #
    # Variable settings- You change these!!!
    #
######################################################
    software_type             => $software_type,

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
    'skip_repeatmodeler'               => '1',                                 # Skip using our repeatmodeler library for the species with repeatmasker, will still run standard repeatmasker
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
	  dust_path                 => $dust_path,
    trf_path                  => $trf_path,
    repeatmasker_path         => $repeatmasker_path,
    red_path                  => $red_path,
    samtools_path             => $samtools_path,
    eponine_java_path         => $eponine_java_path,
    eponine_jar_path          => $eponine_jar_path,
    cpg_path                  => $cpg_path,
    trnascan_path             => $trnascan_path,
    blast                     => $blast,
########################
    # Misc setup info
########################
    'repeatmasker_engine' => 'rmblast',
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
      -input_ids  => [{}],
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
      -rc_name   => '3GB',
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
      -rc_name   => '3GB',
      -flow_into => {
        1 => ['create_softmasked_faidx'],
      },
    },

    {
      -logic_name => 'create_softmasked_faidx',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '3GB',
      -parameters => {
        cmd => 'if [ ! -e "' . $self->o('faidx_softmasked_genome_file') . '.fai" ]; then ' . $self->o('samtools_path') . ' faidx ' . $self->o('faidx_softmasked_genome_file') . ';fi',
      },
    },
  ];
}

sub resource_classes {
  my $self = shift;
  return {
    #inherit other stuff from the base class
    %{ $self->SUPER::resource_classes() },
  }
}

1;
