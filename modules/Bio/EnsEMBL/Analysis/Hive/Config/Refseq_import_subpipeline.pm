
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::Refseq_import_subpipeline;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);

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
    'user_r'           => '',    # read only db user
    'user'             => '',    # write db user
    'password'         => '',    # password for write db user
    'pipe_db_server'   => '',    # host for pipe db
    'dna_db_server'    => '',    # host for dna db
    'pipe_db_port'     => '',    # port for pipeline host
    'dna_db_port'      => '',    # port for dna db host
    'databases_server' => '',
    'databases_port'   => '',

    'release_number' => '' || $self->o('ensembl_release'),
    'production_name' => '',     # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names

    'output_path' => '',         # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses

    'assembly_name'             => '',    # Name (as it appears in the assembly report file)
    'assembly_refseq_accession' => '',    # Versioned GCF accession, e.g. GCF_001857705.1

    ########################
    # Pipe and ref db info
    ########################

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is

    'reference_db_name'   => $self->o('dna_db_name'),
    'reference_db_server' => $self->o('dna_db_server'),
    'reference_db_port'   => $self->o('dna_db_port'),

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . '_core_' . $self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release' => $ENV{ENSEMBL_RELEASE},    # this is the current release version on staging to be able to get the correct database

    ######################################################
    #
    # Mostly constant settings
    #
    ######################################################

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

    'refseq_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_refseq_' . $self->o('release_number'),
      -host   => $self->o('refseq_db_server'),
      -port   => $self->o('refseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    ########################################################
    # URLs for retrieving the INSDC contigs and RefSeq files
    ########################################################
    'ncbi_base_ftp' => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
    'refseq_base_ftp'        => $self->o('ncbi_base_ftp') . '/#expr(substr(#assembly_refseq_accession#, 0, 3))expr#/#expr(substr(#assembly_refseq_accession#, 4, 3))expr#/#expr(substr(#assembly_refseq_accession#, 7, 3))expr#/#expr(substr(#assembly_refseq_accession#, 10, 3))expr#/#assembly_refseq_accession#_#assembly_name#',
    'refseq_import_ftp_path' => $self->o('refseq_base_ftp') . '/#assembly_refseq_accession#_#assembly_name#_genomic.gff.gz',
    'refseq_report_ftp_path' => $self->o('refseq_base_ftp') . '/#assembly_refseq_accession#_#assembly_name#_assembly_report.txt',

    ensembl_analysis_script => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    refseq_import_script_path => catfile( $self->o('ensembl_analysis_script'), 'refseq', 'parse_ncbi_gff3.pl' ),

  };
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [

    {
      -logic_name => 'load_refseq_synonyms',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadRefSeqSynonyms',
      -parameters => {
        'target_db'  => $self->o('reference_db'),
        'output_dir' => $self->o('output_path'),
        'url'        => $self->o('refseq_report_ftp_path'),
      },
      -rc_name   => 'default',
      -input_ids  => [
        {
          assembly_name => $self->o('assembly_name'),
          assembly_refseq_accession => $self->o('assembly_refseq_accession'),
        },
      ],
      -flow_into => {
        1 => ['download_refseq_gff'],
      },
    },

    {
      -logic_name => 'download_refseq_gff',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData',
      -parameters => {
        output_dir => catdir( $self->o('output_path'), 'refseq_import' ),
        url => $self->o('refseq_import_ftp_path'),
        download_method => 'ftp',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_refseq_db'],
      },
    },

    {
      -logic_name => 'create_refseq_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('refseq_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['load_refseq_gff'],
      },
    },

    {
      -logic_name => 'load_refseq_gff',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('refseq_import_script_path') .
          ' -dnahost ' . $self->o( 'dna_db', '-host' ) .
          ' -dnadbname ' . $self->o( 'dna_db', '-dbname' ) .
          ' -dnaport ' . $self->o( 'dna_db', '-port' ) .
          ' -dnauser ' . $self->o( 'dna_db', '-user' ) .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'refseq_db', '-host' ) .
          ' -port ' . $self->o( 'refseq_db', '-port' ) .
          ' -dbname ' . $self->o( 'refseq_db', '-dbname' ) .
          ' -write' .
          ' -file ' . catfile( $self->o('output_path'), 'refseq_import', '#assembly_refseq_accession#_#assembly_name#_genomic.gff' ),
      },
      -rc_name => 'refseq_import',
    },

  ];
}

sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder( 'production', 900 ) },
    'refseq_import' => { LSF => $self->lsf_resource_builder( 'production', 9900 ) },
    }
}

1;
