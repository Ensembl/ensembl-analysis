=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

 Bio::EnsEMBL::Production::Pipeline::PipeConfig::DumpCore_vertebrates_conf;

=head1 DESCRIPTION

=head1 AUTHOR 

 ckong@ebi.ac.uk 

=cut
package pre_release_dumps_conf;

use strict;
use warnings;
use File::Spec;
use Data::Dumper;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
       # inherit other stuff from the base class
       %{ $self->SUPER::default_options() },

       'ftp_base'      => "/ebi/ftp/pub/databases/ensembl/vgp",

       ## General parameters
       'registry'      => $self->o('registry'),
       'release'       => $self->o('release'),
       'pipeline_name' => "ftp_pipeline",
       'email'         => $self->o('ENV', 'USER').'@ebi.ac.uk',
       'tmp_dump_dir'  => '',
       'xrefs'         => 0,

       ## 'job_factory' parameters
       'species'              => [],

       ## dump_gff3 & dump_gtf parameter
       'abinitio'        => 1,
       'gene' => 1,

       ## dump_gtf parameters, e! specific
       'gtftogenepred_exe' => 'gtfToGenePred',
       'genepredcheck_exe' => 'genePredCheck',

       ## dump_gff3 parameters
      'gt_exe'          => 'gt',
      'gff3_tidy'       => $self->o('gt_exe').' gff3 -tidy -sort -retainids -force',
      'gff3_validate'   => $self->o('gt_exe').' gff3validator',

      'feature_type'    => ['Gene', 'Transcript', 'SimpleFeature'], #'RepeatFeature'
      'per_chromosome'  => 1,
      'include_scaffold'=> 1,
      'logic_name'      => [],
      'db_type'        => 'core',
      'out_file_stem'   => undef,

      ## dump_fasta parameters
      'dna_sequence_type_list'  => ['dna'],
      # Do/Don't process these logic names
      'process_logic_names' => [],
      'skip_logic_names'    => [],
      # Previous release FASTA DNA files location
      # Previous release number
      'prev_rel_dir' => '/nfs/ensemblftp/PUBLIC/pub/',
      'previous_release' => (software_version() - 1),

      # create BLAST databases, version 2.2.27+
      'ncbiblast_exe' => 'makeblastdb',
      # convert DNA from fasta to 2bit format
      'blat_exe' => 'faToTwoBit',
      # Previous release FASTA DNA files location
      'prev_rel_dir' => '/nfs/ensemblftp/PUBLIC/pub/',
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
	    @{$self->SUPER::pipeline_create_commands},
      'mkdir -p '.$self->o('tmp_dump_dir'),
    ];
}

# Override the default method, to force an automatic loading of the registry in all workers
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  return
      ' -reg_conf ' . $self->o('registry'),
  ;
}

# these parameter values are visible to all analyses,
# can be overridden by parameters{} and input_id{}
sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
             %{$self->SUPER::pipeline_wide_parameters},  # here we inherit anything from the base class
             'pipeline_name' => $self->o('pipeline_name'), #This must be defined for the beekeeper to work properly
             'base_path'     => $self->o('tmp_dump_dir'),
             'release'       => $self->o('release'),
           };
}

sub resource_classes {
    my $self = shift;
    return {
             'default'          => {'LSF' => '-q production-rh74 -n 4 -M 4000   -R "rusage[mem=4000]"'},
             '32GB'             => {'LSF' => '-q production-rh74 -n 4 -M 32000  -R "rusage[mem=32000]"'},
             '64GB'             => {'LSF' => '-q production-rh74 -n 4 -M 64000  -R "rusage[mem=64000]"'},
             '128GB'            => {'LSF' => '-q production-rh74 -n 4 -M 128000 -R "rusage[mem=128000]"'},
             '256GB'            => {'LSF' => '-q production-rh74 -n 4 -M 256000 -R "rusage[mem=256000]"'},
           }
  }

sub pipeline_analyses {
    my ($self) = @_;

    my $pipeline_flow_1  = ['dump_gtf', 'dump_gff3', 'dump_embl', 'dump_genbank', 'dump_fasta_dna', 'dump_tsv_ena', 'dump_tsv_metadata'];
    my $pipeline_flow_2 = ['copy_sm_genome_to_ftp', 'create_readme', 'cp_dumps_to_ftp', 'fan_add_repeatmodeler'];

    return [

	    {
	      -logic_name     => 'dump_job_factory',
	      -module         => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
              -parameters     => {
				   division => 'vertebrates',
				 },
              -input_ids      => [ {} ],
	      -hive_capacity   => -1,
	      -rc_name        => 'default',
              -max_retry_count => 1,
              -flow_into       => { '2' => 'backbone_job_pipeline'},
            },

	   {
	      -logic_name     => 'backbone_job_pipeline',
              -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
              -hive_capacity  => -1,
	      -rc_name        => 'default',
              -flow_into      => {
				    '1->A' => $pipeline_flow_1,
                                    'A->1' => ['create_ftp_jobs'],
                                 }
           },

### GTF
	   {
              -logic_name     => 'dump_gtf',
              -module         => 'Bio::EnsEMBL::Production::Pipeline::GTF::DumpFile',
              -parameters     => {
                                    gtf_to_genepred => $self->o('gtftogenepred_exe'),
                                    gene_pred_check => $self->o('genepredcheck_exe'),
                                    abinitio        => $self->o('abinitio'),
                                    gene            => $self->o('gene')
                                 },
              -hive_capacity  => 50,
              -rc_name        => 'default',
              -flow_into      => { '-1' => 'dump_gtf_32GB', },
           },

	   {
              -logic_name     => 'dump_gtf_32GB',
              -module         => 'Bio::EnsEMBL::Production::Pipeline::GTF::DumpFile',
              -parameters     => {
                                    gtf_to_genepred => $self->o('gtftogenepred_exe'),
                                    gene_pred_check => $self->o('genepredcheck_exe'),
                                    abinitio        => $self->o('abinitio'),
                                    gene            => $self->o('gene')
                                 },
              -hive_capacity  => 50,
              -rc_name       => '32GB',
              -flow_into      => { '-1' => 'dump_gtf_64GB', },
           },

	   {
              -logic_name     => 'dump_gtf_64GB',
              -module         => 'Bio::EnsEMBL::Production::Pipeline::GTF::DumpFile',
              -parameters     => {
                                    gtf_to_genepred => $self->o('gtftogenepred_exe'),
                                    gene_pred_check => $self->o('genepredcheck_exe'),
                                    abinitio        => $self->o('abinitio'),
                                    gene            => $self->o('gene')
                                 },
              -hive_capacity  => 50,
              -rc_name       => '64GB',
              -flow_into      => { '-1' => 'dump_gtf_128GB', },
           },

	   {
              -logic_name     => 'dump_gtf_128GB',
              -module         => 'Bio::EnsEMBL::Production::Pipeline::GTF::DumpFile',
              -parameters     => {
                                    gtf_to_genepred => $self->o('gtftogenepred_exe'),
                                    gene_pred_check => $self->o('genepredcheck_exe'),
                                    abinitio        => $self->o('abinitio'),
                                    gene => $self->o('gene')
                                 },
              -hive_capacity  => 50,
              -rc_name       => '128GB',
           },

### GFF3
	   {
              -logic_name     => 'dump_gff3',
              -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
              -parameters     => {
                                    feature_type       => $self->o('feature_type'),
                                    per_chromosome     => $self->o('per_chromosome'),
                                    include_scaffold   => $self->o('include_scaffold'),
                                    logic_name         => $self->o('logic_name'),
                                    db_type            => $self->o('db_type'),
                                    abinitio           => $self->o('abinitio'),
                                    gene               => $self->o('gene'),
                                    out_file_stem      => $self->o('out_file_stem'),
                                    xrefs              => $self->o('xrefs'),
                                 },
              -hive_capacity  => 50,
              -rc_name        => 'default',
              -flow_into      => {
                                    '-1' => 'dump_gff3_32GB',
                                    '1'  => 'tidy_gff3',
                                 },
           },

	   {
              -logic_name     => 'dump_gff3_32GB',
              -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
              -parameters     => {
                                    feature_type       => $self->o('feature_type'),
                                    per_chromosome     => $self->o('per_chromosome'),
                                    include_scaffold   => $self->o('include_scaffold'),
                                    logic_name         => $self->o('logic_name'),
                                    db_type            => $self->o('db_type'),
                                    abinitio           => $self->o('abinitio'),
                                    gene               => $self->o('gene'),
                                    out_file_stem      => $self->o('out_file_stem'),
                                    xrefs              => $self->o('xrefs'),
                                 },
              -hive_capacity  => 50,
              -rc_name        => '32GB',
              -flow_into      => {
                                    '-1' => 'dump_gff3_64GB',
                                    '1'  => 'tidy_gff3',
                                  },
           },

	   {
              -logic_name     => 'dump_gff3_64GB',
              -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
              -parameters     => {
                                    feature_type       => $self->o('feature_type'),
                                    per_chromosome     => $self->o('per_chromosome'),
                                    include_scaffold   => $self->o('include_scaffold'),
                                    logic_name         => $self->o('logic_name'),
                                    db_type            => $self->o('db_type'),
                                    abinitio           => $self->o('abinitio'),
                                    gene               => $self->o('gene'),
                                    out_file_stem      => $self->o('out_file_stem'),
                                    xrefs              => $self->o('xrefs'),
                                 },
              -hive_capacity  => 50,
              -rc_name        => '64GB',
              -flow_into      => {
                                    '-1' => 'dump_gff3_128GB',
                                    '1'  => 'tidy_gff3',
                                 },
           },

	   {
              -logic_name     => 'dump_gff3_128GB',
              -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
              -parameters     => {
                                    feature_type       => $self->o('feature_type'),
                                    per_chromosome     => $self->o('per_chromosome'),
                                    include_scaffold   => $self->o('include_scaffold'),
                                    logic_name         => $self->o('logic_name'),
                                    db_type            => $self->o('db_type'),
                                    abinitio           => $self->o('abinitio'),
                                    gene               => $self->o('gene'),
                                    out_file_stem      => $self->o('out_file_stem'),
                                    xrefs              => $self->o('xrefs'),
                                 },
              -hive_capacity  => 50,
              -rc_name        => '128GB',
              -flow_into      => {
                                    '1'  => 'tidy_gff3',
                                 },
           },

### GFF3:post-processing
	   {
	      -logic_name     => 'tidy_gff3',
              -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters     => {
                                    cmd => $self->o('gff3_tidy').' -gzip -o #out_file#.sorted.gz #out_file#',
				 },
              -hive_capacity  => 10,
              -batch_size     => 10,
              -rc_name        => 'default',
              -flow_into      => 'move_gff3',
           },

	   {
              -logic_name     => 'move_gff3',
              -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters     => {
				    cmd => 'mv #out_file#.sorted.gz #out_file#',
				 },
              -hive_capacity  => 10,
              -rc_name        => 'default',
              -meadow_type    => 'LOCAL',
              -flow_into      => 'validate_gff3',
           },

	   {
              -logic_name     => 'validate_gff3',
              -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters     => {
			            cmd => $self->o('gff3_validate').' #out_file#',
			         },
              -hive_capacity  => 10,
              -batch_size     => 10,
              -rc_name        => 'default',
           },

### EMBL
	   {
	      -logic_name    => 'dump_embl',
              -module        => 'Bio::EnsEMBL::Production::Pipeline::Flatfile::DumpFile',
              -parameters    => {
                                    type => 'embl',
				},
              -hive_capacity => 50,
              -rc_name       => 'default',
              -flow_into     => { '-1' => 'dump_embl_32GB', },
           },

	   {
	      -logic_name    => 'dump_embl_32GB',
              -module        => 'Bio::EnsEMBL::Production::Pipeline::Flatfile::DumpFile',
              -parameters    => {
				    type => 'embl',
				},
              -hive_capacity => 50,
              -rc_name       => '32GB',
              -flow_into     => { '-1' => 'dump_embl_64GB', },
           },

	   {
	      -logic_name    => 'dump_embl_64GB',
	      -module        => 'Bio::EnsEMBL::Production::Pipeline::Flatfile::DumpFile',
              -parameters    => {
				    type => 'embl',
				},
              -hive_capacity => 50,
              -rc_name       => '64GB',
              -flow_into     => { '-1' => 'dump_embl_128GB', },
           },

	   {
	      -logic_name    => 'dump_embl_128GB',
	      -module        => 'Bio::EnsEMBL::Production::Pipeline::Flatfile::DumpFile',
	      -parameters    => {
				    type => 'embl',
				},
              -hive_capacity => 50,
              -rc_name       => '128GB',
           },

### GENBANK
	   {
	      -logic_name    => 'dump_genbank',
              -module        => 'Bio::EnsEMBL::Production::Pipeline::Flatfile::DumpFile',
              -parameters    => {
				   type => 'genbank',
				},
              -hive_capacity => 50,
              -rc_name       => 'default',
              -flow_into     => { -1 => 'dump_genbank_32GB', },
           },

	   {
	      -logic_name    => 'dump_genbank_32GB',
              -module        => 'Bio::EnsEMBL::Production::Pipeline::Flatfile::DumpFile',
              -parameters    => {
				    type => 'genbank',
				},
              -hive_capacity => 50,
              -rc_name       => '32GB',
              -flow_into     => { -1 => 'dump_genbank_64GB', },
           },

	   {
	       -logic_name    => 'dump_genbank_64GB',
               -module        => 'Bio::EnsEMBL::Production::Pipeline::Flatfile::DumpFile',
               -parameters    => {
				     type => 'genbank',
				 },
               -hive_capacity => 50,
               -rc_name       => '64GB',
               -flow_into     => { -1 => 'dump_genbank_128GB', },
           },

	   {
	       -logic_name    => 'dump_genbank_128GB',
               -module        => 'Bio::EnsEMBL::Production::Pipeline::Flatfile::DumpFile',
               -parameters    => { type      => 'genbank',},
               -hive_capacity => 50,
               -rc_name       => '128GB',
           },

### FASTA (dna)
	    {
               -logic_name  => 'dump_fasta_dna',
               -module      => 'Bio::EnsEMBL::Production::Pipeline::FASTA::DumpFile',
               -parameters  => {
                                  sequence_type_list  => $self->o('dna_sequence_type_list'),
                               },
               -max_retry_count => 1,
               -hive_capacity   => 10,
               -priority        => 5,
               -rc_name         => 'default',
           },

### TSV
	   {
	       -logic_name    => 'dump_tsv_ena',
	       -module        => 'Bio::EnsEMBL::Production::Pipeline::TSV::DumpFileEna',
	       -hive_capacity => 50,
	       -rc_name       => 'default',
	   },

	   {
	       -logic_name    => 'dump_tsv_metadata',
	       -module        => 'Bio::EnsEMBL::Production::Pipeline::TSV::DumpFileMetadata',
	       -hive_capacity => 50,
	       -rc_name       => 'default',
	   },

### SETUP FTP
	    {
               -logic_name    => 'create_ftp_jobs',
               -module        => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::FTPSetup',
               -rc_name       => 'default',
               -flow_into     => { 1 => 'create_ftp_dirs' },
           },

	   {
               -logic_name    => 'create_ftp_dirs',
	       -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	       -parameters    => {
				    cmd          => 'mkdir -p '.$self->o('ftp_base').'/#species#/#gca#/pre-release_dumps',
				 },
               -rc_name       => 'default',
	       -flow_into     => { 1 => $pipeline_flow_2 },
	   },

	   {
               -logic_name    => 'copy_sm_genome_to_ftp',
               -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
               -parameters    => {
                                    cmd => 'gzip #softmasked_genome_file# '.$self->o('ftp_base').'/#species#/#gca#/#softmasked_genome_file#.gz',
                                 },
               -rc_name       => 'default',
           },

	   {
               -logic_name    => 'create_readme',
               -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
               -parameters    => {
				    free_text    => 'README\n\nThis directory contains pre-release data, i.e. files pertain to the raw annotation from the Ensembl Genebuild team, these are available ahead of full integration in to Ensembl. This is NOT the final annotation that will appear in a future Ensembl release.\n\n-----------------------\nPRE_RELEASE DIRECTORY\n-----------------------\n\n\"pre-release_dumps/\"\n\nPre-release data for #species# can be found here. For more information, see the README found in each data directory.\n\n-----------------------\nSOFTMASKED GENOME FILE\n-----------------------\n\n\"#species#_softmasked_toplevel.fa\"\n\nThe genomic sequence has been screened for repeats using Repeatmasker. For #species#, the Repbase #repbase_library# library was used with Repeatmasker.\n\n-----------------\n',
                                    cmd          => 'cd '.$self->o('ftp_base').'/#species#/#gca#; printf "#free_text#" >README',
                                 },
               -rc_name       => 'default',
           },

	   {
               -logic_name    => 'cp_dumps_to_ftp',
               -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
               -parameters    => {
				    division     => 'vertebrates',
                                    cmd          => 'rsync -rtWv #base_path#/#division#/embl/#species#/ '.$self->o('ftp_base').'/#species#/#gca#/pre-release_dumps/embl; rsync -rtWv #base_path#/#division#/genbank/#species#/ '.$self->o('ftp_base').'/#species#/#gca#/pre-release_dumps/genbank; rsync -rtWv #base_path#/#division#/fasta/#species#/dna/ '.$self->o('ftp_base').'/#species#/#gca#/pre-release_dumps/fasta; rsync -rtWv #base_path#/#division#/gff3/#species#/ '.$self->o('ftp_base').'/#species#/#gca#/pre-release_dumps/gff3; rsync -rtWv #base_path#/#division#/gtf/#species#/ '.$self->o('ftp_base').'/#species#/#gca#/pre-release_dumps/gtf; rsync -rtWv #base_path#/#division#/tsv/#species#/ '.$self->o('ftp_base').'/#species#/#gca#/pre-release_dumps/tsv;',
                                 },
	       -flow_into     => { 1 => 'rm_abinitio_files' },
               -rc_name       => 'default',
           },

	   {
               -logic_name    => 'rm_abinitio_files',
               -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
               -parameters    => {
                                    cmd => 'rm '.$self->o('ftp_base').'/#species#/#gca#/pre-release_dumps/*/*abinitio*',
				 },
               -rc_name       => 'default',
	   },

	   {
	       -logic_name    => 'fan_add_repeatmodeler',
	       -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	       -parameters    => {
			            cmd => 'if [ ! -f "#repeatmodeler_library#" ]; then exit 42; else exit 0;fi',
                                    return_codes_2_branches => {'42' => 2},
				 },
	       -rc_name       => 'default',
	       -flow_into     => {
				    1 => ['cp_repeatmodeler_to_ftp', 'append_repeatmodeler_readme'],
				 },
          },

	  {
               -logic_name    => 'cp_repeatmodeler_to_ftp',
               -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
               -parameters    => {
                                    cmd => 'gzip #repeatmodeler_library# '.$self->o('ftp_base').'/#species#/#gca#/#gca#.repeatmodeler.fa.gz',
                                 },
               -rc_name       => 'default',
          },

	  {
               -logic_name    => 'append_repeatmodeler_readme',
               -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
               -parameters    => {
                                    free_text    => '\n-----------------------------\nCUSTOM REPEATMODELER LIBRARY \n-----------------------------\n\n#gca#.repeatmodeler.fa\n\nA custom repeat library was created using RepeatModeler.\n\n',
                                    cmd          => 'cd '.$self->o('ftp_base').'/#species#/#gca#; printf "#free_text#" >>README',
                                 },
               -rc_name       => 'default',
          },


    ];
}


1;
