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

=head1 NAME

Bam2BigWig_conf

=head1 SYNOPSIS


=head1 DESCRIPTION

The pipeline will need all your BAM files to be in the same directory.
It will create the drectory structure needed by the Production team.
It will link the BAM files in the output directory for easy management.
It will generate the BigWig files, the md5 checksum and the readme file
By default, the file names are:
BAM.bam
BAM.bam.bam
BAM.bam.bw
README.1
md5sum.txt.1

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::Bam2BigWig_conf;

use strict;
use warnings;

use File::Spec::Functions qw(catdir catfile);

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },
    'assembly_name' => '',
    'pipeline_name'                => 'bam2bigwig_'.$self->o('species_name'), # pipeline_name does not like '.' as it is used to create the database name
    'user'                         => '', #!!!!!!!!!!!
    'password'                     => '', #!!!!!!!!!!!
    'port'                         => '', #!!!!!!!!!!!
    'host' => '',
    'species_name' => '',

    'input_dir' => '', # Where your BAM files are
    'output_path' => '', # Where to create the BigWigs and link the BAM files
    'samtools' => catfile($self->o('binary_base'), 'samtools'),
    'bedtools' => catfile($self->o('binary_base'), 'bedtools'),
    'bedGraphToBigWig' => catfile($self->o('binary_base'), 'bedGraphToBigWig'),
  };
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
                         input_dir => $self->o('input_dir'),
                         output_dir => catdir('#output_path#', $self->o('species_name'), $self->o('assembly_name'), 'rnaseq'),
                         output_path => $self->o('output_path'),
    };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'create_directory',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'mkdir -p #output_dir#; which lfs &> /dev/null; if [ $? -eq 0 ]; then lfs getstripe #output_dir# &> /dev/null; if [ $? -eq 0 ];then lfs setstripe -c -1 #output_dir#;fi;fi',
        },
        -rc_name => 'default',
        -input_ids => [{species_name => $self->o('species_name')}],
        -flow_into  => {
          1 => ['create_bam_file_job'],
        },
    },

    {
      -logic_name => 'create_bam_file_job',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputcmd => 'cd #input_dir#; ls *.bam',
        column_names => ['bam_file'],
      },
      -rc_name => 'default',
      -flow_into  => {
        '2->A' => ['create_chromosome_file'],
        'A->1' => ['concat_md5_sum'],
      },
    },

    {
      -logic_name => 'create_chromosome_file',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#samtools# view -H #input_dir#/#bam_file# | grep \@SQ |cut -f2,3 | sed \'s/[SL]N://g\' > #output_dir#/#bam_file#.txt',
          samtools => $self->o('samtools'),
        },
        -rc_name => '3GB',
        -flow_into  => {
          1 => ['bam2bedgraph'],
        },
    },

    {
      -logic_name => 'bam2bedgraph',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#bedtools# genomecov -ibam #input_dir#/#bam_file# -bg -split | LC_COLLATE=C sort -k1,1 -k2,2n > #output_dir#/#bam_file#.bg',
          bedtools => $self->o('bedtools'),
        },
        -rc_name => '3GB',
        -flow_into  => {
          1 => ['bedgrap2bigwig'],
          -1 => ['bam2bedgraph_himem'],
        },
    },

    {
      -logic_name => 'bam2bedgraph_himem',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#bedtools# genomecov -ibam #input_dir#/#bam_file# -bg -split | LC_COLLATE=C sort -k1,1 -k2,2n > #output_dir#/#bam_file#.bg',
          bedtools => $self->o('bedtools'),
        },
        -rc_name => '8GB',
        -flow_into  => {
          1 => ['bedgrap2bigwig_himem'],
        },
    },

    {
      -logic_name => 'bedgrap2bigwig',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#bedGraphToBigWig# #output_dir#/#bam_file#.bg #output_dir#/#bam_file#.txt #output_dir#/#bam_file#.bw',
          bedGraphToBigWig => $self->o('bedGraphToBigWig'),
        },
        -rc_name => '3GB',
        -flow_into  => {
          1 => ['clean_bg_files'],
          -1 => ['bedgrap2bigwig_himem'],
        },
    },

    {
      -logic_name => 'bedgrap2bigwig_himem',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => '#bedGraphToBigWig# #output_dir#/#bam_file#.bg #output_dir#/#bam_file#.txt #output_dir#/#bam_file#.bw',
          bedGraphToBigWig => $self->o('bedGraphToBigWig'),
        },
        -rc_name => '8GB',
        -flow_into  => {
          1 => ['clean_bg_files'],
        },
    },

    {
      -logic_name => 'clean_bg_files',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'rm  #output_dir#/#bam_file#.bg #output_dir#/#bam_file#.txt',
        },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['link_bam'],
        },
    },

    {
      -logic_name => 'link_bam',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cd #output_dir#;ln -s #input_dir#/#bam_file# #bam_file#;ln -s #input_dir#/#bam_file#.bai #bam_file#.bai',
        },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['md5_sum'],
        },
    },

    {
      -logic_name => 'md5_sum',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cd #output_dir#;md5sum #bam_file#.bw #bam_file# #bam_file#.bai > #bam_file#.md5',
        },
        -rc_name => 'default',
    },

    {
      -logic_name => 'concat_md5_sum',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cat #output_dir#/*md5 > #output_dir#/md5sum.txt.1',
        },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['clean_concat_md5_sum'],
        },
    },

    {
      -logic_name => 'clean_concat_md5_sum',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'rm #output_dir#/*md5',
        },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['create_readme'],
        },
    },

    {
      -logic_name => 'create_readme',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cd #output_dir#;FILES=($(ls *.bam));echo "#free_text#" | sed "s/NUM/$((${#FILES[*]}-1))/g;s/ \([a-z]\)\([a-z]\+_\)/ \U\1\E\2/;s/_/ /g" > README.1; IFS=$\'\n\';echo "${FILES[*]}" >> README.1',
          free_text => '"Note\n------\n\n'.
                       'The RNASeq data for #species_name# consists of NUM individual samples and one merged set containing all NUM samples.\n\n'.
                       'All files have an index file (.bai).\n\n'.
                       'Use the md5sum.txt file to check the integrity of the downloaded files.\n\n'.
                       'Files\n-----\n"',
        },
        -rc_name => 'default',
    },
  ];
};

sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder('production-rh7', 500)},
    '3GB' => { LSF => $self->lsf_resource_builder('production-rh7', 3000)},
    '8GB' => { LSF => $self->lsf_resource_builder('production-rh7', 8000)},
    '20GB' => { LSF => $self->lsf_resource_builder('production-rh7', 20000)},
  };
}

1;

