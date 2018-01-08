=head1 LICENSE

Copyright [2017-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic

=head1 SYNOPSIS

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');


=head1 DESCRIPTION

This is the config file for all pipeline sanity checks. It can be used for checking
feature counts versus expected values for a clade

=head1 METHODS

  _master_config_settings: contains all possible parameters

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key, $species) = @_;

  my %config = (
  'default' => {},
  'genome_preparation_checks' => {
     'primates_basic' => {
       # repeats
       'dust'                        => [2000000,'repeat'],
       'repeatmask_repbase_primates' => [3000000,'repeat'],
       'trf'                         => [500000,'repeat'],

       # simple features
       'cpg'                         => [15000,'simple'],
       'trnascan'                    => [300,'simple'],
       'eponine'                     => [30000,'simple'],

       # prediction transcripts
       'genscan'                     => [40000,'prediction transcript'],

       # dna align features
       'unigene'                     => [3000000,'dna align'],
       'vertrna'                     => [3000000,'dna align'],

       # protein align features
       'uniprot'                     => [3000000,'protein align'],
     },
  },

  'gene_db_checks' => {
    'primates_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 125000,
          'genblast_not_best'   => 125000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 20000,
          'primates_pe12_'      => 50000,
          'mammals_pe12_'       => 100000,
        }, # biotypes
      }, # genblast
      'ig_tr' => {
        'logic_names' => {
          'ig_tr_gene'          => 40,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 20,
          'TR_'                 => 20,
        }, # biotypes
      }, # ig_tr
      'projection_coding' => {
        'logic_names' => {
          'project_transcripts' => 30000,
        }, # logic_names
        'biotypes' =>    {
          'projection'          => 30000,
        }, # biotypes
      }, # projection_coding
      'projection_lincrna' => {
        'logic_names' => {
          'project_lincrna' => 2000,
        }, # logic_names
      }, # projection_lincrna
      'projection_pseudogene' => {
        'logic_names' => {
          'project_pseudogene' => 2000,
        }, # logic_names
      }, # projection_pseudogene
      'projection_ig_tr' => {
        'logic_names' => {
          'project_ig_tr' => 50,
        }, # logic_names
      }, # projection_ig_tr
      'realign' => {
        'logic_names' => {
          # Would actually prefer an upper limit on realign as opposed to a lower limit
          'project_transcripts'  => 20000,
          'genblast'             => 1000,
        }, # logic_names
        'biotypes' =>    {
          'realign'             => 20000,
        }, # biotypes
      }, # realign
      'rnaseq_blast' =>  {
        # This one is an issue, logic names, counts are varied and one biotype is a
        # substring of the other. At the moment it's really just a check that theres'
        # some stuff in there
        'logic_names' => {
#          $species.'_merged_rnaseq' => 10000,
        }, # logic_names
        'biotypes' =>    {
          'rnaseq'              => 10000,
         }, # biotypes
      }, # rnaseq_blast
      'layer' => {
        'biotypes' =>    {
          'IG_'                  => 20,
          'TR_'                  => 20,
          'human_pe12_'          => 10000,
          'primates_pe12_'       => 15000,
          'mammals_pe12_'        => 10000,
          'realign_'             => 30000,
          'rnaseq_merged_'       => 5000,
          'rnaseq_tissue_'       => 100,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
          'ensembl'             => 19000,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                  => 20,
          'TR_'                  => 20,
          'human_pe12_'         => 2000,
          'primates_pe12_'      => 2000,
          'mammals_pe12_'       => 1000,
          'realign_'            => 10000,
          'rnaseq_merged_'      => 2000,
        }, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 3000,
        }, # logic_names
        'biotypes' =>    {
          'miRNA'               => 500,
          'misc_RNA'            => 1000,
          'ribozyme'            => 0,
          'rRNA'                => 200,
          'scaRNA'              => 0,
          'snoRNA'              => 200,
          'snRNA'               => 400,
        }, # biotypes
      }, # ncrna
      'final' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 3000,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                  => 20,
          'TR_'                  => 20,
          'realign_'            => 10000,
          'rnaseq_merged_'      => 2000,
          'human_pe12_'         => 2000,
          'primates_pe12_'      => 2000,
          'mammals_pe12_'       => 1000,
          'miRNA'               => 500,
          'misc_RNA'            => 1000,
          'ribozyme'            => 0,
          'rRNA'                => 200,
          'scaRNA'              => 0,
          'snoRNA'              => 200,
          'snRNA'               => 400,
        }, # biotypes
      }, # final
      'core' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 3000,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                  => 20,
          'TR_'                  => 20,
          'protein_coding'       => 25000,
          'pseudogene'           => 50,
          'processed_pseudogene' => 50,
          'miRNA'                => 500,
          'misc_RNA'             => 1000,
          'ribozyme'             => 0,
          'rRNA'                 => 200,
          'scaRNA'               => 0,
          'snoRNA'               => 200,
          'snRNA'                => 400,
        }, # biotypes
      }, # core
    }, # primates_basic
    otherfeatures => {
      'logic_names' => {
        'cdna_alignment'         => 20000,
        'refseq_import'          => 30000,
      },
      'biotypes' => {
        'cdna'                   => 20000,
      },
    }, # otherfeatures
    rnaseq_final => {
      'biotypes' => {
        'protein_coding'         => 5000,
      },
    }, # rnaseq_final
  }, # gene_db_checks

  'final_core_checks' => {
     'primates_basic' => {
       # Fill in
     }, # primates_basic
  }, # final_core_checks

  );
  return $config{$key};
}

1;
