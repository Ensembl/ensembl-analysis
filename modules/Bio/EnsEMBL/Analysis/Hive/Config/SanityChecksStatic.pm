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
     'mammals_basic' => {
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

     'mammals_basic' => {
       # repeats
       'dust'                        => [1392556,'repeat'],
       'repeatmask_repbase_mammals' =>  [1765829,'repeat'],
       'trf'                         => [248970,'repeat'],

       # simple features
       'cpg'                         => [9846,'simple'],
       'trnascan'                    => [328,'simple'],
       'eponine'                     => [14000,'simple'],

       # prediction transcripts
       'genscan'                     => [43449,'prediction transcript'],

       # dna align features
       'unigene'                     => [3306419,'dna align'],
       'vertrna'                     => [4000000,'dna align'],

       # protein align features
       'uniprot'                     => [3005772,'protein align'],
     },
  },

  'gene_db_checks' => {
    'mammals_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 125000,
          'genblast_not_best'   => 125000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 20000,
          'mouse_pe12_'      => 50000,
          'vert_pe12_'      => 50000,
          'mammals_pe12_'       => 100000,
        }, # biotypes
      }, # genblast
      'ig_tr' => {
        'logic_names' => {
          'ig_tr_gene'          => 80,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 100,
          'TR_'                 => 5,
        }, # biotypes
      }, # ig_tr
      'projection_coding' => {
        'logic_names' => {
          'project_transcripts' => 51000,
        }, # logic_names
        'biotypes' =>    {
          'projection'          => 51000,
        }, # biotypes
      }, # projection_coding
      'projection_lincrna' => {
        'logic_names' => {
          'project_lincrna' => 40,
        }, # logic_names
      }, # projection_lincrna
      'projection_pseudogene' => {
        'logic_names' => {
          'project_pseudogene' => 3900,
        }, # logic_names
      }, # projection_pseudogene
      'projection_ig_tr' => {
        'logic_names' => {
          'project_ig_tr' => 80,
        }, # logic_names
      }, # projection_ig_tr
      'realign' => {
        'logic_names' => {
          # Would actually prefer an upper limit on realign as opposed to a lower limit
          'project_transcripts'  => 40000,
          'genblast'             => 9000,
        }, # logic_names
        'biotypes' =>    {
          'realign'             => 50000,
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
          'mouse_pe12_'       => 15000,
          'vert_pe12_'       => 15000,
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
          'mouse_pe12_'      => 2000,
          'vert_pe12_'      => 2000,
          'mammals_pe12_'       => 1000,
          'realign_'            => 10000,
          'rnaseq_merged_'      => 2000,
        }, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 2500,
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
          'mouse_pe12_'      => 2000,
          'vert_pe12_'      => 2000,
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
    }, # mammals_basic
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
          'ig_tr_gene'          => 80,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 100,
          'TR_'                 => 5,
        }, # biotypes
      }, # ig_tr
      'projection_coding' => {
        'logic_names' => {
          'project_transcripts' => 51000,
        }, # logic_names
        'biotypes' =>    {
          'projection'          => 51000,
        }, # biotypes
      }, # projection_coding
      'projection_lincrna' => {
        'logic_names' => {
          'project_lincrna' => 40,
        }, # logic_names
      }, # projection_lincrna
      'projection_pseudogene' => {
        'logic_names' => {
          'project_pseudogene' => 3900,
        }, # logic_names
      }, # projection_pseudogene
      'projection_ig_tr' => {
        'logic_names' => {
          'project_ig_tr' => 80,
        }, # logic_names
      }, # projection_ig_tr
      'realign' => {
        'logic_names' => {
          # Would actually prefer an upper limit on realign as opposed to a lower limit
          'project_transcripts'  => 40000,
          'genblast'             => 9000,
        }, # logic_names
        'biotypes' =>    {
          'realign'             => 50000,
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
    'mammals_basic' => {
	'genblast' => {
	    'logic_names' => {
          'genblast'            => 190000,
          'genblast_not_best'   => 103000,
	    }, # logic_names
	    'biotypes' =>    {
          'human_pe12_'         => 40000,
          'mouse_pe12_'      => 60000,
          'mammals_pe12_'       => 210000,
	  'vert_pe12_'       => 145000,
	    }, # biotypes
	}, # genblast
	'ig_tr' => {
	    'logic_names' => {
          'ig_tr_gene'          => 200,
	    }, # logic_names
	    'biotypes' =>    {
          'IG_'                 => 800,
          'TR_'                 => 40,
	    }, # biotypes
	}, # ig_tr
	'ncrna' => {
	    'logic_names' => {
          'ncrna' => 2800,
	    }, # logic_names
	    'biotypes' =>    {
          'miRNA'               => 200,
          'misc_RNA'            => 300,
          'ribozyme'            => 5,
          'rRNA'                => 190,
          'scaRNA'              => 20,
          'snoRNA'              => 200,
          'snRNA'               => 1000,
	    }, # biotypes
	}, # ncrna
	# mammals_basic
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

 }
);
  return $config{$key};

}
1;
