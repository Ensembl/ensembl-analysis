=head1 LICENSE

Copyright [2017-2024] EMBL-European Bioinformatics Institute

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

     },

     'mammals_basic' => {
       # repeats
       'dust'                        => [1392556,'repeat'],
       'repeatmask_repbase_mammals' =>  [1586621,'repeat'],
       'trf'                         => [248970,'repeat'],

       # simple features
       'cpg'                         => [9846,'simple'],
       'trnascan'                    => [328,'simple'],
       'eponine'                     => [14000,'simple'],

     },

     'rodentia_basic' => {
       # repeats
       'dust'                        => [3362887,'repeat'],
       'repeatmask_repbase_rodentia' =>  [3416232,'repeat'],
       'trf'                         => [1411057,'repeat'],

       # simple features
       'cpg'                         => [9846,'simple'],
       'trnascan'                    => [328,'simple'],
       'eponine'                     => [14000,'simple'],

     },

     'bird_basic' => {
       # repeats
       'dust'                        => [1392556,'repeat'],
       'repeatmask_repbase_mammals' =>  [1765829,'repeat'],
       'trf'                         => [248970,'repeat'],

       # simple features
       'cpg'                         => [9846,'simple'],
       'trnascan'                    => [328,'simple'],
       'eponine'                     => [14000,'simple'],

     },

     'fish_basic' => {
       'dust'                       => [500000, 'repeat'],
       'repeatmask_repbase_teleost' => [50000, 'repeat'],
       'trf'                        => [150000, 'repeat']
     },

    'sharks_basic' => {
       'dust'                       => [500000, 'repeat'],
       'repeatmask_repbase_teleost' => [50000, 'repeat'],
       'trf'                        => [150000, 'repeat']
     },

     'insects_basic' => {
       # repeats
       'dust'                       => [500000, 'repeat'],
       'repeatdetector'             => [200000, 'repeat'],
       'trf'                        => [50000, 'repeat'],

       # simple features
       'cpg'                        => [2000, 'simple'],
       'eponine'                    => [30000, 'simple'],
       'trnascan'                   => [300, 'simple'],
     },

     'reptiles_basic' => {
       # repeats
       'dust'                           => [2000000,'repeat'],
       'repeatmask_repbase_vertebrates' => [400000,'repeat'],
       'trf'                            => [400000,'repeat'],

       # simple features
       'cpg'                         => [20000,'simple'],
       'trnascan'                    => [1000,'simple'],
       'eponine'                     => [40000,'simple'],

     },

    'amphibians_basic' => {
    # repeats
     'dust'                           => [2000000, 'repeat'],
     'repeatmask_repbase_vertebrates' => [400000, 'repeat'],
     'trf'                            => [400000, 'repeat'],

    # simple features
       'cpg'                         => [20000,'simple'],
       'trnascan'                    => [1000,'simple'],
       'eponine'                     => [40000,'simple'],

    },

    'distant_vertebrate' => {
       # repeats
       'dust'                       => [2260000, 'repeat'],
       'repeatdetector'             => [153800, 'repeat'],
       'trf'                        => [504000, 'repeat'],

       # simple features
       'cpg'                        => [2000, 'simple'],
       'eponine'                    => [30000, 'simple'],
       'trnascan'                   => [300, 'simple'],
     },
      },

'gene_db_checks' => {

    'aves_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 10000,
          'genblast_not_best'   => 10000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 5000,
          'mammals_pe12_'       => 5000,
          'aves_pe12_'          => 5000,
        }, # biotypes
      }, # genblast
     'ig_tr' => {
        'logic_names' => {
          'ig_tr_gene'          => 20,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 10,
          'TR_'                 => 10,
        }, # biotypes
      }, # ig_tr
      'projection_coding' => {
        'logic_names' => {
          'project_transcripts' => 20000,
        }, # logic_names
        'biotypes' =>    {
          'projection'          => 10000,
        }, # biotypes
      }, # projection_coding
      'rnaseq_blast' =>  {
        # This one is an issue, logic names, counts are varied and one biotype is a
        # substring of the other. At the moment it's really just a check that theres'
        # some stuff in there
        'logic_names' => {
        }, # logic_names
        'biotypes' =>    {
          'rnaseq'              => 10000,
         }, # biotypes
      }, # rnaseq_blast
      'layer' => {
        'logic_names' =>    {
          'genblast'                => 100,
          'genblast_rnaseq_support' => 1000,
          'project_transcripts' => 2000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'          => 5000,
          'mammals_pe12_'         => 5000,
          'rnaseq_tissue_'       => 10000,
          'aves_pe12_'          => 5000,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
	    'ensembl'             => 15000,
	},
	'biotypes' =>    {
	}, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 350,
        }, # logic_names
        'biotypes' =>    {
          'miRNA'               => 10,
          'rRNA'                => 10,
          'snoRNA'              => 10,
          'snRNA'               => 10,
        }, # biotypes
      }, # ncrna
      'final' => {
        'logic_names' => {
          'ensembl'             => 15000,
          'ncrna'               => 400,
        }, # logic_names
        'biotypes' =>    {

        }, # biotypes
      }, # final
      'core' => {
        'logic_names' => {
          'ensembl'             => 15000,
          'ncrna'               => 400,
        }, # logic_names
        'biotypes' =>    {
          'protein_coding'       => 15000,
        }, # biotypes
      }, # core
      'otherfeatures' => {
        'logic_names' => {
          'refseq_import'          => 20000,
        },
	  'biotypes' =>    {
        }, # biotypes
      }, # otherfeatures
       'rnaseq_final' => {
	   'logic_names' => {
	   }, # logic_names
	   'biotypes' => {
	       'protein_coding'         => 5000,
	   },
      }, # rnaseq_final
    }, # aves_basic

     'reptiles_basic' => {
	      'genblast' => {
	        'logic_names' => {
	          'genblast'            => 10000,
	          'genblast_not_best'   => 10000,
	        }, # logic_names
	        'biotypes' =>    {
	          'human_pe12_'         => 50000,
	          'mammals_pe12_'       => 300000,
			  'amphibians_pe12_'    => 50000,
			  'reptiles_pe12_'      => 1200,
	          'aves_pe12_'          => 50000,
	        }, # biotypes
	      }, # genblast


	     'ig_tr' => {
	        'logic_names' => {
	          'ig_tr_gene'          => 5,
	        }, # logic_names
	        'biotypes' =>    {
	          'IG_'                 => 5,
	          'TR_'                 => 5,
	        }, # biotypes
	      }, # ig_tr


	      'projection_coding' => {
	        'logic_names' => {
			  'cesar'               => 15000,
	          'project_transcripts' => 20000,
	        }, # logic_names
	        'biotypes' =>    {
	          'projection'          => 40000,
	        }, # biotypes
	      }, # projection_coding


	      'rnaseq_blast' =>  {
	        # This one is an issue, logic names, counts are varied and one biotype is a
	        # substring of the other. At the moment it's really just a check that theres'
	        # some stuff in there
	        'logic_names' => {
	        }, # logic_names
	        'biotypes' =>    {
	          'rnaseq'              => 100000,
	         }, # biotypes
	      }, # rnaseq_blast



	      'layer' => {
	        'logic_names' =>    {
	          'genblast'                => 1600,
	          'genblast_rnaseq_support' => 5000,
	          'project_transcripts' => 600,
	        }, # logic_names
	        'biotypes' =>    {
	          'human_pe12_'         => 3000,
	          'mammals_pe12_'       => 800,
	          'rnaseq_tissue_'      => 30000,
			  'reptiles_pe12_'      => 50,
	          'aves_pe12_'          => 300,
	        }, # biotypes
	      }, # layer


	      'genebuilder' => {
	        'logic_names' => {
		    'ensembl'             => 15000,
		},
		'biotypes' =>    {
		     'rnaseq_tissue_'      => 20000,
		}, # biotypes
	      }, # genebuilder


	      'ncrna' => {
	        'logic_names' => {
	          'ncrna' => 700,
	        }, # logic_names
	        'biotypes' =>    {
	          'miRNA'               => 100,
			  'misc_RNA'			=> 10,
			  'ribozyme'            => 5,
	          'rRNA'                => 20,
			  'scaRNA'              => 15,
	          'snoRNA'              => 200,
	          'snRNA'               => 200,
	        }, # biotypes
	      }, # ncrna


	      'final' => {
	        'logic_names' => {
	          'ensembl'             => 17000,
	          'ncrna'               => 700,
	        }, # logic_names
	        'biotypes' =>    {
				 'IG_'                 => 10,
				 'TR_'                 => 10,
				 'human_pe12_'         => 1800,
				 'lncRNA'              => 4000,
				 'aves_pe12_'          => 250,
				 'mammals_pe12_'       => 1000,
				 'miRNA'               => 120,
				 'misc_RNA'            => 10,
				 'rnaseq_tissue_'      => 20000,
				 'ribozyme'            => 10,
				 'rRNA'                => 20,
				 'scaRNA'              => 15,
				 'snoRNA'              => 200,
				 'snRNA'               => 200,
	        }, # biotypes
	      }, # final


	      'core' => {
	        'logic_names' => {
	          'ensembl'             => 20000,
	          'ncrna'               => 700,
	        }, # logic_names
	        'biotypes' =>    {
			  'IG_'                  => 10,
			  'TR_'                  => 10,
			  'lncRNA'               => 4000,
			  'protein_coding'       => 21000,
			  'pseudogene'           => 50,
			  'processed_pseudogene' => 10,
			  'miRNA'                => 80,
			   'misc_RNA'             => 10,
			   'ribozyme'             => 8,
			   'rRNA'                 => 20,
			   'scaRNA'               => 10,
			   'snoRNA'               => 100,
			   'snRNA'                => 100,
	        }, # biotypes
	      }, # core


	      'otherfeatures' => {
	        'logic_names' => {
	          'refseq_import'          => 20000,
	        },
		  'biotypes' =>    {
	        }, # biotypes
	      }, # otherfeatures


	       'rnaseq_final' => {
		   'logic_names' => {
		   }, # logic_names
		   'biotypes' => {
		       'protein_coding'         => 5000,
		   },
	      }, # rnaseq_final
	    }, # reptiles_basic


     'amphibians_basic' => {
	      'genblast' => {
	        'logic_names' => {
	          'genblast'            => 10000,
	          'genblast_not_best'   => 10000,
	        }, # logic_names
	        'biotypes' =>    {
	          'human_pe12_'         => 50000,
	          'fish_pe12_'          => 50000,
			  'amphibians_pe12_'    => 50000,
			  'reptiles_pe12_'      => 1200,
	          'vert_pe12_'          => 50000,
	        }, # biotypes
	      }, # genblast


	     'ig_tr' => {
	        'logic_names' => {
	          'ig_tr_gene'          => 2,
	        }, # logic_names
	        'biotypes' =>    {
	          'IG_'                 => 2,
	          'TR_'                 => 2,
	        }, # biotypes
	      }, # ig_tr


	      'projection_coding' => {
	        'logic_names' => {
			  'cesar'               => 15000,
	          'project_transcripts' => 20000,
	        }, # logic_names
	        'biotypes' =>    {
	          'projection'          => 40000,
	        }, # biotypes
	      }, # projection_coding


	      'rnaseq_blast' =>  {
	        'biotypes' =>    {
	          'rnaseq'              => 100000,
	         }, # biotypes
             'logic_names' => {
             } #logic_names
	      }, # rnaseq_blast

	      'layer' => {
	        'logic_names' =>    {
	          'genblast'                => 1600,
	          'genblast_rnaseq_support' => 5000,
	          'project_transcripts' => 600,
	        }, # logic_names
	        'biotypes' =>    {
	          'human_pe12_'       => 3000,
	          'vert_pe12_'        => 800,
	          'rnaseq_tissue_'    => 30000,
			  'reptiles_pe12_'    => 50,
              'fish_pe12_'        => 50,
	          'amphibians_pe12_'  => 300,
	        }, # biotypes
	      }, # layer


	      'genebuilder' => {
	        'logic_names' => {
		    'ensembl'             => 15000,
		},
		'biotypes' =>    {
		     'rnaseq_tissue_'      => 20000,
		}, # biotypes
	      }, # genebuilder


	      'ncrna' => {
	        'logic_names' => {
	          'ncrna' => 700,
	        }, # logic_names
	        'biotypes' =>    {
	          'miRNA'               => 100,
			  'misc_RNA'			=> 10,
			  'ribozyme'            => 5,
	          'rRNA'                => 20,
			  'scaRNA'              => 15,
	          'snoRNA'              => 200,
	          'snRNA'               => 200,
	        }, # biotypes
	      }, # ncrna


	      'final' => {
	        'logic_names' => {
	          'ensembl'             => 17000,
	          'ncrna'               => 700,
	        }, # logic_names
	        'biotypes' =>    {
				 'IG_'                 => 2,
				 'TR_'                 => 2,
				 'human_pe12_'         => 1800,
				 'lncRNA'              => 1000,
				 'amphibians_pe12_'    => 250,
				 'reptiles_pe12_'      => 100,
                 'fish_pe12_'          => 100,
                 'vert_pe12_'          => 1000,
				 'miRNA'               => 120,
				 'misc_RNA'            => 10,
				 'rnaseq_tissue_'      => 20000,
				 'ribozyme'            => 10,
				 'rRNA'                => 20,
				 'scaRNA'              => 15,
				 'snoRNA'              => 200,
				 'snRNA'               => 200,
	        }, # biotypes
	      }, # final


	      'core' => {
	        'logic_names' => {
	          'ensembl'             => 20000,
	          'ncrna'               => 700,
	        }, # logic_names
	        'biotypes' =>    {
			  'IG_'                  => 2,
			  'TR_'                  => 2,
			  'lncRNA'               => 4000,
			  'protein_coding'       => 18000,
			  'pseudogene'           => 50,
			  'processed_pseudogene' => 10,
			  'miRNA'                => 80,
			   'misc_RNA'             => 10,
			   'ribozyme'             => 8,
			   'rRNA'                 => 20,
			   'scaRNA'               => 10,
			   'snoRNA'               => 100,
			   'snRNA'                => 100,
	        }, # biotypes
	      }, # core


	      'otherfeatures' => {
	        'logic_names' => {
	          'refseq_import'          => 20000,
	        },#logic_names
	      }, # otherfeatures


	       'rnaseq_final' => {
		   'biotypes' => {
		       'protein_coding'         => 5000,
		   },#biotypes
	      }, # rnaseq_final
	    }, # amphibians_basic


    'mammals_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 10000,
          'genblast_not_best'   => 10000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 1000,
	  'mouse_pe12_'      => 1000,
	  'mammals_pe12_'      => 500,
        }, # biotypes
      }, # genblast
      'ig_tr' => {
        'logic_names' => {
          'ig_tr_gene'          => 20,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 10,
          'TR_'                 => 10,
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
      'rnaseq_blast' =>  {
        # This one is an issue, logic names, counts are varied and one biotype is a
        # substring of the other. At the moment it's really just a check that theres'
        # some stuff in there
        'logic_names' => {
        }, # logic_names
        'biotypes' =>    {
          'rnaseq'              => 10000,
         }, # biotypes
      }, # rnaseq_blast
      'layer' => {
        'logic_names' =>    {
          'genblast'                => 100,
          'best_targetted'          => 100,
          'project_transcripts' => 10000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'          => 1000,
          'mouse_pe12_'       => 1000,
          'rnaseq_tissue_'       => 10000,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
          'ensembl'             => 19000,
        }, # logic_names
	 'biotypes' =>    {
        }, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 1000,
        }, # logic_names
        'biotypes' =>    {#these numbers are very low, it's simply a check that something is there
          'miRNA'               => 10,
          'misc_RNA'            => 10,
          'rRNA'                => 10,
          'snoRNA'              => 10,
          'snRNA'               => 10,
        }, # biotypes
      }, # ncrna
      'final' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 600,
        }, # logic_names
        'biotypes' =>    {
        }, # biotypes
      }, # final
      'core' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 600,
        }, # logic_names
        'biotypes' =>    {
          'protein_coding'       => 19000,
        }, # biotypes
      }, # core
      'otherfeatures' => {
        'logic_names' => {
          'refseq_import'          => 20000,
        },
	 'biotypes' =>    {
        }, # biotypes
      }, # otherfeatures
      'rnaseq_final' => {
        'logic_names' => {
        }, # logic_names
	'biotypes' => {
          'protein_coding'         => 5000,
        },
      }, # rnaseq_final
    }, # mammals_basic

    'primates_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 40000,
          'genblast_not_best'   => 60000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 20000,
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
      'rnaseq_blast' =>  {
        # This one is an issue, logic names, counts are varied and one biotype is a
        # substring of the other. At the moment it's really just a check that theres'
        # some stuff in there
        'logic_names' => {

        }, # logic_names
        'biotypes' =>    {
          'rnaseq'              => 10000,
         }, # biotypes
      }, # rnaseq_blast
      'layer' => {
	  # This is pretty variable so we need to think of a sensible plan for the logic_name checks
	  'logic_names' =>    {
          'genblast'                => 100,
          'best_targetted'          => 100,
          'genblast_rnaseq_support' => 1000,
          'project_transcripts' => 10000,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                  => 20,
          'TR_'                  => 20,
	  'human_pe12_'          => 2000,
	  #projection from human is quite relaible in primates so having this set high is reasonable
	  'projection_'          => 20000,
	  #I'm leaving this number low so it's just a check to see if some rnaseq stuff exists
          'rnaseq_tissue_'       => 100,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
          'ensembl'             => 19000,
        }, # logic_names
	 'biotypes' =>    {
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
          'rRNA'                => 50,
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
          'miRNA'               => 500,
          'misc_RNA'            => 1000,
          'ribozyme'            => 0,
          'rRNA'                => 50,
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
          'protein_coding'       => 19000,
          'pseudogene'           => 50,
          'processed_pseudogene' => 0,
          'miRNA'                => 500,
          'misc_RNA'             => 1000,
          'ribozyme'             => 0,
          'rRNA'                 => 50,
          'scaRNA'               => 0,
          'snoRNA'               => 200,
          'snRNA'                => 400,
        }, # biotypes
      }, # core
   'rnaseq_final' => {
       'logic_names' => {
       }, # logic_names
       'biotypes' => {
	   'protein_coding'         => 5000,
       },
      }, # rnaseq_final
    }, # primates_basic

    'fish_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 20000,
          'genblast_not_best'   => 40000,
        }, # logic_names
        'biotypes' =>    {
          'fish_pe12_' => 20000,
          'human_pe12_'=> 10000,
        }, # biotypes
      }, # genblast
      'ig_tr' => {
        'logic_names' => {
          'ig_tr_gene'          => 10,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 10,
          'TR_'                 => 1,
        }, # biotypes
      }, # ig_tr
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
        'logic_names' =>    {
          'genblast'                => 10000,
          'genblast_rnaseq_support' => 20000,
      }, # logic_names
        'biotypes' =>    {
          'IG_'                  => 10,
          'TR_'                  => 1,
          'fish_pe12_'          => 10000,
          'human_pe12_'          => 1000,
          'rnaseq_tissue_'       => 20000,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
          'ensembl'             => 19000,
        }, # logic_names
	    'biotypes' =>    {
		'rnaseq_tissue_'      => 10000,
        }, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 300,
        }, # logic_names
        'biotypes' =>    {
          'miRNA'               => 20,
          'misc_RNA'            => 5,
          'rRNA'                => 10,
          'snoRNA'              => 80,
          'snRNA'               => 80,
        }, # biotypes
      }, # ncrna
      'final' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 500,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 10,
          'TR_'                 => 2,
          'fish_pe12_'          => 5000,
          'human_pe12_'         => 1000,
          'rnaseq_tissue_'      => 20000,
          'lncRNA'              => 2000,
          'miRNA'               => 30,
          'misc_RNA'            => 5,
          'rRNA'                => 50,
          'snoRNA'              => 50,
          'snRNA'               => 100,
        }, # biotypes
      }, # final
      'core' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 400,
        }, # logic_names
        'biotypes' =>    {
          'protein_coding'       => 19000,
          'IG_'                 => 10,
          'TR_'                 => 2,
          'pseudogene'           => 50,
          'miRNA'                => 25,
          'misc_RNA'             => 5,
          'rRNA'                 => 50,
          'snoRNA'               => 50,
          'snRNA'                => 75,
        }, # biotypes
      }, # core
      'otherfeatures' => {
        'logic_names' => {
          'refseq_import'          => 18000,
        },
        'biotypes' =>    {
        }, # biotypes
      }, # otherfeatures
      'rnaseq_final' => {
        'logic_names' => {
       }, # logic_names
        'biotypes' => {
          'protein_coding'      => 50000,
        }, # biotypes
      }, # rnaseq_final
    }, # fish_basic

    'sharks_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 20000,
          'genblast_not_best'   => 40000,
        }, # logic_names
        'biotypes' =>    {
          'fish_pe12_' => 20000,
          'human_pe12_'=> 10000,
        }, # biotypes
      }, # genblast
      'ig_tr' => {
        'logic_names' => {
          'ig_tr_gene'          => 0,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 0,
          'TR_'                 => 0,
        }, # biotypes
      }, # ig_tr
      'rnaseq_blast' =>  {
        # This one is an issue, logic names, counts are varied and one biotype is a
        # substring of the other. At the moment it's really just a check that theres'
        # some stuff in there
        'logic_names' => {
        }, # logic_names
        'biotypes' =>    {
          'rnaseq'              => 10000,
         }, # biotypes
      }, # rnaseq_blast
      'layer' => {
        'logic_names' =>    {
          'genblast'                => 10000,
          'genblast_rnaseq_support' => 20000,
      }, # logic_names
        'biotypes' =>    {
          'IG_'                  => 10,
          'TR_'                  => 1,
          'fish_pe12_'          => 10000,
          'human_pe12_'          => 1000,
          'rnaseq_tissue_'       => 20000,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
          'ensembl'             => 18000,
        }, # logic_names
	'biotypes' =>    {
	    'rnaseq_tissue_'      => 10000,
        }, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 300,
        }, # logic_names
        'biotypes' =>    {
          'miRNA'               => 1,
          'misc_RNA'            => 5,
          'rRNA'                => 10,
          'snoRNA'              => 80,
          'snRNA'               => 80,
        }, # biotypes
      }, # ncrna
      'final' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 500,
        }, # logic_names
        'biotypes' =>    {
          'fish_pe12_'          => 5000,
          'human_pe12_'         => 1000,
          'rnaseq_tissue_'      => 15000,
          'lncRNA'              => 1000,
          'miRNA'               => 1,
          'misc_RNA'            => 5,
          'rRNA'                => 50,
          'snoRNA'              => 50,
          'snRNA'               => 100,
        }, # biotypes
      }, # final
      'core' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 400,
        }, # logic_names
        'biotypes' =>    {
          'protein_coding'       => 18000,
          'pseudogene'           => 20,
          'miRNA'                => 1,
          'misc_RNA'             => 5,
          'rRNA'                 => 50,
          'snoRNA'               => 50,
          'snRNA'                => 75,
        }, # biotypes
      }, # core
      'otherfeatures' => {
        'logic_names' => {
          'refseq_import'          => 18000,
        },
        'biotypes' =>    {
        }, # biotypes
      }, # otherfeatures
      'rnaseq_final' => {
        'logic_names' => {
       }, # logic_names
        'biotypes' => {
          'protein_coding'      => 50000,
        }, # biotypes
      }, # rnaseq_final
    }, # shark_basic


    'rodentia_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 125000,
          'genblast_not_best'   => 100000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 20000,
          'mouse_pe12_'      => 50000,
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
        'logic_names' =>    {
          'genblast'                => 100,
          'genblast_rnaseq_support' => 1000,
          'project_transcripts' => 30000,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                  => 20,
          'TR_'                  => 20,
          'human_pe12_'          => 1000,
          'mouse_pe12_'          => 1000,
          'mammals_pe12_'        => 1000,
          'rnaseq_tissue_'       => 10000,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
          'ensembl'             => 19000,
        }, # logic_names
	'biotypes' =>    {
	}, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 1000,
        }, # logic_names
        'biotypes' =>    {
          'miRNA'               => 100,
          'misc_RNA'            => 100,
          'ribozyme'            => 10,
          'rRNA'                => 200,
          'scaRNA'              => 500,
          'snoRNA'              => 500,
          'snRNA'               => 500,
        }, # biotypes
      }, # ncrna
      'final' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 600,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 20,
          'TR_'                 => 20,
          'human_pe12_'         => 1000,
          'lncRNA'              => 2000,
          'mouse_pe12_'         => 1000,
          'mammals_pe12_'       => 1000,
          'miRNA'               => 100,
          'misc_RNA'            => 10,
          'rnaseq_tissue_'      => 20000,
          'ribozyme'            => 10,
          'rRNA'                => 20,
          'scaRNA'              => 10,
          'snoRNA'              => 100,
          'snRNA'               => 100,
        }, # biotypes
      }, # final
      'core' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 600,
        }, # logic_names
        'biotypes' =>    {
          'IG_'                  => 20,
          'TR_'                  => 20,
          'protein_coding'       => 21000,
          'pseudogene'           => 50,
          'processed_pseudogene' => 10,
          'miRNA'                => 80,
          'misc_RNA'             => 10,
          'ribozyme'             => 10,
          'rRNA'                 => 20,
          'scaRNA'               => 10,
          'snoRNA'               => 100,
          'snRNA'                => 100,
        }, # biotypes
      }, # core
      'otherfeatures' => {
        'logic_names' => {
          'refseq_import'          => 20000,
        },
        'biotypes' =>    {
        }, # biotypes
      }, # otherfeatures
      'rnaseq_final' => {
        'logic_names' => {
       }, # logic_names
        'biotypes' => {
          'protein_coding'         => 5000,
        },
      }, # rnaseq_final
    }, # rodentia_basic
    'distant_vertebrate' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 125000,
          'genblast_not_best'   => 125000,
        }, # logic_names
        'biotypes' =>    {
          'vert_pe12_'         => 100000,
        }, # biotypes
      }, # genblast
      'ig_tr' => {
        'logic_names' => {
          'ig_tr_gene'          => 0,#setting to 0 until we have better ig detection outside of human
        }, # logic_names
        'biotypes' =>    {
          'IG_'                 => 0,
          'TR_'                 => 0,
        }, # biotypes
      }, # ig_tr
      'projection_coding' => {
        'logic_names' => {
          'project_transcripts' => 10000,
        }, # logic_names
        'biotypes' =>    {
          'projection_'          => 10000,
        }, # biotypes
      }, # projection_coding
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
	  'logic_names' =>    {
          'genblast'                => 1000,
          'best_targetted'          => 1000,
          'genblast_rnaseq_support' => 1000,
          'project_transcripts' => 10000,
	  },
	  'biotypes' =>    {
	      'verts_pe12_'        => 20000,
	      'rnaseq_tissue_'       => 100,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
          'ensembl'             => 20000,
        }, # logic_names
        'biotypes' =>    {
        }, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 500,
        }, # logic_names
        'biotypes' =>    {
          'miRNA'               => 5,
          'misc_RNA'            => 1,
          'ribozyme'            => 0,
          'rRNA'                => 20,
          'scaRNA'              => 0,
          'snoRNA'              => 20,
          'snRNA'               => 20,
        }, # biotypes
      }, # ncrna
      'final' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 500,
        }, # logic_names
        'biotypes' =>    {
          'protein_coding'       => 18000,
        }, # biotypes
      }, # final
      'core' => {
        'logic_names' => {
          'ensembl'             => 18000,
          'ncrna'               => 500,
        }, # logic_names
        'biotypes' =>    {
          'protein_coding'       => 18000,
        }, # biotypes
      }, # core
    }, # distant_vertebrate_basic

   }, # gene_db_checks


);
  return $config{$key};
}

1;
