#
# package Bio::EnsEMBL::Analysis::Config::ImportArrays
# 
# Cared for by EnsEMBL (dev@ensembl.org)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Config::ImportArrays

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis:Config::ImportAarrys;

=head1 DESCRIPTION

This contains the configuration for importing arrays from flat files.
It is entirely dependant on the arrays.env environment which can be used 
to set up and run the pipeline in an easy and interactive way. This contains 
all possible configurations which will then be set dynamically by the RunnableDB
for each instance using the input_id as a key into a separate ImportArrays.conf 
file, listed here as ARRAY_FORMAT_FILE.


The layout of the configuration is a set of hashes,
each one keyed by logic name. There is also a DEFAULT hash,
which is used as the default for all logic names (this
was the configuration pattern stolen from Exonerate2Genes,
although in this case it's very unlikely you will need to have
different configs by logic name).

=head1 CONTACT

=cut


# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/Config/ImportArrays.pm,v $
# $Version: $
package Bio::EnsEMBL::Analysis::Config::ImportArrays;

use warnings ;
use strict;
use vars qw( %Config );

# Hash containing config info
# -- one hashnode per logic name, with a 'DEFAULT' logic name provided
#

%Config = 
  (

   #This entire hash is exported as the global $ARRAY_CONFIG var
   #each key will be exported as $ARRAY_CONFIG->{'_CONFIG_'.$key}
   #Dependant on logic name of RunnableDB

   ARRAY_CONFIG => 
   {
    DEFAULT => 
    {
     #These are now defined dynamically or via the ImportArrays.conf file
     # All input probes must be kept in one huge (possibly redundant) fasta file
     #QUERYSEQS            => $ENV{'RAW_FASTA'},
     # The output of this module writes a set of affy probes into the OUTDB.affy_probe table,
     # and also writes the nonredundant probes into this fasta file,
     # with the fasta headers keyed with the affy probes' internal id. 
     #NON_REDUNDANT_PROBE_SEQS => $ENV{'NR_FASTA'},
	 
     # DB containing all affy_arrays, affy_probes and (next step) affy_features
     OUTDB => {
               -dbname          => $ENV{'DB_NAME'},
               -host            => $ENV{'DB_HOST'},
               -port            => $ENV{'DB_PORT'},
               -user            => $ENV{'DB_USER'},
               -pass            => $ENV{'DB_PASS'},
               -species         => $ENV{'SPECIES'}, #Only here until we fix the DBAadptor new method
               -multispecies_db => $ENV{'MULTISPECIES_DB'},
               -species_id      => $ENV{'SPECIES_ID'}
              },


     #Optional, must define if dnadb is not on ensembldb
     #Not used, but will fail if dnadb autoguessing fails
     DNADB => {
               -dbname          => $ENV{'DNADB_NAME'},
               -host            => $ENV{'DNADB_HOST'},
               -port            => $ENV{'DNADB_PORT'},
               -user            => $ENV{'DNADB_USER'},
               -pass            => $ENV{'DNADB_PASS'},
               -species         => $ENV{'SPECIES'},
               -multispecies_db => $ENV{'DNADB_MULTISPECIES_DB'},
               -species_id      => $ENV{'DNADB_SPECIES_ID'}
              },
	 
     #Used for building the format specific NR fasta file
     OUTPUT_DIR           => $ENV{'WORK_DIR'},


     #This defines how to parse the file headers
     IIDREGEXP =>  '^>probe:(\S+):(\S+):(\S+:\S+;).*$', #AFFY
				  
     #We also need a has to define the input field order
     #This will be used to set the relevant hash values
     IFIELDORDER => {
                     #do we need to add fields for class to enable skipping on control probes
                     #here and in regexp
                     #We duplicate the field 0 between array.name and array_chip.design_id
                     #-name       => 2,
                     #-array      => 0,
                     #-array_chip => 0,
                     #-probe_set   => 1,
                    },

     #ISKIPLIST/REGEX
     #ISKIPFIELD



	 

     ARRAY_PARAMS => {
                      #'MG-U74Cv2' => {
                      #				  -name => 'MG-U74Cv2',
                      #				  -vendor => 'AFFY',
                      #				  #-setsize => undef,
                      #				  -format  => 'EXPRESSION',
                      #				  -type    => 'OLIGO', 
                      #  -class => 'AFFY_ST',
                      #				  #-description => '',
                      #				 },

                      # 'MoGene-1_0-st-v1' => {
                      #						 -name => 'MoGene-1_0-st-v1',
                      #						 -vendor => 'AFFY',
                      #						 #-setsize => undef,
                      #						 -format  => 'EXPRESSION',
                      #						 -type    => 'OLIGO',
                      #						 #-description => '',
                      #  -class => 'AFFY_ST',
                      #						},


                     },

	 
    },


    #%{$Config::ArrayMapping::import_arrays},

    IMPORT_AFFY_UTR_ARRAYS => 
    {
     IIDREGEXP => '^>probe:(\S+):(\S+):(\S+:\S+;).*$',
	
     IFIELDORDER => {
                     -name       => 2, -array_chip => 0,
                     -array      => 0, -probe_set   => 1
                    },
 	 
     #Can we remove name from these hashes?

     ARRAY_PARAMS => 
     {

      #Remove all the redundant values and set them in ImportArrays?

      Porcine =>  {
                   -name => 'Porcine',
                   -vendor => 'AFFY',
                   #-setsize => undef,
                   -format  => 'EXPRESSION',
                   -type    => 'OLIGO',
                   #-description => '',
                   -class   => 'AFFY_UTR',
                  },



      #Platypus
      #NOTE: This is not actually AFFY_UTR, but a custom array 
      #which we are treating the same. Remember to manually change
      #array.vendor/class to CUSTOM after RunTranscriptXrefs!

      #IIDREGEXP => '^>probe:(\S+):(\S+):(\S+);.*$',
      'platypus_exon' => {
                          -name => 'platypus_exon',
                          -vendor => 'CUSTOM',
                          #-setsize => undef,
                          -format  => 'EXPRESSION',
                          -type    => 'OLIGO',
                          #-description => '',
                          -class   => 'CUSTOM',
                         },
					  



      #Frog
      'X_tropicalis' => {
                         -name => 'X_tropicalis',
                         -vendor => 'AFFY',
                         -format  => 'EXPRESSION',
                         -type    => 'OLIGO',
                         #-description => '',
                         -class   => 'AFFY_UTR',
                        },



      #Dog
      'Canine_2' => {
                     -name => 'Canine_2',
                     -vendor => 'AFFY',
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     #-description => '',
                     -class   => 'AFFY_UTR',
                    },


      #Macaque
      'Rhesus' => {
                   -name => 'Rhesus',
                   -vendor => 'AFFY',
                   -format  => 'EXPRESSION',
                   -type    => 'OLIGO',
                   #-description => '',
                   -class   => 'AFFY_UTR',
                  },

      #C intestinalis
      'CINT06a520380F' => {
                           -name => 'CINT06a520380F',
                           -vendor => 'AFFY',
                           -format  => 'EXPRESSION',
                           -type    => 'OLIGO',
                           #-description => '',
                           -class   => 'AFFY_UTR',
                          },

      #Cow
      'Bovine' => {
                   -name => 'Bovine',
                   -vendor => 'AFFY',
                   -format  => 'EXPRESSION',
                   -type    => 'OLIGO',
                   #-description => '',
                   -class   => 'AFFY_UTR',
                  },



      #Chicken
      'Chicken' => {
                    -name => 'Chicken',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },


      #C_elegans

      'C_elegans' => {
                      -name => 'C_elegans',
                      -vendor => 'AFFY',
                      -format  => 'EXPRESSION',
                      -type    => 'OLIGO',
                      #-description => '',
                      -class   => 'AFFY_UTR',
                     },


      #Zebra fish
      'Zebrafish' => {
                      -name => 'Zebrafish',
                      -vendor => 'AFFY',
                      -format  => 'EXPRESSION',
                      -type    => 'OLIGO',
                      #-description => '',
                      -class   => 'AFFY_UTR',
                     },


      #Rat

      'RAE230A' => {
                    -name => 'RAE230A',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },


      'RAE230B' => {
                    -name => 'RAE230B',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },
      'Rat230_2' => {
                     -name => 'Rat230_2',
                     -vendor => 'AFFY',
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     #-description => '',
                     -class   => 'AFFY_UTR',
                    },



      'RG-U34A' => {
                    -name => 'RG-U34A',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },
      'RG-U34B' => {
                    -name => 'RG-U34B',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },

      'RG-U34C' => {
                    -name => 'RG-U34C',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },

      'RN-U34' => {
                   -name => 'RN-U34',
                   -vendor => 'AFFY',
                   -format  => 'EXPRESSION',
                   -type    => 'OLIGO',
                   #-description => '',
                   -class   => 'AFFY_UTR',
                  },


      'RT-U34' => {
                   -name => 'RT-U34',
                   -vendor => 'AFFY',
                   -format  => 'EXPRESSION',
                   -type    => 'OLIGO',
                   #-description => '',
                   -class   => 'AFFY_UTR',
                  },


      #Human
	
      'HC-G110' => {
                    -name => 'HC-G110',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },


      'U133_X3P' => {
                     -name => 'U133_X3P',
                     -vendor => 'AFFY',
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     #-description => '',
                     -class   => 'AFFY_UTR',
                    },

	  
      'HuGeneFL' => {-name => 'HuGeneFL',
                     -vendor => 'AFFY',
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     #-description => '',
                     -class   => 'AFFY_UTR'},

      'HG_U95A' => {-name => 'HG-U95A',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR'},

      'HG-U95E' => {-name => 'HG-U95E',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR'},

      'HG-U95D' => {-name => 'HG-U95D',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR'},

      'HG-U95C' => {-name => 'HG-U95C',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR'},

      'HG-U95B' => {-name => 'HG-U95B',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR'},
	  
      'HG_U95Av2' => {-name => 'HG-U95Av2',
                      -vendor => 'AFFY',
                      -format  => 'EXPRESSION',
                      -type    => 'OLIGO',
                      #-description => '',
                      -class   => 'AFFY_UTR'},

      'HG-U133_Plus_2' => {-name => 'HG-U133_Plus_2',
                           -vendor => 'AFFY',
                           -format  => 'EXPRESSION',
                           -type    => 'OLIGO',
                           #-description => '',
                           -class   => 'AFFY_UTR'},

      'HG-U133B' => {-name => 'HG-U133B',
                     -vendor => 'AFFY',
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     #-description => '',
                     -class   => 'AFFY_UTR'},

      'HG-U133A' => {-name => 'HG-U133A',
                     -vendor => 'AFFY',
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     #-description => '',
                     -class   => 'AFFY_UTR'},

      'HG-U133A_2' => {-name => 'HG-U133A_2',
                       -vendor => 'AFFY',
                       -format  => 'EXPRESSION',
                       -type    => 'OLIGO',
                       #-description => '',
                       -class   => 'AFFY_UTR'},

      'HG-Focus' => {-name => 'HG-Focus',
                     -vendor => 'AFFY',
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     #-description => '',
                     -class   => 'AFFY_UTR'},





      #Mouse
      'MG-U74Cv2' => {
                      -name => 'MG-U74Cv2',
                      -vendor => 'AFFY',
                      -format  => 'EXPRESSION',
                      -type    => 'OLIGO',
                      #-description => '',
                      -class   => 'AFFY_UTR',
                     },


      'MG-U74A' => {
                    -name => 'MG-U74A',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
					  
                   },
	  
      'MG-U74Av2' => {
                      -name => 'MG-U74Av2',
                      -vendor => 'AFFY',
                      -format  => 'EXPRESSION',
                      -type    => 'OLIGO',
                      #-description => '',
                      -class   => 'AFFY_UTR',
                     },
	  


      'MG-U74B' => {
                    -name => 'MG-U74B',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },
	  
      'MG-U74Bv2' => {
                      -name => 'MG-U74Bv2',
                      -vendor => 'AFFY',
                      -format  => 'EXPRESSION',
                      -type    => 'OLIGO',
                      #-description => '',
                      -class   => 'AFFY_UTR',
                     },
	  

      'MG-U74C' => {
                    -name => 'MG-U74C',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },

      'MOE430A' => {
                    -name => 'MOE430A',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },

      'MOE430B' => {
                    -name => 'MOE430B',
                    -vendor => 'AFFY',
                    -format  => 'EXPRESSION',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },

      'Mouse430A_2' => {
                        -name => 'Mouse430A_2',
                        -vendor => 'AFFY',
                        -format  => 'EXPRESSION',
                        -type    => 'OLIGO',
                        #-description => '',
                        -class   => 'AFFY_UTR',
                       },

	  
      'Mouse430_2' => {
                       -name => 'Mouse430_2',
                       -vendor => 'AFFY',
                       -format  => 'EXPRESSION',
                       -type    => 'OLIGO',
                       #-description => '',
                       -class   => 'AFFY_UTR',
                      },



	    'Mu11KsubA' => {
                      -name => 'Mu11KsubA',
                      -vendor => 'AFFY',
                      -format  => 'EXPRESSION',
                      -type    => 'OLIGO',
                      #-description => '',
                      -class   => 'AFFY_UTR',
                     },
      'Mu11KsubB' => {
                      -name => 'Mu11KsubB',
                      -vendor => 'AFFY',
                      -format  => 'EXPRESSION',
                      -type    => 'OLIGO',
                      #-description => '',
                      -class   => 'AFFY_UTR',
                     },
	  
      #Drosophila
      'DrosGenome1' => {-name => 'DrosGenome1',
                        -vendor => 'AFFY',
                        -format  => 'EXPRESSION',
                        -type    => 'OLIGO',
                        #-description => '',
                        -class   => 'AFFY_UTR'},
	  
      'Drosophila_2' => {-name => 'Drosophila_2',
                         -vendor => 'AFFY',
                         -format  => 'EXPRESSION',
                         -type    => 'OLIGO',
                         #-description => '',
                         -class   => 'AFFY_UTR'},
	  	

      #Yeast
      'Yeast_2' => {
                    -name => 'Yeast_2',
                    -vendor => 'AFFY',
                    #-setsize => undef,
                    -format  => 'EXPRESSION', #? UTR?
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AFFY_UTR',
                   },

      'YG-S98' => {
                   -name => 'YG-S98',
                   -vendor => 'AFFY',
                   #-setsize => undef,
                   -format  => 'EXPRESSION', #? UTR?
                   -type    => 'OLIGO',
                   #-description => '',
                   -class   => 'AFFY_UTR',
                  },


					
      #EColi
      'E_coli_2' => {
                     -name => 'E_coli_2',
                     -vendor=>'AFFY',
                     -format => 'EXPRESSION',
                     -type=>'OLIGO',
                     -class=>'AFFY_UTR'
                    },	  

      'E_coli_Antisense' => {
                             -name => 'E_coli_Antisense',
                             -vendor => 'AFFY',
                             -format => 'EXPRESSION',
                             -type => 'OLIGO',
                             -class => 'AFFY_UTR'
                            },
      #S_aureus

      'S_aureus' => {
                     -name => 'S_aureus',
                     -vendor=>'AFFY',
                     -format => 'EXPRESSION',
                     -type=>'OLIGO',
                     -class=>'AFFY_UTR'
                    },
					    
      # plants

      'ATH1-121501' => {
                        -name => 'ATH1-121501',
                        -vendor=>'AFFY',
                        -format => 'EXPRESSION',
                        -type=>'OLIGO',
                        -class=>'AFFY_UTR'
                       },
 
      'Rice' => {
                 -name => 'Rice',
                 -vendor=>'AFFY',
                 -format => 'EXPRESSION',
                 -type=>'OLIGO',
                 -class=>'AFFY_UTR'
                },

      'Poplar' => {
                   -name => 'Poplar',
                   -vendor=>'AFFY',
                   -format => 'EXPRESSION',
                   -type=>'OLIGO',
                   -class=>'AFFY_UTR'
                  },

      'Vitis_Vinifera' => {
                           -name => 'Vitis_Vinifera',
                           -vendor=>'AFFY',
                           -format => 'EXPRESSION',
                           -type=>'OLIGO',
                           -class=>'AFFY_UTR'
                          },				     


      #Then add user defined/custom ones here?
      #values %{$ArrayConfig->{ARRAY_PARAMS}}
      #Could write this automatically from env or script?

     },
	 
     INPUT_FORMAT => 'FASTA',
    },

    IMPORT_AFFY_ST_ARRAYS => 
    {
     IIDREGEXP => '^>probe:(\S+):(\S+);\S+:\S+;.*[TranscriptCluster|ProbeSet]ID=([0-9]+);.*$',
	 
     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     -probe_set   => 2,
                    },
	 	 
     ARRAY_PARAMS => {

				

                      #Rat
                      'RaEx-1_0-st-v1' => {
                                           -name => 'RaEx-1_0-st-v1',
                                           -vendor => 'AFFY',
                                           #-setsize => undef,
                                           -format  => 'EXPRESSION',
                                           -type    => 'OLIGO',
                                           #-description => '',
                                           -class   => 'AFFY_ST',
                                          },

					  
                      'RaGene-1_0-st-v1' => {
                                             -name => 'RaGene-1_0-st-v1',
                                             -vendor => 'AFFY',
                                             #-setsize => undef,
                                             -format  => 'EXPRESSION',
                                             -type    => 'OLIGO',
                                             #-description => '',
                                             -class   => 'AFFY_ST',
                                            },


                      #Human
                      'HuGene-1_0-st-v1' => {
                                             -name => 'HuGene-1_0-st-v1',
                                             -vendor => 'AFFY',
                                             #-setsize => undef,
                                             -format  => 'EXPRESSION',
                                             -type    => 'OLIGO',
                                             #-description => '',
                                             -class   => 'AFFY_ST',
                                            },


                      'HuEx-1_0-st-v2' => {
                                           -name => 'HuEx-1_0-st-v2',
                                           -vendor => 'AFFY',
                                           -format  => 'EXPRESSION',
                                           -type    => 'OLIGO',
                                           #-description => '',
                                           -class   => 'AFFY_ST',
                                          }, 


                      #Mouse
                      'MoGene-1_0-st-v1' => {
                                             -name => 'MoGene-1_0-st-v1',
                                             -vendor => 'AFFY',
                                             #-setsize => undef,
                                             -format  => 'EXPRESSION',
                                             -type    => 'OLIGO',
                                             #-description => '',
                                             -class   => 'AFFY_ST',
                                            },


                      'MoEx-1_0-st-v1' => {
                                           -name => 'MoEx-1_0-st-v1',
                                           -vendor => 'AFFY',
                                           -format  => 'EXPRESSION',
                                           -type    => 'OLIGO',
                                           #-description => '',
                                           -class   => 'AFFY_ST',
                                          },
					  


                     },
	 
     INPUT_FORMAT => 'FASTA',
    },

    IMPORT_ILLUMINA_WG_ARRAYS => 
    {
     IIDREGEXP => '^>(\S+):(\S+).*$',
	 
     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     #-probe_set   => 2,#This could be annotation
                    },
	 	 
     ARRAY_PARAMS => {
					  
                      'MouseWG_6_V1' => {
                                         -name => 'MouseWG_6_V1',
                                         -vendor => 'ILLUMINA',
                                         #-setsize => undef,
                                         -format  => 'EXPRESSION',
                                         -type    => 'OLIGO',
                                         #-description => '',
                                         -class   => 'ILLUMINA_WG',
                                        },
					  
					  
                      'MouseWG_6_V2' => {
                                         -name => 'MouseWG_6_V2',
                                         -vendor => 'ILLUMINA',
                                         #-setsize => undef,
                                         -format  => 'EXPRESSION',
                                         -type    => 'OLIGO',
                                         #-description => '',
                                         -class   => 'ILLUMINA_WG',
                                        },
                      #V1 is no longer accesible via website?
                      #Only on ftp site
                      'HumanWG_6_V1' => {
                                         -name => 'HumanWG_6_V1',
                                         -vendor => 'ILLUMINA',
                                         #-setsize => undef,
                                         -format  => 'EXPRESSION',
                                         -type    => 'OLIGO',
                                         #-description => '',
                                         -class   => 'ILLUMINA_WG',
                                        },

                      'HumanHT-12' => {
                                       -name => 'HumanHT-12',
                                       -vendor => 'ILLUMINA',
                                       #-setsize => undef,
                                       -format  => 'EXPRESSION',
                                       -type    => 'OLIGO',
                                       #-description => '',
                                       -class   => 'ILLUMINA_WG',
                                      },




                      'HumanWG_6_V2' => {
                                         -name => 'HumanWG_6_V2',
                                         -vendor => 'ILLUMINA',
                                         #-setsize => undef,
                                         -format  => 'EXPRESSION',
                                         -type    => 'OLIGO',
                                         #-description => '',
                                         -class   => 'ILLUMINA_WG',
                                        },

                      'HumanWG_6_V3' => {
                                         -name => 'HumanWG_6_V3',
                                         -vendor => 'ILLUMINA',
                                         #-setsize => undef,
                                         -format  => 'EXPRESSION',
                                         -type    => 'OLIGO',
                                         #-description => '',
                                         -class   => 'ILLUMINA_WG',
                                        },

                      'RatRef-12' => {
                                      -name => 'RatRef-12',
                                      -vendor => 'ILLUMINA',
                                      #-setsize => undef,
                                      -format  => 'EXPRESSION',
                                      -type    => 'OLIGO',
                                      #-description => '',
                                      -class   => 'ILLUMINA_WG',
                                     },

					  
					  

                     },
	 
     INPUT_FORMAT => 'FASTA',
    },



    IMPORT_ILLUMINA_INFINIUM_ARRAYS => 
    {

     #These are DNA methylation arrays so need a reverse strand mapping?
     #Or are SourceSeqs target seqs? So we don't need to do anything?

     IIDREGEXP => '^>(\S+):(\S+).*$',
	 
     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     #-probe_set   => 2,#This could be annotation
                    },
	 	 
     ARRAY_PARAMS => {
					  
					  
                      'HumanMethylation27' => {
                                               -name => 'HumanMethylation27',
                                               -vendor => 'ILLUMINA',
                                               -format  => 'METHYLATION',
                                               -type    => 'OLIGO',
                                               #-description => '',
                                               -class   => 'ILLUMINA_INFINIUM',
                                               #-version => '1.2',
                                              },

                      'HumanMethylation450' => {
                                                -name => 'HumanMethylation450',
                                                -vendor => 'ILLUMINA',
                                                -format  => 'METHYLATION',
                                                -type    => 'OLIGO',
                                                #-description => '',
                                                -class   => 'ILLUMINA_INFINIUM',
                                                #-version => '1.1',
                                               },
					  
                     },
	 
     INPUT_FORMAT => 'FASTA',
    },



    #CODELINK

    IMPORT_CODELINK_ARRAYS => 
    {
     IIDREGEXP => '^>(\S+):(\S+).*$',
	 
     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,

                     #-probe_set   => 2,#This could be annotation
                    },
	 	 
     ARRAY_PARAMS => {
					  
                      'CODELINK' => {
                                     -name => 'CODELINK',
                                     -vendor => 'CODELINK',
                                     #-setsize => undef,
                                     -format  => 'EXPRESSION',
                                     -type    => 'OLIGO',
                                     #-description => '',
                                     -class   => 'CODELINK',
                                    },
					  
                     },
	 
     INPUT_FORMAT => 'FASTA',
    },

    #AGILENT

    IMPORT_AGILENT_ARRAYS => 
    {
	 
     IIDREGEXP => '^>(\S+):(\S+)\s*(.*)$',
     #IIDREGEXP => '^>(\S+):(.+)', #EG HACK
	 
     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     -description => 2,
                     #-probe_set   => 2,#This could be annotation
                    },
	 	 
     ARRAY_PARAMS => 
     {
					  
					  
      #Danio
      'G2518A' => {
                   -name => 'G2518A',
                   -vendor => 'AGILENT',
                   #-setsize => undef,
                   -format  => 'EXPRESSION',
                   -type    => 'OLIGO',
                   #-description => '',
                   -class   => 'AGILENT',
                  },

      'G2519F' => {
                   -name => 'G2519F',
                   -vendor => 'AGILENT',
                   #-setsize => undef,
                   -format  => 'EXPRESSION',
                   -type    => 'OLIGO',
                   #-description => '',
                   -class   => 'AGILENT',
                  },



			  


      #plant - ara
      'G2519F-015059' => {
                          -name => 'G2519F-015059',
                          -vendor => 'AGILENT',
                          #-setsize => undef,
                          -format  => 'EXPRESSION',
                          -type    => 'OLIGO',
                          #-description => '',
                          -class   => 'AGILENT',
                         },
      'G2519F-021169' => {
                          -name => 'G2519F-021169',
                          -vendor => 'AGILENT',
                          #-setsize => undef,
                          -format  => 'EXPRESSION',
                          -type    => 'OLIGO',
                          #-description => '',
                          -class   => 'AGILENT',
                         },



      #plant - Rice 

      'G2519F-015241' => {
                          -name => 'G2519F-015241',
                          -vendor => 'AGILENT',
                          #-setsize => undef,
                          -format  => 'EXPRESSION',
                          -type    => 'OLIGO',
                          #-description => '',
                          -class   => 'AGILENT',
                         },
      'G4138A-012106' => {
                          -name => 'G4138A-012106',
                          -vendor => 'AGILENT',
                          #-setsize => undef,
                          -format  => 'EXPRESSION',
                          -type    => 'OLIGO',
                          #-description => '',
                          -class   => 'AGILENT',
                         },


					

				
      #human/mouse/rat
      'WholeGenome_4x44k_v1' => {
                                 -name => 'WholeGenome_4x44k_v1',
                                 -vendor => 'AGILENT',
                                 #-setsize => undef,
                                 -format  => 'EXPRESSION',
                                 -type    => 'OLIGO',
                                 #-description => '',
                                 -class   => 'AGILENT',	
                                },

      'WholeGenome_4x44k_v2' => {
                                 -name => 'WholeGenome_4x44k_v2',
                                 -vendor => 'AGILENT',
                                 #-setsize => undef,
                                 -format  => 'EXPRESSION',
                                 -type    => 'OLIGO',
                                 #-description => '',
                                 -class   => 'AGILENT',	
                                },


      'SurePrint_G3_GE_8x60k' => {
                                  -name => 'SurePrint_G3_GE_8x60k',
                                  -vendor => 'AGILENT',
                                  #-setsize => undef,
                                  -format  => 'EXPRESSION',
                                  -type    => 'OLIGO',
                                  #-description => '',
                                  -class   => 'AGILENT',	
                                 },
					  
                    
      #Rabbit only

      'SurePrint_G2519F_4x44k' => {
                                   -name => 'SurePrint_G2519F_4x44k',
                                   -vendor => 'AGILENT',
                                  #-setsize => undef,
                                  -format  => 'EXPRESSION',
                                  -type    => 'OLIGO',
                                  #-description => '',
                                  -class   => 'AGILENT',	
                                  },

      #Rat only
      'WholeGenome_4x44k_v3' => {
                                 -name => 'WholeGenome_4x44k_v3',
                                 -vendor => 'AGILENT',
                                 #-setsize => undef,
                                 -format  => 'EXPRESSION',
                                 -type    => 'OLIGO',
                                 #-description => '',
                                 -class   => 'AGILENT',	
                                },

      #human
      'CGH_44b' => {
                    -name => 'CGH_44b',
                    -vendor => 'AGILENT',
                    #-setsize => undef,
                    -format  => 'CGH',
                    -type    => 'OLIGO',
                    #-description => '',
                    -class   => 'AGILENT',	
                   },

      #Celegans
      #This is mixed array and array_chip config
      'OligoArray_012795' => {
                              -name => 'OligoArray',
                              -vendor => 'AGILENT',
                              #-setsize => undef,
                              -format  => 'EXPRESSION',
                              -type    => 'OLIGO',
                              #-description => '',
                              -class   => 'AGILENT',
                              -design_id => '012795', #array_chip maybe this shoudl just be integrated into the name?
                             },


     },
	 
     INPUT_FORMAT => 'FASTA',
    },

	




    #PHALANX
    #Human
    #ftp://ftp.phalanxbiotech.com/pub/probe_sequences/hoa
    #Mouse
    #ftp://ftp.phalanxbiotech.com/pub/probe_sequences/moa
    IMPORT_PHALANX_ARRAYS => 
    {
     IIDREGEXP => '^>(\S+):(\S+).*$',
	 
     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     #-probe_set   => 2,#This could be annotation
                    },
	 	 
     ARRAY_PARAMS => {
					  
                      'OneArray' => {
                                     -name => 'OneArray',
                                     -vendor => 'PHALANX',
                                     #-setsize => undef,
                                     -format  => 'EXPRESSION',
                                     -type    => 'OLIGO',
                                     #-description => '',
                                     -class   => 'PHALANX',
                                    },
					  
                     },
	 
     INPUT_FORMAT => 'FASTA',
    },


    #LEIDEN

    IMPORT_LEIDEN_ARRAYS => 
    {
     IIDREGEXP => '^>(\S+):(\S+)\s*(.*)$',
	 
     IFIELDORDER => {
                     -name        => 1,
                     -array_chip  => 0,
                     -array       => 0,
                     -description => 2,
                     #-probe_set   => 2,#This could be annotation
                    },
	 	 
     ARRAY_PARAMS => {
                      #Danio
                      'LEIDEN2' => {
                                    -name => 'LEIDEN2',
                                    -vendor => 'LEIDEN',
                                    #-setsize => undef,
                                    -format  => 'EXPRESSION',
                                    -type    => 'OLIGO',
                                    #-description => '',
                                    -class   => 'LEIDEN',
                                   },

					  
                      'LEIDEN3' => {
                                    -name => 'LEIDEN3',
                                    -vendor => 'LEIDEN',
                                    #-setsize => undef,
                                    -format  => 'EXPRESSION',
                                    -type    => 'OLIGO',
                                    #-description => '',
                                    -class   => 'LEIDEN',
                                   },

					  
                     },
	 
     INPUT_FORMAT => 'FASTA',
    },

    #STEMPLE

    IMPORT_STEMPLE_LAB_SANGER_ARRAYS => 
    {
     IIDREGEXP => '^>(\S+):(\S+)\s*(.*)$', #Need to add desc field here
	 
     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     -description   => 2,
                    },
	 	 
     ARRAY_PARAMS => {
                      #Danio
                      'MattArray1' => {
                                       -name => 'MattArray1',
                                       -vendor => 'STEMPLE_LAB_SANGER',
                                       #-setsize => undef,
                                       -format  => 'EXPRESSION',
                                       -type    => 'OLIGO',
                                       #-description => '',
                                       -class   => 'STEMPLE_LAB_SANGER',
                                      },
					  
					  
                      'MattArray2' => {
                                       -name => 'MattArray2',
                                       -vendor => 'STEMPLE_LAB_SANGER',
                                       #-setsize => undef,
                                       -format  => 'EXPRESSION',
                                       -type    => 'OLIGO',
                                       #-description => '',
                                       -class   => 'STEMPLE_LAB_SANGER',
                                      },
					  
					  
                     },
	 
     INPUT_FORMAT => 'FASTA',
    },

    IMPORT_CATMA_ARRAYS =>
    {
     IIDREGEXP => '^>(\S+):(\S+)', #Need to add desc field here

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                    },

     ARRAY_PARAMS => {
                      'CATMA' => {
                                  -name => 'CATMA',
                                  -vendor => 'CATMA',
                                  #-setsize => undef,
                                  -format  => 'EXPRESSION',
                                  -type    => 'OLIGO',
                                  #-description => '',
                                  -class   => 'CATMA',
                                 },
                     },
     INPUT_FORMAT => 'FASTA',
    },





    IMPORT_NSF_ARRAYS =>
    {
     IIDREGEXP => '^>(\S+):(\S+)', #Need to add desc field here

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                    },

     ARRAY_PARAMS => {

                      'BGIYale' => {
                                    -name => 'BGIYale',
                                    -vendor => 'NSF',
                                    #-setsize => undef,
                                    -format  => 'EXPRESSION',
                                    -type    => 'OLIGO',
                                    #-description => '',
                                    -class   => 'NSF',
                                   },


                      'NSF20K' => {
                                   -name => 'NSF20K',
                                   -vendor => 'NSF',
                                   #-setsize => undef,
                                   -format  => 'EXPRESSION',
                                   -type    => 'OLIGO',
                                   #-description => '',
                                   -class   => 'NSF',
                                  },

                      'NSF45K' => {
                                   -name => 'NSF45K',
                                   -vendor => 'NSF',
                                   #-setsize => undef,
                                   -format  => 'EXPRESSION',
                                   -type    => 'OLIGO',
                                   #-description => '',
                                   -class   => 'NSF',
                                  },




                     },
     INPUT_FORMAT => 'FASTA',
    },

   }
  );

sub import {
  my ($callpack) = caller(0);   # Name of the calling package
  my $pack = shift;             # Need to move package off @_

  # Get list of variables supplied, or else everything
  my @vars = @_ ? @_ : keys( %Config );
  return unless @vars;
  
  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
  die $@ if $@;


  foreach (@vars) {
    if ( defined $Config{$_} ) {
      no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$Config{ $_ };
    } else {
	    die "Error: Config: $_ not known\n";
    }
  }
}

1;
