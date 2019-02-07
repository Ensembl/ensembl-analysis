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

Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key) = @_;

  my %config = (
    default => {
      IIDREGEXP           => '(\d+):(\d+)',
      OPTIONS             => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'dna',
      GENOMICSEQS         => '#genome_file#',
      PROGRAM             => '#exonerate_path#',
      SOFT_MASKED_REPEATS => '#repeat_libraries#', # This should be an arrayref
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
    },
    exonerate_cdna => {
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1',
      FILTER  => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage => 50,
          -percent_id => 50,
        },
      },
    },

    exonerate_projection_dna => {
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1 --maxintron 100000 --forcegtag ',
      FILTER  => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage => 50,
          -percent_id => 50,
        },
      },
    },

    exonerate_projection_pseudogene => {
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE'.
                 ' --bestn 5 --maxintron 100000 --minintron 1 --frameshift 10',
      FILTER  => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage => 50,
          -percent_id => 50,
        },
      },
    },

    exonerate_projection_ig_tr_protein => {
      IIDREGEXP => '(\d+):(\d+)',
      OPTIONS   => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1 --maxintron 50000',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'protein',
    },


    exonerate_projection_coding => {
        COVERAGE_BY_ALIGNED => 1,
        OPTIONS => "--model cdna2genome --forwardcoordinates FALSE ".
        "--softmasktarget TRUE --exhaustive FALSE ".
        "--refine region --refineboundary 5000 --forcegtag 1 ".
        "--score 500 --saturatethreshold 100 ".
        "--dnahspthreshold 60 --dnawordlen 15 --bestn 10",
        FILTER => { OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
          PARAMETERS => { -coverage => 50,
            -percent_id => 50,
            -best_in_genome => 0,
            -reject_processed_pseudos => 1,
          },
        },
      QUERYTYPE           => 'dna',
    },

    exonerate_protein => {
      IIDREGEXP => '(\d+):(\d+)',
      OPTIONS   => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'protein',
    },

    exonerate_protein_pseudo => {
      IIDREGEXP => '(\d+):(\d+)',
      OPTIONS   => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1 --maxintron 100000 --minintron 1 --frameshift 10',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'protein',
      FILTER  => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage => 70,
          -percent_id => 25,
        },
      },
    },

    exonerate_protein_recover => {
      IIDREGEXP => '(\d+):(\d+)',
      OPTIONS   => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'protein',
    },

    exonerate_protein_human_patch => {
      IIDREGEXP => '(\d+):(\d+)',
      OPTIONS   => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive TRUE --bestn 5 --maxintron 50000 --minintron 20',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'protein',
      FILTER => {
                  OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                  PARAMETERS                  => {
                    -coverage   => '#exonerate_cov#',
                    -percent_id => '#exonerate_pid#',
                    -best_in_genome           => 1,
                    -reject_processed_pseudos => 1,
                  },
      },
    },

    exonerate_protein_human_patch_non_exhaustive => {
      IIDREGEXP => '(\d+):(\d+)',
      OPTIONS   => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 5 --maxintron 50000 --minintron 20',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'protein',
      FILTER => {
                  OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                  PARAMETERS                  => {
                    -coverage   => '#exonerate_cov#',
                    -percent_id => '#exonerate_pid#',
                    -best_in_genome           => 1,
                    -reject_processed_pseudos => 1,
                  },
      },
    },

    exonerate_cdnaupdate => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 0,
          -coverage => '#exonerate_cov#',
          -percent_id => '#exonerate_pid#',
          -reject_processed_pseudos => 1,
          -verbosity => 2,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 100000 --softmasktarget TRUE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    exonerate_cdnaupdate_loose => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 10,
          -coverage => '#exonerate_cov#',
          -percent_id => '#exonerate_pid#',
          -reject_processed_pseudos => 1,
          -verbosity => 2,
        }
      },
      KILL_TYPE => undef,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --bestn 10 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    exonerate => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 0,
          -coverage => 90,
          -percent_id => 97,
          -reject_processed_pseudos => 1,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    exonerate_2 => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 10,
          -coverage => 90,
          -percent_id => 97,
          -reject_processed_pseudos => 1,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    exonerate_short_intron => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 0,
          -coverage => 50,
          -percent_id => 50,
          -reject_processed_pseudos => 1,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 100000 --softmasktarget TRUE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    exonerate_short_intron_lose => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 0,
          -coverage => 50,
          -percent_id => 50,
          -reject_processed_pseudos => 1,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 200000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    cdna_selection => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 0,
          -coverage => 90,
          -percent_id => 97,
          -reject_processed_pseudos => 0,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 200000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    gene_selection => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 0,
          -coverage => 90,
          -percent_id => 95,
          -reject_processed_pseudos => 1,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 200000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    pacbio_selection => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::PacBioTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 0,
          -coverage => 90,
          -percent_id => 95,
          -reject_processed_pseudos => 1,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 200000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    pacbio_selection_bestn => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::PacBioTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 1,
          -coverage => 90,
          -percent_id => 95,
          -reject_processed_pseudos => 1,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 200000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    cdna2genome => {
        COVERAGE_BY_ALIGNED => 1,
        OPTIONS => '--model cdna2genome --forwardcoordinates FALSE '.
                   '--softmasktarget TRUE --exhaustive FALSE '.
                   '--score 500 --saturatethreshold 100 '.
                   '--dnawordlen 15 --codonwordlen 15 '.
                   '--dnahspthreshold 60 --bestn 10',
        FILTER => {
          OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::Filter::cDNA2GenomeTranscriptFilter',
          PARAMETERS => {
            -coverage => 50,
            -percent_id => 90,
            -best_in_genome => 1,
            -reject_processed_pseudos => 1,
          },
        },
    },
    exonerate_cov_per_sub => {
      IIDREGEXP           => '(\d+):(\d+)',
      OPTIONS             => ' --maxintron 100000 --model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE ',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'dna',
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage   => '#exonerate_cdna_cov#',
          -percent_id => '#exonerate_cdna_pid#',
        },
      },
    },
    cdna_est2genome => {
      IIDREGEXP           => '(\d+):(\d+)',
      OPTIONS             => ' --model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14 --score 500',
      COVERAGE_BY_ALIGNED => 1,
      QUERYTYPE           => 'dna',
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage   => '#exonerate_cdna_cov#',
          -percent_id => '#exonerate_cdna_pid#',
          -best_in_genome => 1,
          -reject_processed_pseudos => 1,
        },
      },
    },
    exonerate_cov_per_bestn_sub => {
      IIDREGEXP           => '(\d+):(\d+)',
      OPTIONS             => ' --model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn #exonerate_bestn#',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'dna',
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage   => '#exonerate_cdna_cov#',
          -percent_id => '#exonerate_cdna_pid#',
        },
      },
    },
    exonerate_cov_per_bestn_loose_sub => {
      IIDREGEXP           => '(\d+):(\d+)',
      OPTIONS             => ' --maxintron 400000 --model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn #exonerate_bestn#',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'dna',
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage   => '#exonerate_cdna_cov#',
          -percent_id => '#exonerate_cdna_pid#',
        },
      },
    },
    pacbio_exonerate => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 0,
          -coverage => 90,
          -percent_id => 95,
          -reject_processed_pseudos => 1,
          -verbosity => 1,
        }
      },
      KILL_TYPE => undef,
      USE_KILL_LIST => 0,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 200000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
    },
    protein_cov_per_bestn_sub => {
      OPTIONS             => ' --model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn #exonerate_bestn#',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'protein',
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage   => '#exonerate_cdna_cov#',
          -percent_id => '#exonerate_cdna_pid#',
        },
      },
    },
    protein_cov_per_bestn_maxintron_sub => {
      OPTIONS             => ' --model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn #exonerate_bestn# --maxintron #exonerate_max_intron#',
      COVERAGE_BY_ALIGNED => 0,
      QUERYTYPE           => 'protein',
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -coverage   => '#exonerate_cdna_cov#',
          -percent_id => '#exonerate_cdna_pid#',
          -best_in_genome => 0,
          -reject_processed_pseudos => 1,
        },
      },
    },
# This is coming from the IgSegBuilder configs: ensembl-config/human/GRCh37_ig/Bio/EnsEMBL/Analysis/Config/Exonerate2Genes.pm
# Thibaut have modified it for the Hive pipeline
    c_segment => {
      QUERYTYPE => 'protein',
      OPTIONS  => "--model protein2genome --softmasktarget FALSE --maxintron 5000 --forcegtag TRUE --bestn 5 ",
      FILTER => {
        OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::IgTranscriptFilter',
        PARAMETERS => {
          -minpercent   => 95,
          -percentrange => 0,
        },
      },
    },

    d_segment =>  {
      QUERYTYPE => 'protein',
      OPTIONS  => "--model protein2genome --softmasktarget TRUE --maxintron 0 --gappedextension TRUE --bestn 1 ",
      FILTER => {
        OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::IgTranscriptFilter',
        PARAMETERS => {
          -minpercent  =>  95,
          -percentrange => 0,
        },
      },
    },

    j_segment =>  {
      QUERYTYPE => 'protein',
      OPTIONS  => "--model protein2genome --softmasktarget TRUE --maxintron 0 --percent 90 --score 0 --bestn 5 ",
      FILTER => {
        OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::IgTranscriptFilter',
        PARAMETERS => {
          -minpercent   => 95,
          -percentrange => 0,
        },
      },
    },

    v_segment => {
      QUERYTYPE => 'protein',
      # minintron set artificially low to 20 due to bug in exonerate;
      # softmasktarget set to FALSE because some leader exons are low-complexity
      OPTIONS  => "--model protein2genome --softmasktarget FALSE --maxintron 1000 --minintron 0 --bestn 5 ",
      FILTER => {
        OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::IgTranscriptFilter',
        PARAMETERS => {
          -minpercent   => 95,
          -percentrange => 0,
        },
      },
    },
  );
  return $config{$key};
}

1;
