=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
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
  );
  return $config{$key};
}

1;
