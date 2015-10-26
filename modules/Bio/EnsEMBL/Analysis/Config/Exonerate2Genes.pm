=head1 NAME

    Bio::EnsEMBL::Analysis::Config::Exonerate2Genes

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;

=head1 DESCRIPTION

    It imports and sets a number of standard global variables into the
    calling package, which are used in many scripts in the human sequence
    analysis system. The variables are first declared using "use vars", so
    that it can be used when "use strict" is in use in the calling script.
    Without arguments all the standard variables are set, and with a list,
    only those variables whose names are provided are set. The module will
    die if a variable which doesn't appear in its C<%Exonerate> hash is
    asked to be set.

    Since The RunnableDB that this config controls can be used for
    inferring transcript structures from (different sets of) EST,
    cDNA and proteins, and several uses may be required in the same
    pipeline, this Config contains one primary config variable,
    EXONERATE_TRANSCRIPT_CONFIG. This is hash keyed off logic name,
    each entry of which is a hash containing the variable that affect
    the behaviour of the RunnableDB. When the RunnableDB instance is
    created, the correct entry is identified by logic name and value for a
    corresponding set of local variables are set.

=head1 LICENSE


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


=head1 CONTACT

    Please email comments or questions to the public Ensembl
    developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

    Questions may also be sent to the Ensembl help desk at
    <http://www.ensembl.org/Help/Contact>.

=cut


package Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;

use strict;
use vars qw( %Config );

# Hash containing config info
%Config = (
  EXONERATE_CONFIG_BY_LOGIC => {
    DEFAULT => {
      COVERAGE_BY_ALIGNED => undef,
      FILTER => undef,
      GENOMICSEQS => '/your/soft/masked/genome.fa',
      IIDREGEXP => undef,
      KILL_TYPE => undef,
      NONREF_REGIONS => 1,
      OPTIONS => undef,
      OUTDB => undef,
      PROGRAM => undef,
      QUERYANNOTATION => undef,
      QUERYSEQS => undef,
      QUERYTYPE => undef,
      SEQFETCHER_OBJECT => 'Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher',
      SEQFETCHER_PARAMS => {
        -db => [
          '/path/to/dir/containing/index/'
        ]
      },
      SOFT_MASKED_REPEATS => [
        'repeatmask',
        'dust'
      ],
      USE_KILL_LIST => 0
    },
    CDNA_logic_name => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 1,
          -coverage => 90,
          -percent_id => 97,
          -reject_processed_pseudos => 1
        }
      },
      KILL_TYPE => undef,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      OUTDB => 'ESTCDNA_DB',
      QUERYSEQS => '/directory/with/cdna/chunks',
      QUERYTYPE => 'dna',
      USE_KILL_LIST => 0
    },
    EST_logicname => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 1,
          -coverage => 90,
          -percent_id => 97,
          -reject_processed_pseudos => 1
        }
      },
      KILL_TYPE => undef,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      OUTDB => 'ESTCDNA_DB',
      QUERYSEQS => '/directory/with/EST/chunks',
      QUERYTYPE => 'dna',
      USE_KILL_LIST => 0
    },
    Protein_logicname => {
      COVERAGE_BY_ALIGNED => 0,
      IIDREGEXP => '(\d+):(\d+)',
      KILL_TYPE => undef,
      OPTIONS => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE  --bestn 1',
      OUTDB => 'EXONERATE_DB',
      QUERYSEQS => '/file/containing/all/proteins.fa',
      QUERYTYPE => 'protein',
      USE_KILL_LIST => 0
    },
    cdna_update => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 1,
          -coverage => 90,
          -percent_id => 97,
          -reject_processed_pseudos => 1,
          -verbosity => 0
        }
      },
      GENOMICSEQS => '/data/blastdb/Ensembl/human/GRCh38/genome/softmasked_dusted/toplevel.with_nonref_and_GRCh38_p5.no_duplicate.softmasked_dusted.fa',
      KILL_TYPE => undef,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      OUTDB => 'ESTCDNA_DB',
      QUERYSEQS => '/lustre/scratch109/sanger/cgg/cdna_update_human83/chunks',
      QUERYTYPE => 'dna',
      USE_KILL_LIST => 0
    },
    cdna_update_2 => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 1,
          -coverage => 90,
          -percent_id => 97,
          -reject_processed_pseudos => 1,
          -verbosity => 1
        }
      },
      GENOMICSEQS => '/data/blastdb/Ensembl/human/GRCh38/genome/softmasked_dusted/toplevel.with_nonref_and_GRCh38_p5.no_duplicate.softmasked_dusted.fa',
      KILL_TYPE => undef,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --bestn 10 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      OUTDB => 'ESTCDNA_DB',
      QUERYSEQS => '/lustre/scratch109/sanger/cgg/cdna_update_human83/chunks2',
      QUERYTYPE => 'dna',
      USE_KILL_LIST => 0
    },
    cdna_update_3 => {
      COVERAGE_BY_ALIGNED => 1,
      FILTER => {
        OBJECT => 'Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter',
        PARAMETERS => {
          -best_in_genome => 1,
          -coverage => 90,
          -percent_id => 97,
          -reject_processed_pseudos => 1,
          -verbosity => 1
        }
      },
      GENOMICSEQS => '/data/blastdb/Ensembl/human/GRCh38/genome/softmasked_dusted/toplevel.with_nonref_and_GRCh38_p5.no_duplicate.softmasked_dusted.fa',
      KILL_TYPE => undef,
      OPTIONS => '--model est2genome --forwardcoordinates FALSE --maxintron 400000 --bestn 10 --softmasktarget FALSE --exhaustive FALSE  --score 500 --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      OUTDB => 'ESTCDNA_DB',
      QUERYSEQS => '/lustre/scratch109/sanger/cgg/cdna_update_human83/chunks3',
      QUERYTYPE => 'dna',
      USE_KILL_LIST => 0
    },
    orthologue_recovery => {
      COVERAGE_BY_ALIGNED => 0,
      IIDREGEXP => '(\d+):(\d+)',
      OPTIONS => '--model protein2genome --forwardcoordinates FALSE  --softmasktarget TRUE --exhaustive FALSE  --bestn 1',
      OUTDB => 'EXONERATE_DB',
      QUERYSEQS => '/path/to/where/you/want/your/seq/dumped/by/OrthologueAnalyis',
      QUERYTYPE => 'protein'
    }
  }
);

sub import {
  my ($callpack) = caller(0);    # Name of the calling package
  my $pack = shift;              # Need to move package off @_

  # Get list of variables supplied, or else everything
  my @vars = @_ ? @_ : keys(%Config);
  return unless @vars;

  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join( ' ', map { '$' . $_ } @vars ) . ")";
  die $@ if $@;

  foreach (@vars) {
    if ( defined $Config{$_} ) {
      no strict 'refs';
      # Exporter does a similar job to the following
      # statement, but for function names, not
      # scalar variables:
      *{"${callpack}::$_"} = \$Config{$_};
    } else {
      die "Error: Config: $_ not known\n";
    }
  }
} ## end sub import

1;
