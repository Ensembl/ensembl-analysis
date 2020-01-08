=head1 NAME

    Bio::EnsEMBL::Analysis::Config::SelenoBuild

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::SelenoBuild;

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
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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


package Bio::EnsEMBL::Analysis::Config::SelenoBuild;

use warnings ;
use strict;
use vars qw( %Config );

# Hash containing config info
%Config = (
  EXONERATE_CONFIG_BY_LOGIC => {
      DEFAULT => {

        # The GENOMICSEQS can be a file (string, A.), a
        # directory with files (string, B.) or an anonymous
        # array containing the names of multiple directories
        # (array of strings, C.). In the latter case the order
        # of appearance will determine the usage of directories
        # and files; new versions should therefore be listed in
        # the first directory.
        #
        # format: A. '/your/soft/masked/genome.fa'
        #      or B. '/your/soft/masked/genomeDIR'
        #      or C. ['/your/soft/masked/genomeDIR_1', '/your/soft/masked/genomeDIR_2']
        GENOMICSEQS         => '/your/soft/masked/genome.fa',
        QUERYTYPE           => undef,
        QUERYSEQS           => undef,
        QUERYANNOTATION     => '/annotation/file/that/will/be/create/during/exonerateForSeleno/run',
        IIDREGEXP           => undef,
        OUTDB               => undef,
        FILTER              => undef,
        COVERAGE_BY_ALIGNED => undef,
        OPTIONS             => undef,

        # Comment out or set to undef if non-reference regions (eg DR52
        # in human) should NOT be fetched, ie should be ignored:
        NONREF_REGIONS => 1,

        # If program not defined, will look in program_file of analysis
        # table. Or will take default 0.8.3 if neither is defined.
        PROGRAM => '/usr/local/ensembl/bin/exonerate-0.9.0',    # /usr/local/ensembl/bin/exonerate
      },

      exonerate_for_seleno_logicname => {
        GENOMICSEQS => '/data/blastdb/Ensembl/Human/GRCh37/genome/softmasked/softmasked_dusted.fa',
        PROGRAM => '/usr/local/ensembl/bin/exonerate-0.9.0',
        # GENOMICSEQS obtained from DEFAULT
        # IIDREGEXP not set; input ids are file names

        QUERYTYPE => 'dna',
        QUERYSEQS => '/directory/with/seleno/sequences/in/embl/format',

        QUERYANNOTATION     => '/annotation/file/that/will/be/create/during/exonerateForSeleno/run',

        # CHANGE! You can now either submit all db-connection parameters
        # for OUTDB (old method ) - or you can just supply a string. The string
        # must be assigned to a database-connection in Databases.pm
        OUTDB => "SELENO_DB",    # HASH-key out of Databases.pm

        FILTER => {
          OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
          PARAMETERS => { -coverage                 => 50,
                          -percent_id               => 50,
                          -best_in_genome           => 1,
                          -reject_processed_pseudos => 1,
          },
        },

        COVERAGE_BY_ALIGNED => 1,
        OPTIONS             => "--model est2genome --forwardcoordinates FALSE "
          . "--softmasktarget TRUE --exhaustive FALSE --score 50 "
          . "--saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 12",
      },
      seleno_build_logicname => {
        GENOMICSEQS => '/data/blastdb/Ensembl/Human/GRCh37/genome/softmasked/softmasked_dusted.fa',
        PROGRAM => '/nfs/users/nfs_j/jb16/cvs_checkout/ensembl-personal/genebuilders/scripts/exonerate.hacked.cdna2genome',

        # GENOMICSEQS obtained from DEFAULT
        # IIDREGEXP not set; input ids are file names

        QUERYTYPE => 'dna',
        QUERYSEQS => '/directory/with/seleno/sequences/in/swissprot/format',

        QUERYANNOTATION     => '/annotation/file/that/was/created/during/exonerateForSeleno/run',

        # CHANGE! You can now either submit all db-connection parameters
        # for OUTDB (old method ) - or you can just supply a string. The string
        # must be assigned to a database-connection in Databases.pm
        OUTDB => "SELENO_DB",    # HASH-key out of Databases.pm

        FILTER => {
          OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
          PARAMETERS => { -coverage                 => 90,
                          -percent_id               => 90,
                          -best_in_genome           => 1,
                          -reject_processed_pseudos => 1,
          },
        },

        COVERAGE_BY_ALIGNED => 1,
        OPTIONS             => "--model cdna2genome --forwardcoordinates FALSE "
          . "--softmasktarget TRUE --exhaustive FALSE --score 50 "
          . "--saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 12",
      },


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
