# Ensembl module for Bio::EnsEMBL::Analysis::Config::Funcgen::Nessie
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Config::Funcgen::Nessie

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Config::Funcgen::Nessie;
  
  use Bio::EnsEMBL::Analysis::Config::Funcgen::Nessie qw(CONFIG);

=head1 DESCRIPTION

This is a module needed to provide configuration for the
Nessie RunnableDBs.

NESSIE_CONFIG is an hash of hashes which contains analysis specific
settings and is keyed on logic_name

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Config::Funcgen::Nessie;

use strict;
use vars qw(%Config);

%Config = 
    (
     NESSIE_CONFIG => {
         DEFAULT => {
             PROGRAM           => $ENV{NE_PROGRAM},
             PROGRAM_FILE      => $ENV{NE_PROGRAM_FILE},
             VERSION           => $ENV{NE_VERSION},
             PARAMETERS        => $ENV{NE_PARAMETERS},
             EXPERIMENT        => $ENV{EXPERIMENT},
             NORM_ANALYSIS     => $ENV{NORM_ANALYSIS},
             RESULT_SET_REGEXP => $ENV{RESULT_SET_REGEXP},
             DATASET_NAME      => $ENV{DATASET_NAME},
             ANALYSIS_WORK_DIR => $ENV{ANALYSIS_WORK_DIR},
         },
         NESSIE => {
             PROGRAM           => $ENV{NE_PROGRAM},
         },
         NESSIE_NG => {
             PROGRAM           => $ENV{NE_PROGRAM},
         },
#         NESSIE_PCRARRAY => {
#             VERSION => '0.13',
#             RESULT_SET_REGEXP => '_BR\d+_TR\d+',
#             PARAMETERS =>
#                 # PCR array based data
#                 ' --distribution="normal" --fit-method="full" '.
#                 # Number of replicates per factor
#                 ' --replicates=6 '.
#                 # One step trainig and annotation run  
#                 ' --onestep="/nfs/acari/graef/bin/train_state.para"',
#             },
#         NESSIE_NG_STD_2 => {
#             VERSION => '0.13',
#             RESULT_SET_REGEXP => '_BR\d+_TR\d+',
#             PARAMETERS =>
#                 # NimbleGen data
#                 # (recommended for standard ChIP-chip fragment length ~300-600)
#                 ' --distribution="weibull" --fit-method="slope" '.
#                 # Number of replicates per factor
#                 ' --replicates=2 '.
#                 # One step trainig and annotation run  
#                 ' --onestep="/nfs/acari/graef/bin/train_state.para"',
#         }
     });


sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_
    
    # Get list of variables supplied, or else all
    my @vars = @_ ? @_ : keys(%Config);
    return unless @vars;
    
    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
        . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;
    
    
    foreach (@vars) {
        if (defined $Config{ $_ }) {
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
