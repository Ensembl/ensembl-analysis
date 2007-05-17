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

### want to make the config inherit from Funcgen config that then 
### contains the experiment info!!!

package Bio::EnsEMBL::Analysis::Config::Funcgen::Nessie;

use strict;
use vars qw(%Config);

my %ftype = (
             'H3K4me1' => ['Histone 3 Lysine 4 Mono-Methylation','HISTONE'],
             'H3K4me2' => ['Histone 3 Lysine 4 Di-Methylation','HISTONE'],
             'H3K4me3' => ['Histone 3 Lysine 4 Tri-Methylation','HISTONE'],
             'H3ac'    => ['Histone 3 Acethylation','HISTONE'],
             'H3K9ac'  => ['Histone 3 Lysine 9 Acethylation','HISTONE'],
             'H4ac'    => ['Histone 4 Acethylation','HISTONE'],
             'CTCF'    => ['CCCTC-binding factor sites','INSULATOR']
             );
my %ctype = (
             'HeLa'    => 'Human Epithelial Carcinoma Cells',
             'GM06990' => 'Human B-Lymphocyte Cells',
             'U2OS'    => 'Human Bone Osteosarcoma Epithelial Cells',
             'IMR90'   => 'Human Primary Fibroblast Cells',
             'U937'    => 'Human Leukemic Monocyte Lymphoma Cells'
             );


%Config = (
           NESSIE_CONFIG => {
           DEFAULT => {
               EFG_EXPERIMENT => 'Experiment Name',
               LOGIC_NAME     => 'Nessie',
               PROGRAM        => 'nessie',
               OPTIONS        => 
                    # PCR array based data
                    '--distribution="normal" --fit-method="full" '.
                    # NimbleGen data
                    # (recommended for standard ChIP-chip fragment length ~300-600)
                    #'--distribution="weibull" --fit-method="slope" ',
                    # Affymetrix data  and NimbleGen
                    # (use with very short fragments hybridised to arrays ~200bp)
                    #'--distribution="normal" --fit-method="slope",
                    # Number of replicates per factor
                    '--replicates=6 '.
                    # One step trainig and annotation run  
                    '--onestep="train_state.para"'
               },
           NESSIE => {
               EFG_EXPERIMENT => {
                   # Sanger Encode
                   'H3K4me1-HeLa' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H3K4me1',
                       FT_CLASS   => $ftype{'H3K4me1'}[1],
                       FT_DESC    => $ftype{'H3K4me1'}[0],
                       CT_NAME    => 'HeLa',
                       CT_DESC    => $ctype{'HeLa'}},
                   'H3K4me2-HeLa' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H3K4me2',
                       FT_CLASS   => $ftype{'H3K4me2'}[1],
                       FT_DESC    => $ftype{'H3K4me2'}[0],
                       CT_NAME    => 'HeLa',
                       CT_DESC    => $ctype{'HeLa'}},
                   'H3K4me3-HeLa' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H3K4me3',
                       FT_CLASS   => $ftype{'H3K4me3'}[1],
                       FT_DESC    => $ftype{'H3K4me3'}[0],
                       CT_NAME    => 'HeLa',
                       CT_DESC    => $ctype{'HeLa'}},
                   'H3ac-HeLa' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H3ac',
                       FT_CLASS   => $ftype{'H3ac'}[1],
                       FT_DESC    => $ftype{'H3ac'}[0],
                       CT_NAME    => 'HeLa',
                       CT_DESC    => $ctype{'HeLa'}},
                   'H4ac-HeLa' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H4ac',
                       FT_CLASS   => $ftype{'H4ac'}[1],
                       FT_DESC    => $ftype{'H4ac'}[0],
                       CT_NAME    => 'HeLa',
                       CT_DESC    => $ctype{'HeLa'}},
                   'H3K4me1-GM06990' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H3K4me1',
                       FT_CLASS   => $ftype{'H3K4me1'}[1],
                       FT_DESC    => $ftype{'H3K4me1'}[0],
                       CT_NAME    => 'GM06990',
                       CT_DESC    => $ctype{'GM06990'}},
                   'H3K4me2-GM06990' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H3K4me2',
                       FT_CLASS   => $ftype{'H3K4me2'}[1],
                       FT_DESC    => $ftype{'H3K4me2'}[0],
                       CT_NAME    => 'GM06990',
                       CT_DESC    => $ctype{'GM06990'}},
                   'H3K4me3-GM06990' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H3K4me3',
                       FT_CLASS   => $ftype{'H3K4me3'}[1],
                       FT_DESC    => $ftype{'H3K4me3'}[0],
                       CT_NAME    => 'GM06990',
                       CT_DESC    => $ctype{'GM06990'}},
                   'H3ac-GM06990' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H3ac',
                       FT_CLASS   => $ftype{'H3ac'}[1],
                       FT_DESC    => $ftype{'H3ac'}[0],
                       CT_NAME    => 'GM06990',
                       CT_DESC    => $ctype{'GM06990'}},
                   'H4ac-GM06990' => {
                       ANALYSIS   => 'SangerPCR',
                       FT_NAME    => 'H4ac',
                       FT_CLASS   => $ftype{'H4ac'}[1],
                       FT_DESC    => $ftype{'H4ac'}[0],
                       CT_NAME    => 'GM06990',
                       CT_DESC    => $ctype{'GM06990'}},

                   # NCMLS
                   'Stunnenberg_all_OID_1963' => {
                       ANALYSIS   => 'VSN_GLOG',
                       FT_NAME    => 'H3K9ac',
                       FT_CLASS   => $ftype{'H3K9ac'}[1],
                       FT_DESC    => $ftype{'H3K9ac'}[0],
                       CT_NAME    => 'U2OS',
                       CT_DESC    => $ctype{'U2OS'}},

                   # CTCF (Bing Ren)
                   'ctcf_ren' => {
                       ANALYSIS   => 'VSN_GLOG',
                       FT_NAME    => 'CTCF',
                       FT_CLASS   => $ftype{'CTCF'}[1],
                       FT_DESC    => $ftype{'CTCF'}[0],
                       CT_NAME    => 'IMR90',
                       CT_DESC    => $ctype{'IMR90'}},
                   
                   #FAIRE (Duke University)
                   #EXPERIMENT => 'HeLaS3_FAIRE_WholeGenome',
 
               }               
           }
       }
           );

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
