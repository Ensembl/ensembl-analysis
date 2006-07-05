# Ensembl module fohr Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases
#
# Copyright (c) 2006 Ensembl
#
# POD documentation - main docs before the code
#

=head1 NAME

Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases 

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases ; 

=head1 DESCRIPTION

Databases.pm is the main configuration file which holds the different 
parameters (usernames, hosts, passwords, ports, database-names) to 
connect to different databases used in the Ensembl-Analysis pipeline. 

It imports and sets a number of standard global variables into the
calling package. Without arguments all the standard variables are set,
and with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Config> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%Config> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables. 

Databases is a pure ripoff of humConf written by James Gilbert.
humConf is based upon ideas from the standard perl Env environment
module.


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases ; 

use strict;
use vars qw(%Config);

%Config= (

           #
           # This config-file is currently only used by 
           #           TranscriptCoalescer.pm  


  DATABASES => { 
                  
                 # The REFERENCE_DB (formely known as GB_DB) holds sequence + repeats + features 
                 # from raw computes (e.g. ab-inito predictions, dna- or protein alignments ) 
                  
                 REFERENCE_DB => 
                                 { 
                                   -dbname => 'sw4_stick_ref',
                                   -host => 'ia64g',
                                   -port => '3306',
                                   -user => 'ensro',
                                   -pass => '',
                                  },
  
 
                 # The GENEWISE_DB holds genes made by FPC_TargettedGenewise or FPC_BlastMiniGenewise 
                 # (TGE_gw or similarity_genewise - genes ) 

                 GENEWISE_DB => 
                                 { 
                                   -dbname => 'sw4_stick_SimGW',
                                   -host => 'ia64f',
                                   -port => '3306',
                                   -user => 'ensro',
                                   -pass => '',
                                  },


                 # The EXONERATE_DB ( formerly GB_cDNA ) holds alignments to cDNA's or 
                 # EST's gene-structtures made by exonerate (Exonerate2Genes.pm) 
                 
                 EXONERATE_DB => 
                                 { 
                                   -dbname => 'sw4_stick_EST',
                                   -host => 'ia64f',
                                   -port => '3306',
                                   -user => 'ensro',
                                   -pass => '',
                                  },


                 # The BLESSED_DB (formerly GB_BLESSED) holds the 'blessed' gene-set ( if there is one ) 

                 BLESSED_DB => 
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
                                 
                   
                 # The UTR_DB (formerly GB_COMB) holds genes made by the UTR-addtion-run 
                 # Combine_Genewises_and_E2Gs.pm writes to UTR_DB 
                   
                 UTR_DB  => 
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
                                 
                   
                 # GENEBUILD_DB (formerly GB_FINALDB) is the Database where 
                 # GeneBuilder.pm writes it's results to 

                 GENEBUILD_DB =>
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
         
                        
                 # PSEUDO_DB holds the pseudo-genes
           
                 PSEUDO_DB => 
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
                                 
                   
                 # COALESCER_DB is the DB where TranscriptCoalescer writes it's results to
                   
                 COALESCER_DB =>
                                 { 
                                   -dbname => 'sw4_stick_EST_GENE2',
                                   -host => 'ecs2',
                                   -port => '3362',
                                   -user => 'ensadmin',
                                   -pass => 'ensembl',
                                  },
                                 
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
