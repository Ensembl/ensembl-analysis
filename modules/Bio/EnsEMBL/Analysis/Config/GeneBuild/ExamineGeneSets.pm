#
# Ensembl configuration file used in 
#
# Bio::EnsEMBL::Analysis::RunnableDB::OrthologueEvaluator
#
# Copyright (c) 2006 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator

=head1 SYNOPSIS

    Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator

=head1 DESCRIPTION


OrthologueEvaluator - Configuration 

This is the main configuration file for OrthologueEvaluator, a perl 
module which uses information from an Ensembl Compara database to 
compare and assess gene predictions. 

The parameters to connect to various databases are defiend in 

  - modules/Bio/EnsEMBL/Analysis/Config/GeneBuild/Databases.pm 
  - modules/Bio/EnsEMBL/Analysis/Config/OrthologueEvaluator.pm 

The general function of this config file is to import  a number of 
standard global variables into the calling package. Without arguments 
all the standard variables are set, and with a list, only those variables 
whose names are provided are set.  The module will die if a variable 
which doesn\'t appear in its C<%Config> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%Config> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Analysis::Config::GeneBuild::ExamineGeneSets; 

use strict;
use vars qw(%Config);

%Config= 
 ( 

   #
   # GENE_DBNAME points to the DB of the species to investigate - this DB
   # should be the same as the one which was used to do the compara run  
   # and is defined in Analysis/Config/Genebuild/Databases.pm 
   #
   GENE_DBNAME => "REFERENCE_DB" ,    
  # GENE_DBNAME => "ORTHOLOGUES_DB" ,    

   # databases where to write recovered genes to ( connection details are defined in Databases.pm) 
  
   # if DO_NOT_READ_THE_EXONERATE2GENES_CONFIG_FILE == 1 the Exonerate2Genes.pm config file will be 
   # ignored (not read ) and the default parameters (hard-coded in the RunnableDB ) will be used 
   # to run Exonerat2Genes will be used. Nevertheless, you have to provide GENOMIC_SEQ and OUT_DBNAME
   # if you decide to use the E2G config these settings will be ignored !
   
   DO_NOT_READ_THE_EXONERATE2GENES_CONFIG_FILE => 1, 

       # this should point to dir || file where dumped genomic seq lives 
       
       GENOMIC_SEQ => "/data/blastdb/Ensembl/Dog/BROADD2/genome/softmasked_dusted/toplevel.fa", 

       # db where you want to write the results to ( should point to a db-connection in Databases.pm)
       # ( you can also supply a hash-ref with the db-connection parameters if you like ) 
       
       OUT_DBNAME  => "ORTHOLOGUE_DB" ,   

       # do you want to run the exonerate jobs now or do you want to setup a post-analysis
       # and upload the input ids ?  
       
       SETUP_POST_ANALYSIS___DO_NOT_RUN_EXONERATES_NOW => 0 ,  

       # the default can be found in Runnable/BaseExonerate.pm ( exonerate 0.8.3 )  
       
       EXONERATE_PROGRAM_FILE =>"/usr/local/ensembl/bin/exonerate-1.0.0", 



   # output-directory where we will write the files to 
   
   OUTPUT_DIR  => "/lustre/scratch1/ensembl/jhv/patches/dog",  

   # WARNING !!!!! Only use the names out of the genome_db table in the compara databaase 
  
 
    SPECIES_TO_COMPARE => [ 'Mus musculus' , 'Homo sapiens'] ,    

    #
    # S T A T I S T I C S
    #


    # this is for module FullStats.pm to plot exon/intron/cds/distributions 
    # of all genes in the databases no matter if there are orthologues or not 
   
       FULL_STATS_SPECIES       =>  [ "Mus musculus","Homo sapiens"],   

    # if you want to limit the statistics-analysis to certain biotypes add them here - 
    # otherwise, all genes will be used
    
       LIMIT_TO_GENE_BIOTYPES   =>  ["protein_coding" ], 

       R_BINARY_LOCATION => "/vol/software/linux-x86_64/R-2.4.0/bin/R",  

       R_OUTPUT_DIR  => "/lustre/work1/ensembl/jhv/project_genestructure_comparison/r_output",  
       #
       # dir with old scripts : 
       # R_SCRIPTS_DIR => "/lustre/work1/ensembl/jhv/project_genestructure_comparison/rscripts_old/", 
       #
       # new scripts are kept in cvs-personal 
       R_SCRIPTS_DIR =>  "/nfs/acari/jhv/cvs_checkout/ensembl-personal/jhv/projects/genestruc_comp/r_scripts",


       # specify here the r-script you like to run ( it should be in the dir R_SCRITPTS_DIR,
       # the data-file you want to use as input and, if your script writs output-files, 
       # specify the name of the output -file
       R_ANALYSIS => {
                      # plot CDS-length per transcript and exons per transcript 
                      "gt_analysis" => {
                                        'r_script'   => "exon_dist.R", 
                                        'data_file'  => "full_exon_report.txt", 
                                        'ouput_file' => "",
                                        'use_sweave' => 0,
                                        },

#                      # plot CDS-length per transcript
#                      "cds_length_trans" => {
#                                        'r_script'      => "trans_cds_length_stats.r",
#                                        'data_file'     => "full_gene_trans_report.txt", 
#                                        'use_sweave'   =>1 ,
#                                        'sweave_script' => "trans_cds_length_stats.Snw",
#                                        }, 
#
#
#                      # plot percentage of exons-per-transcript in genome  
#                      "exon_per_trans" => {
#                                        'r_script'   => "full_gene_trans_stats_relative.r",
#                                        'data_file'  => "full_gene_trans_report.txt",
#                                        'use_sweave'   =>1 ,
#                                        'sweave_script' => "full_gene_trans_stats_relative.Snw",
#                                        },
                     }, 
  );




sub import {
    my ($callpack) = caller(0); 
    my $pack = shift; 
    my @vars = @_ ? @_ : keys(%Config);
    return unless @vars;
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;
    foreach (@vars) {
	if (defined $Config{ $_ }) {
            no strict 'refs';
	    *{"${callpack}::$_"} = \$Config{ $_ };
	} else {
	    die "Error: Config: $_ not known\n";
	}
    }
}
1;
