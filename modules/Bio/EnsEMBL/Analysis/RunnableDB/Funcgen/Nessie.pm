# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Nessie
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie

=head1 SYNOPSIS

  my $nessie = Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie->new
     (
        -input_id => 'chromosome::20:1:100000:1',
        
     );
  $nessie->fetch_input;
  $nessie->run;
  $nessie->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable Nessie which wraps the program Nessie

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the predicted_feature table with respect to the resut_set, data_set 
and feature_set table

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Funcgen::Nessie;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA); 

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Funcgen);



=head2 new

  Arg [1]     : 
  Arg [2]     : 
  Description : Instantiates new Nessie runnabledb
  Returntype  : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie object
  Exceptions  : 
  Example     : 

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->read_and_check_config($NESSIE_CONFIG);
  #print Dumper $self;

  return $self;
}

=head2 fetch_input

  Arg [1]     : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie
  Description : fetch data out of database and create runnable
  Returns     : 1
  Exceptions  : none
  Example     : 

=cut

sub fetch_input {
    my ($self) = @_;
    print "Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie::fetch_input\n";

    my $sa = $self->db->get_SliceAdaptor();
    #warn( Dumper($self->input_id) );
    my $slice = $sa->fetch_by_name($self->input_id);
    #warn( Dumper($sa, $slice) );
    $self->query($slice);

    my $ea = $self->db->get_ExperimentAdaptor();
    my $e = $ea->fetch_by_name($self->experiment);
    
    my $analysis = $self->EFG_EXPERIMENT->{$self->experiment}->{ANALYSIS};
    my $aa = $self->db->get_AnalysisAdaptor();
    my $a = $aa->fetch_by_logic_name($analysis);

    my $rsa = $self->db->get_ResultSetAdaptor();
    my $rsets = $rsa->fetch_all_by_Experiment_Analysis($e, $a);
    #print scalar(@$rsets);
    $self->result_sets($rsets);
    
    my %features = ();
    foreach my $rset (@{$rsets}) {
        #print Dumper $rset->name();

        next if ($rset->name() !~ m/_TR\d+$/);
        print join(" ", $rset->dbID, $rset->name), "\n";

        my $results = $rset->get_ResultFeatures_by_Slice($self->query());
        #print "  ResultFeatures_by_Slice:\t", scalar(@$results), "\n";

        my @features = ();
        my $ft_cnt = 1;
        foreach my $prb_ft (@{$results}) {
            #print join(" ", $self->query()->seq_region_name, 
            #           @$prb_ft, $ft_cnt++), "\n";
            push (@features,
                  [ 'chr'.$self->query()->seq_region_name, @$prb_ft, $ft_cnt++ ]);
        }

        $features{$rset->name} = \@features;

        #print Dumper $results;

    }
    
    $self->probe_features(\%features);
    #print Dumper $self->probe_features();

    # set program options defined in Config
    my %parameters = %{$self->parameters_hash($self->OPTIONS)};
    #print Dumper %parameters;

    if(!$self->analysis->program){
        $self->analysis->program($self->PROGRAM);
    }
    if(!$self->analysis->program_file){
        $self->analysis->program_file($self->PROGRAM);
    }
    
    my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'.$self->LOGIC_NAME;
    $runnable = $runnable->new
        (
         -query => $self->query,
         -program => $self->analysis->program_file,
         -analysis => $self->analysis,
         -features => $self->probe_features,
         %parameters,
         );
    
    $self->runnable($runnable);

    return 1;

}

################################################################################
# Declare and set up config variables
################################################################################

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($NESSIE_CONFIG);

  ##########
  # CHECKS
  ##########

  # check that compulsory options have values
  foreach my $config_var 
      (
       qw(
          PROGRAM
          OPTIONS
          LOGIC_NAME
          EFG_EXPERIMENT
          )
       ){
          if ( not defined $self->$config_var ){
              throw("You must define $config_var in config.");
          }
          
          #print Dumper $self->$config_var;
          
      }
  
  # make sure EFG_EXPERIMENT exists, is a hash, and contains all 
  # compulsory options
  throw("EFG_EXPERIMENT for ".$self->experiment." is not defined.")
      if(! exists $self->EFG_EXPERIMENT->{$self->experiment});
  throw("EFG_EXPERIMENT must be a hash ref not ".$self->EFG_EXPERIMENT.
        " Chipotle::read_and_check_config")
      if(ref($self->EFG_EXPERIMENT) ne 'HASH');
  foreach my $config_var
      (
       qw(
          FT_NAME
          FT_CLASS
          FT_DESC
          CT_NAME
          CT_DESC
          )
       ){
          throw("Must define $config_var in EFG_EXPERIMENT config.")
          if (! defined $self->EFG_EXPERIMENT->{$self->experiment}->{$config_var});
      }

}

#############################################################
###     end of config
#############################################################

1;
