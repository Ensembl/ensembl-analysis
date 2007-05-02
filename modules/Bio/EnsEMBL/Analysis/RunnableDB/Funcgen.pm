# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Funcgen
#
# Copyright (c) 2007 Ensembl
#
=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Fungen

=head1 SYNOPSIS

=head1 DESCRIPTION

This module is the base class for the Fungen Runnabledbs that act as an 
interface between the function genomics database and the Funcgen Runnables 
both fetching input data and writing data back to the databases.

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;

use Bio::EnsEMBL::Analysis::Config::General;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

=head2 fetch_input

  Arg [1]     : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Chipotle
  Description : fetch data out of database and create runnable
  Returns     : 1
  Exceptions  : none
  Example     : 

=cut

sub fetch_input {
    my ($self) = @_;
    
    my $sa = $self->db->get_SliceAdaptor();
    #warn( Dumper($self->input_id) );
    my $slice = $sa->fetch_by_name($self->input_id);
    #warn( Dumper($sa, $slice) );

    $self->query($slice);
    #warn("write file: ", Dumper $self->db);

    my $ea = $self->db->get_ExperimentAdaptor();
    my $e = $ea->fetch_by_name($self->EFG_EXPERIMENT);

    my $aa = $self->db->get_AnalysisAdaptor();
    my $a = $aa->fetch_by_logic_name($self->EFG_ANALYSIS);

    my $rsa = $self->db->get_ResultSetAdaptor();
    my $rsets = $rsa->fetch_all_by_Experiment_Analysis($e, $a);
    
    # maybe good enough for now, but might cause trouble in the case of 
    # result set duplication 
    throw("More than one resultset for Experiment_Analysis!") 
        if (scalar(@{$rsets}) > 1);

    my $rset = $rsets->[0];
    #print Dumper $rset;

    $self->result_set($rset);

    my $pfa = $self->db->get_ProbeFeatureAdaptor();

    my @features = ();

    #print "result_set dbID: ", $rset->dbID(), "\n";
        
    my $prb_ft = $pfa->fetch_all_by_Slice_ExperimentalChips
        ( $slice, $rset->get_ExperimentalChips() );
        
    foreach my $ft (@{$prb_ft}) {
        next if (! defined $ft->get_result_by_ResultSet($rset));
        #print join(' ', $ft->dbID(), $ft->start(), $ft->end()), "\n";
        push @features, [ $ft->probe->get_probename(), 
                          $ft->slice()->seq_region_name(),
                          $ft->start(),
                          $ft->end(), 
                          $ft->get_result_by_ResultSet($rset) ];
        
    }
    $self->probe_features(\@features);

    # set program options defined in Config
    my %parameters = %{$self->parameters_hash($self->OPTIONS)};
    #print Dumper %parameters;

    if(!$self->analysis->program_file){
        $self->analysis->program_file($self->PROGRAM);
    }

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle->new
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

=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Chipotle
  Function  : set analysis and slice on each feature
  Returntype: 1
  Exceptions: none
  Example   : 

=cut

sub write_output{

  my ($self) = @_;

  # needed adaptors
  my $pfa = $self->db->get_PredictedFeatureAdaptor();
  my $dsa = $self->db->get_DataSetAdaptor();
  my $fsa = $self->db->get_FeatureSetAdaptor();
  my $cta = $self->db->get_CellTypeAdaptor();
  my $fta = $self->db->get_FeatureTypeAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();

  # stored results_set from fetch_input
  my $rset = $self->result_set();
  
  #print Dumper $self->analysis();
  
  my $logic_name = $self->analysis->logic_name;
  if (! $logic_name) {
      $logic_name = $self->LOGIC_NAME;
  }
  #print Dumper($logic_name);  

  my $analysis = $aa->fetch_by_logic_name($logic_name);
  if (! $analysis) {

      # add analysis
      $analysis = Bio::EnsEMBL::Analysis->new
          (
           -logic_name      => $logic_name,
           #-db              => 'NULL',
           #-db_version      => 'NULL',
           #-db_file         => 'NULL',
           #-program         => 'NULL',
           #-program_version => 'NULL',
           -program_file    => $self->analysis->program_file,
           #-gff_source      => 'NULL',
           #-gff_feature     => 'NULL',
           -module          => $logic_name,
           #-module_version  => 'NULL',
           -parameters      => $self->OPTIONS,
           #-created         => 'NULL',
           -description     => $self->analysis->description,
           -display_label   => $self->analysis->display_label,
           -displayable     => 1,
           );
      $aa->store($analysis);
  }
  #print Dumper($analysis);
  
  my $ftype = $fta->fetch_by_name($self->EFG_FT_NAME);
  if (! $ftype) {
      
      # add new feature type
      $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new
          (
           -name => $self->EFG_FT_NAME,
           -class => $self->EFG_FT_CLASS,
           #-description => 'NULL'
           );
      $fta->store($ftype);
  }
  #print Dumper($ftype);
  
  my $ctype = $cta->fetch_by_name($self->EFG_CT_NAME);
  if (! $ctype) {
      
      # add new cell type
      $ctype = Bio::EnsEMBL::Funcgen::CellType->new
          (
           -name => $self->EFG_CT_NAME,
           -display_label => $self->EFG_CT_NAME,
           #-description => 'NULL'
           );
      $cta->store($ctype);
  }
  #print Dumper($ctype);
    
  # add new feature set
  my $fset = Bio::EnsEMBL::Funcgen::FeatureSet
      ->new(
            -analysis => $analysis,
            -feature_type => $ftype,
            -cell_type => $ctype
            );
  #print Dumper($fset);
  $fsa->store($fset);

  # add new data set
  my $dset = Bio::EnsEMBL::Funcgen::DataSet->new
      (
       -RESULT_SET => $self->result_set(),
       -FEATURE_SET => $fset,
       );
  #print Dumper($dset);
  $dsa->store($dset);

  my @pf;
  foreach my $ft (@{$self->output}){

      print Dumper $ft;
      my ($peakid, $seqid, $start, $end, 
          $score, # 
          $pval, # P value for that peak
          $avgwinscore # 
          ) = @{$ft};

      my $pf = Bio::EnsEMBL::Funcgen::PredictedFeature->new
          (
           -slice         => $self->query,
           -start         => $start,
           -end           => $end,
           -strand        => 1,
           -display_label => 'enriched_site',
           -score         => $score,
           -feature_set   => $fset,
           );
  
      push @pf, $pf;

  }

  $pfa->store(@pf);
  return 1;

}

=head2 probe_features

  Arg [1]     : Bio::EnsEMBL::Analysis::RunnableDB::Chipotle
  Arg [2]     : arrayref of probe features
  Description : container for probe features
  Returntype  : arrayref
  Exceptions  : throws if no probe feature container is defined
  Example     : 

=cut

sub probe_features {

  my ($self, $features) = @_;

  if($features){
      $self->{'probe_features'} = $features;
  }

  throw("No probe features available in Runnable.") 
      if (!$self->{'probe_features'});

  return $self->{'probe_features'};

}

=head2 result_set

  Arg [1]     : Bio::EnsEMBL::Analysis::RunnableDB::Chipotle
  Arg [2]     : arrayref of result_set object
  Description : container for result_set
  Returntype  : arrayref
  Exceptions  : throws if no result_set container is defined
  Example     : 

=cut

sub result_set {

  my ($self, $rset) = @_;

  if($rset){
      $self->{'result_set'} = $rset;
  }

  throw("No result_set in RunnableDB.") if (!$self->{'result_set'});

  return $self->{'result_set'};

}

1;
