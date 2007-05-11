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

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Arg [3]   : Bio::EnsEMBL::Analysis
  Function  : create a Bio::EnsEMBL::Analysis::RunnableDB object
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB
  Exceptions: throws if not passed either a dbadaptor, input id or
  an analysis object
  Example   : $rdb = $perl_path->new( -analysis => $self->analysis,
                                      -input_id => $self->input_id,
                                      -db => $self->adaptor->db );

=cut



sub new{
  my ($class,@args) = @_;
  my $self = bless {},$class;
  my ($db, $input_id, $experiment, $analysis) = rearrange
    (['DB', 'INPUT_ID', 'EXPERIMENT', 'ANALYSIS'], @args);
  if(!$db || !$analysis || !$experiment || !$input_id){
    throw("Must instantiate RunnableDB with a dbadaptor".
          ", an experiment, an analysis object".
          ", and an input_id.");
  }

  $self->db($db);
  $self->input_id($input_id);
  $self->experiment($experiment);
  $self->analysis($analysis);

  return $self;
}

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
    my $e = $ea->fetch_by_name($self->experiment);
    
    my $analysis = $self->EFG_EXPERIMENT->{$self->experiment}->{ANALYSIS};
    my $aa = $self->db->get_AnalysisAdaptor();
    my $a = $aa->fetch_by_logic_name($analysis);

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

    #print "Getting list of features ";
    foreach my $ft (@{$prb_ft}) {
        next if (! defined $ft->get_result_by_ResultSet($rset));
        #print join(' ', $ft->dbID(), $ft->start(), $ft->end()), "\n";
        push @features, [ $ft->probe->get_probename(), 
                          $ft->slice()->seq_region_name(),
                          $ft->start(),
                          $ft->end(), 
                          $ft->get_result_by_ResultSet($rset) ];
        #print ".";
    }
    #print "done!\n";
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
  ### must implement check if analysis with logic_name also has same options
  if (! defined $analysis) {

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
  
  # get userdefined experiment meta data
  my $e_conf = $self->EFG_EXPERIMENT->{$self->experiment};

  my $ftype = $fta->fetch_by_name($e_conf->{FT_NAME});
  if (! defined $ftype) {
      
      # add new feature type
      $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new
          (
           -name => $self->EFG_FT_NAME,
           -class => $self->EFG_FT_CLASS,
           -description => $self->EFG_FT_DESC
           );
      $ftype = $fta->store($ftype);
  }
  #print Dumper($ftype);
  
  my $ctype = $cta->fetch_by_name($e_conf->{CT_NAME});
  if (! defined $ctype) {
      
      # add new cell type
      $ctype = Bio::EnsEMBL::Funcgen::CellType->new
          (
           -name => $self->EFG_CT_NAME,
           -display_label => $self->EFG_CT_NAME,
           -description => $self->EFG_CT_DESC
           );
      $ctype = $cta->store($ctype);
  }
  #print Dumper($ctype);

  my $fset_name = $self->LOGIC_NAME.'_'.$e_conf->{FT_NAME}.
      '_'.$e_conf->{CT_NAME};
  my $fset = $fsa->fetch_all_by_name($fset_name);
  if (! @{$fset}) {
      # add new feature set
      $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
          (
           -analysis => $analysis,
           -feature_type => $ftype,
           -cell_type => $ctype,
           -name => $fset_name
           );
      $fset = $fsa->store($fset);
  } else {
      throw("More than one feature_set is currently not supported.")
          if (scalar(@{$fset}) > 1);
  }
  $fset = $fset->[0];
  
  my $dset = $dsa->fetch_all_by_FeatureSet($fset);
  
  my $add_dset = 1;
  foreach my $ds (@{$dset}) {
      foreach my $rs (@{$ds->get_ResultSets}){
          print Dumper $rs->dbID;
          
          if ($rs->dbID == $self->result_set()->dbID) {
              warn("Data_set ".$ds->dbID.
                   " linking result_set ".$self->result_set()->dbID.
                   " with feature_set ".$fset->dbID." already exists");
              $add_dset=0;
              last;
          }
      }
  }
  if ($add_dset) {
      # add new data set
      $dset = Bio::EnsEMBL::Funcgen::DataSet->new
          (
           -result_set => $self->result_set(),
           -feature_set => $fset,
           #-name => $fset_name
           );
      $dsa->store($dset);
  }

  my ($transfer, $slice);
  if($self->query->start != 1 || $self->query->strand != 1) {
      my $sa = $self->db->get_SliceAdaptor();
      $slice = $sa->fetch_by_region($self->query->coord_system->name(),
                                    $self->query->seq_region_name(),
                                    undef, #start
                                    undef, #end
                                    undef, #strand
                                    $self->query->coord_system->version());
      $transfer = 1;
  } else {
      $slice = $self->query;
  }

  my $pf = $pfa->fetch_all_by_Slice_FeatureSet($slice, $fset);


  if (@$pf) {

      throw("database already contains ".scalar(@$pf)." predicted features ".
            " of this feature_set on this slice. Not importing!");

  } else {
      my @pf;
      foreach my $ft (@{$self->output}){
          
          #print Dumper $ft;
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
          # make sure feature coords are relative to start of entire seq_region
          if ($transfer) {
              #warn("original pf:\t", join("\t", $pf->start, $pf->end), "\n");
              $pf = $pf->transfer($slice);
              #warn("remapped pf:\t", join("\t", $pf->start, $pf->end), "\n");
          }

          push(@pf, $pf);

      }
      $pfa->store(@pf);

  }

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

sub experiment{
  my $self = shift;
  $self->{'experiment'} = shift if(@_);
  return $self->{'experiment'};
}

1;
