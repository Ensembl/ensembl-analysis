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

    #print "Funcgen::new\n";
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($experiment, $result_set_analysis) = rearrange
        ( ['EXPERIMENT', 'RESULT_SET_ANALYSIS'], @args);
    #print Dumper ($experiment, $result_set_analysis);
    
    throw("Must instantiate RunnableDB with an experiment name and ".
          "result set analysis.") if (!$experiment || !$result_set_analysis);

    my $ea = $self->db->get_ExperimentAdaptor();
    my $e = $ea->fetch_by_name($experiment);
    throw("Can't fetch experiment with name ".$experiment) if (!$e);
    $self->experiment($e);

    my $aa = $self->db->get_AnalysisAdaptor();
    my $a = $aa->fetch_by_logic_name($result_set_analysis);
    throw("Can't fetch result set analysis ".$result_set_analysis) if (!$a);
    $self->result_set_analysis($a);

    my $sa = $self->db->get_SliceAdaptor();
    #warn( Dumper($self->input_id) );
    my $slice = $sa->fetch_by_name($self->input_id);
    throw("Can't fetch slice ".$self->input_id) if (!$slice);
    $self->query($slice);

    #print Dumper ($self);
    return $self;
}

### 

sub experiment {
    my $self = shift;
    $self->{'experiment'} = shift if(@_);
    return $self->{'experiment'};
}

sub result_set_analysis {
    my $self = shift;
    $self->{'result_set_analysis'} = shift if(@_);
    return $self->{'result_set_analysis'};
}

sub result_sets {
    my ($self, $rsets) = @_;
    $self->{'result_sets'} = $rsets if ($rsets);
    throw("No result_set in RunnableDB.") 
        if (!$self->{'result_sets'});
    return $self->{'result_sets'};
}

sub feature_type {
    my $self = shift;
    $self->{'feature_type'} = shift if(@_);
    return $self->{'feature_type'};
}

sub cell_type {
    my $self = shift;
    $self->{'cell_type'} = shift if(@_);
    return $self->{'cell_type'};
}

################################################################################
### Declare and set up config variables
################################################################################

sub read_and_check_config {

    #print "FunGen::read_and_check_config\n";

    my ($self, $config) = @_;

    $self->SUPER::read_and_check_config($config);

    ### check if analysis with this logic_name already exists in the database
    ### and has same program version and parameters
 
    my $aa = $self->db->get_AnalysisAdaptor();
    my $logic_name = $self->analysis->logic_name();
    my $analysis = $aa->fetch_by_logic_name($logic_name);
    if (defined $analysis) {
        print "Analysis with logic_name '$logic_name' already exists!\n";
     
        my @options = split(/\s+--?/, $self->PARAMETERS);
        my @options1 = split(/\s+--?/, $analysis->parameters);
        #print Dumper(@options, @options1);

        my %seen;
        foreach my $o (@options, @options1) {
            next if ($o =~ m/^$/ || $o eq 'NULL');
            my ($k, $v) = split /[= ]/, $o;
            $v =~ s/[\"\']//g;
            $seen{$k}{$v}++;
        }
        #print Dumper %seen;
        
      CHECK:
        foreach my $k (keys %seen) {
            foreach my $v (keys %{$seen{$k}}) {
                if ($seen{$k}{$v} < 2) {
                    throw('Analysis with logic name \''.$logic_name.'\' already '.
                          'exists, but has different options! Must choose different '.
                          'logic_name for analysis.');
                }
            }
        }
    } else {

        $self->analysis->parameters($self->PARAMETERS);
        $self->analysis->program($self->PROGRAM);
        $self->analysis->program_file($self->PROGRAM_FILE);
        $self->analysis->program_version($self->VERSION);
        $self->analysis->parameters($self->PARAMETERS);
        $self->analysis->displayable(1);
        #print Dumper $self->analysis;

    }

    ### check feature_type and cell_type info for experiment

    my $eca = $self->db->get_ExperimentalChipAdaptor();
    my $echips = $eca->fetch_all_by_Experiment($self->experiment());
    
    my (%ftype, $ftype, %ctype, $ctype);
    foreach my $echip (@{$echips}) {
        $ftype = $echip->feature_type();
        throw("No feature type defined for experimental chip (unique id: ".
              $echip->unique_id.")\n") if (! defined $ftype);
        $ftype{$ftype->dbID()}++;

        $ctype = $echip->cell_type();
        throw("No cell type defined for experimental chip (unique id: ".
              $echip->unique_id.")\n") if (! defined $ctype);
        $ctype{$ctype->dbID()}++;
    }
    
    #print Dumper %ftype;
    throw("Experiment '".$self->experiment."' has more than one feature type")
        if (keys(%ftype) > 1);
    $self->feature_type($ftype);

    #print Dumper %ctype;
    throw("Experiment '".$self->experiment."' has more than one cell type")
        if (keys(%ctype) > 1);
    $self->cell_type($ctype);

}

### container for CONFIG parameter

sub PROGRAM {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_PROGRAM'} = $value;
  }

  if ( exists( $self->{'_CONFIG_PROGRAM'} ) ) {
    return $self->{'_CONFIG_PROGRAM'};
  } else {
    return undef;
  }
}

sub PROGRAM_FILE {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_PROGRAM_FILE'} = $value;
  }

  if ( exists( $self->{'_CONFIG_PROGRAM_FILE'} ) ) {
    return $self->{'_CONFIG_PROGRAM_FILE'};
  } else {
    return undef;
  }
}

sub VERSION {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_VERSION'} = $value;
  }

  if ( exists( $self->{'_CONFIG_VERSION'} ) ) {
    return $self->{'_CONFIG_VERSION'};
  } else {
    return undef;
  }
}

sub PARAMETERS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_PARAMETERS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_PARAMETERS'} ) ) {
    return $self->{'_CONFIG_PARAMETERS'};
  } else {
    return undef;
  }
}

sub RESULT_SET_REGEXP {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_RESULT_SET_REGEXP'} = $value;
  }

  if ( exists( $self->{'_CONFIG_RESULT_SET_REGEXP'} ) ) {
    return $self->{'_CONFIG_RESULT_SET_REGEXP'};
  } else {
    return undef;
  }
}

################################################################################
### end of config
################################################################################

=head2 fetch_input

  Arg [1]     : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen
  Description : fetch data out of database and create runnable
  Returns     : 1
  Exceptions  : none
  Example     : 

=cut

sub fetch_input {
    my ($self) = @_;
    #print "Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::fetch_input\n";

    #warn("write file: ", Dumper $self->db);

    my $rsets = $self->fetch_ResultSets();
    #print Dumper $rsets;
    $self->result_sets($rsets);
    print "No. of result sets: ", scalar(@$rsets), "\n";

    my %probe_features = ();
    foreach my $rset (@{$rsets}) {
        #print Dumper $rset->name();

        print join(" ", $rset->dbID, $rset->name), "\n";
        
        my $probe_features = $rset->get_ResultFeatures_by_Slice($self->query());
        print "No. of ResultFeatures_by_Slice:\t", scalar(@$probe_features), "\n";

        throw("No result_features on slice ".$self->query()->name())
            if scalar(@$probe_features) == 0;

        my @probe_features = ();
        my $ft_cnt = 1;
        foreach my $prb_ft (@{$probe_features}) {
            #print join(" ", $self->query()->seq_region_name, 
            #           @$prb_ft, $ft_cnt++), "\n";
            push (@probe_features,
                  [ $self->query()->seq_region_name, @$prb_ft, $ft_cnt++ ]);
        }

        $probe_features{$rset->name} = \@probe_features;

    }
    
    #$self->probe_features(\%features);
    #print Dumper $self->probe_features();

    my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'.$self->analysis->module;
    $runnable = $runnable->new
        (
         -query => $self->query,
         -program => $self->analysis->program_file,
         -analysis => $self->analysis,
         -probe_features => \%probe_features,
         );
    
    $self->runnable($runnable);

    return 1;

}

#sub fetch_input_alt {
#    my ($self) = @_;
#    my $pfa = $self->db->get_ProbeFeatureAdaptor();
#    my @features = ();
#    #print "result_set dbID: ", $rset->dbID(), "\n";
#    my $prb_ft = $pfa->fetch_all_by_Slice_ExperimentalChips
#        ( $self->query(), $rset->get_ExperimentalChips() );
#    #print Dumper $prb_ft;
#    #print "Getting list of features ";
#    foreach my $ft (@{$prb_ft}) {
#        next if (! defined $ft->get_result_by_ResultSet($rset));
#        #print join(' ', $ft->dbID(), $ft->start(), $ft->end()), "\n";
#        push @features, [ $ft->probe->get_probename(), 
#                          $ft->slice()->seq_region_name(),
#                          $ft->start(),
#                          $ft->end(), 
#                          $ft->get_result_by_ResultSet($rset) ];
#        #print ".";
#    }
#    #print "done!\n";
#    
#    throw("EXIT: No probe results features for result set!") if (! @features);
#    warn("No. of features:\t", scalar(@features));
#
#    my %features = ();
#    $features{$rset->name()} = \@features;
#    $self->probe_features(\%features);
#
#    # set program options defined in Config
#    my %parameters = %{$self->parameters_hash($self->OPTIONS)};
#    #print Dumper %parameters;
#
#    my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'.$self->LOGIC_NAME;
#    $runnable = $runnable->new
#        (
#         -query => $self->query,
#         -program => $self->analysis->program_file,
#         -analysis => $self->analysis,
#         -features => $self->probe_features,
#         %parameters,
#         );
#    
#    $self->runnable($runnable);
#    return 1;
#}

=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Chipotle
  Function  : set analysis and slice on each feature
  Returntype: 1
  Exceptions: none
  Example   : 

=cut

sub write_output{

    print "RunnableDB::Funcgen::write_output\n";
    my ($self) = @_;

    my $analysis_dbID = $self->analysis()->dbID();
    if (! defined $analysis_dbID) {

        my $aa = $self->db->get_AnalysisAdaptor();
        print "Storing new analysis '".$self->analysis->logic_name()."'\n";
        $analysis_dbID = $aa->store($self->analysis());

    }
    
    ### get feature_set/data_set name by determine the longest common 
    ### substring of the result_set names
    my @names;
    map { push @names, $_->name } @{$self->result_sets()};
    #print Dumper @names;
    my %hash = (); my $i = 0;
    while (keys(%hash) != 1) {
        %hash = ();
        foreach my $n (@names) {
            $n =~ s/(_[^_]+){$i}$//;
            print $n, "\n";
            $hash{$n}++;
        }
        $i++;
    }
    my $set_name = $self->analysis->logic_name.'_'.(keys(%hash))[0];
    #print Dumper $set_name;
    
    ### implement check if this is a subset of the chosen replicates by 
    ### comparing number of keys in %hash with return value from still to be implemented method

    throw("Need to implement/update how to store fset/dset/rset via new() and add_ResultSet() ".
          "of DataSet and add contains_ResultSet ... (see comments!)")
    # Logic: 1) Check whether feature_set already exists in data_set otherwise create 2) for each 
    # ResultSet check whether it is already stored with FeatureSet, otherwise add_ResultSet to
    # FeatureSet
    
    ### feature_set
    my $fsa = $self->db->get_FeatureSetAdaptor();
    my $fset = $fsa->fetch_by_name($set_name);
    if (defined $fset) {
        
        warn("Feature set '".$set_name."' already exists.");
        #print Dumper $fset->name;

    } else {
    
        $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
            (
             -analysis => $self->analysis,
             -feature_type => $self->feature_type,
             -cell_type => $self->cell_type,
             -name => $set_name
             );
        #print Dumper $fset;
        
        $fset = $fsa->store($fset)->[0];
    }
    #print Dumper $fset;

    ### data set
    my $dsa = $self->db->get_DataSetAdaptor();
    my $dsets = $dsa->fetch_all_by_FeatureSet($fset);
    #print Dumper $dsets;

    my %rset_ids;
    foreach my $ds (@{$dsets}) {
        foreach my $rs (@{$ds->get_ResultSets}) {
            $rset_ids{$rs->dbID}++;
        }
    }
    
    foreach my $rset (@{$self->result_sets()}) {
        next if (exists $rset_ids{$rset->dbID});
        
        # add new data set
        my $dset = Bio::EnsEMBL::Funcgen::DataSet->new
            (
             -result_set => $rset,
             -feature_set => $fset,
             -name => $set_name
             );
        $dsa->store($dset);
        
    }

    ### predicted features
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
    
    my $pfa = $self->db->get_PredictedFeatureAdaptor();
    my $pf = $pfa->fetch_all_by_Slice_FeatureSet($self->query, $fset);
    if (@$pf) {
        
        throw("NOT IMPORTING ".scalar(@{$self->output})." predicted features! Slice ".
              join(':', $self->query->seq_region_name,$self->query->start,$self->query->end).
              " already contains ".scalar(@$pf)." predicted features of feature set ".
              $fset->dbID.".");
        
    } else {
        my @pf;
        foreach my $ft (@{$self->output}){
            
            #print Dumper $ft;
            my ($seqid, $start, $end, $score) = @{$ft};
            
            my $pf = Bio::EnsEMBL::Funcgen::PredictedFeature->new
                (
                 -slice         => $self->query,
                 -start         => $start,
                 -end           => $end,
                 -strand        => 0,
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

1;
