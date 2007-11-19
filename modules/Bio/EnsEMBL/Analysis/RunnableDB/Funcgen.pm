# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Funcgen
#
# Copyright (c) 2007 Ensembl
#
=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Fungen

=head1 SYNOPSIS

=head1 DESCRIPTION

This module is the base class for the Fungen Runnabledbs that act as an 
interface between the functional genomics database and the Funcgen Runnables 
both fetching input data and writing data back to the databases.

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

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
    print "Analysis::RunnableDB::Funcgen::new\n";
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    #warn( Dumper($self->input_id) );
    my $slice = $self->adaptor('Slice')->fetch_by_name($self->input_id);
    #print Dumper $slice;
    throw("Can't fetch slice ".$self->input_id) if (!$slice);
    $self->query($slice);
    
    #print Dumper ($self);
    return $self;
}

### 

# generic ensembl db adapter method
sub adaptor {
    my ($self, $name) = @_;
    my $adaptor = $name.'Adaptor';
    my $method = 'get_'.$adaptor;

    $self->{$adaptor} = $self->db->$method
        if (! $self->{$adaptor} );
    throw("No ".$adaptor."Adaptor in RunnableDB.") 
        if (!$self->{$adaptor});
    return $self->{$adaptor};
}

sub experiment {
    my $self = shift;
    $self->{'experiment'} = shift if(@_);
    return $self->{'experiment'};
}

sub norm_analysis {
    my $self = shift;
    $self->{'norm_analysis'} = shift if(@_);
    return $self->{'norm_analysis'};
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

sub DataSet {
    my ($self, $set) = @_;
    $self->{'DataSet'} = $set if ($set);
    throw("No DataSet in RunnableDB.") 
        if (!$self->{'DataSet'});
    return $self->{'DataSet'};
}

sub FeatureSet {
    my ($self, $set) = @_;
    $self->{'FeatureSet'} = $set if ($set);
    throw("No FeatureSet in RunnableDB.") 
        if (!$self->{'FeatureSet'});
    return $self->{'FeatureSet'};
}

sub ResultSets {
    my ($self, $sets) = @_;
    $self->{'ResultSets'} = $sets if ($sets);
    throw("No ResultSets in RunnableDB.") 
        if (!$self->{'ResultSets'});
    return $self->{'ResultSets'};
}

################################################################################
### Declare and set up config variables
################################################################################

sub read_and_check_config {

    print "Analysis::RunnableDB::Funcgen::read_and_check_config\n";

    my ($self, $config) = @_;

    #print Dumper $config;

    $self->SUPER::read_and_check_config($config);

    # make sure we have a valid analysis object
    $self->check_Analysis();
    
    # make sure we have valid feature and data set objects incl. supporting sets
    $self->check_Sets();
    
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

sub EXPERIMENT {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_EXPERIMENT'} = $value;
    }

    if ( exists( $self->{'_CONFIG_EXPERIMENT'} ) ) {
        return $self->{'_CONFIG_EXPERIMENT'};
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

sub NORM_ANALYSIS {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_NORM_ANALYSIS'} = $value;
    }

    if ( exists( $self->{'_CONFIG_NORM_ANALYSIS'} ) ) {
        return $self->{'_CONFIG_NORM_ANALYSIS'};
    } else {
        return undef;
    }
}

sub DATASET_NAME {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_DATASET_NAME'} = $value;
    }

    if ( exists( $self->{'_CONFIG_DATASET_NAME'} ) ) {
        return $self->{'_CONFIG_DATASET_NAME'};
    } else {
        return undef;
    }
}

sub ANALYSIS_WORK_DIR {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_ANALYSIS_WORK_DIR'} = $value;
    }

    if ( exists( $self->{'_CONFIG_ANALYSIS_WORK_DIR'} ) ) {
        return $self->{'_CONFIG_ANALYSIS_WORK_DIR'};
    } else {
        return undef;
    }
}

sub check_Analysis {

    my ($self) = @_;

    # get/set experiment
    my $e = $self->adaptor('Experiment')->fetch_by_name($self->EXPERIMENT);
    throw("Can't fetch experiment with name ".$self->EXPERIMENT) if (!$e);
    $self->experiment($e);

    # get/set norm_analysis
    #print Dumper $ENV{NORM_ANALYSIS};
    my $a = $self->adaptor('Analysis')->fetch_by_logic_name($self->NORM_ANALYSIS);
    throw("Can't fetch result set analysis ".$self->NORM_ANALYSIS) if (!$a);
    $self->norm_analysis($a);

    ### check if analysis with this logic_name already exists in the database
    ### and has same program version and parameters

    my $logic_name = $self->analysis->logic_name();

    my $analysis = Bio::EnsEMBL::Analysis->new
        (
         -logic_name => $logic_name,
         -module => $self->analysis->module,
         -program => $self->PROGRAM,
         -program_file => $self->PROGRAM_FILE,
         -program_version => $self->VERSION,
         -parameters => $self->PARAMETERS,
         -displayable => 1
         );

    #print Dumper ($self->PROGRAM, $self->analysis->program);
    #print Dumper ($self->PROGRAM_FILE, $self->analysis->program_file);
    #print Dumper ($self->VERSION, $self->analysis->program_version);
    #print Dumper ($self->PARAMETERS, $self->analysis->parameters);
    
    if (defined $self->analysis->dbID) {
        
        ### analysis compare
        # returns  1 if this analysis is special case of given analysis
        # returns  0 if they are equal
        # returns -1 if they are completely different
        throw('Analysis with logic name \''.$logic_name.'\' already '.
              'exists, but has different options! Must choose different '.
              'logic_name for analysis and edit config hash of $self accordingly.') 
            if ($self->analysis->compare($analysis));

    } else {

        warn("New Analysis with logic name $logic_name.");
        #print Dumper $analysis;
        my $dbID = $self->adaptor('Analysis')->store($analysis);
        $self->analysis($self->adaptor('Analysis')->fetch_by_dbID($dbID));

    }

}

sub check_Sets {

    my ($self) = @_;

    ### get result sets to process and check that result_set, data_set, 
    ### feature_set are storable

    my $rsets = $self->fetch_ResultSets();

    ### get feature_set/data_set name by determining the longest common 
    ### prefix of the result_set names
    #print Dumper @{$self->ResultSets()};
    my @names = map { $_->name } @{$self->ResultSets()};
    #print Dumper @names;

    my %hash = ();
    foreach my $n (@names) {
        my @chars = split(//, $n);

        my $string = '';
        foreach my $c (@chars) {
            $string .= $c;
            #print Dumper $string;
            $hash{$string}++;
        }
    }
    my $lcp = '';

    foreach (sort keys %hash) {
        last if ($hash{$_} < scalar(@$rsets));
        #print Dumper( $_, $hash{$_}, scalar(@$rsets) ) ;
        $lcp = $_;
    } 
    $lcp =~ s/_$//;

    my $set_name = $lcp.'_'.$self->DATASET_NAME;
    print 'Set name: ', $set_name, "\n";
    #print Dumper $self->analysis;

    my $fset = $self->adaptor('FeatureSet')->fetch_by_name($set_name);
    #print Dumper $fset;

    if ( ! defined $fset ) {

        $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
            (
             -analysis => $self->analysis,
             -feature_type => $self->feature_type,
             -cell_type => $self->cell_type,
             -name => $set_name,
             -type => 'annotated'
             );
        #print Dumper $fset;
        
        warn("Storing new feature set \'$set_name\'");
        eval { $fset = ${$self->adaptor('FeatureSet')->store($fset)}[0] };
        throw("Coudn't store feature set \'$set_name\': $!") if ($@);

    } else {

        warn("Feature set with name $set_name already exists.");

    }
    # save FeatureSet
    $self->FeatureSet($fset);

    my $dset = $self->adaptor('DataSet')->fetch_by_name($set_name);
    #print Dumper $dset;

    if ( ! defined $dset ) {

        $dset = Bio::EnsEMBL::Funcgen::DataSet->new
            (
             -SUPPORTING_SETS     => $self->ResultSets,
             -FEATURE_SET         => $fset,
             -DISPLAYABLE         => 1,
             -NAME                => $set_name,
             -SUPPORTING_SET_TYPE => 'result',
             );
        #print Dumper $dset;

        warn("Storing new feature set \'$set_name\'");
        eval { $dset = ${$self->adaptor('DataSet')->store($dset)}[0] };
        throw("Coudn't store data set \'$set_name\': $!") if ($@);

    } else {

        warn("Data set with name $set_name already exists.");

        # need to check whether ResultSets and supporting_sets are the same and 
        # possibly add ResultSet to supporting_sets

        my $ssets = $dset->get_supporting_sets();

        my %ssets_dbIDs = ();
        map { $ssets_dbIDs{$_->dbID}='' } (@{$ssets});
        map { 
            $dset->add_supporting_sets([ $_ ]) if (! exists $ssets_dbIDs{$_->dbID}); 
        } @{$self->ResultSets} 

    }
    # save DataSet
    $self->DataSet($dset);

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
    print "Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::fetch_input\n";
    
    #warn("write file: ", Dumper $self->db);
    
    my %result_features = ();
    foreach my $rset (@{$self->ResultSets}) {
        #print Dumper $rset->name();
        
        print join(" ", $rset->dbID, $rset->name), "\n";
        
            my $result_features = $rset->get_ResultFeatures_by_Slice($self->query());
            print "No. of ResultFeatures_by_Slice:\t", scalar(@$result_features), "\n";
            
            throw("No result_features on slice ".$self->query()->name())
                if scalar(@$result_features) == 0;
            
            my @result_features = ();
            my $ft_cnt = 1;
            foreach my $prb_ft (@{$result_features}) {
                #print join(" ", $self->query()->seq_region_name, 
                #           @$prb_ft, $ft_cnt++), "\n";
                push (@result_features,
                      [ $self->query()->seq_region_name, @$prb_ft, $ft_cnt++ ]);
            }

            $result_features{$rset->name} = \@result_features;

        }
        
        #print Dumper %result_features;
        
        #$self->result_features(\%features);
        #print Dumper $self->result_features();

        my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'
            .$self->analysis->module;
        $runnable = $runnable->new
            (
             -query => $self->query,
             -program => $self->analysis->program_file,
             -analysis => $self->analysis,
             -result_features => \%result_features,
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

        print "RunnableDB::Funcgen::write_output\n";
        my ($self) = @_;

        # store analysis, feature set and data set
        # analysis was alredy been stored while checking config in read_and_check_config
        #$self->adaptor('Analysis')->store($self->analysis());
        if (! defined $self->FeatureSet->dbID) {
            $self->adaptor('FeatureSet')->store($self->FeatureSet());
        }
        if (! defined $self->DataSet->dbID) {
            $self->adaptor('DataSet')->store($self->DataSet());
        }

        ### annotated features

        my ($transfer, $slice);
        if($self->query->start != 1 || $self->query->strand != 1) {
            my $sa = $self->adaptor('Slice');
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
        
#    my $af = $self->adaptor('AnnotatedFeature')->fetch_all_by_Slice_FeatureSet(
#		$self->query, $fset);

        my $fset = $self->FeatureSet;	
        my $fs_id = $fset->dbID();
        my $constraint = qq( af.feature_set_id = $fs_id );

        my $af = $self->adaptor('AnnotatedFeature')->fetch_all_by_Slice_FeatureSet
            ($slice, $fset, $self->analysis->logic_name);
        
        print 'No. of annotated features already stored: '.scalar(@$af)."\n";
        print 'No. of annotated features to be stored: '.scalar(@{$self->output})."\n";

        if (@$af) {
            
            throw("NOT IMPORTING ".scalar(@{$self->output})." annotated features! Slice ".
                  join(':', $self->query->seq_region_name,$self->query->start,$self->query->end).
                  " already contains ".scalar(@$af)." annotated features of feature set ".
                  $fset->dbID.".");
            
        } else {
            my @af;
            foreach my $ft (@{$self->output}){
                
                #print Dumper $ft;
                my ($seqid, $start, $end, $score) = @{$ft};
                
                #print Dumper ($seqid, $start, $end, $score);
                my $af = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
                    (
                     -slice         => $self->query,
                     -start         => $start,
                     -end           => $end,
                     -strand        => 0,
                     -display_label => $self->analysis->logic_name,
                     -score         => $score,
                     -feature_set   => $fset,
                     );
                
                # make sure feature coords are relative to start of entire seq_region
                if ($transfer) {
                    #warn("original af:\t", join("\t", $af->start, $af->end), "\n");
                    $af = $af->transfer($slice);
                    #warn("remapped af:\t", join("\t", $af->start, $af->end), "\n");
                }
                push(@af, $af);
            }

            $self->adaptor('AnnotatedFeature')->store(@af);
        }
        return 1;
    }


=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : fetch and set ResultSets of interest
  Returntype: 1
  Exceptions: none
  Example   : 

=cut


sub fetch_ResultSets
{
    print "Analysis::RunnableDB::Funcgen::fetch_ResultSets\n";

    my $self = shift;
    #print Dumper $self;

    my $rsa = $self->adaptor('ResultSet');

    ### Can we speed this up???
    my $rsets = $rsa->fetch_all_by_Experiment_Analysis
        ($self->experiment, $self->norm_analysis);

    print "No. of available ResultSets: ", scalar(@$rsets), "\n";

    my @rsets = ();

    #print Dumper( $self->experiment->name, $self->RESULT_SET_REGEXP);
    my $regex = $self->RESULT_SET_REGEXP;
    foreach my $rset (@{$rsets}) {
        #print Dumper $rset->name();
        next if ($rset->name() !~ m/$regex/);
        push(@rsets, $rset);

        # check feature_type
        if (! defined $self->feature_type ) {
            $self->feature_type($rset->feature_type);
        } else {
            throw("replicates differ in feature types")
                if ($self->feature_type->dbID != $rset->feature_type->dbID);
        }

        # check cell_type
        if ( ! defined $self->cell_type() ) {
            $self->cell_type($rset->cell_type);
        } else {
            throw("replicates differ in cell types")
                if ($self->cell_type->dbID != $rset->cell_type->dbID);
        }

    }

    if (!@rsets) {
        
        my $rset_list = join(' ', map { $_->name } @{$rsets});
        throw ("RESULT_SET_REGEXP doesn't match any the following result set:\n$rset_list");
    }
    
    $self->ResultSets(\@rsets);

    print "Selected result sets: ", join(', ', map { $_->name } @rsets), 
    ' (in total ', scalar(@rsets), ")\n";

    return \@rsets;

}

1;
