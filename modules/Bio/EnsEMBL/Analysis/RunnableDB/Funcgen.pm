# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Fungen

=head1 SYNOPSIS

=head1 DESCRIPTION

This module is the base class for the Fungen Runnabledbs that act as an 
interface between the functional genomics database and the Funcgen Runnables 
both fetching input data and writing data back to the databases.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw( parse_config );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning stack_trace_dump );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


# To do
# 1 Done Remove all the ucfirst input_id_types, lc all! This has also been changed in the env 

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
    #print stack_trace_dump();
    my $self = $class->SUPER::new(@args);
    return $self;

}


################################################################################
### Declare and set up config variables
################################################################################

sub read_and_check_config {

    print "Analysis::RunnableDB::Funcgen::read_and_check_config\n";

    my ($self, $config) = @_;

    # get config for logic_name
    warn("MODULE:  \t".$self->analysis->module);
    warn("LOGIC_NAME:\t".$self->analysis->logic_name);
    warn("INPUT_ID:\t".$self->input_id);
    warn("INPUT_ID_TYPE:\t".$self->analysis->input_id_type);

    warn("Reading config for '".$self->analysis->logic_name."'");
    parse_config($self, $config, $self->analysis->logic_name);
    #print Dumper $self;

    # first get default and eFG analysis config
    #warn("Reading config for '$logic_name'");
    #parse_config($self, $config, $logic_name);
    # then get config for pipeline analysis (individual experiment settings,
    # like RESULT_SET_REGEXP)
    

	#NJ Must define dnadb first and pass to efgdb during creation
	#This is to avoid errors when autoguessing the dnadb from ensembldb
	if($self->DNADB->{-dbname}){
	  $self->dnadb(Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{ $self->DNADB })); 
	}


    # Make sure we have the correct DB adaptors!!!
    $self->efgdb(Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
				 (
				  %{ $self->EFGDB },
				  -dnadb => $self->dnadb,
				 ));

   
    # Set analysis
    my $efg_analysis = new Bio::EnsEMBL::Analysis( -logic_name => $self->analysis->logic_name );
    $self->efg_analysis($efg_analysis);

    # extract experiment name from input_id

    # Set experiment and query
    $self->check_InputId();

    # Make sure we deal with the correct analysis object!!!
    # -- moved to child's news method --
    #$self->analysis->logic_name($logic_name);

    # Set normalization method (only for ChIP-chip on slices)
    unless ($self->analysis->input_id_type eq "file") {
        my $m = $self->efgdb->get_AnalysisAdaptor->fetch_by_logic_name($self->NORM_METHOD);
        $self->norm_method($m) 
            or throw("Can't fetch analysis object for norm method ".$self->NORM_METHOD);
        #print Dumper $self->norm_method;
    }

    # make sure we have valid feature and data set objects incl. supporting sets
    # -- moved to child's news method --
    #$self->check_Sets();
    

    # set score factor (hack to detect depletion)
    #if ($self->experiment->name =~ m/H3-Core|H3K9me2/) {
    #    $self->SCORE_FACTOR(-1);
    #} else {
        $self->SCORE_FACTOR(1);
    #}

}

sub check_InputId {

    print "Analysis::RunnableDB::Funcgen::check_InputId\n";

    my ($self) = @_;
    my @input = split(':', $self->input_id);

    # set experiment
    my $ename = shift @input;
    my $e = $self->efgdb->get_ExperimentAdaptor->fetch_by_name($ename);
    $self->experiment($e) 
        or throw("Can't fetch experiment with name ".$ename);
    #print Dumper $self->experiment;
    warn("EXPERIMENT:\t".$self->experiment->name);

    # set query
    my $input_id_type = $self->analysis->input_id_type();

    my $sa = $self->efgdb->get_SliceAdaptor;

    my @slices = ();
    if ($input_id_type eq 'slice') {

        @slices = ( $sa->fetch_by_name(join(':', @input)) );
        throw("Can't fetch slice with name ".join(':', @input)) 
            unless (@slices);
        
    } elsif ($input_id_type eq 'array') {

        @slices = @{$sa->fetch_all('toplevel')};
        throw("Can't fetch toplevel slices")
            unless (@slices);
        
    } else {

        throw("Input_id type '$input_id_type' not implemented.");

    }

    $self->query(\@slices);

}

sub check_Analysis {

    my ($self) = @_;

    # add/update efg_analysis config
    $self->efg_analysis->module($self->MODULE);
    $self->efg_analysis->program($self->PROGRAM);
    $self->efg_analysis->program_file($self->PROGRAM_FILE);
    $self->efg_analysis->program_version($self->VERSION);
    $self->efg_analysis->parameters($self->PARAMETERS);
    $self->efg_analysis->displayable(1);

    my $logic_name = $self->efg_analysis->logic_name;
    #print Dumper $self->efg_analysis, $logic_name;

    ### check if analysis with this logic_name already exists in the database
    ### and has got same settings
    
    my $aa = $self->efgdb->get_AnalysisAdaptor;
    my $analysis = $aa->fetch_by_logic_name($logic_name);
    
    if ( ! defined $analysis ) { # NEW

        warn("Storing new analysis with logic name $logic_name.");
        $aa->store($self->efg_analysis);

    } elsif ( $self->efg_analysis->compare($analysis) ) { # UPDATE

        ### analysis compare
        # returns  1 if this analysis is special case of given analysis
        # returns  0 if they are equal
        # returns -1 if they are completely different

        warn('Analysis with logic name \''.$logic_name.'\' already '.
             'exists, but has different options! Updating analysis ...');

        $self->efg_analysis->dbID($analysis->dbID);
        $self->efg_analysis->adaptor($self->efgdb->get_AnalysisAdaptor);
        $aa->update($self->efg_analysis);
 
    } else { # EXISTS

        warn('Analysis with logic name \''.$logic_name.'\' already '.
             'exists.');

    }

    $self->efg_analysis($aa->fetch_by_logic_name($logic_name));
	
}

sub check_Sets {
  #This should use Utils::Helper::define_and_validate_sets

    my ($self) = @_;

    ### get result sets to process and check that result_set, data_set, 
    ### feature_set are storable

    my $rsets = $self->fetch_ResultSets();

    ### get feature_set/data_set name by determining the longest common 
    ### prefix of the result_set names
    #print Dumper @{$self->ResultSets()};
    my @names = map { $_->name; } @{$self->ResultSets()};
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
    
    my $set_name = $self->efg_analysis->logic_name.'_'.$lcp; #.'_'.$self->DATASET_NAME;
    print 'Set name: ', $set_name, "\n";
    
    my $fsa = $self->efgdb->get_FeatureSetAdaptor();
    my $fset = $fsa->fetch_by_name($set_name);
    #print Dumper $fset;

    if ( ! defined $fset ) {

        $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
            (
             -analysis      => $self->efg_analysis,
             -feature_type  => $self->feature_type,
             -cell_type     => $self->cell_type,
             -name          => $set_name,
             -feature_class => 'annotated'
             );
        #print Dumper $fset;
        
        warn("Storing new feature set \'$set_name\'");
        eval { 
            ($fset) = @{$fsa->store($fset)};
        };
        throw("Coudn't store feature set \'$set_name\': $!") if ($@);

    } else {

        warn("Feature set with name $set_name already exists.");

    }
    # save FeatureSet
    $self->feature_set($fset);

    my $dsa = $self->efgdb->get_DataSetAdaptor;
    my $dset = $dsa->fetch_by_name($set_name);
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

        warn("Storing new data set \'$set_name\'");
        eval { 
            ($dset) = @{$dsa->store($dset)} 
        };
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
    $self->data_set($dset);

}

=head2 fetch_ResultSets

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

    # fetch all available result_sets for given experiment and norm_method
    my $rsa = $self->efgdb->get_ResultSetAdaptor;
    my $rsets = $rsa->fetch_all_by_Experiment_Analysis
        ($self->experiment, $self->norm_method);
    print "No. of available ResultSets: ", scalar(@$rsets), "\n";

    # select only those result_sets that match regular expression
    my @rsets = ();
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
    my $norf;

    # make sure the cache directory exists
    my $cachedir = $self->ANALYSIS_WORK_DIR.'/cache';
    if ( ! -d $cachedir) {
        my $retval = system("mkdir -p $cachedir");
        throw("Couldn't create cache directory $cachedir") 
            unless ($retval == 0);
    }


    warn("INPUT_ID_TYPE: ". $self->analysis->input_id_type);
    
    my $input_id_type = $self->analysis->input_id_type;
    my $query_name;

    foreach my $rset (@{$self->ResultSets}) {

        print join(" ", $rset->dbID, $rset->name), "\n";

        if ($input_id_type eq 'slice') {
            
            $query_name = $self->query->[0]->name();

        } elsif ($input_id_type eq 'array') {
            
            my @echips = @{$rset->get_ExperimentalChips()};
            warn("WARNING: Result set '".$rset->name."' comprises more than one experimental chip, ".
                 "which will result in very large datafiles!")
                if (scalar @echips > 1);
            $query_name = 'ec_' . join (':', map { $_->unique_id } @echips);

        } else {

            throw("No idea how to deal with input_id type ".$input_id_type);

        }



        my $datfile = $cachedir.'/'.$query_name.'.'.$rset->name.'.dat';
        warn('datafile: '.$datfile);

        my @result_features = ();

        warn("SCORE_FACTOR: ".$self->SCORE_FACTOR);

        unless ( -e $datfile ) { # dat file already dumped; skipping sql

            foreach my $slice ( @{$self->query} ) {
            
                #my $result_features = $rset->get_ResultFeatures_by_Slice($self->query(),'DISPLAYABLE',1);
                #my $result_features = $rset->get_ResultFeatures_by_Slice($self->query(),'DISPLAYABLE');
                my $result_features = $rset->get_ResultFeatures_by_Slice($slice);

                print "No. of ResultFeatures_by_Slice:\t", scalar(@$result_features), "\n";
        
                if (scalar(@$result_features) == 0) {
                    warn("No result_features on slice ".$slice->name());
                    next;
                }

                my $ft_cnt = 1;
                
                foreach my $rft (sort {$a->start <=> $b->start} @{$result_features}) {
                    #print '<', join("><", 
                    #           $self->query()->seq_region_name, 
                    #           $rft->start, 
                    #           $rft->end, 
                    #           $rft->score*($self->SCORE_FACTOR), 
                    #           $rft->probe->get_probename), ">\n";
                    push (@result_features, 
                          [ 
                            $slice->seq_region_name,
                            $rft->start, 
                            $rft->end, 
                            $rft->score*$self->SCORE_FACTOR,
                            $ft_cnt++,
                            #$rft->probe->get_probename,
                            ]
                          );
                }
            }

            open(RF, "> $datfile")
                or throw("Can't open file $datfile");
            
            map {
                print RF join("\t", @$_),"\n";
            } @result_features;
            
            close RF
                or throw("Can't close file $datfile");
            
        } else {

            print "Using cached ResultFeatures:\t", $datfile, "\n";

            open(CACHE, $datfile)
                or throw("Can't open file cache $datfile");
            while (<CACHE>) {
                chomp;
                #print;

                my @values = split(/\t/);
                push (@result_features, 
                      [ 
                        @values
                        ]
                      );
            }           
            close CACHE;
            
        }

        
        warn("No. of ResultFeatures_by_Slice:\t".scalar(@result_features));

        throw ("No of result_features doesn't match") 
            if (@result_features && $norf && $norf != scalar(@result_features));

        $result_features{$rset->name} = \@result_features;
        
        $norf = scalar(@result_features);
        
    }

    if (scalar(keys %result_features) == 0) {

        warn('No ResultFeatures on slice');
        return 1;

    }
    
    #print Dumper %result_features;
    
    #$self->result_features(\%features);
    #print Dumper $self->result_features();
    #warn('EFG ANALYSIS PARAMETERS: '.$self->efg_analysis->parameters);

    my %parameters_hash = %{$self->parameters_hash($self->efg_analysis->parameters)};

    $parameters_hash{-result_features} = \%result_features;

    my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'
        .$self->efg_analysis->module;
    $runnable = $runnable->new
        (
         -query => $self->query,
         -program => $self->efg_analysis->program_file,
         -analysis => $self->efg_analysis,
         -workdir => $self->ANALYSIS_WORK_DIR,
         %parameters_hash
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
    
    print "Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::write_output\n";
    my ($self) = @_;

    if (scalar(@{$self->output}) == 0) {
        warn("No features to annotate on slice ".$self->query->name.
             " for experiment ".$self->experiment->name()."!");
        return 1;
    }
    
    # store analysis, feature set and data set
    # analysis was alredy been stored while checking config in read_and_check_config
    #$self->adaptor('Analysis')->store($self->analysis());
    if (! defined $self->feature_set->dbID) {
        $self->efgdb->get_FeatureSetAdaptor->store($self->feature_set());
    }
    if (! defined $self->data_set->dbID) {
        $self->efgdb->get_DataSetAdaptor->store($self->data_set());
    }
    
    ### annotated features

    #my $af = $self->adaptor('AnnotatedFeature')->fetch_all_by_Slice_FeatureSet
    #    ($self->query, $fset);

    my $fset = $self->feature_set;	
    my $fs_id = $fset->dbID();

    my @stored_af = ();
    foreach my $s (@{$self->query}) {

        push @stored_af, 
        @{$self->efgdb->get_AnnotatedFeatureAdaptor->fetch_all_by_Slice_FeatureSets($s, [$fset])};

    }

    warn('No. of annotated features already stored: '.scalar(@stored_af));
    warn('No. of annotated features to be stored:  '.scalar(@{$self->output}));

    if (@stored_af) {
        
        warn("NOT IMPORTING ".scalar(@{$self->output})." annotated features! Slice(s) ".
             join (', ', map {$_->name} @{$self->query})." already contains ".scalar(@stored_af).
             " annotated features of feature set '".$fset->name."' (dbId ".$fset->dbID.').');
        return 1;
        
    } else {
        
        # we need to transfer annotated features onto the toplevel slice, if the input
        # input slice is not a toplevel slice (starting with 1). This applies only if 
        # the input id is of type 'Slice'.
        
        my ($transfer, $slice, $tl_slice);
        my $input_id_type = $self->analysis->input_id_type;
        
        if ($input_id_type eq 'slice' && 
            ($self->query->[0]->start != 1 || $self->query->[0]->strand != 1)) {
            my $sa = $self->efgdb->get_SliceAdaptor;
            $tl_slice = $sa->fetch_by_region($self->query->[0]->coord_system->name(),
                                             $self->query->[0]->seq_region_name(),
                                             undef, #start
                                             undef, #end
                                             undef, #strand
                                             $self->query->[0]->coord_system->version());
            $transfer = 1;
            warn('TRANSFER: 1');
        }
        
        my %af_slice = ();
        
        if ($input_id_type eq 'array') {
            map { $af_slice{$_->seq_region_name} = $_ } @{$self->query};
        }
        
        my @af = ();
        foreach my $ft (@{$self->output}){
            
            #print Dumper $ft;
            my ($sr_name, $start, $end, $score) = @{$ft};
            
            $slice = $input_id_type eq 'slice' ? $self->query->[0] : $af_slice{"$sr_name"};
            
            #print Dumper ($sr_name, $start, $end, $score);
            my $af = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
                (
                 -slice         => $slice,
                 -start         => $start,
                 -end           => $end,
                 -strand        => 0,
                 -display_label => $self->efg_analysis->logic_name,
                 -score         => $score,
                 -feature_set   => $fset,
                 );
            
            # make sure feature coords are relative to start of entire seq_region
            if ($transfer) {
                #warn("original af:\t", join("\t", $af->start, $af->end), "\n");
                $af = $af->transfer($tl_slice);
                #warn("remapped af:\t", join("\t", $af->start, $af->end), "\n");
            }
            push(@af, $af);
        }
        
        $self->efgdb->get_AnnotatedFeatureAdaptor->store(@af);
    }
    return 1;
}


### container for config hash parameter

sub MODULE {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_MODULE'} = $value;
    }

    if ( exists( $self->{'_CONFIG_MODULE'} ) ) {
        return $self->{'_CONFIG_MODULE'};
    } else {
        return undef;
    }
}

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

sub SCORE_FACTOR {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_SCORE_FACTOR'} = $value;
    }

    if ( exists( $self->{'_CONFIG_SCORE_FACTOR'} ) ) {
        return $self->{'_CONFIG_SCORE_FACTOR'};
    } else {
        return undef;
    }
}


### 

=head2 dnadb

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : container for dbadaptor
  Returntype: Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions: throws if not passed a Bio::EnsEMBL::DBSQL::DBConnection
    object
  Example   : 

=cut


sub dnadb{
  my ($self, $db) = @_;
  if($db){
      throw("Must pass RunnableDB:db a Bio::EnsEMBL::DBSQL::DBAdaptor ".
            "not a ".$db) 
          unless($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
      $self->{'dnadb'} = $db;
  }
  return $self->{'dnadb'};
}

=head2 efgdb

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Function  : container for dbadaptor
  Returntype: Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Exceptions: throws if not passed a Bio::EnsEMBL::Funcgen::DBSQL::DBConnection
    object
  Example   : 

=cut


sub efgdb{
  my ($self, $db) = @_;
  if($db){
      throw("Must pass RunnableDB:db a Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor ".
            "not a ".$db) 
          unless($db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'));
      $self->{'efgdb'} = $db;
  }
  return $self->{'efgdb'};
}

sub efg_analysis {
    my ($self, $a) = @_;
    if ($a) {
      throw("Must pass a Bio::EnsEMBL::Analysis not a ".$a)
          unless($a->isa('Bio::EnsEMBL::Analysis'));         
      $self->{'efg_analysis'} = $a;      
    }
    return $self->{'efg_analysis'};
}

sub query {

    my ($self, $query) = @_;

    if ( $query ) {
        
        throw("Must pass RunnableDB:Funcgen:query a array ref not a ".
              ref($query)) unless (ref($query) eq 'ARRAY');

        map { 
            throw($_->name . " is not a Bio::EnsEMBL::Slice")
                unless ($_->isa("Bio::EnsEMBL::Slice"));
        } @$query;

        $self->{'query'} = $query;
    }

    return $self->{'query'};

}

sub experiment {
    my ($self, $e) = @_;
    if ($e) {
      throw("Must pass a Bio::EnsEMBL::Funcgen::Experiment not a ".$e)
          unless($e->isa('Bio::EnsEMBL::Funcgen::Experiment'));         
      $self->{'experiment'} = $e;      
    }
    return $self->{'experiment'};
}

sub norm_method {
    my ($self, $m) = @_;
    if ($m) {
      throw("Must pass a Bio::EnsEMBL::Ananysis not a ".$m)
          unless($m->isa('Bio::EnsEMBL::Analysis'));         
      $self->{'norm_method'} = $m;
    }
    return $self->{'norm_method'};
}

sub feature_type {
    my ($self, $ft) = @_;
    if ($ft) {
      throw("Must pass a Bio::EnsEMBL::Funcgen::FeatureType not a ".$ft)
          unless($ft->isa('Bio::EnsEMBL::Funcgen::FeatureType'));         
      $self->{'feature_type'} = $ft;
    }
    return $self->{'feature_type'};
}

sub cell_type {
    my ($self, $ct) = @_;
    if ($ct) {
      throw("Must pass a Bio::EnsEMBL::Funcgen::CellType not a ".$ct)
          unless($ct->isa('Bio::EnsEMBL::Funcgen::CellType'));         
      $self->{'cell_type'} = $ct;
    }
    return $self->{'cell_type'};
}

sub feature_set {
    my ($self, $fs) = @_;
    if ($fs) {
      throw("Must pass a Bio::EnsEMBL::Funcgen::FeatureSet not a ".$fs)
          unless($fs->isa('Bio::EnsEMBL::Funcgen::FeatureSet'));         
      $self->{'feature_set'} = $fs;
    }
    return $self->{'feature_set'};
}

sub data_set {
    my ($self, $ds) = @_;
    if ($ds) {
      throw("Must pass a Bio::EnsEMBL::Funcgen::DataSet not a ".$ds)
          unless($ds->isa('Bio::EnsEMBL::Funcgen::DataSet'));         
      $self->{'data_set'} = $ds;
    }
    return $self->{'data_set'};
}




sub ResultSets {
    my ($self, $sets) = @_;
    $self->{'ResultSets'} = $sets if ($sets);
    throw("No ResultSets in RunnableDB.") 
        if (!$self->{'ResultSets'});
    return $self->{'ResultSets'};
}



sub NORM_METHOD {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_NORM_METHOD'} = $value;
    }

    if ( exists( $self->{'_CONFIG_NORM_METHOD'} ) ) {
        return $self->{'_CONFIG_NORM_METHOD'};
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

sub EFGDB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_EFGDB'} = $value;
  }

  if ( exists( $self->{'_CONFIG_EFGDB'} ) ) {
    return $self->{'_CONFIG_EFGDB'};
  } else {
    return undef;
  }
}

sub DNADB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_DNADB'} = $value;
  }

  if ( exists( $self->{'_CONFIG_DNADB'} ) ) {
    return $self->{'_CONFIG_DNADB'};
  } else {
    return undef;
  }
}



1;
