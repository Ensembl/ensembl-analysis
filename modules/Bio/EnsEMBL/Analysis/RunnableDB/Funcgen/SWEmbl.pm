# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl
#
# Copyright (c) 2008 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl->new
     (
         -db       => $db,
         -input_id => 'chromosome::20:1:100000:1',
         -analysis => $analysis,
     );
  $runnable->fetch_input;
  $runnable->run;
  $runnable->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl functional genomics 
database and the Runnable SWEmbl which wraps the ChIP-Seq peak caller SWEmbl.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Funcgen::SWEmbl;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA); 

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Funcgen);

=head2 new

  Arg [1]     : 
  Arg [2]     : 
  Description : Instantiates new SWEmbl runnabledb
  Returntype  : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl object
  Exceptions  : 
  Example     : 

=cut

sub new {

    print "Analysis::RunnableDB::Funcgen::SWEmbl::new\n";
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);


    $self->read_and_check_config($CONFIG);

    # make sure we have the correct analysis object
    $self->check_Analysis();

    # make sure we can store the correct feature_set, data_sets, and result_sets
    $self->check_Sets();

    return $self;
	
}

sub check_Experiment {

    warn("Analysis::RunnableDB::Funcgen::SWEmbl::check_Experiment");

    my ($self) = @_;

    # extract experiment name and query file from input_id
    warn("INPUT_ID:\t".$self->input_id);
    my ($experiment, $file) = split(':', $self->input_id);

    
    my $infile = $self->ANALYSIS_WORK_DIR.'/infiles/'.$file;
    warn('INFILE: '.$self->ANALYSIS_WORK_DIR.'/infiles/'.$file);
    $self->query($infile);

    my $ea = $self->efgdb->get_ExperimentAdaptor();
    my $e = $ea->fetch_by_name($experiment);

    my @date = (localtime)[5,4,3];
    $date[0] += 1900; $date[1]++;
    #print Dumper @date;
    
    if (! defined $e) {

        warn("Experiment NOT defined");

        my $exp = Bio::EnsEMBL::Funcgen::Experiment->new
            (
             -NAME => $experiment,
             -GROUP => "$ENV{EFG_GROUP}",
             -DATE => join('-', @date),
             -PRIMARY_DESIGN_TYPE => 'binding_site_identification',
             -ADAPTOR => $ea,
             );

        my ($g_dbid) = $self->efgdb->fetch_group_details($exp->group());

        if (! $g_dbid) {
            
            warn("Group specified does, not exist. Importing (group, location, contact)");
            my $sql = "INSERT INTO experimental_group (name, location, contact) ".
                "values ('$ENV{EFG_GROUP}','$ENV{EFG_LOCATION}','$ENV{EFG_CONTACT}')";
            eval {
                $self->efgdb->dbc->do($sql);
            };
            throw("Couldn't import group information. Double-check that the environment variables ".
                  "EFG_GROUP, EFG_LOCATION, and EFG_CONTACT are set.") if ($@);
            
        }

        ($e) =  @{$ea->store($exp)};

    }

    $self->experiment($e) 
        or throw("Can't fetch experiment with name ".$experiment);
    #print Dumper $self->experiment;

}
   
sub check_Sets {

    warn("Analysis::RunnableDB::Funcgen::SWEmbl::check_Experiment");

    my ($self) = @_;

    my $set_name = $self->experiment->name();
    my ($ct_name, $ft_name) = split(/_/, $set_name);
    ###
    #### Hack to import LMI methylation data again ####
    ###
    $set_name .= '_e'.$ENV{VERSION};

    my $feature_type = $self->efgdb->get_FeatureTypeAdaptor()->fetch_by_name($ft_name)
        or throw("Feature type $ft_name does not exit");
    $self->feature_type($feature_type);

    my $cell_type = $self->efgdb->get_CellTypeAdaptor()->fetch_by_name($ct_name)
        or throw("Cell type $ct_name does not exit");
    $self->cell_type($cell_type);

    my $esa = $self->efgdb->get_ExperimentalSetAdaptor();
    my $eset = $esa->fetch_by_name($set_name);
    
    if (! defined $eset){

        warn("ExperimentalSet NOT defined");

        $eset = Bio::EnsEMBL::Funcgen::ExperimentalSet->new
            (
             -name         => $set_name,
             -experiment   => $self->experiment(),
             -feature_type => $feature_type,
             -cell_type    => $cell_type,
             -vendor       => 'SOLEXA',
             -format       => 'SEQUENCING',
             #-analysis     => $self->feature_analysis,
             );
        
        warn("Storing new experimental set \'$set_name\'");
        eval { 
            ($eset)  = @{$esa->store($eset)};
        };
        throw("Coudn't store experimental set \'$set_name\': $!") if ($@);
    }

    my $fsa = $self->efgdb->get_FeatureSetAdaptor();
    my $fset = $fsa->fetch_by_name($set_name);
    #print Dumper $fset;

    if ( ! defined $fset ) {

        $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
            (
             -analysis => $self->efg_analysis,
             -feature_type => $self->feature_type,
             -cell_type => $self->cell_type,
             -name => $set_name,
             -type => 'annotated'
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
             -SUPPORTING_SETS     => [$eset],
             -FEATURE_SET         => $fset,
             -DISPLAYABLE         => 1,
             -NAME                => $set_name,
             -SUPPORTING_SET_TYPE => 'experimental',
             );
        #print Dumper $dset;

        warn("Storing new data set \'$set_name\'");
        eval { 
            ($dset) = @{$dsa->store($dset)}
           };
        throw("Coudn't store data set \'$set_name\': $!") if ($@);

    } else {

        warn("Data set with name $set_name already exists.");

        # need to check whether ExperimentalSets and supporting_sets are the same and 
        # possibly add ExperimentalSet to supporting_sets

        my $ssets = $dset->get_supporting_sets();

        my %ssets_dbIDs = ();
        map { $ssets_dbIDs{$_->dbID}='' } (@{$ssets});
        $dset->add_supporting_sets([ $eset ]) if (! exists $ssets_dbIDs{$eset->dbID}); 

    }
    # save DataSet
    $self->data_set($dset);

}

sub fetch_input {
    
    my ($self) = @_;
    #warn ("Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl::fetch_input");
 
    # make sure the cache directory exists
    my $cachedir = $self->ANALYSIS_WORK_DIR.'/cache';
    if ( ! -d $cachedir) {
        my $retval = system("mkdir -p $cachedir");
        throw("Couldn't create cache directory $cachedir") unless ($retval == 0);
    }

    warn('infile: '.$self->query);

    ### Checking for gzip and bed format has been moved to create_input_id stage!!!
    #my $open = ($self->is_gzip($self->query))? 'gzip -dc' : 'cat';
    #my $cmd = $open.' '.$self->query.' |';
    #warn $cmd;
    ### bed file format ?
    #open(INFILE, "$cmd")
    #    or throw("Can't open gzipped infile ".$self->query);
    #while (<INFILE>) {
    #
    #    my @line = split("\t");
    #    throw("Infile does not have 6 or more columns. We expect bed format: CHROM START END NAME SCORE STRAND.") 
    #        if scalar @line < 6;
    #    
    #    throw ("1st column must contain name of seq_region (e.g. chr1 or 1)")
    #        unless ($line[0] =~ m/^(chr[MTXY\d]+)$/);
    #    throw ("2nd and 3rd column must contain start and end respectively")
    #        unless ($line[1] =~ m/^\d+$/ && $line[2] =~ m/^\d+$/);
    #    throw ("6th column must define strand (either '+' or '-')")
    #        unless ($line[5] =~ m/^[+-]$/);
    #    last;
    #}
    #close INFILE;
    
    (my $cachefile = $self->query) =~ s,.*/([^/]+)$,$1,;
    my $datfile = $cachedir.'/'.$cachefile;
    warn('datafile: '.$datfile);
    
    unless (-e $datfile) { # dat file already dumped sorted; skipping dump

        # make sure reads are dumped in sorted order
        my $sort = "gzip -dc ".$self->query.
            " | sort -k1,1 -k2,2n -k3,3n | gzip -c > $datfile |";
        warn("Executing  $sort");
        open(SORT, "$sort")
            or throw("Can't open and sort gzipped file ".$self->query);
        while (<SORT>) {warn($_)}
        close SORT;
        
    }

    $self->query($datfile);

    my %parameters_hash = %{$self->parameters_hash($self->efg_analysis->parameters)};

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

sub query {
  my $self = shift;
  $self->{'query'} = shift if(@_);

  throw("file ".$self->{'query'}. " doesn't exist") 
      unless ( -e $self->{'query'});

  return $self->{'query'};
}

sub write_output{
    
    print "Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl::write_output\n";
    my ($self) = @_;

    if (scalar(@{$self->output}) == 0) {
        warn("No features to annotate on slice ".$self->query.
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
    
    my $fset = $self->feature_set;	
    my $fs_id = $fset->dbID();
    
    my $af = $self->efgdb->get_AnnotatedFeatureAdaptor->fetch_all_by_FeatureSets(($fset));
    
    warn('No. of annotated features already stored: '.scalar(@$af).' ('.$self->query.' '.$fset->name.')');
    warn('No. of annotated features to be stored: '.scalar(@{$self->output}).' ('.$self->query.')');

    if (@$af) {
        
        warn("NOT IMPORTING ".scalar(@{$self->output})." annotated features! File ".
             $self->query." already has been processed; contains ".scalar(@$af).
             " annotated features of feature set ".$fset->dbID.".");
        return 1;
        
    } else {
        my @af;

        my $sa = $self->efgdb->get_SliceAdaptor();
        my %slice;
        
        foreach my $ft (@{$self->output}){
            
            #print Dumper $ft;
            my ($seqid, $start, $end, $score) = @{$ft};
            #print Dumper ($seqid, $start, $end, $score);
            
            # skip mito calls
            next if ($seqid =~ m/^M/);

            unless (exists $slice{"$seqid"}) {
                
                $slice{"$seqid"} = $sa->fetch_by_region('chromosome', $seqid);
                
            }
     
            my $af = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
                (
                 -slice         => $slice{"$seqid"},
                 -start         => $start,
                 -end           => $end,
                 -strand        => 0,
                 -display_label => $self->efg_analysis->logic_name,
                 -score         => $score,
                 -feature_set   => $fset,
                 );
            
            push(@af, $af);

        }
        
        $self->efgdb->get_AnnotatedFeatureAdaptor->store(@af);
    }
    return 1;
    
}


1;
