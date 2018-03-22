# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
1
=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( is_gzipped get_file_format );
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Funcgen);

# To do
# 1 Handle sam/bam files properly i.e. do not need sort them and create dat files or cache them?


=head2 new

  Arg [1]     : 
  Arg [2]     : 
  Description : Instantiates new SWEmbl runnabledb
  Returntype  : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl object
  Exceptions  : 
  Example     : 

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

	#my $file = '/nfs/users/nfs_n/nj1/scratch/peaks/rfbug_homo_sapiens_funcgen_58_37c/peaks_out/SWEmbl_DNase/CD4_DNase1_*';
	#my @tmp = `ls $file`;
	#throw("FOUND @tmp") if @tmp;

    $self->read_and_check_config($CONFIG);

    # make sure we have the correct analysis object
    $self->check_Analysis();

    # make sure we can store the correct feature_set, data_sets, and result_sets
    $self->check_Sets();

    return $self;
	
}

sub check_InputId {
    my ($self) = @_;
    my ($ename, $file) = split(':', $self->input_id);

    # set experiment
    my $ea = $self->efgdb->get_ExperimentAdaptor;
    my $e = $ea->fetch_by_name($ename);

    my @date = (localtime)[5,4,3];
    $date[0] += 1900; $date[1]++;
    
    if (! defined $e) {

        warn("Experiment NOT defined");

        my $exp = Bio::EnsEMBL::Funcgen::Experiment->new
            (
             -NAME => $ename,
             -GROUP => "$ENV{EFG_GROUP}",
             -DATE => join('-', @date),
             -PRIMARY_DESIGN_TYPE => 'binding_site_identification',
             -ADAPTOR => $ea,
             );

        my ($g_dbid) = $self->efgdb->fetch_group_details($exp->group());
		
        if (! $g_dbid) {
            
	    #This is causing concurrency issues... move this to a InnoDB table?... create GroupAdaptor class??
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
        or throw("Can't fetch experiment with name ".$ename);

    my $infile = $self->ANALYSIS_WORK_DIR.'/infiles/'.$self->input_id;
    warn('INFILE: '.$infile);
    $self->query($infile);

}
   
sub check_Sets {

    warn("Analysis::RunnableDB::Funcgen::SWEmbl::check_Experiment - should use EFGUtils::define_and_validate_sets");
	#Also, resolve wrt to Funcgen::check_Sets

    my $self = shift;

	my $iset_name = $self->experiment->name;
    my $set_name  = $iset_name.'_'.$self->analysis->logic_name;
	#Remove logic name from set_name?
	#If we need to discrimate between different analysis runs then we should
	#Or at least move to end of name?
	#Set a different exp_suffix via CreateAlignIDs
	#Will this fail or overwrite if it is already found?
	#Using define_and_validate_set will allow rollback of previously called features(given an env var in the config)


    my ($ct_name, $ft_name) = split(/_/, $iset_name);

    my $feature_type = $self->efgdb->get_FeatureTypeAdaptor()->fetch_by_name($ft_name)
        or throw("Feature type '$ft_name' does not exist");
    $self->feature_type($feature_type);

    my $cell_type = $self->efgdb->get_CellTypeAdaptor()->fetch_by_name($ct_name)
        or throw("Cell type '$ct_name' does not exist");
    $self->cell_type($cell_type);

    my $isa = $self->efgdb->get_InputSetAdaptor();
    my $iset = $isa->fetch_by_name($iset_name);

	#Strip the Experiment name off the input_id to give the
	#file name for the input_subset.name
	my $file_name;
	($file_name = $self->input_id) =~ s/.*://;
	
    
    if (! defined $iset){

	  $iset = Bio::EnsEMBL::Funcgen::InputSet->new
            (
             -name         => $iset_name,
             -experiment   => $self->experiment(),
             -feature_type => $feature_type,
             -cell_type    => $cell_type,
             -vendor       => 'SOLEXA',
             -format       => 'SEQUENCING',
			 -feature_class => 'result'
             #-analysis     => $self->feature_analysis,
			);
	  print "Storing new InputSet:\t$iset_name\n";
	  ($iset)  = @{$isa->store($iset)};
	  $iset->add_new_subset($file_name);
	  $iset->adaptor->store_InputSubsets($iset->get_InputSubsets);
	}
	else{
	  #We only expect one subset here
	  print "InputSet already exists:\t$iset_name\n";
	  my @issets = @{$iset->get_InputSubsets};

	  if(scalar(@issets) > 1){
		throw("InputSet $iset_name has more than one InputSubset:\t".join("\t", (map $_->name, @issets)));
	  }
	  elsif((scalar(@issets) == 1) &&
			($issets[0]->name ne $file_name)){
		throw("InputSet $iset_name already has an InputSubset(".$issets[0]->name.") which does not match ".$file_name);
	  }
	  elsif(scalar(@issets) == 0){#we can just add this InputSubset
		$iset->add_new_subset($file_name);
		$iset->adaptor->store_InputSubsets($iset->get_InputSubsets);
	  }
	}



    my $fsa = $self->efgdb->get_FeatureSetAdaptor();
    my $fset = $fsa->fetch_by_name($set_name);


    if ( ! defined $fset ) {

        $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
            (
             -analysis      => $self->efg_analysis,
             -feature_type  => $self->feature_type,
             -cell_type     => $self->cell_type,
             -name          => $set_name,
             -feature_class => 'annotated'
             );

		print "Storing new FeatureSet:\t$set_name\n";
		($fset) = @{$fsa->store($fset)};
		
    } 
	else {
	  print "FeatureSet already exists:\t$set_name\n";
    }
    # save FeatureSet
    $self->feature_set($fset);

    my $dsa = $self->efgdb->get_DataSetAdaptor;
    my $dset = $dsa->fetch_by_name($set_name);


    if ( ! defined $dset ) {

	  $dset = Bio::EnsEMBL::Funcgen::DataSet->new
		(
		 -SUPPORTING_SETS     => [$iset],
		 -FEATURE_SET         => $fset,
		 -DISPLAYABLE         => 1,
		 -NAME                => $set_name,
             -SUPPORTING_SET_TYPE => 'input',
		);
	  
	  print "Storing new DataSet:\t$set_name\n";
	  ($dset) = @{$dsa->store($dset)}
	} 
	else {

	  print "DataSet already exists:\t$set_name\n";
	  
	  # need to check whether InputSets and supporting_sets are the same and 
	  # possibly add InputSet to supporting_sets

	  my $ssets = $dset->get_supporting_sets();

	  my %ssets_dbIDs = ();
	  map { $ssets_dbIDs{$_->dbID}='' } (@{$ssets});
	  $dset->add_supporting_sets([ $iset ]) if (! exists $ssets_dbIDs{$iset->dbID}); 
	  
    }
 
    $self->data_set($dset);
}

sub fetch_input {
  my ($self) = @_;
   
  # make sure the cache directory exists
  my $cachedir = $self->ANALYSIS_WORK_DIR.'/cache';

  if ( ! -d $cachedir) {
	system("mkdir -p $cachedir") && throw("Couldn't create cache directory $cachedir");
  }

  print "Input file:\t".$self->query."\n";
    
  (my $cachefile = $self->query) =~ s,.*/([^/]+)$,$1,;
  $cachefile = $cachedir.'/'.$cachefile;
  print "Filtered file:\t".$cachefile."\n";
  

  

  

  #Should really set file format here, so we don't have to test in the Runnable
  #Also required in the default Runnable::Funcgen::parse_results
  #This should beset in the file parser!
  my $file_format = &get_file_format($self->query);
  $self->{'input_format'} =  $file_format;#HACK!!!!! Write input_format method in RunnableDB::Funcgen.pm?
  my $sort;


  unless (-e $cachefile) { # dat file already dumped sorted; skipping dump
	#What if file is incomplete??? We need an md5 here!

	# make sure reads are dumped in sorted order
	#Explain this sort key?
	#Change this to simply pipe the gzip sort as the input file for SWEmbl
	#Need to check is compressed also
	#Implement Bed(other) parser here?
	#Why bother dumping, why don't we just use this as this sort as the input file handle?
	
	#These sorts should be in the Bed/Sam file parsers
	#Implement file parsers in Runnables(Need to separate parsers from importer functions first)
	#Or wait and use BioPerl? Then just have simple wrappers to include our our config/methods sat on top 
	#of the BioPerl parsers

	if($file_format eq 'bed'){
	  $sort = "gzip -dc ".$self->query.
		' | grep -vE "^MT" | sort -k1,1 -k2,2n -k3,3n | gzip -c > '.$cachefile.' |';
	}
	elsif($file_format eq 'sam'){
	  #We need to preserve the header line beginning with @
	  #But remove and unmapped reads with bitwise flag of 4 as 2nd field

	  #Quick fix is to just rmeove the headers for now
	  #Could really do with removing them, then sorting and adding back in.


	  warn("This does not sort sam by end position! Hence will not work for paired ends reads. Need to implement samtools sort on bam files");
	  #No other bits set for unmapped?
	  #Not necessarily true!
	  #Need to bitwise and this value!
	  $sort = "gzip -dc ".$self->query.' | grep -vE \'^@\' | grep -vE "^[^[:space:]]+[[:blank:]]4[[:blank:]].*$" | grep -vE "^[^[:space:]]+[[:blank:]][^[:space:]]+[[:blank:]][^[:space:]]+\:[^[:space:]]+\:MT\:" |sort -k3,3 -k4,4n | gzip -c > '.$cachefile.' |';
	}
	else{
	  throw("$file_format file format not supported");
	}

	warn("Executing  $sort");
	open(SORT, "$sort") or throw("Can't open and sort gzipped file ".$self->query);
	while (<SORT>) {warn($_)}
	close SORT;
  }

  $self->query($cachefile);
  my %parameters_hash = %{$self->parameters_hash($self->efg_analysis->parameters)};

  my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'.$self->efg_analysis->module;
  $runnable = $runnable->new(
							 -query      => $self->query,
							 -program    => $self->efg_analysis->program_file,
							 -analysis   => $self->efg_analysis,
							 -workdir    => $self->ANALYSIS_WORK_DIR,
							 #-input_format => '$file_format',#? WOuld need to add to Runnable::Funcgen.pm
							 %parameters_hash
							);
    
  $self->runnable($runnable);
    
  return;   
}

sub query {
  my $self = shift;
  $self->{'query'} = shift if(@_);

  throw("file ".$self->{'query'}. " doesn't exist") 
      unless ( -e $self->{'query'});

  return $self->{'query'};
}

sub write_output{
  #Is this not the same as Funcgen::write_output?
  #There is only one way to write AnnotatedFeatures!


  
    print "Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl::write_output\n";
    my ($self) = @_;

	#Where is this getting set?

	#The following warns should be throws!
	#But we need to up the retry count so that we don't run the job again.


    if (scalar(@{$self->output}) == 0) {
	  #Surely this should be throw?
        warn("No features called from ".$self->query.
             " for experiment ".$self->experiment->name()."!");
        return 1;
    }

    # store analysis, feature set and data set
    # analysis was alredy been stored while checking config in read_and_check_config
    #$self->adaptor('Analysis')->store($self->analysis());

    #This seems to be causing conflict... 
    if (! defined $self->feature_set->dbID) {
        $self->efgdb->get_FeatureSetAdaptor->store($self->feature_set());
    }
    if (! defined $self->data_set->dbID) {
        $self->efgdb->get_DataSetAdaptor->store($self->data_set());
    }
    
    ### annotated features
    
    my $fset = $self->feature_set;	
    my $fs_id = $fset->dbID();
    
    my $af = $self->efgdb->get_AnnotatedFeatureAdaptor->fetch_all_by_FeatureSets([$fset]);
    
    warn('No. of annotated features already stored: '.scalar(@$af).' ('.$self->query.' '.$fset->name.')');
    warn('No. of annotated features to be stored: '.scalar(@{$self->output}).' ('.$self->query.')');

    #??!!! Why is this done here and not before SWEmbl is actually run?
    #Maybe so we can rerun with different params without storing?
    #This should be done with different logic name and specify a no write env var(or pipeline opt)?
    

    if (@$af) {
      #Surely this should be throw?
      warn("NOT IMPORTING ".scalar(@{$self->output})." annotated features! File ".
	   $self->query." already has been processed; contains ".scalar(@$af).
	   " annotated features of feature set ".$fset->dbID.".");
      return 1;
      
    } else {
        my @af;

        my $sa = $self->efgdb->get_SliceAdaptor();
        my %slice;
        
        foreach my $ft (@{$self->output}){

	  my ($seqid, $start, $end, $score, $summit) = @{$ft};

	  #warn "Debugging Region: ".$seqid." Start: ".$start." End: ".$end." Score: ".$score." Summit: ".$summit."\n";
	  
	  #Hack mostly to ignore header (should be in parse output somewhere?)
	  if(! (($start =~ /^-?\d+$/) && ($end =~ /^\d+$/))){  
	    warn "Feature being ignored due to incorect positions: Region:".$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
	    
	    next;
	  }
		  

	  $summit = int($summit);#Round down
	  
	  # skip mito calls
	  #remove this when we have pre-processing step to filter alignments using blacklist?
	  next if ($seqid =~ m/^M/);
	  
	  unless (exists $slice{"$seqid"}) {
	    $slice{"$seqid"} = $sa->fetch_by_region(undef, $seqid);
	  }

	  #Sometimes there are naming issues with the slices... e.g. special contigs... which are not "valid" slices in ENSEMBL
	  if(!$slice{"$seqid"}){
	    warn "Feature being ignored due to inexistent slice: Region:".$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
	    next;
	  }
	  
	  my $af = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
	    (
	     -slice         => $slice{"$seqid"},
	     -start         => $start,
	     -end           => $end,
	     -strand        => 0,
	     -score         => $score,
		 -summit        => $summit,
	     -feature_set   => $fset,
	    );
	  
	  push(@af, $af);
	  
	}
        
        $self->efgdb->get_AnnotatedFeatureAdaptor->store(@af);
    }

	#This should all be done in Funcgen.pm?
	#For all RunnableDBs

	$fset->adaptor->set_imported_states_by_Set($fset);
	#Question about whether we need an imported status on the Input(Sub)Set here too?
	#We need to idnetify which subsets have been imported as we lose this info once
	#a file has been imported
	#however, what defines imported?
	#on which assembly
	#using which analysis

	#IMPORTED status only denotes that some import has been completed succesfully
	#IMPORTED_CSVERSION denotes which assembly
	#So IMPORTED is only appropriate for none CS related data i.e. array designs
	#where as IMPORTED_CSVERSION should be used for all coord dependant data
	#This does not account for analysis specific imports
	#i.e. want to reimport an inputset using a different analysis
	#Will the current code just ignore files which have been imported on another analysis?
	


    return 1;
    
}


#=head2 run
#  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::SWEmbl
#  Function  : overrides the general RunnableDB run just to pass an extra info to the runnable(s)... 
#  Returntype: array ref
#  Exceptions: none
#=cut

sub run{
  my ($self) = @_;

  foreach my $runnable(@{$self->runnable}){
    $runnable->has_control($self->HAS_CONTROL);
    $runnable->run;
    $self->output($runnable->output);
  }
  return $self->{'output'};
}

#Maybe pass it to more generic Funcgen? (other runnables may need this)
sub HAS_CONTROL {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_HAS_CONTROL'} = $value;
    }

    if ( exists( $self->{'_CONFIG_HAS_CONTROL'} ) ) {
        return $self->{'_CONFIG_HAS_CONTROL'};
    } else {
        return undef;
    }
}

1;
