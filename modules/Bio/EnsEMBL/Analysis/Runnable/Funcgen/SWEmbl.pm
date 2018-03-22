# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'program.pl',
  );
  $runnable->run;
  my @features = @{$runnable->output};

=head1 DESCRIPTION

SWEmbl expects bed or maq mapview files as input and predicts features which 
can be stored in the annotated_feature table in the eFG database

=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (get_file_format is_gzipped);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Funcgen);

=head2 run

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::SWEmbl
    Arg [2]     : filename
    Description : 
    Returntype  : 
    Exceptions  : 
    Example     : 

=cut

sub run {

    my ($self, $dir) = @_;

    $self->run_analysis;


    print "Parsing results ... ";
    $self->parse_results();#This is using the generic Runnable::Funcgen::parse_results!
    print "done!\n";

}


=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::SWEmbl
  Arg [2]     : string, program name
  Usage       : 
  Description : 
  Returns     : 
  Exceptions  : 

=cut

sub run_analysis {
        
    my ($self, $program) = @_;
    
    if(!$program){
        $program = $self->program;
    }

    throw($program." is not executable") if(! ($program && -x $program));

    my @fields = (1,2,6,9);
	#removed seq_region_name 0 as this is now matched by the regex 
	#in Runnable::Funcgen dependant no the input format
	#Could infact just extend the regex to do the output field assignment instead of the split

    $self->output_fields(\@fields);
	
	#Override default random file naming as we want to be able to look at the output easily
	#This will not make the directory if it does not exist!!
	#(my $resultsfile = $self->infile());

	#Need to checkdir exists anyway and also set to logic_name, not module!
	#This needs moving to Runnable::Funcgen?
	#And handling feature_set name where appropriate?
	
	
	my $results_dir = $self->workdir.'/'.$self->analysis->logic_name;
	my $suffix = &get_file_format($self->query);
	my $results_file = $self->query;
	$results_file =~ s/.*\///;	
	$suffix .= '.gz' if &is_gzipped($self->query, 1);
	$results_file =~ s/${suffix}$/dat/;
	$results_file = $results_dir.'/'.$results_file;
	

	if(! -d $results_dir){
	  system("mkdir -p $results_dir") && throw("Could not create results dir:\t$results_dir");
	}
	elsif(-e $results_file){
	  #This will stop us overwriting
	  throw("Results file already exists, please move away before running:\t$results_file");
	}


    $self->resultsfile($results_file);
    print "RESULTS FILE: ".$results_file."\n";
    
	my %fswitches = (
					 bed   => '-B',
					 sam   => '-S',
					 maq   => '-M',
					 eland => '-E',
					 #bam   => '-?',
					);
	my $format_switch = $fswitches{&get_file_format($self->query)};
	my $zip_switch = (&is_gzipped($self->query, 1)) ? '-z' : '';
		
	throw("Could not identify valid SWEmbl input format") if(! $format_switch);
	
    my $command = $self->program . " $format_switch $zip_switch -i " . $self->query . ' ' . 
      $self->options . ' -o ' . $results_file;

    if($self->has_control){ 
      #Identify the file name of the control based on the query name...
      #This assuming the control was added to the cache folder and conforms to the naming policy...
      if($self->query =~ /^(.*\/[^_]+)_[^_]+_(.+):.*$/){ 
	my $control_file = $1."_control_".$2.".samse.".&get_file_format($self->query).".gz";
	if(!(-e $control_file)){ throw("Could not find control file ".$control_file); }
	$command = $command." -r ".$control_file;
      } else {
	throw("Could not infer the name of the control file: Check the naming policy on ".$self->query);
      }
    }
    
    print "Running analysis:\t$command\n";
    system($command) && throw("FAILED to run $command: ", $@);
	print "Finished analysis:\t$command\n";

    # Here maybe post-process the output file to be sure it's ok for the following steps??
	
}


sub query{
  my $self = shift;
  $self->{'query'} = shift if(@_);

  throw("file ".$self->{'query'}. " doesn't exist") if (! -e $self->{'query'});

  return $self->{'query'};
}


sub config_file {
    my $self = shift;
    $self->{'config_file'} = shift if(@_);
    return $self->{'config_file'};
}

sub has_control {
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
