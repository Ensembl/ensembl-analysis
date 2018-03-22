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

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::EponineTSS - 

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::EponineTSS->new
  (
   -query => $slice,
   -program => 'java',
  );
  $runnable->run;
  my @simple_features = @{$runnable->output};

=head1 DESCRIPTION

EponineTSS expects to run the program EponineTSS and produces SimpleFeature
which can be stored in the simple_feature table in the core database

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::EponineTSS;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EponineTSS
  Arg [2]   : string, path to the eponine jar
  Arg [3]   : int, threshold
  Function  : create a new Bio::EnsEMBL::Analysis::Runnable::EponineTSS
  Returntype: Bio::EnsEMBL::Analysis::Runnable::EponineTSS
  Exceptions: none
  Example   :

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($epojar, $threshold) = rearrange(['EPOJAR', 'THRESHOLD'], @args);

  if(!$self->program){
    $self->program('java');
  }
  
  if ($epojar) {
  	$self->epojar($epojar);
  } else {
    $self->epojar('eponine-scan.jar');	
  }
  	
  if ($threshold) {
    $self->threshold($threshold);
  } else {
    $self->threshold(50);
  }

  return $self;
}



=head2 epojar

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EponineTSS
  Arg [2]   : string , path to the eponine jar
  Function  : container
  Returntype: string
  Exceptions: none (find file will throw if it cant locate the file)
  Example   :

=cut


sub epojar{
  my $self = shift;
  my $epojar = shift;
  if($epojar){
    if(! -e $epojar){
      my $temp = $self->find_file($epojar);
      $epojar = $temp;
    }
    $self->{'epojar'} = $epojar;
  }
  return $self->{'epojar'};
}



=head2 threshold

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EponineTSS
  Arg [2]   : int, threshold
  Function  : container
  Returntype: int
  Exceptions:
  Example   :

=cut


sub threshold{
  my $self = shift;
  $self->{'threshold'} = shift if(@_);
  return $self->{'threshold'};
}



=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EponineTSS
  Arg [2]   : string, program name
  Function  : construct and open commandline for running
  eponine
  Returntype: none
  Exceptions: throws if fail to run program
  Example   :

=cut



sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable EponineTSS::run_analysis ")
    unless($program && -x $program);
  my $command = $program." ";
  $command .= $self->options." " if($self->options);
  $command .= "-jar ".$self->epojar." -seq ".$self->queryfile.
    " -threshold ".$self->threshold." > ".$self->resultsfile;
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EponineTSS
  Arg [2]   : string, filename
  Function  : to open and parse the results file
  Returntype: none
  Exceptions: throws on failure to open or close the results file
  or if the results file doesnt exist
  Example   :

=cut




sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  my $ff = $self->feature_factory;
  if(!-e $results){
    throw("Can't parse an no existance results file ".$results.
          " EponineTSS:parse_results");
  }
  my @output;
  open(EponineTSS, $results) or throw("FAILED to open ".$results.
                                      " EponineTSS:parse_results");
 LINE:while(<EponineTSS>){
    if (! /^\#/){ #ignore introductory lines
      chomp;
      my @element = split;
      my ($name, $start, $end, $score, $temp_strand) =
        @element[0, 3, 4, 5, 6];
      my $strand = 1;
      # Eponine is reporting coordinate with a one base offset when predicting on the reverse strand
      if($temp_strand eq '-'){
        $strand = -1;
        $start--;
        $end--;
      }
      # This should not happen anymore as we fix the offset but I will leave it here for some time
      if ($self->query->seq_region_length < $start) {
          warning("WRONG feature $start $end $score $temp_strand as length is ".$self->query->seq_region_length."\n");
          next if ($start == $end);
      }
      $score = $self->trunc_float_3($score);
      my $sf = $ff->create_simple_feature($start, $end, $strand, $score,
                                          '', $name, $self->query);
      push(@output, $sf);
    }
  }
  $self->output(\@output);
  close(EponineTSS) or throw("FAILED to close ".$results.
                             " EponineTSS:parse_results");
}



=head2 trunc_float_3

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::EponineTSS
  Arg [2]   : int, number to be truncated
  Function  : truncate a floating point number to three decimal places
  Returntype:
  Exceptions:
  Example   :

=cut



sub trunc_float_3 {
    my ($self, $arg) = @_;

    # deal only with valid numbers
    # and only need cases of the form [+/-]xx.yyyyy
    return $arg unless $arg =~ /^[+-]?\d*\.\d+$/;

    return 0.001 * int (1000 * $arg);
}
