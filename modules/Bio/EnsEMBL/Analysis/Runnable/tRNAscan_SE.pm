=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE - 

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE->new
  (
   -query => $slice,
   -program => 'trnascan-SE',
  );
  $runnable->run;
  my @simple_features = @{$runnable->output};

=head1 DESCRIPTION

tRNAscan_SE expects to run the program tRNAscan-SE and produces 
SimpleFeature which can be stored in the simple_feature table in the 
core database

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE
  Function  : produce a tRNAscan_SE runnable and set the options
  and program if undefined
  Returntype: Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE
  Exceptions: 
  Example   : 

=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  ######################
  #SETTING THE DEFAULTS#
  ######################
  if(!$self->options){
    $self->options('-q');
  }
  if(!$self->program){
    $self->program('tRNAscan-SE');
  }
  ######################
  return $self;
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::tRNAscan_SE
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
  if(!-e $results){
    throw("Can't parse an no existance results file ".$results.
          " tRNAscan_SE:parse_results");
  }
  my $ff = $self->feature_factory;
  my @output;
  open(CPG, $results) or throw("FAILED to open ".$results.
                               " tRNAscan_SE:parse_results");
 LINE:while(<CPG>){
    next LINE if(/^Sequence/ ||/^Name/ || /^---/); 
    #ignore introductory lines
    my @element = split (/\s+/, $_); 
    my ($name, $start, $end, $display_label, $score) 
      = @element[0, 2, 3, 4, 8];
    my $strand = 1;
    if($start > $end){
      $strand = -1;
      my $temp_end = $start;
      $start = $end;
      $end = $temp_end;
    }

    unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
      warning("The start and end were not integers, something went wrong. Skipping");
      next;
   }

    my $sf = $ff->create_simple_feature($start, $end, $strand, $score,
                                        $display_label, 
                                        $name, $self->query);
    push(@output, $sf)
  }
  $self->output(\@output);
  close(CPG) or throw("FAILED to close ".$results.
                      " tRNAscan_SE:parse_results");
}


1;
