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

Bio::EnsEMBL::Analysis::Runnable::Snap - 

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::Snap->new(
      -query => $slice,
      -program => 'snap',
     );
  $runnable->run;
  my @predictions = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run the genefinder gene predictor and then parse the results
into prediction transcripts

=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Runnable::Snap;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;
use Bio::SeqIO;
use Bio::Seq;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->program('snap') if(!$self->program);
  ######################

  throw("Must defined a matrix file like ".
        "/usr/local/ensembl/Zoe/HMM/worm with the -matrix option ")
    if(!$self->matrix);
  $self->unaltered_slice($self->query);
  my $new_seq = Bio::Seq->new(
                              -seq => $self->query->seq,
                              -id => $self->query->name,
                             );
  $self->query($new_seq);
  return $self;
}




=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Snap
  Arg [2]   : string, program name
  Function  : create and open a commandline for the program genefinder
  Returntype: none
  Exceptions: throws if the program in not executable or the system
  command fails to execute 
  Example   : 

=cut

sub run_analysis{
  my ($self, $program) = @_;

  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable Snap::run_analysis ") 
    unless($program && -x $program);

  my $command = $self->program." ".$self->matrix." ".$self->queryfile.
    " > ".$self->resultsfile;
  print STDERR "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
}



=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Snap
  Arg [2]   : string, resultsfile name
  Function  : parse the results file into prediction exons then
  collate them into prediction transcripts and calculate their
  phases
  Returntype: none 
  Exceptions: throws if cant open or close results file or the parsing
  doesnt work
  Example   : 

=cut



sub parse_results{
  my ($self, $results) = @_;

  if(!$results){
    $results = $self->resultsfile;
  }

  $self->query($self->unaltered_slice);
  my $ff = $self->feature_factory;
  my $exon_count;
  my $current_trans;

  open(OUT, "<".$results) or throw("FAILED to open ".$results.
                                   "Snap:parse_results");

  while(<OUT>){
    /^\>/ and next;

    my @element = split;
    throw("Unable to parse Snap output ".@element." in output ".
          "array expecting 9") if(@element != 9);
    if (!$current_trans || 
        $current_trans ne $element[8]) {
	    $exon_count = 0;
	    $current_trans = $element[8];
    }
    $exon_count++;
    my $name = $current_trans.".".$exon_count;
    my $start = $element[1];
    my $end = $element[2];
    my $score = $element[4];
    throw("strand wrongly formated $element[6] not + or -")
      unless ($element[3] eq '+' || $element[3] eq '-');
    my $strand = 1;
    $strand = -1 if($element[3] eq '-');
    
    my $start_phase = (3 - $element[5]) % 3;

    my $exon = $ff->create_prediction_exon($start, 
                                           $end, 
                                           $strand, 
                                           $score, 
                                           0, 
                                           $start_phase,
                                           $name, 
                                           $self->query, 
                                           $self->analysis);
    $self->exon_groups($current_trans, $exon);
  }
  $self->create_transcripts;

  close(OUT) or throw("FAILED to close ".$results.
                      "Snap:parse_results");

}


=head2 unaltered_slice

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Snap
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : holder for the given slice as the one passed
  to Snap has to have desc removed otherwise parser doesnt
  work properly
  Returntype: Bio::EnsEMBL::Slice
  Exceptions: 
  Example   : 

=cut



sub unaltered_slice{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'unaltered_slice'} = $arg;
  }
  return $self->{'unaltered_slice'};
}

1;
