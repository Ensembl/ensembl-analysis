# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Snap
#
# Copyright (c) 2005 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Snap

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

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

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
    " -gff -aa ".$self->protfile." > ".$self->resultsfile;
  print "Running analysis ".$command."\n";
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

   my $in  = Bio::SeqIO->new ( 
                              '-format' => 'Fasta' , 
                              -file => $self->protfile
                             );
  
  my %peptides;
  while(my $seq = $in->next_seq) {
    $peptides{$seq->id} = $seq->seq;
  }
  $self->peptides(\%peptides);
  open(OUT, "<".$results) or throw("FAILED to open ".$results.
                                   "Snap:parse_results");
  $self->query($self->unaltered_slice);
  my $ff = $self->feature_factory;
  my $exon_count;
  my $current_trans;
  while(<OUT>){
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
    my $start = $element[3];
    my $end = $element[4];
    my $score = $element[5];
    throw("strand wrongly formated $element[6] not + or -")
      unless ($element[6] eq '+' || $element[6] eq '-');
    my $strand = 1;
    $strand = -1 if($element[6] eq '-');
    my $exon = $ff->create_prediction_exon($start, $end, $strand, 
                                           $score, 0, 0, 
                                           $name, $self->query, 
                                           $self->analysis);
    $self->exon_groups($current_trans, $exon);
  }
  $self->create_transcripts;
  $self->calculate_phases;
  close(OUT) or throw("FAILED to close ".$results.
                      "Snap:parse_results");

}


=head2 protfile

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Snap
  Arg [2]   : string, filename
  Function  : accessor method for prot file name, will create
  one if one is requested but not defined
  Returntype: string, filename
  Exceptions: 
  Example   : 

=cut



sub protfile{
  my ($self, $filename) = @_;

  if($filename){
    $self->{'protfile'} = $filename;
  }
  if(!$self->{'protfile'}){
    $self->{'protfile'} = $self->create_filename('snap', 'prot');
  }
  $self->files_to_delete($self->{'protfile'});
  return $self->{'protfile'};
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
