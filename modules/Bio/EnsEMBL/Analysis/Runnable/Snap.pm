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



sub run_analysis{
  my ($self, $program) = @_;

  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable BaseAbInitio::run_analysis ") 
    unless($program && -x $program);

  my $command = $self->program." ".$self->matrix." ".$self->queryfile.
    " -gff -aa ".$self->protfile." > ".$self->resultsfile;
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
}




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


sub unaltered_slice{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'unaltered_slice'} = $arg;
  }
  return $self->{'unaltered_slice'};
}

1;
