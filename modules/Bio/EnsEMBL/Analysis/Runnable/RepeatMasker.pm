# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::RepeatMasker
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::RepeatMasker

=head1 SYNOPSIS

  my $repeat_masker = Bio::EnsEMBL::Analysis::Runnable::RepeatMasker->
  new(
      -query => 'slice',
      -program => 'repeatmasker',
      -options => '-low'
     );
  $repeat_masker->run;
  my @repeats = @{$repeat_masker->output};

=head1 DESCRIPTION

RepeatMasker expects to run the program RepeatMasker
and produce RepeatFeatures which can be stored in the repeat_feature
and repeat_consensus tables in the core database


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Runnable::RepeatMasker;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::RepeatConsensus;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::RepeatMasker
  Function  : constructs a commandline and runs the program passed
  in, the generic method in Runnable isnt used as RepeatMasker doesnt
  fit this module
  Returntype: none
  Exceptions: throws if run failed because sysetm doesnt
  return 0 or the output file doesnt exist
  Example   : 

=cut


sub run_analysis{
  my ($self) = @_;
  my $cmd = $self->program." ";
  $cmd .= $self->options." " if($self->options);
  $cmd .= $self->queryfile;
  print "Running analysis ".$cmd."\n";
  system($cmd) == 0 or throw("FAILED to run ".$cmd.
                             " RepeatMasker:run_analysis");
  foreach my $file(glob $self->queryfile."*"){
    $self->files_to_delete($file);
  }
  if(! -e $self->resultsfile){
    throw("FAILED to run repeatmasker on ".$self->queryfile." ".
          $self->resultsfile." has not been produced ".
          "RepeatMasker:run_analysis");
  }
}



=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::RepeatMasker
  Arg [2]   : string, filename
  Function  : open and parse the results file into repeat
  features
  Returntype: none 
  Exceptions: throws on failure to open or close output file
  Example   : 

=cut



sub parse_results{
  my ($self, $results) =  @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  open(OUT, "<".$results) or throw("FAILED to open ".$results.
                                   "RepeatMasker:parse_results");
  REPEAT:while(<OUT>){
      if (/no repetitive sequences detected/) {
        print "RepeatMasker didn't find any repetitive sequences";
        last REPEAT; #no repeats found no point carrying on
      }
      my @columns;
      if(/\d+/){ #ignoring introductory lines
        chomp;
        @columns = split;
        pop @columns if $columns[-1] eq '*';
        if (@columns != 15) {
          throw("Can't parse repeatmasker output unexpected number ".
                "of columns in the output ".@columns." in ".$_." ".
                "RepeatMasker:parse_results");
        }
        my ($score, $query_name, $query_start, $query_end, $strand,
            $repeat_name, $repeat_class) 
          =  @columns[0, 4, 5, 6, 8, 9, 10];
        my $start_column;
        if($strand eq '+'){ 
          $start_column = 11;
          $strand = 1;
        }
        if($strand eq 'C'){
          $start_column = 13;
          $strand = -1;
        }#the location of the start and end inside the repeat 
        #is different depending on the strand
        my $repeat_start = $columns[$start_column];
        my $repeat_end = $columns[12];
        if($repeat_end < $repeat_start){
          my $temp_end = $repeat_start;
          $repeat_start = $repeat_end;
          $repeat_end = $temp_end;
        }
        my $rc = $self->_get_consensus($repeat_name, $repeat_class);
        my $rf = Bio::EnsEMBL::RepeatFeature->new;
        $rf->score            ($score);
        $rf->start            ($query_start);
        $rf->end              ($query_end);
        $rf->strand           ($strand);
        $rf->hstart           ($repeat_start);
        $rf->hend             ($repeat_end);
        $rf->repeat_consensus ($rc);
        $rf->slice($self->query);
        $self->output([$rf]);
      }
    }
  close(OUT) or throw("FAILED to close ".$results.
                      "RepeatMasker:parse_results");
}


=head2 _get_consensus

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::RepeatMasker
  Arg [2]   : string, name of repeat
  Arg [3]   : string, name of class
  Function  : create a new Bio::EnsEMBL::RepeatConsensus object
  Returntype: Bio::EnsEMBL::RepeatConsensus
  Exceptions: 
  Example   : 

=cut


sub _get_consensus {
    my ($self, $name, $class) = @_;
    my $cons;

    unless ($cons = $self->{'_consensi'}{$name}) {
      $cons = new Bio::EnsEMBL::RepeatConsensus;
      $cons->name($name);
      $cons->repeat_class($class);
      $self->{'_consensi'}{$name} = $cons;
    }

    return $cons;

}
