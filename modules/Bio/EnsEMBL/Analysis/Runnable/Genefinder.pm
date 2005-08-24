package Bio::EnsEMBL::Analysis::Runnable::Genefinder;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($tablenamefile, $intronpenalty, $exonpenalty) = 
    rearrange(['TABLENAMEFILE', 
               'INTRONPENALTY', 
               'EXONPENALTY'], @args);

  
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->tablenamefile('/usr/local/ensembl/data/nemtables/sanger.tablenamefile.cod');
  $self->intronpenalty('/usr/local/ensembl/data/nemtables/intron_penalty.lookup');
  $self->exonpenalty('/usr/local/ensembl/data/nemtables/exon_penalty.lookup');
  $self->program('genefinder') if(!$self->program);
  ######################
  
  $self->tablenamefile($tablenamefile);
  $self->intronpenalty($intronpenalty);
  $self->exonpenalty($exonpenalty);

  return $self;
}


sub run_analysis{
  my ($self, $program) = @_;

  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable Genefinder::run_analysis ") 
    unless($program && -x $program);

  my $err_string = $self->filecheck;
  throw("Filecheck found missing files $err_string") if($err_string);

  my $command = $program .' -tablenamefile '.$self->tablenamefile.
    ' -intronPenaltyFile '.$self->intronpenalty.' -exonPenaltyFile '.
      $self->exonpenalty.' '.$self->queryfile.' > '.
        $self->resultsfile;
  print "Running analysis ".$command."\n";
  my $value = system($command); #or throw("FAILED to run ".$command);
  throw("FAILED to run ".$command) if(!-s $self->resultsfile);
}



sub filecheck{
    my ($self) = @_;
    my $err_string;
    if(!-e $self->tablenamefile){
      $err_string = "Tablenamefile ".$self->tablenamefile." does not ".
        "exist\n";
    }
    if(!-e $self->intronpenalty){
       $err_string .= "IntronPenalty ".$self->intronpenalty." does not".
         " exists\n";
    }
    if(!-e $self->exonpenalty){
      $err_string .= "ExonPenalty ".$self->exonpenalty." does not ".
        "exist\n";
    }

    return $err_string;
}

sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  
  open(OUT, "<".$results) or throw("FAILED to open ".$results.
                                   "Genefinder:parse_results");
  my $ff = $self->feature_factory;
  my $prefix;
  my @forward_lines;
  my @reverse_lines;
  while(<OUT>){ # this loop sorts the lines in to foward 
    # and reverse genes
    chomp;
    if (/Sequence: (\S+)/) {
      $prefix = $1; 
    }
    if($_ =~ /\.\./){
      if($_ =~ /\*/){
        push(@reverse_lines, $_);
      }else{
        push(@forward_lines, $_);
      }
    }
  }
  
  my $gene_count = $self->parse_lines(\@forward_lines, $prefix, $ff, 
                                   0, 1);

  
  $gene_count = $self->parse_lines(\@reverse_lines, $prefix, $ff, 
                                   $gene_count, -1);
 
  close(OUT) or throw("FAILED to close ".$results.
                      "Genefinder:parse_results");

  $self->create_transcripts;
  $self->calculate_phases;
}

#accessor methods


sub calculate_phases{
  my ($self) = @_;

  my @phases = (0, 1, 2);
  my @transcripts = @{$self->output};
  my @checked_transcripts;
  my $ff = $self->feature_factory;
  $self->clean_output;
 TRANS:foreach my $transcript(@transcripts){
  PHASE:foreach my $phase(@phases){
      my $temp_exons = $self->set_phases($phase,
                                         $transcript->get_all_Exons,
                                        );
      my $new = $ff->create_prediction_transcript($temp_exons, 
                                                  $self->query,
                                                  $self->analysis);
      if($new->translate->seq !~ /\*/){
        push(@checked_transcripts, $new);
        next TRANS;
      }
    }
  }
  if(scalar(@checked_transcripts) != scalar(@transcripts)){
    throw("One of the transcripts doesn't translate as there are ".
          @transcripts." from genefinder but only ".
          @checked_transcripts." transcripts which translate");
  }
  $self->output(\@checked_transcripts);
}

sub parse_lines{
  my ($self, $lines, $prefix, $ff, $gene_count, $strand) = @_;

  my $phase = 0;
  foreach my $line(@$lines){
    print $line."\n";
    $line =~ s/\[TSL: \S+\s+\S+\]//;
    $line =~ s/start//;
    $line =~ s/end//; #end is stripped off
    $line =~ s/\[ -?(\d+\.\d+) \]//;
    my $score = $1;
    if($line =~ s/\*// && $strand != -1){
      throw("Transcript is reverse strand but it isn't ".
            "being set to that");
    }
    $line =~ s/\*//;
    if($strand == -1){
      my @elements = split /\s+/, $line;
      $line = '';
      foreach my $v (reverse(@elements)){
        $line .= $v." ";
      } 
    }
    my @values = split /\s+/, $line;
 
    if(!$values[0]){
      my $first = shift @values;  
    }
    if($values[0] =~ /U(\d)/){

      $phase = $1;
      my $phase_variable = shift @values; 

    }

    my $count = 0;
    my $exonname = $prefix.".".$gene_count;

    foreach my $coord(@values){

      if($coord =~ /U\d/){
        next;
      }

      my ($start, $end) = split /\.\./, $coord;

      my $end_phase = ($phase + ($end-$start) + 1)%3;
      my $exon = $ff->create_prediction_exon($start, $end, $strand, 
                                             $score, 0, $phase, 
                                             $exonname, $self->query,
                                             $self->analysis);
      $self->exon_groups($exonname, $exon);
      $phase =  $end_phase;
      my $count++;
      if($values[$count] && $values[$count] =~ /U(\d)/){
        if($phase != $1){
          warning(" phase ".$phase." and continuation phase ".$1." aren't the same may be issues with translation\n");
        }
      }
    }
    $gene_count++;
  }
  return $gene_count;
}


sub tablenamefile{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'tablenamefile'} = $arg;
  }
  return $self->{'tablenamefile'};
}


sub exonpenalty{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'exonpenalty'} = $arg;
  }
  return $self->{'exonpenalty'};
}

sub intronpenalty{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'intronpenalty'} = $arg;
  }
  return $self->{'intronpenalty'};
}


1;
