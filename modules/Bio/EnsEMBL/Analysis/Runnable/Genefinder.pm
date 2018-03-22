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

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Genefinder - 

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::Genefinder->new(
      -query => $slice,
      -program => 'genefinder',
     );
  $runnable->run;
  my @predictions = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run the genefinder gene predictor and then parse the results
into prediction transcripts

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::Genefinder;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Arg [2]   : string, tablename file
  Arg [3]   : string, intron penalties file
  Arg [4]   : string, exon penalties file
  Function  : create a Genefinder runnable
  Returntype: Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Exceptions: 
  Example   : 

=cut



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



=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Arg [2]   : string, program name
  Function  : create and open a commandline for the program genefinder
  Returntype: none
  Exceptions: throws if the program in not executable or the system
  command fails to execute or the appropriate files are not found
  Example   : 

=cut

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



=head2 filecheck

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Function  : make sure the config files exist
  Returntype: string, any errors found
  Exceptions: 
  Example   : 

=cut



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


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
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

# accessor methods



=head2 calculate_phases

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Function  : works out which phase to make the exons to get 
  a complete cds
  Returntype: none
  Exceptions: throws if the number of transcripts on the way out isnt
  the same as the number on the way in
  Example   : 

=cut


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



=head2 parse_lines

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Arg [2]   : array, of lines from file
  Arg [3]   : string
  Arg [4]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [5]   : int, starting count of genes
  Arg [6]   : int, 1 or -1 for strand
  Function  : parse lines into PredictionExon objects
  Returntype: int, new gene count (incremented for each line parsed)
  Exceptions: throws if * is found and strand isnt -1 and gives a
  warning about phase inconsistencies
  Example   : 

=cut


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




=head2 tablenamefile

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Arg [2]   : string
  Function  : these are all acessor methods for config filenames
  Returntype: string
  Exceptions: 
  Example   : 

=cut



sub tablenamefile{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'tablenamefile'} = $arg;
  }
  return $self->{'tablenamefile'};
}


=head2 exonpenalty

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Arg [2]   : string
  Function  : these are all acessor methods for config filenames
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub exonpenalty{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'exonpenalty'} = $arg;
  }
  return $self->{'exonpenalty'};
}

=head2 intronpenalty

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genefinder
  Arg [2]   : string
  Function  : these are all acessor methods for config filenames
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub intronpenalty{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'intronpenalty'} = $arg;
  }
  return $self->{'intronpenalty'};
}


1;
