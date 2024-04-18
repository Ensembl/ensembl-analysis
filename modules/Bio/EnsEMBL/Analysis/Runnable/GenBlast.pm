# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::GenBlast

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::GenBlast->new(
      -query => $slice,
      -program => 'genblast',
     );
  $runnable->run;
  my @predictions = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run the genBlast gene prediction program that is based on
protein homology (http://genome.sfu.ca/genblast/) against a set of
proteins.  The resulting output file is parsed into prediction
transcripts.

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Runnable::GenBlast;

use strict;
use warnings;

use File::Basename;
use File::Spec::Functions qw(tmpdir);

use Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio);


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::GenBlast
  Function  : create a Bio::EnsEMBL::Analysis::Runnable::GenBlast runnable
  Returntype: Bio::EnsEMBL::Analysis::Runnable::GenBlast
  Exceptions: none
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($database,$ref_slices) = rearrange([qw(DATABASE REFSLICES)], @args);
  $self->database($database) if defined $database;
  $self->genome_slices($ref_slices) if defined $ref_slices;

  throw("You must supply a database") if not $self->database; 
  throw("You must supply a query") if not $self->query;
  throw("You must supply a hash of reference slices") if not $self->genome_slices;

  return $self;
}

############################################################
#
# Analysis methods
#
############################################################

=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : string, directory
  Function  : a generic run method. This checks the directory specifed
  to run it, write the query sequence to file, marks the query sequence
  file and results file for deletion, runs the analysis parses the 
  results and deletes any files
  Returntype: 1
  Exceptions: throws if no query sequence is specified
  Example   : 

=cut


sub run{
  my ($self, $dir) = @_;
  $self->workdir($dir) if($dir);
  throw("Can't run ".$self." without a query sequence") 
    unless($self->query);
  $self->checkdir();
  $self->run_analysis();
  $self->parse_results;
  return 1;
}



=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, program name
  Function  : create and open a commandline for one
  of the ab initio gene finding programs
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

  throw($program." is not executable GenBlast::run_analysis ") 
    unless($program && -x $program);

  my $workdir = tmpdir();

  # set up environment variables
  # we want the path of the program_file
  my $dir = dirname($program);
  $ENV{GBLAST_PATH} = $dir;
  $ENV{path} = "" if (!defined $ENV{path});
  $ENV{path} = "(".$ENV{path}.":".$dir.")";

  # link to alignscore.txt if not already linked
  chdir $workdir;
  my $ln_cmd = "ln -s ".$dir."/alignscore.txt alignscore.txt";
  my $value = system($ln_cmd) unless (-e "alignscore.txt"); 

  # genBlast sticks "_1.1c_2.3_s1_0_16_1" on the end of the output
  # file for some reason - it will probably change in future
  # versions of genBlast.  
  my $outfile_suffix = "_1.1c_2.3_s1_0_16_1";
  my $outfile_glob_prefix = $self->query . $outfile_suffix;

  # if there are old files around, need to get rid of them
  foreach my $oldfile (glob("${outfile_glob_prefix}*")) {
    unlink $oldfile;
  }

  my $command = $program .
  ' -p genblastg '.
  ' -q '.$self->query.
  ' -t '.$self->database.
  ' -o '.$self->query.
  ' -cdna -pro   '.$self->options;

  $self->resultsfile($self->query. $outfile_suffix. ".gff");

  system($command) == 0 or throw("FAILED to run ".$command);

  foreach my $file (glob("${outfile_glob_prefix}*")) {
    $self->files_to_delete($file);
  }
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genscan
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
                                   "GenBlast:parse_results");
  my (%transcripts, @transcripts);

  LINE:while(<OUT>){
    chomp;
    if(/^#/){
      next LINE;
    }
    # ##gff-version 3
    # ##sequence-region       1524908_group1  1       6746
    # 1524908 genBlastG       transcript      72538   75301   81.7114 +       .       ID=2RSSE.1-R1-1-A1;Name=2RSSE.1
    # 1524908 genBlastG       coding_exon     72538   72623   .       +       .       ID=2RSSE.1-R1-1-A1-E1;Parent=2RSSE.1-R1-1-A1
    # 1524908 genBlastG       coding_exon     73276   73336   .       +       .       ID=2RSSE.1-R1-1-A1-E2;Parent=2RSSE.1-R1-1-A1
    # 1524908 genBlastG       coding_exon     73694   73855   .       +       .       ID=2RSSE.1-R1-1-A1-E3;Parent=2RSSE.1-R1-1-A1
    # 1524908 genBlastG       coding_exon     74260   74372   .       +       .       ID=2RSSE.1-R1-1-A1-E4;Parent=2RSSE.1-R1-1-A1
    # 1524908 genBlastG       coding_exon     74629   74800   .       +       .       ID=2RSSE.1-R1-1-A1-E5;Parent=2RSSE.1-R1-1-A1
    # 1524908 genBlastG       coding_exon     75092   75301   .       +       .       ID=2RSSE.1-R1-1-A1-E6;Parent=2RSSE.1-R1-1-A1

    if(/transcript|coding_exon/i){
      my @elements = split;
      if(@elements != 9){
        throw("Can't parse ".$_." splits into wrong number of elements ".
              "GenBlast:parse_results");
      }
      my ($chromosome, $type, $start, $end, $score, $strand, $other) =  @elements[0, 2, 3, 4, 5, 6, 8];


      if ($type eq 'transcript') {
        my ($group, $hitname) = ($other =~ /ID=(\S+?);Name=(\S+)/);
        $transcripts{$group}->{score} = $score;
        $transcripts{$group}->{hitname} = $hitname;
      } elsif ($type eq 'coding_exon') {
        my ($group) = ($other =~ /Parent=(\S+)/);

        if (not exists $self->genome_slices->{$chromosome}) {
          throw("No slice supplied to runnable with for $chromosome");
        }
          

        my $prediction_exon = Bio::EnsEMBL::PredictionExon->new(-start => $start,
                                                                -end   => $end,
                                                                -strand => $strand eq '-' ? -1 : 1,
                                                                -analysis => $self->analysis,
                                                                -score => $score,
                                                                -slice => $self->genome_slices->{$chromosome});
        push @{$transcripts{$group}->{exons}}, $prediction_exon;
      }
      
    }
  }
  close(OUT) or throw("FAILED to close ".$results.
                      "GenBlast:parse_results");
  
  foreach my $tid (keys %transcripts) {
    my @exons = sort { $a->start <=> $b->start } @{$transcripts{$tid}->{exons}};
    map { $_->score($transcripts{$tid}->{score}) } @exons;
    $self->set_phases(0, \@exons);

    my $tran = Bio::EnsEMBL::PredictionTranscript->new(-exons => \@exons,
                                                       -slice => $exons[0]->slice,
                                                       -analysis => $self->analysis,
                                                       -display_label => $transcripts{$tid}->{hitname});

    push @transcripts, $tran;

    my $pep = $tran->translate->seq;
    if ($pep =~ /\*/) {
      printf STDERR "Bad translation for $tid : $pep\n";
    }
  }
  
  $self->clean_output;
  $self->output(\@transcripts);

}






############################################################
#
# get/set methods
#
############################################################

=head2 query

    Title   :   query
    Usage   :   $self->query($seq)
    Function:   Get/set method for query.  If set with a Bio::Seq object it
                will get written to the local tmp directory
    Returns :   filename
    Args    :   Bio::PrimarySeqI, or filename

=cut

sub query {
  my ($self, $val) = @_;
  
  if (defined $val) {
    if (not ref($val)) {   
      throw("[$val] : file does not exist\n") unless -e $val;
    } elsif (not $val->isa("Bio::PrimarySeqI")) {
      throw("[$val] is neither a Bio::Seq not a file");
    }
    $self->{_query} = $val;
  }
  
  return $self->{_query}
}

=head2 database
  
    Title   :   database
    Usage   :   $self->database($seq)
    Function:   Get/set method for database.  If set with a Bio::Seq object it
                will get written to the local tmp directory
    Returns :   filename
    Args    :   Bio::PrimarySeqI, or filename

=cut

sub database {
  my ($self, $val) = @_;
  
  if (defined $val) {
    if (not ref($val)) {   
      throw("[$val] : file does not exist\n") unless -e $val;
    } else {
      if (ref($val) eq 'ARRAY') {
        foreach my $el (@$val) {
          throw("All elements of given database array should be Bio::PrimarySeqs")
        }
      } elsif (not $val->isa("Bio::PrimarySeq")) {
        throw("[$val] is neither a file nor array of Bio::Seq");
      } else {
        $val = [$val];
      }
    }
    $self->{_database} = $val;
  }
  
  return $self->{_database};
}


sub genome_slices {
  my ($self, $val) = @_;
  
  if (defined $val) {
    $self->{_genome_slices} = $val;
  }

  return $self->{_genome_slices}
}


1;
