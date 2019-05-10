=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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

package Bio::EnsEMBL::Analysis::Runnable::GenBlastGene;

use strict;
use warnings;
use feature 'say';


use File::Spec::Functions qw(splitpath catfile);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_timer);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent('Bio::EnsEMBL::Analysis::Runnable');

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

  my ($database,$genblast_program,$max_rank,$genblast_pid,$database_adaptor) = rearrange([qw(DATABASE
                                                                                             GENBLAST_PROGRAM
                                                                                             MAX_RANK
                                                                                             GENBLAST_PID
                                                                                             DATABASE_ADAPTOR
                                                                                         )], @args);
  $self->database($database) if defined $database;
  # Allows the specification of exonerate or genewise instead of genblastg. Will default to genblastg if undef
  $self->genblast_program($genblast_program || 'genblastg');
  $self->max_rank($max_rank) if defined $max_rank;
  $self->genblast_pid($genblast_pid) if defined $genblast_pid;
  $self->database_adaptor($database_adaptor) if defined $database_adaptor;
  throw("You must supply a database") if not $self->database;
  throw("You must supply a database adaptor to fetch the genomic sequence from for supporting features") if not $self->database_adaptor;

  # Default max rank
  unless($self->max_rank()) {
    $self->max_rank(5);
    warning("No max_rank parameter provided so defaulting to: ".$self->max_rank());
  }

  # Default pid
  unless($self->genblast_pid()) {
    $self->genblast_pid(50);
    warning("No genblast_pid parameter provided so defaulting to: ".$self->genblast_pid());
  }

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
  $self->checkdir();
  $self->write_seq_file if ($self->query and ref($self->query) eq 'ARRAY');
  throw("Can't run ".$self." without a query sequence")
    unless($self->queryfile);
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

  # link to alignscore.txt if not already linked
  my (undef, $directory, $progname) = splitpath($program);
  if ($directory) {
    chdir $directory;
  }

  # if there are old files around, need to get rid of them
  my $genblast_program = $self->genblast_program;
  unless($genblast_program) {
    $genblast_program = "genblastg";
  }

  my $outfile = $self->resultsfile;
  my $command = $program .
  ' -p '.$genblast_program.
  ' -q '.$self->queryfile.
  ' -t '.$self->database.
  ' -o '.$outfile.
  ' -g T -pid -r '.$self->max_rank.' '.$self->options;

  execute_with_timer($command, $self->timer);

  my @files = glob qq(${outfile}*.gff);
  if (@files == 1) {
    $self->resultsfile($files[0]);
  }
  else {
    throw("More than one gff file has been created: ".join(',', @files));
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

  my $dba = $self->database_adaptor();
  my $slice_adaptor = $dba->get_SliceAdaptor();

  open(OUT, "<".$results) or throw("FAILED to open ".$results."\nGenBlast:parse_results");
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
      my $slice = $slice_adaptor->fetch_by_name($chromosome);

      if ($type eq 'transcript') {
        #ID=Q502Q5.1-R1-1-A1;Name=Q502Q5.1;PID=67.05;Coverage=99.36;Note=PID:67.05-Cover:99.36
        my ($group, $hitname, $pid, $cov) = ($other =~ /ID=(\S+?);Name=([^;]+);PID=([^;]+);Coverage=([^;]+);/);
        $group =~ /^$hitname\-R(\d+)\-/;
        my $rank = $1;
        $transcripts{$group}->{score} = $score;
        $transcripts{$group}->{hitname} = $hitname;
        $transcripts{$group}->{pid} = $pid;
        $transcripts{$group}->{cov} = $cov;
        $transcripts{$group}->{rank} = $rank;
        $transcripts{$group}->{slice} = $slice;
      } elsif ($type eq 'coding_exon') {
        my ($group) = ($other =~ /Parent=(\S+)/);
#        if (not exists $self->genome_slices->{$chromosome}) {
#          throw("No slice supplied to runnable with for $chromosome");
#        }

        my $exon = Bio::EnsEMBL::Exon->new(-start => $start,
                                           -end   => $end,
                                           -strand => $strand eq '-' ? -1 : 1,
                                           -analysis => $self->analysis,
                                           -slice => $transcripts{$group}->{slice});

        push @{$transcripts{$group}->{exons}}, $exon;
      }
    }
  }
  close(OUT) or throw("FAILED to close ".$results.
                      "GenBlast:parse_results");

  foreach my $tid (keys %transcripts) {

    # Skip transcripts that fail the pid cut-off
    unless($transcripts{$tid}->{pid} >= $self->genblast_pid()) {
      warning("Skipping transcript because of percent identity cut-off:\n".
              "Name: ".$tid."\nPercent identity: ".$transcripts{$tid}->{pid}."\nCut-off: ".$self->genblast_pid());
      next;
    }

    my @exons = sort { $a->start <=> $b->start } @{$transcripts{$tid}->{exons}};


    my $tran = Bio::EnsEMBL::Transcript->new(-exons => \@exons,
                                             -analysis => $self->analysis,
                                             -stable_id => $transcripts{$tid}->{hitname},
                                             -slice => $transcripts{$tid}->{slice});

    $tran->{'accession'} = $transcripts{$tid}->{hitname};
    $tran->{'pid'} = $transcripts{$tid}->{pid};
    $tran->{'cov'} = $transcripts{$tid}->{cov};
    $tran->{'rank'} = $transcripts{$tid}->{rank};
    $tran->{'genblast_id'} = $tid;
    $tran->{'slice_name'} = $transcripts{$tid}->{slice_name};

    # Reverse the exons for negative strand to calc the translations
    my $strand = $exons[0]->strand;
    if($strand == -1) {
      @exons = sort { $b->start <=> $a->start } @{$transcripts{$tid}->{exons}};
    }

    my $start_exon = $exons[0];
    my $end_exon = $exons[scalar(@exons)-1];
    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($start_exon);
    $translation->start(1);
    $translation->end_Exon($end_exon);
    $translation->end($end_exon->length());
    $tran->translation($translation);

    # Set the phases
    calculate_exon_phases($tran, 0);

    # Set the exon and transcript supporting features
    $self->set_supporting_features($tran);

    push @transcripts, $tran;

#    my $pep = $tran->translate->seq;
#    if ($pep =~ /\*/) {
#      printf STDERR "Bad translation for $tid : $pep\n";
#    }
  }

  $self->clean_output;
  $self->output(\@transcripts);

}

sub set_supporting_features {
  my ($self,$transcript) = @_;

  # This sets supporting features for all exons. It does this by extracting the final alignment from the report file
  # and comparing all non-split codons in the exons to their corresponding positions in the alignment.
  # For each non-split codon, it makes a feature pair where the start and end are the codon start and end and the
  # hit start and end is the position of the corresponding query sequence amino acid in the UNGAPPED version of the
  # query sequence. So it 1) Finds the codon in the alignment 2) Finds the index of the corresponding query amino acid in the
  # unaligned query
  # Once codons that have support have had the corresponding hit start and end worked out, it then joins then together
  # It does this in a simple way, as the proto feature pairs are already in codon by codon order it just looks at consecutive
  # ones and checks that they are from adjacent codons. If they are it then checks if the hit end of the first is directly
  # beside the hit start of the second. If they are then the left proto-feature pair is deleted and the right one is extended
  # to encompass the left one. This repeats for as long as possible to build maximum lenght feature pairs
  # In the cases where there is a gap in the query, there will be no proto SF for the corresponding codon
  # In the case where there is a gap in the target (genome), the hit start and hit ends will not be side by side and thus not
  # joined together
  # By following this an exon might have several feature pairs, some might be side by side and unjoined, some might be on either
  # side of codons with no feature pairs. Often there is just one long feature pair covering all non-split codons in the exon
  # All of the feature pairs on an exon are then passing into a single object (DnaPepAlignFeature). By doing this the API will
  # automatically work out the cigar strings for each exon. By looking at gaps between the feature pairs and at feature pairs
  # that are adjacent but the hit end/start are not adjacent it will model the indels across the exon
  # For the transcript supporting feature you need to add every individual feature pair into an array ref and add to a DnaPepAlignFeature

  # This will store all feature pairs in the end so that they can be added as a transcript supporting feature
  my $all_exon_supporting_features = [];

  # This is the id from the gff file that is unique to each gene. Use this to find the alignment in the report file
  my $genblast_id = $transcript->{'genblast_id'};

  my $report_file = $self->resultsfile;
  $report_file =~ s/\.gff$//;

  unless(-e $report_file) {
    throw("Tried to find the report file but couldn't. Path checked:\n".$report_file);
  }


  my @report_array = ();
  my $query_seq;
  my $target_seq;
  my $gene_info;
#  my $transcript_percent_id = $transcript->{'pid'};
#  my $transcript_coverage = $transcript->{'cov'};
  my $transcript_rank = $transcript->{'rank'};

  my $found = 0;
  open(GENBLAST_REPORT,$report_file);
  while(<GENBLAST_REPORT>) {
    my $line = $_;
    chomp $line;

    push(@report_array,$line);
    if(scalar(@report_array) > 5) {
      shift(@report_array);
    }

    if($line =~ /^Gene\:ID\=$genblast_id\|/) {
      $found = 1;
      $query_seq = $report_array[0];
      $query_seq =~ s/^query\://;
      $target_seq = $report_array[2];
      $target_seq =~ s/^targt\://;
      last;
    }
  }
  close(GENBLAST_REPORT);

  unless($query_seq && $target_seq) {
    throw("Did not retrieve query and target sequence from the report file");
  }

  my $coverage;
  my $percent_id;
  ($query_seq,$target_seq,$coverage,$percent_id) = $self->realign_translation($query_seq,$target_seq);

  # Add a stop to the alignment seqs. Basically this will allow me to ignore the final codon (which is a stop)
  $query_seq .= '*';
  $target_seq .= '*';

  say "";
  say "------------------------------------------------";
  say "NEW TRANSCRIPT";
  say "------------------------------------------------";

  my $codon_index= 0;
  my $exons = $transcript->get_all_Exons();
  my $i=0;
  for($i=0; $i<scalar(@{$exons}); $i++) {
    my $proto_supporting_features = [];
    my $exon = $$exons[$i];
    my $exon_seq = $exon->seq->seq();
    my @nucleotide_array = split('',$exon_seq);
    my $phase = $exon->phase();
    my $end_phase = $exon->end_phase();
    my $start_index = 0;

    say "";
    say "------------------------------------------------";
    say "NEW EXON";
    say "------------------------------------------------";
    say "FM2 ESTART: ".$exon->start;
    say "FM2 EEND: ".$exon->end;
    say "FM2 ESTRAND: ".$exon->strand;
    say "FM2 EPHASE: ".$exon->phase;
    say "FM2 EENDPHASE: ".$exon->end_phase;

    # If the phase is not 0 then the first codon is a split one. The phase is then number of bases missing from
    # the codon. so if you take the phase from 3 you get the number of bases in the split codon
    if($phase) {
      $start_index += (3 - $phase);
    }

    my @unjoined_feature_pairs = ();
    for(my $k=$start_index; $k<scalar(@nucleotide_array); $k+=3) {

     # Ending on a split codon, so increase the index and finish the loop
     if($k+2 >= scalar(@nucleotide_array)) {
        $codon_index++;
        say "FM2 SKIPPING SPLIT CODON";
        last;
      }

      my $target_start;
      my $target_end;
      if($exon->strand == 1) {
         $target_start = $exon->start + $k;
         $target_end = $target_start + 2;
      } else {
        $target_start = $exon->end - $k - 2;
        $target_end = $target_start + 2;
      }

      # Now need to look at the alignment index for this codon. If there is a gap in the query sequence at that index
      # then the codon should be skipped. If it isn't a gap at that codon position then
      my $codon_alignment_index = $self->find_codon_alignment_index($codon_index,$target_seq);

      my $query_char = substr($query_seq,$codon_alignment_index,1);
      my $target_char = substr($target_seq,$codon_alignment_index,1);
      if($query_char eq '-') {
        say "FM2 TSTART: ".$target_start;
        say "FM2 TEND: ".$target_end;
        say "FM2 CODON INDEX: ".$codon_index;
        say "FM2 CODON ALIGNMENT INDEX: ".$codon_alignment_index;
        say "FM2 CODON CHARS: '".$nucleotide_array[$k].$nucleotide_array[$k+1].$nucleotide_array[$k+2]."'";
        say "FM2 ALIGNMENT CHARS: '".$query_char."'='".$target_char."'";
        say "FM2 SKIPPING CODON BECAUSE OF ALIGMENT GAP";
        $codon_index++;
        next;
      } elsif(($codon_alignment_index == length($query_seq)-1) && $query_char eq '*' && $target_char eq '*') {
        say "FM2 LAST CODON IS STOP SO SKIPPING";
        next;
      }

      my $hit_start = $self->find_hit_start($codon_alignment_index,$query_seq);
      my $hit_end = $hit_start;

      say "FM2 TSTART: ".$target_start;
      say "FM2 TEND: ".$target_end;
      say "FM2 HSTART: ".$hit_start;
      say "FM2 HEND: ".$hit_end;
      say "FM2 CODON INDEX: ".$codon_index;
      say "FM2 CODON ALIGNMENT INDEX: ".$codon_alignment_index;
      say "FM2 CODON CHARS: '".$nucleotide_array[$k].$nucleotide_array[$k+1].$nucleotide_array[$k+2]."'";
      say "FM2 ALIGNMENT CHARS: '".$query_char."'---'".$target_char."'";
      $codon_index++;

      push(@{$proto_supporting_features},{'tstart' => $target_start,
                                          'tend'   => $target_end,
                                          'hstart' => $hit_start,
                                          'hend'   =>$hit_end});
    }

    say "\nQUERY:\n".$query_seq;
    say "\nTARGET:\n".$target_seq;

   my $joined_supporting_features = $self->join_supporting_features($proto_supporting_features,$exon->strand);
   my $exon_feature_pairs = [];
   foreach my $joined_feature (@{$joined_supporting_features}) {
     my $feature_pair = Bio::EnsEMBL::FeaturePair->new(
                                                        -start      => $joined_feature->{'tstart'},
                                                        -end        => $joined_feature->{'tend'},
                                                        -strand     => $exon->strand,
                                                        -hseqname   => $transcript->{'accession'},
                                                        -hstart     => $joined_feature->{'hstart'},
                                                        -hend       => $joined_feature->{'hend'},
                                                        -hcoverage  => $coverage,
                                                        -percent_id => $percent_id,
                                                        -slice      => $exon->slice,
                                                        -analysis   => $transcript->analysis);
     say "FM2 ADD SUPPORTING EVIDENCE START: ".$feature_pair->start;
     say "FM2 ADD SUPPORTING EVIDENCE END: ".$feature_pair->end;
     push(@{$exon_feature_pairs},$feature_pair);
     push(@{$all_exon_supporting_features},$feature_pair);
   }

    if(scalar(@{$exon_feature_pairs})) {
      my $final_exon_supporting_features = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => $exon_feature_pairs, -align_type => 'ensembl');
       $exon->add_supporting_features($final_exon_supporting_features);
    } else {
      warning("No supporting features added for exon.\nExon start: ".$exon->start."\nExon end: ".$exon->end);
    }

  }

  my $transcript_supporting_features = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => $all_exon_supporting_features, -align_type => 'ensembl');
  $transcript->add_supporting_features($transcript_supporting_features);

}

sub find_codon_alignment_index {
  my ($self,$codon_index,$align_seq) = @_;

  my $align_index = -1;
  my $char_count = 0;
  for(my $j=0; $j<length($align_seq); $j++) {
    my $char = substr($align_seq,$j,1);

    if($char eq '-') {
      next;
    }

    if($char_count == $codon_index) {
      $align_index = $j;
      last;
    }

    $char_count++;

  }

  # For the last codon
  if($align_index == -1 && $char_count == $codon_index) {
    $align_index = length($align_seq) - 1;
  }

  unless($align_index >= 0) {
    throw("Did not find the alignment index for the codon");
  }

  return($align_index);
}

sub find_hit_start {
  my ($self,$alignment_index,$align_seq) = @_;

  my $sub_seq = substr($align_seq,0,$alignment_index + 1);
  $sub_seq =~ s/\-//g;

  my $hit_start = length($sub_seq);

  return($hit_start);
}

sub join_supporting_features {
  my ($self,$proto_supporting_features,$strand) = @_;

  # At this point all supporting features for the exon have been built on a per codon basis. Before making
  # a feature pair we want to combine adjacent supporting features. If both the genomic end of the left
  # feature matches the genomic start - 1 of the next then they can be joined if the hit end of the first
  # is hit start - 1 of the next

  my $joined_supporting_features = [];
  if($strand == 1) {
    for(my $i=0; $i<scalar(@{$proto_supporting_features})-1; $i++) {
      my $left_proto = $$proto_supporting_features[$i];
      my $right_proto = $$proto_supporting_features[$i+1];

      if($left_proto->{'tend'} == ($right_proto->{'tstart'} - 1)) {
        if($left_proto->{'hend'} == ($right_proto->{'hstart'} - 1)) {
          # If this is the case then the codons and hits are contiguous and so they can be joined
          $right_proto->{'tstart'} = $left_proto->{'tstart'};
          $right_proto->{'hstart'} = $left_proto->{'hstart'};
          $$proto_supporting_features[$i] = 0;
          $$proto_supporting_features[$i+1] = $right_proto;
        }
      }
    }
  } else {
    for(my $i=scalar(@{$proto_supporting_features})-1; $i>0; $i--) {
      my $left_proto = $$proto_supporting_features[$i];
      my $right_proto = $$proto_supporting_features[$i-1];

      say "FM2 PROTO LTS: ".$left_proto->{'tstart'};
      say "FM2 PROTO LTE: ".$left_proto->{'tend'}; 
      say "FM2 PROTO LHS: ".$left_proto->{'hstart'};
      say "FM2 PROTO LHE: ".$left_proto->{'hend'};
      say "FM2 PROTO RTS: ".$right_proto->{'tstart'};
      say "FM2 PROTO RTE: ".$right_proto->{'tend'};
      say "FM2 PROTO RHS: ".$right_proto->{'hstart'};
      say "FM2 PROTO RHE: ".$right_proto->{'hend'};

      if($left_proto->{'tend'} == ($right_proto->{'tstart'} - 1)) {
        if($left_proto->{'hend'} == ($right_proto->{'hstart'} + 1)) {
          # If this is the case then the codons and hits are contiguous and so they can be joined
          $right_proto->{'tstart'} = $left_proto->{'tstart'};
          $right_proto->{'hstart'} = $left_proto->{'hstart'};
          $$proto_supporting_features[$i] = 0;
          $$proto_supporting_features[$i-1] = $right_proto;
        }
      }
    }
  }
  foreach my $proto_sf (@{$proto_supporting_features}) {
    if($proto_sf) {
      say "FM2 JOINED TSTART: ".$proto_sf->{'tstart'};
      say "FM2 JOINED TEND: ".$proto_sf->{'tend'};
      say "FM2 JOINED HSTART: ".$proto_sf->{'hstart'};
      say "FM2 JOINED HEND: ".$proto_sf->{'hend'};

      # If it's the negative strand then swap the start and end of the hit
      if($strand == -1) {
        my $temp = $proto_sf->{'hstart'};
        $proto_sf->{'hstart'} = $proto_sf->{'hend'};
        $proto_sf->{'hend'} = $temp;
      }
      push(@{$joined_supporting_features},$proto_sf);
    }
  }

  return($joined_supporting_features);

}

sub project_to_alignment {
  my ($self,$codon_index,$query_seq,$target_seq) = @_;

  my $ungapped_index=0;
  for(my $i=0; $i<length($query_seq); $i++) {
    my $query_char = substr($query_seq,$i,1);
    if($query_char ne "-") {
      if($ungapped_index == $codon_index) {
        return($i);
      }
      $ungapped_index++;
    }
  }
}


sub realign_translation {
  my ($self,$query_seq,$target_seq) = @_;

  $query_seq =~ s/\-//g;
  $target_seq =~ s/\-//g;

  my $align_input_file = $self->resultsfile;
  $align_input_file =~ s/\.gff$/\.prealn/;

  my $align_output_file = $self->resultsfile;
  $align_output_file =~ s/\.gff$/\.aln/;

  open(INPUT,">".$align_input_file);
  say INPUT ">query";
  say INPUT $query_seq;
  say INPUT ">target";
  say INPUT $target_seq;
  close INPUT;

  my $align_program_path = 'muscle';

  my $cmd = $align_program_path." -in ".$align_input_file." -out ".$align_output_file;
  my $result = system($cmd);

  if($result) {
    throw("Got a non-zero exit code from alignment. Commandline used:\n".$cmd);
  }

  my $file = "";
  open(ALIGN,$align_output_file);
  while(<ALIGN>) {
    $file .= $_;
  }
  close ALIGN;

  $file =~ /\>.+\n(([^>]+\n)+)\>.+\n(([^>]+\n)+)/;
  my $aligned_query_seq = $1;
  my $aligned_target_seq = $3;

  $aligned_query_seq =~ s/\n//g;
  $aligned_target_seq =~ s/\n//g;

  `rm $align_input_file`;
  `rm $align_output_file`;

  # Work out coverage
  my $coverage;
  my $temp = $aligned_target_seq;
  my $target_gap_count = $temp =~ s/\-//g;
  my $ungapped_query_seq = $aligned_query_seq;
  $ungapped_query_seq  =~ s/\-//g;

  if(length($ungapped_query_seq) == 0) {
    $coverage = 0;
  } else {
    $coverage = 100 - (($target_gap_count/length($ungapped_query_seq)) * 100);
  }

  # Work out precent identity
  my $match_count = 0;
  my $aligned_positions = 0;
  for(my $j=0; $j<length($aligned_query_seq); $j++) {
    my $char_query = substr($aligned_query_seq,$j,1);
    my $char_target = substr($aligned_target_seq,$j,1);
    if($char_query eq '-' || $char_target  eq '-') {
      next;
    }
    if($char_query eq $char_target) {
      $match_count++;
    }
    $aligned_positions++;
  }

  unless($aligned_positions) {
    throw("Pairwise alignment between the query sequence and the translation shows zero aligned positions. Something has gone wrong");
  }

  my $percent_id = ($match_count / $aligned_positions) * 100;

  return($aligned_query_seq,$aligned_target_seq,$coverage,$percent_id);

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

  return $self->{_genome_slices};
}

sub genblast_program {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_genblast_program} = $val;
  }

  return $self->{_genblast_program};
}

sub max_rank {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_max_rank} = $val;
  }

  return $self->{_max_rank};
}

sub genblast_pid {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_genblast_pid} = $val;
  }

  return $self->{_genblast_pid};
}

sub database_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_database_adaptor} = $val;
  }

  return $self->{_database_adaptor};
}


1;
