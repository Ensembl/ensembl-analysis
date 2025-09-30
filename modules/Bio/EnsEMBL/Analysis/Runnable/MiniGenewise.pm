=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::MiniGenewise - 

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::Runnable::MiniGenewise->new(
      -genomic  => $genseq,
      -features => $features)

    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION


=head1 METHODS

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Runnable::MiniGenewise;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::Runnable::Genewise;
use Bio::EnsEMBL::Analysis::Tools::MiniSeq;
use Bio::EnsEMBL::Analysis::Tools::PairAlign;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(create_Exon);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils 
  qw(create_feature_from_gapped_pieces);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(Transcript_info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

use parent ('Bio::EnsEMBL::Analysis::Runnable');



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my( $protein, $features, $minimum_intron, $terminal_padding, $exon_padding,
      $max_iterate_dist, $genewise_options, $calculate_coverage_and_pid,$best_in_genome)
    = rearrange([qw(PROTEIN FEATURES MINIMUM_INTRON TERMINAL_PADDING EXON_PADDING
                    MAX_SPLIT_ITERATE_DIST GENEWISE_OPTIONS CALCULATE_COVERAGE_AND_PID BEST_IN_GENOME
               )], @args);

  ####SETTING_DEFAULTS###
  $self->minimum_intron(1000);
  $self->exon_padding(200);
  $self->terminal_padding(20000);
  $self->max_split_iterate_distance(0);
  #######################
  $self->protein_sequence($protein);
  $self->features($features);
  $self->minimum_intron($minimum_intron);
  $self->exon_padding($exon_padding);
  $self->terminal_padding($terminal_padding);
  $self->max_split_iterate_distance($max_iterate_dist);
  $self->genewise_options($genewise_options);
  $self->calculate_coverage_and_pid($calculate_coverage_and_pid);
  $self->best_in_genome_transcript($best_in_genome);

  throw("MiniGenewise needs a query sequence") unless($self->query);
  throw("MiniGenewise needs a peptide sequence") unless($self->protein_sequence);
  throw("MiniGenewise needs an array of features")
    unless($self->features && scalar(@{$self->features}));

  return $self;
}



sub best_in_genome_transcript {
   my ($self,$val) = @_;

   if(defined($val) && $val==1) {
     $self->{'_best_in_genome_transcript'} = 1;
   } elsif(defined($val) && $val==0) {
     $self->{'_best_in_genome_transcript'} = 0;
   }

   return($self->{'_best_in_genome_transcript'});
 }

#ACCESSOR METHODS


=head2 genewise_options

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : various things depending on the accessor
  Function  : to store the passed in variable
  Returntype: various
  Exceptions: methods will throw if the wrong type of variable is passed in
  Example   : 

=cut



sub genewise_options {
  my ($self,$arg) = @_;

  $self->{'_genewise_options'} = {} if(!$self->{'_genewise_options'});

  if ($arg) {
    throw("Must pass genewise options a hash ref not a $arg") 
      unless(ref($arg) eq "HASH");
    $self->{'_genewise_options'} = $arg;
  }

  return $self->{'_genewise_options'};
}


=head2 features

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : various things depending on the accessor
  Function  : to store the passed in variable
  Returntype: various
  Exceptions: methods will throw if the wrong type of variable is passed in
  Example   : 

=cut

sub features {
  my ($self,$arg) = @_;

  $self->{features} = [] if(!$self->{features});
  if($arg){
    throw($arg." must be an arrayref of Bio::EnsEMBL::FeaturePair objects")
      unless(ref($arg) eq "ARRAY" && $arg->[0]->isa("Bio::EnsEMBL::FeaturePair"));
    $self->{features} = $arg;
  }
  return $self->{features};
}


=head2 minimum_intron

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : various things depending on the accessor
  Function  : to store the passed in variable
  Returntype: various
  Exceptions: methods will throw if the wrong type of variable is passed in
  Example   : 

=cut
sub minimum_intron {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{'_minimum_intron'} = $arg;
  }

  return $self->{'_minimum_intron'};
}


=head2 exon_padding

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : various things depending on the accessor
  Function  : to store the passed in variable
  Returntype: various
  Exceptions: methods will throw if the wrong type of variable is passed in
  Example   : 

=cut
sub exon_padding {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{'_padding'} = $arg;
  }

  return $self->{'_padding'};
}

=head2 terminal_padding

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : various things depending on the accessor
  Function  : to store the passed in variable
  Returntype: various
  Exceptions: methods will throw if the wrong type of variable is passed in
  Example   : 

=cut
sub terminal_padding {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{'_terminal_padding'} = $arg;
  }
  return $self->{'_terminal_padding'};
}


=head2 max_split_iterate_distance

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : various things depending on the accessor
  Function  : to store the passed in variable
  Returntype: various
  Exceptions: methods will throw if the wrong type of variable is passed in
  Example   : 

=cut
sub max_split_iterate_distance {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{'_iterate'} = $arg;
  }

  return $self->{'_iterate'};
}


=head2 protein_sequence

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : various things depending on the accessor
  Function  : to store the passed in variable
  Returntype: various
  Exceptions: methods will throw if the wrong type of variable is passed in
  Example   : 

=cut
sub protein_sequence {
  my( $self, $value ) = @_;

  if ($value) {
    throw("Input isn't a Bio::PrimarySeq but ".$value) 
      unless($value->isa("Bio::PrimarySeqI"));
    $self->{'_protein_sequence'} = $value;
  }
  return $self->{'_protein_sequence'};
}

=head2 miniseq

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : various things depending on the accessor
  Function  : to store the passed in variable
  Returntype: various
  Exceptions: methods will throw if the wrong type of variable is passed in
  Example   : 

=cut
sub miniseq{
  my ($self, $value) = @_;
  if($value){
    $self->{mini_seq} = $value;
  }
  return $self->{mini_seq};
}

#FUNCTIONALITY


sub run{
  my ($self) = @_;

  my $must_try_again = 1;
  my $features = $self->features;
  my $gw_output;
  while($must_try_again){
    $must_try_again = 0;

    $gw_output = $self->run_Genewise($features);
    if($gw_output eq "FAILED"){
      return;
    }
  
    logger_info("Have ".@$gw_output." output from the genewise run");
    
    if($self->max_split_iterate_distance){
    
      my @bridge_features;
      
      foreach my $g(@$gw_output){

        foreach my $e($g->get_all_Exons){

          my @pieces = @{$self->miniseq->convert_SeqFeature($e)};
        
          if (@pieces > 1) {

            @pieces = sort {$a->start <=> $b->start} @pieces;

            my $can_merge = 1;
            for(my $i=1; $i < @pieces; $i++) {
              my $dist = $pieces[$i]->start - $pieces[$i-1]->end - 1;
              if ($dist > $self->max_split_iterate_distance) {
                $can_merge = 0;
                last;
              }
            }

            if ($can_merge) {
              push @bridge_features, Bio::EnsEMBL::Feature
                ->new(-start => $pieces[0]->start,
                      -end   => $pieces[-1]->end);
            }
          }
        }
      }
      if(@bridge_features){
        push(@$features, @bridge_features);
        $must_try_again = 1;
        warning($self->protein_sequence->id." has a predicted feature which splits ".
                " when mapped back to the genomic, trying again");
      }
    }
  }

  my $output = [];
  foreach my $gene(@$gw_output){
    my $new_gene = $self->convert_to_genomic($gene);
    warning("Failed to convert ".$gene." to genomic coordinates, skipping") 
      unless($new_gene);

    if($self->calculate_coverage_and_pid && $new_gene) {
      $self->realign_translation($new_gene);
    }

    if(defined($self->best_in_genome_transcript()) && $self->best_in_genome_transcript() == 0 && $new_gene) {
      my $analysis = $self->analysis();
      $analysis->logic_name($analysis->logic_name()."_not_best");
      $new_gene->analysis($analysis);
    }

    push(@$output, $new_gene) if($new_gene);
  }
  $self->output($output);
}



sub run_Genewise{
  my ($self, $features) = @_;

  if(!$features){
    $features = $self->features;
  }
  # foreach my $f(@$features){
  #  print "MINIGENEWISE ".$f->start." ".$f->end." ".$f->strand." ".$f->hseqname."\n";
  #}
  my $miniseq_object = $self->make_MiniSeq($features);
  my $miniseq = $miniseq_object->get_cDNA_sequence;
  my %params = %{$self->genewise_options} if($self->genewise_options);
  my $gw = new Bio::EnsEMBL::Analysis::Runnable::Genewise
      (
       -query   => $miniseq,
       -protein  => $self->protein_sequence,
       -reverse  => $self->is_reversed,
       -analysis => $self->analysis,
       %params
      );
  #throw("Minigenewise Can't run without context in sequence") 
  #  if(!$miniseq =~ /[CATG]{3}/);
  eval{
    $gw->run;
  };
  if($@){
    throw("Failed ".$gw." run $@"); # warning changed to throw() to prevent silent failing
  }
  my @output = @{$gw->output};

  return \@output;
}


sub convert_to_genomic{
  #my ($self, $gene) = @_;
  my ($self, $transcript) = @_;

  my $strand = 1;
  $strand = -1 if($self->is_reversed);
  my @converted_exons;
  my @all_supp_features;
  
  # need to convert all the exons and all the supporting features
  my @exons = @{$transcript->get_all_Exons};
  my $nexon = scalar(@exons);
  
 TEXON:for(my $i=0; $i < $nexon; $i++) {
    my $exon = $exons[$i];
    $exon->strand($strand);
    
    # need to convert whole exon back to genomic coordinates
   
    my @genomics = @{$self->miniseq->convert_SeqFeature($exon)};
 
    if(!@genomics){
      warning(id($exon)." ".coord_string($exon)." didn't produce any features ".
              "when converted to genomic coordinates, skipping");
      return undef;
    }
    if (@genomics > 1) {
      # all hell will break loose as the sub alignments will probably 
      # not map cheerfully and we may start introducing in frame stops ...
      # for now, ignore this feature.
      warning("Warning : feature converts into " . scalar(@genomics).
              " features; ");
      if ($i == $nexon - 1) {
        warning("Skipping last exon $i");
        next TEXON;
      } elsif ($exon->phase == $exons[$i+1]->phase) {
        warning("Skipping compatible exon $i");
        next TEXON;
      } else {
        warning("Cannot skip exon without adjusting phase; skipping gene");
        return undef;
      }
    }  else {
      my $genomic_exon = create_Exon($genomics[0]->start, $genomics[0]->end,
                                     $exon->phase, $exon->end_phase, $strand,
                                     $self->analysis, undef, undef, $self->query);
      # also need to convert each of the sub alignments back to genomic coordinates
      foreach my $sf (@{$exon->get_all_supporting_features}) {
        my @features;
        my @ungapped = $sf->ungapped_features;
        foreach my $aln (@ungapped) {
          my @alns = @{$self->miniseq->convert_PepFeaturePair($aln)};
          if ($#alns > 0) {
            warning("Warning : sub_align feature converts into > 1 features " . 
                    scalar(@alns));
          }
          my $align = create_feature_from_gapped_pieces
            ("Bio::EnsEMBL::DnaPepAlignFeature", \@alns, 100, undef, 
             undef, $self->query, $self->protein_sequence->id, $self->analysis);
          push @features,$align;
        }
        push @all_supp_features,@features;
        
        my $gapped = create_feature_from_gapped_pieces
          ("Bio::EnsEMBL::DnaPepAlignFeature", \@features, 100, undef, 
           undef, $self->query, $self->protein_sequence->id, $self->analysis);
        $genomic_exon->add_supporting_features($gapped);
      }
      push(@converted_exons,$genomic_exon);
    }
  }
  
  # make a new transcript from @converted_exons
  my $converted_transcript  = new Bio::EnsEMBL::Transcript;
  my $converted_translation = new Bio::EnsEMBL::Translation;
  $converted_transcript->translation($converted_translation);
  
  if (scalar(@all_supp_features)) {
    my $daf = create_feature_from_gapped_pieces
      ("Bio::EnsEMBL::DnaPepAlignFeature", \@all_supp_features, 100, undef, 
       undef, $self->query, $self->protein_sequence->id, $self->analysis);
    $converted_transcript->add_supporting_features($daf);
  }
  
  if ($#converted_exons < 0) {
    warning("Odd.  No exons found in MiniGenewise running on ".
            $self->query->seq_region_name." ".$self->protein_sequence->id);
    return undef;
  } else {
    if ($converted_exons[0]->strand == -1) {
      @converted_exons = sort {$b->start <=> $a->start} @converted_exons;
    } else {
      @converted_exons = sort {$a->start <=> $b->start} @converted_exons;
    }
    
    $converted_translation->start_Exon($converted_exons[0]);
    $converted_translation->end_Exon  ($converted_exons[$#converted_exons]);
    
    # phase is relative to the 5' end of the transcript (start translation)
    if ($converted_exons[0]->phase == 0) {
      $converted_translation->start(1);
    } elsif ($converted_exons[0]->phase == 1) {
      $converted_translation->start(3);
    } elsif ($converted_exons[0]->phase == 2) {
      $converted_translation->start(2);
    }
    
    $converted_translation->end  ($converted_exons[$#converted_exons]->end -
                                  $converted_exons[$#converted_exons]->start + 1);
    foreach my $exon(@converted_exons){
      $converted_transcript->add_Exon($exon);
    }
  }
  
  $converted_transcript->slice($self->query);

  return $converted_transcript;
}


=head2 make_MiniSeq

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniGenewise
  Arg [2]   : arrayref of Bio::EnsEMBL::FeaturePairs
  Function  : create MiniSeq object from the given feature pairs
  Returntype: Bio::EnsEMBL::Analysis::Tools::MiniSeq
  Exceptions: 
  Example   : 

=cut


sub make_MiniSeq {
  my ($self, $features) = @_;

  my $pairaln  = new Bio::EnsEMBL::Analysis::Tools::PairAlign;
  my $ff = $self->feature_factory;
  my @genomic_features;
		
  my $prevend     = 0;
  my $count       = 0;
  my $mingap      = $self->minimum_intron;
  my $seqname     = $features->[0]->seqname;

  my @features = sort {$a->start <=> $b->start} @$features;

 FEAT:foreach my $f (@features) {
    my $start = $f->start - $self->exon_padding;
    my $end   = $f->end   + $self->exon_padding;

    if ($start < 1) { 
      $start = 1;
    }

    if ($end   > $self->query->length) {
      $end = $self->query->length;
    }

    my $gap     =    ($start - $prevend);

    # Extend the region is the gap between features is
    # below a certain size - otherwise start a new region

    if ($count > 0 && ($gap < $mingap)) {

      if ($end < $prevend) { 
        $end = $prevend;
      }
      
      $genomic_features[$#genomic_features]->end($end);
      $prevend     = $end;
	
    } else {
	
      my $newfeature = new Bio::EnsEMBL::Feature;
	
      $newfeature->seqname   ($f->seqname);
      $newfeature->start     ($start);
      $newfeature->end       ($end);
      $newfeature->strand    (1);
      $newfeature->slice($self->query);
      push(@genomic_features,$newfeature);
	
      $prevend     = $end;
    }
    $count++;
  }

  # make a forward strand sequence, but tell genewise to run reversed if the 
  # features are on the reverse strand - handled by _is_reversed

  @genomic_features = sort {$a->start <=> $b->start } @genomic_features;

  # pad the termini of the sequence with 20K (configurable?) of genomic 
  # sequence - to catch small terminal exons that are missed by blast
  my $adjusted_start = $genomic_features[0]->start;
  $adjusted_start -= $self->terminal_padding;
  if($adjusted_start < 1 ) {
    $adjusted_start = 1;
  }
  $genomic_features[0]->start($adjusted_start);

  my $adjusted_end = $genomic_features[$#genomic_features]->end;
  $adjusted_end += $self->terminal_padding;
  if($adjusted_end > $self->query->length) {
    $adjusted_end = $self->query->length;
  }
  $genomic_features[$#genomic_features]->end($adjusted_end);

  my $current_coord = 1;
  foreach my $f (@genomic_features) {
    my $cdna_start = $current_coord;
    my $cdna_end   = $current_coord + ($f->end - $f->start);
    
    my $fp  = $ff->create_feature_pair($f->start, $f->end, $f->strand, undef,
                                       $cdna_start, $cdna_end, 1, 
                                       $f->seqname.".cdna", undef, undef, 
                                       $self->query->seq_region_name, $self->query, 
                                       $self->analysis);
    $pairaln->addFeaturePair($fp);
    $current_coord = $cdna_end+1;
  }

  my $miniseq = new Bio::EnsEMBL::Analysis::Tools::MiniSeq(-id        => 'test',
                                                           -pairalign => $pairaln);

  $self->miniseq($miniseq);
  return $miniseq;

}



=head2 is_reversed

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::MiniSeq
  Function  : Calculates if most of the features are on the
  forward or reverse strand
  Returntype: boolean
  Exceptions: warns if features are of mixed strand
  Example   : 

=cut



sub is_reversed {
  my ($self) = @_;

  if (!($self->{_reverse})) {
    my $strand = 0;
    my $fcount = 0;
    my $rcount = 0;

    foreach my $f (@{$self->features}) {
      if ($f->strand == 1){
	$fcount++;
      } elsif ($f->strand == -1) {
	$rcount++;
      }
    }

    if ($fcount > $rcount) {
      $self->{_reverse} = 0;
    } else {
      $self->{_reverse} = 1;
    }
    if($fcount && $rcount){
      warning("Forward and reverse strand features see to have been ".
              "passed to MiniGenewise together");
    }
  }
  return $self->{_reverse};
}


sub calculate_coverage_and_pid {
  my ($self, $value) = @_;
  if($value){
    $self->{_calculate_coverage_and_pid} = $value;
  }
  return $self->{_calculate_coverage_and_pid};
}


sub realign_translation {
  my ($self,$transcript) = @_;

  my $query_seq = $self->protein_sequence->seq();
#  my $transcripts = $gene->get_all_Transcripts();
#  my $transcript = ${$transcripts}[0];
  my $translation = $transcript->translate->seq();

  my $align_input_file = create_file_name('genewise_align', 'fa');
  my $align_output_file = create_file_name('genewise_align', 'aln');

  open(INPUT,">".$align_input_file);
  say INPUT ">query";
  say INPUT $query_seq;
  say INPUT ">target";
  say INPUT $translation;
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

  unless($file =~ /\>.+\n(([^>]+\n)+)\>.+\n(([^>]+\n)+)/) {
    throw("Could not parse the alignment file for the alignment sequences. Alignment file: ".$align_output_file);
  }

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

  # Get all exons and transcript supporting features
  my $transcript_supporting_features = $transcript->get_all_supporting_features();
  my $exons = $transcript->get_all_Exons();

  # Now clean these out
  $transcript->flush_Exons();
  $transcript->flush_supporting_features();

  # Loop through the TSFs and add the coverage and pid, then add back into transcript
  foreach my $transcript_supporting_feature (@{$transcript_supporting_features}) {
    $transcript_supporting_feature->hcoverage($coverage);
    $transcript_supporting_feature->percent_id($percent_id);
    $transcript->add_supporting_features($transcript_supporting_feature);
  }

  # Loop through exons, get supporting features for each, flush existing SF, add coverage and pid, add back to exon, add exon to transcript
  foreach my $exon (@{$exons}) {
    my $exon_supporting_features = $exon->get_all_supporting_features();
    $exon->flush_supporting_features();
    foreach my $exon_supporting_feature (@{$exon_supporting_features}) {
      $exon_supporting_feature->hcoverage($coverage);
      $exon_supporting_feature->percent_id($percent_id);
      $exon->add_supporting_features($exon_supporting_feature);
    }
    $transcript->add_Exon($exon);
  }

}

1;
