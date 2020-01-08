# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise

=head1 SYNOPSIS

my $obj = Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise->
    new(
        -max_exon_length => '20000',
        -multi_exon_min_coverage => '25',
        -single_exon_min_coverage => '80',
        -max_intron_length => '200000',
        -min_split_coverage => 95,
        -max_low_complexity => 101,
        );

 my ($accepted, $rejected) = $obj->filter_genes($genes);


=head1 DESCRIPTION

This module is designed to assess gene structures on several factors including 
coverage of evidence, size of introns and presence of complete ORF

The limits of these variables are defined in the constructure

The filter method then return 2 sets of genes, the first is the accepted set
the second is the rejected set

=head1 CONTACT

http://lists.ensembl.org/mailman/listinfo/dev

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 
  qw(convert_to_genes are_strands_consistent are_phases_consistent 
     is_not_folded has_no_unwanted_evidence all_exons_are_valid 
     intron_lengths_all_less_than_maximum list_evidence 
     evidence_coverage_greater_than_minimum Transcript_info
     split_Transcript tidy_split_transcripts evidence_coverage);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ConfigDependent::TranscriptUtils qw(low_complexity_less_than_maximum) ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(exon_length_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(clone_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id seq_region_coord_string 
                                                     lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils 
  qw(validate_Translation_coords contains_internal_stops print_Translation print_peptide);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Attribute;
use vars qw (@ISA);

@ISA = qw();


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise
  Function  : to create a BlastMiniGenewise filter object
  Returntype: Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise
  Exceptions: there are several args which musst be defined otherwise it
  will throw
  Example   : 

=cut



sub new{
  my ($class,@args) = @_;
  my $self = bless {},$class;
  my ($max_exon_length, $split_genes, $slice, $unwanted_evidence, 
      $max_intron_length, $min_coverage, $max_low_complexity, $seqfetcher,
      $single_exon_min_coverage, $min_split_coverage) 
    = rearrange([qw(MAX_EXON_LENGTH SPLIT_GENES SLICE UNWANTED_EVIDENCE 
                    MAX_INTRON_LENGTH MULTI_EXON_MIN_COVERAGE MAX_LOW_COMPLEXITY 
                    SEQFETCHER SINGLE_EXON_MIN_COVERAGE MIN_SPLIT_COVERAGE)], @args);
  #Setting defaults
  #$self->max_exon_length(20000);
  #  $self->split_multi_transcript_genes(0);
  #  $self->max_intron_length(200000);
  #  $self->min_coverage(25);
  #  $self->min_split_coverage(95);
  #  $self->single_exon_min_coverage(80);
  #  $self->max_low_complexity(101);
  #decided against setting defaults to avoid accidentally using the wrong
  #value
  #################

  $self->max_exon_length($max_exon_length);
  $self->split_multi_transcript_genes($split_genes);
  $self->slice($slice);
  $self->unwanted_evidence($unwanted_evidence);
  $self->min_coverage($min_coverage);
  $self->single_exon_min_coverage($single_exon_min_coverage);
  $self->max_low_complexity($max_low_complexity);
  $self->min_split_coverage($min_split_coverage);
  $self->seqfetcher($seqfetcher);
  $self->max_intron_length($max_intron_length);

  throw("Filter::BlastMiniGenewise needs a seqfetcher object to run") 
    if(!$self->seqfetcher);
  throw("Filter::BlastMiniGenewise needs a max exon length ") 
    if(! defined($self->max_exon_length));
  throw("Filter::BlastMiniGenewise needs a max intron length ") 
    if(! defined($self->max_intron_length));
  throw("Filter::BlastMiniGenewise needs a multi exon min coverage") 
    if(! defined($self->min_coverage));
  throw("Filter::BlastMiniGenewise needs a single exon min coverage")
    if(!defined($self->single_exon_min_coverage));
  throw("Filter::BlastMiniGenewise needs a max low complexity value ")
    if(!defined($self->max_low_complexity));
  throw("Filter::BlastMiniGenewise needs a min split coverage") 
    if(! defined($self->min_split_coverage));
  
  return $self;
}


#accessor methods


=head2 max_exon_length 

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise
  Arg [2]   : int/hash ref/array ref/object
  Function  : to store the given value
  Returntype: the give value
  Exceptions: some with throw if given the wrong thing
  Example   : 

=cut



sub max_exon_length{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{max_exon_length} = $arg;
  }
  return $self->{max_exon_length};
}

sub max_intron_length{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{max_intron_length} = $arg;
  }
  return $self->{max_intron_length};
}


sub max_low_complexity{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{max_low_complexity} = $arg;
  }
  return $self->{max_low_complexity};
}

sub min_coverage{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{min_coverage} = $arg;
  }
  return $self->{min_coverage};
}

sub min_split_coverage{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{min_split_coverage} = $arg;
  }
  return $self->{min_split_coverage};
}


sub single_exon_min_coverage{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{single_exon_min_coverage} = $arg;
  }
  return $self->{single_exon_min_coverage};
}

sub split_multi_transcript_genes{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{split_multi_transcript_gene} = $arg;
  }
  return $self->{split_multi_transcript_gene};
}

sub slice{
  my ($self, $arg) = @_;
  if(defined $arg){
    throw($arg." should be a Bio::EnsEMBL::Slice") 
      unless($arg->isa("Bio::EnsEMBL::Slice"));
    $self->{slice} = $arg;
  }
  return $self->{slice};
}

sub unwanted_evidence{
  my ($self, $arg) = @_;
  if(defined $arg){
    throw($arg." should be a hashref") 
      unless(ref($arg) eq "HASH");
    $self->{unwanted_evidence} = $arg;
  }
  return $self->{unwanted_evidence};
}


sub seqfetcher{
  my ($self, $arg) = @_;
  if($arg){
    throw("Filter::BlastMiniGenewise ".
          $arg." must have a method get_Seq_by_acc") 
      unless($arg->can("get_Seq_by_acc"));
    $self->{seqfetcher} = $arg;
  }
  return $self->{seqfetcher};
}


sub accepted_genes{
  my ($self, $arg) = @_;
  $self->{accepted_genes} = [] if(!$self->{accepted_genes});
  if($arg){
    my $temp;
    if(ref($arg) ne "ARRAY"){
      throw($arg." must be a Bio::EnsEMBL::Gene object") 
        unless($arg->isa("Bio::EnsEMBL::Gene"));
      $temp = [$arg];
    }else{
      $temp = $arg;
    }
    push(@{$self->{accepted_genes}}, @$temp);
  }
  return $self->{accepted_genes};
}

sub rejected_genes{
  my ($self, $arg) = @_;
  $self->{rejected_genes} = [] if(!$self->{rejected_genes});
  if($arg){
    my $temp;
    if(ref($arg) ne "ARRAY"){
      throw($arg." must be a Bio::EnsEMBL::Gene object") 
        unless($arg->isa("Bio::EnsEMBL::Gene"));
      $temp = [$arg];
    }else{
      $temp = $arg;
    }
    push(@{$self->{rejected_genes}}, @$temp);
  }
  return $self->{rejected_genes};
}


sub split_genes{
  my ($self, $arg) = @_;
  $self->{split_genes} = [] if(!$self->{split_genes});
  if($arg){
    my $temp;
    if(ref($arg) ne "ARRAY"){
      throw($arg." must be a Bio::EnsEMBL::Gene object") 
        unless($arg->isa("Bio::EnsEMBL::Gene"));
      $temp = [$arg];
    }else{
      throw($arg->[0]." must be a Bio::EnsEMBL::Gene object") 
        unless($arg->[0]->isa("Bio::EnsEMBL::Gene"));
      $temp = $arg;
    }
    push(@{$self->{split_genes}}, @$temp);
  }
  return $self->{split_genes};
}



=head2 filter_genes

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise
  Arg [2]   : Arrayref of Bio::EnsEMBL::Gene objects
  Function  : filter genes on defined characteristics
  Returntype: 2 arrayrefs, first accepted genes, then rejected
  Exceptions: 
  Example   : 

=cut



sub filter_genes{
  my ($self, $genes) = @_;

  my @accepted;
 GENE:foreach my $gene(@$genes, @{$self->split_genes}){
    my @transcripts = @{$gene->get_all_Transcripts};
    if(@transcripts > 1){
      my $cloned_gene = clone_gene($gene);
      warning("Filter::BlastMiniGenewise only works on single transcript ".
              "structures. Will convert ".$gene." into to one gene one ".
              "transcript structures");
      my $split_genes = convert_to_genes($cloned_gene->get_all_Transcripts);
      push(@$genes, @$split_genes) if($self->split_multi_transcript_genes);
      $self->split_genes($gene);
      next GENE;
    }
    my $is_valid = $self->validate_Transcript($transcripts[0]);
    my $hit_name = ${$transcripts[0]->get_all_supporting_features}[0]->hseqname;
    if($is_valid){
      logger_info("Have accepted ".id($gene)." $hit_name ".seq_region_coord_string($gene));
      push(@accepted, $gene);
    }else{
      logger_info("Have rejected ".id($gene)." $hit_name ".seq_region_coord_string($gene));
      $self->rejected_genes($gene);
    }
  }
  ACCEPTED:foreach my $accepted(@accepted){
      my @transcripts = @{$accepted->get_all_Transcripts};
      my $evidence = $self->get_Transcript_supporting_evidence($transcripts[0]);
      if(evidence_coverage_greater_than_minimum($transcripts[0], $evidence, ($self->min_split_coverage -1))){
        $self->accepted_genes($accepted);
        next ACCEPTED;
      }
      if(@{$transcripts[0]->get_all_Exons} == 1){
        $self->accepted_genes($accepted);
        next ACCEPTED;
      }  
 
      my $split_transcripts = split_Transcript($transcripts[0], $self->max_intron_length);
      my $acceptable_transcripts = tidy_split_transcripts($transcripts[0], $split_transcripts);
      my $genes = convert_to_genes($acceptable_transcripts, $accepted->analysis, $accepted->biotype);
      foreach my $gene(@$genes){
        my $attrib = $self->get_Attribute("split_tscript");
        foreach my $transcript(@{$gene->get_all_Transcripts}){
          $transcript->add_Attributes($attrib);
        }
        $self->accepted_genes($gene);
      }
    }
  return $self->accepted_genes, $self->rejected_genes;
}


=head2 validate_Transcript

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise
  Arg [2]   : Bio::EnsEMBL::Transcript
  Function  : validate the transcript on the basis of several factors
  Returntype: boolean, 1 if valid 0 if not
  Exceptions: 
  Example   : 

=cut



sub validate_Transcript{
  my ($self, $transcript) = @_;
  my $hit_name = ${$transcript->get_all_supporting_features}[0]->hseqname;
  #print "VALIDATING ".Transcript_info($transcript)."\n";
  my $slice = $self->slice;
  $slice = $transcript->slice if(!$slice);
  my $is_valid = 0;
  #basic transcript validation
  unless(are_strands_consistent($transcript)){
    $is_valid++;
    my $attrib = $self->get_Attribute("incons_strands");
    $transcript->add_Attributes($attrib);
  }
  #print "IS VALID is ".$is_valid." after strand consistency\n";
  unless(are_phases_consistent($transcript)){
    $is_valid++;
    my $attrib = $self->get_Attribute("incons_phases");
    $transcript->add_Attributes($attrib);
  }
  #print "IS VALID is ".$is_valid." phase consistency\n";
  unless(is_not_folded($transcript)){
    $is_valid++;
    my $attrib = $self->get_Attribute("is_folded");
    $transcript->add_Attributes($attrib);
  }
  #print "IS VALID is ".$is_valid." folded \n";
  unless(has_no_unwanted_evidence($transcript, 
                                  $self->unwanted_evidence)){
    $is_valid++;
    my $attrib = $self->get_Attribute("unwanted_evidence");
    $transcript->add_Attributes($attrib);
  }
  #print "IS VALID is ".$is_valid." after unwanted evidence\n";
  #$is_valid++ unless(all_exons_are_valid($transcript, $self->max_exon_length));
 EXON:foreach my $exon(@{$transcript->get_all_Exons}){
    if(exon_length_less_than_maximum($exon, $self->max_exon_length)){
      next EXON;
    }else{
      $is_valid++;
      last EXON;
      my $attrib = $self->get_Attribute("exon_too_long");
      $transcript->add_Attributes($attrib);
    }
  }
  #print "IS VALID is ".$is_valid." after exon validation\n";
  #basic translation validation
  #print_peptide($transcript);
  if(contains_internal_stops($transcript)){
    warning(Transcript_info($transcript)." $hit_name contains internal stop codons");
    $is_valid++;
    my $attrib = $self->get_Attribute("contains_stops");
    $transcript->add_Attributes($attrib);          
  }
  #print "IS VALID is ".$is_valid." after contains internal stops\n";
  unless(validate_Translation_coords($transcript)){
    $is_valid++;
    my $attrib = $self->get_Attribute("borked_coords");
    $transcript->add_Attributes($attrib);          
  }
  #print "IS VALID is ".$is_valid." after translation coord validation\n";
  #intron checks
  #$is_valid++ 
  #  unless(intron_lengths_all_less_than_maximum($transcript,                                                $self->max_intron_length));
  #print "IS VALID is ".$is_valid." intron length check\n";
  #checking low complexity
  
  unless(low_complexity_less_than_maximum($transcript, 
                                          $self->max_low_complexity)){
    $is_valid++;
    my $attrib = $self->get_Attribute("low_complex");
    $transcript->add_Attributes($attrib)
  }
  #print "IS VALID is ".$is_valid." low complexity check\n";
  #evidence coverage
  #note the coverage passed into the method is min - 1 as the method returns true if
  #the coverage is greater than the given value to allow the utils conventions to 
  #stand but in the old code the reverse was true it returned false if coverage
  #is less than the given so in order for this to produce the same results I needed
  #to add the min - 1 convention
  my $evidence = $self->get_Transcript_supporting_evidence($transcript);
  if(@{$transcript->get_all_Exons} >= 2){
    my $coverage = evidence_coverage($transcript, $evidence);
    unless(evidence_coverage_greater_than_minimum($transcript, $evidence, 
                                                  ($self->min_coverage - 1))){
      $is_valid++;
      my $attrib = $self->get_Attribute("evi_coverage");
      $transcript->add_Attributes($attrib);          
    }
  }else{
    unless(evidence_coverage_greater_than_minimum($transcript, $evidence, 
                                                  ($self->single_exon_min_coverage
                                                   -1))){
      $is_valid++;
      my $attrib = $self->get_Attribute("evi_coverage");
      $transcript->add_Attributes($attrib);          
    }
  }
  #print "IS VALID is ".$is_valid." after evidence coverage check\n";
  warning(Transcript_info($transcript)." $hit_name failed ".$is_valid." tests out of 9 ".
          "returning 0") if($is_valid >= 1);
  return 0 if($is_valid >= 1);
  return 1;
}

sub get_Transcript_supporting_evidence{
  my ($self, $transcript, $seqfetcher) = @_;

  $seqfetcher = $self->seqfetcher if(!$seqfetcher);
  throw("Can't get ".$transcript." transcrpts supporting features without ".
        "a seqfetcher object") if(!$seqfetcher);
  my $ids = list_evidence($transcript);
  my $sequence;
  foreach my $id(@$ids){
    $sequence = $seqfetcher->get_Seq_by_acc($id);
    last if($sequence);
  }
  return $sequence;
}

sub create_Attribute{
  my ($self, $reason, $description) = @_;
  $description = $reason if(!$description);
  my $attrib = Bio::EnsEMBL::Attribute->new(
                                            -code => $reason,
                                            -name => $reason,
                                            -description => $description,
                                            -value => $description,
                                           );
  return $attrib;
}

sub get_Attribute{
  my ($self, $code, $description) = @_;
  if(!$self->{'attribute_hash'}){
    $self->{'attribute_hash'} = {};
  }
  my $attrib = $self->{'attribute_hash'}->{$code};
  if(!$attrib){
    $attrib = $self->create_Attribute($code, $description);
    $self->{'attribute_hash'}->{$code} = $attrib;
  }
  return $attrib;
}

1;
