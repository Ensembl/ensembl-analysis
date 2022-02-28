# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

  Bio::EnsEMBL::Analysis::Tools::FeatureFactory

=head1 SYNOPSIS

  my $featurefactory = new Bio::EnsEMBL::Analysis::Tools::FeatureFactory;

  my $feature_pair = $featurefactory->create_feature_pair(1, 13, -1, 300, 
                                                          45, 56, 1 Q12732,
                                                          95, 2.3e-35, '',
                                                          $slice,
                                                          $analysis); 
  $featurefactory->validate($feature_pair);

=head1 DESCRIPTION

This is a utilities module which provides methods for feature creation
and feature validation for various feature types

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut


package Bio::EnsEMBL::Analysis::Tools::FeatureFactory;


use strict;
use warnings;


use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::RepeatConsensus;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::PredictionExon;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MarkerFeature;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Programs;


use vars qw (@ISA);

@ISA = qw();


sub new{
  my ($class,@args) = @_;
  my $self = bless {},$class;
  return $self;
}

#feature creation methods#



=head2 create_simple_feature

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : int, start
  Arg [3]   : int, end,
  Arg [4]   : int, stand (must be 0, 1 or -1)
  Arg [5]   : int score,
  Arg [6]   : string, display label,
  Arg [7]   : string, sequence name
  Arg [8]   : Bio::EnsEMBL::Slice
  Arg [9]   : Bio::EnsEMBL::Analysis
  Function  : creata a Bio::EnsEMBL::SimpleFeature
  Returntype: Bio::EnsEMBL::SimpleFeature
  Exceptions: 
  Example   : 

=cut



sub create_simple_feature{
  my ($self, $start, $end, $strand, $score, $display_label,
      $seqname, $slice, $analysis) = @_;
  my $simple_feature = Bio::EnsEMBL::SimpleFeature->new
    (
     -start => $start,
     -end => $end,
     -strand => $strand,
     -score => $score,
     -display_label => $display_label,
     -seqname => $seqname,
     -slice => $slice,
     -analysis => $analysis,
    );
  return $simple_feature;
}


=head2 create_simple_feature

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : string, name
  Arg [3]   : string, class
  Arg [4]   : string, type
  Arg [5]   : string, consensus sequence
  Arg [6]   : int length
  Function  : creata a Bio::EnsEMBL::RepeatConsensus
  Returntype: Bio::EnsEMBL::RepeatConsensus
  Exceptions: 
  Example   : 

=cut

sub create_repeat_consensus{
  my ($self, $name, $class, $type, $consensus_seq, $length) = @_;
  my $repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new
    (
     -name => $name,
     -length => $length,
     -repeat_class => $class,
     -repeat_consensus => $consensus_seq,
     -repeat_type => $type,
    );
  return $repeat_consensus;
}

=head2 create_simple_feature

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : int, start
  Arg [3]   : int, end,
  Arg [4]   : int, stand (must be 0, 1 or -1)
  Arg [5]   : int score,
  Arg [6]   : int, repeat_start
  Arg [7]   : int, repeat_end
  Arg [8]   : Bio::EnsEMBL::RepeatConsensus
  Arg [9]   : Bio::EnsEMBL::Slice
  Arg [10]   : Bio::EnsEMBL::Analysis
  Function  : creata a Bio::EnsEMBL::RepeatFeature
  Returntype: Bio::EnsEMBL::RepeatFeature
  Exceptions: 
  Example   : 

=cut

sub create_repeat_feature{
  my ($self, $start, $end, $strand, $score, $repeat_start, $repeat_end,
      $repeat_consensus, $seqname, $slice, $analysis) = @_;

  my $repeat_feature = Bio::EnsEMBL::RepeatFeature->new
    (
     -start            => $start,
     -end              => $end,
     -strand           => $strand,
     -slice            => $slice,
     -analysis         => $analysis,
     -repeat_consensus => $repeat_consensus,
     -hstart           => $repeat_start,
     -hend             => $repeat_end,
     -score            => $score,
     -seqname          => $seqname, 
    );
  return $repeat_feature;
}


=head2 create_feature_pair

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : int, start
  Arg [3]   : int, end
  Arg [4]   : int, strand must be 0, 1 or -1
  Arg [5]   : int, score
  Arg [6]   : int, hstart
  Arg [7]   : int, hend
  Arg [8]   : int, hstrand
  Arg [9]   : string, hseqname
  Arg [10]  : int, percent id
  Arg [11]  : int, p value
  Arg [12]  : string, seqname
  Arg [13]   : Bio::EnsEMBL::Slice
  Arg [14]  : Bio::EnsEMBL::Analysis
  Function  : creates a Bio::EnsEMBL::FeaturePair
  Returntype: Bio::EnsEMBL::FeaturePair
  Exceptions: 
  Example   : 

=cut


sub create_feature_pair {
    my ($self, $start, $end, $strand, $score, $hstart, $hend, 
        $hstrand, $hseqname, $percent_id, $p_value, $seqname,
        $slice, $analysis, $positive_matches, $identical_matches) = @_;
    my $fp = Bio::EnsEMBL::FeaturePair->new(
                                            -start    => $start,
                                            -end      => $end,
                                            -strand   => $strand,
                                            -hstart   => $hstart,
                                            -hend     => $hend,
                                            -hstrand  => $hstrand,
                                            -percent_id => $percent_id,
                                            -score    => $score,
                                            -p_value  => $p_value,
                                            -hseqname => $hseqname,
                                            -analysis => $analysis,
                                           );

    $fp->seqname($seqname);
    $fp->slice($slice);
    $fp->positive_matches($positive_matches);
    $fp->identical_matches($identical_matches);
    return $fp;
}

=head2 create_misc_feature

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : int, start,
  Arg [3]   : int, end
  Arg [4]   : int, strand
  Arg [5]   : Bio::EnsEMBL::Slice

=cut

sub create_misc_feature {

  my ($self, $start, $end, $strand, $slice) = @_;
  
  my $mf = Bio::EnsEMBL::MiscFeature->new
    (
     -start => $start,
     -end => $end,
     -strand => $strand,
     -slice => $slice,
    );
  return $mf;
}


=head2 add_misc_feature_attribute
  
  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : Bio::EnsEMBL::MiscFeature
  Arg [3]   : string, code
  Arg [4]   : string, name
  Arg [5]   : string, description
  Arg [6]   : string, value

=cut

sub add_misc_feature_attribute {

  my ($self, $mf, $code, $name, $description, $value) = @_;

  $mf->add_Attribute ( Bio::EnsEMBL::Attribute->new
    (-CODE   => $code,
     -NAME   => $name,
     -DESCRIPTION => $description,
     -VALUE  => $value,
    )
                     ); 
}
=head2 add_misc_set
  
  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : Bio::EnsEMBL::MiscFeature
  Arg [3]   : string, code
  Arg [4]   : string, name
  Arg [5]   : string, description
  Arg [6]   : string, value

=cut

sub add_misc_set {

  my ($self, $mf, $code, $name, $description, $longest_feature) = @_;

  $mf->add_MiscSet ( Bio::EnsEMBL::MiscSet->new
    (-CODE   => $code,
     -NAME   => $name,
     -DESCRIPTION => $description,
     -VALUE  => $longest_feature,
    )
                     ); 
}
  
=head2 create_prediction_exons

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : int, start,
  Arg [3]   : int, end
  Arg [4]   : int, strand
  Arg [5]   : int, score
  Arg [6]   : float, p value
  Arg [7]   : int, phase
  Arg [8]   : string, seqname
  Arg [9]   : Bio::EnsEMBL::Slice
  Arg [10]  : Bio::EnsEMBL::Analysis
  Function  : create a Bio::EnsEMBL::PredictionExon
  Returntype: Bio::EnsEMBL::PredictionExon
  Exceptions: 
  Example   : 

=cut



sub create_prediction_exon{
  my ($self, $start, $end, $strand, $score, $pvalue, $phase, $seqname,
      $slice, $analysis) = @_;
  my $exon = Bio::EnsEMBL::PredictionExon->new
    (
     -start => $start,
     -end => $end,
     -strand => $strand,
     -score => $score,
     -p_value => $pvalue,
     -phase => $phase,
     -slice => $slice,
     -seqname => $seqname,
     -analysis => $analysis,
     );
  return $exon;
}



=head2 create_prediction_transcript

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : arrayref, array of Bio::EnsEMBL::PredictionExons
  Arg [3]   : Bio::EnsEMBL::Slice
  Arg [4]   : Bio::EnsEMBL::Analysis
  Function  : 
  Returntype: 
  Exceptions: 
  Example   : 

=cut



sub create_prediction_transcript{
   my ($self, $exons, $slice, $analysis) = @_;
   my $transcript = Bio::EnsEMBL::PredictionTranscript->new
     (
      -exons => $exons,
      -slice => $slice,
      -analysis => $analysis,
     );
   return $transcript;
}



=head2 create_marker

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : int, database id
  Function  : create a marker with a specific database id
  Returntype: Bio::EnsEMBL::Map::Marker
  Exceptions: 
  Example   : 

=cut


sub create_marker{
  my ($self, $dbID) = @_;
  my $m = Bio::EnsEMBL::Map::Marker->new();
  $m->dbID($dbID);
  return $m;
}


=head2 create_marker_feature

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : int, start
  Arg [3]   : int, end
  Arg [4]   : int, strand
  Arg [5]   : Bio::EnsEMBL::Map::Marker
  Arg [6]   : string, seqname
  Arg [7]   : Bio::EnsEMBL::Slice
  Arg [8]  : Bio::EnsEMBL::Analysis
  Function  : 
  Returntype: Bio::EnsEMBL::Map::MarkerFeature
  Exceptions: 
  Example   : 

=cut


sub create_marker_feature{
  my ($self, $start, $end, $strand, $marker,
     $seqname, $slice, $analysis) = @_;
  my $mf = Bio::EnsEMBL::Map::MarkerFeature->new();
  $mf->start($start);
  $mf->end($end);
  $mf->strand($strand);
  $mf->marker($marker);
  $mf->seqname($seqname);
  $mf->slice($slice);
  $mf->analysis($analysis);
  return $mf;
}


#validation methods#

=head2 validate

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Arg [2]   : Bio::EnsEMBL::Feature
  Function  : validates feature
  Returntype: Bio::EnsEMBL::Feature
  Exceptions: throws if no slice or analysis is defined
  if the start, end or strand arent defined, if start or end are
  less than one or if start is greater than end
  Example   : 

=cut


sub validate{
  my ($self, $feature) = @_;
  # print STDERR "validating: ".$feature->start."-".$feature->end.":".$feature->hseqname."::".$feature->hstart."-".$feature->hend."\n";
  my @error_messages;
  if(!$feature){
    throw("Can't validate a feature without a feature ".
          "FeatureFactory::validate");
  }
  if(!($feature->isa('Bio::EnsEMBL::Feature'))){
    throw("Wrong type ".$feature." must be a Bio::EnsEMBL::Feature ".
          "object FeatureFactory::validate");
  }
  if(not defined $feature->slice){
    my $string = "No slice defined";
    push(@error_messages, $string);
  }
  if(not defined $feature->analysis){
    my $string = "No analysis defined";
    push(@error_messages, $string);
  }
  if(not defined $feature->start){
    my $string = "No start defined";
    push(@error_messages, $string); 
  }
  if(not defined $feature->end){
    my $string = "No end defined";
    push(@error_messages, $string); 
  }
  if(not defined $feature->strand){
    my $string = "No strand defined";
    push(@error_messages, $string); 
  }
  if($feature->start > $feature->end){
    my $string = "Start is greater than end ".$feature->start." ".
      $feature->end;
    push(@error_messages, $string); 
  }
  if(@error_messages > 0){
    print STDERR join("\n", @error_messages);
    throw("Invalid feature ".$feature." FeatureFactory:validate");
  }
}


sub validate_prediction_transcript{
  my ($self, $pt, $attach_to_exons) = @_;
  if(!$pt){
    throw("Can't validate a prediction transcript without a ".
          "prediction transcript ".
          "FeatureFactory:validate_prediction_transcript");
  }
  if(!($pt->isa('Bio::EnsEMBL::PredictionTranscript'))){
    throw("Wrong type ".$pt." must be a Bio::EnsEMBL::PredictionTranscript".
          "FeatureFactory:validate_prediction_transcript");
  }
  my @exons = @{$pt->get_all_Exons};
  if(@exons == 0){
    throw("problem ".$pt." has no exons");
  }
  foreach my $e(@exons){
    if($attach_to_exons){
      $e->slice($pt->slice) if(!$e->slice);
      $e->analysis($pt->analysis) if(!$e->analysis);
    }
    $self->validate($e);
  }
  $self->validate($pt);
  my $tseq;
  eval{
    $tseq = $pt->translate;
  };
  if(!$tseq){
    throw($pt." translate didn't return a sequence $@");
  }
  if ($tseq->seq =~ /\*/) {
    my $msg = $pt." doesn't have a valid translation ".$tseq->seq;
    $msg .= "\n$@" if($@);
    throw($msg);
  }
  return 1;
}

1;
