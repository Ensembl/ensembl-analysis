package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id);
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw (print_Evidence clone_Evidence Evidence_info create_feature 
              create_feature_from_gapped_pieces);


sub print_Evidence{
  my ($feature, $indent) = @_;
  print Evidence_info($feature, $indent)."\n";
}



sub Evidence_info{
 my ($feature, $indent) = @_;
  my $coord_string = coord_string($feature);
  my $id = id($feature);
  my $tag;
  if($feature->isa('Bio::EnsEMBL::DnaDnaAlignFeature')){
    $tag = "DNA EVIDENCE";
  }elsif($feature->isa('Bio::EnsEMBL::DnaPepAlignFeature')){
    $tag = "PROTEIN EVIDENCE";
  }else{
    throw("ExidenceUtils:print_Evidence Don't know what to do ">
          "with a ".$feature);
  }
  my $score = $feature->score || ".";
  my $percent_id = $feature->percent_id || ".";
  my $p_value = $feature->p_value || ".";
  $indent = '' if(!$indent);
  return $indent.$tag.": ".$coord_string." ".$score." ".
    $feature->hseqname." ".$feature->hstart." ".
      $feature->hend." ".$feature->hstrand." ".
        $percent_id." ".$p_value." ".
          $feature->cigar_string;
}


sub clone_Evidence{
  my ($feature) = @_;
  my $feature_string;
  if($feature->isa('Bio::EnsEMBL::DnaDnaAlignFeature')){
    $feature_string = 'Bio::EnsEMBL::DnaDnaAlignFeature';
  }elsif($feature->isa('Bio::EnsEMBL::DnaPepAlignFeature')){
    $feature_string = 'Bio::EnsEMBL::DnaPepAlignFeature';
  }else{
    throw("ExidenceUtils:clone_Evidence Don't know what to do ".
          "with a ".$feature);
  }

  my $newfeature = create_feature($feature_string, 
                                  $feature->start,
                                  $feature->end,
                                  $feature->strand,
                                  $feature->hstart,
                                  $feature->hend,
                                  $feature->hstrand,
                                  $feature->percent_id,
                                  $feature->p_value,
                                  $feature->score,
                                  $feature->hseqname,
                                  $feature->cigar_string,
                                  $feature->analysis,
                                  $feature->external_db_id,
                                  $feature->hcoverage,
                                  $feature->slice);
  $newfeature->dbID($feature->dbID);
  return $newfeature;

}


sub create_feature{
  my ($feature_string, $start, $end, $strand, $hstart, $hend, 
      $hstrand, $percent_id, $p_value, $score, $hseqname,
      $cigar_string, $analysis, $external_db_id, $hcoverage, $slice) = @_;

  my $feature = $feature_string->
    new(
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
        -external_db_id => $external_db_id,
        -hcoverage => $hcoverage,
        -cigar_string => $cigar_string,
       );
  $feature->slice($slice);
  $feature->hseqname($hseqname);
  return $feature;
}


sub create_feature_from_gapped_pieces{
  my ($feature_string, $pieces, $score, $percent_id, $p_value, $slice, $hseqname, 
      $analysis, $external_db_id, $hcoverage) = @_;
  
  my $feature = $feature_string->new(
                                     -features => $pieces,
                                     -percent_id => $percent_id,
                                     -score    => $score,
                                     -p_value  => $p_value,
                                     -hseqname => $hseqname,
                                     -analysis => $analysis,
                                     -external_db_id => $external_db_id,
                                     -hcoverage => $hcoverage,
                                    );
  $feature->slice($slice);
  $feature->hseqname($hseqname);
  $feature->analysis($analysis);
  return $feature;
}

#Decided against wrapping the get gene/transcript by evidence 
#methods as don't really need to the scripts themselves can do that
#themselves!


1;
