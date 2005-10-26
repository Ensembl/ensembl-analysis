package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;

use strict;
use warnings;
use Exporter;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(coord_string id);


sub coord_string{
  my $feature = shift;
  throw("Must be passed a feature") if(!$feature);
  my $string = $feature->start." ".$feature->end." ".$feature->strand." ".$feature->slice->seq_region_name;
  return $string;
}



sub id{
  my $feature = shift;
  my $id;
  if($feature->can('stable_id')){
    $id = $feature->stable_id;
  }elsif($feature->dbID){
    $id = $feature->dbID;
  }else{
    $id = 'no-id';
  }
  if($feature->can('biotype') && $feature->biotype){
    $id .= $feature->biotype
  }
  return $id;
}

1;
