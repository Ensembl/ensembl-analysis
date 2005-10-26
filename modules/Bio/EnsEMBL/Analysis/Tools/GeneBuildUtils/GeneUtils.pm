
package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;

use strict;
use warnings;
use Exporter;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(print_Gene clone_Gene print_just_Gene);

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(print_Transcript clone_Transcript);

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw (id coord_string);
use Bio::EnsEMBL::Gene;


sub print_Gene{
  my $gene = shift;
  print_just_Gene($gene);
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    my $indent = "\t";
    print_Transcript($transcript, $indent);
  }
}


sub clone_Gene{
  my ($gene) = @_;
  my $newgene = Bio::EnsEMBL::Gene->new();
  $newgene->dbID($gene->dbID);
  $newgene->biotype($gene->biotype);
  $newgene->analysis($gene->analysis);
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    my $newtranscript = clone_Transcript($transcript);
    $newgene->add_Transcript($newtranscript);
  }
  $newgene->slice($newgene->slice);
  return $newgene;
}

sub print_just_Gene{
  my ($gene) = @_;
  my $coord_string = coord_string($gene);
  my $id = id($gene);
  print "GENE: ".$id." ".$coord_string." ".$gene->analysis->logic_name."\n";
}

1;
