package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePseudopipeSetEnv;

use warnings;
use strict;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);

sub param_defaults {
    return {
      exon_mask_fields => '2 3',
            }
}

sub fetch_input {
  my($self) = @_;
  $self->create_analysis;
  $self->param_required('output_path');
  $self->param_required('tfasty_path');

}

sub write_output {
  my ($self) = @_;

  my $output_path = $self->param('output_path');
  my $tfasty_path = $self->param('tfasty_path');
  my $exon_mask_fields =  $self->param('exon_mask_fields');

# create output dir if it does not exist
  if (not -e $self->param('output_path')) {
     run_command("mkdir -p ".$self->param('output_path'),"Create output path.",0);
   }

  my $strand_array = [['P','plus'], ['M','minus']];

  foreach my $strand (@$strand_array){
    my $pseudo_dir = $self->param('output_path') ."/". $$strand[1];
    if (not -e $pseudo_dir) {
     run_command("mkdir -p ".$pseudo_dir,"Create ".$$strand[1]." path.",0);
     run_command("mkdir -p ".$pseudo_dir."/log","Create ".$$strand[1]." path.",0);
    }
    my $out_text = "dataDir=".$output_path."\nexport BlastoutSortedTemplate=".$output_path."/pep/%s_".$$strand[0]."_blastHits.sorted\nexport ChromosomeFastaTemplate=".$output_path."/dna/%s.fa\nexport ExonMaskTemplate=".$output_path."/mysql/%s_exLocs\nexport ExonMaskFields=\'".$exon_mask_fields."\'\nexport FastaProgram=".$tfasty_path."\nexport ProteinQueryFile=".$output_path."/pep/allPep.fa";
    open(ENVOUT, '>', $pseudo_dir.'/setenvPipelineVars');
    print ENVOUT $out_text;
  }
}


1;
