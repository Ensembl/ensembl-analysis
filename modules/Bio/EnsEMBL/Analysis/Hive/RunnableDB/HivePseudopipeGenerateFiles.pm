package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePseudopipeGenerateFiles;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

sub fetch_input {
  my($self) = @_;
  $self->create_analysis;

  my $dna_dba = $self->get_database_by_name('dna_db');

  my $slice = $dna_dba->get_SliceAdaptor->fetch_by_name($self->input_id);
  $self->param('_input_slice',$slice);
}

sub write_output {
  my ($self) = @_;
  my $slice = $self->param('_input_slice');
  my $output_path = $self->param('output_path');

# create peptide file
  my $pepfile = "allPep.fa";
  my $pep_dir = $self->param('output_path') ."/pep";
  if (not -e $pep_dir) {
     run_command("mkdir -p ".$pep_dir,"Create pep path.",0);
   }
  open( pepOUT, ">>", $pep_dir . "/" . $pepfile );

# for each slice
  say "Slice: ".$slice->name;

# print exon locations to file per slice
  my @exons = @{$slice->get_all_Exons};
  my $mysql_dir = $self->param('output_path') ."/mysql";
  if (not -e $mysql_dir) {
    run_command("mkdir -p ".$mysql_dir,"Create mysql path.",0);
  }
  open( exOUT, ">", $mysql_dir . "/" . $slice->name."_exLocs" );

  foreach my $exon (@exons) {
    my $seqname = $exon->seqname();
    my $length = $exon->length();
    my $start = $exon->start();
    my $end = $exon->end();

    say exOUT $seqname ."\t". $length ."\t". $start ."\t". $end;
  }

# print slice dna to file per slice
  my $dna_dir = $self->param('output_path') ."/dna";
  if (not -e $dna_dir) {
    run_command("mkdir -p ".$dna_dir,"Create dna path.",0);
  }

  open( dnaOUT, ">", $dna_dir . "/" . $slice->name . ".dna.fa" );
  my $dna_seq = $slice->seq;

  say dnaOUT ">".$slice->name."\n".$dna_seq;

# print all peptides to allPep.fa
  my @genes = @{$slice->get_all_Genes};
  foreach my $gene (@genes) {
    my $transcript = $gene->canonical_transcript;

    if ( $transcript->translation ) {
      my $seq_reg_name = $transcript->seqname();
      my $gene_stable_id = $transcript->get_Gene->stable_id();
      my $transcript_stable_id = $transcript->stable_id();
      my $pep_stable_id = $transcript->translation->stable_id;
      my $pep_seq = $transcript->translation->seq();

      say pepOUT ">".$pep_stable_id." pep:novel ".$seq_reg_name." gene:".$gene_stable_id." transcript:".$transcript_stable_id."\n".$pep_seq;
    }
  }
  close exOUT;
  close dnaOUT;
  close pepOUT;
}


1;
