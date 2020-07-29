package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar2Introns;

use warnings;
use strict;
use feature 'say';
use POSIX;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::IntronSupportingEvidence;

sub param_defaults {
  my ($self) = @_;

  return {
	  %{$self->SUPER::param_defaults},
    small_intron_size => 0,
    min_intron_depth => 0,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: 
 Returntype : None
 Exceptions : 

=cut

sub fetch_input {
  my ($self) = @_;

  my $source_db = $self->get_database_by_name('source_db');
  $self->hrdb_set_con($source_db,'source_db');
  $self->hrdb_set_con($source_db,'target_db');

  my $output_dba = $self->hrdb_get_dba($self->param('intron_db'));
  $self->param('_output_dba',$output_dba);

  if($self->param('use_genome_flatfile')) {
    say "Ingoring dna table and using fasta file for sequence fetching";
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
    } else {
    say "Attaching dna db";
    my $dna_dba = $self->get_database_by_name('dna_db');
    $source_db->dnadb($dna_dba);
  }

  my $junction_files;
  if($self->param('star_junctions_dir')) {
    my $star_junctions_dir = $self->param('star_junctions_dir');
    unless(-e $star_junctions_dir) {
      $self->throw("The STAR junction dir param was provided but the dir does not exist Path provided:\n".$star_junctions_dir);
    }
    $junction_files = [glob($star_junctions_dir . '/*_SJ.out.tab')];
  }
  unless(scalar(@$junction_files)) {
    $self->throw("No junction files found, something has gone wrong");
  }
  $self->param('_junction_files',$junction_files);

  my $slice = $source_db->get_SliceAdaptor->fetch_by_name($self->input_id);
  my $exons = $slice->get_all_Exons();

  $self->param('_exons_on_slice',$exons);

}

sub run {
  my ($self) = @_;

  my $out_dbc = $self->param('_output_dba')->dbc;

  my $junction_files = $self->param('_junction_files');
  my $small_intron_size = $self->param('small_intron_size');
  my $min_intron_depth = $self->param('min_intron_depth');

  my $analysis_id=9000;
  foreach my $junction_file (@$junction_files) {
    unless(-e $junction_file) {
      $self->throw("The STAR splice junction file in the input id does not exist. Path provided:\n".$junction_file);
    }

    print "Parsing junction file: ".$junction_file."\n";

    open(IN,$junction_file);
    while(<IN>) {
      my $line = $_;
      my @elements = split("\t",$line);
      my $seq_region_name = $elements[0];
      my $start = $elements[1];
      my $end = $elements[2];
      my $strand = $elements[3];
      if($strand == 1) {
        next;
      } elsif($strand == 2) {
        $strand = -1;
      }
      my $is_canonical = $elements[4];
      if($is_canonical > 0) {
        $is_canonical = 1;
      } elsif($strand == 0) {
        $strand = 0;
      }
      my $unique_map_count = $elements[6];
      my $multi_map_count = POSIX::ceil($elements[7] / 2);
      my $depth = $unique_map_count + $multi_map_count;

      my $hit_seq_region_name = $seq_region_name;
      $hit_seq_region_name =~ s/\./\*/;
      my $hit_name = $seq_region_name.":".$start.":".$end.":".$strand;
      if ($is_canonical){
        $hit_name = $hit_name.":canon";
      }

      unless((($end - $start + 1) >= $small_intron_size) && $depth >= $min_intron_depth) {
        next;
      }

      my $id_sql = "SELECT seq_region_id from seq_region where name=?";
      my $sth = $out_dbc->prepare($id_sql);
      $sth->bind_param(1,$seq_region_name);
      $sth->execute();
      my $seq_region_id = $sth->fetchrow_array();

#      say 'INSERT INTO intron_supporting_evidence (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, score, score_type, is_splice_canonical) values ('.$seq_region_id.', '.$start.', '.$end.', '.$strand.', "'.$hit_name.'", '.$depth.', "DEPTH", '.$is_canonical.');';

      my $sth_write = $out_dbc->prepare('INSERT INTO intron_supporting_evidence (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, score, score_type, is_splice_canonical) values ('.$analysis_id.', '.$seq_region_id.', '.$start.', '.$end.', '.$strand.', "'.$hit_name.'", '.$depth.', "DEPTH", '.$is_canonical.');');
      eval {
        $sth_write->execute();
      };
      if($@){
        say "failed";
      }

# to create new intron need 5' and 3' exons
#      my $new_intron = Bio::EnsEMBL::Intron->new(
#        -SEQ_REGION_NAME => $seq_region_name,
#        -START           => $start,
#        -END             => $end,
#        -STRAND          => $strand,
#        -IS_SPLICE_CANONICAL=> $is_canonical,
#        -SCORE           => $depth,
#        -SCORE_TYPE      => "DEPTH",
#        -HIT_NAME        => $hit_name,
#      );

#      $self->output([$new_intron]);
   }
    close IN;
    $analysis_id++;
  } # End foreach my $junction_file
}

sub write_output {
  my ($self) = @_;

#  my $introns = $self->output;
#  my $out_dba = $self->hrdb_get_con('output_db');

#  my $intron_adaptor = $out_dba->get_IntronAdaptor;
#  foreach my $intron ( @{$introns} ) {
#    empty_Intron($intron);
#    $intron_adaptor->store($intron);
#  }

  return 1;
}
1;
