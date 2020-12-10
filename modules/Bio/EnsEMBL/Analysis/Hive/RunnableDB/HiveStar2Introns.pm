package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar2Introns;

use warnings;
use strict;
use feature 'say';
use POSIX;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::IntronSupportingEvidence;
use Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses;

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
#  $self->hrdb_set_con($source_db,'target_db');

  my $output_dba = $self->hrdb_get_dba($self->param('intron_db'));
  $self->param('_output_dba',$output_dba);

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
}

sub run {
  my ($self) = @_;

  my $out_dbc = $self->param('_output_dba')->dbc;

  my $junction_files = $self->param('_junction_files');
  my $small_intron_size = $self->param('small_intron_size');
  my $min_intron_depth = $self->param('min_intron_depth');

  my $ise_table = {};
  my @analyses;
  foreach my $junction_file (@$junction_files) {
    unless(-e $junction_file) {
      $self->throw("The STAR splice junction file in the input id does not exist. Path provided:\n".$junction_file);
    }

    print "Parsing junction file: ".$junction_file."\n";

    my $csv_data_adaptor = $self->db->get_NakedTableAdaptor;
    $csv_data_adaptor->table_name('csv_data');
    my $srr = (split'\/', $junction_file)[-1];
    $srr =~ s/_.*//;
    my $results = $csv_data_adaptor->fetch_all();
    my $logic_name;
    for my $row (@$results){
      if ($row->{$self->param('sample_id_column')} eq $srr){
	$logic_name = $row->{$self->param('sample_column')}
      }
    }
    if ($logic_name){
      $logic_name .= "_ise";
      my $analysis = {'-logic_name' => $logic_name, '-module' => 'Star2Introns'};
      push (@analyses, Bio::EnsEMBL::Analysis->new(%$analysis));
    }
    else{
      $self->throw("Sample name does not exist for ".$srr);
    }
    $self->param('_analyses', \@analyses);

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

      if (exists ($ise_table->{$hit_name})){
        $ise_table->{$hit_name}->{'score'} += $depth;
      }
      else{
        $ise_table->{$hit_name}->{'seq_region_id'} = $seq_region_id;
        $ise_table->{$hit_name}->{'seq_region_start'} = $start;
        $ise_table->{$hit_name}->{'seq_region_end'} = $end;
        $ise_table->{$hit_name}->{'seq_region_strand'} = $strand;
        $ise_table->{$hit_name}->{'score'} = $depth;
        $ise_table->{$hit_name}->{'is_splice_canonical'} = $is_canonical;
        $ise_table->{$hit_name}->{'logic_name'} = $logic_name;
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
    } # End while
    close IN;
  } # End foreach my $junction_file

  $self->param('_ise_table', $ise_table);
}

sub write_output {
  my ($self) = @_;

  my $out_dbc = $self->param('_output_dba')->dbc;
  my $ise_table = $self->param('_ise_table');
  my $analyses = $self->param('_analyses');

  my $target_adaptor = $self->param('_output_dba')->get_AnalysisAdaptor;
  foreach my $analysis (@{$analyses}) {
    $target_adaptor->store($analysis);
  }

  foreach my $hit_name ( keys $ise_table){
    my $logic_name = $ise_table->{$hit_name}->{'logic_name'};
    my $analysis_id = $target_adaptor->fetch_by_logic_name($logic_name)->dbID();
    my $seq_region_id = $ise_table->{$hit_name}->{'seq_region_id'};
    my $start = $ise_table->{$hit_name}->{'seq_region_start'};
    my $end = $ise_table->{$hit_name}->{'seq_region_end'};
    my $strand = $ise_table->{$hit_name}->{'seq_region_strand'};
    my $depth = $ise_table->{$hit_name}->{'score'};
    my $is_canonical = $ise_table->{$hit_name}->{'is_splice_canonical'};

    my $sth_write = $out_dbc->prepare('INSERT INTO intron_supporting_evidence (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, score, score_type, is_splice_canonical) values ('.$analysis_id.', '.$seq_region_id.','.$start.', '.$end.', '.$strand.', "'.$hit_name.'", '.$depth.', "DEPTH", '.$is_canonical.');');
    eval {
      $sth_write->execute();
    };
    if($@){
      say 'INSERT INTO intron_supporting_evidence (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, score, score_type, is_splice_canonical) values ('.$analysis_id.', '.$seq_region_id.','.$start.', '.$end.', '.$strand.', "'.$hit_name.'", '.$depth.', "DEPTH", '.$is_canonical.');';
      say "Failed to insert intron supporting evidence ".$hit_name."\n";
    }
  }
  return 1;
}
1;
