=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar2Introns

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar2Introns;

use warnings;
use strict;

use POSIX;
use Bio::EnsEMBL::Analysis;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Defaults parameters
                small_intron_size => 0,
                min_intron_depth => 0,
                table_name => 'csv_data',
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    small_intron_size => 0,
    min_intron_depth => 0,
    table_name => 'csv_data',
    has_merge_set => undef,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Check that sample_id_column and sample_column are set, check the databases
              and fetch the files containing the splce sites. It will also create the set
              analysis and cache the seq_regions
 Returntype : None
 Exceptions : Throws if no files where found

=cut

sub fetch_input {
  my ($self) = @_;

  $self->param_required('sample_id_column');
  $self->param_required('sample_column');
  my $output_dba = $self->get_database_by_name('intron_db');
  $self->hrdb_set_con($output_dba, 'target_db');

  my $junction_files;
  my $star_junctions_dir = $self->param_required('star_junctions_dir');
  $self->throw("The STAR junction dir provided does not exist:\n".$star_junctions_dir) unless (-d $star_junctions_dir);
  $junction_files = [glob($star_junctions_dir . '/*_SJ.out.tab')];

  if (scalar(@$junction_files)) {
    $self->param('_junction_files',$junction_files);
    my $csv_data_adaptor = $self->db->get_NakedTableAdaptor;
    $csv_data_adaptor->table_name($self->param('table_name'));
    my %analyses;
    my %unique_analyses;
    foreach my $row (@{$csv_data_adaptor->fetch_all()}) {
      my $logic_name = $self->param_required('species').'_'.$row->{$self->param('sample_column')}.'_rnaseq_daf';
      if (!exists $unique_analyses{$logic_name}) {
        $unique_analyses{$logic_name} = Bio::EnsEMBL::Analysis->new(-logic_name => $logic_name, -module => 'Star2Introns');
      }
      $analyses{$row->{$self->param('sample_id_column')}} = $unique_analyses{$logic_name};
    }
    #if (scalar(keys %unique_analyses) > 1) {
    #  my $logic_name = $self->param_required('species').'_merged_rnaseq_daf';
    #  $analyses{$logic_name} = Bio::EnsEMBL::Analysis->new(-logic_name => $logic_name, -module => 'Star2Introns');
    #  $unique_analyses{$logic_name} = $analyses{$logic_name};
    #  $self->param('has_merge_set', $logic_name);
    #}
    $self->param('_analyses', \%analyses);
    $self->param('_unique_analyses', \%unique_analyses);
    my $slice_adaptor = $output_dba->get_SliceAdaptor;
    $slice_adaptor->fetch_all('toplevel');
    $self->param('_slice_adaptor', $slice_adaptor);
  }
  else {
    $self->throw("No junction files found, something has gone wrong");
  }
}


=head2 run

 Arg [1]    : None
 Description: Read the splice junctions files and create DnaDnaAlign objects from it
 Returntype : None
 Exceptions : Throws if cannot open or close the files
              Throws if the accession is not found in the CSV table 'table_name'

=cut

sub run {
  my ($self) = @_;

  my $junction_files = $self->param('_junction_files');
  my $analyses = $self->param('_analyses');
  my %daf_table;
  foreach my $junction_file (@$junction_files) {
    $self->say_with_header("Parsing junction file: $junction_file");
    my $srr = (split'\/', $junction_file)[-1];
    $srr =~ s/_.*//;
    $self->throw("Sample name does not exist for ".$srr) unless (exists $analyses->{$srr});

    open(IN,$junction_file) or $self->throw("Coule not open $junction_file");
    while(<IN>) {
      my ($seq_region_name, $start, $end, $strand, $is_canonical, undef, $unique_map_count, $multi_map_count) = split("\t", $_);
      my $depth = $unique_map_count + POSIX::ceil($multi_map_count/2);
      my $intron_id = "$start:$end:$strand:$is_canonical";
      if (exists $daf_table{$seq_region_name}->{$intron_id}) {
        $daf_table{$seq_region_name}->{$intron_id}->{$analyses->{$srr}->logic_name} += $depth;
      }
      else{
        $daf_table{$seq_region_name}->{$intron_id}->{$analyses->{$srr}->logic_name} = $depth;
      }
    }
    close(IN) or $self->throw("Could not close $junction_file");
  }
  $self->param('_daf_table', \%daf_table);
}


=head2 write_output

 Arg [1]    : None
 Description: Store the splice site as DnaDnaAlignFeature
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;
  my $out_dbc = $self->hrdb_get_con('target_db')->dbc;
  my $daf_table = $self->param('_daf_table');
  my $analyses = $self->param('_unique_analyses');
  my $small_intron_size = $self->param('small_intron_size');
  my $min_intron_depth = $self->param('min_intron_depth');
  #my $merged_logic_name = $self->param('has_merge_set');

  my $target_adaptor = $self->hrdb_get_con('target_db')->get_AnalysisAdaptor;
  foreach my $analysis (values %$analyses) {
    $target_adaptor->store($analysis);
  }
  my $slice_adaptor = $self->param('_slice_adaptor');

  my $max_intron_size = 0;
  my $num_dafs = 0;
  $out_dbc->do("ALTER TABLE dna_align_feature DISABLE KEYS");
  my $sth_write = $out_dbc->prepare('INSERT INTO dna_align_feature (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, hit_start, hit_end, hit_strand, analysis_id, score, cigar_line, align_type) values (?, ?, ?, ?, ?, 1, ?, 1, ?, ?, ?, "ensembl")');
  foreach my $seq_region ( keys %$daf_table) {
    my $seq_region_id = $slice_adaptor->fetch_by_region('toplevel', $seq_region)->get_seq_region_id;
    foreach my $item (keys %{$daf_table->{$seq_region}}) {
      my ($start, $end, $strand, $is_canonical) = split(':', $item);
      #my $merged_depth = 0;
      my $diff = $end-$start;
      if ($diff > $max_intron_size) {
        $max_intron_size = $diff;
      }
      $is_canonical = $is_canonical > 0 ? 'canon' : 'non canon';
      foreach my $logic_name (keys %{$daf_table->{$seq_region}->{$item}}) {
        if ((($diff + 1) >= $small_intron_size) && $daf_table->{$seq_region}->{$item}->{$logic_name} >= $min_intron_depth) {
          $num_dafs++;
          $sth_write->bind_param(1, $seq_region_id);
          $sth_write->bind_param(2, $start);
          $sth_write->bind_param(3, $end);
          $sth_write->bind_param(4, $strand == 2 ? -1 : 1);
          $sth_write->bind_param(5, "$num_dafs:$is_canonical");
          $sth_write->bind_param(6, $diff);
          $sth_write->bind_param(7, $analyses->{$logic_name}->dbID);
          $sth_write->bind_param(8, $daf_table->{$seq_region}->{$item}->{$logic_name});
          $sth_write->bind_param(9, ($diff+1).'M');
          $sth_write->execute();
          #$merged_depth += $daf_table->{$seq_region}->{$item}->{$logic_name};
        }
        else {
          $self->say_with_header("Low threshold for $seq_region $item ".$analyses->{$logic_name}->logic_name.' '.$daf_table->{$seq_region}->{$item}->{$logic_name});
        }
      }
      #if ($merged_logic_name and $merged_depth) {
      #  $sth_write->bind_param(1, $seq_region_id);
      #  $sth_write->bind_param(2, $start);
      #  $sth_write->bind_param(3, $end);
      #  $sth_write->bind_param(4, $strand == 2 ? -1 : 1);
      #  $sth_write->bind_param(5, "$num_dafs:$is_canonical");
      #  $sth_write->bind_param(6, $diff);
      #  $sth_write->bind_param(7, $analyses->{$merged_logic_name}->dbID);
      #  $sth_write->bind_param(8, $merged_depth);
      #  $sth_write->bind_param(9, ($diff+1).'M');
      #  $sth_write->execute();
      #}
    }
  }
  $out_dbc->do("ALTER TABLE dna_align_feature ENABLE KEYS");
  my $sth_meta_coord = $out_dbc->prepare('INSERT INTO meta_coord (table_name, coord_system_id, max_length) VALUES ("dna_align_feature", 1, '.($max_intron_size+1).')');
  $sth_meta_coord->execute;
}


=head2 pre_cleanup

 Arg [1]    : None
 Description: Check indexes are 'on', truncate the dna_align_feature table and clean the meta_coord table
 Returntype : None
 Exceptions : None

=cut

sub pre_cleanup {
  my ($self) = @_;

  my $db = $self->get_database_by_name('intron_db');
  $self->warning('Truncating dna_align_feature for '.$db->dbc->dbname.'@'.$db->dbc->host);
  $db->dbc->do("ALTER TABLE dna_align_feature ENABLE KEYS");
  $db->dbc->do("TRUNCATE dna_align_feature");
  $db->dbc->do("DELETE FROM meta_coord WHERE table_name = 'dna_align_feature' AND coord_system_id = 1");
}

1;
