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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Set the default parameters
               methods_list => [],
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    methods_list => [],
    skip_analysis => 0,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Uses the 'group' parameter to decide which check need to be done on the 'input_db' database
              meta_table, analyze_tables and gene_sanity are done by default. If you specify an arrayref
              in 'methods_list' it overrides the defaults
              'core_handover' group add 'supporting_evidence_sanity'
              'otherfeatures_handover' adds no other test
              'rnaseq_handover' adds 'rnaseq_analysis_sanity', 'intron_supporting_evidence_sanity',
              'data_file_sanity' and 'dna_align_feature_sanity'
 Returntype : None
 Exceptions : Throws if the group is not known

=cut

sub fetch_input {
  my $self = shift;

  if($self->param('skip_analysis')) {
    $self->complete_early('Skip analysis flag is enabled, so no healthcheck will occur');
  }

  my $group = $self->param_required('group');
  my $hc_db = $self->get_database_by_name('input_db');
  $self->hrdb_set_con($hc_db, 'hc_db');
  my $methods_list = $self->param('methods_list');
  if (@$methods_list == 0) {
    push(@$methods_list, qw(
      meta_table
      analyze_tables
      gene_sanity
    ));
    if ($group eq 'core_handover') {
      push(@$methods_list, qw(
        supporting_evidence_sanity
        coding_supporting_evidence_presence
      ));
    }
    elsif ($group eq 'otherfeatures_handover') {
      push(@$methods_list, qw(
      ));
    }
    elsif ($group eq 'rnaseq_handover') {
      push(@$methods_list, qw(
        rnaseq_analysis_sanity
        data_file_sanity
        intron_supporting_evidence_sanity
        dna_align_feature_sanity
      ));
    }
    elsif ($group eq 'protein_cdna') {
      push(@$methods_list, qw(
        supporting_evidence_sanity
        supporting_evidence_presence
      ));
    }
    else {
      $self->throw("You are trying to use an unknown group '$group'");
    }
  }
  else {
    $self->say_with_header('Using specified list of checks: '.join(', ', @$methods_list));
  }
}


=head2 run

 Arg [1]    : None
 Description: Run the tests prepared in fetch_input. If a test fails, it stores the name of the test
              in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my $self = shift;

  foreach my $method_name (@{$self->param('methods_list')}) {
    $self->$method_name;
  }
}


=head2 write_output

 Arg [1]    : None
 Description: If at least of of the analyses fails, it prepares a red/green
              table to be displayed in html in GuiHive. Specifying -debug 1
              in the worker option will print the problems for each failing
              analysis
 Returntype : None
 Exceptions : Throws if at least one of the test fails

=cut

sub write_output {
  my $self = shift;

  my $failed_analyses = $self->output;
  if (@$failed_analyses) {
    my $html_msg;
    my %failed;
    foreach my $test (@{$self->output}) {
      $failed{$test} = 1;
    }
    foreach my $test (@{$self->param('methods_list')}) {
      my $color = exists $failed{$test} ? 'style="background-color:red";' : 'style="background-color:green";';
      $html_msg .= "<div $color><p>$test</p></div>";
    }
    $self->throw($html_msg);
  }
}


=head2 meta_table

 Arg [1]    : None
 Description: Check that meta keys species.url and species.production_name are present
 Returntype : None
 Exceptions : None

=cut

sub meta_table {
  my ($self) = @_;

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $meta_adaptor = $hc_db->get_MetaContainerAdaptor;
  my $meta_results = $meta_adaptor->list_value_by_key('species.production_name');
  my $species_production_name;
  my $failed = 0;
  if (@$meta_results) {
    if (@$meta_results == 1) {
      $species_production_name = $meta_results->[0];
      $self->param('production_name', $species_production_name);
    }
    else {
      $self->say_with_header('You have more than 1 meta key for "species.production_name":'.join(' ', @$meta_results));
    }
  }
  else {
    $failed = 1;
    $self->say_with_header('You are missing "species.production_name"');
  }
  $meta_results = $meta_adaptor->list_value_by_key('species.url');
  my $species_url;
  if (@$meta_results) {
    if (@$meta_results == 1) {
      $species_url = $meta_results->[0];
      if ($species_production_name) {
        $species_production_name =~ s/^(\w)/\U$1\E/;
        $species_production_name =~ s/gca(\d+)v(\d+)$/GCA_$1.$2/;
        $self->say_with_header("$species_production_name and $species_url are different, one of them could be wrong")
          unless ($species_production_name eq $species_url);
      }
    }
    else {
      $self->say_with_header('You have more than 1 meta key for "species.url":'.join(' ', @$meta_results));
    }
  }
  else {
    $failed = 1;
    $self->say_with_header('You are missing "species.url"');
  }
  if ($failed) {
    $self->output(['meta_table']);
  }
}


=head2 analyze_tables

 Arg [1]    : None
 Description: Run 'ANALYZE TABLE' on all tables of the database
 Returntype : None
 Exceptions : None

=cut

sub analyze_tables {
  my ($self) = @_;

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $dbc = $hc_db->dbc;
  my $sth = $dbc->db_handle->table_info(undef, undef, undef, 'TABLE');
  foreach my $row (@{$sth->fetchall_arrayref}) {
    $dbc->do('ANALYZE TABLE '.$row->[2]);
  }
}


=head2 coding_supporting_evidence_presence

 Arg [1]    : None
 Description: Checks that all coding transcript has a transcript supporting evidence
 Returntype : None
 Exceptions : None

=cut

sub coding_supporting_evidence_presence {
  my ($self) = @_;

  my @sql_queries = (
		     'SELECT t.* FROM transcript t LEFT JOIN transcript_supporting_feature tsf ON t.transcript_id = tsf.transcript_id LEFT JOIN analysis a ON t.analysis_id = a.analysis_id WHERE t.biotype NOT LIKE "%RNA" AND t.biotype NOT LIKE "ribozyme" AND t.biotype NOT LIKE "IG%" AND t.biotype NOT LIKE "TR%" AND tsf.transcript_id IS NULL AND a.logic_name NOT LIKE "mt_genbank_import" AND a.logic_name NOT LIKE "ncrna"',
  );

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $dbc = $hc_db->dbc;
  my $max_lines = 5;
  my $failed = 0;
  foreach my $query (@sql_queries) {
    my $sth = $dbc->prepare($query);
    $sth->execute;
    my $count = 0;
    my $msg = '';
    foreach my $row (@{$sth->fetchall_arrayref}) {
      $msg .= join("\t", @$row) if ($count++ < $max_lines);
    }
    if ($count) {
      $self->say_with_header("$count rows where it should be 0. You should rerun this query: $query");
      $failed = 1;
    }
  }
  if ($failed) {
    $self->output(['coding_supporting_evidence_presence']);
  }
}


=head2 supporting_evidence_presence

 Arg [1]    : None
 Description: Check that all transcripts and all exons have a supporting evidence
              This is useful for protein based alignemnt or cDNAs based alignment.
              Some transcriptomic sets will not have them
 Returntype : None
 Exceptions : None

=cut

sub supporting_evidence_presence {
  my ($self) = @_;

  my @sql_queries = (
    'SELECT t.* FROM transcript t LEFT JOIN transcript_supporting_feature tsf ON t.transcript_id = tsf.transcript_id WHERE tsf.transcript_id IS NULL',
    'SELECT e.* FROM exon e LEFT JOIN supporting_feature sf ON e.exon_id = sf.exon_id WHERE sf.exon_id IS NULL',
  );

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $dbc = $hc_db->dbc;
  my $max_lines = 5;
  my $failed = 0;
  foreach my $query (@sql_queries) {
    my $sth = $dbc->prepare($query);
    $sth->execute;
    my $count = 0;
    my $msg = '';
    my $small_exon = 0;
    foreach my $row (@{$sth->fetchall_arrayref}) {
      if ($row->[3]-$row->[2]+1 > 3) {
        $msg .= join("\t", @$row) if ($count++ < $max_lines);
      }
      else {
        $small_exon = 1;
      }
    }
    if ($count) {
      if ($small_exon) {
        $self->say_with_header("$count exons smaller than 3bp do not have a supporting evidence but it is OK");
      }
      else {
        $self->say_with_header("$count rows where it should be 0. You should rerun this query: $query");
        $failed = 1;
      }
    }
  }
  if ($failed) {
    $self->output(['supporting_evidence_presence']);
  }
}


=head2 supporting_evidence_sanity

 Arg [1]    : None
 Description: Run basic SQL queries to check that supporting evidences are correctly linked
 Returntype : None
 Exceptions : None

=cut

sub supporting_evidence_sanity {
  my ($self) = @_;

  my @sql_queries = (
    'SELECT t.transcript_id, t.seq_region_id, daf.seq_region_id FROM transcript t LEFT JOIN transcript_supporting_feature tsf ON t.transcript_id = tsf.transcript_id LEFT JOIN dna_align_feature daf ON daf.dna_align_feature_id = tsf.feature_id WHERE tsf.feature_type = "dna_align_feature" AND t.seq_region_id != daf.seq_region_id',
    'SELECT t.transcript_id, t.seq_region_id, paf.seq_region_id FROM transcript t LEFT JOIN transcript_supporting_feature tsf ON t.transcript_id = tsf.transcript_id LEFT JOIN protein_align_feature paf ON paf.protein_align_feature_id = tsf.feature_id WHERE tsf.feature_type = "protein_align_feature" AND t.seq_region_id != paf.seq_region_id',
    'SELECT t.transcript_id, t.seq_region_id, ise.seq_region_id FROM transcript t LEFT JOIN transcript_intron_supporting_evidence tsf ON t.transcript_id = tsf.transcript_id LEFT JOIN intron_supporting_evidence ise ON ise.intron_supporting_evidence_id = tsf.intron_supporting_evidence_id WHERE t.seq_region_id != ise.seq_region_id',
    'SELECT t.exon_id, t.seq_region_start, t.seq_region_end, daf.seq_region_start, daf.seq_region_end FROM exon t LEFT JOIN supporting_feature sf ON t.exon_id = sf.exon_id LEFT JOIN dna_align_feature daf ON daf.dna_align_feature_id = sf.feature_id WHERE sf.feature_type = "dna_align_feature" AND NOT (t.seq_region_start <= daf.seq_region_end AND t.seq_region_end >= daf.seq_region_start)',
    'SELECT t.exon_id, t.seq_region_start, t.seq_region_end, paf.seq_region_start, paf.seq_region_end FROM exon t LEFT JOIN supporting_feature sf ON t.exon_id = sf.exon_id LEFT JOIN protein_align_feature paf ON paf.protein_align_feature_id = sf.feature_id WHERE sf.feature_type = "protein_align_feature" AND NOT (t.seq_region_start <= paf.seq_region_end AND t.seq_region_end >= paf.seq_region_start)',
  );

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $dbc = $hc_db->dbc;
  my $max_lines = 5;
  my $failed = 0;
  foreach my $query (@sql_queries) {
    my $sth = $dbc->prepare($query);
    $sth->execute;
    my $count = 0;
    my $msg = '';
    foreach my $row (@{$sth->fetchall_arrayref}) {
      $msg .= join("\t", @$row) if ($count++ < $max_lines);
    }
    if ($count) {
      $self->say_with_header("$count rows where it should be 0. You should rerun this query: $query");
      $failed = 1;
    }
  }
  if ($failed) {
    $self->output(['supporting_evidence_sanity']);
  }
}


=head2 intron_supporting_evidence_sanity

 Arg [1]    : None
 Description: Check that all intron supporting evidences are link to a transcript and that all anaslyses
              are *_ise
 Returntype : None
 Exceptions : None

=cut

sub intron_supporting_evidence_sanity {
  my ($self) = @_;

  my @sql_queries = (
    'SELECT tise.* FROM transcript_intron_supporting_evidence tise LEFT JOIN transcript t ON t.transcript_id = tise.transcript_id WHERE t.transcript_id IS NULL',
  );

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $dbc = $hc_db->dbc;
  my $max_lines = 5;
  my $failed = 0;
  foreach my $query (@sql_queries) {
    my $sth = $dbc->prepare($query);
    $sth->execute;
    my $count = 0;
    my $msg = '';
    foreach my $row (@{$sth->fetchall_arrayref}) {
      $msg .= join("\t", @$row) if ($count++ < $max_lines);
    }
    if ($count) {
      $self->say_with_header("$count rows where it should be 0. You should rerun this query: $query");
      $failed = 1;
    }
  }
  my $sth = $dbc->prepare('SELECT a.logic_name FROM intron_supporting_evidence ise LEFT JOIN analysis a ON a.analysis_id = ise.analysis_id GROUP BY ise.analysis_id');
  $sth->execute;
  foreach my $row (@{$sth->fetchall_arrayref}) {
    if ($row->[0] !~ /ise$/) {
      $failed = 1;
      $self->say_with_header('Logic name "'.$row->[0].'" should not be in your intron_supporting_evidence table');
    }
  }
  my $iseadaptor = $hc_db->get_IntronSupportingEvidenceAdaptor;
  foreach my $analysis (@{$hc_db->get_AnalysisAdaptor->fetch_all}) {
    next unless ($analysis->logic_name =~ /ise$/);
    if (!$iseadaptor->generic_count('analysis_id = '.$analysis->dbID)) {
      $failed = 1;
      $self->say_with_header('Logic name "'.$analysis->logic_name.'" does not have data in the intron_supporting_evidence table');
    }
  }
  if ($failed) {
    $self->output(['intron_supporting_evidence_sanity']);
  }
}


=head2 dna_align_feature_sanity

 Arg [1]    : None
 Description: Check that:
                1. all the dna_align_feature rows are linked to an analysis whose logic name ends with "_daf"
                2. all the analysis rows whose logic name ends with "_daf" are linked to at least 1 row in the dna_align_feature table
              and output 'dna_align_feature_sanity' if the above checks fail. 
 Returntype : None
 Exceptions : None

=cut

sub dna_align_feature_sanity {
  my ($self) = @_;

  my @sql_queries = (
    'SELECT DISTINCT a.logic_name FROM analysis a RIGHT JOIN dna_align_feature daf ON a.analysis_id=daf.analysis_id WHERE a.logic_name NOT LIKE "%_daf" OR a.logic_name IS NULL', # check 1 (see description)
    'SELECT logic_name FROM analysis a LEFT JOIN dna_align_feature daf ON a.analysis_id=daf.analysis_id WHERE a.logic_name LIKE "%_daf" AND dna_align_feature_id IS NULL', # check 2 (see description)
  );

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $dbc = $hc_db->dbc();
  my $failed = 0;

  foreach my $query (@sql_queries) {
    my $sth = $dbc->prepare($query);
    $sth->execute();
    foreach my $row (@{$sth->fetchall_arrayref()}) {
      $row->[0] = "NULL" if (!($row->[0]));
      my $error_msg = "Analysis whose logic_name is ".$row->[0]." is not linked to the right set of dna_align_feature rows or viceversa.\n";
      $self->say_with_header($error_msg);
      $failed++;
    }
  }

  if ($failed) {
    $self->output(['dna_align_feature_sanity']);
  }
}


=head2 gene_sanity

 Arg [1]    : None
 Description: Run a basic set of SQL queries to check that all gene/transcript/exons are linked correctly
              Check the gene/transcripts/exons are on the same seq_region and overlapping each other
 Returntype : None
 Exceptions : None

=cut

sub gene_sanity {
  my ($self) = @_;

  my @sql_queries = (
    'SELECT g.* FROM gene g LEFT JOIN transcript t ON t.gene_id = g.gene_id WHERE t.gene_id IS NULL',
    'SELECT g.gene_id, g.seq_region_id, t.seq_region_id, t.transcript_id FROM gene g, transcript t WHERE t.gene_id = g.gene_id AND t.seq_region_id != g.seq_region_id',
    'SELECT g.gene_id, g.seq_region_id, g.seq_region_start, g.seq_region_end, t.seq_region_id, t.seq_region_start, t.seq_region_end, t.transcript_id FROM gene g, transcript t WHERE t.gene_id = g.gene_id AND t.seq_region_id = g.seq_region_id AND NOT (t.seq_region_start <= g.seq_region_end AND t.seq_region_end >= g.seq_region_start)',
    'SELECT t.* FROM transcript t LEFT JOIN exon_transcript et ON t.transcript_id = et.transcript_id WHERE et.transcript_id IS NULL',
    'SELECT et.* FROM exon_transcript et LEFT JOIN exon e ON et.exon_id = e.exon_id WHERE e.exon_id IS NULL',
    'SELECT tln.* FROM translation tln LEFT JOIN transcript t  ON tln.transcript_id = t.transcript_id WHERE t.transcript_id IS NULL',
    'SELECT t.* FROM transcript t LEFT JOIN translation tln ON t.transcript_id = tln.transcript_id WHERE t.biotype = "protein_coding" AND tln.transcript_id IS NULL',
    'SELECT tln.start_exon_id FROM translation tln LEFT JOIN exon_transcript et ON et.exon_id = tln.start_exon_id WHERE et.exon_id IS NULL',
    'SELECT tln.end_exon_id FROM translation tln LEFT JOIN exon_transcript et ON et.exon_id = tln.end_exon_id WHERE et.exon_id IS NULL',
    'SELECT t.transcript_id as TID, t.seq_region_id as TSR, e.seq_region_id as ESR, e.exon_id as EID FROM transcript t, exon_transcript et, exon e WHERE t.transcript_id = et.transcript_id AND et.exon_id = e.exon_id AND t.seq_region_id != e.seq_region_id',
    'SELECT t.transcript_id as TID, t.seq_region_id as TSR, e.seq_region_id as ESR, e.exon_id as EID FROM transcript t, exon_transcript et, exon e WHERE t.transcript_id = et.transcript_id AND et.exon_id = e.exon_id AND t.seq_region_id = e.seq_region_id AND NOT (t.seq_region_start <= e.seq_region_end AND t.seq_region_end >= e.seq_region_start)',
  );

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $dbc = $hc_db->dbc;
  my $max_lines = 5;
  my $failed = 0;
  foreach my $query (@sql_queries) {
    my $sth = $dbc->prepare($query);
    $sth->execute;
    my $count = 0;
    my $msg = '';
    foreach my $row (@{$sth->fetchall_arrayref}) {
      $msg .= join("\t", @$row) if ($count++ < $max_lines);
    }
    if ($count) {
      $self->say_with_header("$count rows where it should be 0. You should rerun this query: $query");
      $failed = 1;
    }
  }
  if ($failed) {
    $self->output(['gene_sanity']);
  }
}


=head2 rnaseq_analysis_sanity

 Arg [1]    : None
 Description: Check that 'other_protein' is present, check that all the other logic_names have 'rnaseq'
              in the name and that there is four similar analyses for each samples. The type is expected
              to be bam, daf, ise, gene.
 Returntype : None
 Exceptions : None

=cut

sub rnaseq_analysis_sanity {
  my ($self) = @_;

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $analysis_adaptor = $hc_db->get_AnalysisAdaptor;
  my %bases;
  my %types;
  my $failed = 1;
  my $total = 0;
  my $lc_count = 0;
  foreach my $analysis (@{$analysis_adaptor->fetch_all}) {
    ++$total;
    my $ln = $analysis->logic_name;
    if ($ln eq 'other_protein') {
      $failed = 0;
    }
    elsif ($ln =~ /rnaseq/) {
      my ($base, $type) = $ln =~ /^(\w+)_([^_]+)$/;
      ++$bases{$base};
      ++$types{$type};
    }
    else {
      $failed = 1;
      $self->say_with_header('You should not have this analysis: '.$ln);
    }
    if ($ln eq lc($ln)) {
      ++$lc_count;
    }
    else {
      $self->say_with_header("$ln has uppercase letters");
    }
  }
  if ($total != $lc_count) {
    $failed = 1;
    $self->say_with_header('Some of your analyses have uppercase characters');
  }
  if (4 != keys %types) {
    $failed = 1;
    $self->say_with_header('Some of your analyses are missing either gene, daf, bam or ise');
  }
  foreach my $base (keys %bases) {
    if ($bases{$base} != 4) {
      $failed = 1;
      $self->say_with_header('There is a problem with your analysis '.$base);
      $self->say_with_header('Check if gene, daf, bam and ise analyses are present.');
      $self->say_with_header('Total number of analyses in the database: '.$total);
    }
  }
  if ($failed) {
    $self->output(['rnaseq_analysis_sanity']);
  }
}


=head2 data_file_sanity

 Arg [1]    : None
 Description: Check that there is data in the data_file table
              and check that no logic_name is used for more than one file.
 Returntype : None
 Exceptions : None

=cut

sub data_file_sanity {
  my ($self) = @_;

  my $hc_db = $self->hrdb_get_con('hc_db');
  my $datafile_adaptor = $hc_db->get_DataFileAdaptor;
  my $analysis_adaptor = $hc_db->get_AnalysisAdaptor;
  my %bases;
  my %types;
  my $failed = 0;
  my $count = 0;
  foreach my $analysis (@{$analysis_adaptor->fetch_all}) {
    if ($analysis->logic_name =~ /_rnaseq_bam$/) {
      my $df = $datafile_adaptor->fetch_all_by_logic_name($analysis->logic_name);
      if (@$df == 0) {
        $failed = 1;
        $self->say_with_header('You are missing data_file for '.$analysis->logic_name);
      }
      elsif (@$df == 1) {
        ++$count;
      }
      else {
        $failed = 1;
        $self->say_with_header('You have more than 1 data_file for '.$analysis->logic_name);
      }
    }
  }
  $failed = 1 unless ($count);
  if ($failed) {
    $self->output(['data_file_sanity']);
  }
}

1;
