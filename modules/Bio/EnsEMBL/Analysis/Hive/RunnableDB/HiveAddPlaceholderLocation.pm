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

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddPlaceholderLocation;

use warnings;
use strict;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

 Arg [1]    : None
 Description: Set core db connection
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;

  my $core_dba = $self->hrdb_get_dba($self->param_required('input_db'));
  $self->hrdb_set_con($core_dba, 'input_db');
  my $support_test = $core_dba->dbc->prepare('select count(*) from supporting_feature');
  $support_test->execute;

  my $test = $support_test->fetchrow_array;
  $self->{'has_supporting'} = 0;
  if ($test > 0) {
    $self->{'has_supporting'} = 1;
  }
}

=head2 run

 Arg [1]    : None
 Description: Choose a transcript for the placeholder location
              (checks 10 longest seq_regions for transcript with high support,
              i.e coverage>=99 and percent_id>=75. If no supporting evidence
              available, uses largest gene by number of exons.)
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $core_dba = $self->hrdb_get_con('input_db');
  my $sa = $core_dba->get_SliceAdaptor;
  my $sth_longest = $core_dba->dbc->prepare('select seq_region_id from seq_region order by length desc limit 10');
  $sth_longest->execute;
  my $sample_transcript;

  if ($self->{'has_supporting'} == 1) {

    LOOP: while (my $seq_region_id = $sth_longest->fetchrow_array) {
      my $region = $sa->fetch_by_seq_region_id($seq_region_id);
      TRANSCRIPT:foreach my $transcript (@{ $region->get_all_Transcripts_by_type('protein_coding') }) {
        my $supporting_features = $transcript->get_all_supporting_features;
        foreach my $support (@$supporting_features) {
          if ($support->hcoverage() >= 99 && $support->percent_id() >= 75) {
            $sample_transcript=$transcript;
	          last LOOP;
          }
          else {
            next TRANSCRIPT;
          }
        }
      }
    }#end while
  }

  else {
    my $most_exons = 0;
    my $selected_gene;

    while (my $seq_region_id = $sth_longest->fetchrow_array) {
      my $region = $sa->fetch_by_seq_region_id($seq_region_id);
      foreach my $gene (@{$region->get_all_Genes_by_type('protein_coding')}) {
        my $exons = $gene->get_all_Exons;
        if (scalar(@$exons) > $most_exons) {
          $most_exons = scalar(@$exons);
    	    $selected_gene = $gene;
        }
      }
    }

    $sample_transcript = $selected_gene->canonical_transcript;
  }

  if ($sample_transcript) {
    my $db_name = $core_dba->dbc->dbname;
    my $sample_gene = $sample_transcript->get_Gene;
    my $sample_coord = $sample_gene->seq_region_name().':'.$sample_gene->seq_region_start().'-'.$sample_gene->seq_region_end();

    my @output = (['sample.location_param', $sample_coord],
                  ['sample.location_text', $sample_coord],
                  ['sample.gene_param', $sample_gene->stable_id()],
                  ['sample.gene_text', $sample_gene->stable_id()],
                  ['sample.transcript_param', $sample_transcript->stable_id()],
                  ['sample.transcript_text', $sample_transcript->stable_id()]);

    $self->output(\@output);
  }
  else {
    self->throw("Could not obtain a suitable placeholder transcript.")
  }
}

=head2 write_output

 Arg [1]    : None
 Description: Executes mysql commands to add sample location info to the meta table
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;
  my $core_dba = $self->hrdb_get_con('input_db');
  my $meta_container = $core_dba->get_MetaContainer();
  foreach my $meta_pair ( @{$self->output} ){
    my $meta_key = $meta_pair->[0];
    my $meta_value = $meta_pair->[1];
    $meta_container->store_key_value($meta_key, $meta_value);
  }
}

1;
