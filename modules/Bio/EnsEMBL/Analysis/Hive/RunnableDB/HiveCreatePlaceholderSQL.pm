=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreatePlaceholderSQL;

use warnings;
use strict;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

 Arg [1]    : None
 Description: Implements fetch_input() interface method of Bio::EnsEMBL::Hive::Process that is used to read in parameters and load data.
              Here we have nothing to do.
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
}

=head2 run

 Arg [1]    : None
 Description: Creates an sql command add placeholder location information to the meta table
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $core_dba = $self->hrdb_get_dba($self->param_required('input_db'));
  my $sa = $core_dba->get_SliceAdaptor;

  my $sth_longest = $core_dba->dbc->prepare('select seq_region_id from seq_region order by length desc limit 10');
  $sth_longest->execute;

  my $sample_transcript;
 LOOP: while (my $seq_region_id = $sth_longest->fetchrow_array) {
    my $region = $sa->fetch_by_seq_region_id($seq_region_id);

  TRANSCRIPT:foreach my $transcript (@{ $region->get_all_Transcripts_by_type('protein_coding') }){
      print Dumper $transcript;
      my $supporting_features = $transcript->get_all_supporting_features;
      foreach my $support (@$supporting_features){
        if ($support->hcoverage() >= 99 && $support->percent_id() >= 75){
          $sample_transcript=$transcript;
          last LOOP;
        }
        else{
          next TRANSCRIPT;
        }
      }
    }
  }#end while

  if ($sample_transcript){
    my $db_name = $core_dba->dbc->dbname;
    my $sample_gene = $sample_transcript->get_Gene;
    my $sample_coord = $sample_gene->seq_region_name().':'.$sample_gene->seq_region_start().'-'.$sample_gene->seq_region_end();

    my @sql = ("INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.location_param', '".$sample_coord."')","INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.location_text', '".$sample_coord."')","INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.gene_param', '".$sample_gene->stable_id()."')","INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.gene_text', 'ensembl_gene')","INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.transcript_param', '".$sample_transcript->stable_id()."')","INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.transcript_text', 'ensembl_transcript')","INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, 'sample.search_text', 'ensembl_gene')");

    my $output_hash = {};
    $output_hash->{'sql_command'} = \@sql;
    $self->dataflow_output_id($output_hash, $self->param('_branch_to_flow_to'));

  }
}

=head2 write_output

 Arg [1]    : None
 Description: Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
              Here we have nothing to do.
 Returntype : None
 Exceptions : None

=cut

sub write_output {
}

1;
