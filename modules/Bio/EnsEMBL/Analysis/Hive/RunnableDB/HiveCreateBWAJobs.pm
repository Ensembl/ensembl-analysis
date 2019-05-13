=head1 LICENSE
 
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute
 
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateBWAJobs;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');


=head2 fetch_input

 Arg [1]    : None
 Description: Creates input id based on a custom table 'csvfile_table' in the hive database
              It will generate the parameters for BWA based on the data for each file
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
    my $self = shift;

    my $table_adaptor = $self->db->get_NakedTableAdaptor;
    $table_adaptor->table_name($self->param('csvfile_table'));
    my $results = $table_adaptor->fetch_all($self->param('sample_id_column').' = "'.$self->param($self->param('sample_id_column')).'" AND '.$self->param('sample_column').' = "'.$self->param($self->param('sample_column')).'"',
        ($self->param_is_defined('one_per_key') ? $self->param('one_per_key') : undef),
        ($self->param_is_defined('key_list') ? $self->param('key_list') : undef),
        ($self->param_is_defined('value_column') ? $self->param('value_column') : undef));
    my @output_ids;
    my $column_hash = $self->param('column_names');
    foreach my $result (@$results) {
        my @row;
        foreach my $key (@$column_hash) {
            push(@row, $result->{$key});
        }
        my $short_read_aligner_options = '-n '.(int($result->{read_length}/2)).' -i '.$result->{read_length};
        $short_read_aligner_options .= ' -t '. $self->param('use_threading') if ($self->param('use_threading'));
        $short_read_aligner_options .= ' -I ' if ($result->{'is_13plus'});
        push(@row, $short_read_aligner_options);
        push(@output_ids, \@row);
    }
    push(@$column_hash, 'short_read_aligner_options');
    $self->param('column_names', $column_hash);
    $self->param('inputlist', \@output_ids);
}

1;
