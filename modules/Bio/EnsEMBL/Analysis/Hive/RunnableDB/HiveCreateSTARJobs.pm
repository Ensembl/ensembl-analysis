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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateSTARJobs;

use strict;
use warnings;

use File::Spec;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');


=head2 param_defaults

 Description: It allows the definition of default parameters for all inherting module.
              These are the default values:
               bunzip => 'bzcat',
               gunzip => 'zcat',
 Returntype : Hashref, containing all default parameters
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
      %{$self->SUPER::param_defaults},
      bunzip => 'bzcat',
      gunzip => 'zcat',
  }
}
=head2 fetch_input

 Arg [1]    : None
 Description: Creates input id based on a custom table 'csvfile_table' in the hive database
              It will generate the parameters for STAR based on the data for each file
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : Throws if there is too many read mates in the set

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
  my %mates;
  my $column_hash = $self->param('column_names');
  my @row;
  foreach my $result (@$results) {
    my $sample_id_column = $result->{$self->param('sample_id_column')};
    my $short_read_aligner_options = $self->param_is_defined('base_options') ? $self->param('base_options') : '';
    $short_read_aligner_options .= ' --runThreadN '.($self->param_is_defined('use_threading') ? $self->param('use_threading') : 1);
    $short_read_aligner_options .= ' --outQSconversionAdd -31 ' if ($result->{'is_13plus'});
    my $decompress;
    if ($result->{filename} =~ /gz$/) {
      $decompress = $self->param('gunzip');
    }
    elsif ($result->{filename} =~ /bz2$/) {
      $decompress = $self->param('bunzip');
    }
    if ($result->{is_mate_1}) {
      foreach my $key (@$column_hash) {
        push(@row, $result->{$key});
      }
# if I'm mate 1 I can store the data in @row and I store a ref to @row in %mates if I have not seen
# mate 2 yet. Otherwise I can add the filename of mate 2 as fastqpair
      if (exists $mates{$sample_id_column}) {
        $self->throw('You have more read mates than expected') if (ref($mates{$sample_id_column}) eq 'ARRAY');
        push(@row, $mates{$sample_id_column}, File::Spec->catfile($self->param('wide_output_dir'), $sample_id_column), $decompress, $short_read_aligner_options);
        push(@output_ids, \@row);
      }
      else {
        $mates{$sample_id_column} = \@row;
      }
    }
    else {
# If I'm mate 2, if I've already seen mate 1 I can add the filename to the arrayref in %mates. Otherwise
# I store filename in %mates
      if (exists $mates{$sample_id_column}) {
        $self->throw('You have more read mates than expected') if ($result->{is_paired} == 0 or ref($mates{$sample_id_column}) ne 'ARRAY');
        push(@{$mates{$sample_id_column}}, $result->{filename}, File::Spec->catfile($self->param('wide_output_dir'), $sample_id_column.'_'), $decompress, $short_read_aligner_options);
        push(@output_ids, $mates{$sample_id_column});
      }
      else {
        $mates{$sample_id_column} = $result->{filename};
      }
    }
  }
  push(@$column_hash, 'fastqpair', 'output_dir', 'decompress', 'short_read_aligner_options');
  $self->param('column_names', $column_hash);
  $self->param('inputlist', \@output_ids);
}

1;
