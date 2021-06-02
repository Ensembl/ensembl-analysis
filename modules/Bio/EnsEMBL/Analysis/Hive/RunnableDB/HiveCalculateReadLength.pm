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

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCalculateReadLength;

use warnings;
use strict;

use File::Basename;
use File::Spec::Functions qw(catfile file_name_is_absolute);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my ($self) = @_;

  $self->param_required('read_length_table');
  my $filename = $self->input_id;
  if (!file_name_is_absolute($self->input_id)) {
    $filename = catfile($self->param('input_dir'), $self->input_id);
  }

  $self->throw("$filename does not exist")
    unless(-e $filename);
  $self->param('absolute_filepath', $filename);
}

sub run {
  my ($self) = @_;

  my $read_length = 0;
  my $count = 0;

  open(PH, 'zcat '.$self->param('absolute_filepath').' |')
    or $self->throw('Could not open '.$self->param('absolute_filepath'));
  while (my $line = <PH>) {
    if ($count%4 == 1){
      $read_length = length($line) if (length($line) > $read_length);
    }
    ++$count;
  }
  close(PH) or $self->throw('Could not close '.$self->param('absolute_filepath'));

  $self->say_with_header($self->param('absolute_filepath')."  READ LENGTH: $read_length");
  $self->output([$read_length]);
}


sub write_output {
  my ($self) = @_;

  my $table_adaptor = $self->db->get_NakedTableAdaptor;
  $table_adaptor->table_name($self->param('read_length_table'));

  my $filename = $self->input_id;
  if (file_name_is_absolute($filename)) {
    $filename = basename($self->input_id);
  }

  $table_adaptor->store([{'fastq'       => $filename,
                          'read_length' => $self->output->[0]}]);
}

1;
