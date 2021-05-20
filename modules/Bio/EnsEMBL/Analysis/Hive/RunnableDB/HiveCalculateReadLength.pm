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

use File::Spec::Functions qw(catfile);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my ($self) = @_;

  $self->param_required('read_length_table');
  $self->throw(catfile($self->param('input_dir'), $self->param('iid')).' does not exist')
    unless(-e catfile($self->param_required('input_dir'), $self->param_required('iid')));
}

sub run {
  my ($self) = @_;

  my $read_length = 0;
  my $count = 0;

  open(PH, 'zcat '.catfile($self->param('input_dir'), $self->param('iid')).' |')
    or $self->throw('Could not open '.catfile($self->param('input_dir'), $self->param('iid')));
  while (my $line = <PH>) {
    if ($count%4 == 1){
      $read_length = length($line) if (length($line) > $read_length);
    }
    ++$count;
  }
  close(PH) or $self->throw('Could not close '.catfile($self->param('input_dir'), $self->param('iid')));

  $self->say_with_header(catfile($self->param('input_dir'), $self->param('iid'))."  READ LENGTH: $read_length");
  $self->output([$read_length]);
}


sub write_output {
  my ($self) = @_;

  my $table_adaptor = $self->db->get_NakedTableAdaptor;
  $table_adaptor->table_name($self->param('read_length_table'));

  $table_adaptor->store([{'fastq'       => $self->param('iid'),
                          'read_length' => $self->output->[0]}]);
}

1;
