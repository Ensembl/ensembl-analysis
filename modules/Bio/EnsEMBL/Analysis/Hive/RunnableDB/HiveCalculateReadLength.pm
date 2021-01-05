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
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 write_output

 Arg [1]    : None
 Description: Calculate the read length of a FASTQ file and insert the values
              in the table 'read_length_table' of the Hive DB
 Returntype : None
 Exceptions : Throws if the file does not exist
              Throws if the read length is 0

=cut

sub write_output {
  my ($self) = @_;
  my $fastq = $self->param('iid');
  my $path = $self->param('input_dir');
  my @output_ids;

  if(!-e $path.'/'.$fastq) {
    $self->throw("$path/$fastq does not exist");
  }

  my $read_length_cmd="zcat $path/$fastq| awk \'{if(NR%4==2) print length(\$1)}\' | sort -n | uniq -c";
  my @read_lengths=`$read_length_cmd`;

  my $max = 0;
  my $read_length=0;

  foreach my $length_line (@read_lengths){
    $length_line =~ s/^\s+//;
    my @line_array = split(/ /, $length_line);
    if ($line_array[0] > $max){
      $max = $line_array[0];
      $read_length = $line_array[1];
    }
  }
  $self->throw('Read length is 0, something went wrong') unless ($read_length);

  say $fastq."  READ LENGTH: ".$read_length;

  my $table_adaptor = $self->db->get_NakedTableAdaptor;
  $table_adaptor->table_name($self->param('read_length_table'));

  my $insert_row = [{'fastq'         => $fastq,
                      'read_length'  => $read_length}];

  $table_adaptor->store($insert_row);
}
1;
