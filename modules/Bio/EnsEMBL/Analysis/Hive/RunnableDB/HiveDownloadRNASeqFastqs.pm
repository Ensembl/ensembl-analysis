=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs;

use warnings;
use strict;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub write_output {
  my ($self) = @_;
  my $ftp_base_url = $self->param('ftp_base_url');
  my $fastq = $self->param('iid');
  my $path = $self->param('input_dir');
  my $srr;

  if(-e $path.'/'.$fastq) {
    $self->complete_early('Input file already exists, will not download');
  }

  if ($fastq =~ m/_/){
    $srr = (split /_/, $fastq)[0];
  }
  else{
    $srr = (split /\./, $fastq)[0];
  }
  my $first = substr $srr, 0, 6;
  my $second = '00'.(substr $srr, -1, 1);

  my $res = $self->run_system_command(['wget', '-qq', "$ftp_base_url/$first/$second/$srr/$fastq",  '-P', $path]);
  if ($res) {
    $res >>= 8;
    if ($res == 8) {
      $res = $self->run_system_command(['wget', '-qq', "$ftp_base_url/$first/$srr/$fastq",  '-P', $path]);
      if ($res) {
        $res >>= 8;
        $self->throw("Could not download file $fastq error code is $res");
      }
    }
    elsif ($res) {
      $self->throw("wget died with error code $res");
    }
  }

  unless(-e $path.'/'.$fastq) {
    $self->throw("Did not find the fastq file on the expected path. Path:\n".$path."/".$fastq);
  }

}

1;
