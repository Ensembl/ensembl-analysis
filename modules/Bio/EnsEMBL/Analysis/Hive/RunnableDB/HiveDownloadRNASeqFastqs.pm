=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
  if ($fastq =~ m/_/){
    $srr = (split /_/, $fastq)[0];
  }
  else{
    $srr = (split /\./, $fastq)[0];
  }
  my $first = substr $srr, 0, 6;
  my $second = '00'.(substr $srr, -1, 1);

  my $cmd = 'wget '. $ftp_base_url.'/'.$first.'/'.$second.'/'.$srr.'/'.$fastq.' -P '.$path;
  my $cmd2 = 'wget '. $ftp_base_url.'/'.$first.'/'.$srr.'/'.$fastq.' -P '.$path;

  my $trigger = 0;
  my @err = `$cmd 2>&1`;
  foreach my $line (@err){
    if ($line =~ m/No such directory/){
      $trigger = 1;
    }
  }

  if ($trigger){
    my @err2 = `$cmd2 2>&1`;
    foreach my $line2 (@err2){
      if ($line2 =~ m/No such directory/){
        say "fucked mate";
      }
    }
  }

}

1;
