# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Sam2Bam




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::Sam2Bam->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses samtools to convert SAM files from a directory into a
single sorted indexed merged BAM file

=head1 CONTACT

Post general queries to B<http://lists.ensembl.org/mailman/listinfo/dev>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSam2Bam;

use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::Runnable::Sam2Bam;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

 Arg [1]    : None
 Description: It will create a BAM file from the SAM files given as input via 'filename' using samtools
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;

  if (@{$self->param_required('filename')} > 0) {
    my $program = $self->param('samtools');
    $self->throw("Samtools program not defined in analysis \n") unless (defined $program);
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Sam2Bam->new
      (
       -analysis => $self->create_analysis,
       -header   => $self->param('headerfile'),
       -program  => $program,
       -samfiles => $self->param('filename'),
       -bamfile  => $self->param('intron_bam_file'),
       -genome   => $self->param('genome_file'),
      );
    $self->runnable($runnable);
  }
  else {
    $self->input_id->autoflow(0);
    $self->complete_early('You do not have any SAM files');
  }
}


=head2 run

  Arg [1]   : None
  Function  : Overrides run as we want to disconnect from the database as the job will last at least 30 minutes
  Returntype: None
  Exceptions: None

=cut

sub run {
  my ($self) = @_;
  $self->throw("Can't run - no runnable objects") unless ( $self->runnable );
  $self->dbc->disconnect_if_idle() if ($self->param('disconnect_jobs'));
  my ($runnable) = @{$self->runnable};
  $runnable->run;
}


=head2 write_output

  Arg [1]   : None
  Function  : It does nothing as the BAM file is written and checked in the run method
  Returntype: Integer, 1
  Exceptions: None

=cut

# override write output as we have nothing for the db
sub write_output {
    my ($self) = @_;

    return 1;
}

1;
