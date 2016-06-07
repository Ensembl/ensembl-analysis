# Copyright [1999-2016] the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveIndexGenome;

use strict;
use warnings;

use File::Spec;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

 Arg [1]    : None
 Description: Create the command to execute genome indexing with STAR
 Returntype : None
 Exceptions : Throw if splitpath does not return an existing directory name

=cut

sub fetch_input {
  my ($self) = @_;

  my (undef, $dirname, $file) = File::Spec->splitpath($self->param('wide_genome_file'));
  $self->throw("File::Spec->splitpath failed, $dirname does not exist") unless (-e $dirname);
  if (-e "$dirname/SA") {
    $self->complete_early($self->param('wide_genome_file').'is already indexed!');
  }
  else {
    my @command = ($self->param('wide_short_read_aligner'), '--runMode genomeGenerate');
    push(@command, '--runThreadN', $self->param('use_threading'))
      if ($self->param_is_defined('use_threading') and $self->param('use_threading') > 0);
    push(@command, '--genomeDir', $dirname);
    push(@command, '--genomeFastaFiles', $file);
    push(@command, '--sjdbGTFfile', $self->param('annotation_gtf'))
      if ($self->param_is_defined('annotation_gtf'));
    push(@command, '--sjdbOverhang', $self->param('read_length')-1);
      if ($self->param_is_defined('read_length') and $self->param('read_length') > 1);
    push(@command, $self->param('extra_options'))
      if ($self->param_is_defined('extra_options'));
    $self->param('commandline', \@command);
  }
}


=head2 run

 Arg [1]    : None
 Description: Run the STAR command, it will generate the indexes
 Returntype : None
 Exceptions : Throws if STAR fails

=cut

sub run {
  my ($self) = @_;

  $self->throw('Could not execute: '.join(' ', @{$self->param('commandline')}))
    if (system(@{$self->param('commandline')}));
}


=head2 write_output

 Arg [1]    : None
 Description: Return 1 to override SUPER method
 Returntype : Integer 1
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  return 1;
}

1;
