#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFastqDownloadJobs;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

=head2 fetch_input

 Arg [1]    : 
 Description: 

 Returntype : None
 Exceptions : None

=cut


sub write_output {
    my $self = shift;
    my $inputfile = $self->param('inputfile');
    my @fastq_list = `cut -d\$'\t' -f4 $inputfile`;
    my @output_ids;
    foreach my $fastq (@fastq_list){
      chomp $fastq;
      if ($fastq ne ""){
	push(@output_ids, {iid => $fastq})
      }
    }
  $self->dataflow_output_id(\@output_ids, $self->param('fan_branch_code'));
  }

1;
