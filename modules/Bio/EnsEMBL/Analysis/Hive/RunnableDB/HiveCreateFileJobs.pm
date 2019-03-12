#!/usr/bin/env perl

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateFileJobs;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub fetch_input {
    my $self = shift;

    my $list_cmd = "ls ".$self->param('file_path')."/".$self->param('file_regex');
    my @files_list = `$list_cmd`;

    my @output_ids;

    foreach my $file (@files_list){
     chomp $file;
     $file = (split '/', $file)[-1];

     if ($self->param('remove_regex')){
       $file =~ s/$self->param('remove_regex')//;
     }

     my @params;
     push(@params, $file);
     push(@params, $self->param('additional'));

     push (@output_ids, \@params);
    }
    my $column_hash = $self->param('column_names');
    $self->param('column_names', $column_hash);
    $self->param('inputlist', \@output_ids);
}

1;
