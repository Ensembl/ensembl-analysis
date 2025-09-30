#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStoreUnmappedcDNAs;

use strict;
use warnings;


use Bio::EnsEMBL::UnmappedObject;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub fetch_input {
  my $self = shift;

  $self->create_analysis;
  my $db = $self->get_database_by_name('target_db');
  $self->hrdb_set_con($db, 'target_db');
  return 1;
}

sub run {
  my ($self) = shift;

  return 1;
}

sub write_output {
  my $self = shift;

  my $unmapped_adaptor = $self->hrdb_get_con('target_db')->get_UnmappedObjectAdaptor;
  foreach my $iid (@{$self->param('iid')}) {
    $unmapped_adaptor->store(Bio::EnsEMBL::UnmappedObject->new(
     -type => 'cDNA',
     -identifier => $iid,
     -summary => 'No output from Exonerate',
     -full_desc => 'Exonerate returned no hits using standard parameters plus options --maxintron 400000 and --softmasktarget FALSE',
     -analysis => $self->analysis,
    ));
  }
}


1;
