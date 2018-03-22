#!/usr/env perl
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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


use strict;
use warnings;

use Getopt::Long;

use Bio::SeqIO;


# Connection to the target DB
my $embl;
my $fasta;
my $do_vertrna = 0;
my $print_seq = 1;
my $print_description = 0;
my $acc_display = 'accession.version';

&GetOptions (
    'fasta=s'       => \$fasta,
    'embl=s'        => \$embl,
    'acc_display=s' => \$acc_display,
    'vertrna!'      => \$do_vertrna,
    'desc!'         => \$print_description,
    );

my $eh = new Bio::SeqIO(-format => 'embl', -file => $embl);
my $fh = new Bio::SeqIO(-format => 'fasta', -file => '>'.$fasta);
$fh->preferred_id_type($acc_display);
while (my $seq = $eh->next_seq) {
  if ($do_vertrna) {
    next if ($seq->molecule ne 'mRNA');
    $print_seq = 0;
    foreach my $feat ($seq->get_SeqFeatures) {
      if ($feat->primary_tag eq 'source') {
        $print_seq = 1 if (grep('mRNA', $feat->get_tag_values('mol_type')));
        if ($seq->division eq 'SYN') {
          $print_seq = 0 if ($seq->species->db_handle->get_taxonids('Craniata') == 0);
        }
        last;
      }
    }
  }
  if ($print_seq) {
    $seq->desc(' ') unless ($print_description);
    $fh->write_seq($seq);
  }
}
