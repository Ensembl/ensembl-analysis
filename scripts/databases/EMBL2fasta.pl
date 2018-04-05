#!/usr/env perl
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


use strict;
use warnings;

use Getopt::Long;

use Bio::SeqIO;

use Bio::EnsEMBL::KillList::DBSQL::DBAdaptor;
use Bio::EnsEMBL::KillList::KillObject;
use Bio::EnsEMBL::KillList::Sequence;

my $embl;
my $fasta;
my $do_vertrna = 0;
my $print_seq = 1;
my $print_description = 0;
my $acc_display = 'accession.version';
my $kill_list_host;
my $kill_list_port;
my $kill_list_user;
my $kill_list_pass;
my $kill_list_dbname = 'gb_kill_list';
my $kill_protein = 0;
my $format = 'embl';
my $reason_why = 'ensembl_tr_1';
my $username = $ENV{USER};

&GetOptions (
    'fasta=s'       => \$fasta,
    'embl=s'        => \$embl,
    'format=s'      => \$format,
    'acc_display=s' => \$acc_display,
    'vertrna!'      => \$do_vertrna,
    'desc!'         => \$print_description,
    'kill!'         => \$kill_protein,
    'kill_list_host=s' => \$kill_list_host,
    'kill_list_port=s' => \$kill_list_port,
    'kill_list_user=s' => \$kill_list_user,
    'kill_list_pass=s' => \$kill_list_pass,
    'kill_list_dbname=s' => \$kill_list_dbname,
    'reason_why=s' => \$reason_why,
    'user=s'       => \$username,
    );

if (!$fasta) {
  $fasta = $embl;
  $fasta =~ s/\.\w+$/\.fasta/;
}

my $reason;
my %killlist;
my $ko_adaptor;
my $user;
my $species_adaptor;
my %species;
if ($kill_protein) {
  my $killlist_db = Bio::EnsEMBL::KillList::DBSQL::DBAdaptor->new(
    -host => $kill_list_host,
    -port => $kill_list_port,
    -user => $kill_list_user,
    -pass => $kill_list_pass,
    -dbname => $kill_list_dbname,
  );
  $reason = $killlist_db->get_ReasonAdaptor->fetch_by_why($reason_why);
  die("Could not find reason $reason_why") unless ($reason);
  $ko_adaptor = $killlist_db->get_KillObjectAdaptor;
  $species_adaptor = $killlist_db->get_SpeciesAdaptor;
  $user = $killlist_db->get_UserAdaptor->fetch_by_user_name($username);
  die("could not find user '$username'") unless ($user);
  foreach my $kill_object (@{$ko_adaptor->fetch_all_by_reasonID($reason->dbID)}) {
    $killlist{$kill_object->version} = 1;
  }
}

my $eh = Bio::SeqIO->new(-format => $format, -file => $embl);
my $fh = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$fasta);
$fh->preferred_id_type($acc_display);
while (my $seq = $eh->next_seq) {
  if ($seq->molecule and $seq->molecule eq 'mRNA') {
    if ($do_vertrna) {
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
  }
  elsif (!$seq->molecule and $seq->primary_seq->alphabet eq 'protein') {
    my $annot_collection = $seq->annotation;
    if ($annot_collection) {
      my @objects = $annot_collection->get_Annotations('seq_update');
      my $sequence_version = $objects[0]->value || 1;
      $seq->version($sequence_version); # BioPerl uses the entry accession...
      if ($kill_protein and $seq->namespace eq 'TrEMBL') {
        my @evidence = $annot_collection->get_Annotations('evidence');
        if (@evidence and substr($evidence[0]->value, 0, 1) == 1) {
          foreach my $annotation ($annot_collection->get_Annotations('comment')) {
            if ($annotation->text =~ /CAUTION:.+Ensembl\s+automatic\s+analysis\s+pipeline.+preliminary\s+data/s) {
              next if (exists $killlist{$seq->accession.'.'.$sequence_version});
              my $taxon;
              my $tax_id = $seq->species->ncbi_taxid;
              if (exists $species{$tax_id}) {
                $taxon = $species{$tax_id};
              }
              else {
                $taxon = $species_adaptor->fetch_by_dbID($tax_id);
                die("Could not find taxon id '$tax_id'") unless ($taxon);
                $species{$tax_id} = $taxon;
              }
              my $killed = Bio::EnsEMBL::KillList::KillObject->new(
                  -mol_type => 'protein',
                  -accession => $seq->accession,
                  -version => $seq->accession.'.'.$sequence_version,
                  -description => $seq->description,
                  -external_db_id => 2200,
                  -reasons => [$reason],
                  -taxon => $taxon,
                  -user => $user,
                  -sequence => Bio::EnsEMBL::KillList::Sequence->new(
                      -sequence => $seq->seq,
                    ),
                  );
              $killed->flush_Analyses_allowed;
              $ko_adaptor->store($killed, 0, 1);
            }
          }
        }
      }
    }
  }
  if ($print_seq) {
    $seq->desc(' ') unless ($print_description);
    $fh->write_seq($seq);
  }
}
