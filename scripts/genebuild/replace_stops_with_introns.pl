=head1 LICENSE

 Copyright [2022] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 DESCRIPTION

 This script replaces the stops with introns in non-'polymorphic_pseudogene' transcripts whose translations contain stops in a given database.

=cut

use warnings;
use strict;
use feature 'say';
use Digest::MD5 qw(md5);

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(clone_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(replace_stops_with_introns);

my $dbname = '';
my $dbuser   = '';
my $dbhost = '';
my $dbport;
my $dbpass;

my $options = GetOptions ("dbuser=s" => \$dbuser,
                          "dbhost=s" => \$dbhost,
                          "dbport=i" => \$dbport,
                          "dbname=s" => \$dbname,
                          "dbpass=s" => \$dbpass);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dbport,
  -user    => $dbuser,
  -pass    => $dbpass,
  -host    => $dbhost,
  -dbname  => $dbname);

my $gene_adaptor = $db->get_GeneAdaptor();

foreach my $old_gene (@{$gene_adaptor->fetch_all()}) {
  my $max_stops = 999;
  my $num_transcripts = 0;
  my $num_transcripts_changed = 0;
  my $transcripts = $old_gene->get_all_Transcripts();
  my $gene = clone_Gene($old_gene);
  print STDERR "Gene: ".$old_gene->dbID()." ".$old_gene->stable_id()."\n";
  $gene->flush_Transcripts();
  TRANSCRIPT: foreach my $transcript (@$transcripts) {
    if ($transcript->biotype() ne 'polymorphic_pseudogene' and
        $transcript->translation() and
        $transcript->translation()->seq() =~ /\*/) {
      # we replace the stop codons with introns
      print STDERR "old translation: ".$transcript->translation()->seq()."\n";
      $transcript = replace_stops_with_introns($transcript,$max_stops);
      print STDERR "new translation: ".$transcript->translation()->seq()."\n";
      if ($transcript and $transcript->translate()->seq() !~ /\*/) {
        $gene->add_Transcript($transcript);
	$num_transcripts++;
	$num_transcripts_changed++;
        print STDERR "The transcript ".$transcript->dbID()." ".$transcript->stable_id()." has been changed.\n";
      } elsif ($transcript and !$transcript->translate()) {
        print STDERR "Transcript ".$transcript->dbID()." ".$transcript->stable_id()." (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().",".$transcript->seq_region_end().",".$transcript->seq_region_strand().",".$transcript->seq_region_name().") does not translate after replacing a maximum of $max_stops stops. Discarded.\n";
      } elsif ($transcript) {
        print STDERR "Transcript ".$transcript->dbID()." ".$transcript->stable_id()." (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().",".$transcript->seq_region_end().",".$transcript->seq_region_strand().",".$transcript->seq_region_name().") does not translate after replacing a maximum of $max_stops stops. Discarded.\n";
      } else {
        print STDERR "No transcript defined after replacing stops.\n";
      }
    } else {
      # we skip the transcript by adding it to the gene as it is
      print STDERR "Transcript ".$transcript->dbID()." ".$transcript->stable_id()." does not contain any stop or it is a polymorphic_pseudogene. Kept as it is.\n";
      $gene->add_Transcript($transcript);
      $num_transcripts++;
    }
  }
  
  if ($num_transcripts and $num_transcripts_changed) {
    print STDERR "Removing old gene: ".$old_gene->dbID()." ".$old_gene->stable_id()."\n";
    $gene_adaptor->remove($old_gene);
    $gene->dbID(undef);
    $gene_adaptor->store($gene);
    print STDERR "Stored new gene: ".$gene->dbID()." ".$gene->stable_id()."\n\n";
  } elsif ($num_transcripts) {
    print STDERR "Gene kept: ".$gene->dbID()." ".$gene->stable_id()."\n\n";
  } else {
    print STDERR "Gene does not have any transcript. Not storing.\n\n";
  }
}

