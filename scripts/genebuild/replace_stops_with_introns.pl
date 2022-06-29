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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(clone_Gene empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(replace_stops_with_introns);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(setup_fasta_db);

my $dbname = '';
my $dbuser   = '';
my $dbhost = '';
my $dbport;
my $dbpass;
my $output_dir;
my $fasta_file;
my $write = 0;

my $options = GetOptions ("dbuser=s"     => \$dbuser,
                          "dbhost=s"     => \$dbhost,
                          "dbport=i"     => \$dbport,
                          "dbname=s"     => \$dbname,
                          "dbpass=s"     => \$dbpass,
                          "fasta_file=s" => \$fasta_file,
                          "write!"       => \$write,
                          "output_dir=s" => \$output_dir);

if ($fasta_file) {
  setup_fasta_db;
}
my $log_file = $output_dir."/replace_stops_with_introns_".$dbname.".log";
open(LOG,">".$log_file) or die("could not open '$log_file' for writing");

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dbport,
  -user    => $dbuser,
  -pass    => $dbpass,
  -host    => $dbhost,
  -dbname  => $dbname);

if ($fasta_file) {
  $db->get_SequenceAdaptor->fasta($fasta_file);
}
my $slice_adaptor = $db->get_SliceAdaptor;
my @genes_to_process;
my $count = 0;
foreach my $slice (@{$slice_adaptor->fetch_all('toplevel', undef, undef, 1)}) {
  foreach my $old_gene (@{$slice->get_all_Genes(undef, undef, 1)}) {
    ++$count;
    my $max_stops = 999999;
    my $num_transcripts = 0;
    my $num_transcripts_changed = 0;
    my $transcripts = $old_gene->get_all_Transcripts();
    my $gene = clone_Gene($old_gene, 1);
    say LOG "Gene: ".$old_gene->dbID()." ".$old_gene->stable_id();
    $gene->flush_Transcripts();
    TRANSCRIPT: foreach my $transcript (@$transcripts) {
      my $translation = $transcript->translation;
      if ($translation and $transcript->biotype() ne 'polymorphic_pseudogene' and index($translation->seq(),'*') >= 0) {
        # we replace the stop codons with introns
        say LOG "old translation: ".$translation->seq();
        my $new_transcript = replace_stops_with_introns($transcript,$max_stops);
        say LOG "new translation: ".$new_transcript->translation()->seq();
        if ($new_transcript and index($new_transcript->translate()->seq(), '*') == -1) {
          $gene->add_Transcript($new_transcript);
          if ($old_gene->canonical_transcript == $transcript) {
            $gene->canonical_transcript($new_transcript);
          }
          $num_transcripts++;
          $num_transcripts_changed++;
          say LOG "The transcript ", $new_transcript->dbID(), " ", $new_transcript->stable_id(), " has been changed";
        } elsif ($new_transcript and !$new_transcript->translate()) {
          say LOG "Transcript ", $new_transcript->dbID(), " ", $new_transcript->stable_id(), " (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (", $new_transcript->seq_region_start(), ",", $new_transcript->seq_region_end(), ",", $new_transcript->seq_region_strand(), ",", $new_transcript->seq_region_name(), ") does not translate after replacing a maximum of $max_stops stops. Discarded.";
        } elsif ($new_transcript) {
          say LOG "Transcript ", $new_transcript->dbID(), " ", $new_transcript->stable_id(), " (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (", $new_transcript->seq_region_start(), ",", $new_transcript->seq_region_end(), ",", $new_transcript->seq_region_strand(), ",", $new_transcript->seq_region_name(), ") does not translate after replacing a maximum of $max_stops stops. Discarded!";
        } else {
          say LOG "No transcript defined after replacing stops.";
        }
      } else {
        # we skip the transcript by adding it to the gene as it is
        say LOG "Transcript ", $transcript->dbID(), " ", $transcript->stable_id(), " does not contain any stop or it is a polymorphic_pseudogene. Kept as it is.";
        $gene->add_Transcript($transcript);
        $num_transcripts++;
      }
    }
    
    if ($num_transcripts and $num_transcripts_changed) {
      say LOG "Updating gene: ", $old_gene->dbID(), " ", $old_gene->stable_id();
      empty_Gene($gene);
      push(@genes_to_process, $old_gene, $gene);
    } elsif ($num_transcripts) {
      say LOG "Gene kept: ", $gene->dbID(), " ", $gene->stable_id();
    } else {
      say LOG "Gene does not have any transcript. Not storing.";
    }
    if ($count%1000 == 0) {
      say "Processed $count genes at ".localtime;
    }
  }
}

my $gene_adaptor = $db->get_GeneAdaptor;
foreach my $gene (@genes_to_process) {
  if ($gene->dbID) {
    say LOG "Removing old gene: ", $gene->dbID(), " ", $gene->stable_id();
    $gene_adaptor->remove($gene);
  }
  else {
    $gene_adaptor->store($gene);
    say LOG "Stored new gene: ", $gene->dbID(), " ", $gene->stable_id();
  }
}
close LOG or die("could not close '$log_file'");

1;
