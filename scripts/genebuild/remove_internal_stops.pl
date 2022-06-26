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

 A script to compare genes and transcripts in a reference database to those in a target database and generate mapping stats. Works under the assumption that the
 mapping process used the Minimap2Remap module and also that pipeline in general as it bases the stats off values in the gene and transcript descriptions, which
 are produced by that module and the mapping pipeline in general. Assuming the data are compatible it will provide stats in terms of mapping percentages for genes
 and transcripts across various biotypes. It also utilises the information about transcript coverage and percent id to provide stats on those

=cut

use warnings;
use strict;
use feature 'say';
use Digest::MD5 qw(md5);

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(align_nucleotide_seqs);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(replace_stops_with_introns);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(clone_Gene empty_Gene);
use Getopt::Long qw(:config no_ignore_case);

my $coord_system = 'toplevel';

my $query_dbname = '';
my $query_user   = '';
my $query_pass   = '';
my $query_host = '';
my $query_port;

my $options = GetOptions ("user=s"              => \$query_user,
                          "host=s"              => \$query_host,
                          "port=i"              => \$query_port,
                          "dbname=s"            => \$query_dbname,
                          "pass=s"            => \$query_pass);

my $query_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $query_port,
  -user    => $query_user,
  -host    => $query_host,
  -dbname  => $query_dbname,
  -pass    => $query_pass);


my $target_gene_adaptor = $query_db->get_GeneAdaptor();

my $target_genes = $target_gene_adaptor->fetch_all();

my $genes_to_remove = [];
my $genes_to_store = [];
foreach my $gene (@$target_genes) {
  unless($gene->get_Biotype->biotype_group() eq 'coding') {
    next;
  }

  if($gene->biotype eq 'polymorphic_pseudogene') {
    next;
  }

  my $updated_transcripts = [];
  my $gene_updated = 0;
  say "Checking gene ".$gene->stable_id();
  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@$transcripts) {
    if($transcript->translation()) {
      my $translation_seq = $transcript->translate->seq;
      my $num_internal_stops = $translation_seq =~ tr/\*/\*/;
      if($num_internal_stops) {
        say "Transcript ".$transcript->stable_id()." has ".$num_internal_stops." internal stops in the translation seq";
        say "  Translation:\n  ".$translation_seq;
        my $no_internal_stop_transcript = replace_stops_with_introns($transcript, $num_internal_stops);
        if($no_internal_stop_transcript) {
          my $transcript_attribs = $transcript->get_all_Attributes();
          $transcript = $no_internal_stop_transcript;
          foreach my $attrib (@$transcript_attribs) {
            $transcript->add_Attributes($attrib);
          }
          say "  Created a new transcript with the internal stops removed, translation seq:\n  ".$transcript->translate->seq();
          $gene_updated = 1;
        } else {
          warning('Could not replace internal stops for '.$transcript->stable_id());
        }
      }
    } # end if($transcript->translation()
    push($updated_transcripts,$transcript);
  } # end foreach my $transcript

  if($gene_updated) {
    my $new_gene = clone_Gene($gene);
    $new_gene->flush_Transcripts();
    foreach my $transcript (@$updated_transcripts) {
      $new_gene->add_Transcript($transcript);
    }
    push(@$genes_to_remove,$gene);
    push(@$genes_to_store,$new_gene);
  }
} # end foreach my $gene (@$target_genes)

unless(scalar(@$genes_to_store) == scalar(@$genes_to_remove)) {
  throw("There is a difference in terms of the number of genes to store (".scalar(@$genes_to_store).") and the number of genes to remove ".scalar(@$genes_to_remove).", this shouldn't happen");
}

say "Storing ".scalar(@$genes_to_store)." updated genes";
foreach my $gene (@$genes_to_store) {
  empty_Gene($gene);
  $target_gene_adaptor->store($gene);
}

say "Removing ".scalar(@$genes_to_remove)." genes";
foreach my $gene (@$genes_to_remove) {
  $target_gene_adaptor->remove($gene);
}


exit;
