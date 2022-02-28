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

=head1 NAME

  find_and_remove_duplicates.pl

=head1 DESCRIPTION

This script takes database options and finds and removes duplicate genes.

=cut 

use warnings;
use strict;
use feature 'say';


use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception;

my $host;
my $port=3306;
my $dbname;
my $user;
my $pass;
my $dnahost;
my $dnaport;
my $dnadbname;
my $dnauser;
my $biotype;
my $delete_duplicate_transcripts = 0;

GetOptions( 'dbhost:s'                     => \$host,
            'dbport:n'                     => \$port,
            'dbname:s'                     => \$dbname,
            'dbuser:s'                     => \$user,
            'dbpass:s'                     => \$pass,
            'dnadbhost:s'                  => \$dnahost,
            'dnadbport:s'                  => \$dnaport,
            'dnadbuser:s'                  => \$dnauser,
            'dnadbname:s'                  => \$dnadbname,
            'biotype:s'                    => \$biotype,
            'delete_duplicate_transcripts:s' => \$delete_duplicate_transcripts);

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                              -user   => $user,
                                              -port   => $port,
                                              -dbname => $dbname,
                                              -pass   => $pass, );



my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                -port    => $dnaport,
                                                -user    => $dnauser,
                                                -host    => $dnahost,
                                                -dbname  => $dnadbname);

$dba->dnadb($dnadb);

my $gene_adaptor = $dba->get_GeneAdaptor;
my $transcript_adaptor = $dba->get_TranscriptAdaptor;
my $slice_adaptor = $dba->get_SliceAdaptor;
my $slices = $slice_adaptor->fetch_all('toplevel');

foreach my $slice (@{$slices}) {
  my $genes = $gene_adaptor->fetch_all_by_Slice($slice,undef,undef,undef,$biotype);
  my $gene_strings;

  foreach my $gene (@{$genes}) {
    my $gene_string = "";# $gene->start.":".$gene->end.":".$gene->seq_region_name;
    my $transcripts = $gene->get_all_Transcripts();
    my $transcript_strings;
    foreach my $transcript (@{$transcripts}) {
      my $transcript_string = $transcript->biotype.":".$transcript->start.":".$transcript->end.":".$transcript->seq_region_name;
      my $exons = $transcript->get_all_Exons();
      my $exon_string = generate_exon_string($exons);
      $transcript_string .= ":".$exon_string;
      if($transcript->translation) {
        $transcript_string .= ":".$transcript->translation->seq;
      }

      if($transcript_strings->{$transcript_string}) {
        say "Found a duplicate transcript within a gene:";
        say "Duplicate transcript id: ".$transcript->dbID;
        say "Duplicate transcript start: ".$transcript->start;
        say "Duplicate transcript end: ".$transcript->end;
        say "Duplicate strand: ".$transcript->strand;
        say "Duplicate name: ".$transcript->seq_region_name;
	if($delete_duplicate_transcripts) {
          $transcript_adaptor->remove($transcript);
        }
      } else {
        $transcript_strings->{$transcript_string} = 1;
      }

      $gene_string .= $transcript_string
    } # end foreach my $transcript

    if($gene_strings->{$gene_string}) {
      say "Found a duplicate gene:";
      say "Duplicate id: ".$gene->dbID;
      say "Duplicate start: ".$gene->start;
      say "Duplicate end: ".$gene->end;
      say "Duplicate strand: ".$gene->strand;
      say "Duplicate name: ".$gene->seq_region_name;

      eval {
        say "Deleting duplicate: ".$gene->dbID();
        $gene_adaptor->remove($gene);
      }; if($@){
        say "Couldn't remove gene ".$gene->dbID.":\n($@)\n";
      }

    } else {
      $gene_strings->{$gene_string} = 1;
    }

  } # end foreach my $gene

} # end foreach my $slice

exit;


sub generate_exon_string {
  my ($exon_array) = @_;

  my $exon_string = "";
  foreach my $exon (@{$exon_array}) {
    my $start = $exon->start();
    my $end = $exon->end();
    $exon_string .= $start."..".$end.":";
  }

  return($exon_string);
}

