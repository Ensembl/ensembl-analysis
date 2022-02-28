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

=cut

  This is a script for adding CCDS models as
  supporting features to a new database.

=cut


use strict;
use warnings;
use Data::Dumper;

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

my $host   = '';
my $user   = 'ensro';
my $dbname = '';
my $port   = 3306;
my $pass   = '';
my $path;

my $file;
my $verbose = 0;
my %transcripts = ();

GetOptions(
  'host|dbhost|h:s'   => \$host,
  'user|dbuser|u:s'   => \$user,
  'dbname|db|D:s' => \$dbname,
  'path:s'   => \$path,
  'port|dbport|P:s'   => \$port,
  'pass|dbpass|p:s'   => \$pass,
  'file:s'   => \$file,
  'verbose!' => \$verbose,
);

open (IN, "<$file") or die "Can't open file $file\n";

while (<IN>){
  chomp;
  my @columns = split(' ');
  $transcripts{$columns[1]} = $columns[0]; # hash having keys = transcript stable IDs and values = CCDS IDs 
}
close IN;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host     => $host,
  -user     => $user,
  -port     => $port,
  -pass     => $pass,
  -dbname   => $dbname,
  -path     => $path
);

my $transcript_adaptor = $db->get_TranscriptAdaptor();
my $tsf_adaptor = $db->get_TranscriptSupportingFeatureAdaptor;

my $undef = 0;
foreach my $transcript_stable_id (keys %transcripts){
  print "Modifying transcript: $transcript_stable_id matching CCDS: ".$transcripts{$transcript_stable_id}."\n";
  my  $transcript =
  $transcript_adaptor->fetch_by_stable_id($transcript_stable_id);

  if (!defined( $transcript ) ) {
    $undef++;
    print STDERR "Can't get $transcript_stable_id\n";
    next;
  }
  if ($verbose) {
    print "Transcript analysis: ", $transcript->analysis->dbID,"\n";
  }

  my @exon_features;
  my @features;

  my @exons = @{$transcript->get_all_translateable_Exons};
  if ($verbose) {
    print "Biotype: " . $transcript->biotype
        . "\tExons: " . scalar(@exons) . ".\n";
  }

  if ($exons[0]->strand == 1){
    @exons =sort {$a->start <=> $b->start} @exons;
  }else{
    @exons =sort {$b->start <=> $a->start} @exons;
  }

  my $exon_start = 1;
  my $exon_length= 0;

  foreach my $exon(@exons){
    if ($verbose) {
      print "Exon: ",$exon->start," end ",$exon->end," strand ",$exon->strand,"\t";
    }

    my $fp = new Bio::EnsEMBL::FeaturePair();
    $fp->start   ($exon->start);
    $fp->end     ($exon->end);
    $fp->strand  ($exon->strand);
    $fp->seqname ($exon->slice->seq_region_name);
    $fp->hseqname($transcripts{$transcript_stable_id});
    $fp->hstart  ($exon_start);
    $exon_length = $exon->end-$exon->start;
    $fp->hend    ($exon_start+$exon_length);
    if ($verbose) {
      print "Feature start: ", $fp->hstart, " end: ", $fp->hend,"\n";
    }
    $exon_start+=$exon_length+1;
    $fp->hstrand (1);

    push (@exon_features, $fp);
  }
  my $tsf = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@exon_features);
  $tsf->seqname($transcript->slice->seq_region_name);
  $tsf->slice($transcript->slice);
  $tsf->score(100);
  $tsf->analysis($transcript->analysis);

  push (@features, $tsf);

  $tsf_adaptor->store($transcript->dbID(), \@features);
}
print "Number of stable_ids that could not be processed: $undef\n";
