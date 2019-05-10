#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::IO::Parser::GFF3;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;

# Connection to the target DB
my $host   = '';
my $port   = '3306';
my $user   = 'ensadmin';
my $pass   = '';
my $dbname = '';
my $file;
my $write;
my $logic_name = 'ncrna';

my $transcript_type = 'miRNA_primary_transcript';
my $coord_system = 'toplevel';
my $biotype = 'miRNA';
my $source = 'mirbase';
my $batch_size;

&GetOptions (
            'host=s'       => \$host,
            'port=s'       => \$port,
            'user=s'       => \$user,
            'pass=s'       => \$pass,
            'dbname=s'     => \$dbname,
            'file=s'       => \$file,
            'batch_size=s' => \$batch_size,
            'source=s'     => \$source,
            'biotype=s'    => \$biotype,
            'logic_name=s' => \$logic_name,
            'write!'       => \$write,
        );

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host => $host,
    -user => $user,
    -port => $port,
    -pass => $pass,
    -dbname => $dbname
    );

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
if (!defined $analysis) {
    $analysis = Bio::EnsEMBL::Analysis->new(
        -logic_name => $logic_name,
        );
}

my $slice_adaptor = $db->get_SliceAdaptor;
my $gff_file = Bio::EnsEMBL::IO::Parser::GFF3->open($file);
my @genes;
my $gene_count = 0;
while($gff_file->next) {
    next unless ($gff_file->get_type eq $transcript_type);
    my $name = $gff_file->get_seqname;
    $name =~ s/^chr//;
    my $slice = $slice_adaptor->fetch_by_region($coord_system, $name);
    my $exon = Bio::EnsEMBL::Exon->new(
        -start => $gff_file->get_start,
        -end => $gff_file->get_end,
        -strand => $gff_file->get_strand,
        -slice => $slice,
        -phase => -1, # It should always be -1
        -end_phase => -1, # It should always be -1
        );
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->add_Exon($exon);
    $transcript->analysis($analysis);
    $transcript->biotype($biotype);
    $transcript->source($source);
    $transcript->stable_id($gff_file->get_attribute_by_name('Name'));
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->add_Transcript($transcript);
    $gene->analysis($analysis);
    $gene->biotype($biotype);
    $gene->source($source);
    $gene->stable_id($gff_file->get_attribute_by_name('ID'));
    push(@genes, $gene);
    $gene_count++;
    if (defined $batch_size and $gene_count == $batch_size) {
        store_genes(\@genes);
        undef @genes;
        $batch_size = 0;
    }
}

store_genes(\@genes);

sub store_genes {
    my $genes = shift;

    my $gene_adaptor = $db->get_GeneAdaptor;
    foreach my $gene (@$genes) {
        $gene_adaptor->store($gene);
    }
}
