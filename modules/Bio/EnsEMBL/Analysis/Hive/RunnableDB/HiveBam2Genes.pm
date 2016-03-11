# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Solexa2Genes

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
my $refine_genes = Bio::EnsEMBL::Analysis::RunnableDB::Solexa2Genes->new (
        -db      => $db,
        -input_id   => $input_id
        -analysis   => $analysis );
$refine_genes->fetch_input();
$refine_genes->run();
$refine_genes->write_output(); #writes to DB


=head1 DESCRIPTION


The databases containing the various features to combine is defined in
Bio::EnsEMBL::Analysis::Config::Databases and the configuration for the
module is defined in Bio::EnsEMBL::Analysis::Config::GeneBuild::Bam2Genes

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes;

use warnings ;
use strict;

use Bio::DB::Sam;
use Bio::EnsEMBL::Analysis::Runnable::Bam2Genes;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

Title        :   fetch_input
Usage        :   $self->fetch_input
Returns      :   nothing
Args         :   none

=cut

sub fetch_input {
    my ($self) = @_;

    my $reference_db = $self->get_database_by_name('dna_db');
    my $slice_adaptor = $reference_db->get_SliceAdaptor;

    my $id = $self->input_id;
    $self->param('slice', $self->fetch_sequence($id, $reference_db));
    my $sam = Bio::DB::Sam->new(
            -bam => $self->param('alignment_bam_file'),
            -expand_flags => 1,
            );
    $self->throw('Bam file ' . $self->param('alignment_bam_file') . "  not found \n") unless $sam;
    $self->runnable(Bio::EnsEMBL::Analysis::Runnable::Bam2Genes->new(
                -analysis => $self->create_analysis,
                -query   => $self->param('slice'),
                -bamfile => $sam,
                -min_length => $self->param('min_length'),
                -min_exons  =>  $self->param('min_exons'),
                -paired => $self->param('paired'),
                -pairing_regex => $self->param('pairing_regex'),
                -max_intron_length => $self->param('max_intron_length'),
                -min_single_exon_length => $self->param('min_single_exon_length'),
                -min_span   => $self->param('min_span'),
                -use_ucsc_naming   => $self->param('wide_use_ucsc_naming'),
                ));
}

sub write_output{
    my ($self) = @_;

    my $outdb = $self->get_database_by_name('output_db');
    my $gene_adaptor = $outdb->get_GeneAdaptor;

    my $fails = 0;
    my $total = 0;

    $gene_adaptor->dbc->disconnect_when_inactive(0);
    foreach my $gene ( @{$self->output} ) {
        eval {
            $gene_adaptor->store($gene);
        };
        if ($@){
            $self->warning("Unable to store gene!!\n$@");
            print STDERR $gene->start ." " . $gene->end ."\n";
            $fails++;
        }
        $total++;
    }
    $self->throw("Not all genes could be written successfully ($fails fails out of $total)") if ($fails);
    $gene_adaptor->dbc->disconnect_when_inactive(1);
    print STDERR "$total genes written after filtering\n";
}

1;
