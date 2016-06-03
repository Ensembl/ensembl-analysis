# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Runnable::Bam2Genes;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(convert_to_ucsc_name);

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

    my $slice = $self->fetch_sequence($self->input_id, $reference_db);
    my $sam = Bio::DB::Sam->new(
            -bam => $self->param('alignment_bam_file'),
            -expand_flags => 1,
            );
    $self->throw('Bam file ' . $self->param('alignment_bam_file') . "  not found \n") unless $sam;
    $self->create_analysis;
    my ($exon_clusters, $cluster_link, $full_slice) = $self->exon_cluster($slice, $sam);
    $self->runnable(Bio::EnsEMBL::Analysis::Runnable::Bam2Genes->new(
                -analysis => $self->analysis,
                -query   => $full_slice,
                -exon_clusters => $exon_clusters,
                -cluster_data => $cluster_link,
                -min_length => $self->param('min_length'),
                -min_exons  =>  $self->param('min_exons'),
                -paired => $self->param('paired'),
                -max_intron_length => $self->param('max_intron_length'),
                -min_single_exon_length => $self->param('min_single_exon_length'),
                -min_span   => $self->param('min_span'),
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


=head2 exon_cluster
  Arg[0]      : Bio::EnsEMBL::Slice object representing the input_id
  Arg[1]      : Bio::DB::Sam object containing the aligned short reads
  Usage       : $self->exon_cluster($slice, $bam)
  Returns     : Array ref Bio::EnsEMBL::Exon, a hash ref representing the links between exons, a Bio::EnsEMBL::Slice
                object representing the full region covered by the input_id
  Description : clusters individual reads into blocks representing exons
                uses pair information to link blocks into transcripts
                filters out poorly supported blocks
=cut

sub exon_cluster {
    my ($self, $slice, $bam) = @_;
    print STDERR "CLUSTER EXON\n";
    my $seq_region_name = $self->param('wide_use_ucsc_naming') ? convert_to_ucsc_name($slice->seq_region_name, $slice) : $slice->seq_region_name;
    my $region = $seq_region_name.':'.$slice->start.'-'.$slice->end;
    # BWA has been run on whole genome. If the slice is not starting at 1, the Core API
    # will shift the coordinate of our clusters which is wrong
    my $full_slice = $slice->start == 1 ? $slice : $slice->seq_region_Slice;
    my %exon_clusters;
    my @exon_clusters;
    my $cluster_data;
    my $cluster_count = 0;
    my $read_count = 0;
    # I can't give parameters to $bam->fetch() so it's easier to create the callback
    # inside the method. Maybe with the low level method it's better.
    my $_process_reads = sub {
        my $read = shift;
        ++$read_count;
        # It seems we always get the unmmapped mate, so we need to remove it
        return if ($read->get_tag_values('XS') or $read->get_tag_values('UNMAPPED'));
        my $name = $read->query->name;
        my $start  = $read->start;
        my $end    = $read->end;
        my $hstart = $read->query->start;
        my $hend   = $read->query->end;
        my $paired = $read->get_tag_values('MAP_PAIR');
        # ignore spliced reads
        # make exon clusters and store the names of the reads and associated cluster number
        for (my $index = @exon_clusters; $index > 0; $index--) {
            my $exon_cluster = $exon_clusters[$index-1];
            if ( $start <= $exon_cluster->end+1 &&  $end >= $exon_cluster->start-1 ) {
                # Expand the exon_cluster
                $exon_cluster->start($start) if $start < $exon_cluster->start;
                $exon_cluster->end($end)     if $end   > $exon_cluster->end;
                $exon_cluster->score($exon_cluster->score + 1);
                # only store the connection data if it is paired in mapping
                if ($paired) {
                    if (exists $cluster_data->{$name}->{$exon_cluster->hseqname}) {
                        $cluster_data->{$name}->{$exon_cluster->hseqname}++;
                        delete $cluster_data->{$name} if ($cluster_data->{$name}->{$exon_cluster->hseqname} == 2);
                    }
                    else {
                        $cluster_data->{$name}->{$exon_cluster->hseqname} = 1;
                    }
                }
                # only allow it to be a part of a single cluster
                return;
            }
        }
        # start a new cluster if there is no overlap
        $cluster_count++;
        # make a feature representing the cluster
        my $feat = Bio::EnsEMBL::FeaturePair->new
            (
             -start      => $start,
             -end        => $end,
             -strand     => -1,
             -slice      => $full_slice,
             -hstart     => $hstart,
             -hend       => $hend,
             -hstrand    => 1,
             -score      => 1,
             -percent_id => 100,
             -hseqname   => "C". $cluster_count,
             -analysis   => $self->analysis,
            );
        # store the clusters in a hash with a unique identifier
        push(@exon_clusters, $feat);
        # store the key within the feature
        if ($paired) {
            if (exists $cluster_data->{$name}->{$feat->hseqname}) {
                $cluster_data->{$name}->{$feat->hseqname}++;
                delete $cluster_data->{$name} if ($cluster_data->{$name}->{$feat->hseqname} == 2);
            }
            else {
                $cluster_data->{$name}->{$feat->hseqname} = 1;
            }
        }
    };
    $bam->fetch($region, $_process_reads);
    print STDERR "Processed $read_count reads\n";
    %exon_clusters = map {$_->hseqname => $_} @exon_clusters;
    return \%exon_clusters, $cluster_data, $full_slice;
}

1;
