=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Bam2Introns

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::RunnableDB::Bam2Introns->new(
    );

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

Takes an input id of a rough transcript and fetches reads associated with that model
from a BAM file, then runs a spliced exonerate alignment against genomic DNA or the
rough transcript. Writes output as SAM files.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRecoveryIntrons;

use warnings ;
use strict;

use Bio::SeqIO;
use Bio::EnsEMBL::Analysis::Runnable::Bam2Introns;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns');

# fetch input
# get the transcript in question
# get all the reads that overlap
# do we need to hold them in memory or can we stream them into the file?
sub fetch_input {
  my ($self) = @_;

  my $fh = Bio::SeqIO->new(-format => 'fasta', -file => $self->param('iid'));
  my %seq_hash;
  while(my $seq = $fh->next_seq) {
    $seq_hash{$seq->display_id} = $seq;
  }
  $self->param('seq_hash', \%seq_hash);
  # set uo the runnable
  my $program = $self->param('program_file');
  $program = "/software/ensembl/genebuild/bin/exonerate64-0.9.0" unless $program;

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Bam2Introns->new(
     -analysis     => $self->create_analysis,
     -program      => $program,
     -basic_options => $self->get_aligner_options,
     -target_file => $self->param('genome_file'),
     -query_file => $self->param('iid'),
     -percent_id   => $self->param('percent_id'),
     -coverage     => $self->param('coverage'),
     -missmatch     => $self->param('missmatch'),
    );
  $self->runnable($runnable);
}

sub process_features {
    my ($self, $flist) = @_;

# first do all the standard processing, adding a slice and analysis etc.
# unless we are projecting in which case we dont really nead a slice

    my $slice_adaptor = $self->get_database_by_name('reference_db')->get_SliceAdaptor;
    my @dafs;
    foreach my $f (@$flist) {
        my @features;
# get as ungapped features
        foreach my $ugf ( $f->ungapped_features ) {
            my $out_slice = $slice_adaptor->fetch_by_name($ugf->seqname);
            $ugf->slice($out_slice);
        }
        push(@dafs, @{$self->build_dna_align_features($f, \@features)});
    }
    return \@dafs;
}

1;
