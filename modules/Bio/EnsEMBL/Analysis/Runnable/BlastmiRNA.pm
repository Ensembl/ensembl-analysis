=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::BlastmiRNA - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::BlastmiRNA->new
     (
      -query    => $slice,
      -program  => 'wublastn',
      -database => 'embl_vertrna',
      -options  => 'hitdist=40 -cpus=1',
      -parser   => $bplitewrapper,
      -filter   => $featurefilter,
     );
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

Modified blast runnable for specific use with miRNAs
Use for running BLASTN of genomic vs miRNAs prior to 
miRNA analysis.
Keeps the coverage in the dna_align_feature score field.
Also clusters overlapping hits and picks the one with the lowest 
evalue to represent that cluster.


=cut


package Bio::EnsEMBL::Analysis::Runnable::BlastmiRNA;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Blast);

=head2 run

  Arg       : string, directory
  Function  : a generic run method. This checks the directory specifed
            : to run it, write the query sequence to file, marks the query sequence
            : file and results file for deletion, runs the analysis parses the 
            : results and deletes any files
  Returntype: 1
  Exceptions: throws if no query sequence is specified
  Example   : 

=cut

sub run{
  my ($self, $dir) = @_;
  $self->workdir($dir) if($dir);
  throw("Can't run ".$self." without a query sequence") 
    unless($self->query);
  $self->checkdir();
  my $filename = $self->write_seq_file();
  $self->files_to_delete($filename);
  $self->files_to_delete($self->resultsfile);
  $self->run_analysis;
  $self->parse_results;
  $self->delete_files;
  return 1;
}

=head2 parse_results

  Arg [1]   : Scalar coverage cutoff value
  Function  : overrides the parse results method in runnable to allow for 
            : coverage calculations
  Returntype: none
  Exceptions: none
  Example   : $self->parse_results(80);

=cut

sub parse_results{
  my ($self,$coverage_cutoff) = @_;
  my $slice = $self->query;
  my $analysis = $self->analysis;
  my $results = $self->results_files;
  my @daf_coverage_results;
  my $filtered_output;
  my $bplite = $self->parser->get_parsers($results);
  foreach my $blast (@{$bplite}){
      while( my $subject = $blast->nextSbjct){
	 while (my $hsp = $subject->nextHSP) {
	   my @daf_results;
	   my $hsp_length = $hsp->length."\t";
	   my $subject_length = $subject->{'LENGTH'};
	   $subject_length = 1 unless ($subject_length);
	   my $coverage = $hsp_length/$subject_length*100;
	   $coverage =~ s/\.\d+//;
	   if ($coverage_cutoff){
	     next unless($coverage > $coverage_cutoff);
	   }
	   $subject->name =~ /^\s*\S+\s+(\w+)/;
	   my $name = $1;
	   push @daf_results, $self->parser->split_hsp($hsp,$name);
	   # add coverage into daf score
	   foreach my $daf(@daf_results){
             #####################################################################
	     # swaps strands over if blast aligns +strand genomic to -ve strand in 
	     # RFAM file. RFAM file sequences are the correct orientation whereas 
             # the genomic can be either, with unspliced DNA it becomes
	     # impossibe to tell automatically
	     if ($daf->hstrand == -1 && $daf->strand == 1){
	       $daf->strand(-1);
	       $daf->hstrand(1);
	     }
	     $daf->score($coverage);
	     $daf->external_db_id(3300);
             $daf->slice($slice);
             $daf->analysis($analysis);
	     push  @daf_coverage_results, $daf;
	   }
	 }
      }
    }
  return undef unless  ( @daf_coverage_results);
  my $output = $self->cluster(\@daf_coverage_results);
  $self->output($output);
}

=head2 cluster

  Arg [1]   : Array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Clusters overlapping blast hits and chooses the one to represent the cluster.
  Returntype: Array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Exceptions: None
  Example   : my $output = $self->cluster(\@dna_align_features);

=cut

sub cluster{
  my ($self,$dafs_ref)=@_;
  my @dafs = @$dafs_ref; 
  @dafs = sort{$a->p_value <=> $b->p_value} @dafs;
  my $start =0;
  my @representative_sequences;
 DAFS: foreach my $daf (@dafs){
    my @cluster;
    $start++;
    next DAFS unless($daf);	
    push @cluster,$daf;
  MATCHES:  for (my $index = $start; $index <= $#dafs ; $index++){
      next MATCHES unless ($dafs[$index]);
      if ($daf->end >= $dafs[$index]->start() && $daf->start() <= $dafs[$index]->end){
	push @cluster,$dafs[$index];
	$dafs[$index] = undef;
      }
    }
    # want to pick identical full length hits by preference
    foreach my $daf ( @cluster ){
      if($daf->score >= 100 && $daf->percent_id == 100){
	push @representative_sequences, $daf;
	next DAFS;
      }
    }
    # otherwise sort by e_value
    @cluster = sort{$b->p_value <=> $a->p_value} @cluster;
    push @representative_sequences, pop @cluster;
  }
  return \@representative_sequences;
}

1;
