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

Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene2x - 

=head1 SYNOPSIS

my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB->new(
						-db => $db_adaptor,
						-input_id => $slice_id,		
					        -analysis => $analysis,
								       );

$runnabledb->fetch_input();
$runnabledb->run();
my @array = @{$runnabledb->output};
$runnabledb->write_output();

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::Pseudogene2x.pm 

Specidfic behaviour for Pseudogene identification in 2x gene-builds


=head1 METHODS


=head1 APPENDIX


=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene2x;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::DBConnection;

use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning 
                                      stack_trace);



use Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB;
use Bio::EnsEMBL::Analysis::Runnable::Pseudogene2x;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB ); 


sub fetch_input {
  my( $self) = @_;

  $self->SUPER::fetch_input;

  my $genes_slice = $self->query;

  my @seq_level_features;
  foreach my $bit (@{$genes_slice->project('seqlevel')}) {
    push @seq_level_features, Bio::EnsEMBL::Feature->new(-start => $bit->from_start + $genes_slice->start - 1,
                                                         -end   => $bit->from_end   + $genes_slice->start - 1,
                                                         -strand => 1);
  }

  foreach my $run (@{$self->runnable}) {
    $run->seqlevel(\@seq_level_features);
  }
}


sub make_runnable {
  my ($self) = @_;

  my $run = Bio::EnsEMBL::Analysis::Runnable::Pseudogene2x->new
      ( 
        -analysis => $self->analysis,
        -genes => $self->genes,
        -repeat_features => $self->repeat_blocks,
        -PS_REPEAT_TYPES              => $self->PS_REPEAT_TYPES,
		    -PS_FRAMESHIFT_INTRON_LENGTH  => $self->PS_FRAMESHIFT_INTRON_LENGTH,
				-PS_MAX_INTRON_LENGTH         => $self->PS_MAX_INTRON_LENGTH,
				-PS_MAX_INTRON_COVERAGE       => $self->PS_MAX_INTRON_COVERAGE,
				-PS_MAX_EXON_COVERAGE         => $self->PS_MAX_EXON_COVERAGE,
				-PS_NUM_FRAMESHIFT_INTRONS    => $self->PS_NUM_FRAMESHIFT_INTRONS,
				-PS_NUM_REAL_INTRONS          => $self->PS_NUM_REAL_INTRONS,
				-SINGLE_EXON                  => $self->SINGLE_EXON,
				-INDETERMINATE                => $self->INDETERMINATE,
				-PS_MIN_EXONS                 => $self->PS_MIN_EXONS,
				-PS_MULTI_EXON_DIR            => $self->PS_MULTI_EXON_DIR,
				-BLESSED_BIOTYPES             => $self->BLESSED_BIOTYPES,
				-PS_PSEUDO_TYPE               => $self->PS_PSEUDO_TYPE,
				-PS_BIOTYPE                   => $self->PS_BIOTYPE,
				-PS_REPEAT_TYPE              => $self->PS_REPEAT_TYPE,
				-DEBUG                        => $self->DEBUG,
																																					
      );
  $self->runnable($run);
}


sub get_all_repeat_blocks {
  my ($self,$repeat_ref) = @_;
  my @repeat_blocks;
  my @repeats = @{$repeat_ref};
  @repeats = sort {$a->start <=> $b->start} @repeats;
  my $curblock = undef;

 REPLOOP: foreach my $repeat (@repeats) {
    my $rc = $repeat->repeat_consensus;
    if ($rc->repeat_class !~ /LINE/ && 
        $rc->repeat_class !~ /LTR/ && 
        $rc->repeat_class !~ /SINE/ &&
        $rc->repeat_class !~ /Unknown/) {  
      next REPLOOP;
    }
    if ($repeat->start <= 0) { 
      $repeat->start(1); 
    }
    if (defined($curblock) && $curblock->end >= $repeat->start) {
      if ($repeat->end > $curblock->end) { 
	$curblock->end($repeat->end); 
      }
    } else {
      $curblock = Bio::EnsEMBL::Feature->new(
                                             -START => $repeat->start,
                                             -END => $repeat->end, 
                                             -STRAND => $repeat->strand
                                             );
      push (@repeat_blocks,$curblock);
    }
  }
    @repeat_blocks = sort {$a->start <=> $b->start} @repeat_blocks;
  return\@repeat_blocks;
}

1;
