# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene2x.pm

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


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 APPENDIX



=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene2x;

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
        '-analysis' => $self->analysis,
        '-genes' => $self->genes,
        '-repeat_features' => $self->repeat_blocks,
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
