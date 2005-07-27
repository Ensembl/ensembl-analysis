
=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Protein::Hmmpfam

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Pipeline::Runnable::ProteinAnnotation::Pfam->new 
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => '/data/Pfam_ls;/data/Pfam_fs');
  $hmm->run;
  my @results = @{$hmm->output};

=head1 DESCRIPTION

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Pfam;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);



sub multiprotein{
  my ($self) = @_;
  return 1;
}


##################################
sub run_analysis {
  my ($self) = @_;

  my @other_runnables;
  foreach my $db_file (split(/;/, $self->analysis->db_file)) {
    my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam->
        new(-query => $self->query,
            -analysis => $self->analysis,
            -database => $db_file,
            -options  => $self->options);
    $run->run;
    push @other_runnables, $run;
  }
    
  $self->runnables_by_priority(\@other_runnables);
}


#################################
sub parse_results {
  my ($self) = @_;

  # the filter gets rid "secondary tier" hits (e.g. from
  # fs matches) that overlap with "primary tier" matches.
  
  my (%all_hits_by_seqid, @all_hits);

  foreach my $run (@{$self->runnables_by_priority}) {
    my %these_hits_by_seqid;

    foreach my $hit (@{$run->output}) {
      push @{$these_hits_by_seqid{$hit->seqname}}, $hit;
    }

    foreach my $seqid (keys %these_hits_by_seqid) {
      my @kept_hits;

      foreach my $hit (@{$these_hits_by_seqid{$seqid}}) {
        my $overlap = 0;
        if (exists($all_hits_by_seqid{$seqid})) {
          foreach my $ohit (@{$all_hits_by_seqid{$seqid}}) {
            if ($hit->start <= $ohit->end and $hit->end >= $ohit->start) {
              $overlap = 1;
              last;
            }
          }
        }

        if (not $overlap) {
          push @kept_hits, $hit;
        }
      }

      push @{$all_hits_by_seqid{$seqid}}, @kept_hits;
    }
  }

  # finally, strip off those version numbers from the accessions
  foreach my $seqid (keys %all_hits_by_seqid) {
    foreach my $hit (@{$all_hits_by_seqid{$seqid}}) {
      my $acc = $hit->hseqname;
      $acc =~ s/\.\d+$//;
      $hit->hseqname($acc);

      push @all_hits, $hit;
    }
  }

  $self->output(\@all_hits);
}



################################
sub runnables_by_priority {
  my ($self, $runs) = @_;

  if (defined $runs) {
    $self->{_runs} = $runs;
  }

  return $self->{_runs};
}


1;
