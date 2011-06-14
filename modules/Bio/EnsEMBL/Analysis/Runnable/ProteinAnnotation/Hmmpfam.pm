
=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Protein::Hmmpfam

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam->new 
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => '/data/Pfam_ls');

  $hmm->run;
  my @results = @{$hmm->output};

=head1 DESCRIPTION

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam;

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;


@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);


sub multiprotein{
  my ($self) = @_;
  return 1;
}



sub run_analysis {
  my ($self) = @_;
  
  my $options = "";
  if (defined($self->parameters)) {
    $options = $self->parameters ;
  }
  if (defined($self->options)) {
    $options .= $self->options;
  }
  if ($options !~ /\-\-acc/) {
    $options .= ' --acc';
  }
  if ($options !~ /\-\-cpu/) {
    $options .= ' --cpu 1';
  }

  my $cmd = $self->program 
      . ' ' . $options
      . ' ' . $self->database 
      . ' ' . $self->queryfile 
      .' > '. $self->resultsfile;

  print STDERR "Running:$cmd\n";
  
  throw ("Error running ".$self->program." on ".
         $self->queryfile." against ".$self->database)unless 
         ((system ($cmd)) == 0);
}


=head2 parse_results

 Title    :  parse_results
 Usage    :  $self->parse_results ($filename)
 Function :  parses program output to give a set of features
 Example  :
 Returns  : 
 Args     : filename (optional, can be filename, filehandle or pipe, not implemented)
 Throws   :

=cut

sub parse_results {
  my ($self) = @_;

  my (@hits, $id, $hid);

  my $f = $self->resultsfile;

  if (-e $f and not -z $f) {
    my $fh;
    open ($fh, "<$f") or throw ("Error opening $f");
    
    while (<$fh>) {
      chomp;
      #last if /^Alignments of top-scoring domains/;
      #next if (/^Model/ || /^\-/ || /^$/);
      if (/^Query:\s+(\S+)/) {
        $id = $1;
        next;
      }
      if (/^>> (\S+)/) {
        $hid = $1 ;
      }

      if (my ($score,
              $evalue,
              $hstart,
              $hend,
              $start,
              $end) = /^\s+\d+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)/) {

        my $percentIdentity = 0;
        $evalue = sprintf ("%.3e", $evalue);


        my $fp = $self->create_protein_feature($start, 
                                               $end, 
                                               $score, 
                                               $id, 
                                               $hstart, 
                                               $hend, 
                                               $hid,
                                               $self->analysis, 
                                               $evalue,
                                               $percentIdentity);
        push @hits, $fp;
      }
    }
  }
  $self->output(\@hits);
}

1;
