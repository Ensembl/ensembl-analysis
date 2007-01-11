
=head1 NAME

Prints - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prints;

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

  # because Prints reports matches off the end of the sequence, 
  # we need the sequence lengths to be able to trim them back
      
  my $command =  $self->program ." " . 
      $self->database . " " . 
      $self->queryfile . " " .
      "-fjR  > " . 
      $self->resultsfile;

  throw("Failed during prints run $!\n") unless 
      system($command) == 0 ;
}


sub parse_results {
  my ($self) = @_;

  my %seqlengths;
  
  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -file   => $self->queryfile);
  while (my $seq = $seqio->next_seq) {
    $seqlengths{$seq->id} = $seq->length;
  }
  $seqio->close;
    
  my ($fh);

  my $resfile = $self->resultsfile;
  
  if (-e $resfile) {    
    if (-z $resfile) {  
      # No hits found
      return; 
    }       
    else {
      open ($fh, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");
    }
  }

  my (@fps, %printsac, $seq_id);
  
  while (<$fh>) {
    my $line = $_;
    chomp $line;
    # Pattern match the Sn; field which should contain the SequenceId and Accession
    
    if ($line =~ s/^Sn;//) { # We have identified a Sn; line so there should be the following:	    
      #ENSP00000003603 Gene:ENSG00000000003 Query:AL035608 Contig:AL035608.00001 Chr:chrX basepair:97227305
      ($seq_id) = $line =~ /^\s*(\w+)/;
    }
        
    if ($line =~ s/^1TBH//) {
      my ($id) = $line =~ /^\s*(\w+)/;
      my ($ac) = $line =~ /(PR\w+);?\s*$/;
      $printsac{$id} = $ac;
    }
    
    if ($line =~ s/^3TB//) {
      if ($line =~ s/^[HN]//) {
        my ($num,$temp1,$tot1);
        # Grab these lines
        #  1433ZETA        1  of  6  88.19   1328    1.00e-16  ELTVEERNLLSVAYKNVIGARRASWRIITS   30   35   36   48

        $line =~ s/^\s+//;        
        my @elements = split /\s+/, $line; 
        
        # Name each of the elements in the array
        my ($fingerprintName,
            $motifNumber,
            $temp,
            $tot,
            $percentageIdentity,
            $profileScore,
            $pvalue,
            $subsequence,
            $motifLength,
            $lowestMotifPosition,
            $matchPosition,
            $highestMotifPosition) = @elements;
        
        # If protein is 10,000+ residues (i.e. titin), then last two elements are merged, e.g.:
        # VEGFRECEPTOR    5  of  6  39.20   406     1.05e-05  LIVRNARKENAGKYTLVL                                      18   374  13564392
        if (!defined $highestMotifPosition) {
          # First five characters of $matchPosition is actual $matchPosition
          $highestMotifPosition = $matchPosition;
          $matchPosition = substr($highestMotifPosition, 0, 5, '');
        }
        
        my $start = $matchPosition;
        my $end = $matchPosition + $motifLength - 1;
        
        # Prints sometimes reports -1 for the motif position start;
        $lowestMotifPosition = 1 if $lowestMotifPosition == -1;
        $highestMotifPosition = $lowestMotifPosition + $motifLength -1 
            if $highestMotifPosition == -1;

        # It's possible that near the ends of the sequence,
        # fragment matches are found. These can result in a 
        # start < 1 or end > seq_length. Make a policy decision
        # prune these back in the sequence, but not the hit
        $start = 1 if $start < 1;
        $end = $seqlengths{$seq_id} if $end > $seqlengths{$seq_id};
    
        my $fp = $self->create_protein_feature($start, $end, 
                                               $profileScore,
                                               $seq_id, 
                                               $lowestMotifPosition, 
                                               $highestMotifPosition, 
                                               $printsac{$fingerprintName},
                                               $self->analysis, 
                                               $pvalue, 
                                               $percentageIdentity);
        push @fps, $fp;
      }
    }
  }
  close($fh);

  $self->output(\@fps);
}


1;

