=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

# Maybe it's not our copyright...
#
#
#
=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::ProteinAnnotation::PIRSF

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::Runnable::ProteinAnnotation::PIRSF->
    new ( 
    -query      => $query,
    -analysis   => $analysis,
    -database   => $db,
    -datfile    => $datfile,
    );
    

=head1 DESCRIPTION

=head1 CONTACT

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);

###################
# analysis methods
###################

=head2 run_program

 Title    : run_program
 Usage    : $self->program
 Function : makes the system call to program
 Example  :
 Returns  : 
 Args     :
 Throws   :

=cut


sub run_analysis {
    my ($self) = @_;

    # run program
    print STDERR "running ".$self->program." against ".$self->database."\n";
    print STDERR "FILENAME: ".$self->queryfile."\n";

    $self->resultsfile() ;
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
                                -file   => $self->queryfile);
    while (my $seq = $seqio->next_seq) {
      my $filename = $self->create_filename($seq->id, 'fa') ;
      $filename = $self->write_seq_file($seq, $filename) ;
      my $cmd = $self->program .' '.
  	        $self->analysis->parameters .' '.
  	        '-i ' . $filename.' '.
  		'>> ' . $self->resultsfile;
      print STDERR "$cmd\n";   
      $self->throw ("Error running PIRSF ".$self->program." on ".$self->queryfile) 
       unless ((system ($cmd)) == 0);
    }
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

    my $filehandle;
    my $resfile = $self->resultsfile();
    my @fps;
		

    if (-e $resfile) {
      if (-z $resfile) {  
	print STDERR "PIRSF didn't find any hits\n";
	return; 
      } else {
	open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
      }
    }


 # Example output:
#Query sequence: R186.4  matches PIRSF004870: Molybdenum cofactor molybdenum incorporation protein MoeA
#Full-length Match:
#Model      Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
#
#PIRSF004870   1/1       9   390 ..     1   443 []   298.3    2e-87
#

   my $id;
    while (<CPGOUT>) {
      chomp;
      if (/^Query sequence:\s+(\S+)/) {
	$id = $1;
      }

      if (my ($hid,
              $start,
              $end,
              $hstart,
              $hend,
              $score,
              $evalue) = /^(\S+)\s+\d+\/\d+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {

        $evalue = sprintf ("%.3e", $evalue);
	my $percentIdentity = 0;

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
	push @fps, $fp;
      }
    }
    close (CPGOUT); 
    $self->output(\@fps);
}

1;
