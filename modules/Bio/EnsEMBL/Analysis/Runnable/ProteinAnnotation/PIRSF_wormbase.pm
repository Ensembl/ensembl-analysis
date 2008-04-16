package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF_wormbase;

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
 
    my $cmd = $self->program .' '.
	        $self->analysis->parameters .' '.
	        '-i ' . $self->queryfile.' '.
		'> ' . $self->resultsfile;
    print STDERR "$cmd\n";   
    $self->throw ("Error running PIRSF_wormbase ".$self->program." on ".$self->queryfile) 
     unless ((system ($cmd)) == 0);
    
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

      print "$_\n";
      if (my ($hid, $start, $end, $hstart, $hend, $score, $evalue) = /(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {

	print "matched\n";
	my $percentIdentity = 0;

	my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);
	push @fps, $fp;
      }
    }
    close (CPGOUT); 
    $self->output(\@fps);
}

1;
