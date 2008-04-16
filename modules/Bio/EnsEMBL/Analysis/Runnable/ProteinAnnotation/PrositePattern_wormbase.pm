package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositePattern_wormbase;

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
	        $self->queryfile.' '.
		'> ' . $self->resultsfile;
    print STDERR "$cmd\n";   
    $self->throw ("Error running PrositePattern_wormbase ".$self->program." on ".$self->filename) 
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
	print STDERR "PrositePattern didn't find any hits\n";
	return; 
      } else {
	open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
      }
    }


 # Example output:
#CE20125        PS00031 NUCLEAR_RECEPTOR        30      56      T

   my $id;
    while (<CPGOUT>) {
      chomp;

      print "$_\n";
      if (my ($id, $hid, $start, $end, $status) = /^(\S+)\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)/) {

	my $evalue;
	my $score = 0;
	my $hstart = 0;
	my $hend = 0;

	if ($status eq 'T') {
	  $evalue = '8e-5';
	} else {
# what evalue should we give to the hits with status = '?'
	  $evalue = '0.01';
	}

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
