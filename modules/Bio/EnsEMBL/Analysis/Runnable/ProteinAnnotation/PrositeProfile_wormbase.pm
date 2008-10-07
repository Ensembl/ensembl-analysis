package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositeProfile_wormbase;

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
    print STDERR "running ".$self->program." ".$self->analysis->parameters ." against ".$self->database."\n";

    print STDERR "FILENAME: ".$self->queryfile."\n";
 
    my $cmd = $self->program .' '.
	        $self->analysis->parameters .' '.
	        $self->queryfile.' '.
		'> ' . $self->resultsfile;
    print STDERR "$cmd\n";   
    $self->throw ("Error running PrositeProfile_wormbase ".$self->program." on ".$self->queryfile) 
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
	print STDERR "PrositeProfile didn't find any hits\n";
	return; 
      } else {
	open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
      }
    }


 # Example output for sequence ID 2345:
#2345    ps_scan|v1.4    PS50082 329     362     9.205   .       .       Name "WD_REPEATS_2" ; Level 0 ; RawScore 244 ; FeatureFrom 1 ; FeatureTo -9 ; Sequence "MKSYFGGLLCLAWSPDARYIVTGGEDDLITVYSV--------" ; KnownFalsePos 1
#2345    ps_scan|v1.4    PS50294 285     385     12.947  .       .       Name "WD_REPEATS_REGION" ; Level 0 ; RawScore 389 ; FeatureFrom 1 ; FeatureTo -28 ; Sequence "WAVGSGTLHEFAFSPSddTKLLATVSQDGFLRIFNYHTMELLAYMKSYFGGLLCLAWSPDARYIVTGGEDDLITVYSVVEKRVVCRGQGHRSWISKVAFDP---------------------------" ; KnownFalsePos 0



    my $id;
    while (<CPGOUT>) {
      chomp;

      print "$_\n";
      my $score = 0;
      my $hstart = 1;
      my $hend = -1;

      # is the first column the sequence ID or the file name? - check this
      my ($id, $hid, $start, $end, $evalue ) = /\s*(\S+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+([\d\.\-\e]+)\s+\S+\s+\S+/;

#Calculate the length of the match using the values given for the sequence.
      $hend = $end - $start + 1;

      print "matched\n";
#      my $evalue = 0.01;
      my $percentIdentity = 0;
      
      my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);
      push @fps, $fp;
    }
    close (CPGOUT); 
    $self->output(\@fps);
}

1;
