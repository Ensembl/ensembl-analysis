package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Pfam;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;
use Env qw(PATH);


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
	      '-o ' . $self->resultsfile .' '.
              $self->database . ' ' .
	      $self->queryfile;
    print STDERR "$cmd\n";   

    $PATH = '/software/worm/bin/hmmer/:/software/worm/iprscan/bin/binaries/:'.$PATH; #for hmmpfam
    $PATH = '/software/worm/iprscan/bin/binaries/blast/:'.$PATH; #for blastall

    $self->throw ("Error running Pfam ".$self->program." on ".$self->queryfile) 
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
	print STDERR "Pfam didn't find any hits\n";
	return; 
      } else {
	open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
      }
    }


 # Example output for sequence ID 2345:
#2345   279  319 PF00400.23     1   38    13.7     0.008  WD40
#2345   323  361 PF00400.23     1   38    22.0    0.0023  WD40

# now
# 12311     49   288 PF00149.20      1   124 ls    90.9   4.4e-24       Metallophos         
# 12311    294   517 PF04152.6       1   240 ls   310.9   2.6e-90    Mre11_DNA_bind

# 11845      1   167 PF06653.3       1   160 ls   -12.3   0.00016           DUF1164 
 
    my $id ;
    my $hid ;
    while (<CPGOUT>) {
      chomp;

      if (/^Query:\s+(\S+)/) {
        $id = $1;
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
