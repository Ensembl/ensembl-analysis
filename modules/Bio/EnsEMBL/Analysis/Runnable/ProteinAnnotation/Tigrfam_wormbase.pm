package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Tigrfam_wormbase;

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
	      $self->database .' '.
	      $self->queryfile .' '.
              '> ' . $self->resultsfile;
    print STDERR "$cmd\n";   

    $self->throw ("Error running Tigrfam_wormbase ".$self->program." on ".$self->queryfile) 
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
	print STDERR "Tigrfam didn't find any hits\n";
	return; 
      } else {
	open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
      }
    }



#     my $id;
#     while (<CPGOUT>) {
#       chomp;

#       print "$_\n";

#       last if /^Alignments of top-scoring domains/;
#       next if (/^Model/ || /^\-/ || /^$/);
#       if (/^Query sequence:\s+(\S+)/) {
# 	$id = $1;
#       }

#       if (my ($hid, $start, $end, $hstart, $hend, $score, $evalue) = /^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {


# 	print "matched\n";
# 	$evalue = sprintf ("%.3e", $evalue);
# 	my $percentIdentity = 0;
      
# 	my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);
# 	push @fps, $fp;
#       }
#     }
#     close (CPGOUT); 
#     $self->output(\@fps);


  my $id;
  my $hid;
  while (my $line = <CPGOUT>) {
    chomp;
    if ($line =~ /Query:\s+(\S+)/) {$id = $1}
    if ($line =~ />>\s+(\S+)/) {$hid = $1} 
#   1 ?   10.5   0.0   2.1e-05       2.1      46      81 ..      23      58 ..      20      88 .. 0.90
    if (my ($evalue, $hstart, $hend, $start, $end) = ($line =~ /^\s+\d+\s+\?\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\.\.\s+(\d+)\s+(\d+)/)) {
      $evalue *= 1;           # force evalue to be a float not a string
      my $score = 0;
      my $percentIdentity = 0;
      my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);

      push @fps,$fp;      
    }
	

  }
  close (CPGOUT); 
  $self->output(\@fps);
}

1;
