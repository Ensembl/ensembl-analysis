# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Smart_wormbase;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);

###################
# analysis methods
###################

=head2 run_analysis

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

    $self->throw ("Error running Smart_wormbase ".$self->program." on ".$self->queryfile) 
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
      print STDERR "Smart didn't find any hits\n";
      return; 
    } else {
      open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
    }
  }



  my $id;
  while (<CPGOUT>) {
    chomp;

    print "$_\n";

    last if /^Alignments of top-scoring domains/;
    next if (/^Model/ || /^\-/ || /^$/);
    if (/^Query sequence:\s+(\S+)/) {
      $id = $1;
    }

    if (my ($hid, $start, $end, $hstart, $hend, $score, $evalue) = /^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {


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
