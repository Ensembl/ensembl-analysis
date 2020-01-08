# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositePattern_wormbase;

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
	        $self->queryfile.' '.
		'> ' . $self->resultsfile;
    print STDERR "$cmd\n";   
    $self->throw ("Error running PrositePattern_wormbase ".$self->program." on ".$self->queryfile) 
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
        #sp|Q9BQN2|MT1P3_HUMAN   ps_scan|v1.57   PS00203 13      31      .       .       .       Name "METALLOTHIONEIN_VRT" ; LevelTag "(0)" ; Sequence "CtCAgsCkCkeCkCtsCkK" ; SequenceDescription "sp|Q9BQN2|MT1P3_HUMAN Putative metallothionein MT1P3 OS=Homo sapiens GN=MT1P3 PE=5 SV=1" ; KnownFalsePos 0
        #sp|A0RXF8|HEM1_CENSY    ps_scan|v1.57   PS00436 4       15      .       .       .       Name "PEROXIDASE_2" ; LevelTag "(-1)" ; Sequence "GLinARVtFHNS" ; SequenceDescription "sp|A0RXF8|HEM1_CENSY Glutamyl-tRNA reductase OS=Cenarchaeum symbiosumGN=hemA PE=3 SV=1" ; KnownFalsePos 13
        #sp|A0RXF8|HEM1_CENSY    ps_scan|v1.57   PS00747 97      120     .       .       .       Name "GLUTR" ; LevelTag "(0)"; Sequence "HLlrLTSGLDSmVVGEeQILGQVK" ; SequenceDescription "sp|A0RXF8|HEM1_CENSY Glutamyl-tRNA reductase OS=Cenarchaeum symbiosum GN=hemA PE=3 SV=1" ; KnownFalsePos 0
        ###  etc.


    my $id;
    while (<CPGOUT>) {
      chomp;

      print "$_\n";
      my $score = 0;
      my $hstart = 1;
      my $hend = -1;

      
      my ($id, $hid, $start, $end ) = /\s*(\S+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+\S+/;

#Calculate the length of the match using the values given for the sequence.
      $hend = $end - $start + 1;

      print "matched\n";
      my $evalue = 0.01;
      my $percentIdentity = 0;
      
      my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);
      push @fps, $fp;
    }
    close (CPGOUT); 
    $self->output(\@fps);
}

1;



