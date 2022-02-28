# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Superfamily_wormbase;

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
    print STDERR "running ".$self->program." against ".$self->queryfile."\n";

    print STDERR "FILENAME: ".$self->queryfile."\n";
 
    my $cmd = 	$self->program .' '.
	        $self->analysis->parameters .' '.
	        '-i ' . $self->queryfile.' '.
		'-o ' . $self->resultsfile;
    print STDERR "$cmd\n";   
    $self->throw ("Error running Superfamily_wormbase ".$self->program." on ".$self->queryfile) 
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
	print STDERR "Superfamily didn't find any hits\n";
	return; 
      } else {
	open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
      }
    }


# example output lines:
#R186.4  63882   0019507 9       181     2.6e-37 msssgllkPATLDDVFQKLEDLCKLFpPQEKTVNVTSLkTGRILAEDIITEYDIPAQRTSIVDGFAIIVNQLGTKREIVGLSTAVTPYNAELISNECVRITIGGVVPDGADTVVPIENVALlkEEKCIEVLRKPKEGDNIREVGSEAKTGEILLKDGHHLDTMSITLLHALGISQVEIYKKprvcvlsigsdlnsnkmygsfnrsqllelfqsqgftaidagsstehiteveekirtaasfacvlvtvggaqvirevaktlkfkfeiqdvdstpgnftvstgkidetpvvlsifpeyhvsswiganlfvspilramegqnsetshrfkaeltqpisktsetrflrarsevskgnlistplgcedifgansilevksntcfsagdvvdlrfa      Molybdenum cofactor biosynthesis protein MoeA, N-terminal and linker domains
#//
#R186.4  53218   0017938 180     353     3.9e-14 msssgllkpatlddvfqkledlcklfppqektvnvtslktgrilaediiteydipaqrtsivdgfaiivnqlgtkreivglstavtpynaelisnecvritiggvvpdgadtvvpienvallkeekcievlrkpkegdnirevgseaktgeillkdghhldtmsitllhalgisqveiyKKPRVCVLSIGSDLNSNKMYgSFNRSQLLELFQSQGFTAIDAGSSTEHITEVEEKIRTAASFACVLVTVGGAQVIREVAKTLKFKFEIQDVDSTPGNFTVSTGKIDetPVVLSIFPEYHVSSWIGANLFVSPILRAMEGQNSETSHRFKAELTQPISKTSETRFLRARSEVSKGnlistplgcedifgansilevksntcfsagdvvdlrfa      Molybdenum cofactor biosynthesis proteins
#//

    my $id;
    while (<CPGOUT>) {
      chomp;
      next if (/^\/\//);

      print "$_\n";
      if (my ($id, $hid, $start, $end, $evalue) = /^(\S+)\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)/) {
	print "matched\n";
	my $score = 0;
	my $hstart = 0;
	my $hend = 0;
	my $hid = "SSF$hid";
	my $percentIdentity = 0;
	$evalue *= 1;		# force evalue to be a float not a string

	print "output= $start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity\n";
	my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);
	
	push @fps, $fp;
      }
    }
    close (CPGOUT); 
    $self->output(\@fps);
    for ( @fps ) {
       print "value : " . $_->p_value() . "\n" ;  
    } 
}

1;
