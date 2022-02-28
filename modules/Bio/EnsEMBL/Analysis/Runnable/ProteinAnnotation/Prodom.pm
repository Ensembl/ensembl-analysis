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

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prodom

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prodom->new 
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => 'sf_hmm');

  $hmm->workdir ($workdir);
  $hmm->run;
  my @results = $hmm->output;

=head1 DESCRIPTION

  Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs hmmpfam.
  The resulting output file is parsed to produce a set of Bio::EnsEMBL::FeaturePairs.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prodom;

use warnings ;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);

sub run_analysis {
  my ($self) = @_;

  # run program
  print STDERR "running ".$self->program." against ".$self->database."\n";
   
  my @dbfiles = split(/;/,$self->analysis->db_file); # hmmm ... only the first bit is used

  print STDERR "FILENAME: ".$self->queryfile."\n";
  
                                                                                                                   
  my $cmd = $self->program .' '.
	$self->options .' '.
	'-P ' . '/software/worm/iprscan/bin/binaries/blast/ ' .
	'-p ' . 'blastp' . ' '.
	'-d ' . $dbfiles[0] .' '.
	'-s ' . $self->queryfile . ' '.
	'-t ' . '/tmp' . ' '.
	'-h 0 -f '.
	' >' . $self->resultsfile;

  print STDERR "$cmd\n";   
    $self->throw ("Error running ".$self->program." on ".$self->filename." against ".$dbfiles[0]) 
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
	if (-z $self->resultsfile) {  
	    print STDERR "Prodom didn't find any hits\n";
	    return; }       
	else {
	    open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
	    }
    }

#5_Translation_id_5_gene_5      1     40 //  pd_PD000035;sp_TLL_DROVI_O16845;      33    102 // S=135    E=2e-08  //  (1021)  RECEPTOR NUCLEAR TRANSCRIPTION DNA-BINDING REGULATION ZINC-FINGER METAL-BINDING ZINC HORMONE FAMILY    Length = 70 
# Example output:
#CE27771    203    318 //  pd_PD002460;sp_Q19242_CAEEL_Q19242;     199    344 // S=122    E=6e-06  //  (189)  BIOSYNTHESIS MOLYBDOPTERIN MOEA COFACTOR MOLYBDENUM ENZYME PROBABLE LONG GEPHYRIN DOMAIN        Length = 146
#//

  while (<CPGOUT>) {
    chomp;
	
    if (my ($id, $start, $end, $hid, $hstart, $hend, $score, $evalue) = /^(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+pd_(\S+);\S+\s+(\d+)\s+(\d+)\s+\S+\s+S=(\S+)\s+E=(\S+)/) {
      my $fp= $self->create_protein_feature($start,$end,$score,$id,$hstart,$hend,$hid,$self->analysis, sprintf("%.3e", $evalue),0);
	    
      push @fps,$fp;
    }
  }
  close CPGOUT;
  $self->output(\@fps);
}

1;
