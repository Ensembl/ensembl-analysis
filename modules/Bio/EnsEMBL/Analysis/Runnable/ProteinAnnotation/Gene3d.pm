# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Genes3d

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Genes3d->new 
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => 'Pfam');

  $hmm->workdir ($workdir);
  $hmm->run;
  my @results = $hmm->output;

=head1 DESCRIPTION

  Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs hmmpfam.
  The resulting output file is parsed to produce a set of Bio::EnsEMBL::FeaturePairs.

=head1 CONTACT

   Marc Sohrmann: http://www.ensembl.org/Help/Contact

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Gene3d;

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

    # some of these options require HMMER 2.2g (August 2001)
    
    my @dbfiles = split(/;/,$self->analysis->db_file);

    print STDERR "FILENAME: ".$self->queryfile."\n";
 
# /software/worm/iprscan/bin/binaries/hmmer3/hmmscan --domE 1 -Z 100000 /data/blastdb/Worms/interpro_scan/iprscan/data/gene3d.lib inputfile > resultsfile

    my $cmd = $self->program .' '. 
	$self->analysis->parameters .' '.
	  $self->database .' '.
	    $self->queryfile.' > '.
	      $self->resultsfile;
    print STDERR "$cmd\n";   
    $self->throw ("Error running ".$self->program." on ".$self->queryfile." against ".$dbfiles[0]) 
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
      print STDERR "Genes3D didn't find any hits\n";
      return; }       
    else {
      open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
    }
  }

   
  my %seen_before;




#Domain annotation for each model (and alignments):
#>> 1ukxA00
#   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
# ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#   1 ?   10.5   0.0   2.1e-05       2.1      46      81 ..      23      58 ..      20      88 .. 0.90
#   2 ?   -2.7   0.0      0.24   2.4e+04      67      83 ..     136     152 ..     108     168 .. 0.77
#
#  Alignments for each domain:
#  == domain 1    score: 10.5 bits;  conditional E-value: 2.1e-05
#                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX RF
#      1ukxA00 46 sfkikikpesdeeeeqkvsltLqvelpetYPdeaPe 81
#                  fk+++++  de++ + + + L+++l ++ Pde+P
#  temp_gene_1 23 AFKFHLSTLVDEDDPTPFDFALSFQLKPELPDELPA 58
#                 5999999999************************96 PP
#
#  == domain 2    score: -2.7 bits;  conditional E-value: 0.24
#                  XXXXXXXXXXXXXXXXX RF
#      1ukxA00  67 LqvelpetYPdeaPeie 83
#                  + +  +  YPd+   ++
#  temp_gene_1 136 MLIYKTRYYPDTKCLVK 152
#                  55556677888887765 PP






#First parse what comes from the ls mode matches. Every match in that case is taken
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
