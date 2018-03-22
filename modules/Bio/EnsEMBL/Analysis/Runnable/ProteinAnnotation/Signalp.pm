# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

 Bio::EnsEMBL::Analysis::Runnable::Protein::Signalp

=head1 SYNOPSIS

 my $signalp = Bio::EnsEMBL::Analysis::Runnable::Protein::Signalp->new(
    -query => $seq,
    -analysis => $ana,
 );
 $signalp->run;
 my @results = @{$signalp->output};

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Signalp;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);


# We now leave the length of sequence to analyse at the default 70
## we over-ride write_seq_file because Signalp works best
## if only the first 50 a.a. of the protein is supplied 
#
#sub write_seq_file{
#  my ($self, $seq, $filename) = @_;
#  if(!$seq){
#    $seq = $self->query;
#  }
#  if(!$filename){
#    $filename = $self->queryfile;
#  }
#  my $seqout = Bio::SeqIO->new(
#                               -file => ">".$filename,
#                               -format => 'Fasta',
#                              );
#  eval{
#    my $init_seq = substr($seq->seq, 0, 50);
#    my $newseqobj = Bio::Seq->new(-id => $seq->display_id,
#                                  -seq => $init_seq);
#    $seqout->write_seq($newseqobj);
#  };
#  if($@){
#    throw("FAILED to write $seq to $filename Runnable:write_seq_file");
#  }
#  return $filename;
#}


sub run_analysis {
  my ($self) = @_;
 
  throw("Error running ".$self->program." on ".$self->queryfile) 
      unless ((system($self->program . " -t euk -f short "
                      .$self->queryfile 
                      . " > ".$self->resultsfile)) == 0); 
}


sub parse_results {
  my ($self) = @_;

  my ($fh);

  my $resfile = $self->resultsfile;
  
  if (-e $resfile) {
    # it's a filename
    if (-z $resfile) {  
      # no hits
     print STDERR "BLERG\n"; 

      return;
    }else {
      open ($fh, "<$resfile") or throw ("Error opening $resfile");
    }
  }else {
    # it'a a filehandle
    $fh = $resfile;
  }
 
  # parse
  my (@fps, $id, $end, $yesno);
  while (<$fh>) {
    print STDERR "$_" if $ENV{DEBUG};
    chomp;
    next if (/^#/);
    ($id, $end, $yesno) = (/^(\S+)\s+[\d\.]+\s+\d+\s+[\d\.]+\s+(\d+)\s+[\d\.]+\s+\d+\s+[\d\.]+\s+[\d+\.]+\s+(\S)/);
    if (!defined $id) {
      # throw an exception
      $self->throw ("parsing problem in ".$self->program);
    }
    if ($yesno eq 'Y') {
      my $fp = $self->create_protein_feature(1, $end-1, 0, $id, 0, 0,
					     'Sigp', $self->analysis,
					     0, 0);
      push @fps, $fp;
    }
    
  }
  close($fh);

  $self->output(\@fps);
}



1;
