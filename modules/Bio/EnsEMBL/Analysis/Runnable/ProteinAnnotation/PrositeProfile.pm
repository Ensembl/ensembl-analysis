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


package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositeProfile;
use warnings ;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object


use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);


sub multiprotein{
  my ($self) = @_;
  return 0;
}


sub run_analysis {
  my ($self) = @_;
  
  throw("Failed during Profile run $!\n") unless 
    (system ($self->program . ' -f ' . $self->queryfile. ' ' .
             $self->database . ' > ' .$self->resultsfile) == 0) ;
 
}


sub parse_results {
  my ($self,$seqid) = @_;
  
  my ($fh);
  my $resfile = $self->resultsfile;
  
  if (-e $resfile) {	
    if (-z $resfile) {  
      return; 
    } else {
      open ($fh, "<$resfile") or throw("Error opening ", $resfile,);
    }
  }
  
  my (@pfs);
  while (<$fh>) {
    if (/^\s*(\S+)\s+(\d+)\s*pos\.\s+(\d+)\s+\-\s+(\d+)\s+(\w+)\|/) {
      my ($sc, $rsc, $st, $en, $acc) = ($1, $2, $3, $4, $5);
      my $fp = $self->create_protein_feature($st,
                                             $en,
                                             $sc,
                                             $seqid,
                                             0, 0,
                                             $acc,
                                             $self->analysis,
                                             0, 0);
      push @pfs, $fp;
    }
  }
  
  $self->output(\@pfs);  
}


