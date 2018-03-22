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

  Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Seg

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Seg->new ( -db      => $db,
	    	                                                    -input_id   => $input_id,
                                                                    -analysis   => $analysis,
                                                                  );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->output;
  $seg->write_output; # writes features to to DB

 NB: The input_id can either be a peptide id or the location for a protein file. 

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg
  to add functionality to read and write to databases.
  The query sequence is provided through the input_id.
  The appropriate Bio::EnsEMBL::Analysis object
  must be passed for extraction of parameters.

=head1 CONTACT

  Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

  The rest of the documentation details each of the object methods. 
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Seg;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);

sub fetch_input {
  my ($self, @args) = @_;

  $self->SUPER::fetch_input(@args);

  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg->new(-query     => $self->query,
                                                                          -analysis  => $self->analysis);
  $self->runnable($run);
}


1;
