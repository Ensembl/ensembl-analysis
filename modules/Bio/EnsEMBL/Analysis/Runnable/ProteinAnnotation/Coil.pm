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

 Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Coil

=head1 SYNOPSIS

 my $seqstream = Bio::SeqIO->new ( -file => $queryfile,
                                   -fmt => 'Fasta',
                                 );
 $seq = $seqstream->next_seq;

 my $ncoils = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Coil->new(
    -analysis => $ana,
    -query    => $query);
 $ncoils->run;
 my @list = @{$ncoils->run}

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Coil;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);

  if (not $self->database) {
    throw("You must supply the Coils runnable with a database directory");
  } elsif (not -d $self->database) {
    throw($self->database . " is not a *valid* database directory");
  } else {
    $ENV{'COILSDIR'} = $self->database;
  }

  return $self;
}


sub run_analysis {
  my ($self) = @_;

  my $command = $self->program." -f < ".$self->queryfile." > ".$self->resultsfile;
  
  # run program
  throw ("Error running ".$self->program." on ".$self->queryfile) 
      unless (system($command) == 0);
}


sub parse_results {
  my ($self) = @_;

  my ($fh, @pfs);

  my $resfile = $self->resultsfile;
    
  if (-e $resfile) {
    # it's a filename
    if (-z $resfile) {  
      return;
    }       
    else {
      open($fh, "<$resfile") or throw ("Error opening $resfile");
    }
  }
  else {
    # it'a a filehandle
    $fh = $resfile;
  }
  
  # parse
  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -fh  => $fh);
  while(my $seq = $seqio->next_seq) {
    my $pep = $seq->seq;

    while ($pep =~ m/(x+)/g) {
      my $end = pos($pep);
      my $start = $end - length($1) + 1;

      my $fp = $self->create_protein_feature($start, $end, 0, $seq->id, 0, 0,
                                             'ncoils', $self->analysis,
                                             0, 0);
      push @pfs, $fp;
    }
  }    

  $self->output(\@pfs);
}


1;
