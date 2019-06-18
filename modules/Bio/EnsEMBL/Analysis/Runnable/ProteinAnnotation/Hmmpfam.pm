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
=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Protein::Hmmpfam

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Analysis::Runnable::Protein::Hmmpfam->new
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => '/data/Pfam_ls');

  $hmm->run;
  my @results = @{$hmm->output};

=head1 DESCRIPTION

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam;

use warnings ;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;


@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);



sub run_analysis {
  my ($self) = @_;
  
  my $cmd = $self->program 
      . ' ' . $self->analysis->parameters
      . ' ' . $self->database 
      . ' ' . $self->queryfile 
      .' > '. $self->resultsfile;

  print STDERR "Running:$cmd\n";
  
  throw ("Error running ".$self->program." on ".
         $self->queryfile." against ".$self->database)unless 
         ((system ($cmd)) == 0);
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

  my $resfile = $self->resultsfile();
  my @hits;

  if (-e $resfile) {
    if (-z $resfile) {
      print STDERR "Tigrfam didn't find any hits\n";
      return;
    } else {
      open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
    }
  }

  my $id;
  my $hid;

  while (<CPGOUT>) {
      chomp;

## Query line identifier with Hmmer3
      if (/^Query:\s+(\S+)/) {
        $id = $1;
        next;
      }

## Query line identifier with Hmmer2
      if (/^Query sequence:\s+(\S+)/) {
        $id = $1;
      }


# With Hmmer3, hit name is on a separate line
      if (/^>> (\w+)/) {
        $hid = $1 ;
      }

## Result line with Hmmer2
      if (my ($hid,
              $start,
              $end,
              $hstart,
              $hend,
              $score,
              $evalue) = /^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {

        $evalue = sprintf ("%.3e", $evalue);
        my $percentIdentity = 0;

        # Remove the version at the end of the Pfam accession number if any
        # e.g. 'PF00001.13' => 'PF00001'
        $hid =~ s/\.\d+$//;

        my $fp = $self->create_protein_feature($start,
                                               $end,
                                               $score,
                                               $id,
                                               $hstart,
                                               $hend,
                                               $hid,
                                               $self->analysis,
                                               $evalue,
                                               $percentIdentity);
        push @hits, $fp;
      }


## Result line with Hmmer3
      if (my ($score,
              $evalue,
              $hstart,
              $hend,
              $start,
              $end) = /^\s+\d+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)/) {

        $evalue = sprintf ("%.3e", $evalue);
        my $percentIdentity = 0;

        my $fp = $self->create_protein_feature($start, 
                                               $end, 
                                               $score, 
                                               $id, 
                                               $hstart, 
                                               $hend, 
                                               $hid,
                                               $self->analysis, 
                                               $evalue,
                                               $percentIdentity);
        push @hits, $fp;
      }
  }
  close (CPGOUT);
  $self->output(\@hits);
}

1;
