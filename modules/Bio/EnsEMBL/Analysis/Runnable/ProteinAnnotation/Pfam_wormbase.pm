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

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Pfam_wormbase;

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

    $self->throw ("Error running Pfam_wormbase ".$self->program." on ".$self->queryfile) 
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
	print STDERR "Pfam didn't find any hits\n";
	return; 
      } else {
	open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
      }
    }



    my %seen_before;

#    my $id;
#    while (<CPGOUT>) {
#      chomp;

#      print "$_\n";

#      last if /^Alignments of top-scoring domains/;
#      next if (/^Model/ || /^\-/ || /^$/);
#      if (/^Query sequence:\s+(\S+)/) {
#	$id = $1;
#      }

#      if (my ($hid, $start, $end, $hstart, $hend, $score, $evalue) = /^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {

#	# remove .fs or .ls at end of Pfam ID
#	$hid =~ s/\.fs$//;
#	$hid =~ s/\.ls$//;

#	print "matched\n";
#	$evalue = sprintf ("%.3e", $evalue);
#	my $percentIdentity = 0;
      
#	# only output the unique hits
#	if (!exists $seen_before{"${hid}_${start}_${end}"}) {
#	  my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);
#	  push @fps, $fp;
#	  $seen_before{"${hid}_${start}_${end}"} = 1;
#	  print "Outputting: $hid $start $end $evalue\n";
#	}
#      }
#    }
#    close (CPGOUT); 
#    $self->output(\@fps);


    my $mstatdir = $_[0];
    my $cutoff   = $_[1];
    my $pre      = $_[2];
    my $model;
    my $seqid;
    my $start;
    my $end;
    my $score;
    my $E;
    my $evalue;
    my $hid;
    my $hstart;
    my $hend;
      


    my $id;
    while (<CPGOUT>) {
      last if /^Internal pipeline statistics summary/;

      if(/^Query\:\s+(\S+)/) {
	#Query:       C13D9.8
	$id = $1;
      }

      if (/^>>\s+(\S+)\./) {
	$hid = $1;
      }

      #Domain annotation for each model:
      #>> PF01699.17  Sodium/calcium exchanger protein
      #   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
      # ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
      #   1 !   82.8   6.0   1.1e-27   1.3e-23       1     138 [.      67     200 ..      67     201 .. 0.96
      #   2 !   75.8  10.1   1.6e-25   1.9e-21       1     138 [.     502     639 ..     502     640 .. 0.88

      
      if (/^\s*\d+\s+\!\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S\S\s+\d+\s+\d+\s+\S\S\s+\S+/) {
	$score    = $1;
	$evalue   = $2;
	$hstart   = $3;
	$hend     = $4;
	$start    = $5;
	$end      = $6;

	$evalue = sprintf ("%.3e", $evalue);
	my $percentIdentity = 0;
	  
	  # only output the unique hits
	  if (!exists $seen_before{"${hid}_${start}_${end}"}) {
	    my $fp= $self->create_protein_feature($start, $end, $score, $id, $hstart, $hend, $hid, $self->analysis, $evalue, $percentIdentity);
	    push @fps, $fp;
	    $seen_before{"${hid}_${start}_${end}"} = 1;
	  }
	}
    }
    close (CPGOUT); 
    $self->output(\@fps);
}

1;
