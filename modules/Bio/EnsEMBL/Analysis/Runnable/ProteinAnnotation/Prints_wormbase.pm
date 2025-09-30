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
package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prints_wormbase;

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

###########################
### The E-value calculation
###########################

# The reported P-value of any fingerprint result, is the product of the
# p-values for each motif. The motif p-values represent the probabilty
# that a comparison between the motif and a random sequence would
# achieve a score greater than or equal to the score attributed to the
# match between your query sequence and the motif.
#
# However, when a database search is performed this value must be
# adjusted for the multiple comparisons made therein. An E-value can be
# computed which represents the Expected numberof occurences of
# sequences scoring greater than or equal to the query's score.
#
# There are a number of estimates of this parameter;
#
# E = pD, which states that the probability multiplied by the sample
# database size (D equals the number of sequences), provides the number
# of occurrences, the assumption made here is that all sequences in a
# database are equally likely to be related to a particular motif.
#
# E = pN/n, is the same calculation, but the size of the database is
# calculated by taking the total length of the database in residues and
# dividing it by the width of the scoring entity (in this case the motif
# or fingerprint). This calculation assumes that the chance of a
# sequence being related to a motif (or fingerprint) is proportional to
# length (rather than being equally likely). Incidentally this latter
# approach is the method of choice for calculating BLAST database
# E-values (NCBI).
#
# FingerPRINTScan has now been modified to conform to the E = pN/n
# method of database E-value calculation, however P-value calculations
# remain unchanged
#
# The database from which D and N are calculated is (for example):
#
# SWISS-PROT 37 (77,977 sequences and 28,777,979 residues) and
# TrEMBL 9 (179,066 sequences and 55,577,465 residues)
#
# Therefore in all calculations D = 257,043 and N = 84,355,444.
#
# The current default values used in the /usr/local/ensembl/bin/FingerPRINTScan
# script are:
#         -E #1 #2  E-value calculation parameters.
#              (where #1 is the number of sequences in the primary database (default 80000))
#              (where #2 is the number of resides   in the primary database (default 2.96103e+07))
#              ( the default values are based on SWISS-PROT 38)
#
# The database may get updated out of synch of the scripts, so you may
# need to set '-E #1 #2' explicitly in the parameters= line of the config file.



sub run_analysis {
    my ($self) = @_;

    # run program
    print STDERR "running ".$self->program." against ".$self->database."\n";

    print STDERR "FILENAME: ".$self->queryfile."\n";
 
    my $cmd = $self->program .' '.
		$self->database .' '.
                $self->queryfile. ' '.
	        $self->analysis->parameters .' '.
		'> ' . $self->resultsfile;
    print STDERR "$cmd\n";   
    $self->throw ("Error running Prints_wormbase ".$self->program." on ".$self->filename) 
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
      print STDERR "Prints didn't find any hits\n";
      return; 
    } else {
      open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n"); # 
    }
  }


  my $id;
  my $line;
  
  my $sequenceId;
  my @features;
  my %evalue;
  my %printsac;

  while (<CPGOUT>) {
    $line = $_;

    print STDERR "Next line: $line";
    chomp $line;
    # Pattern match the Sn; field which should contain the SequenceId and Accession
    
    if ($line =~ s/^Sn;//) { # We have identified a Sn; line so there should be the following:

      #ENSP00000003603 Gene:ENSG00000000003 Query:AL035608 Contig:AL035608.00001 Chr:chrX basepair:97227305
      ($sequenceId) = $line =~ /^\s*(\w+)/;
    }


    if ($line =~ s/^1TBH//) {
      my  ($id) = $line =~ /^\s*(\w+)/;
      my ($ac) = $line =~ /(PR\w+)[;\.]*\s*$/;
      $printsac{$id} = $ac;
      print STDERR "1TBH line = $line\n";
      print STDERR "In 1TBH data, printsac{$id} =  $ac\n";
    }

# get the evalues of each of the fingerprint names
    if ($line =~ /^2TBH/) {
      print STDERR "got a 2TBH line: $line\n";
      my @line = split /\s+/, $line;
      $evalue{$line[1]} = $line[9];
      print STDERR "hash evalue of $line[1] = $line[9]\n";
    }

    if ($line =~ s/^3TB//) {
      if ($line =~ s/^[HN]//) {
	my ($num,$temp1,$tot1) = "";
	# Grab these lines
	#       1433ZETA        1  of  6  88.19   1328    1.00e-16  ELTVEERNLLSVAYKNVIGARRASWRIITS  30   35   36   48
	# split line on space, hence strip off all leading spaces first.
	$line =~ s/^\s+//;
	print STDERR "line = $line\n";
	# Place all elements of list into an array
	my @elements = split /\s+/, $line;

	# Name each of the elements in the array
	my ($fingerprintName,$motifNumber,$temp,$tot,$percentageIdentity,$profileScore,$pvalue,$subsequence,$motifLength,$lowestMotifPosition,$matchPosition,$highestMotifPosition) = @elements;
	print STDERR "fingerprintName=$fingerprintName\n";
	my $start = $matchPosition;
	#
	# If the match to the pattern lies at the end of the protein we might get padding of the subsequence with #'s, and the
	# end position will be bigger than the actual end of the protein. So we'll strip the #'s off the end, adjust the
	# motif length accordingly, and only then derive the match end.
	my $hash_substring;
	my $end;
	
	if($subsequence =~ /(\#+)$/){
	  $hash_substring = $1;
	  $end = $matchPosition + $motifLength - 1 - length($hash_substring);
	} else {
	  $end = $matchPosition + $motifLength - 1;
	}
	
# if we don't have a valid hit for this PRINTS ID, then ignore teh rest
	if (! exists $printsac{$fingerprintName}) {next;} 

	my $print =  $printsac{$fingerprintName};
	my $feat = "$print,$start,$end,$percentageIdentity,$profileScore,$evalue{$fingerprintName}";
	print STDERR "features= $feat\n";
	
	print "matched\n";
	my $hstart = 0;
	my $hend = 0;
	my $evalue = $evalue{$fingerprintName};

	print STDERR "writing to database\n";
	my $fp= $self->create_protein_feature($start, $end, $profileScore, $sequenceId, $hstart, $hend, $print, $self->analysis, $evalue, $percentageIdentity);
	push @fps, $fp;
      }
    }
  }

  close (CPGOUT); 
  $self->output(\@fps);
}

1;
