# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prints;

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

  print STDERR "running ".$self->program." against ".$self->database."\n";
  print STDERR "FILENAME: ".$self->queryfile."\n";

  # because Prints reports matches off the end of the sequence, 
  # we need the sequence lengths to be able to trim them back
      
  my $command =  $self->program ." " . 
      $self->database . " " . 
      $self->queryfile . " " .
      $self->analysis->parameters . " " . 
      "> " . $self->resultsfile;

  print STDERR "$command\n";
  throw("Failed during prints run $!\n") unless 
      system($command) == 0 ;
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

  my %seqlengths;
  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -file   => $self->queryfile);
  while (my $seq = $seqio->next_seq) {
    $seqlengths{$seq->id} = $seq->length;
  }
  $seqio->close;
    
  my ($fh);

  my $resfile = $self->resultsfile;
  
  if (-e $resfile) {    
    if (-z $resfile) {  
      # No hits found
      print STDERR "Prints didn't find any hits\n";
      return; 
    }       
    else {
      open ($fh, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");
    }
  }

  my @fps ;
  my %printsac ;
  my $seq_id ;
  my %evalue;

  
  while (<$fh>) {
    my $line = $_;
    chomp $line;
    # Pattern match the Sn; field which should contain the SequenceId and Accession
    
    if ($line =~ s/^Sn;//) { # We have identified a Sn; line so there should be the following:	    
      #ENSP00000003603 Gene:ENSG00000000003 Query:AL035608 Contig:AL035608.00001 Chr:chrX basepair:97227305
      ($seq_id) = $line =~ /^\s*(\w+)/;
    }
        
    if ($line =~ s/^1TBH//) {
      my ($id) = $line =~ /^\s*(\w+)/;
      my ($ac) = $line =~ /(PR\w+)[;\.]*\s*$/;
#      print STDERR "In 1TBH data, printsac{$id} =  $ac\n";
      $printsac{$id} = $ac;
    }

    if ($line =~ /^2TBH/) {
      my @line = split /\s+/, $line;
      $evalue{$line[1]} = $line[9];
#      print STDERR "hash evalue of $line[1] = $line[9]\n";
    }
    
    if ($line =~ s/^3TB//) {
      if ($line =~ s/^[HN]//) {
        my ($num,$temp1,$tot1);
        # Grab these lines
        #  1433ZETA        1  of  6  88.19   1328    1.00e-16  ELTVEERNLLSVAYKNVIGARRASWRIITS   30   35   36   48

        $line =~ s/^\s+//;        
        my @elements = split /\s+/, $line; 
        
        # Name each of the elements in the array
        my ($fingerprintName,
            $motifNumber,
            $temp,
            $tot,
            $percentageIdentity,
            $profileScore,
            $pvalue,
            $subsequence,
            $motifLength,
            $lowestMotifPosition,
            $matchPosition,
            $highestMotifPosition) = @elements;
        
        # If protein is 10,000+ residues (i.e. titin), then last two elements are merged, e.g.:
        # VEGFRECEPTOR    5  of  6  39.20   406     1.05e-05  LIVRNARKENAGKYTLVL                                      18   374  13564392
        if (!defined $highestMotifPosition) {
          # First five characters of $matchPosition is actual $matchPosition
          $highestMotifPosition = $matchPosition;
          $matchPosition = substr($highestMotifPosition, 0, 5, '');
        }
        
        my $start = $matchPosition;
        my $end = $matchPosition + $motifLength - 1;
        
        # Prints sometimes reports -1 for the motif position start;
        $lowestMotifPosition = 1 if $lowestMotifPosition == -1;
        $highestMotifPosition = $lowestMotifPosition + $motifLength -1 
            if $highestMotifPosition == -1;

        # if we don't have a valid hit for this PRINTS ID, then ignore teh rest
        if (! exists $printsac{$fingerprintName}) {next;}
        my $feat = "$printsac{$fingerprintName},$start,$end,$percentageIdentity,$profileScore,$evalue{$fingerprintName}";
#        print STDERR "features= $feat\n" ;

        # It's possible that near the ends of the sequence,
        # fragment matches are found. These can result in a 
        # start < 1 or end > seq_length. Make a policy decision
        # prune these back in the sequence, but not the hit
        $start = 1 if $start < 1;
        $end = $seqlengths{$seq_id} if $end > $seqlengths{$seq_id};
        my $evalue = $evalue{$fingerprintName};
    
        my $fp = $self->create_protein_feature($start, $end, 
                                               $profileScore,
                                               $seq_id, 
                                               $lowestMotifPosition, 
                                               $highestMotifPosition, 
                                               $printsac{$fingerprintName},
                                               $self->analysis, 
                                               $evalue, 
                                               $percentageIdentity);
        push @fps, $fp;
      }
    }
  }
  close($fh);

  $self->output(\@fps);
}


1;

