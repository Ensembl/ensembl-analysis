=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Fgenesh - 

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::Fgenes->new
  (
   -query => $slice,
   -program => 'genscan',
   -matrix => 'HumanIso.smat',
   -analysis => $analysis,
  );
$runnable->run;
my @predictions = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run the fgenesh gene predictor and then parse the results
into prediction transcripts

this is a bare bones module which inherits most of its functionality from
Bio::EnsEMBL::Analysis::Runnable::Genscan

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::Fgenesh;



use vars qw(@ISA);
use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio);

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Fgenesh
  Function  : create a new Bio::EnsEMBL::Analysis::Runnable::Fgenesh object
  and setting up default values
  Returntype: Bio::EnsEMBL::Analysis::Runnable::Fgenesh
  Exceptions: 
  Example   : 

=cut


sub new {
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);

  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->program('fgenesh') if(!$self->program);
  $self->matrix('hum.dat') if(!$self->matrix);
  ######################

  return $self;
}




#FGENESH-1.0 Prediction of potential genes in genomic DNA
#Time:	Thu Sep	 2 10:57:44 2004.
#Seq name: HS120G22 AL031847.17 Human DNA sequence from clone RP1-120G22 on chromosome 1p36.21-36.33
#length of sequence  166518bp  G+C content: 56 Isochore: 3
#number of predicted genes 18 in +chain 13 in -chain 5
#number of predicted exons 97 in +chain 72 in -chain 25

#  Gn S	 Type	Start	      End	  Score		  ORF		    Len
#  -- -	 ----	-----	      ---	  -----		  ---		    ---
#   1 +	 CDSf	  827 -	    849	   6.70	    827 -	847	  21
#   1 +	 CDSi	 4024 -	   4166	  -1.97	   4025 - 4165  141
#   1 +	 CDSi	 5726 -	   5809	   0.50	   5728 - 5808	81
#   1 +	 CDSi	 7096 -	   7154	   7.92	   7098 - 7154	57
#   1 +	 CDSl	 7429 -	   7740	   7.18	   7429 - 7740  312
#   1 +	 PolA	 7762		   2.02

#   2 +	 TSS	 9875		  -2.59
#   2 +	 CDSf	10848 -	  10925	  11.04	  10848 - 10925	78
#   2 +	 CDSi	14545 -	  14719	   5.83	  14545 - 14718 174
#   2 +	 CDSi	14883 -	  14999	  20.17	  14885 - 14998 114
#   2 +	 CDSl	15087 -	  15169	  -2.37	  15089 - 15169	81
#   2 +	 PolA	15355		  -2.78


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Fgenesh
  Arg [2]   : string, filename
  Function  : parse the output from Fgenesh into prediction transcripts
  Returntype: none
  Exceptions: throws if cannot open or close results file
  Example   : 

=cut

sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  open(OUT, "<".$results) or throw("FAILED to open ".$results.
                                   "Genscan:parse_results");

  if (<OUT> =~ m|no reliable predictions|i ){
    print STDERR "No genes predicted\n";
    return;
  }
  while (<OUT>) {
    chomp;
    my $flag = 0;
    if (/^\s*-[-\s]+$/) {
      GENES : while (<OUT>) {
        my @lines;
        until (/^$/) {
          if (/CDSl|CDSi|CDSf|CDSo/i) {
            my @element = split;
            push @lines, \@element; 
            throw("Unable to parse fgenesh ouput (".@element.
                  ") Line: $_\n") unless (scalar(@element) == 11);
          }elsif (/Predicted protein/i) {
            $flag = 1;
            last GENES ;
          }
          #print "Getting next line\n";
          $_ = <OUT>; 
        }
        if ($lines[0]->[1] eq '+') {
          @lines = sort {$a->[3] <=> $b->[3]} @lines;
        } elsif ($lines[0]->[1] eq '-') {
          @lines = reverse sort {$a->[3] <=> $b->[3]} @lines;
        }
        my $exon_num=1;
         foreach my $line (@lines) {
           my ($start, $end, $strand, $score, $pvalue, $phase, $seqname);
           $seqname = $line->[0]+($exon_num/100);
           if ($line->[1] eq '+') {
             $start = $line->[3];
             $end = $line->[5];
             $strand = 1;
             $phase = (3-($line->[7]-$line->[3]))% 3;
           } elsif ($line->[1] eq '-') {
             $start = $line->[3];
             $end = $line->[5];
             $strand = -1;
             $phase = (3-($line->[5]-$line->[9]))% 3;
           }
           $score = $line->[6];
           my $exon = $self->feature_factory->create_prediction_exon
             ($start, $end, $strand, $score, '0', $phase, $seqname, 
              $self->query, $self->analysis);
           $self->exon_groups($line->[0], $exon);
           $exon_num++;
         }
        if ( $flag==1 ) {
          last;
        }
      }
    }
  }
  $self->create_transcripts();
  close(OUT) or throw("FAILED to close ".$results.
                      "Fgenesh:parse_results");
}

1;
