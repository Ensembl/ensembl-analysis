# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Fgenes
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Fgenes

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::Fgenes->new(
      -query => $slice,
      -program => 'genscan',
      -matrix => 'HumanIso.smat',
     );
  $runnable->run;
  my @predictions = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run the fgenesh gene predictor and then parse the results
into prediction transcripts

this is a bare bones module which inherits most of its functionality from
Bio::EnsEMBL::Analysis::Runnable::Genscan

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Fgenesh;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::Genscan;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Genscan);


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
  $self->program('fgenesh') if($self->program eq 'genscan'
                              || !$self->program);
  $self->matrix('hum.dat') if($self->matrix eq 'HumanIso.smat'
                              || $self->matrix eq
                              '/usr/local/ensembl/data/HumanIso.smat'
                              || !$self->matrix);
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
  open(FH, $results) or throw("FAILED to open ".$results.
                              "Fgenesh:parse_results");
  my $ff = $self->feature_factory;
  my $exon_count = 1;
  my $gene_num = 1;
 LINE:while(<FH>){
    chomp;
    if(m|no reliable predictions|i){
      print "No Genes prediction\n";
      last LINE;
    }
    if (/CDSl|CDSi|CDSf|CDSo/i){
      my @line = split;
      if($line[0] != $gene_num){
        $gene_num = $line[0];
        $exon_count = 1;
      }
      my $exon_name = $gene_num.".".$exon_count;
      my ($start, $end, $strand, $phase);
      if($line[1] eq '+'){
        $start = $line[3];
        $end = $line[5];
        $strand = 1;
        $phase = (3-($line[7]-$line[3]))%3;
      }elsif($line[1] eq '-'){
        $start = $line[3];
        $end = $line[5];
        $strand = -1;
        $phase = (3-($line[5]-$line[9]))%3;
      }
      my $score = $line[6];
      my $exon = $ff->create_prediction_exon($start, $end, $strand, 
                                             $score, 0, $phase, 
                                             $exon_name, $self->query,
                                             $self->analysis);
      $self->exon_groups($gene_num, $exon);
    }elsif(/Predicted protein/i){
      last LINE;
    }
  }
  close(FH) or throw("FAILED to close ".$results.
                     "Fgenesh:parse_results");
  $self->create_transcripts;
}

1;
