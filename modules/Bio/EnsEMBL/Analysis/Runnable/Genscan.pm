# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Genscan
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Genscan

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::Genscan->new(
      -query => $slice,
      -program => 'genscan',
      -matrix => 'HumanIso.smat',
     );
  $runnable->run;
  my @predictions = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run the genscan gene predictor and then parse the results
into prediction transcripts

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Genscan;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genscan
  Arg [2]   : string, matrix file
  Function  : create a Bio::EnsEMBL::Analysis::Runnable::Genscan runnable
  Returntype: Bio::EnsEMBL::Analysis::Runnable::Genscan
  Exceptions: none
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($matrix) = rearrange(['MATRIX'], @args);

  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->program('genscan') if(!$self->program);
  $self->matrix('HumanIso.smat') if(!$self->matrix);
  ######################

  return $self;
}





=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genscan
  Arg [2]   : string, resultsfile name
  Function  : parse the results file into prediction exons then
  collate them into prediction transcripts and calculate their
  phases
  Returntype: none 
  Exceptions: throws if cant open or close results file or the parsing
  doesnt work
  Example   : 

=cut


sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  open(OUT, "<".$results) or throw("FAILED to open ".$results.
                                   "Genscan:parse_results");
  my $ff = $self->feature_factory;
 LINE:while(<OUT>){
    chomp;
    if(m|NO EXONS/GENES PREDICTED IN SEQUENCE|i){
      print "No genes predicted\n";
      last LINE;
    }
    # 3.00 Prom +  29940  29979   40                              -7.05
    # 3.01 Init +  33401  33538  138	 1  0  100   65	  110 0.621  10.09
    # 3.02 Intr +  35584  35734  151	 0  1	-7   81	  116 0.301   0.21
    # 3.03 Intr +  35979  36114  136	 1  1	45   44	   95 0.260  -0.49
    # 3.04 Intr +  37050  37124   75	 0  0	97   92	   12 0.130   0.21
    # 3.05 Intr +  43698  44101  404	 0  2  100   85	   74 0.522   1.45
    # 3.06 Term +  44406  44596  191	 1  2  112   48	   99 0.422   4.93
    # 3.07 PlyA +  47675  47680    6                               1.05
    # 6.01 Sngl - 100426  99728  699	 1  0	76   36	  273 0.699  16.88
    if(/init|intr|term|sngl/i){
      my @elements = split;
      if(@elements != 13){
        throw("Can't parse ".$_." splits into wrong number of elements ".
              "Genscan:parse_results");
      }
      my ($name, $strand, $start, $end, $pvalue, $score)
        = @elements[0, 2, 3, 4, 11, 12];
      
      if($strand eq '+'){
        $strand = 1;
      }else{
        $strand = -1;
        my $temp_start = $end;
        $end = $start;
        $start = $temp_start;
      }
      my ($group, $exon_name) = split(/\./, $name);
      
      my $exon = $ff->create_prediction_exon($start, $end, $strand, 
                                             $score, $pvalue, 0, 
                                             $name, $self->query, 
                                             $self->analysis);
      
      $self->exon_groups($group, $exon);
    }elsif(/predicted peptide/i){
      last LINE;
    }else{
      next LINE;
    }
  }
  my $peptide;
  my $group;
  my $hash = {};
 PEP:while(<OUT>){
    chomp;
    if(/predicted peptide/i){
      next;
    }
    if(/^>/){
      if($peptide){
        $peptide =~ s/\s+//;
        $hash->{$group} = $peptide;
      }
      $peptide = undef;
      my @values = split(/\|/, $_);
      $values[1] =~ /GENSCAN_predicted_peptide_(\d+)/;
      $group = $1;
    }else{
      $peptide .= $_;
    }
  }
  $hash->{$group} = $peptide;
  $self->peptides($hash);
  close(OUT) or throw("FAILED to close ".$results.
                      "Genscan:parse_results");
  $self->create_transcripts;
  $self->calculate_phases;
}



=head2 calculate_phases

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genscan
  Function  : starting with phase 0 create a transcript whose start
  phase is either 0, 1, 2 and keep the one whose translation matches the
  translation from genscan
  Returntype: none
  Exceptions: 
  Example   : 

=cut



sub calculate_phases{
  my ($self) = @_;
  my @phases = (0, 1, 2);
  my @output;
  my $ff = $self->feature_factory;
  my $peptides = $self->peptides;
 TRANS:foreach my $trans(@{$self->output}){
    my @exons = @{$trans->get_all_Exons};
    foreach my $phase(@phases){
      my @temp_exons = @{$self->set_phases($phase, \@exons)};
      my $new = $ff->create_prediction_transcript(\@temp_exons, 
                                                  $self->query,
                                                  $self->analysis);
      my $pep = $new->translate->seq;
      my $peptide = $peptides->{$trans->seqname};
      my ($ensembl, $genscan) = $self->subsitute_x_codes($pep, 
                                                         $peptide); 
      $ensembl =~ s/^x//i;
      $ensembl =~ s/x$//i;
      $genscan =~ s/^x//i;
      $genscan =~ s/x$//i;
      if ($ensembl =~ /$genscan/){
        push(@output, $new);
        next TRANS;
      }
    }
    throw("Failed to find translation for ".$trans." ".$exons[0]->seqname)
  }
  $self->clean_output;
  $self->output(\@output);
}



=head2 set_phases

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genscan
  Arg [2]   : number, start phase,
  Arg [3]   : arrayref of Bio::EnsEMBL::PredictionExons
  Function  : starting with the given phase sets of the phase of
  each exon on the basis of the end phase of the last. This is done
  after ordering them on the basis of there strand
  Returntype: arrayref of
  Exceptions: 
  Example   : 

=cut


sub set_phases{
  my ($self, $start_phase, $exons) = @_;
  if(@$exons == 0){
    throw("Can't set phases if have no exons ".$exons." ".@$exons);
  }
  if ($exons->[0]->strand == 1) {
    @$exons = sort {$a->start <=> $b->start} @$exons;
  } else {
    @$exons = sort {$b->start <=> $a->start} @$exons;
  }
  foreach my $e(@$exons){
    $e->phase($start_phase);
    $start_phase = ($e->phase + $e->length)%3;
  }
  return $exons;
}



=head2 subsitute_x_codes

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Genscan
  Arg [2]   : string, ensembl produced peptides
  Arg [3]   : string, genscan predicted peptide
  Function  : makes sure x's and length's of peps are the same
  to allow fair comparison
  Returntype: string, string
  Exceptions: 
  Example   : 

=cut



sub subsitute_x_codes{
  my ($self, $ensembl_pep, $genscan_pep) = @_;
  my $x = 0;
  my $ens_len = length($ensembl_pep);
  my $gen_len = length($genscan_pep);
  if($ens_len == ($gen_len+1)){
    chop($ensembl_pep);
  }
  if($gen_len == ($ens_len+1)){
    chop($genscan_pep);
  }
  $ens_len = length($ensembl_pep);
  $gen_len = length($genscan_pep);
  while (($x = index($ensembl_pep, 'X', $x)) != -1) {
		substr($genscan_pep, $x, 1) = 'X'
      if length($genscan_pep) >= length($ensembl_pep);
		$x++;
  }
  
  $x = 0;
  while (($x = index($genscan_pep, 'X', $x)) != -1) {
		substr($ensembl_pep, $x, 1) = 'X'
      if length($ensembl_pep) >= length($genscan_pep);
		$x++;
  }
  return $ensembl_pep, $genscan_pep;
}

1;
