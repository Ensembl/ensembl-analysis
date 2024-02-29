=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Genscan - 

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

=head1 METHODS

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
      if($pvalue !~ /\d+/){
        warning("Genscan has reported ".$pvalue." rather ".
                "than a number setting p value to 0");
        $pvalue = 0;
      }

      if($pvalue =~ /\+/)
      {
	  warning("Genscan has reported ".$pvalue." removing \+");
	  $pvalue =~ s/\+//;
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

  my $group;
  my $hash = {};
  PEP:while(<OUT>){
    chomp;
    if(/predicted peptide/i){
      next;
    }
    if(/^>/){
      my @values = split(/\|/, $_);
      ($group) = ($values[1] =~ /GENSCAN_predicted_peptide_(\d+)/); 
    }elsif(/(\S+)/ and defined $group) {
      $hash->{$group} .= $1;
    }
  }

  $self->peptides($hash);
  close(OUT) or throw("FAILED to close ".$results.
                      "Genscan:parse_results");
  $self->create_transcripts;
  $self->calculate_phases;
}


1;
