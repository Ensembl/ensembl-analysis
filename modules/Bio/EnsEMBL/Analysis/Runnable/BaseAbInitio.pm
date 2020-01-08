=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

use parent qw(Bio::EnsEMBL::Analysis::Runnable);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($matrix) = rearrange(['MATRIX'], @args);
  $self->matrix($matrix);

  return $self;
}




#containters



=head2 matrix

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, matrix file
  Function  : container for path to matrix file, when passed a filename it
  will use the find_file method to check for its existance
  Returntype: string
  Exceptions: 
  Example   : 

=cut


sub matrix{
  my ($self, $matrix) = @_;
  if($matrix){
    my $file = $self->find_file($matrix);
    $self->{'matrix'} = $file;
  }
  return $self->{'matrix'};
}


=head2 peptides

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, peptide sequence from genscan results
  Function  : container for the peptides sequences from genscan results
  Returntype: arrayref
  Exceptions: 
  Example   : 

=cut



sub peptides{
  my ($self, $peptides) = @_;
  if($peptides){
    $self->{'peptides'} = $peptides;
  }
  return $self->{'peptides'};
}


=head2 exon_groups

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, group name
  Arg [3]   : Bio::EnsEMBL::PredictionExon
  Function  : stores exons in arrayrefs in a hash keyed on group names
  Returntype: hashref
  Exceptions: throws if passed an exon which isnt a 
  Bio::EnsEMBL::PredictionExon
  Example   : 

=cut


sub exon_groups{
  my ($self, $group, $exon) = @_;
  if(!$self->{'exon_groups'}){
    $self->{'exon_groups'} = {};
  }
  if($group && $exon){   
    if(!$exon->isa('Bio::EnsEMBL::PredictionExon')){
      throw("must be passed an exon object ".
            "Bio::EnsEMBL::PredictionExon not an $exon ".
            "BaseAbInitio:exon_groups");
    }
    if(!$self->{'exon_groups'}->{$group}){
      $self->{'exon_groups'}->{$group} = [];
    }
    push(@{$self->{'exon_groups'}->{$group}}, $exon);
  }
  return $self->{'exon_groups'};
}


#utility methods


=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, program name
  Function  : create and open a commandline for one
  of the ab initio gene finding programs
  Returntype: none
  Exceptions: throws if the program in not executable or the system
  command fails to execute
  Example   : 

=cut


sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  throw($program.' is not executable '.ref($self).'::run_analysis ')
    unless($program && -x $program);
  my $command = $program." ".$self->matrix." ".$self->queryfile." > ".
    $self->resultsfile;
  print "Running analysis ".$command."\n";
  execute_with_wait($command);
}

=head2 create_transcripts

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Function  : use the groups of exons to create prediction transcripts
  Returntype: none
  Exceptions: none
  Example   : 

=cut



sub create_transcripts{
  my ($self) = @_;
  my @transcripts;
  my $ff = $self->feature_factory;
  my %exon_groups = %{$self->exon_groups};
  foreach my $group(keys(%exon_groups)){
    my @exons = @{$exon_groups{$group}};
    my $transcript = $ff->create_prediction_transcript(\@exons, $self->query,
                                                       $self->analysis);
    $transcript->seqname($group);
    push(@transcripts, $transcript);
  }
  $self->output(\@transcripts);
}

=head2 calculate_phases

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Function  : works out which phase to make the exons to get 
  a complete cds
  Returntype: none
  Exceptions: throws if it cant find a translation 
  for a transcript
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

  $genscan_pep =~ s/\*$//;

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

  my $x = 0;
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
