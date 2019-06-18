
package Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus;
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

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus->new
  (
   -query => $slice,
   -program => 'augustus',
   -species => 'human',
   -analysis => $analysis,
  );
$runnable->run;
my @predictions = @{$runnable->output};


=head1 DESCRIPTION

Wrapper to run the augustus gene predictor and then parse the results
into prediction transcripts

this is a bare bones module which inherits most of its functionality from
Bio::EnsEMBL::Analysis::Runnable::Genscan

=head1 CONTACT

Post questions to : anacode-people@sanger.ac.uk

=cut

use vars qw(@ISA);
use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio);

=head2 new

  Returntype : Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus
=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($species) = rearrange(['SPECIES'], @args);
  $self->species($species);
  
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->program('/software/anacode/bin/augustus') if(!$self->program);
  $self->species('human') if(!$self->species);
  ######################

  return $self;
}

sub species {
  my ($self,$species) = @_;
  if($species){
    $self->{'species'} = $species;
  }
  return $self->{'species'};
}

sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable Augustus::run_analysis ")
    unless($program && -x $program);

  my $command = $program." --species=".$self->species." ";
  $command .= $self->options." " if($self->options);
  $command .= $self->queryfile." > ".$self->resultsfile;
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
}

## Predicted genes for sequence number 1 on both strands
#### gene g1
#seqname   source     feature    start   end   score   strand   frame    transcript and gene name
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        gene    108106  121071  1       -       .       g1
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        transcript      111161  121071  0.5     -       .       g1.t1
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        stop_codon      111161  111163  .       -       0       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        terminal        111161  111357  0.56    -       2       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        intron  111358  113826  1       -       .       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        internal        113827  114097  1       -       0       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        intron  114098  120939  0.99    -       .       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        initial 120940  121071  0.92    -       0       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        start_codon     121069  121071  .       -       0       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        CDS     111164  111357  0.56    -       2       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        CDS     113827  114097  1       -       0       transcript_id "g1.t1"; gene_id "g1";
#contig::BX901927.7.1.152042:1:152042:1 contig BX901927.7.1.152042       AUGUSTUS        CDS     120940  121071  0.92    -       0       transcript_id "g1.t1"; gene_id "g1";
## protein sequence = [MFFQFGPSIEQQASVMLNIMEEYDWYIFSIVTTYYPGHQDFVNRIRSTVDNSFVGWELEEVLLLDMSVDDGDSKIQNQ
## MKKLQSPVILLYCTKEEATTIFEVAHSVGLTGYGYTWIVPSLVAGDTDNVPNVFPTGLISVSYDEWDYGLEARVRDAVAIIAMATSTMMLDRGPHTLL
## KSGCHGAPDKKGSKSGNPNEVLR]


=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus
  Arg [2]   : string, filename
  Function  : parse the output from Augustus into prediction transcripts
  Returntype: none
  Exceptions: throws if cannot open or close results file
  Example   :

=cut

sub parse_results{
  my ($self, $results) = @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  open(OUT, "<".$results) or throw("FAILED to open ".$results."Augustus:parse_results");

  my $genes;
  my $verbose = 0;
  # parse output
  while (<OUT>) {
    chomp;
    if(/^#/) { next; }
    my @element = split /\t/;
    if($element[2] && $element[2] eq 'CDS')
    {
		my ($transcript_id,$gene_id) = $element[8] =~ /transcript_id "(.*)"; gene_id "(.*)";/;
		my @exon = @element[3..7]; #($start,$end,$score,$strand,$phase)
		$genes->{$gene_id}->{$transcript_id} ||= [];
		push @{$genes->{$gene_id}->{$transcript_id}}, \@exon;
    }
  }

  # create exon objects
  foreach my $gene (keys %$genes) {
  	print STDOUT "Gene $gene\n" if($verbose);
    foreach my $transcript (keys %{$genes->{$gene}}) {
    	print STDOUT "\tTranscript $transcript\n" if($verbose);
    	my $exons = $genes->{$gene}->{$transcript};
    	my $strand = $exons->[0]->[3];
    	if($strand eq '+') {
    		$strand = 1;
    		@$exons = sort {$a->[0] <=> $b->[0]} @$exons;
    	}else{
    		$strand = -1;
    		@$exons = reverse sort {$a->[0] <=> $b->[0]} @$exons;
    	}
    	my $exon_num=1;
    	foreach my $exon (@$exons) {
    		my $start	= $exon->[0];
    		my $end		= $exon->[1];
    		my $score	= $exon->[2];
    		my $phase	= $exon->[4];
    		my $e_name	= $transcript.".e".$exon_num; # g1.t1.e1
    		my $e = $self->feature_factory->create_prediction_exon
             ($start, $end, $strand, $score, '0', $phase, $e_name, $self->query, $self->analysis);
             print "\t\t\t".join("\t",$e_name,$start, $end, $strand, $score, $phase,"\n" ) if($verbose);
           	$self->exon_groups($transcript, $e);
           	$exon_num++;
    	}
    }
  }

  # create transcript objects
  $self->create_transcripts();
  close(OUT) or throw("FAILED to close ".$results."Augustus:parse_results");

}

1;
