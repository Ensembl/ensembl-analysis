=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME
  
Bio::EnsEMBL::Analysis::Runnable::GeneWise - aligns protein hits to genomic sequence

=head1 SYNOPSIS

    # Initialise the genewise module  
    my $genewise = new  Bio::EnsEMBL::Analysis::Runnable::Genewise(
        'genomic'  => $genomic,
	'protein'  => $protein,
        'memory'   => $memory);

   $genewise->run;
    
   my @genes = $genewise->output

=head1 DESCRIPTION


=head1 METHODS

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::Genewise;

use warnings ;
use strict;  
use vars   qw(@ISA);

# Object preamble

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::SeqIO;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(create_Exon);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils 
  qw(create_feature create_feature_from_gapped_pieces);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(Transcript_info);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  #print "IN CONSTRUCTOR\n";
  my ($protein, $memory, $reverse,$endbias, $gap, 
      $ext, $subs, $matrix, $verbose, $hmm, $splice_model) = rearrange([qw(
                                                            PROTEIN 
                                                            MEMORY 
                                                            REVERSE 
                                                            ENDBIAS 
                                                            GAP 
                                                            EXTENSION 
                                                            SUBS 
                                                            MATRIX
                                                            VERBOSE
                                                            HMM
																														SPLICE_MODEL)], @args);
  
	####SETTING DEFAULTS####
  $self->program('genewise') unless($self->program);
  $self->memory(100000);
  $self->gap(12);
  $self->extension(2);
  $self->subs(0.0000001);
  $self->matrix('BLOSUM62.bla');
  $self->options('-quiet') unless($self->options);
  #$self->endbias(1);
	#$self->splice_model(0);
  ########################
 
  $self->protein($protein);
  $self->hmm($hmm);
  throw("Must have a hmm or a protein") if(!$self->protein && !$self->hmm);
  throw("Must not have both a hmm ".$self->hmm." and a protein ".$self->protein) 
    if($self->protein && $self->hmm);
  throw("Must have a query sequence defined")
    if(!$self->query);
  $self->reverse($reverse);             
  $self->endbias($endbias) if(defined $endbias);
  $self->memory ($memory);
  $self->gap($gap);
  $self->extension($ext);
  $self->subs($subs);
  $self->matrix($matrix);
	$self->splice_model($splice_model) if($splice_model);
  $self->verbose($verbose);
  #print "RUNNING On ".$self->query->id." protein ".$self->protein->id."\n";
  return $self;
}


sub run{
  my ($self, $dir) = @_;
  $self->workdir($dir) if($dir);
  $self->checkdir();
  #throw("Sequence appears to have no context") unless($self->query->seq =~ /[CATG]{3}/);
  my $seqfile = $self->write_seq_file;
  $self->files_to_delete($seqfile);
  if($self->protein){
    my $protfile = $self->create_filename("pep", "fa");
    $self->write_seq_file($self->protein, $protfile);
    $self->protfile($protfile);
    $self->files_to_delete($protfile);
  }
  $self->resultsfile($self->create_filename("results", "out"));
  $self->files_to_delete($self->resultsfile);
  $self->run_analysis; 
   # Warning ! 
   # you can't run Genewise.pm straight away, like Bio::EnsEMBL::Genewise->run() , as it has dependecies 
   # about Featuretypes / FeaturePairs - run BlastMiniGenewise or Minigenewise first ! 
   #
  $self->parse_results;
  $self->delete_files;
}


sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  $self->check_environment;
  throw($program." is not executable Runnable::Genewise::run_analysis ") 
    unless($program && -x $program);
  my $options = " -genesf -kbyte ".$self->memory." -ext ".
    $self->extension." -gap ".$self->gap." -subs ".$self->subs." -m ".
      $self->matrix." ".$self->options;

  if($self->endbias == 1){
    $options =~ s/-init endbias//;
    $options =~ s/-splice flat//;
    $options =~ s/-splice_gtag//;
    
		if (!$self->splice_model || $self->splice_model == 0){
		  print "USING STANDARD SPLICE MODEL\n";
		  $options .= " -init endbias -splice_gtag ";
    }elsif ($self->splice_model == 1){
		  print "USING ALTERNATIVE SPLICE MODEL\n";
      $options .= " -nosplice_gtag ";
		}
  }
  if (($self->reverse) && $self->reverse == 1) {
    $options .= " -trev ";
  }
  my $command = $self->program." ";
  if($self->protein){
    $command .= $self->protfile." ".$self->queryfile." ".$options;
  }elsif($self->hmm){
    $command .= $options." -hmmer ".$self->hmm." ".$self->queryfile;
  }else{
    throw("Seem to be without either a hmm or a protein sequence");
  }
  $command .= " > ".$self->resultsfile;
  logger_info($command);
  print "MMG COMMAND for ".$self->protein->id." ".$self->query->id." ".
  $command."\n";
  system($command) == 0 or throw("FAILED to run protein " . $self->protein->id . " : $command "); 
}




=head2 parse_genewise_output

  Arg [1]   : $columns_ref, reference to array of strings
  Arg [2]   : $curr_gene_ref, reference to scalar
  Arg [3]   : $curr_exon_ref, reference to Bio::EnsEMBL::Feature
  Arg [4]   : $exons_ref, reference to array of Bio::EnsEMBL::Feature
  Function  : Parses genewise output lines and adds exons/supporting features to exons_ref
  Returntype: 
  Exceptions: 
  Caller    :
  Example   :

=cut


sub parse_results{
  my ($self, $results) =  @_;
  if(!$results){
    $results = $self->resultsfile;
  }
  open(GW, $results) or throw("Failed to open ".$results);
  my $genesf_exons = $self->parse_genewise_output(\*GW);
  close(GW) or throw("Failed to close ".$results);
  my $transcript = $self->make_transcript($genesf_exons);
  if(defined $transcript){
    $self->output([$transcript]);
  }else{
    warning("No valid transcript made, cannot make gene for ".$self->protein->id); 
  }
}


sub parse_genewise_output {
  my ($self, $fh) = @_;

  my ($hit_name, @genes);

  while(<$fh>) {
    $self->verbose and print;
    chomp;

    if (/^Query\s+\S+:\s+(\S+)/) {
      $hit_name = $1;
      next;
    }

    my @l = split;

    next unless (defined $l[0] && ($l[0] eq 'Gene' || $l[0] eq 'Exon' || $l[0] eq 'Supporting'));

    if($l[0] eq 'Gene' and /^Gene\s+\d+$/){
      push @genes, {
        exons => [],
      };
    }    
    elsif($l[0] eq 'Exon'){
      my $start  = $l[1];
      my $end    = $l[2];
      my $phase  = $l[4];
      my $strand = 1;
      
      if($l[1] > $l[2]){
        $strand = -1;
        $start  = $l[2];
        $end    = $l[1];
      }
      
      my $exon_length = $end - $start + 1;
      
      # end phase is the number of bases at the end of the exon which do not
      # fall in a codon and it coincides with the phase of the following exon.
      
      my $end_phase   = ( $exon_length + $phase ) %3;
      
      my $exon = new Bio::EnsEMBL::Exon;

      $exon->start    ($start);
      $exon->end      ($end);
      $exon->strand   ($strand);
      $exon->phase    ($phase);
      $exon->end_phase($end_phase);      

      if (not exists $genes[-1]->{strand}) {
        $genes[-1]->{strand} = $exon->strand;
      }

      push @{$genes[-1]->{exons}}, {
        ex => $exon,
        sf => [],
      };

    }    
    elsif($l[0] eq 'Supporting') {
      
      my $gstart = $l[1];
      my $gend   = $l[2];
      my $pstart = $l[3];
      my $pend   = $l[4];
      
      my $strand = 1;
      
      if ($gstart > $gend){
        $gstart = $l[2];
        $gend   = $l[1];
        $strand = -1;
      }
      
      if($strand != $genes[-1]->{strand}) {
        warning("incompatible strands between exon and supporting feature - cannot add suppfeat\n");
        return;
      }
      
      if($pstart > $pend){
        warning("Protein start greater than end! Skipping this suppfeat\n");
        return;
      }
      
      my $fp = create_feature("Bio::EnsEMBL::FeaturePair", $gstart, $gend, $strand,
                              $pstart, $pend, 1, undef, undef, undef, $hit_name,
                              undef, $self->analysis);

      if (not exists $genes[-1]->{start}) {
        $genes[-1]->{start}  = $fp->start;
        $genes[-1]->{end}    = $fp->end;
        $genes[-1]->{hstart} = $fp->hstart;
        $genes[-1]->{hend} = $fp->hend;
      } else {
        if ($genes[-1]->{start} > $fp->start) {
          $genes[-1]->{start} = $fp->start;
        }
        if ($genes[-1]->{end} < $fp->end) {
          $genes[-1]->{end} = $fp->end;
        }
        if ($genes[-1]->{hstart} > $fp->hstart) {
          $genes[-1]->{hstart} = $fp->hstart;
        }
        if ($genes[-1]->{hend} < $fp->hend) {
          $genes[-1]->{hend} = $fp->hend;
        }
      }

      push @{$genes[-1]->{exons}->[-1]->{sf}}, $fp;
    }
  }

  my @kept_genes;
  foreach my $gene (@genes) {
    
    if (not @kept_genes) {
      push @kept_genes, $gene;
    } else {
      # only keep this gene if all exons are "consistent" with kept
      # exons and supporting features so far

      my @test_genes = (@kept_genes, $gene);
      @test_genes = sort { $a->{hstart} <=> $b->{hstart} } @test_genes;

      my $bad = 0;

      for(my $i=1; not $bad and $i < @test_genes; $i++) {
        my $this = $test_genes[$i];
        my $prev = $test_genes[$i-1];

        if ($this->{strand} != $prev->{strand}) {
          $bad = 1;
        } elsif ($this->{hstart} <= $prev->{hend}) {
          $bad = 1;
        } elsif ($this->{strand} > 0 and $this->{start} <= $prev->{end}) {
          $bad = 1;
        } elsif ($this->{strand} < 0 and $this->{end} >= $prev->{start}) {
          $bad = 1;
        }
      }
      
      if (not $bad) {
        push @kept_genes, $gene;
      }
    }
  }
  
  my @exons;
  foreach my $g (@kept_genes) {
    foreach my $entry (@{$g->{exons}}) {
      my ($exon, @sfs) = ($entry->{ex}, @{$entry->{sf}});
          
      if (@sfs) {
        my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@sfs, -align_type => 'ensembl');
        $align->seqname($self->query->id);
        $align->score(100);
        $exon->add_supporting_features($align);    
      }
      push @exons, $exon;
    }
  }
    
  return \@exons;
}
=head2 make_transcript

  Arg [1]   : $exons_ref, reference to array of Bio::EnsEMBL::Feature
  Function  : Turns array of exons into transcript & validates it.
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: 
  Caller    :
  Example   :

=cut

################################
sub make_transcript{
  my ($self, $exonsref) = @_;
  my @exons = @$exonsref;

  my $transcript   = Bio::EnsEMBL::Transcript->new;
  my $translation  = Bio::EnsEMBL::Translation->new;
  $transcript->translation($translation);

  if ($#exons < 0) {
    warning("Odd.  No exons found\n");
    return undef;
  }

  else {

    if ($exons[0]->strand == -1) {
      @exons = sort {$b->start <=> $a->start} @exons;
    } else {
      @exons = sort {$a->start <=> $b->start} @exons;
    }

    $translation->start_Exon($exons[0]);
    $translation->end_Exon  ($exons[$#exons]);

    # phase is relative to the 5' end of the transcript (start translation)
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }

    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);

    my ($min_start, $max_end, $min_hstart, $max_hend, $total_hcoverage, $hit_name);
    foreach my $exon(@exons){

      my ($sf) = @{$exon->get_all_supporting_features};

      if (defined $sf) {
        $hit_name = $sf->hseqname;
        if (not defined $min_start or $sf->start < $min_start) {
          $min_start = $sf->start;
        }
        if (not defined $max_end or $sf->hend > $max_end) {
          $max_end = $sf->end;
        }
        if (not defined $min_hstart or $sf->hstart < $min_hstart) {
          $min_hstart = $sf->hstart;
        }
        if (not defined $max_hend or $sf->hend > $max_hend) {
          $max_hend = $sf->hend;
        }
        $total_hcoverage += $sf->hend - $sf->hstart + 1;
      }
      else {
        warning('No supporting evidence for '.$exon->start.' '.$exon->end.' '.$exon->strand."\n");
      }


      eval {$transcript->add_Exon($exon)};
      if ($@) {
        warning($@);
        return;
      }
    }
    my $tsf = Bio::EnsEMBL::FeaturePair->new();
    $tsf->start($min_start);
    $tsf->end($max_end);
    $tsf->strand($transcript->strand);
    $tsf->seqname($self->query->id);
    $tsf->hseqname($hit_name);
    $tsf->hstart($min_hstart);
    $tsf->hend    ($max_hend);
    $tsf->hstrand (1);   
    $tsf->hcoverage(100 * ($total_hcoverage / $self->target_length));
    
    $transcript->add_supporting_features($tsf);
  }
  
  return $transcript;
}
# These all set/get or initializing methods




sub protfile{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'protfile'} = $arg;
  }
  return $self->{'protfile'};
}

sub reverse {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_reverse'} = $arg;
  }
  return $self->{'_reverse'};
}

sub endbias {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_endbias'} = $arg;
  }
  
  if (!defined($self->{'_endbias'})) {
    $self->{'_endbias'} = 0;
  }
  
  return $self->{'_endbias'};
}

sub memory {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_memory'} = $arg;
  }
  
  return $self->{'_memory'};
}

sub gap {
  my ($self,$arg) = @_;
  
  if(!$self->{'_gap'}){
    $self->{'_gap'} = undef;
  }
  if (defined($arg)) {
    $self->{'_gap'} = $arg;
  }
  
  return $self->{'_gap'};
}

sub extension {
  my ($self,$arg) = @_;
  
  if(!$self->{'_ext'}){
    $self->{'_ext'} = undef;
  }
  if (defined($arg)) {
    $self->{'_ext'} = $arg;
  }
  
  return $self->{'_ext'};
}

sub subs{
  my ($self,$arg) = @_;
  
  if(!$self->{'_subs'}){
    $self->{'_subs'} = undef;
  }
  if (defined($arg)) {
    $self->{'_subs'} = $arg;
  }
  
  return $self->{'_subs'};
}

sub matrix{
  my ($self,$arg) = @_;
  
  if(!$self->{'_matrix'}){
    $self->{'_matrix'} = undef;
  }
  if (defined($arg)) {
    $self->{'_matrix'} = $arg;
  }
  
    return $self->{'_matrix'};
}

sub splice_model{
  my ($self, $arg) = @_;

	if(!$self->{'_splice_model'}){
	  $self->{'_splice_model'} = undef;
	}
	if ($arg){
    $self->{'_splice_model'} = $arg;
	}
	return $self->{'_splice_model'};
}
sub verbose {
  my ($self,$arg) = @_;
  
  if(not exists $self->{'_verbose'}){
    $self->{'_verbose'} = 0;
  }
  if (defined($arg)) {
    $self->{'_verbose'} = $arg;
  }
  
  return $self->{'_verbose'};
}

sub protein {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    throw("[$arg] is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI") unless
      ($arg->isa("Bio::SeqI") || 
       $arg->isa("Bio::Seq")  || 
       $arg->isa("Bio::PrimarySeqI"));
    
    $self->{'_protein'} = $arg;
    $self->target_length($arg->length);
  }
  return $self->{'_protein'};
}

sub hmm {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    open(HMM, $arg) or throw("[$arg] is not a readable file");
    my $model_len;
    while(<HMM>) {
      /^LENG\s+(\d+)$/ and do {
        $model_len = $1;
        last;
      }
    }
    close(HMM);
    throw("[$arg] is not a HMM file in HMMER format") if not defined $model_len;
    $self->target_length($model_len);
    
    $self->{'_hmm'} = $arg;
  }
  return $self->{'_hmm'};
}


sub target_length {
  my ($self, $arg) = @_;
  
  if (defined $arg) {
    $self->{'_target_length'} = $arg;
  }
  return $self->{'_target_length'};
}

sub check_environment {
  my($self,$genefile) = @_;
  if (! -d $ENV{WISECONFIGDIR}) {
    throw("No WISECONFIGDIR ["  . $ENV{WISECONFIGDIR} . "]");
  }
}

1;
