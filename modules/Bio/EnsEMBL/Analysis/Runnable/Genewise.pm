=pod

=head1 NAME
  
Bio::EnsEMBL::Pipeline::Runnable::GeneWise - aligns protein hits to genomic sequence

=head1 SYNOPSIS

    # Initialise the genewise module  
    my $genewise = new  Bio::EnsEMBL::Pipeline::Runnable::Genewise  (
        'genomic'  => $genomic,
	'protein'  => $protein,
        'memory'   => $memory);

   $genewise->run;
    
   my @genes = $genewise->output

=head1 DESCRIPTION


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::Genewise;

use strict;  
use vars   qw(@ISA);

# Object preamble

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::SeqIO;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Genewise;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($genomic, $protein, $memory,$reverse,$endbias, $gap, 
      $ext, $subs, $matrix, $verbose, $options) = rearrange([qw(DNA 
                                                                PROTEIN 
                                                                MEMORY 
                                                                REVERSE 
                                                                ENDBIAS 
                                                                GAP 
                                                                EXTENSION 
                                                                SUBS 
                                                                MATRIX
                                                                VERBOSE
                                                                OPTIONS)], @args);
  ####SETTING DEFAULTS####
  $self->program($GB_GENEWISE_EXE) unless($self->program);
  $self->memory($GB_GENEWISE_MEMORY);
  $self->gap($GB_GENEWISE_GAP);
  $self->extension($GB_GENEWISE_EXTENSION);
  $self->subs($GB_GENEWISE_SUBS);
  $self->matrix($GB_GENEWISE_MATRIX);
  $self->options($GB_GENEWISE_OPTIONS) unless($self->options);
  $self->verbose($GB_GENEWISE_VERBOSE);
  ########################
 
  $self->dna($genomic) || throw("No genomic sequence entered for blastwise");
  $self->protein($protein) || throw("No protein sequence entered for blastwise");
  throw("Must have a query sequence defined and it must be a Bio::EnsEMBL::Slice")
    if(!$self->query || !($self->query->isa("Bio::EnsEMBL::Slice")));

  $self->reverse($reverse)   if ($reverse);             
  $self->endbias($endbias)   if ($endbias);
  $self->memory ($memory);
  $self->gap($gap);
  $self->extension($ext);
  $self->subs($subs);
  $self->matrix($matrix);
  $self->verbose($verbose);

  return $self;
}


sub run{
  my ($self, $dir) = @_;
  $self->workdir($dir) if($dir);
  $self->checkdir();
  my $seqfile = $self->create_filename("dna", "fa");
  $self->write_seq_file($self->dna, $seqfile);
  $self->seqfile($seqfile);
  print "Have seqfile ".$seqfile."\n";
  $self->files_to_delete($seqfile);
  my $protfile = $self->create_filename("pep", "fa");
  $self->write_seq_file($self->protein, $protfile);
  $self->protfile($protfile);
  print "Have prot file ".$protfile."\n";
  $self->files_to_delete($protfile);
  $self->resultsfile($self->create_filename("results", "out"));
  $self->files_to_delete($self->resultsfile);
  $self->run_analysis;
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
  my $command = $self->program." ".$self->protfile." ".$self->seqfile.
    " -genesf -kbyte ".$self->memory." -ext ".$self->extension.
      " -gap ".$self->gap." -subs ".$self->subs." -m ".$self->matrix.
        " ".$self->options." > ".$self->resultsfile;
  logger_info($command);
  print $command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
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
  my @genesf_exons = $self->parse_genewise_output(\*GW);
  close(GW) or throw("Failed to close ".$results);
  my $transcript = $self->make_transcript(\@genesf_exons);
  if(defined $transcript){
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->biotype("genewise");
    $gene->add_Transcript($transcript);
    $self->output([$gene]);
  }else{
    logger_info("No valid transcript made, cannot make gene"); 
  }
}

sub parse_genewise_output {
  my ($self, $fh) = @_;

  my ($last_geneid, @exons_and_sfs);

  while(<$fh>) {
    print if($self->verbose);
    chomp;    

    my @l = split;

    next unless (defined $l[0] && ($l[0] eq 'Gene' || $l[0] eq 'Exon' || $l[0] eq 'Supporting'));

    if($l[0] eq 'Gene'){

      # flag a frameshift - ultimately we will do something clever here but for now ...
      # frameshift is here defined as two or more genes produced by Genewise from a single protein
      
      if (/^(Gene\s+\d+)$/){
        if(!defined($last_geneid)){
          $last_geneid = $1;
        } 
        elsif ($1 ne $last_geneid) {
          warning("frameshift!\n");
        }
      }
    }    
    elsif($l[0] eq 'Exon'){
      my $start  = $l[1];
      my $end    = $l[2];
      my $phase  = $l[4];
      my $strand;
      
      if ($l[1] == $l[2]) {
        # can't determine strand; defer
        $strand = 0;
      } else {
        if($l[1] > $l[2]){
          $strand = -1;
          $start  = $l[2];
          $end    = $l[1];
        } else {
          $strand = 1;
        }
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
      $exon->slice($self->query);

      push @exons_and_sfs, { exon => $exon,
                             sfs  => [],
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
      
      if($strand != $exons_and_sfs[-1]->{exon}->strand){
        warning("incompatible strands between exon and supporting feature - cannot add suppfeat\n");
        return;
      }
      
      if($pstart > $pend){
        warning("Protein start greater than end! Skipping this suppfeat\n");
        return;
      }
      
      my $fp = new Bio::EnsEMBL::FeaturePair();
      $fp->start   ($gstart);
      $fp->end     ($gend);
      $fp->strand  ($strand);
      $fp->seqname ($self->dna->id);
      $fp->hseqname($self->protein->id);
      $fp->hstart  ($pstart);
      $fp->hend    ($pend);
      $fp->hstrand (1);
      push(@{$exons_and_sfs[-1]->{sfs}}, $fp);

    }
  }

  my ($consensus_strand, @strandless_exons);

  foreach my $entry (@exons_and_sfs) {
    my ($exon, @sfs) = ($entry->{exon}, @{$entry->{sfs}});

    if ($exon->strand != 0 and not defined $consensus_strand) {
      $consensus_strand = $exon->strand;
    } else {
      push @strandless_exons, $exon;
    }
          
    if (@sfs) {
      my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@sfs);
      $align->seqname($self->query->seq_region_name);
      $align->slice($self->query);
      $align->score(100);
      $exon->add_supporting_features($align);    
    }
  }

  map { $_->strand($consensus_strand) } @strandless_exons;

  return map { $_->{exon} } @exons_and_sfs;
}

=head2 make_transcript

  Arg [1]   : $exons_ref, reference to array of Bio::EnsEMBL::Feature
  Function  : Turns array of exons into transcript & validates it.
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: 
  Caller    :
  Example   :

=cut

sub make_transcript{
  my ($self, $exonsref) = @_;
  my @exons = @$exonsref;

  my $transcript   = Bio::EnsEMBL::Transcript->new;
  my $translation  = Bio::EnsEMBL::Translation->new;
  $transcript->translation($translation);

  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
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

    my ($min_start, $max_end, $min_hstart, $max_hend, $total_hcoverage);
    foreach my $exon(@exons){

      my ($sf) = @{$exon->get_all_supporting_features};

      if (defined $sf) {
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

      $transcript->add_Exon($exon);

    }

    my $tsf = Bio::EnsEMBL::FeaturePair->new();
    $tsf->start($min_start);
    $tsf->end($max_end);
    $tsf->strand($transcript->strand);
    $tsf->seqname($self->dna->id);
    $tsf->hseqname($self->protein->id);
    $tsf->hstart($min_hstart);
    $tsf->hend    ($max_hend);
    $tsf->hstrand (1);   
    $tsf->hcoverage(100 * ($total_hcoverage / $self->protein->length));
    
    $transcript->add_supporting_features($tsf);
  }
  
  $transcript->slice($self->query);
  return $transcript;
}

# These all set/get or initializing methods


sub seqfile{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'seqfile'} = $arg;
  }
  return $self->{'seqfile'};
}

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

    return $self->{'_memory'} || 100000;
}

sub gap {
    my ($self,$arg) = @_;

    if(!$self->{'_gap'}){
      $self->{'_gap'} = undef;
    }
    if (defined($arg)) {
      $self->{'_gap'} = $arg;
    }

    return $self->{'_gap'} || 12;
}

sub extension {
    my ($self,$arg) = @_;

    if(!$self->{'_ext'}){
      $self->{'_ext'} = undef;
    }
    if (defined($arg)) {
      $self->{'_ext'} = $arg;
    }

    return $self->{'_ext'} || 2;
}

sub subs{
    my ($self,$arg) = @_;

    if(!$self->{'_subs'}){
      $self->{'_subs'} = undef;
    }
    if (defined($arg)) {
      $self->{'_subs'} = $arg;
    }

    return $self->{'_subs'} || 0.0000001;
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


sub dna {
  my ($self,$arg) = @_;
  
  if ($arg) {
    throw("Genomic sequence input is not a Bio::PrimarySeqI") unless
      ($arg->isa("Bio::PrimarySeqI"));
    $self->{'_genomic'} = $arg;
  }
  return $self->{'_genomic'};
}

sub protein {
    my ($self,$arg) = @_;

    if (defined($arg)) {
			throw("[$arg] is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI") unless
				($arg->isa("Bio::SeqI") || 
				 $arg->isa("Bio::Seq")  || 
				 $arg->isa("Bio::PrimarySeqI"));
			
			$self->{'_protein'} = $arg;
    }
    return $self->{'_protein'};
}

sub slice {
    my ($self,$arg) = @_;

    if (defined($arg)) {
			throw("slice sequence input is not a Bio::EnsEMBL::Slice") unless
				($arg->isa("Bio::EnsEMBL::Slice"));
			
			$self->{'_slice'} = $arg;
    }
    return $self->{'_slice'};
}



sub check_environment {
  my($self,$genefile) = @_;

  if (! -d $ENV{WISECONFIGDIR}) {
    throw("No WISECONFIGDIR ["  . $ENV{WISECONFIGDIR} . "]");
  }

}

1;
