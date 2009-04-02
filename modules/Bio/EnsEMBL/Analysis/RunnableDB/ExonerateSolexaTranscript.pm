
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaTranscript




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaTranscript->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
  );

$runnableDB->fetch_input();
$runnableDB->run();
$runnableDB->write_output(); #writes to DB

=head1 DESCRIPTION

Extends Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa to allow
reads to be aligned against transcript sequences and then projected
onto the genome

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaTranscript;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::ExonerateSolexa;
use Bio::EnsEMBL::FeaturePair;
use vars qw(@ISA);

@ISA =  ("Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa");


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  return $self;
}


sub fetch_input {
  my ($self) = @_;
  # do all the normal stuff
  $self->SUPER::fetch_input();
  # then get all the transcripts and exons
  my $trans_db = $self->get_dbadaptor($self->TRANSDB);
  my $trans_adaptor = $trans_db->get_TranscriptAdaptor;
  my %trans_by_id;
  my $biotype = $self->BIOTYPE;
  my @trans;
  # fetch genes, transcripts
  if ( $biotype ){
    push @trans , @{$trans_adaptor->generic_fetch("biotype = \"$biotype\"")};
  } else {
    push @trans , @{$trans_adaptor->fetch_all(undef,undef,undef)};
  }
  
  foreach my $trans ( @trans ) {
    $trans_by_id{$trans->display_id} = $trans;
  }
  $self->transcripts(\%trans_by_id);
  return ;
}


sub filter_solexa {
  my ($self,$features_ref) = @_;
  my @features = @$features_ref;
  my @filtered_features;
  #�allow no more than MISSMATCH missmatches and
  foreach my $feat ( @features ) {
   # check missmatches
    if ( $self->MISSMATCH ) {
      my $aligned_length = abs($feat->hend - $feat->hstart) +1;
      #print $feat->hseqname ." ";
      #print $feat->percent_id . " " ;
      #print " aligned length = $aligned_length ";
      my $matches = $aligned_length *  $feat->percent_id / 100;
      #print " matches $matches  ";
      my $missmatches = ( $aligned_length - $matches) / $aligned_length * 100;
      #print " missmatches $missmatches \n";
      next if $missmatches > $self->MISSMATCH;
      #print "accepting\n";
      
    }
    push @filtered_features, $feat;
  }
  return \@filtered_features;
}

sub run {
  my ($self) = @_;

  throw("Can't run - no runnable objects") unless ( $self->runnable );

  my ($runnable) = @{$self->runnable};

  $runnable->run;
  my $features = $runnable->output;

  if ($self->MISSMATCH) {
    $features = $self->filter_solexa($features);
  }

  # Pair features together if they come from paired end reads
  if ( $self->PAIREDEND ) {
    $features = $self->pair_features($features);
  }

  my $genomic_features = $self->process_features($features);

  $self->output($genomic_features);
}

=head2 process_features

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Uses cdna2genome to convert the ungapped align 
              features into genomic coords
  Returntype: 1
  Exceptions: Throws if the transcript the read has aligned to is
              not stored in the transcript hash

=cut

sub process_features {
  my ( $self, $flist ) = @_;

  # first do all the standard processing, adding a slice and analysis etc.
  # unless we are projecting in which case we dont really nead a slice 
  unless ($self->PROJECT) {
    $self->SUPER::process_features($flist);
  }

  my $slice_adaptor = $self->db->get_SliceAdaptor;
  my @dafs;

FEATURE:  foreach my $f (@$flist) {
    
    my $trans_id;
    my $accept = 0;
    
    $accept = 1 unless  $self->INTRON_OVERLAP or $self->INTRON_MODELS ;

    if ($f->seqname =~ /\S+\:\S+\:(\S+)\:\d+\:\d+:\d+/ ) {
      $trans_id = $1;
    } else {
      $trans_id = $f->seqname;
    }

    if ( not exists $self->transcripts->{$trans_id} ) {
      $self->throw("Transcript $trans_id not found\n");
    }

    my $trans = $self->transcripts->{$trans_id};

    if ( $self->INTRON_OVERLAP ) {
      # only consider introns overlapping exons by at least x base pairs 
      # on both sides of the intron
      # make a hash with the positions of the exon boundaries
      # check first if you have already stored it though, no need to 
      # make more than is needed
	unless ( $self->exon_starts($trans->display_id) && $self->exon_ends($trans->display_id) ){
	my %exon_starts;
	my %exon_ends;
	foreach my $e ( @{$trans->get_all_Exons} ) {
	  $exon_starts{$e->cdna_start($trans)} = 1;
	  $exon_ends{$e->cdna_end($trans)} = 1;
	}
	$self->exon_starts($trans->display_id,\%exon_starts);
	$self->exon_ends($trans->display_id,\%exon_ends);
      }
      my $ee = undef;
      my $es = undef;
      # scan along the read and look for overlaps with the exon boundaries
      for ( my $i = $f->start ; $i <= $f->end ; $i++ ){
	# exon end position within the feature;
	$ee = $i - $f->start if $self->exon_ends($trans->display_id)->{$i};
	# next exon start position within the feature
	$es = $i - $f->start if $self->exon_starts($trans->display_id)->{$i};
      }
      # allow reads covering introns
      $accept = 1 if ( $es && $ee && $ee >= $self->INTRON_OVERLAP 
		       && $es <= $f->length - $self->INTRON_OVERLAP );

      # disallow unspliced reads spanning introns  unless the intron is size 0
      # then allow the unspliced alignment to  bridge the gap
      # also allow any spliced modles even if they dont lie in an intron
      
      if ( $self->INTRON_MODELS ) {
	unless ( $f->{"_intron"} or ( $es && $ee && ($es - $ee == 1) ) ) {
	  $accept = 0;
	}
      }
    }

    # restrict just to splice models where read is within an exon
    if ( $self->INTRON_MODELS ) {
      $accept = 1 if $f->{"_intron"};
    }

    next FEATURE unless $accept;

    if ( $self->PROJECT ) {
      my @mapper_objs;
      my @features;
      my $start = 1;
      my $cannonical = 1;
      my $end = $f->length;
      my $out_slice = $slice_adaptor->fetch_by_name($trans->slice->name);
      # get as ungapped features
      foreach my $ugf ( $f->ungapped_features ) {
	# Project onto the genome
	foreach my $obj ($trans->cdna2genomic($ugf->start, $ugf->end)){
	  if( $obj->isa("Bio::EnsEMBL::Mapper::Coordinate")){
	    # make into feature pairs?
	    my $qstrand = $f->strand  * $trans->strand;
	    my $hstrand = $f->hstrand * $trans->strand;
	    my $fp;
	    $fp = Bio::EnsEMBL::FeaturePair->new
	      (-start    => $obj->start,
	       -end      => $obj->end,
	       -strand   => $qstrand,
	       -slice    => $trans->slice,
	       -hstart   => 1,
	       -hend     => $obj->length,
	       -hstrand  => $hstrand,
	       -percent_id => $f->percent_id,
	       -score    => $f->score,
	       -hseqname => $f->hseqname,
	       -hcoverage => $f->hcoverage,
	       -p_value   => $f->p_value,
	      );
	    push @features, $fp->transfer($out_slice);
	  }
	}
      }
      @features = sort { $a->start <=> $b->start } @features;
      # if we have a spliced alignment check to see if it's non-cannonical
      # if so tag it so we can tell later on

      if ( $f->{'_intron'} ) {
	print $f->hseqname . "  ";
	if ( scalar(@features) == 2 ) {
	  my $left_splice = $slice_adaptor->fetch_by_region('toplevel',
							    $out_slice->seq_region_name,
							    $features[0]->end+1,
							    $features[0]->end+2,
							    $features[0]->strand
							   );
	  my $right_splice = $slice_adaptor->fetch_by_region('toplevel',
							     $out_slice->seq_region_name,
							     $features[1]->start-2,
							     $features[1]->start-1,
							     $features[0]->strand
							    );	
	  if ( $left_splice->seq eq 'NN' && $right_splice->seq eq 'NN' ) {
	    warn("Cannot find dna sequence for " . $f->display_id .
		 " this is used in detetcting non cannonical splices\n");
	  } else {
	    #�is it cannonical
	    if ( $features[0]->strand  == 1 ) {
	      print "Splice type " . $left_splice->seq ."- ".  $right_splice->seq ." ";
	      #�is it GTAG?
	      unless ( $left_splice->seq eq 'GT' && $right_splice->seq eq 'AG' ) {
		$cannonical = 0;
	      }
	    } else {
	      print "Splice type " . $right_splice->seq ."- ".  $left_splice->seq ." ";
	      #�is it GTAG?
	      unless ( $right_splice->seq eq 'GT' && $left_splice->seq eq 'AG' ) {
		$cannonical = 0;
	      }
	    }
	  }
	}
      }
      print "Making feat " . $f->hseqname . "\n";
      my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@features);
      # corect for hstart end bug
      $feat->hstart($f->hstart);
      $feat->hend($f->hend);
      $feat->analysis($self->analysis);
      unless ( $cannonical ) {
	print "Non cannonical \n";
	#�mark it as non cannonical splice
	$feat->hseqname($feat->hseqname.":NC");
      } else {
	print "Cannonical \n";
      }
      # dont store the same feature twice because it aligns to a different transcript in the same gene.
      my $unique_id = $feat->seq_region_name .":" .
	$feat->start .":" .
	  $feat->end .":" .	
	    $feat->strand .":" .
	      $feat->hseqname;
      unless ($self->stored_features($unique_id)){
	push @dafs,$feat;
	# keep tabs on it so you don't store it again.
	$self->stored_features($unique_id,1);
      }
    } else {
      # just write the features on the transcript coord system 
      # assming there is one ( useful for writing reads that support introns 
      # on the transcrtipt coord system
      push @dafs,$f;
    }
  }
  return \@dafs;
}


##########################################################################

sub transcripts {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_transcripts'} = $value;
  }

  if (exists($self->{'_transcripts'})) {
    return $self->{'_transcripts'};
  } else {
    return undef;
  }
}


sub exon_starts {
  my ($self,$key,$value) = @_;

  return undef unless defined ($key);

  if (defined $value) {
    $self->{'_exon_starts'}{$key} = $value;
  }

  if (exists($self->{'_exon_starts'}{$key})) {
    return $self->{'_exon_starts'}{$key};
  } else {
    return undef;
  }
}

sub exon_ends {
  my ($self,$key,$value) = @_;

  return undef unless defined ($key);

  if (defined $value) {
    $self->{'_exon_ends'}{$key} = $value;
  }

  if (exists($self->{'_exon_ends'}{$key})) {
    return $self->{'_exon_ends'}{$key};
  } else {
    return undef;
  }
}

sub stored_features {
  my ($self,$key,$value) = @_;

  return undef unless defined ($key);

  if (defined $value) {
    $self->{'_stored_features'}{$key} = $value;
  }

  if (exists($self->{'_stored_features'}{$key})) {
    return $self->{'_stored_features'}{$key};
  } else {
    return undef;
  }
}
1;
