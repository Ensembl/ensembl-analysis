
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
  my @genes;
  my %trans_by_id;
  my $biotype = $self->BIOTYPE;
  my @slices = @{$trans_db->get_SliceAdaptor->fetch_all('toplevel')};
  foreach my $slice (@slices){
    # fetch genes, transcripts and exons
    if ( $biotype ){
      push @genes , @{$slice->get_all_Genes_by_type($biotype,undef,1)};
    } else {
      push @genes , @{$slice->get_all_Genes(undef,undef,1)};
    }

  }

  # store them in a hash so they can be retreived quickly

  if (scalar(@genes) > 0 ){
    foreach my $gene ( @genes ) {
      foreach my $trans ( @{$gene->get_all_Transcripts} ) {
	$trans_by_id{$trans->display_id} = $trans;
      }
    }
    $self->transcripts(\%trans_by_id);
  }
  return ;
}


sub run {
  my ($self) = @_;

  throw("Can't run - no runnable objects") unless ( $self->runnable );

  my ($runnable) = @{$self->runnable};

  $runnable->run;
  my $features = $runnable->output;

  if ($self->filter) {
    $features = $self->filter->filter_results($features);
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

  my $slice_adaptor = $self->db->get_SliceAdaptor;
  my @dafs;

FEATURE:  foreach my $f (@$flist) {

    my $trans_id = $f->seqname;

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
      for ( my $ i = $f->start ; $i <= $f->end ; $i++ ){
	# exon end position within the feature;
	$ee = $i - $f->start if $self->exon_ends($trans->display_id)->{$i};
	# next exon start position within the feature
	$es = $i - $f->start if $self->exon_starts($trans->display_id)->{$i};
      }
      next FEATURE unless ( $es && $ee && $ee >= $self->INTRON_OVERLAP 
			    && $es <= $f->length - $self->INTRON_OVERLAP );
    }

    my @mapper_objs;
    my @features;
    my $start = 1;
    my $end = $f->length;
    my $out_slice = $slice_adaptor->fetch_by_name($trans->slice->name);
    # get as ungapped features
    foreach my $ugf ( $f->ungapped_features ) {
      # Project onto the genome
      push @mapper_objs,$trans->cdna2genomic($ugf->start, $ugf->end);
    }

    # Convert all these mapper objects into a dna_align_feature on the new coord system
    foreach my $obj ( sort { $a->start <=> $b->start }@mapper_objs ){
      if( $obj->isa("Bio::EnsEMBL::Mapper::Coordinate")){
	# make into feature pairs?
	my $fp;
	if ( $f->hstrand ==  1 ){
	  $fp = Bio::EnsEMBL::FeaturePair->new
	    (-start    => $obj->start,
	     -end      => $obj->end,
	     -strand   => 1,
	     -slice    => $trans->slice,
	     -hstart   => $start,
	     -hend     => $start+$obj->length-1,
	     -hstrand  => $f->hstrand,
	     -percent_id => $f->percent_id,
	     -score    => $f->score,
	     -hseqname => $f->hseqname,
	     -hcoverage => $f->hcoverage,
	     -p_value   => $f->p_value,
	    );
	  $start += $obj->length;
	} else {
	  $fp = Bio::EnsEMBL::FeaturePair->new
	    (-start    => $obj->start,
	     -end      => $obj->end,
	     -strand   => 1,
	     -slice    => $trans->slice,
	     -hstart   => $end-$obj->length+1,
	     -hend     => $end,
	     -hstrand  => $f->hstrand,
	     -percent_id => $f->percent_id,
	     -score    => $f->score,
	     -hseqname => $f->hseqname,
	     -hcoverage => $f->hcoverage,
	     -p_value   => $f->p_value,
	    );
	  $end -= $obj->length;
	}
	push @features, $fp->transfer($out_slice);
      }
    }
    my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@features);
    $feat->analysis($self->analysis);
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
