
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

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Finished::MiniEst2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::Runnable::Finished::MiniEst2Genome->new('-genomic'    => $genseq,
								    '-features'   => $features,
								    '-seqcache' => $seqcache,
								    '-analysis'   => $analysis,
								   )

    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Modified by Sindhu K. Pillai B<email> sp1@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

Modified by Sindhu K. Pillai B<email> sp1@sanger.ac.uk
=cut

package Bio::EnsEMBL::Analysis::Runnable::Finished::MiniEst2Genome;

use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::Tools::MiniSeq;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::DB::RandomAccessI;
use Bio::PrimarySeqI;
use Bio::SeqIO;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Est2Genome;
use Bio::EnsEMBL::DnaDnaAlignFeature;

use base 'Bio::EnsEMBL::Analysis::Runnable';

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_fplist'} = []; #create key to an array of feature pairs

  my( $genomic, $features, $seqcache, $analysis ) = rearrange([qw(GENOMIC
								 FEATURES
								 SEQCACHE
								 ANALYSIS)], @args);

  $self->throw("No genomic sequence input")
    unless defined($genomic);
  $self->throw("[$genomic] is not a Bio::PrimarySeqI")
    unless $genomic->isa("Bio::PrimarySeqI");
  $self->genomic_sequence($genomic) if defined($genomic);

  $self->throw("No seqcache provided")
    unless defined($seqcache);
  $self->seq_cache($seqcache) if defined($seqcache);
  $self->analysis($analysis) if defined $analysis;

  if (defined($features)) {
    if (ref($features) eq "ARRAY") {
      my @f = @$features;
      foreach my $f (@f) {
	$self->addFeature($f);
      }
    } else {
      $self->throw("[$features] is not an array ref.");
    }
  }
  return $self; # success - we hope!
}

=head2 addFeature

    Title   :   addFeature
    Usage   :   $self->addFeature($f)
    Function:   Adds a feature to the object for realigning
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   Bio::EnsEMBL::FeaturePair

=cut

sub addFeature {

    my( $self, $value ) = @_;

    if(!defined($self->{'_features'})) {
	$self->{'_features'} = [];
    }

    if ($value) {
        $value->isa("Bio::EnsEMBL::FeaturePair") || $self->throw("Input isn't a Bio::EnsEMBL::FeaturePair");
	push(@{$self->{'_features'}},$value);
    }
}


=head2 get_all_FeaturesById

    Title   :   get_all_FeaturesById
    Usage   :   $hash = $self->get_all_FeaturesById;
    Function:   Returns a ref to a hash of features.
                The keys to the hash are distinct feature ids
    Returns :   ref to hash of Bio::EnsEMBL::FeaturePair
    Args    :   none

=cut

sub get_all_FeaturesById {
    my( $self) = @_;
    my  %idhash;

    FEAT: foreach my $f ($self->get_all_Features) {
    if (!(defined($f->hseqname))) {
	warning("No hit name for " . $f->seqname . "\n");
	    next FEAT;
	}
	if (defined($idhash{$f->hseqname})) {
	    push(@{$idhash{$f->hseqname}},$f);
	} else {
	    $idhash{$f->hseqname} = [];
	    push(@{$idhash{$f->hseqname}},$f);
	}

    }

    return (\%idhash);
}


=head2 get_all_Features

    Title   :   get_all_Features
    Usage   :   @f = $self->get_all_Features;
    Function:   Returns the array of features
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   none

=cut


sub get_all_Features {
    my( $self, $value ) = @_;
    return (@{$self->{'_features'}});
}


=head2 get_all_FeatureIds

  Title   : get_all_FeatureIds
  Usage   : my @ids = get_all_FeatureIds
  Function: Returns an array of all distinct feature hids
  Returns : @string
  Args    : none

=cut

sub get_all_FeatureIds {
    my ($self) = @_;

    my %idhash;

    foreach my $f ($self->get_all_Features) {
	if (defined($f->hseqname)) {
	    $idhash{$f->hseqname} = 1;
	} else {
	    warning("No sequence name defined for feature. " . $f->seqname . "\n");
	}
    }

    return keys %idhash;
}

=head2 make_miniseq

  Title   : make_miniseq
  Usage   :
  Function: makes a mini genomic from the genomic sequence and features list
  Returns :
  Args    :

=cut

sub make_miniseq {

    my ($self,@features) = @_;
    my $seqname = $features[0]->seqname;
    @features = sort {$a->start <=> $b->start} @features;
    my $count  = 0;
    my $mingap = $self->minimum_intron;

    my $pairaln = new Bio::EnsEMBL::Analysis::Tools::PairAlign;

    my @genomic_features;

    my $prevend     = 0;
    my $prevcdnaend = 0;
  FEAT: foreach my $f (@features) {

      my $start = $f->start;
      my $end   = $f->end;
      $start = $f->start - $self->exon_padding;
      $end   = $f->end   + $self->exon_padding;

      if ($start < 1) { $start = 1;}
      if ($end   > $self->genomic_sequence->length) {$end = $self->genomic_sequence->length;}

      my $gap     =    ($start - $prevend);

      if ($count > 0 && ($gap < $mingap)) {
	# STRANDS!!!!!
	  if ($end < $prevend) { $end = $prevend;}
	  $genomic_features[$#genomic_features]->end($end);
	  $prevend     = $end;
	  $prevcdnaend = $f->hend;

      } else {

	    my $newfeature = new Bio::EnsEMBL::SeqFeature;

	    $newfeature->seqname ($f->hseqname);
	    $newfeature->start     ($start);
	    $newfeature->end       ($end);
	    $newfeature->strand    (1);
# ???	    $newfeature->strand    ($strand);
	    $newfeature->attach_seq($self->genomic_sequence);

	    push(@genomic_features,$newfeature);
	    $prevend = $end;
	    $prevcdnaend = $f->hend;

	}
	$count++;
    }

    # Now we make the cDNA features
    # but presumably only if we actually HAVE any ...
    return unless scalar(@genomic_features);

    my $current_coord = 1;

    # make a forward strand sequence, est2genome runs -reverse
    @genomic_features = sort {$a->start <=> $b->start } @genomic_features;

    foreach my $f (@genomic_features) {
	$f->strand(1);
	my $cdna_start = $current_coord;
	my $cdna_end   = $current_coord + ($f->end - $f->start);

	my $tmp = new Bio::EnsEMBL::SeqFeature(
					       -seqname => $f->seqname.'.cDNA',
					       -start => $cdna_start,
					       -end   => $cdna_end,
					       -strand => 1);

	my $fp  = new Bio::EnsEMBL::FeaturePair(-feature1 => $f,
						-feature2 => $tmp);

	$pairaln->addFeaturePair($fp);

#	$self->print_FeaturePair($fp);

	$current_coord = $cdna_end+1;
    }

    #changed id from 'Genomic' to seqname
    my $miniseq = new Bio::EnsEMBL::Analysis::Tools::MiniSeq(-id        => $seqname,
						      -pairalign => $pairaln);

    my $newgenomic = $miniseq->get_cDNA_sequence->seq;
    $newgenomic =~ s/(.{72})/$1\n/g;
#    print ("New genomic sequence is " . $newgenomic. "\n");
    return $miniseq;

}


=head2 minimum_intron

  Title   : minimum_intron
  Usage   :
  Function: Defines minimum intron size for miniseq
  Returns :
  Args    :

=cut

sub minimum_intron {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_minimum_intron'} = $arg;
    }

    return $self->{'_minimum_intron'} || 1000;
}

=head2 exon_padding

  Title   : exon_padding
  Usage   :
  Function: Defines exon padding extent for miniseq
  Returns :
  Args    :

=cut

sub exon_padding {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_padding'} = $arg;
    }

#    return $self->{'_padding'} || 100;
    return $self->{'_padding'} || 1000;

}

=head2 print_FeaturePair

  Title   : print_FeaturePair
  Usage   :
  Function: for debugging
  Returns :
  Args    :

=cut

sub print_FeaturePair {
    my ($self,$nf) = @_;
    #changed $nf->id to $nf->seqname
    print(STDERR "FeaturePair is " . $nf->seqname    . "\t" .
	  $nf->start . "\t" .
	  $nf->end   . "\t(" .
	  $nf->strand . ")\t" .
	  $nf->hseqname  . "\t" .
	  $nf->hstart   . "\t" .
	  $nf->hend     . "\t(" .
	  $nf->hstrand  . ")\n");
}




=head2 run_blaste2g

  Title   : run_blaste2g
  Usage   : $self->run_blaste2g()
  Function: Runs est2genome on a MiniSeq
  Returns : none
  Args    :

=cut

sub run_blaste2g {
    my ( $self, $est, $features,$analysis ) = @_;
    my $count = @$features;
    my $miniseq = $self->make_miniseq(@$features);
    my $hseq    = $self->seq_cache->{$est} or throw("Can't find sequence id '$est' in seq_cache");
	my $g = $miniseq->get_cDNA_sequence;

    my $eg = new Bio::EnsEMBL::Analysis::Runnable::Finished::Est2Genome(
        -genomic => $g,
        -est     => $hseq,
		-analysis => $analysis,
    );

    $eg->run;

    foreach my $fp ( @{$eg->output} ) {

      my @converted = @{$miniseq->convert_FeaturePair($fp)};

      if ( @converted > 1 ) {
	warning "feature converts into '" . scalar(@converted) . "' > 1 features - ignoring\n";
      } else {
	# convert_FeaturePair zaps strand and hseqname,0
	# so we put them back here.
	my $new = $converted[0];

#            $new->seqname( $fp->hseqname );  #
#            $new->strand( $fp->strand );     #
#            $new->hseqname( $fp->hseqname ); #
#            $new->percent_id( $fp->percent_id ); #
#            $new->feature2->percent_id( $fp->percent_id ); #
# reverse logic here, change the $fp passed to equal the converted [miniseq]
# change seqname, hend, hstart, gsf_start, gsf_end
# hend & hstart appear to only be calculated from the length of the end & start
# producing the same length for the hit as the contig.  not good....
#	    $fp->hend(   $new->hend   );
#	    $fp->hstart( $new->hstart );
	$fp->end(    $new->end   );
	$fp->start(  $new->start );
#            $self->add_output(@converted);
#	    print STDERR "CIGAR_STRING : " . $fp->cigar_string . "\n";
#	    print STDERR $fp->gffstring . "\n";
	$self->add_output($fp);

      }
    }
  }


=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}

=head2 seq_cache

    Title   :   seq_cache
    Usage   :   $self->seq_cache($hash)
    Function:   Get/set method for Seq Cache
    Returns :   HASH
    Args    :   HASH

=cut

sub seq_cache {
    my( $self, $hash ) = @_;
    if ($hash) {
        $self->{'_seq_cache'} = $hash;
    }
    return $self->{'_seq_cache'};
}

=head2 analysis

    Title   :   analysis
    Usage   :   $self->analysis($analysis)
    Function:   Get/set method for analysis
    Returns :   Bio::EnsEMBL::Analysis object
    Args    :   Bio::EnsEMBL::Analysis object

=cut

sub analysis {
    my( $self, $value ) = @_;
    if ($value) {
        $value->isa("Bio::EnsEMBL::Analysis") || throw("[$value] isn't a Bio::EnsEMBL::Analysis");
        $self->{'analysis'} = $value;
    }
    return $self->{'_analysis'};
}


sub find_extras {
  my ($self,@features) = @_;
  my @output = $self->output;
  my @new;
 FEAT: foreach my $f (@features) {
    my $found = 0;
    if (($f->end - $f->start) < 50) {
      next FEAT;
    }
    #	print ("New feature\n");
    #$self->print_FeaturePair($f);
    foreach my $out (@output) {
      foreach my $sf ($out->sub_SeqFeature) {

	if (!($f->end < $out->start || $f->start >$out->end)) {
	  $found = 1;
	}
      }
    }
    if ($found == 0) {
      push(@new,$f);
    }
  }
  return \@new;
}

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
    my ($self) = @_;
    if (!defined($self->{'_output'})) {
	$self->{'_output'} = [];
    }
    return \@{$self->{'_output'}};
}




=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs est2genome on MiniSeq representation of genomic sequence for each EST
  Returns : none
  Args    :

=cut


sub run {

    my ($self) = @_;
    my ($esthash) = $self->get_all_FeaturesById;
    my @ests    = keys %$esthash;
    my $number_of_errors = 0;

    foreach my $est (@ests) {
        my $features = $esthash->{$est};
        my @exons;
        next unless (ref($features) eq "ARRAY");
        next unless (scalar(@$features) >= 1);
        $self->run_blaste2g($est, $features,$self->{analysis});
    }

    return 1;
}

sub add_output{
    my($self, @feat_pairs) = @_;
    push( @{ $self->{'_output'} } , @feat_pairs);

}

1;

__END__

