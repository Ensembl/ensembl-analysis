package Bio::EnsEMBL::Analysis::Tools::Pmatch::First_PMF;
use Bio::EnsEMBL::Pipeline::Tools::Pmatch::ContigHit;
use Bio::EnsEMBL::Pipeline::Tools::Pmatch::ProteinHit;
use Bio::EnsEMBL::Pipeline::Tools::Pmatch::CoordPair;
use Bio::EnsEMBL::Pipeline::Tools::Pmatch::MergedHit;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw();

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($plengths, $maxintronlen, $min_coverage) = 
                 rearrange(['PROTEIN_LENGTHS','MAX_INTRON_LENGTH', 'MIN_COVERAGE'], @args);

  #Setting a default
  $self->min_coverage(25);

  #print STDERR "Creating FIRST_PMF with @args\n";
  throw("No protlengths data") unless defined $plengths;
  $self->protein_lengths($plengths);

  throw("No maximum intron length") unless defined $maxintronlen;
  $self->max_intron_len($maxintronlen);
  $self->min_coverage($min_coverage);
  return $self;

}

=head2 make_coord_pair

 Title   : make_coord_pair
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_coord_pair {
  my ($self) = @_;
  #print STDERR "Making coord pair with ".$_."\n";
  my @cols = split();
  # clean up the refseq IDs
  my $protid = $cols[5];
  # alter this as necessary 
  if($protid =~ /\w+\|\w+\|\w+\|([\w\.]+)\|/) {
    $protid = $1;
  }
  
  # sort out strand hmmm if strand = -1 then shouldn't we switch qstart & qend? Currently don't ...
  my $strand = 1;
  if ( $cols[3] < $cols[2]  ) { $strand=-1; }
  
  my $cp = new Bio::EnsEMBL::Pipeline::Tools::Pmatch::CoordPair(
							 -query  => $cols[1],
							 -target => $protid,
							 -qstart => $cols[2],
							 -qend   => $cols[3],
							 -tstart => $cols[6],
							 -tend   => $cols[7],
							 -percent=> $cols[8],
							 -strand => $strand,
			);

  # where to put the CoordPair?
  # find the relevant ProteinHit, or make a new one
  my %proteins;
  my $ph = $proteins{$protid};

  if (!defined $ph) {
     #make a new one and add it into %proteins
    $ph = new Bio::EnsEMBL::Pipeline::Tools::Pmatch::ProteinHit(-id=>$protid);
    $proteins{$protid} = $ph;
  }
  $self->proteins(\%proteins);
  # now find the relevant ContigHit, or make a new one
  my $ch = $ph->get_ContigHit($cols[1]);
  if (!defined $ch) {
     # make a new one and add it into $ph
    $ch = new Bio::EnsEMBL::Pipeline::Tools::Pmatch::ContigHit(-id => $cols[1]);
    $ph->add_ContigHit($ch);
  }
  
  # now add the CoordPair
  $ch->add_CoordPair($cp);

}

=head2 merge_hits

 Title   : merge_hits
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub merge_hits {
  my ($self) = @_;
  my @merged;

  # merge the hits together.
  my $protref = $self->proteins;

  foreach my $p_hit(values %$protref) {
    my $allhits = $self->make_mergelist($p_hit);
    
    my $chosen = $self->prune($allhits);
    #print STDERR "\nNo hits good enough for " . $p_hit->id() . "\n"
    #  unless scalar(@chosen);
    push(@merged,@$chosen);
  }
  #print STDERR "PMF have ".@merged." ".$merged[0]." results\n";
  return \@merged;
  
}

=head2 make_mergelist

 Title   : make_mergelist
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_mergelist {
  my ($self, $ph) = @_;

  my @hits = ();
  foreach my $ch($ph->each_ContigHit()) {
    # forward strand
    my @cps = $ch->each_ForwardPair();
    # sort!
    @cps = sort {$a->qstart <=> $b->qstart} @cps;
    #reverse strand
    my @mps = $ch->each_ReversePair();
    @mps = sort {$b->qstart <=> $a->qstart} @mps;
    push (@cps,@mps);
    
    # deal with forward & reverse separately?
    my $first = shift(@cps);
    my $mh = $self->new_merged_hit($first);
    
    push(@hits,$mh);

    my $prev = $hits[$#hits];
    my $prev_cp = $first;

    CP: 
    foreach my $cp(@cps) {

      # need to compare with the last entry in @merged

      # first the strand
      my $strand = 1;
      if ($cp->qend() < $cp->qstart()) { $strand = -1; }

      # does this CoordPair extend the current hit?
      # need a fudge factor - pmatch could overlap them by 1 or 2 ... or 3
      if( $strand == $prev->strand &&
	 ( ($cp->tstart >=  $prev_cp->tend) ||
	 (($prev_cp->tend -  $cp->tstart) <= 3)) &&
	  # Steve's fix - Added qstart/qend condition (*strand should allow it to work on
	  #     either strand)
	  # no overlap currently allowed
          $cp->qstart*$strand >= $prev_cp->qend*$strand 
	)
	{
	  #extend existing MergedHit
	  my $coverage = $cp->tend - $cp->tstart + 1;
	  $coverage += $prev->coverage();

	  # compensate for any overlap 
	  my $diff = $prev_cp->tend -  $cp->tstart;
	  if($diff >=0) {
	    $diff++;
	    $coverage -= $diff;
	  }

	  $prev->coverage($coverage);
	  $prev->add_CoordPair($cp);
          $prev_cp = $cp;
	}
      else {
	# make a new MergedHit
	my $mh = $self->new_merged_hit($cp);
	push(@hits,$mh);	

        $prev = $hits[$#hits];
        $prev_cp = $cp;
      }
    }
  }

  # my $times = join " ", times;
  # print "Before extend " . $times . "\n";
  # extend those merged hits

  # print "Number of hits = " . scalar(@hits) . "\n";

  @hits = @{$self->extend_hits(\@hits)};

  # $times = join " ", times;
  # print "After extend " . $times . "\n";
  # print "Number of hits = " . scalar(@hits) . "\n";

  # sort out coverage
  my $plengths = $self->protein_lengths;
  my $protlen = $$plengths{$ph->id};
  warn "No length for " . $ph->id . "\n" unless $protlen;
  return unless $protlen; 
  foreach my $hit(@hits) {
    my $percent = $hit->coverage;
    $percent *= 100;
    $percent /= $protlen;
    $percent=sprintf("%.1f",$percent);
    $hit->coverage($percent);
  }

  return \@hits;
}

=head2 new_merged_hit

 Title   : new_merged_hit
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_merged_hit {
  my ($self,$cp) = @_;
  my $coverage = $cp->tend - $cp->tstart + 1;
  my $mh = new Bio::EnsEMBL::Pipeline::Tools::Pmatch::MergedHit( -query    =>  $cp->query(),
			  -target   =>  $cp->target(),
			  -strand   =>  $cp->strand(),
			  -coverage =>  $coverage,
			);
  $mh->add_CoordPair($cp);
  return $mh;
}

sub proteins{
  my ($self, $arg) = @_;
  if(!$self->{'proteins'}){
    $self->{'proteins'} = {};
  }
  if($arg){
    $self->{'proteins'} = $arg;
  }
  return $self->{'proteins'};
}

sub protein_lengths {
  my ($self, $plengths) = @_;
  
  # $plengths is a hash reference
  if(defined($plengths)){
    if (ref($plengths) eq "HASH") {
      $self->{_plengths} = $plengths;
    } 
    else {
      throw("[$plengths] is not a hash ref.");
    }
  }
     
  return $self->{_plengths};

}



sub max_intron_len {
  my ($self, $maxintronlen) = @_;
  
  if(defined($maxintronlen)){
    if ($maxintronlen =~ /\d+/) {
      $self->{_maxintronlen} = $maxintronlen;
    } 
    else {
      throw("[$maxintronlen] is not numeric.");
    }
  }
  return $self->{_maxintronlen};
  
}

sub min_coverage {
  my ($self, $mincoverage) = @_;
  
  if(defined($mincoverage)){
    if ($mincoverage =~ /\d+/) {
      $self->{_mincoverage} = $mincoverage;
    } 
    else {
      throw("[$mincoverage] is not numeric.");
    }
  }
  return $self->{_mincoverage};
  
}

sub output{
  my ($self) = @_;
  if (!defined($self->{_output})) {
    $self->{_output} = [];
  } 
  return $self->{_output};
}

sub extend_hits {
  my ($self, $hits) = @_;

  # we want to do essentially what we did to create the merged hits but we can skip over
  # intervening hits if we need to.
  my @fhits;
  my @rhits;
  my @newhits;
  #separate the forward and reverse strand hits into 2 separate arrays
  foreach my $mh(@$hits){
    if($mh->strand == 1){
      push (@fhits, $mh);
    }
    else {
      push (@rhits, $mh);
    }
  }

  @fhits = sort {$a->qstart <=> $b->qstart} @fhits; #sort the forward strand 
  #lowest to highest
  @rhits = sort {$b->qstart <=> $a->qstart} @rhits;
  #and the reverse strand highest to lowest
  #This is so both sets of features sit in 5" to 3" order

  while(scalar(@fhits)){ #This is a recursive process
    my $hit = shift(@fhits); #We take current most 5" hit on the array
    push (@newhits, $hit); #and add it to the new array

    # can we link it up to a subsequent hit?
    foreach my $sh(@fhits){ #For each of the subsequenct hits
      die ("Argh!") unless $hit->strand == $sh->strand; #die if its not on the
      #same strand
      last if ($hit->qend+$self->max_intron_len < $sh->qstart);
      if ($sh->tstart > $hit->tend &&
	  $sh->tend   > $hit->tend &&
	  abs($sh->tstart - $hit->tend) <= 3 &&
	  # qstart/qend condition - no overlap currently allowed
	  $sh->qstart >= $hit->qend
	 ) {
	# add the coord pairs from $sh into $hit
        $hit->subsume_MergedHit($sh);
      }
    }
  }

  
# same for rhits
  while(scalar(@rhits)){
    my $hit = shift(@rhits);
    push (@newhits, $hit);
    

    # can we link it up to a subsequent hit?
    foreach my $sh(@rhits){
      die ("Argh!") unless $hit->strand == $sh->strand;
# hmmmm On minus strand qstart is currently > qend
# ie $sh->qend <= $sh->qstart <= $hit->qend <= $hit->qstart
#      last if ($hit->qstart-$self->max_intron_len > $sh->qend);
      last if ($hit->qend-$self->max_intron_len > $sh->qstart);
      if ($sh->tstart > $hit->tend &&
	  $sh->tend   > $hit->tend &&
	  abs($sh->tstart - $hit->tend) <= 3 &&
	  # qstart/qend condition - no overlap currently allowed
	  $hit->qend >= $sh->qstart
	 ) {
	# add the coord pairs from $sh into $hit
        $hit->subsume_MergedHit($sh);
      }
    }
  }

  # return extended hits
  return (\@newhits);  
}

=head2 prune

 Title   : prune
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub prune {
  my ($self, $all) = @_;
  my @chosen = ();
  # reject hits with < 25% coverage
  my $lower_threshold = $self->min_coverage;

  # sort by descending order of coverage
  my @all = sort {$b->coverage <=> $a->coverage} @$all;

  my $first = shift(@all);
  if($first->coverage < $lower_threshold){
    # we're done
    return \@chosen;
  }
  
  push (@chosen,$first);

  # don't select any hits that have coverage less than 2% below that of the first hit, be it 100 or 99.9 or ...
  my $curr_pc = $first->coverage() - 2;

 PRUNE:
  foreach my $hit(@all) {
    last PRUNE if $hit->coverage < $lower_threshold;
    last PRUNE if $hit->coverage < $curr_pc;
    push (@chosen,$hit);
  }

  return \@chosen;

}

1;
