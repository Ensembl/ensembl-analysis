###############################################################################
# Bio::EnsEMBL::Analysis::Tools::BPlite::Sbjct
###############################################################################
#
# The original BPlite.pm module has been written by Ian Korf !
# see http://sapiens.wustl.edu/~ikorf
#
# You may distribute this module under the same terms as perl itself

package Bio::EnsEMBL::Analysis::Tools::BPlite::Sbjct;

use warnings ;
use strict;

use Scalar::Util qw(looks_like_number);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::BPlite::HSP; # we want to use HSP
#use overload '""' => 'name';
use vars qw(@ISA);

@ISA = qw();

sub new {
    my ($caller, @args) = @_;
    
    my $class = ref($caller) || $caller;
    my $self = bless({}, $class);

    
    ($self->{'NAME'},$self->{'LENGTH'},$self->{'FH'},
     $self->{'LASTLINE'},$self->{'PARENT'}) =
	 rearrange([qw(NAME
			       LENGTH
			       FH
			       LASTLINE
			       PARENT
			       )],@args);
    
    $self->{'HSP_ALL_PARSED'} = 0;
    
  return $self;
}

=head2 name

 Title    : name
 Usage    : $name = $obj->name();
 Function : gets/sets the name of the Sbjct
 Example  :
 Returns  : name of the Sbjct
 Args     :

=cut

sub name {
  my ($self,$name) = @_;
  if($name) {
    $self->{'NAME'} = $name;
  }
  return($self->{'NAME'});
}

=head2 query_name

 Title    : query_name
 Usage    : $query_name = $obj->query_name();
 Function : returns the name of the query sequence
 Example  :
 Returns  : name of the query
 Args     :

=cut

sub query_name {
 my ($self) = @_;
 return($self->{'PARENT'}->{'QUERY'});
}


=head2 slice_name

 Title    : slice_name
 Usage    : $slice_name = $obj->slice_name();
 Function : getter/setter for the name of the slice
 Example  :
 Returns  : name of the slice
 Args     :

=cut

sub slice_name {
  my ($self,$slice_name) = @_;
  if($slice_name) {
    $self->{'SLICE_NAME'} = $slice_name;
  }
  return($self->{'SLICE_NAME'});
}



=head2 nextFeaturePair

 Title    : nextFeaturePair
 Usage    : $name = $obj->nextFeaturePair();
 Function : same as the nextHSP function 
 Example  : 
 Returns  : next FeaturePair 
 Args     :

=cut

sub nextFeaturePair {shift->nextHSP}; # just another name

=head2 nextHSP

 Title    : nextHSP
 Usage    : $hsp = $obj->nextHSP();
 Function : returns the next available High Scoring Pair
 Example  : 
 Returns  : Bio::Tools::HSP  or null if finished
 Args     :

=cut

sub nextHSP {
  my ($self) = @_;
  return undef if $self->{'HSP_ALL_PARSED'};
  
  ############################
  # get and parse scorelines #
  ############################

  my $scoreline = $self->{'LASTLINE'};
  my $FH = $self->{'FH'};
  my $nextline = <$FH>;
  return undef if not defined $nextline;
  $scoreline .= $nextline;
  my ($score, $bits);
  if ($scoreline =~ /\d bits\)/) {
    ($score, $bits) = $scoreline =~
      /Score = (\d+) \((\S+) bits\)/; # WU-BLAST
  }
  else {
    ($bits, $score) = $scoreline =~
      /Score =\s+(\S+) bits \((\d+)/; # NCBI-BLAST
  }
  
  my ($match, $length, $percent) = $scoreline =~ /Identities = (\d+)\/(\d+)(?:\s*\((\d+)\%\))?/;
  my ($positive) = $scoreline =~ /Positives = (\d+)/;
  my $frame = '0';
  $positive = $match if not defined $positive;
  my ($p)        = $scoreline =~ /[Sum ]*P[\(\d+\)]* = ([^\,\s]+)/;
  if (not defined $p) {(undef, $p) = $scoreline =~ /Expect(\(\d+\))? =\s+([^\,\s]+)/}
  $p = 0 if $p eq '0.';         # (DBD::)MySQL doesn't like '0.' for zero. 0.0 or 0 is fine though.
  $p = "1$p" if $p =~ /^(e-?\d+)$/; # e-169 => 1e-169
  if (defined($p) and not looks_like_number($p)) {
      my $name = $self->name;
      throw("Hit '$name': P value [$p] not numeric in '$scoreline'");
  }
  throw("Unable to parse '$scoreline'") if not defined $score;
  
  #######################
  # get alignment lines #
  #######################
  my @hspline;
  while(<$FH>) {
    if ($_ =~ /^WARNING:|^NOTE:/) {
      while(<$FH>) {last if $_ !~ /\S/}
    }
    elsif ($_ !~ /\S/)            {next}
    elsif ($_ =~ /Strand HSP/)    {next} # WU-BLAST non-data
    elsif ($_ =~ /^\s*Strand/)    {next} # NCBI-BLAST non-data
    elsif ($_ =~ /^\s*Score/)     {$self->{'LASTLINE'} = $_; last}
    elsif ($_ =~ /^>|^Parameters|^\s+Database:|^CPU\stime|^Lambda/)   {
      $self->{'LASTLINE'} = $_;
      $self->{'PARENT'}->{'LASTLINE'} = $_;
      $self->{'HSP_ALL_PARSED'} = 1;
      last;
    }
    elsif( $_ =~ /^\s*Frame\s*=\s*([-\+]\d+)/ ) {
	$frame = $1;
    }

    else {
      push @hspline, $_;           #      store the query line
      $nextline = <$FH> ;
# Skip "pattern" line when parsing PHIBLAST reports, otherwise store the alignment line
      my $l1 = ($nextline =~ /^\s*pattern/) ? <$FH> : $nextline;
#     my $l1 =  push @hspline, $l1; # grab/store the alignment line
      push @hspline, $l1; # store the alignment line
      my $l2 = <$FH>; push @hspline, $l2; # grab/store the sbjct line
    }
  }
  
  #########################
  # parse alignment lines #
  #########################
  my ($ql, $sl, $as) = ("", "", "");
  my ($qb, $qe, $sb, $se) = (0,0,0,0);
  my (@QL, @SL, @AS); # for better memory management
  
  for(my $i=0;$i<@hspline;$i+=3) {
    # warn $hspline[$i], $hspline[$i+2];
    $hspline[$i]   =~ /^Query:?\s+(\d+)\s*(\S+)\s+(\d+)/;
    $ql = $2; $qb = $1 unless $qb; $qe = $3;
    
    my $offset = index($hspline[$i], $ql);
    $as = substr($hspline[$i+1], $offset, CORE::length($ql));
    
    $hspline[$i+2] =~ /^Sbjct:?\s+(\d+)\s*(\S+)\s+(\d+)/;
    $sl = $2; $sb = $1 unless $sb; $se = $3;
    
    push @QL, $ql; push @SL, $sl; push @AS, $as;
  }
  
  ##################
  # the HSP object #
  ##################
  $ql = join("", @QL);
  $sl = join("", @SL);
  $as = join("", @AS);
  my $hsp = new Bio::EnsEMBL::Analysis::Tools::BPlite::HSP('-score'=>$score, 
					'-bits'=>$bits, 
					'-match'=>$match,
					'-positive'=>$positive, 
					'-percent'=>$percent,
					'-p'=>$p,
					'-queryBegin'=>$qb, 
					'-queryEnd'=>$qe, 
					'-sbjctBegin'=>$sb,
					'-sbjctEnd'=>$se, 
					'-querySeq'=>$ql, 
					'-sbjctSeq'=>$sl,
					'-homologySeq'=>$as, 
					'-queryName'=>$self->{'PARENT'}->query,
					'-sbjctName'=>$self->{'NAME'},
					'-queryLength'=>$self->{'PARENT'}->qlength,
					'-sbjctLength'=>$self->{'LENGTH'},
					'-frame'   => $frame);
  return $hsp;
}

1;
