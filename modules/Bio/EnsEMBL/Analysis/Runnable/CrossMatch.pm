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

Bio::EnsEMBL::Analysis::Runnable::CrossMatch - 

=head1 DESCRIPTION

Module to provide object-oriented access to Phil Green's B<cross_match>
Smith-Waterman alignment program.

=head1 SYNOPSIS

	use CrossMatch;
	
	# Create a factory object
	$matcher = CrossMatch::Factory->new( '/nfs/disk2001/this_dir',
	                                     '/home/jgrg/that_dir' );
	$matcher->minMatch( 80 );

        # process full crossmatch alignments
        $matcher->alignments;

	# Match two fasta fomat sequence files, generating
        # a match object
	$cm = $matcher->crossMatch( 'dJ334P19.seq', 'cB49C12.aa' );

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::CrossMatch;
use warnings ;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(info verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::RunnableI

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

# new() is written here 

sub new {

  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  bless $self,$class;
  
  $self->{'_fp_array'} =[];
  
  my ($seq1,
      $seq2,
      $program,
      $workdir,
      $minscore,
      $minmatch,
      $masklevel) = rearrange([qw(
				  SEQ1 
				  SEQ2 
				  PROGRAM
				  WORKDIR 
				  MINSCORE
				  MINMATCH 
				  MASKLEVEL)],@args);
  
  if( !defined $seq1 || !defined $seq2 ) {
    $self->throw("Must pass in both seq1 and seq1 args");
  }
  
  my $queryfile = $self->queryfile();
  my $targetfile = $queryfile."target";
  
  $self->write_seq_file($seq1,$queryfile);
  $self->write_seq_file($seq2,$targetfile); print "$seq2 and $targetfile\n";
  
  if( $workdir) { 
    $self->workdir($workdir); 
  } else {
    $self->workdir("/tmp");
  }

  $minmatch = 30 unless (defined $minmatch); # $minmatch = $minscore compatatible with cvs version 1.7
  $minscore = $minmatch unless (defined $minscore); 
  $masklevel = 101 unless (defined $masklevel);
  
  if( defined $minmatch ) {
    $self->minmatch($minmatch);
  }
  if( defined $minscore ) {
    $self->minmatch($minscore);
  }
  if( defined $masklevel ) {
    $self->masklevel($masklevel);
  }
  
  #my $path = $self->locate_executable("cross_match");
  
  #$self->program($path);

  my $options = "-minmatch $minmatch -minscore $minscore -masklevel $masklevel -alignments $queryfile $targetfile ";
  
  if ($options){
    $self->options($options);
  }
    
  # set stuff in self from @args
  return $self;
}

sub run {
  my ($self, $dir) = @_;

  if (!$dir) {
    $dir = $self->workdir;
  }
  $self->checkdir($dir);
  my $filename = $self->write_seq_file();
  $self->files_to_delete($filename);
  $self->files_to_delete($filename."target");
  $self->files_to_delete($filename.".log");
  $self->files_to_delete($self->resultsfile);
  $self->run_analysis();
  $self->parse_results;
  $self->delete_files;
  return 1;
}

=head2 run_analysis

Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::ExonerateArray
Arg [2]   : string, program name
Function  : constructs a commandline and runs the program passed
  in, the generic method in Runnable isnt used as ExonerateArray doesnt
  fit this module
  Returntype: none
  Exceptions: throws if run failed because sysetm doesnt
  return 0 or the output file doesnt exist
 Example   :

=cut


sub run_analysis {

  my ($self, $program) = @_;
  
  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable CrossMatch::run_analysis ")
    unless($program && -x $program);
  my $cmd = $self->program." ";
  $cmd .= $self->options." " if($self->options);
  $cmd .= ">".$self->resultsfile;
  print "Running analysis ".$cmd."\n";
  
  system($cmd) == 0 or throw("FAILED to run ".$cmd." CrossMatch::run_analysis ");

  if(! -e $self->resultsfile){
    throw("FAILED to run CrossMatch on ".$self->queryfile." ".
          $self->resultsfile." has not been produced ".
          "CrossMatch:run_analysis");
  }
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::ExonerateArray
  Arg [2]   : string, filename
  Function  : open and parse the results file into misc_features
  features
  Returntype: none
  Exceptions: throws on failure to open or close output file
  Example   :

=cut


sub parse_results{

  my ($self, $results) =  @_;
  if(!$results){
    $results = $self->resultsfile;
  }

  my (@res,$ialign,@align);

  open( CROSS_MATCH, $results ) || throw("FAILED to open ".$results." CrossMatch::parse_results");
  
  while (<CROSS_MATCH>){
    
    info ($_) ;
    
    # process alignment lines if requested
    # print STDERR "Processing....$_";
    
    if(/^(\w*)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)$/){  ###this is alignment line
      if($2 eq $res[4] && $ialign==1){
	$align[0].=$4;
	$align[2]=$1;
	$ialign=2; 
      }elsif(($2 eq $res[8] || $2 eq $res[9]) && $ialign==2){
	$align[1].=$4;
	$align[3]=$1;
	$ialign=1;
      }else{
	die "alignment parsing error in Crossmatch.pm\n  $_";
      }
    }
    
    # this is used to exclude all except alignment summary lines
    next unless /\(\d+\).*\(\d+\)/;
    
    # Remove parentheses and asterisks
    tr/()\*//d;
    
    # save previous hit, with alignments
    &_save_hit($self,\@res,\@align);
    
    @res = split;
    $ialign=1;
  }

  #print STDERR "About to process...\n";
  &_save_hit($self,\@res,\@align);
  #print STDERR "saved hits \n";
      
  undef @res;
  undef @align;

  # Check exit status of cross_match
  unless (close CROSS_MATCH) {
    my $error = $! ? "Error from cross_match: $!"
      : "Error: cross_match exited status $?";
    warning( "$error\nCommand: 'cross_command'");
  }
  
# Pre-rearranged fields
# 0    1    2    3    4        5     6     7       8 9        10      11    12

# 98   0.00 0.00 0.00 130N4    1     104   84065   W 92M18    68800   68903 0
# 98   0.00 0.00 0.00 92M18    68800 68903 0       W 130N4    1       104   84065
# 8251 0.00 0.00 0.00 130N4    1     84169 0       W 130N4    1       84169 0
# 103  0.00 0.00 0.00 CFAT5    20771 20874 (0)     W A1280    1       104   (22149) | W Bs
# 103  0.00 0.00 0.00 CFAT5    20771 20874 (0)     C A1280.RC (0)     22253 22150   | C Ar Br
# 103  0.00 0.00 0.00 CFAT5.RC 1     104   (20770) C A1280    (22149) 104   1       | C As Bs
# 26355  1.37 0.42 0.36  Em:AP000350    133977 162328 (0)  C Em:AP000352   (125738) 28369     1

# 120 16.53 0.00 0.00  bK363A12    32474 32715 (50)  C cE129H9   (4343) 32827 32586  
# 100  0.00 0.00 0.00  bK363A12    32666 32765 (0)    cE129H9        1   100 (37070) *
#  * indicates that there is a higher-scoring match whose domain partly includes the domain of this match.


#process alignments above score
  foreach my $fp ( $self->fp($self->minmatch) ) {
    my ($seq1,$seq2,$score,$start,$end,$hstart,$hend) = split(/:/,$fp);
  
    my ($strand,$hstrand,$swap);

    if ( $start > $end ) {
      $strand = -1;
      $swap = $start;
      $start = $end;
      $end = $swap;
    } else {
      $strand = 1;
    }

    if ( $hstart > $hend ) {
      $hstrand = -1;
      $swap = $hstart;
      $hstart = $hend;
      $hend = $swap;
    } else {
      $hstrand = 1;
    }


    $fp = Bio::EnsEMBL::FeaturePair->new();
    #print STDERR "Processing FP with $start-$end to $hstart-$hend\n";

    $fp->start($start);
    $fp->end($end);
    $fp->strand($strand);
    $fp->seqname($seq1);
    $fp->hstart($hstart);
    $fp->hend($hend);
    $fp->hstrand($hstrand);
    $fp->hseqname($seq2);
    $fp->score($score);
       
    $self->add_fp($fp);
  }
}

=head2 output

Title   : output
Usage   :
Function:
Example :
Returns : 
Args    :


=cut

sub output{
  my ($self,@args) = @_;
  
  return @{$self->{'_fp_array'}};
}



=head2 add_fp

Title   : add_fp
Usage   :
Function:
Example :
Returns : 
Args    :


=cut

sub add_fp{
  my ($self,@args) = @_;
  
  push(@{$self->{'_fp_array'}},@args);
}

=head2 minmatch

  Title   : minmatch
    Usage   : $obj->minmatch($newval)
  Function: could set and return the minmatch value
    Example : 
    Returns : value of minmatch option used by crossmatch
    Args    : newvalue (optional)


=cut

sub minmatch {
  my ($obj,$value) = @_;
  if ( defined $value) {
    $obj->{'minmatch'} = $value;
  }
  return $obj->{'minmatch'};

}

=head2 minscore

  Title   : score
    Usage   : $obj->score($newval)
  Function: could set and return the score value
    Example : 
    Returns : value of minscore option used by crossmatch
    Args    : newvalue (optional)


=cut
sub minscore {
  my ($obj,$value) = @_;
  if ( defined $value) {
    $obj->{'minscore'} = $value;
  }
  return $obj->{'minscore'};

}

=head2 masklevel

  Title   : masklevel
    Usage   : $obj->masklevel($newval)
  Function: could set and return the masklevel value
    Example : 
    Returns : value of masklevel option used by crossmatch
    Args    : newvalue (optional)


=cut

sub masklevel {
  my ($obj,$value) = @_;
  if ( defined $value) {
    $obj->{'masklevel'} = $value;
  }
  return $obj->{'masklevel'};

}

sub _save_hit{
  my($self,$rares,$raalign)=@_;

  # only save if something to save
  if (@$rares) {
    
    # parsing to deal different numbers of items
    # if 13 then second sequence is complement
    # if 12 then its not, so insert 'W'
    my @nres;
    if (@$rares == 13) {
      @nres=( @$rares[0..9, 12, 11, 10]       );
    } elsif (@$rares == 12) {
      @nres=( @$rares[0..7], 'W', @$rares[8..11] );
    }
    
    # alignment is stored in a directional fashion: x->y maps to a->b where y>x
    my $raw_align;
    if (@$raalign) {
      # reverse if required
      my $st1=$nres[5];
      my $en1=$nres[6];
      my $st2;
      my $en2;
      my $dirl;
      if ($$raalign[2] eq 'C') {
	$$raalign[3]='C';
	$$raalign[0]=reverse($$raalign[0]);
	$$raalign[1]=reverse($$raalign[1]);
	$dirl='C';
	$st2=$nres[11];
	$en2=$nres[10];
      } else {
	$st2=$nres[10];
	$en2=$nres[11];
      }

      # check length
      my $l1=length($$raalign[0]);
      my $l2=length($$raalign[1]);
      if ($l1!=$l2) {
	print "Lengths of alignment are different [$l1,$l2]\n$$raalign[0]\n$$raalign[1]\n";
	print join(',',@nres)."\n";
	throw( "failed at CrossMatch :_save_hit");
      }

      # walk along sequence in blocks
      $raw_align="$st1:$dirl$st2";
      my $seq1=$$raalign[0];
      my $seq2=$$raalign[1];
      {
	my($s1a,$s1b,$s1c,$s2a,$s2b,$s2c);
	if ($seq1=~/^([^\-]+)(\-+)(\S+)$/) {
	  ($s1a,$s1b,$s1c)=($1,$2,$3);
	} else {
	  $s1a=$seq1;
	  $s1b=$s1c='';
	}
	if ($seq2=~/^([^\-]+)(\-+)(\S+)$/) {
	  ($s2a,$s2b,$s2c)=($1,$2,$3);
	} else {
	  $s2a=$seq2;
	  $s2b=$s2c='';
	}
	#print STDERR "heads are ",substr($s1a,0,10),";",substr($s1b,0,10),":",substr($s2a,0,10),";",substr($s2b,0,10),"\n";

	# escape if no more gaps
	next if(length($s1c)==0 && length($s2c)==0);
	# do shortest first
	my $lab;
	if ( length($s1a.$s1b) == 0 ||
	     length($s2a.$s2b) == 0 ) {
	  print STDERR "Dodgy alignment processing catch! Bugging out\n";
	  next;
	}
	if (length($s1a.$s1b)<length($s2a.$s2b)) {
	  # update seq1
	  $lab=length($s1a.$s1b);
	  $st1+=length($s1a);
	  #print STDERR "st1 is $st1 with $lab\n";
	  $seq1=$s1c;
	  #print STDERR "New head is ",substr($seq1,0,10),"\n";
	  # update seq2
	  $seq2=~s/^\S{$lab}//;
	  #print STDERR "new $lab.. seq2 head is ",substr($seq2,0,10),"\n";
	} else {
	  # update seq2
	  $lab=length($s2a);
	  $seq2=$s2c;
	  # update seq1
	  my $l2ab=length($s2a.$s2b);
	  $seq1=~s/^\S{$l2ab}//;
	  $st1+=$l2ab;
	}
	if ($dirl eq 'C') {
	  $st2-=$lab;
	} else {
	  $st2+=$lab;
	}
	$raw_align.=",$st1:$dirl$st2";
	redo;
      }
      $raw_align.=",$en1:$dirl$en2";

    }
    $self->hit(@nres,$raw_align);

    # clear data
    @$rares=();
    @$raalign=();
  }
}
    
# Store data from a hit

sub hit {
  my $self = shift;
  
  if (@_ == 14) {
    push( @{$self->{'hit'}}, [@_]);
    push( @{$self->{'active'}}, $#{$self->{'hit'}} );
  } else {
    warning( "Bad number of elements (", scalar @_, ") in '@_'");
  }
}

# Create access methods to access the data from the matches

BEGIN {
    my( %fields );
    {
        my $i = 0;
        %fields = map {$_, $i++} qw( score pSub pDel pIns
                                     aName aStart aEnd aRemain
                                     strand
                                     bName bStart bEnd bRemain raw_align);
    }

    foreach my $func (keys %fields) {
        no strict 'refs';

        my $i = $fields{ $func };

        *$func = sub {
            my( $match, @rows ) = @_;

            if (wantarray) {
                # Extract the requested values
                unless (@rows) {
                    # Get a list of all the row indices
                    @rows = @{$match->{'active'}}; 
                }

               # Make a vertical slice through the hits
                return map $match->{'hit'}[$_][$i], @rows;
            } else {
                # Return just one value in scalar context
                if (defined $rows[0]) {
                    # Get field from row requested
                    return $match->{'hit'}[$rows[0]][$i];
                } else {
                    # Just get value from first row in active list
                    my $row = $match->{'active'}[0]; 
                    return $match->{'hit'}[$row][$i];
                }
            }
        }
    }
}




# output's full featurepair list
sub fp {
  my( $self, $minscore ) = @_;
  # loop over all rows
  my @hits; 
  foreach my $row (@{$self->{'active'}}) {
    
    my $score=$self->score($row); 
    next if($score<$minscore);
    my $aname=$self->aName($row);
    my $bname=$self->bName($row);
    my @sted=&_trans_fp($self->raw_align($row));
    foreach my $sted2 (@sted) {
      push(@hits,join(":",$aname,$bname,$score,@$sted2));
    }
  }
  return @hits;
}

sub _trans_fp{
  my($align)=@_;
  my($st1p,$st2p);
  $st1p=$st2p=0;
  my @hits;
  foreach my $pair (split(',',$align)) {
    my $dirl;
    my($st1,$st2)=split(':',$pair);
    if ($st2=~/C(\d+)/) {
      $st2=$1;
      $dirl='C';
    }

    # only output if previous saved
    if ($st1p) {
      # calculate ends
      my $l1=$st1-$st1p;
      my $l;
      if ($dirl) {
	$l=$st2p-$st2;
      } else {
	$l=$st2-$st2p;
      }
      # if exact match, must be end, so need full length
      if ($l1==$l) {
	$l++;
      } elsif ($l1<$l) {
	$l=$l1;
      }
      $l1=$l-1;
      my $ed1p=$st1p+$l1;
      my $ed2p;
      if ($dirl) {
	$ed2p=$st2p-$l1;
      } else {
	$ed2p=$st2p+$l1;
      }
      push(@hits,[$st1p,$ed1p,$st2p,$ed2p]);
    }
    ($st1p,$st2p)=($st1,$st2);
  }
  # if higher than end then no match
  return @hits;
}

1;


  __END__

=head1 NAME - CrossMatch.pm

=head1 DESCRIPTION

Module to provide object-oriented access to Phil Green's B<cross_match>
Smith-Waterman alignment program.

=head1 SYNOPSIS

	use CrossMatch;
	
	# Create a factory object
	$matcher = CrossMatch::Factory->new( '/nfs/disk2001/this_dir',
	                                     '/home/jgrg/that_dir' );
	$matcher->minMatch( 80 );

        # process full crossmatch alignments
        $matcher->alignments;

	# Match two fasta fomat sequence files, generating
        # a match object
	$cm = $matcher->crossMatch( 'dJ334P19.seq', 'cB49C12.aa' );

=head1 FUNCTIONS

=over 4

=item new

Create a new CrossMatch factory object.  Optional arguments to new is a
list of directories to search (see B<dir> below).

=item crossMatch

Do a cross_match on the two B<fasta formatted> sequence files supplied
as agruments, returning a B<CrossMatch> object.  The B<CrossMatch>
object records the parameters used in the match, and the match data
returned my cross_match (if any).

I<Note:> The two fasta files may each contain multiple sequences.

=item dir

Change the list of directories in which B<crossMatch> searches for
files.  The list of directories supplied to B<dir> completely replaces
the existing list.  Directories are searched in the order supplied.

This list defaults to the current working directory when the
B<CrossMatch::Factory> object is created.

=item extn

A convenience function which allows you to supply a list of filename
extensions to be considered by B<crossMatch> when finding a file to
pass to cross_match.  For example:

	$matcher->extn( 'seq', 'aa' );
	$cm = $matcher->crossMatch( 'dJ334P19', 'cB49C12' );

B<crossMatch> will look for files "dJ334P19.seq", "dJ334P19.aa".  If
no B<extn>s had been set, then only "dJ334P19" would have been
searched for.

=item minMatch

Set the B<minmatch> parameter passed to cross_match.  Defaults to 50.

=item minScore

Set the B<minscore> parameter passed to cross_match.  Defaults to 50.

=item maskLevel

Set the B<masklevel> parameter passed to
cross_match.  Defaults to 101, which displays all
overlapping matches.

=item alignments

Causes the full crossmatch alignments to be parsed and stored in the
Crossmatch object.  These can be accessed by the a2b and fp methods.

=back

=head1 CrossMatch

Data stored in B<CrossMatch> objects, generated by a
B<CrossMatch::Factory>, can be retrieved with a variety of queries.
All the matches are stored internally in an array, in the order in
which they were generated by cross_match.  You can apply filters and
sorts (sorts are not yet implemented) to the object, which re-order or
hide matches.  These operations save their changes by altering the
active list, which is simply an array of indices which refer to
elements in the array of matches.  The data access funtions all
operate through this active list, and the active list can be reset to
show all the matches found by calling the B<unfilter> method.

=over 4

=item SYNOPSIS

	my $numberOfHits = $cm->count;
	foreach my $hit ($cm->list) {
	    print $cm->score($hit);
        }

	$firstName = $cm->aName;
	@allscores = $cm->score();
	@someSores = $cm->score(0,3,4);


=item FUNCTIONS

=over 4

=item count

Returns a count of the number of hits in the active list

=item list

Returns a list of array indices for the current list.

=item filter

	$cm->filter(\&CrossMatch::endHits);

Takes a reference to a subroutine (or an anonymous subroutine), and
alters the active list to contain only those hits for which the
subroutine returns true.  The subroutine should expect to be passed a
B<CrossMatch> object, and an integer corresponding to the index of a
match.  Applying a series of filters to the same object removes
successively more objects from the active list.

=item unfilter

Resets the active list to include all the hits, in the order in which
they were generated by cross_match.  Returns a list of the active
indices.

=item FILTERS

=over 4

=item endHits

A filter which returns true if one of the sequnces in a hit has its
end in the hit.

=back

=back

=item PARAMETERS

The parameters used by cross_match when doing the match can be
retrieved from the resulting Match object.  The two sequences matched
are labelled I<a> and I<b>, where I<a> is the sequence from the first
file passed to B<crossMatch>, and I<b>, the sequence from the second
file.  For example:

	$path_to_a = $cm->aFile;

retrieves the path of the file supplied as the first argument to
cross_match.

=over 4

=item aFile bFile

The full paths to the files used for sequences I<a> and I<b>.

=item minM

The minmatch parameter used.

=item minS

The minscore parameter used.

=item mask

The masklevel parameter used.

=back

=item DATA FIELDS

Syntax:

	$field = $match->FIELD();
	$field = $match->FIELD(INDEX);
	@fields = $match->FIELD();
	@fields = $match->FIELD(LIST);

Examples:

	$firstScore = $match->score(0);
	@aEnds = $match->aEnd();

These methods provide access to all the data fields for each hit found
by cross_match.  Each of the fields is described below.  In scalar
context, a single data field is returned, which is the field from the
first row in the active list if no I<INDEX> is given.  In array
context it returns a list with the I<FIELD> from all the hits in the
active list when called without arguments, or from the hits specified
by a I<LIST> of indices.

=over 4

=item score

The Smith-Waterman score of a hit, adjusted for the complexity of the
matching sequence.

=item pSub

The percent substitutions in the hit.

=item pDel

The percent deletions in the hit (sequence I<a> relative to sequence I<b>).

=item pIns

The percent insertions in the hit (sequence I<a> relative to sequence I<b>).

=item aName bName

ID of the first sequence (sequence I<a>) and second seqeunce (I<b>) in
the match, respecitvely.

=item aStart bStart

Start of hit in I<a>, or I<b>, respectively.

=item aEnd bEnd

End of hit in I<a>, or I<b>, respectively.

=item aRemain

Number of bases in I<a> after the end of the hit.  ("0" if the hit
extends to the end of I<a>.)

=item bRemain

The equivalent for sequence I<b> of aRemain, but note that if the
strand of I<b> which matches is the reverse strand, then this is the
number of bases left in I<b> after the hit, back towards the beginning
of I<b>.

=item strand

The strand on sequence I<b> which was matched.  This is "B<W>" for the
forward (or parallel or B<Watson>) strand, and "B<C>" for the reverse
(or anti-parallel or B<Crick>) strand.

=item a2b

Returns coordinate in sequence I<b> equivalent to value in sequence
I<a> passed to method.  If no corresponding base, returns I<-1>.

=item fp

Returns an array of strings contining full list of ungapped alignment
fragment coordinates, filtered by score value passed to method.

=back

=item QUERYING BY SEQUENCE NAME

These methods allow access to the B<start>, B<end>, and B<remain>
fields for those occasions where you know the name of your sequence,
but don't necessarily know if it is sequence I<a> or I<b>.  Syntax and
behaviour is the same as for the I<DATA FIELDS> functions, but the
first argument to the function is the name of the sequence you want to
get data about.

I<Note:> These methods perform a filtering operation, reducing the
active list to only those hits which match the name given (in either
the I<a> or I<b> columns).  You'll therefore need to B<unfilter> your
match if you need data from hits which have now become hidden.

If the sequence name is found in both columns I<a> and I<b>, then a
warning is printed, and the I<a> column is chosen.

For example, suppose you have matched two fasta files containing only
the sequences cN75H12 and bK1109B5 (in that order), then the following
calls retrieve the same data:

	@ends = $match->aEnd();
	@ends = $match->end('cN75H12');

	$start = $match->bStart(0);
	$start = $match->start('bK1109B5', 0);

A warning is printed to STDERR if a match contains hits from the name
supplied in both columns, and only hits from the a column are
returned.

=over 4

=item start end remain

Access the aStart or bStart, aEnd or bEnd, and aRemain or bRemain
fields, depending upon wether the name supplied matches the aName or
bName fields respectively.

=back

=head1 BUGS


=head1 AUTHOR

B<James Gilbert> email jgrg@sanger.ac.uk

B<Tim Hubbard> email th@sanger.ac.uk (alignment processing extensions)

=cut





