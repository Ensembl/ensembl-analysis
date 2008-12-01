
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ExonerateProbe

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::ExonerateProbe->new(
     -query_seqs     => \@q_seqs,
     -query_type     => 'dna',
     -target_seqs    => \@t_seqs,
     -options        => $options,
    );

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module handles a specific use of the Exonerate (G. Slater) program, to
align probes to a target genome. (The resulting alignments will be stored in an 
ansembl Funcgen db as Bio::EnsEMBL::ProbeFeature objects.)

NOTE: the ProbeFeature objects refer to Probe id's, and they in turn 
refer to ArrayChip and Array id's. Hence, Arrays, ArrayChips and Probes
should be pre-loaded into the ensembl db: there are separate RunnableDB
/RunnableDB's to do this from the Affymetrix data sets r use the EFG Importer 
to load other arrays e.g. Nimblegen or Sanger. 
This runnable just creates fake Probes in order to create reasonable-looking
affy features???????

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateProbe;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;
#use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseExonerate);

sub new {
  my ( $class, @args ) = @_;
  #my $self = $class->SUPER::new(@args);

  #slightly rearranged order as we want to pass some different defaults to BaseExonerate

  my ($max_mismatches, $mapping_type, $basic_options) = rearrange(['max_mismatches', 'mapping_type', 'basic_options'], @args);
  
  my %basic_opts;

  #We could change mapping_type to same_strand hits filter
  #As this is the only thing we're using it for here


  if(! defined $basic_options){
	#parse result depends on the output format options
	#only override if you intend overload or rewrite the parse_results method.
	#Now let's reset the default BaseExonerate options to remove vulgar and add scores
	#RESULT: 9382225 0 25 + ENSMUST00000060050 1076 1051 - 116 96.00 25 1215 . 1 scores:0:5:5:5:-4:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:0:
	$basic_opts{'-basic_options'} = "--showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %em scores:{%Ps:}\\n\" ";
  }

  my $self = $class->SUPER::new(@args, %basic_opts);

  if(! defined $max_mismatches){
	throw("Must provide a -max_mismatches parameter e.g. 0, 1 or 2");
  }

  if(! (defined $mapping_type && $mapping_type =~ /genomic|transcript/)){
	throw("Must provide a valid -mapping_type parameter e.g. genomic or transcript");
  }

  #reset to the new default exonerate
  if (not $self->program) {
	#This is an architecture specific build path
    $self->program('/lustre/work1/ensembl/gs2/local/x86_64/bin/exonerate');
  }
 


  #if(! (defined $filter_method && (ref($filter_method) ne 'CODE' || $self->can($filter_method)))){
#	throw('You must pass a -filter_method name or CODEREF to filter the ProbeFeatures');
#  }

  ##Set code ref or pointer to internal method
  #We will always have explicitly pass self to the coderef
  

  #$self->{'filter_method'} = (ref($filter_method) eq 'CODE') ? $filter_method : $self->can($filter_method);
  $self->{'max_mismatches'} = $max_mismatches;
  $self->{'mapping_type'}   = $mapping_type;
  


  return $self;
}


sub max_mismatches{
  return $_[0]->{max_mismatches};
}

sub mapping_type{
  return $_[0]->{mapping_type};
}


#
# Implementation of method in abstract superclass
#
sub parse_results {
  my ( $self, $fh ) = @_;
  
  my @features;
  #my $filter_method = $self->{'match_rules'};
  #my $filter_method = $self->filter_method;
  #No, now uses code ref to allow definition in config

  #print "Parsing results from fh ".Data::Dumper::Dumper($fh)."\n";

  
  my ($tag, $probe_id, $q_start, $q_end, $q_strand, 
	 $t_id, $t_start, $t_end, $t_strand, $score, $tscore,
	 $perc_id, $q_length, $t_length, $mismatch_count, $scores,
	$match_length, $align_mismatch, $total_mismatches, $tmp);

  my $max_mismatches = $self->max_mismatches;
  my $mapping_type   = $self->mapping_type;


  while (<$fh>){
    #print STDERR $_ if $self->_verbose;
    next unless /^RESULT:/;
    chomp;
    
	#--showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %g %V\\n\"	
	#--bestn 100 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 14 --dnawordlimit 11 
	
	#Vulgar blocks are also report in in-between coords! So need to add 1(only to start???)
	
	#Shows  the  alignments  in "vulgar" format.  Vulgar is Verbose Useful Labelled Gapped Alignment Report, This format also starts with the same 9 fields as sugar output (see above), and is followed by a series of <label, query_length, target_length> triplets.  The label may be one of the following:

#              M      Match
#              G      Gap
#              N      Non-equivalenced region
#              5      5' splice site
#              3      3' splice site
#              I      Intron
#              S      Split codon
#              F      Frameshift

	#Should always be M probe_length target_length(which would be equal to probe length for ungapped alignment))

	#Can we have a way of dynamically linking the ryo string to a particular method?
	#We shouldn't need different formats!!!

	#"RESULT: %S %pi %ql %tl %g %em scores:{%Ps:} %V \n"
	#RESULT: 9380677 0 25 + ENSMUST00000095075 3120 3145 + 125 100.00 25 4006 . 0 M 25 25  :0:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:0

	#RESULT: 9380676 0 25 + ENSMUST00000095075 3047 3072 + 116 96.00 25 4006 . 1 M 25 25  :0:5:5:5:5:-4:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:0

	#we don't have a M 24 24 here, but that would be the other possibility.

	#So we use %em or %V to figure out if we don't have a perfect match, than we use %em or %V and the info string to figure out where the mismatches are
	#< 25M effectively means mismatches at the start of end of a query sequence


    ($tag, $probe_id, $q_start, $q_end, $q_strand, 
	 $t_id, $t_start, $t_end, $t_strand, $score, 
	 $perc_id, $q_length, $t_length, $mismatch_count, $scores) = split;
	
	#$gene_orientation always .?
	#, $match_type, $query_match_length, $target_match_length, @rest_of_vulgar #vulgar blocks	
	#query strand should always be +?
	
	#Do we really need vulgar? 
	#We never get anything than M N N for ungapped alignments?
	
	#if(@rest_of_vulgar){
	#  throw ("There is more than a simple match ('M') vulgar output:\n".
	#		 "$match_type $query_match_length $target_match_length @rest_of_vulgar\n".
	#		 "for affy tag $probe_id mapped to region $t_id \n");
	#}
	
	#if(!($match_type eq 'M')){
	#  throw "I have received a starting Vulgar symbol $match_type which is not a match!"; 
	#}
  


	#1 25 + ENSMUST00000115314 1593 1617 + 120 100.00 25 3173 . 0 M 24 24
	#0 25 + ENSMUST00000109902 2318 2343 + 116 96.00 25 2449 . 1 M 25 25 

	#because of the 'in-between' coordinates.???
	$t_start += 1;	
	#What about t_end?
	#Not all seem to be inbetween:
	#0 25 -represents 1-25 of a 25mer
	#Only seems to affect start!
	#what about query? We also need to add to query if we were using it.
	
	
	if(!($probe_id =~ /\d+/)){
	  throw "Probe headers MUST be the internal db ids of the Probes for this parser to work!\n";
	}
  


	#We could filter here for -ve strand on transcripts, unless we want to capture this infor for antisense transcription?
	#But this may not follow cdna, so genomic alignment is fine for this.
	#Would -ve to -ve alignments be valid? Yes, but they will never be report over +ve to +ve?
	#SHouldn't this be filtered in the RunnableDB ProbeAlign and leave this to be generic?

	if($mapping_type eq 'transcript' && $t_strand eq '-'){

	  if($q_strand eq '-'){
		warn "Found -ve to -ve mapping for probe $probe_id, need to account for this";
	  }
	  else{
		next;
	  }
	}


	#RESULT: 9379131 0 25 + chromosome:NCBIM37:16:1:98319150:1 33212985 33213010 + 116 96.00 25 98319150 1 scores:0:5:5:5:5:5:5:5:5:5:5:-4:5:5:5:5:5:5:5:5:5:5:5:5:5:5:0: at /nfs/acari/nj1/src/ensembl-analysi


	#filter_emthod is FILTER_METHOD coderef set in the RunnableDB analysis config hash

	#We don't need this filter method any more as we can simply us the query length - match_length + mismatch count!!!
	#Will query start end always have full query length?
	#How will we know full query length if we have different sized probes?

	#my $mismatch_count = $self->filter_method($query_match_length, $q_length, $score);#, $cigar_line);


	$match_length = $q_end - $q_start;
	$align_mismatch = $q_length - $match_length;

	#We also need to alter $t_end or $t_start if we have an align_mismatch!?
	#Should we be reporting these align mismatches at the end?
	#We simply need to know them for filtering
	#We have no way of knowing the full query length if we don't pass the full length
	#cigar line to the Runnable DB
	#Storing the hanging mismatches is not incorrect, just non-sandard
	#We would have to pass the probe length as part of the ProbeFeature
	#to enable storage of query % ID
	#Leave for now as we have no way of passing the query start stop to ProbeAlign
	#hence we would lose where the flanking mismatch was.

	$total_mismatches = $mismatch_count + $align_mismatch;
	my @soft_cigar_line;
	
	#5' unaligned
	if($align_mismatch){
	  push @soft_cigar_line, $q_start.'m' if $q_start;#set this to the value of start if not 0
	  #We want to subtract from start if +ve hit
	  #else we want to add to end if -ve strand, end is actually start in ensembl terms
	  $t_start = ($t_strand eq '+') ? ($t_start - $q_start) : ($t_end + $q_start);
	}

	#mismatches
	if($mismatch_count){
	  my @scores = split/:/, $scores;
	  #remove scores:0
	  shift @scores;
	  shift @scores;
	  #remove last :0
	  pop @scores;
	  #Is this quicker that a splice?

	  my $cnt        = 0;
	  my $prev_score = 5;#always starts with a match otherwise it wouldn't have been reported.

	  while($tscore = shift @scores){
		
		if($tscore == $prev_score){
		  $cnt += 1;
		}else{
		  $tmp = ($prev_score == 5) ?  $cnt.'M' :  $cnt.'m';
		  push @soft_cigar_line, $tmp;
		  $cnt = 1;
		  $prev_score       = $tscore;
		}
	  }

	  #handle last align length
	  $tmp = ($prev_score == 5) ? $cnt.'M' : $cnt.'m';
	  push @soft_cigar_line, $tmp;
	}

	#3' unaligned
	if($align_mismatch != $q_start){
	  #Add end mismatch if
	  #not accounted for by 5' mismatch
	  my $three_mismatch = ($q_length - $q_end);
	  push @soft_cigar_line, $three_mismatch.'m';

	  #Either add to end for +ve or subtract from ensembl start(which is end) for -ve
	  $t_end = ($t_strand eq '+') ? ($t_end + $three_mismatch) : ($t_start - $three_mismatch);
	}



	

	#Do we need to account for this?!!
	#i.e. do we need to swap the t_start t_end?
	#9381894 0 25 + ENSMUST00000039541 2799 2774 - 
	#We don't want -ve strand matches from transcripts?!
	#Only valid -ve strand matches will be not cross introns
	#So we are not interested in them for transcript mapping.
	#Only for genomic mapping

	if($total_mismatches <= $max_mismatches){

	  if($q_strand eq '+'){
		if($t_strand eq '+'){
		  $t_strand = 1;

		  #Alter start/end if we have flanking mismatches here or above?


		}elsif($t_strand eq '-'){
		  $t_strand = -1;
		  ($t_start, $t_end) = reverse($t_start, $t_end);

		  #Reverse the cigar line
		  #As match scores are reported wrt to strand of hit
		  #And we want cigar line wrt to +ve strand of hit
		  @soft_cigar_line = reverse(@soft_cigar_line);

		}else{
		  throw "Unrecognised target strand symbol: $t_strand\n";
		}
	  }
	  elsif($q_strand eq '-'){

		throw('We have found a -ve query match');
		#Exonerate only reports +ve query strand matches??? Really?
		#Are we even interested in -ve query strand matches?

		#if($t_strand eq '-'){
		#  $t_strand = 1;
		#}elsif($t_strand eq '+'){
		#  $t_strand = -1;
		#}else{
		#  throw "unrecognised target strand symbol: $t_strand\n";
		#}
	  }else{    
		throw "unrecognised query strand symbol: $q_strand\n";
	  }

	  push @features, new Bio::EnsEMBL::Funcgen::ProbeFeature
		(
		 #-probe => $probe,
		 -probe_id => $probe_id,
		 -start => $t_start,
		 -end => $t_end,
		 -strand => $t_strand,
		 -mismatchcount => $total_mismatches,
		 -cigar_line => join(':', @soft_cigar_line) || $match_length.'M',
		 -seqname => $t_id,
		);

	  # attach the slice name onto the feature: let the runnabledb
	  # sort out whether it's valid.
	  #$feature->seqname($t_id);
  

	}
	else{
      warn "Feature from probe :$probe_id doesnt match well enough\n";
    }
  }

  return \@features;
}

#we need the full unmatched query length to calculat the mismatches!
#Can't assume will always be 25 or 50
#This should be calculated based on % id?

#This can now coderef'd in the analysis config hash
#Altho this obfuscates this module slighlty, it gather all element which need to be 
#configured/added into one place, making this truly generic/extendable

sub filter_mismatches{
  
  my ($self, $query_match_length, $q_length, $score) = @_;

  my $mismatch;
  
  #score is + 5 for a match, -4 for mismatch: hence -9 from total score for one mismatch.
  my $full_score     = $q_length * 5;
  # Should add this as a config var?
  
  #my $mismatch_length = $q_length - 1;
  #my $mismatch_score = $mismatch_length * 5;
  
  #Can we not rewrite this to use the allowed mismatches conf var

  #We could just use the mismatch output instead of calculating all this nonsense!
  #Then we just build the soft cigar line from the info string and the match length query start!


  if($query_match_length == $q_length){
	
	if($score == $full_score){
	  $mismatch = 0;
	}
  }

  if(! defined $mismatch){
	
	my $max_mismatch = $self->max_mismatches;
	
	for my $i(1..$max_mismatch){
	  
	  my $mismatch_length = $q_length - $i;
	  my $mismatch_score = $mismatch_length * 5;
	  
	  if($query_match_length == $q_length){
		
		if ($score == ($full_score - ($i*9))) {
		  $mismatch = $i;
		}
	  }
	  elsif($query_match_length == $mismatch_length){
		$mismatch = $i if ($score == $mismatch_score);
	  }
	}
  }

  return $mismatch;
}



#This simply contains/returns the coderef which has been set to an internal method or sub specified in the RunnableDB config
#And passed to ExonerateProbe via the -filter_method param.

sub filter_method{
  my $self = shift;
  return $self->{'filter_method'};
}

#
# Create affy features attached to 'fake' affy array's and affy probe
# objects: these are identified by name only, and are then persisted by
# runnabledb: the runnabledb can sort out whether the attached probe/array
# is valid (exists in the db).
#
sub make_feature{
  my ($self, @args) = @_;
 
  #huh? Isn't this assign self to $tag?
 
  my (
    $tag, $q_id, $q_start, $q_end, $q_strand, 
    $t_id, $t_start, $t_end, $t_strand, $score, 
    $q_length, $query_match_length
  ) = @_;


  #because of the 'in-between' coordinates.???
  $t_start += 1;
  #What about t_end and q_start/end?


  my $probe_internal_id = $q_id;
  my $mismatch_count;
  
  if(!($probe_internal_id =~ /\d+/)){
    throw "Probe headers MUST be the internal db ids of the Probes for this parser to work!\n";
  }
  


  #we need to full unmatched query length to calculat the mismatches!
  #Can't assume will always be

  #=pod

  #If we miss by more than one, don't create the feature.
  #This is based on 25mers
  #we need to be able to adapt this for a % hit

  #This should be calculated based on % id?
  #score is + 5 for a match -4 for mismatch: hence -9 from total score for one mismatch.

  #if($query_match_length == $q_length){
 

    #if($score == 125){#25mers
#    if($score == 250){ #50mers
  #    $mismatch_count = 0;
    #}elsif ($score == 116) {
      #}elsif ($score == 241) {#50mer 98%
    #  $mismatch_count = 1;
      #}elsif ($score == 232) {#50mer 96%
      #  $mismatch_count = 2;  
    #} else {
 #     return undef;
    #}
  #}elsif($query_match_length == $q_length -1){
  #  if ($score == 120) {
      #if ($score == 245) {#50mer 100% over 49 bases
  #    $mismatch_count = 1;
      #}elsif ($score == 236) {#50mer just under 98% over 49 bases
      #  $mismatch_count = 2;  
  #  }
  #  else {
  #    return undef;
  #  }
  #}
    #elsif($query_match_length == $q_length -2){#50mers only!
    #  if($score == 240){
    #    $mismatch_count=2;
    #  }
    #}else{
    #  return undef;
    #}

  #=cut

  ##NATH HACK
    #100% ID for tiling arrays
    if(($query_match_length == $q_length) && $score == 250){
      $mismatch_count = 0;
    } else {
      return undef;
    }
  #98 % =  241 for 50 mers
  
  
  
  #Create a dummy probe object to pass out. It must be checked and
  #replaced by the parent runnabledb before being stored.

  #my $probe = 
  #  new Bio::EnsEMBL::Funcgen::Probe(
  #    -dbID =>$q_id,
  #    -name => "replace this probe name",
  #    -arrayname => "replace this array name",
  #  );
   
  # Everything is 'flipped' into the forward strand of the probe -
  # so a hit on the reverse strand of the probe (q_strand = '-1')
  # is altered:  q_strand = '+1' and t_strand => -1 x t_strand. 
  if($q_strand eq '+'){
    if($t_strand eq '+'){
      $t_strand = 1;
    }elsif($t_strand eq '-'){
      $t_strand = -1;
    }else{
      throw "unrecognised target strand symbol: $t_strand\n";
    }
  }elsif($q_strand eq '-'){
    if($t_strand eq '-'){
      $t_strand = 1;
    }elsif($t_strand eq '+'){
      $t_strand = -1;
    }else{
      throw "unrecognised target strand symbol: $t_strand\n";
    }
  }else{    
      throw "unrecognised query strand symbol: $q_strand\n";
  }
  

  my $feature =
    new Bio::EnsEMBL::Funcgen::ProbeFeature(
											#-probe => $probe,
											-probe_id => $q_id,
											-start => $t_start,
											-end => $t_end,
											-strand => $t_strand,
											-mismatchcount => $mismatch_count,
										   );

  # attach the slice name onto the feature: let the runnabledb
  # sort out whether it's valid.
  $feature->seqname($t_id);
  
  return $feature;
}


1;

