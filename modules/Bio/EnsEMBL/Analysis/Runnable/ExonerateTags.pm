=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::ExonerateTags - 

=head1 SYNOPSIS

 my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::ExonerateTags->new(
     -QUERYFILE     => $ditag_file,
     -TARGETFILE    => $genome,
     -TAGTYPE       => $tagtype,
    );
 $runnable->run();
 my @results = @{ $runnable->output() };

=head1 DESCRIPTION

This module handles a specific use of the Exonerate program (G. Slater),
to align CAGE, GIS/GSC ditags to genomic sequences.

=head1 TODO

When removing duplicates from the input sequences, this should be
added to the tag_count of the remaining ditags.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateTags;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Map::Ditag;
use Bio::EnsEMBL::Map::DitagFeature;
use Bio::EnsEMBL::Analysis::Config::ExonerateTags;
use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseExonerate);


=head2 new

  Args       : various
  Description: Runnable constructor
  Returntype : Bio::EnsEMBL::Analysis::Runnable::ExonerateTags
  Caller     : general

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $query_seqs, $query_file, $program, $analysis, $maxmismatch, $options, $maxdistance,
       $mindistance, $delete_queryfile, $specificoptions, $splitseqs, $maxhitsallowed,
       $keep_order) = rearrange(
       [ 'QUERYSEQS', 'QUERY_FILE' ,'PROGRAM', 'ANALYSIS', 'MAXMISMATCH', 'OPTIONS',
         'MAXDISTANCE', 'MINDISTANCE', 'DELETEQUERYFILE', 'SPECOPTIONS', 'SPLITSEQS',
	 'MAXHITS', 'KEEPORDER' ],
       @args );

  $self->query_seqs($query_seqs)            if $query_seqs;
  $self->query_file($query_file)            if $query_file;
  $self->options($options)                  if $options;
  $self->files_to_delete($self->query_file) if $delete_queryfile;
  $self->specificoptions($specificoptions)  if $specificoptions;
  $self->splitseqs($splitseqs)              if($splitseqs);
  $self->keep_order($keep_order)            if($keep_order);

  if( $program ) {
    $self->program($program);
  }
  else{
    $self->program('/usr/local/ensembl/bin/exonerate-1.0.0');
  }
  if( $analysis ) {
    $self->analysis($analysis);
  }
  else{
    $self->analysis('ExonerateTags');
  }
  if( $maxdistance ) {
    $self->maxdistance($maxdistance);
  }
  else{
    $self->maxdistance('600000');
  }
  if( $mindistance ) {
    $self->mindistance($mindistance);
  }
  else{
    $self->mindistance('100');
  }
  if( $maxmismatch ) {
    $self->maxmismatch($maxmismatch);
  }
  else {
    $self->maxmismatch('1');
  }
  if( $maxhitsallowed ) {
    $self->maxhitsallowed($maxhitsallowed);
  }
  else {
    $self->maxhitsallowed('4');
  }

  $self->{'_verbose'} = 0;

  return $self;
}


=head2 run

  Args       : various
  Description: Runnable runner
  Returntype : none (results stored as $self->output)
  Caller     : general

=cut

sub run {
  my ($self) = @_;
  my @output = ();
  my $result_features;

  if ( $self->query_seqs ) {
    # Write query sequences to file if necessary
    my $query_file = $self->workdir . "/ditagexonerate_q.$$";
    my $seqout = Bio::SeqIO->new( '-format' => 'fasta',
                                  '-file'   => ">$query_file" );
    foreach my $seq ( @{ $self->query_seqs } ) {
      $seqout->write_seq($seq);
    }
    $self->query_file($query_file);
    # register the file for deletion
    $self->files_to_delete($self->query_file);
  }
  print "query file: ".$self->query_file()."\n" if($self->_verbose);

  if ($self->target_seqs) {
    # Write query sequences to file if necessary
    my $target_file = $self->workdir . "/ditagexonerate_t.$$";
    my $seqout = Bio::SeqIO->new(
                                 '-format'   => 'fasta',
                                 '-file'     => ">$target_file"
                                 );
    foreach my $seq ( @{$self->target_seqs} ) {
      $seqout->write_seq($seq);
    }
    $self->target_file($target_file);
    $self->files_to_delete($target_file);
  }
  print "target file(s): ".$self->target_file()."\n" if($self->_verbose);

  # Build exonerate command
  my $command =
    $self->program . " " . $self->options     .
       " --querytype "   . $self->query_type  .
       " --targettype "  . $self->target_type .
       " --query "       . $self->query_file  .
       " --target "      . $self->target_file .
       " " . $self->specificoptions;
  print "Exonerate command : $command\n";

  # Execute exonerate command
  my $exo_fh;
  open( $exo_fh, "$command |" ) or throw("Error opening exonerate command: $? $!");

  #parse n save the output
  my $tmp_output = $self->parse_results( $exo_fh );

  close( $exo_fh ) or throw ("Error closing exonerate command: $? $!");

  #clean up
  $self->delete_files;

  #filter the hits
  if($self->splitseqs()) {
    $result_features = $self->filter_split_hits($tmp_output);
  }
  else {
    $result_features = $self->filter_nonsplit_hits($tmp_output);
  }

  #create proper ditagFeature objects
  my $output = $self->make_features($result_features);

  $self->output($output);
}


=head2 run_previous

  Args       : Exonerate output file
  Description: For testing and for ever-failing jobs, Tags can be aligned manually with exonerate
               and the outpu file fed back into the system with this function.
  Returntype : none, stores into output hash.

=cut

sub run_previous {
  my ($self, $exo_fh) = @_;
  my ($result_features, $tmp_output, $output);

  $self->splitseqs(1);
  $tmp_output        = $self->parse_results( $exo_fh );
  if($self->splitseqs()) {
    $result_features = $self->filter_split_hits($tmp_output);
  }
  else {
    $result_features = $self->filter_nonsplit_hits($tmp_output);
  }
  $output            = $self->make_features($result_features);

  $self->output($output);
}


=head2 parse_results

  Args       : open file handle with filtered exonerate output lines
  Description: Create temporary feature objects to store exonerate hits
  Returntype : Hash ref with temporary feature objects sorted into left and right.

=cut

sub parse_results {
  my ( $self, $fh ) = @_;

  my %alltags = ();
  my $c = 0;
  my ($q_name, $q_side, $position);

  while (my $hits = <$fh>) {

    #parse the result lines
    next unless $hits =~ /^RESULT:/;
    chomp $hits;
    my (
        $null,     $q_id,             $q_start,    $q_end,
        $q_strand, $t_id,             $t_start,    $t_end,
        $t_strand, $score,            $perc_id,    $q_length,
        $t_length, $gene_orientation, @vulgar_blocks
       ) = split(' ', $hits);

    #get cigar string
    my $cigar_string        = '';
    my $match_type          = '';
    my $query_match_length  = '';
    my $target_match_length = '';
    while (@vulgar_blocks){
      throw("Something funny has happened to the input vulgar string." .
	    " Expecting components in multiples of three, but only have [" .
	    scalar @vulgar_blocks . "] items left to process.")
	unless scalar @vulgar_blocks >= 3;
      $match_type          = shift @vulgar_blocks;
      $query_match_length  = shift @vulgar_blocks;
      $target_match_length = shift @vulgar_blocks;
      if ($match_type eq "G"){
	if ($query_match_length == 0){
	    $match_type = "D";
            $query_match_length = $target_match_length;
        }elsif ($target_match_length == 0){
            $match_type = "I";
	}
      }
      $cigar_string .= $query_match_length.$match_type;
    }

    #decipher the header line
    if($self->splitseqs()){

      $q_id =~ /^(\S+)_([RL])_([0-9]+)/;
      $q_name   = $1;
      $q_side   = $2;
      $position = $3;

      if(!exists($alltags{$q_name})) {
        $alltags{$q_name}{'L'} = [];
        $alltags{$q_name}{'R'} = [];
      }
    }
    else{
      $q_id =~ /^(\S+)_F/;
      $q_name   = $1;
      $q_side   = "F";
      $position = 0;

      if(!exists($alltags{$q_name})){
        $alltags{$q_name}{'F'} = [];
      }
    }
    if(!$q_name){
      warning "\nproblem with header of $hits.\n";
      next;
    }

    #store a temporary (extended ditag) object
    my %tempobject = (
                      'slice_name'    => $t_id,
                      'start'         => $t_start,
                      'end'           => $t_end,
                      'strand'        => $q_strand,
                      'hit_start'     => $q_start,
                      'hit_end'       => $q_end,
                      'hit_strand'    => $t_strand,
                      'perc'          => $perc_id,
                      'length'        => $q_length,
                      'score'         => $score,
                      'ditag_id'      => $q_name,
                      'ditag_side'    => $q_side,
                      'cigar_line'    => $cigar_string,
                      'position'      => $position,
                      'ali_length'    => $query_match_length,
		      'ditag_pair_id' => 1,
                     );

    #store by its name/id and side: $alltags{$q_name}{$q_side}->@taghits
    push(@{$alltags{$q_name}{$q_side}}, \%tempobject);
    $c++;
  }
  print "\ntotal hits:\t$c\n";

  return \%alltags;
}


=head2 filter_split_hits

  Args       : ref to hash with paired hash with result lines from exonerate run (L&R)
  Description: filter the Exonerate hits in pairs. We only want tags, that
               are pairs in the same genomic region and that have almost perfect identidy.
  Returntype : Arrayref with temporary feature objects in pairs

=cut

sub filter_split_hits {
  my ( $self, $alltags ) = @_;
  my (@tokeep, @tokeep_1, @tokeep_2);
  my $strnd      = 0;
  my $chrom      = 0;
  my $dist       = 0;
  my $unspec     = 0;
  my $scorer     = 0;
  my $maxlength  = 0;
  my $mismatch   = 0;
  my $topscore   = 0;
  my $toomany    = 0;
  my $nopair     = 0;
  my $initial    = 0;
  my %combscores = ();
  my $order      = 0;
  my ($start_1, $start_2);

  #object: $alltags{$q_name}{$q_side}->@taghits

 ALLTAGS:
  foreach my $tagsides (keys %$alltags){
    my @taglist_L = @{$$alltags{$tagsides}{"L"}};
    my @taglist_R = @{$$alltags{$tagsides}{"R"}};

    $initial += (scalar @taglist_L) + (scalar @taglist_R);
    print "Working on split hits of ditag $tagsides.\n".
          "Initial hits: $initial\n" if($self->{'_verbose'});

    #throw away ditags with too many hits
    if( (scalar @taglist_L > 100) or (scalar @taglist_R > 100) ){
      print "Not storing tags: too many hits (".(scalar @taglist_L)."/".
	    (scalar @taglist_R).").\n" if($self->{'_verbose'});
      $toomany++;
      next ALLTAGS;
    }

    #throw away ditags without hits for both tags
    if( (!scalar @taglist_L) or (!scalar @taglist_R) ){
      print "Not storing tags: one match only.\n" if($self->{'_verbose'});
      $nopair += (scalar @taglist_L) + (scalar @taglist_R);
      next ALLTAGS;
    }

    #compare hits
   OUTER:
    for(my $i=0; $i<(scalar @taglist_L); $i++){

      #discard hits with high mismatches
      if( $taglist_L[$i]->{'ali_length'} < ($taglist_L[$i]->{'length'} - $self->maxmismatch()) ){
        print "high mismatch\n" if($self->{'_verbose'});
        $mismatch++;
        next OUTER;
      }

     INNER:
      for(my $j=0; $j<(scalar @taglist_R); $j++){

        #discard hits with high mismatches
        if($taglist_R[$j]->{'ali_length'} < ($taglist_R[$j]->{'length'} - $self->maxmismatch()) ){
          print "high mismatch\n" if($self->{'_verbose'});
          $mismatch++;
          next INNER;
        }

        #check chromosome
        if($taglist_L[$i]->{'slice_name'} ne $taglist_R[$j]->{'slice_name'}){
          print "chromosomes don t match.\n" if($self->{'_verbose'});
          $chrom++;
          next INNER;
        }

        #check strands
        if($taglist_L[$i]->{'strand'} ne $taglist_R[$j]->{'strand'}){
          print "strands don t match.\n" if($self->{'_verbose'});
          $strnd++;
          next INNER;
        }

        #look at locations
        $start_1 = $taglist_L[$i]->{'start'};
        $start_2 = $taglist_R[$j]->{'start'};

	#if applicable discard pairs where order of start/end is wrong
	if( ($self->keep_order) and
	    ($taglist_L[$i]->{'strand'} eq "+" and ($start_1 > $start_2)) or
	    ($taglist_L[$i]->{'strand'} eq "-" and ($start_1 < $start_2)) ) {
		print "start_end_order is wrong.\n" if($self->{'_verbose'});
		$order++;
		next INNER;
	}

	#check distance: too long
        if(abs($start_1 - $start_2) > $self->maxdistance){
          print "distance too large: ".(abs($start_1 - $start_2)).".\n"
	    if($self->{'_verbose'});
          $dist++;
          next INNER;
        }
	#check distance: too short
        if(abs($start_1 - $start_2) < $self->mindistance){
          print "distance too small: ".(abs($start_1 - $start_2)).".\n"
	    if($self->{'_verbose'});
          $dist++;
          next INNER;
        }

        #save this pair
        print "keep ($i/$j): ".$taglist_L[$i]->{'ditag_id'}." / ".
               $taglist_R[$j]->{'ditag_id'}."\n" if($self->{'_verbose'});
        push(@tokeep_1, $taglist_L[$i]);
        push(@tokeep_2, $taglist_R[$j]);
      }
    }

    #only 1-2 pairs should be left
    #check quality, keep only the ones with the best score of both tags
    %combscores = ();
    $topscore = 0;
    if(scalar @tokeep_1 > 1){
      #keep the top combined scorers
      for(my $i=0; $i<(scalar @tokeep_1); $i++){
        my $combscore = $tokeep_1[$i]->{'score'} + $tokeep_2[$i]->{'score'};
        if($topscore < $combscore){
          $topscore = $combscore;
          $combscores{$topscore} = [ $i ];
        }
        elsif($topscore == $combscore){
          push(@{ $combscores{$topscore} }, "$i");
        }
      }
      my (@bestscore_1, @bestscore_2);
      my $index = 1;
      for(my $k=0; $k<(scalar @{$combscores{$topscore}}); $k++){
	#clone the tag objects, add an id for the pair
	my %savethis = %{$tokeep_1[$k]};
	my %savethat = %{$tokeep_2[$k]};
	$savethis{'ditag_pair_id'} = $index;
        $savethat{'ditag_pair_id'} = $index;
	$index++;
        push(@bestscore_1, \%savethis);
        push(@bestscore_2, \%savethat);
      }
      @tokeep_1 = @bestscore_1;
      @tokeep_2 = @bestscore_2;
    }

    #if still many identical hits: throw all away because of unspecificity.
    if(scalar @tokeep_1 > $self->maxhitsallowed){
      @tokeep_1 = ();
      @tokeep_2 = ();
      print "removed unspecific hits.\n" if($self->{'_verbose'});

    }
    else {
      #join both resultlists
      push(@tokeep, @tokeep_1);
      push(@tokeep, @tokeep_2);
      @tokeep_1 = ();
      @tokeep_2 = ();
    }

  }

  print "\n\nSTATS\tInitial_tags     : $initial".
        "\nSTATS\ttoo_many         : $toomany".
	"\nSTATS\tno_pair          : $nopair".
	"\nSTATS\tKeeping_tags     : ".(scalar @tokeep)."\n";

  return \@tokeep;
}


=head2 filter_nonsplit_hits

  Args       : ref to hash with paired hash with result lines from exonerate run
  Description: filter the single Exonerate hits. We only want tags, that
               that have almost perfect identidy.
  Returntype : Arrayref with temporary feature objects

=cut

sub filter_nonsplit_hits {
  my ( $self, $alltags ) = @_;
  my (@tokeep, @tokeep_1);
  my $unspec     = 0;
  my $maxscore   = 0;
  my $maxlength  = 0;
  my $mismatch   = 0;
  my $topscore   = 0;
  my %combscores = ();

  #object: $alltags{$q_name}{"F"}->@taghits

 ALLTAGS:
  foreach my $tagsides (keys %$alltags){
    my @taglist_F = @{$$alltags{$tagsides}{"F"}};

    print "Working on non-split hits of ditag $tagsides.\n".
          " have ".(scalar @taglist_F)."\n" if($self->{'_verbose'});

    #throw way ditags with too many hits
    if(scalar @taglist_F > 100){
      print "Not storing tags: too many hits (".(scalar @taglist_F).").\n"
	if($self->{'_verbose'});
      next ALLTAGS;
    }

    #compare hits
    for(my $j=0; $j<(scalar @taglist_F); $j++){

      #discard hits with high mismatches
      if($taglist_F[$j]->{'ali_length'} < ($taglist_F[$j]->{'length'} - $self->maxmismatch()) ){
        print "high mismatch\n" if($self->{'_verbose'});
        $mismatch++;
        next;
      }

      #save this hit
      print "keep ($j): ".$taglist_F[$j]->{'ditag_id'}."\n" if($self->{'_verbose'});
      push(@tokeep_1, $taglist_F[$j]);
      #kepp track of max score value
      if($taglist_F[$j]->{'score'} > $maxscore){
	$maxscore = $taglist_F[$j]->{'score'};
      }

    }

    #only ~4 should be left
    #check quality, keep only the ones with the best score
    if(scalar @tokeep_1 > 1){
      #sort hits by score
      @tokeep_1 = sort {$b->{'score'} <=> $a->{'score'}} @tokeep_1;
      #keep the top scorers only
      $combscores{$maxscore} = [ ];
      for(my $i=0; $i<(scalar @tokeep_1); $i++){
        if($maxscore == $tokeep_1[$i]->{'score'}){
          push(@{ $combscores{$maxscore} }, $i);
        }
        else{
          last;
        }
      }
      my @bestscore_1;
      for(my $k=0; $k<scalar @{$combscores{$maxscore}}; $k++){
	#set pair-id to "0" as there is no pair
	$tokeep_1[$k]->{'ditag_pair_id'} = -1;
        push(@bestscore_1, $tokeep_1[$k]);
      }
      @tokeep_1 = @bestscore_1;
    }

    #if still many identical hits: throw all away because of unspecificity.
    if(scalar @tokeep_1 > $self->maxhitsallowed){
      $unspec += (scalar @tokeep_1);
      print "removed unspecific hits.\n" if($self->{'_verbose'});

    }
    else {
      push(@tokeep, @tokeep_1);
    }
    @tokeep_1 = ();

  }

  print "\n\nSTATS\tKeeping tags    : ".(scalar @tokeep).
        "\nSTATS\tlow identity    : $mismatch".
        "\nSTATS\ttoo_many_hits   : $unspec\n";

  return \@tokeep;
}


=head2 make_features

  Args       : array ref with temporary feature objects
  Description: Create DitagFeature objects for good exonerate hits
  Returntype : Bio::EnsEMBL::Map:DitagFeature

=cut

sub make_features {
  my ( $self, $tempfeatures ) = @_;
  my @ditagFeatures = ();

  foreach my $tempfeature (@$tempfeatures) {

    #set strand symbols
    my $q_strand = $tempfeature->{'strand'};
    my $t_strand = $tempfeature->{'hit_strand'};
    if    ( $q_strand eq '+' ) { $q_strand = 1 }
    elsif ( $q_strand eq '-' ) { $q_strand = -1 }
    else { throw "unrecognised query strand symbol: $q_strand\n" }
    if    ( $t_strand eq '+' ) { $t_strand = 1 }
    elsif ( $t_strand eq '-' ) { $t_strand = -1 }
    else { throw "unrecognised target strand symbol: $t_strand\n" }

    #work out real position
    my $real_tstart = $tempfeature->{'start'};
    my $real_tend   = $tempfeature->{'end'};
    my $real_qstart = $tempfeature->{'hit_start'} + $tempfeature->{'position'};
    my $real_qend   = $tempfeature->{'hit_end'}   + $tempfeature->{'position'};

    # for reverse strand matches, Exonerate reports end => start 
    if ($real_tstart > $real_tend) {
      ($real_tstart, $real_tend) = ($real_tend, $real_tstart);
    }
    if ($real_qstart > $real_qend) {
      ($real_qstart, $real_qend) = ($real_qend, $real_qstart);
    }
    #start at 1 not 0
    $real_tstart++;
    $real_qstart++;

    my $feature = Bio::EnsEMBL::Map::DitagFeature->new(
                                                       -slice         => $tempfeature->{'slice_name'},
                                                       -start         => $real_tstart,
                                                       -end           => $real_tend,
                                                       -strand        => $q_strand, 
                                                       -hit_start     => $real_qstart,
                                                       -hit_end       => $real_qend,
                                                       -hit_strand    => $t_strand,
                                                       -ditag_id      => $tempfeature->{'ditag_id'},
						       -ditag_pair_id => $tempfeature->{'ditag_pair_id'},
                                                       -ditag_side    => $tempfeature->{'ditag_side'},
                                                       -cigar_line    => $tempfeature->{'cigar_line'},
                                                       -analysis      => $self->analysis,
                                                      );
    push(@ditagFeatures, $feature);
  }
  print "\nhave Ditags:".(scalar @ditagFeatures)."\n";

  return \@ditagFeatures;
}


=head2 maxdistance

 Args       :  (optional) maximum distance allowed between the two end of a ditag
 Description:  Getter/Setter for maxdistance
 Returntype :  int

=cut

sub maxdistance {
  my $self = shift;
  $self->{'maxdistance'} = shift if (@_);
  return $self->{'maxdistance'};
}


=head2 mindistance

 Args       :  (optional) minimum distance allowed between the two end of a ditag
 Description:  Getter/Setter for mindistance
 Returntype :  int

=cut

sub mindistance {
  my $self = shift;
  $self->{'mindistance'} = shift if (@_);
  return $self->{'mindistance'};
}

=head2 maxmismatch

 Args       :  (optional) maximum number of mismatches allowed for exonerate hit
 Description:  Getter/Setter for maxmismatch
 Returntype :  int

=cut

sub maxmismatch {
  my $self = shift;
  $self->{'maxmismatch'} = shift if (@_);
  return $self->{'maxmismatch'};
}


=head2 maxhitsallowed

 Args       :  (optional) maximum number of hits of equally good quality
               before all are filtered out
 Description:  Getter/Setter for maxhitsallowed
 Returntype :  int

=cut

sub maxhitsallowed {
  my $self = shift;
  $self->{'maxhitsallowed'} = shift if (@_);
  return $self->{'maxhitsallowed'};
}


=head2 specificoptions

 Args       :  (optional) further options for exonerate
 Description:  Getter/Setter for specific exonerate options
 Returntype :  int

=cut

sub specificoptions {
  my $self = shift;
  $self->{'specificoptions'} = shift if (@_);
  return $self->{'specificoptions'};
}


=head2 splitseqs

 Args       :  (optional) 1 or 0
 Description:  Getter/Setter for split sequence indicator
 Returntype :  int

=cut

sub splitseqs {
  my $self = shift;
  $self->{'splitseqs'} = shift if (@_);
  return $self->{'splitseqs'};
}

=head2 repeatnumber

 Args       :  (optional) 1 or 0
 Description:  Getter/Setter for number of allowed repeated bases
 Returntype :  int

=cut

sub repeatnumber {
  my $self = shift;
  $self->{'repeatnumber'} = shift if (@_);
  return $self->{'repeatnumber'};
}

=head2 keep_order

 Args       :  (optional) 1 or 0
 Description:  Getter/Setter for option to keep the order of start/end tag when checking alignments
 Returntype :  true or false

=cut

sub keep_order {
  my $self = shift;
  $self->{'keep_order'} = shift if (@_);
  return $self->{'keep_order'};
}

1;
