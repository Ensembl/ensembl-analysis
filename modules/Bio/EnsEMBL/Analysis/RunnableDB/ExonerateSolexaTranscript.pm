=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaTranscript - 

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

=head1 METHODS


=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexaTranscript;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::ExonerateSolexa;
use Bio::EnsEMBL::FeaturePair;
use Bio::SeqIO;
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
    print "MSG: fetching transcripts of biotype $biotype !\n"; 
    push @trans , @{$trans_adaptor->generic_fetch("biotype = \"$biotype\"")};
  } else {
    push @trans , @{$trans_adaptor->fetch_all(undef,undef,undef)};
  }
  print scalar(@trans) . " transcripts fetched \n";   

  foreach my $trans ( @trans ) {
    $trans_by_id{$trans->display_id} = $trans;
  }
  $self->transcripts(\%trans_by_id);

  # fetch the sequence of the query seqs if we want to write them as SAM format
  if ( defined $self->OUT_SAM_DIR ) {
    #open the file and read the sequences into a hash
    my $query_seqs = Bio::SeqIO->new( -file => $self->QUERYSEQS."/".$self->input_id,
				      -format => 'fasta'
				    );
    $self->throw("Cannot open " .  $self->QUERYSEQS ."/" .$self->input_id  . 
		 " for reading querys into hash \n") unless $query_seqs;
    while ( my $seq = $query_seqs->next_seq ) {
      my $id = $seq->display_id;
      $id =~ s/:A$/:a/;
      $id =~ s/:B$/:b/;
      $self->seq_hash($id,$seq);
    }
  }
  # close unused connections to transcript db
  $self->get_dbadaptor($self->TRANSDB)->disconnect_when_inactive(1);
  $self->get_dbadaptor($self->main_reference_db)->disconnect_when_inactive(1);
  $self->get_dbadaptor($self->main_reference_db,'pipeline')->disconnect_when_inactive(1);
  return ;
}

sub filter_solexa {
  my ($self,$features_ref) = @_;
  my @features = @$features_ref;
  my @filtered_features;
  # allow no more than MISSMATCH missmatches and
  foreach my $feat ( @features ) {
    # restrict just to splice models where read is within an exon
    if ( $self->INTRON_MODELS ) {
      next unless $feat->{"_intron"};
    }
    # check missmatches
    if ( $self->MISSMATCH ) {
      my $aligned_length = abs($feat->hend - $feat->hstart) +1;
      my $matches = $aligned_length *  $feat->percent_id / 100;
      my $missmatches = ( $aligned_length - $matches) / $aligned_length * 100;
      next if $missmatches > $self->MISSMATCH;      
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
  # force a recconect so we can get the dna and write output
  $self->get_dbadaptor($self->main_reference_db)->dbc->connect 
    unless  $self->get_dbadaptor($self->main_reference_db)->dbc->connected;
  $self->get_dbadaptor($self->main_reference_db,'pipeline')->dbc->connect
    unless  $self->get_dbadaptor($self->main_reference_db,'pipeline')->dbc->connected;
  $self->get_dbadaptor($self->TRANSDB)->disconnect_when_inactive(0);
  $self->get_dbadaptor($self->main_reference_db)->disconnect_when_inactive(0);
  if ($self->MISSMATCH or  $self->INTRON_MODELS ) {
    $features = $self->filter_solexa($features);
  }

  # Pair features together if they come from paired end reads
  # dont pair intron reads as the strands get screwed up
  if ( $self->PAIREDEND && !$self->INTRON_MODELS ) {
    $features = $self->pair_features($features);
  }

  my $genomic_features = $self->process_features($features);
  $self->output($genomic_features);
}

# overide write output if BAm files are to be used

=head2 write_output

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Overrides the write_output method from the superclass
              to allow Writing in SAM format
  Returntype: 1
  Exceptions: Throws if the feature cannot be stored

=cut

sub write_output {
  my ( $self ) = @_;
  my @output = @{$self->output};
  print "Got " .  scalar(@output) ." genomic features \n";
  unless (  $self->OUT_SAM_DIR ) {
    # write to the db
    $self->SUPER::write_output();
    return;
  } else {
    # write to file
    my $filename = $self->input_id;
    $filename =~ s/.fa$//;
    open ( SAM ,">".$self->OUT_SAM_DIR."/$filename.sam" ) or
      $self->throw("Cannot open file for writing " . $self->OUT_SAM_DIR."./$filename.sam\n");
    # write header
    my $refdb = $self->get_dbadaptor("REFERENCE_DB");
    my $sa = $refdb->get_SliceAdaptor;
    print SAM "\@HD\tVN:1.0\n";
    foreach my $slice ( @{$sa->fetch_all('toplevel')} ) {
      print SAM "\@SQ\tSN:" . $slice->seq_region_name ."\tLN:" . $slice->length ."\n"; 
    }
    foreach my $feature ( @output ) {
      my $line = $self->convert_to_sam($feature);
      print SAM $line if $line;
    }
  }
}

sub convert_to_sam {
  my ( $self,$feature ) = @_;
  my $line;
  # I have a feeling the hit strands are all backwards in the db
  # I will have a go at reversing them 
  #$feature->hstrand($feature->hstrand * -1);
  # strip off the :NC of the end of the seq names
  my $tmp = $feature->hseqname;
  $tmp   =~ s/:NC$//;
  $feature->hseqname($tmp);
  my $seq = "*";
  my $type;
  my $flag = 1;
  my $paired = 1;
  # reverses cigar line if feature is reversed
  $feature = $self->check_cigar($feature);
  my $cigar = $feature->cigar_string;
  # N isued to represent splicing in sam format
  $cigar =~ s/I/N/g;
  # is the feature paired in the alignment
  if ( $feature->hseqname =~ /(\S+):(a:aa|b:bb|A:aa|B:bb)/ ) {
    # paired but not joined
    $paired = 0;	 
    $feature->hseqname($1);
    $flag +=2 ;
    if ( $2 =~ /a/ ) {
      $flag += 64;
      $type = "a";
    } elsif ( $2 =~ /b/ ) {
      $type = "b";
      $flag += 128;
    }
  } elsif ( $feature->hseqname =~ /(\S+):(a|b|A|B|a3p|b3p)/ ) { 
    $paired = 0;
    $feature->hseqname($1);
    if ( $2 =~ /(a|A)/ ) {
      $flag += 64;
      $type = "a";
    } elsif ( $2 =~ /(b|B)/ ) {
      $type = "b";
      $flag += 128;
    }
  } else {
    $flag =3 ;
  }
  # strand wrt the reference
  # do we know the strand?
  my $strand = $feature->strand ;
  if ( $strand == -1 ) {
    $flag +=16;
  }
  # strand of the mate = -1 so we add nothing
  my $feature_seq = $self->seq_hash($feature->hseqname.":$type");
#  print "SEQ " . $feature->hseqname . " " . $feature_seq->seq ,"\n";
  $self->throw("cannot find sequence for " . $feature->hseqname .":$type\n")
    unless $feature_seq;
  my $length = length($feature_seq->seq);
  # add soft clip info to the cigar
  # cigar reversed if the hit strand is -1;
  if ( ( $feature->hstrand == -1 &&  $feature->strand == -1) or ( $feature->hstrand == 1 &&  $feature->strand == 1 ) ) {
    $cigar = $cigar . ($feature->hstart - 1) ."S"	if  $feature->hstart > 1;
    $cigar = ( $length - $feature->hend ) ."S$cigar"	if  $feature->hend < $length;
  } else {
    $cigar = ( $feature->hstart -1 ) ."S$cigar"	if  $feature->hstart > 1;
    $cigar = $cigar . ( $length - $feature->hend ) ."S"	if  $feature->hend < $length;
  }
  my $check = $self->check_cigar_length($cigar);
  if ( $check != $length ) {
    print STDERR "Losing " .  $feature->hseqname . " ($type) I get cigar length $check rather than $length, $cigar\n";
    return ;
  }
  if ( $feature->hstrand == 1 &&  $feature->strand == 1 ) {
    # need the reverse compliment
    $seq = $feature_seq->revcom->seq;
  } elsif( $feature->strand == -1 &&  $feature->hstrand == -1 )   {
    # need the reverse compliment
    $seq = $feature_seq->revcom->seq;   
  }else {
    $seq = $feature_seq->seq;
  }
  # if the feature is aligned to the -ve strand we want the complement of he sequence as well
  if ( $feature->strand == -1 ) {
    $seq =~ tr/ACGTacgt/TGCAtgca/;
  }

  # strand of the mate = -1 so we add nothing
  $line .= $feature->hseqname ."\t"; 
  $line .= "$flag\t";
 # $line .= $feature->hstrand . "\t" . $feature->strand . "\t";
  $line .= $feature->seq_region_name ."\t";
  $line .=  $feature->start ."\t";
  $line .= "0\t";
  $line .=  "$cigar\t";
  if ( $paired ) {
    # what is the position of the other pair?
    $line .= $feature->seq_region_name ."\t" . $feature->hstart ."\t" . $feature->hend . "\t$seq\t*\n";
  } else {
    $line .= "*\t0\t0\t$seq\t*\n";
  }
  return $line;
}

sub check_cigar {
  my ($self,$feat) = @_;
  my @ugfs;
  my $string = $feat->cigar_string;
  my @pieces = ( $string =~ /(\d*[MDI])/g );
  # if  strand is -1 strand cigar line is reversed
  @pieces = reverse @pieces if     $feat->strand == -1 ;
  my $cigar;
  foreach my $piece (@pieces) {
    my ($length,$match) = ( $piece =~ /^(\d*)(M|I|D)/ );
    if( $length eq "" ) { $length = 1 }
    $cigar.= $length;
    $cigar.= $match;
  }
  $feat->cigar_string($cigar);
  return $feat;
}

sub check_cigar_length{
  my ($self,$cigar) = @_;
  my @pieces = ( $cigar =~ /(\d*[SMDN])/g );
  my $total_length;
  foreach my $piece (@pieces) {
    my ($length,$match) = ( $piece =~ /^(\d*)(S|M|D|N)/ );
    if( $length eq "" ) { $length = 1 };
    $total_length += $length unless $match eq "N";
  }
  return $total_length;
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
    print "FEATURE " . $f->hseqname ."\n";
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

    if ( $self->INTRON_OVERLAP && !$self->INTRON_MODELS ) {
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
    }

    if ( $self->INTRON_MODELS ) {
      #we have already filtered everything else out
      $accept = 1;
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
      # if so tag it so we can tell later on (only if you are writing to db)
      # if we are writing to SAM we cannot store the info in the same way
      if ( $f->{'_intron'} and not $self->OUT_SAM_DIR) {
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
	    # is it cannonical
	    if ( $features[0]->strand  == 1 ) {
	      print "Splice type " . $left_splice->seq ."- ".  $right_splice->seq ." ";
	      # is it GTAG?
	      unless ( $left_splice->seq eq 'GT' && $right_splice->seq eq 'AG' ) {
		$cannonical = 0;
	      }
	    } else {
	      print "Splice type " . $right_splice->seq ."- ".  $left_splice->seq ." ";
	      # is it GTAG?
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
	# mark it as non cannonical splice
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

sub seq_hash {
  my ($self,$key,$value) = @_;

  return undef unless defined ($key);

  if (defined $value) {
    $self->{'_seq_hash'}{$key} = $value;
  }

  if (exists($self->{'_seq_hash'}{$key})) {
    return $self->{'_seq_hash'}{$key};
  } else {
    return undef;
  }
}
1;
