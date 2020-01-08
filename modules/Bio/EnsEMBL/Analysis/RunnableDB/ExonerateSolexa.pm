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

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa - 

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
  );

$runnableDB->fetch_input();
$runnableDB->run();
$runnableDB->write_output(); #writes to DB

=head1 DESCRIPTION

Extends Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature to allow
use warnings ;
use of compressed dna align features, useful when aligning millions of short 
Solexa reads

=head1 METHODS


=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::ExonerateSolexa;
use Bio::EnsEMBL::Analysis::Tools::Utilities;


use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature);


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->db->disconnect_when_inactive(1); 
  $self->read_and_check_config($EXONERATE_SOLEXA_CONFIG_BY_LOGIC);  
  $self->db->disconnect_when_inactive(1); 
  return $self;
}


sub run {
  my ($self) = @_;
 # do all the normal stuff 
 
  $self->db->disconnect_when_inactive(1);  
  $self->SUPER::run();  
  print "run finished 1 \n"; 
  $self->db->disconnect_when_inactive(1);  
  # filter the results carefully to only allow strict matches
  my $filtered_features = $self->filter_solexa($self->output);
  $self->db->disconnect_when_inactive(1);  
  # Pair features together if they come from paired end reads
  if ( $self->PAIREDEND ) {
    my $paired_features = $self->pair_features($filtered_features);
    # clear the output array to load in the modified features
    $self->{'output'} = [];
    $self->output($paired_features);
  }
  $self->db->disconnect_when_inactive(1);  
}

sub filter_solexa {
  my ($self,$features_ref) = @_;
  my @features = @$features_ref;
  my @filtered_features;
  # allow no more than MISSMATCH missmatches and
  # dont allow gapping...
  foreach my $feat ( @features ) {
    # dont allow gapping
    if ( $feat->cigar_string =~ /[ID]/ ) {
      # gapped reject it
      next
    }
   # check missmatches
    if ( $self->MISSMATCH ) {
      my $aligned_length = abs($feat->hend - $feat->hstart) +1;
     # print $feat->hseqname ." ";
     # print $feat->percent_id . " " ;
     # print " aligned length = $aligned_length ";
      my $matches = $aligned_length *  $feat->percent_id / 100;
     # print " matches $matches  ";
      my $missmatches = ( $aligned_length - $matches) / $aligned_length * 100;
     # print " missmatches $missmatches \n";
      next if $missmatches > $self->MISSMATCH;
     # print "accepting\n";
    }
   push @filtered_features, $feat;
  }
  return \@filtered_features;
}


sub pair_features {
  my ($self,$features_ref) = @_;
  my @features = @$features_ref;
  my $paired_features;
  my @selected_reads;
  foreach my $feat ( @features ) {
     
    if ( $feat->hseqname =~ /(\S+):((a|b|A|B|HA|HB))$/ ) {
      my $id = $1;
      my $suffix = lc($2);
      $suffix =~ s/^h//;
      push @{$paired_features->{$id}->{$suffix}},$feat;
    } else {
      $self->throw("Read name not recognised " . $feat->hseqname . "\n");
    }
  }

  foreach my $pair ( keys %$paired_features ) {
    my ( @as, @bs ) ;
    @as =  @{$paired_features->{$pair}->{'a'}} if $paired_features->{$pair}->{'a'};
    @bs =  @{$paired_features->{$pair}->{'b'}} if $paired_features->{$pair}->{'b'};
    my @potential_pairs;
    # Pick the potential paired reads bases on chr positon
    foreach my $afeat ( @as ) {
    READ:  foreach my $bfeat ( @bs ) {
	if ( $afeat->seqname eq $bfeat->seqname ) {
	  # They lie on the same chr within PAIREDGAP bases of each other
	  # allow them to overlap
	  if ( ( $afeat->end   > $bfeat->start && 
		 $afeat->start - $bfeat->end <= $self->PAIREDGAP ) or 
	       ( $bfeat->end   > $afeat->start && 
		 $bfeat->start - $afeat->end <= $self->PAIREDGAP ) or
	       # they overlap
	       ( $afeat->start   <= $bfeat->end && 
		 $afeat->end     >= $bfeat->start )) {
	    # they should be oriented pointing towards each other on
	    # opposite strands...
	    if ( $afeat->strand == 1 && $bfeat->strand == 1
	       && $afeat->hstrand != $bfeat->hstrand ) {
	      # are they in the right order
	      # on the forward strand a should come before b
	      # and vice versa
	      next READ unless 
		( 
		 ($afeat->hstrand == 1  && $afeat->start < $bfeat->start) or 
		 ($afeat->hstrand == -1 && $bfeat->start < $afeat->start) 
		);
	    }

	    # storing the pair using the score as the index
	    # so I can easily pull out the highest scoring pair
	    # if they score the SAME then choose the pair that is
	    # closest together 
	    my $combined_score = $afeat->score + $bfeat->score;
	    # work out the lengths the pairs are separated by
	    if ( $potential_pairs[$combined_score] ) {
	      my $new_length;
	      my $old_length;
	      my $old_afeat = $potential_pairs[$combined_score]->[0];
	      my $old_bfeat = $potential_pairs[$combined_score]->[1];
	      if ( $afeat->start > $bfeat->end ) {
		$new_length = $afeat->start - $bfeat->end;
	      } else {
		$new_length = $bfeat->start - $afeat->end;
	      }
	      if ( $old_afeat->start > $old_bfeat->end ) {
		$old_length = $old_afeat->start - $old_bfeat->end;
	      } else {
		$old_length = $old_bfeat->start - $old_afeat->end;
	      }
	      # only replace the pair if the new pair with the same combined score
	      # are closer together ( hopefully this might stop clusters getting
	      # tangled up
	      if ( $old_length > $new_length ) {
		$potential_pairs[$combined_score] = [($afeat,$bfeat)] ;
	      }
	    } else {
	      $potential_pairs[$combined_score] = [($afeat,$bfeat)] ;
	    }
	  }
	}
      }
    }
    if ( scalar(@potential_pairs) >= 1 ){
      # lets join together the highest scoring pair of reads to make one new feature
      push @selected_reads, @{$self->merge_pair(pop @potential_pairs)} if scalar(@potential_pairs) >= 1;
    } else {
      if ( scalar(@as) >= 1 ) {
	@as = sort { $a->score <=> $b->score } @as ;
	my $selected = pop @as;
	if ( $selected->hseqname =~ /(\S+):A$/ ) {
	  $selected->hseqname($1 .":a3p");
	}
	if ( $selected->hseqname =~ /(\S+):HA$/ ) {
	  $selected->hseqname($1 .":ha3p");
	}	
	push @selected_reads, $selected;
      }
      if ( scalar(@bs) >= 1 ) {    
	@bs = sort { $a->score <=> $b->score } @bs ;
	my $selected = pop @bs;
	if ( $selected->hseqname =~ /(\S+):B$/ ) {
	  $selected->hseqname($1 .":b3p");
	}
	if ( $selected->hseqname =~ /(\S+):HB$/ ) {
	  $selected->hseqname($1 .":hb3p");
	}
	push @selected_reads, $selected;
      }
    }
  }
  return \@selected_reads;
}

sub merge_pair {
  my ( $self, $read_pair) = @_;
  my @ugfs;
  my @features;
  if ( $self->PAIR ) {
    # Use the combined scores and identity  of the reads for the new read
    my $score = $read_pair->[0]->score +  $read_pair->[1]->score;
    my $pid = ( $read_pair->[0]->percent_id + $read_pair->[1]->percent_id ) /2;
  #  print "SELECTED "  . $read_pair->[0]->seq_region_name , " " .       $read_pair->[0]->start . " " .	$read_pair->[0]->end ." " .	  $read_pair->[0]->strand . " " .	    $read_pair->[0]->hseqname ."\n";
  #  print "SELECTED "  . $read_pair->[1]->seq_region_name , " " .       $read_pair->[1]->start . " " .	$read_pair->[1]->end ." " .	  $read_pair->[1]->strand . " " .	    $read_pair->[1]->hseqname ."\n";   
    # Because the features are split you can get the second half of
    # the read aligning before the first which can screw up the 
    # hit start ends as you might expect
    # I am  going to flip a and b over
    if ( $read_pair->[1]->end <= $read_pair->[0]->start ) {
      my $tmp = $read_pair->[0];
      $read_pair->[0] = $read_pair->[1];
      $read_pair->[1] = $tmp;
      
    }
    # sort out the hit names so they are consistant
    # we need a 3 prime read to be labeled on 
    # both reads or the merging code will
    # fall over with inconsistant hit names
    for ( my $i = 0 ; $i <=1 ; $i++ ) {
      my $read = $read_pair->[$i];
     if ( $read->hseqname =~ /(\S+):((a|b))$/ ) {
     	$read->hseqname($1);
     }
      if ( $read->hseqname =~ /(\S+):(A|B)$/ ) {
     	$read_pair->[0]->hseqname($1 .":3p");
     	$read_pair->[1]->hseqname($1 .":3p");
        $read = $read_pair->[$i];
     }    
      if ( $read->hseqname =~ /(\S+):(HA|HB)$/ ) {
     	$read_pair->[0]->hseqname($1 .":h3p");
     	$read_pair->[1]->hseqname($1 .":h3p");
        $read = $read_pair->[$i];
     }
   }
   
   for ( my $i = 0 ; $i <=1 ; $i++ ) {
      my $read = $read_pair->[$i];
     foreach my $ugf ( $read->ungapped_features ) {
	# need to modify the 'b' read hit starts ends to make it
	# appear to come from one long read
	# neeed to flip it onto the forward strand or else 
	# the dna align feature will fall over
	my $hstart = $ugf->hstart ;
	my $hend   = $ugf->hend ;
	if ( $i == 1 ) {
	  if ( $ugf->hstrand == -1 ) {
	    $ugf->hstrand(1);
	  }
	} else {
	  $ugf->hstrand(1) if $ugf->hstrand == -1;
	}
	$ugf->score($score);
	$ugf->percent_id($pid);
	push @ugfs, $ugf;
      }
    }
    # if they overlap just make 1 feature that combines the two reads
    # @ugfs = sort { $a->start <=> $b->start } @ugfs;
    # they overlap - merge them
    if ( scalar(@ugfs) == 2 ) {
      if (  $ugfs[0]->end   >= $ugfs[1]->start && 
	    $ugfs[0]->start <= $ugfs[1]->end ) {
	$ugfs[0]->start($ugfs[1]->start) if
	  $ugfs[1]->start < $ugfs[0]->start;
	$ugfs[0]->end($ugfs[1]->end) if
	  $ugfs[1]->end > $ugfs[0]->end;
	$ugfs[0]->hstart( 1);
	$ugfs[0]->hend( ($ugfs[0]->end -$ugfs[0]->start )+1);
	# remove the last read just let through the modified read
	pop(@ugfs);
      }
    } else {
      # dont try to merge complex pairs where there is a 
      # difficult cigar line
      $read_pair->[0]->hseqname($read_pair->[0]->hseqname . ":aa");
      $read_pair->[1]->hseqname($read_pair->[1]->hseqname . ":bb");
      push @features,($read_pair->[0],$read_pair->[1]);   
      return \@features;
    }
    
    my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@ugfs);
    $feat->analysis($self->analysis);
    push @features,$feat;
  } else {  
    $read_pair->[0]->hseqname($read_pair->[0]->hseqname . ":aa");
    $read_pair->[1]->hseqname($read_pair->[1]->hseqname . ":bb");
    push @features,($read_pair->[0],$read_pair->[1]);
  }
  # Otherwise leave them as they are but change the a and b suffix to indictate that they have
  # been paired but are not joined together.
  return \@features;
}

=head2 write_output

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Overrides the write_output method from the superclass
              to allow use of compressed DnaAlignFeatures
  Returntype: 1
  Exceptions: Throws if the feature cannot be stored

=cut

sub write_output {
  my ( $self, @output ) = @_;
  # Flag set to 1 = return a pipeline adaptor
  my $outdb; 

   
  my $fa;
  if ( $self->COMPRESSION ) {
    # Flag set to 1 = return a pipeline adaptor 
    # Remember you need to have the pipeline_tables in your OUT_DB if you like to use compression  
    $outdb = $self->get_dbadaptor($self->OUT_DB, 'pipeline');
    $outdb->disconnect_when_inactive(1);
    $fa = $outdb->get_CompressedDnaAlignFeatureAdaptor;
  } else {
    if ( defined $self->OUT_DBS ) { 
       
      # randomly pick an output db from an array of possible dbs
      my @dbs =  @{$self->OUT_DBS};
      my $number = int(rand(scalar(@dbs))) ;
      # dont let it fall off the end of the array
      $number-- if $number == scalar(@dbs) ;
      $outdb = $self->get_dbadaptor($dbs[$number]);
      print STDERR "Picking db  " . $outdb->dbc->dbname ." out of ". scalar(@dbs) . " possible output databases\n";

    } else {
      # return a NORMAL adaptor NOT a pipeline adaptor and dno not attach dna_db 
      $outdb = $self->get_dbadaptor($self->OUT_DB,undef,1);  
    } 
    $outdb->disconnect_when_inactive(1);
    $fa = $outdb->get_DnaAlignFeatureAdaptor;
  }

    print "writing output\n";  
    # calculate the hcoverage and use the evalue feild to store the
    # depth of the features
    
   if ( $self->COMPRESSION ) { 
     foreach my $f (@{$self->output}){
      $f->p_value(1) ;
      eval{
        $fa->store($f);
      };
      if ($@) {
        $self->throw("Unable to store DnaAlignFeature\n $@");
      }
     }
   } else { 
      print "no compression\n";  
      eval{
        $fa->store(@{$self->output});
      };
      if ($@) {
       $self->throw("Unable to store DnaAlignFeature\n $@");
     }
  }
  #$outdb->disconnect_if_idle();
}

###########################################################
# containers

sub INTRON_MODELS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_INTRON_MODELS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_INTRON_MODELS'})) {
    return $self->{'_CONFIG_INTRON_MODELS'};
  } else {
    return undef;
  }
}

sub INTRON_OVERLAP {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_INTRON_OVERLAP'} = $value;
  }
  
  if (exists($self->{'_CONFIG_INTRON_OVERLAP'})) {
    return $self->{'_CONFIG_INTRON_OVERLAP'};
  } else {
    return undef;
  }
}

sub TRANSCRIPT_BIOTYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_TRANS_BIOTYPE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_TRANS_BIOTYPE'})) {
    return $self->{'_CONFIG_TRANS_BIOTYPE'};
  } else {
    return undef;
  }
}

sub TRANSDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_TRANSDB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_TRANSDB'})) {
    return $self->{'_CONFIG_TRANSDB'};
  } else {
    return undef;
  }
}

sub OUT_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUT_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OUT_DB'})) {
    return $self->{'_CONFIG_OUT_DB'};
  } else {
    return undef;
  }
} 

sub COMPRESSION {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_COMPRESSION'} = $value;
  }
  
  if (exists($self->{'_CONFIG_COMPRESSION'})) {
    return $self->{'_CONFIG_COMPRESSION'};
  } else {
    return undef;
  }
}

sub PAIREDEND {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PAIREDEND'} = $value;
  }
  
  if (exists($self->{'_CONFIG_PAIREDEND'})) {
    return $self->{'_CONFIG_PAIREDEND'};
  } else {
    return undef;
  }
}

sub PAIR {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PAIR'} = $value;
  }
  
  if (exists($self->{'_CONFIG_PAIR'})) {
    return $self->{'_CONFIG_PAIR'};
  } else {
    return undef;
  }
}

sub PROJECT {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PROJECT'} = $value;
  }
  
  if (exists($self->{'_CONFIG_PROJECT'})) {
    return $self->{'_CONFIG_PROJECT'};
  } else {
    return undef;
  }
}

sub PAIREDGAP {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PAIREDGAP'} = $value;
  }
  
  if (exists($self->{'_CONFIG_PAIREDGAP'})) {
    return $self->{'_CONFIG_PAIREDGAP'};
  } else {
    return undef;
  }
}

sub MISSMATCH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MISSMATCH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MISSMATCH'})) {
    return $self->{'_CONFIG_MISSMATCH'};
  } else {
    return undef;
  }
}

sub OUT_DBS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUT_DBS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OUT_DBS'}) and join(@{$self->{'_CONFIG_OUT_DBS'}}) ne '') {
    return $self->{'_CONFIG_OUT_DBS'};
  } else {
    return undef;
  }
}

1;
