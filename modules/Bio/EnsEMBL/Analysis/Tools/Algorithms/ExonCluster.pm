=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut


package Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster;


use warnings ;
use Bio::EnsEMBL::Transcript ; 
use vars qw(@ISA);
use Carp;
use strict;
use Bio::EnsEMBL::Utils::Exception qw (warning throw) ; 

sub new {
  my ($class,@args) = @_;

  my $self = bless {},ref($class) || $class;

  $self->{_start}    = undef;
  $self->{_end}      = undef;
  $self->{_exonhash} = undef;
  $self->{_exonidhash} = undef;
  $self->{_transcripthash} = undef;
  $self->{_transcriptidhash} = undef;
  $self->{_internal_index} = 0;
  $self->{_exon_2_biotype} = undef ;   # holds biotype of exon
  $self->{_all_exons_in_cluster} = {} ; # holds all exons ( key is exon itself ) 

  if (@args) {
    $self->throw("Constructor does not expect any parameters");
  }

  return $self;
}

sub start{
  my ($self,$start) = @_;
  if ($start){
    $self->throw( "$start is not an integer") unless $start =~/^[-+]?\d+$/;
    $self->{'_start'} = $start;
  }
  return $self->{'_start'};
}

sub end{
  my ($self,$end) = @_;
  if ($end){
    $self->throw( "$end is not an integer") unless $end =~/^[-+]?\d+$/;
    $self->{'_end'} = $end;
  }
  return $self->{'_end'};
}

sub length{
  my $self = shift @_;
  if (@_){
    $self->confess( ref($self)."->length() is read-only");
  }
  return ( $self->{'_end'} - $self->{'_start'} + 1 );
}

sub strand{
  my ($self,$strand) = @_;
  if ($strand){
    $self->{'_strand'} = $strand;
  }
  return $self->{'_strand'};
}

sub type {
  my ($self,$ignore_strand, @args) = @_;
  if (@args){
    $self->throw("type is a get only method");
  }
  if (!defined($self->{'_type'})) {
    $self->_determine_type($ignore_strand);
  }
  return $self->{'_type'};
}

sub merge {
  my ($self,$cluster, $ignore_strand) = @_;

  my %transhash =  %{ $cluster->each_transcripts_exons } ; 
  foreach my $transref (keys %transhash) {
    foreach my $exon (@{$transhash{$transref}}) {
      $self->add_exon($exon,$cluster->transcript_from_ref($transref), $ignore_strand);
    }
  }
  $self->{'_type'} = undef;
}

sub merge_new_exon {
  my ($self,$cluster, $ignore_strand) = @_;

  my %transhash =  %{ $cluster->each_transcripts_exons };  
  foreach my $transref (keys %transhash) {
    foreach my $exon (@{$transhash{$transref}}) {
      $self->add_exon_if_not_present($exon,$cluster->transcript_from_ref($transref), $ignore_strand);
    }
  }
  $self->{'_type'} = undef;
}

sub add_exon {
  my ($self,$exon,$transcript, $ignore_strand) = @_;

  if (!$self->contains_exon($exon)) {
    $self->_add_new_exon($exon, $ignore_strand);
    # print "Added exon " . $exon->dbID . "\n";
  }
  $self->_add_transcript_reference($exon,$transcript);
  $self->_add_exon_biotype($exon,$transcript) ; 
  
  if ($self->{_internal_index} != 1) {
    #print "Got more than one exon in cluster\n";
  }
  $self->{'_type'} = undef;
}

# this is a variation on the method above.
# calls contains_exon_with_dbid_and_dbname
# instead of contains_exon
sub add_exon_if_not_present {
  my ($self,$exon,$transcript, $ignore_strand) = @_;

  if (!$self->contains_exon_with_dbid_and_dbname($exon)) {
    $self->_add_new_exon($exon,$ignore_strand);
    # print "Added exon " . $exon->dbID . "\n";
  } else {
    # do not return here or you will not add
    # the transcript_reference to _transcripthash
  }
  $self->_add_transcript_reference($exon,$transcript);
  $self->_add_exon_biotype($exon,$transcript) ;

  if ($self->{_internal_index} != 1) {
    #print "Got more than one exon in cluster\n";
  }
  $self->{'_type'} = undef;
}

sub _add_new_exon {
  my ($self,$exon, $ignore_strand) = @_;

  if (!defined($self->start) || $exon->start < $self->start) {
    $self->start($exon->start);
  }
  if (!defined($self->end) || $exon->end > $self->end) {
    $self->end($exon->end);
  }
  if (!$ignore_strand) {
    if (!defined($self->strand)) {
      $self->strand($exon->strand);
    } elsif ($self->strand != $exon->strand) {
    #  confess("Trying to add exon with strand ". $exon->strand . " to cluster with strand " . $self->strand);
      carp("Trying to add exon with strand ". $exon->strand . " to cluster with strand " . $self->strand);
    }
  }
  $self->{_exonhash}{"$exon"} = $self->{_internal_index}++;
  #$self->{_exonidhash}{$exon->stable_id} = $exon;
  if  (exists $self->{_exonidhash}{$exon->dbID.$exon->adaptor->db->dbc->dbname}) { 
    $self->throw("Error : there seem to be exons with the same dbID and dbname in the databases ".$exon->dbID." ".$exon->adaptor->db->dbname) ; 
  }
  $self->{_exonidhash}{$exon->dbID.$exon->adaptor->db->dbc->dbname} = $exon;
  $self->{_all_exons_in_cluster}{$exon}=$exon ; 
}


=head2 _add_transcript_reference

Name : _add_transcript_reference($exon, $transcript) 
Arg[0] :  Bio::EnsEMBL::Exon
Arg[0] : Bio::EnsEMBL::Transcript
Function : called by add_exon , buids relation href{$transcript}=\@exons

=cut 

sub _add_transcript_reference {
  my ($self,$exon,$transcript) = @_;
  # if there's not already a reference to transcript stored make an arrayref 
  if (!$self->contains_transcript($transcript)) {
    $self->{_transcripthash}{"$transcript"} = [];
  }
  # store exons of transcript (key: transcript) 
  push @{$self->{_transcripthash}{"$transcript"}}, $exon;
  $self->{_transcriptrefhash}{"$transcript"} = $transcript;
}



=head2 transcript_from_ref

Name transcript_from_ref
Arg :      String "Bio::EnsEMBL::Transcript(HASHXXXXXX)" used as key in _transcriptrefhash
Function : Returns for a given STRING the reference to a Bio::EnsEMBL::Transcript-Object
Returnval: Bio::EnsEMBL::Transcript

=cut

sub transcript_from_ref {
  my ($self,$ref) = @_;
  return $self->{_transcriptrefhash}{$ref}
}


sub contains_exon {
  my ($self,$exon) = @_;
  return  (exists $self->{_exonhash}{"$exon"});
}

# this is a variation on the method above. contains_exon looks in the 
# _exonhash which keys on the memory location of the Exon objects. This means 
# that, when we fetch the same exon from the same db at two different times,
# the _exonhash will not recognise that they are the same because their location
# in memory is different. contains_exon only works when the Exon object has not 
# gone out of scope.
sub contains_exon_with_dbid_and_dbname {
  my ($self,$exon) = @_;
  return (exists $self->{_exonidhash}{"".$exon->dbID.$exon->adaptor->db->dbc->dbname});
}

=head2 get_all_Exons_in_ExonCluster

Name      : get_all_Exons_in_ExonCluster () 
Arg[0]    : Bio::EnsEMBL::Analysis::Runnable::Condense_EST::ExonCluster  
Function  : returns an Array of all Exons (all types) which belong to an ExonCluster 
           if Arg[1] is supplied the list of returned Bio::EnsEMBL::Exon Objects will 
           be filtered and only exons >= the supplied rank will be returned 
Returnval : Arrayref to Array of  Bio::EnsEMBL::Exon Objects

=cut

sub get_all_Exons_in_ExonCluster{ 
  my ($self ) = @_ ; 

  my @all_exons_in_cluster = values %{$self->{_all_exons_in_cluster}} ;  
  return \@all_exons_in_cluster ; 
}  



=head2 get_all_Exons_of_EvidenceSet

   Name      : get_all_Exons_of_EvidenceSet
   Arg       : String describing Name of evidence_set
   Function  : returns an Array of all Exons (all types) which belong to an ExonCluster 
               and belong to the specified Evidence_set. 
               be filtered and only exons >= the supplied rank will be returned 
   Returnval : Arrayref to Array of  Bio::EnsEMBL::Exon Objects

=cut


sub get_all_Exons_of_EvidenceSet{
  my ($self, $setname) = @_ ;

  my @result ;  
  for my $e (@{ $self->get_all_Exons_in_ExonCluster }) {
    if ($e->ev_set =~m/$setname/){
      push @result, $e ; 
    }
  }
  @result = sort { $b->seq_region_length <=> $a->seq_region_length} @result ; 
  return \@result ; 
}




sub get_all_Exons {
  my $self =shift;
  # return array of exons in cluster stored by their id 
  return [values %{$self->{_exonidhash}} ] ; 
}


sub contains_transcript {
  my ($self,$transcript) = @_;
  return  (exists $self->{_transcripthash}{"$transcript"});
}


# returns transcsripts which exactly hit the exon
sub transcripts_containing_exon {
  my ($self,$exon) = @_;
  my @transcripts;
  my %transhash =  %{ $self->each_transcripts_exons } ; 
 TRANS:
  foreach my $trans (keys %transhash) {
    foreach my $e (@{$transhash{$trans}}) {
      if ($e == $exon) {
        push @transcripts, $self->transcript_from_ref($trans); 
        next TRANS;
      }
    }
  }
  return \@transcripts;
}


sub get_transcripts_having_this_Exon_in_ExonCluster {
  my ($self,$exon) = @_;
  # self = exoncluster, exon = exon
  my @transcripts;
  my %transhash =  %{ $self->each_transcripts_exons }  ;
 
  #print "\n\nEC-start: ".$self->start."\tEC-end: " .$self->end . "\n" ;
 TRANS:foreach my $trans (keys %transhash) {
    #print "New Transcript:\n";
    unless ( $trans =~m/PredictionTranscript/ ) { 
    
      # we only processs "real" transcripts, no predictionTranscripts from ab-initio
      foreach my $ex_to_test (@{$transhash{$trans}}) {
        #print " testing exon: start: " . $ex_to_test->start . "\t" . $ex_to_test->end ."\n";
        if($ex_to_test->stable_id eq $exon->stable_id){
          #print "  pushing \t".$self->transcript_from_ref($trans)->dbID."\n";
          push @transcripts, $self->transcript_from_ref($trans); 
          next TRANS;
        }
      }
    }
  }
  return \@transcripts;
}


sub get_prediction_transcripts_which_have_exon_in_ExonCluster {
  my ($self,$exon) = @_;
  # self = exoncluster, exon = exon
  my @transcripts;
  my %transhash =  %{ $self->each_transcripts_exons } ; 
 
  #print "\n\nEC-start: ".$self->start."\tEC-end: " .$self->end . "\n" ;
 TRANS:foreach my $trans (keys %transhash) {
    #print "New Transcript:\n";
    if  ( $trans=~m/PredictionTranscript/) { 
      # we only processs "real" transcripts, no predictionTranscripts from ab-initio
      foreach my $ex_to_test (@{$transhash{$trans}}) {
        #print " testing exon: start: " . $ex_to_test->start . "\t" . $ex_to_test->end ."\n";
        if($ex_to_test->start >= $self->start &&  $ex_to_test->end <= $self->end){
          #print "  pushing \t".$self->transcript_from_ref($trans)->dbID."\n";
          push @transcripts, $self->transcript_from_ref($trans); 
          next TRANS;
        }
      }
    }
  }
  return \@transcripts;
}


sub check_if_ExonCluster_has_est_evidence {
  my ( $self,$ev_sets ) = @_ ; 

  my @est_bioytpes  = @{ $$ev_sets{'est'} }; 
  my %est_bt ; 
  @est_bt{@est_bioytpes} = () ;  
  my %transhash =  %{ $self->each_transcripts_exons } ; 
  my $cluster_has_real_evidence = 0 ; 
  my @rank_of_trans ; 

  TRANS: for my $trans_refname (keys %transhash) {
    my $trans = $self->transcript_from_ref( $trans_refname ); 
    my $trans_biotype = $trans->biotype() ; 

    if (exists $est_bt{$trans->biotype}) { 
      #print "ec holds exon which belongs to ev_set \'est\' ".$trans->biotype."\n" ;  
      return 1 ; 
    }
  }
  # print "ec holds no gene which belong to ev-set \'est\' \n" ;  
  return 0 ; 
}


# returns array of all transcripts which share the same exon 
#

sub unique_exon_combinations {
  my $self = shift;
  my $ignore_strand = shift;

  my %unique_combs;

  my %transhash =  %{ $self->each_transcripts_exons } ; 
  foreach my $trans (keys %transhash) {
    my $keystr;
    if (!$ignore_strand) {
      @{$transhash{$trans}} = sort _by_stranded_start @{$transhash{$trans}};
    } else {
      @{$transhash{$trans}} = sort _sort_by_forward_start @{$transhash{$trans}};
    }

    foreach my $exon (@{$transhash{$trans}}) {
      $keystr .= ":" . $self->{_exonhash}{$exon};
    }
    if (! exists($unique_combs{$keystr})) {
      $unique_combs{$keystr}{exons} = [@{$transhash{$trans}}];
    } 
    push @{$unique_combs{$keystr}{transcripts}}, $self->transcript_from_ref($trans); 
  }
  return [values %unique_combs]; 
}

sub each_transcripts_exons {
  my $self = shift;

  return $self->{_transcripthash};
}

sub _determine_type {
  my ($self, $ignore_strand) = @_;

  my @combs = @{ $self->unique_exon_combinations($ignore_strand) } ; 
  
  if (scalar(@combs) == 1) {
    $self->{'_type'} = 0;
    return 0;
  }

  my $maxexoncount = 0;
  foreach my $comb (@combs) {
    if (scalar(@{$comb->{exons}} > $maxexoncount)) {
      $maxexoncount = scalar(@{$comb->{exons}});
    }
  }

# All clusters have one exon which all have the same coordinates
# In this case the difference should be in the phase or the 
# translation start position - if there's no difference in these then
# this is an exon duplication.

  if ($maxexoncount == 1) {
    my $failed = 0;
    my $first_start;
    my $first_end;
    my $first_phase;
    my $subtype;
    foreach my $comb (@combs) {
      if (!defined($first_start)) {
        $first_start = $comb->{exons}[0]->start;
        $first_end   = $comb->{exons}[0]->end;
        $first_phase = $comb->{exons}[0]->phase;
      }
      if ($comb->{exons}[0]->start != $first_start ||
          $comb->{exons}[0]->end   != $first_end) {
        $failed = 1;
        last;
      } else {
        if ($comb->{exons}[0]->phase != $first_phase) {
          $subtype .= "p";
        } else {
          my $countterm = 0;
          foreach my $trans (@{$comb->{transcripts}}) {
            if ($trans->translation->start_Exon == $comb->{exons}[0] ||
                $trans->translation->end_Exon == $comb->{exons}[0]) {
              $countterm++;
            }
          }

          if ($countterm == scalar(@{$comb->{transcripts}})) {
            $subtype .= "t";
          } elsif ($countterm) {
            print "WARNING: Exon which isn't always the terminal exon in a translation\n";
          }
        }
      }
    }

    if (!$failed) {
      #print "Type 1 $subtype ";
      #print_clust_summary(@combs);
      $self->{'_type'} = 1;
      return 1;
    }
  }


# One or more clusters have ONLY one exon in all its transcripts
# Other clusters have longer exon and can have surrounding exons

  if ($maxexoncount == 1) {
# Foreach comb determine exon extents and number of exons in
# transcripts
    my @allsingle;
    my @notallsingle;
    foreach my $comb (@combs) {
      my $nsingle = 0;
      foreach my $trans (@{$comb->{transcripts}}) {
        if (scalar(@{$trans->get_all_Exons}) == 1) {
          $nsingle++;
        }
      }
      if ($nsingle == scalar(@{$comb->{transcripts}})) {
        push @allsingle,$comb;
      } else {
        push @notallsingle,$comb;
      }
    }

    my $failed = 0;
    if (@allsingle) {
      my $notsingle_start = undef;
      my $notsingle_end   = undef;
      SINGLE:
      foreach my $comb (@allsingle) {
        foreach my $comp_comb (@notallsingle) {
          if ($comb->{exons}[0]->start < $comp_comb->{exons}[0]->start ||
              $comb->{exons}[0]->end   > $comp_comb->{exons}[0]->end) {
            $failed = 1;
            last SINGLE;
          } elsif (!defined($notsingle_start)) {
            $notsingle_start = $comp_comb->{exons}[0]->start;
            $notsingle_end   = $comp_comb->{exons}[0]->end;
          } elsif ($notsingle_start != $comp_comb->{exons}[0]->start ||
                   $notsingle_end   != $comp_comb->{exons}[0]->end) {
            $failed = 1;
            last SINGLE;
          }
        }
      }
    } else {
      $failed = 1;
    }

    if (!$failed) {
      #print "Type 2 ";
      #print_clust_summary(@combs);
      $self->{'_type'} = 2;
      return 2;
    }
  }

# One cluster has a single exon which is the terminal exon in all its
# transcripts
# Other clusters has single (possibly) longer exon which has splice side
# boundary conserved. It may have extra exons to non spliced side
  if ($maxexoncount == 1) {
    if ($ignore_strand) {
      throw("This method is not designed to handle cases where we ignore strand");
    }
    my @allterminal;
    my @notallterminal;

    my $failed = 0;
    my $first_end = undef;
    COMB:
    foreach my $comb (@combs) {
      my $nterminal = 0;
      my $end = undef;
      foreach my $trans (@{$comb->{transcripts}}) {
        if (scalar(@{$trans->get_all_Exons}) == 1) {
          $nterminal++;
          print "Single exon transcript\n";
        } elsif ($comb->{exons}[0]->strand == 1) {
          if ($comb->{exons}[0] == $trans->start_Exon &&
                   (!defined($end) || ($end == 0))) {
            $end = 0;
            if (!defined($first_end)) {
              $first_end = $end;
            }
            $nterminal++;
          } elsif ($comb->{exons}[0] == $trans->end_Exon &&
                   (!defined($end) || $end == 1)) {
            $end = 1;
            if (!defined($first_end)) {
              $first_end = $end;
            }
            $nterminal++;
          } else {
            print "Failed forward strand condition\n";
          }
        } else {
          if ($comb->{exons}[0] == $trans->start_Exon &&
                   (!defined($end) || ($end == 1))) {
            $end = 1;
            if (!defined($first_end)) {
              $first_end = $end;
            }
            $nterminal++;
          } elsif ($comb->{exons}[0] == $trans->end_Exon &&
                   (!defined($end) || $end == 0)) {
            $end = 0;
            if (!defined($first_end)) {
              $first_end = $end;
            }
            $nterminal++;
          } else {
            print "Failed reverse strand condition\n";
          }
        }
        if (defined($end) && $end != $first_end) {
          $failed = 1;  
          last COMB;
        }
      }
      if ($nterminal == scalar(@{$comb->{transcripts}})) {
        push @allterminal,$comb;
      } else {
        push @notallterminal,$comb;
      }
    }

    if (!defined($first_end)) {
      $failed = 1;
    }

    if (!$failed) {
      if (scalar(@allterminal)) {
# first find (possibly) conserved end and max extent of 
# non conserved end
        my $conspos; 
        if ($first_end == 0) {
          $conspos = $allterminal[0]->{exons}[0]->end;
        } else {
          $conspos = $allterminal[0]->{exons}[0]->start;
        }

        my $nonconspos = undef;
        foreach my $comb (@allterminal) {
          if ($first_end == 0) {
            if (!defined($nonconspos) || $comb->{exons}[0]->start < $nonconspos) {
              $nonconspos = $comb->{exons}[0]->start;
            }
          } else {
            if (!defined($nonconspos) || $comb->{exons}[0]->end > $nonconspos) {
              $nonconspos = $comb->{exons}[0]->end;
            }
          }
        }

# now check that non of the combs with extra exons on the non conserved side
# 
        TERMINAL:
        foreach my $comb (@combs) {
          if ($first_end == 0 && $comb->{exons}[0]->end != $conspos ) {
            $failed = 1;
            last TERMINAL; 
          } elsif ($first_end == 1 && $comb->{exons}[0]->start != $conspos) {
            $failed = 1;
            last TERMINAL; 
          }
        }
      } else {
        $failed = 1;
      }
    }

    if (!$failed) {
      #print "Type 3 ";
      #print_clust_summary(@combs);
      $self->{'_type'} = 3;
      return 3;
    }
  }

# All have a single exon but do not match any of above cases
  if ($maxexoncount == 1) {
    #print "Type 4 ";
    #print_clust_summary(@combs);
    $self->{'_type'} = 4;
    return 4;
  } else {

# Clusters with multiple exons

#Ones which have very short introns (probably frameshift introns)
    foreach my $comb (@combs) {
      if (scalar(@{$comb->{exons}}) > 1) {
        my $prev_exon = undef; 
        my @exons = sort {$a->start <=> $b->start} @{$comb->{exons}};
        foreach my $exon (@exons) {
          if (defined($prev_exon)) {
            if ($exon->start - $prev_exon->end > 8) {
              #print "Type 6 (has multiple exon clusters) ";
              #print_clust_summary(@combs);
              $self->{'_type'} = 6;
              return 6;
            }
          }
          $prev_exon = $exon;
        }
      }
    }

    #print "Type 5 (has multiple frameshift exons) ";
    #print_clust_summary(@combs);
    $self->{'_type'} = 5;
    return 5;
  }
}
sub print_clust_summary {
  my (@combs) = @_;

  foreach my $comb (@combs) {
    print ": "; 
    foreach my $exon (@{$comb->{exons}}) {
      print $exon->dbID. " ";
    } 
    print ": "; 
  }
  print "\n";
}

sub _sort_by_forward_start {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;

  # we ignore strand
  $alow = $a->start;
  $ahigh = $a->end;
  $blow = $b->start;
  $bhigh = $b->end;

  if ($alow != $blow) {
    return $alow <=> $blow;
  } else {
    return $ahigh <=> $bhigh;
  }
}

sub _by_stranded_start {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;

  if ($a->strand != $b->strand) {
    confess("Mixed strand in sort comparison routine - Bad!");
  }

  if ($a->strand == 1) {
    $alow = $a->start;
    $ahigh = $a->end;
  } else {
    $blow = $a->start;
    $bhigh = $a->end;
  }

  if ($b->strand == 1) {
    $blow = $b->start;
    $bhigh = $b->end;
  } else {
    $alow = $b->start;
    $ahigh = $b->end;
  }

  if ($a->strand == 1) {
    if ($alow != $blow) {
      return $alow <=> $blow;
    } else {
      return $ahigh <=> $bhigh;
    }
  } else {
    if ($ahigh != $bhigh) {
      return $ahigh <=> $bhigh;
    } else {
      return $alow <=> $blow;
    }
  }
}

sub _add_exon_biotype{
  my ( $self, $exon, $transcript ) = @_  ; 
  unless (exists $self->{_exon_2_biotype}{$exon}) {
    $self->{_exon_2_biotype}{$exon}=$transcript->biotype ; 
  } else {
    my $old_bt = $self->{_exon_2_biotype}{$exon}; 
    my $new_bt = $transcript->biotype ;   
    if ($old_bt eq $new_bt) { 
      $self->{_exon_2_biotype}{$exon}=$transcript->biotype ; 
    } else { 
      $self->{_exon_2_biotype}{$exon}=$transcript->biotype ; 
      warning("Changing biotype of exon : $old_bt changed to $new_bt\n" ) ;  # jhv 
    }
  }
}

# could also do this with my $ex = ${$ec->get_all_Exons_in_ExonCluster}{$exon}}
# print $ex->biotype() ; 
#
sub get_biotype_of_Exon{
  my ($self,$exon) = @_ ; 
  return $self->{_exon_2_biotype}{$exon}; 
}   

1;
