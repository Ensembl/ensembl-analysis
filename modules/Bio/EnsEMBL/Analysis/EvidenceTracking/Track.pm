=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::EvidenceTracking::Track - Tracking system

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::EvidenceTracking::Track;

  sub fetch_input{
      ...
      my $track = Bio::EnsEMBL::Analysis::EvidenceTracking::Track->new(
                -runnabledb     => $blastminigenewise,
                -tracking       => 1,
                -evidence_names => $ra_evidences);
      ...
      $self->track($track);
  }

  sub run {
      ...
      $self->track->update($supporting_evidence);
      ...
  }

  sub write_output {
      ...
      $self->track->write_tracks;
      ...
  }

=head1 DESCRIPTION
  
  The track module should be instanciated in the fetch_input method of
  the RunnableDB. It will populate the evidence set with the sequence
  that will be fetch. Then update the evidence set whenever you like.
  Finally, write the tracks during the write_output method of the
  RunnableDB.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Track;

use strict;

use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence;
use Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;


=head2 new

 Arg [1]    : $runnabledb, a Bio::EnsEMBL::Analysis::RunnableDB object
 Arg [2]    : $default_track, a Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack object
 Arg [3]    : $ra_evidences, listref of string
 Example    : $track = Bio::EnsEMBL::Analysis::EvidenceTracking::Track->new(
                -runnabledb     => $blastminigenewise,
                -tracking       => 1,
                -evidence_names => $ra_evidences);
 Description: Constructor
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::Track
 Exceptions : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($runnabledb, $tracking, $ra_evidence_names, $file) =
          rearrange([qw(RUNNABLEDB TRACKING EVIDENCE_NAMES FILE)],@args);

  if (!$tracking) {
      $self->tracking(0);
      return $self;
  }
  $self->tracking(1);
   if ($file) {
     open(RF, 'grep \> '.$file.' | ') or throw('Could not grep the sequences name for '.$file);
     my @seq_names = <RF>;
     close(RF);
     chomp @seq_names;
     s/^>(\S+).*$/$1/ for @seq_names;
     $ra_evidence_names = \@seq_names;
   }

  $self->db($runnabledb->db);
  my $analysisrun_adaptor = $self->db->get_AnalysisRunAdaptor;
  my $reason_adaptor = $self->db->get_ReasonAdaptor;
  $self->_reasons($reason_adaptor->fetch_all);
  my $analysis = $runnabledb->analysis;
  my $analysis_id = $analysis->dbID;
  my $analysis_run;
  $analysis_run = $analysisrun_adaptor->get_current_analysis_run_by_analysis_id($analysis_id);
  if (!$analysis_run) {
      $analysis_run = Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun->new(
        -analysis_id => 0,
        -analysis_run_id => 0
        );
  }
  $self->analysis_run($analysis_run);
  $self->input_id($runnabledb->input_id);

  foreach my $name (@{$ra_evidence_names}) {
      print STDERR "NAME: $name\n";
      my $input_seq = $self->get_input_seq($name);
      my $evidence = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
                                  -input_seq  => $input_seq,
                                  -is_aligned => 'u',
                                  );
      my $evidence_track = Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
          -evidence     => $evidence,
          -analysis_run => $analysis_run,
          -reason       => $self->get_reason_by_code('Unknown'),
          -input_id     => $runnabledb->input_id
          );
      $self->add_track($name, $evidence_track);
  }
  return $self; # success - we hope!
}

=head2 db

 Example    : $self->db($db);
 Description: Get the DBAdaptor for the evidence tracking system database
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor object


=cut

sub db {
    my $self = shift;
    my $db = shift if (@_);

    if ( !exists $self->{'db'}) {
        throw('Not a object!') unless ($db and $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
        my $evidencetracking_db = new Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor (
            -dbconn => $db->dbc
            );
        $self->{'db'} = $evidencetracking_db;
    }
    return $self->{'db'};
}

=head2 update

 Arg [1]    : $support, supporting evidence
 Arg [2]    : $reason_id, int (optional)
 Example    : $self->update($feature, $reason_id);
 Description: It update the evidence with the supporting
              evidence provided. Add the position on the slice
              change the reason if it has been provided
 Returntype : void
 Exceptions : 


=cut

sub update {
    my $self = shift;
    return undef unless $self->tracking;
    my $support = shift;
    my $reason_id = shift;
    print STDERR "Entry: ", Dumper($support);
    
    my $name;
    if (ref($support) eq '') {
        $support = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
            -input_seq  => $self->get_input_seq($support),
            -is_aligned => 'n'
            );
    }
    elsif ($support->isa('Bio::EnsEMBL::Analysis::Tools::Pmatch::MergedHit')) {
        $support = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
            -input_seq         => $self->get_input_seq($support->target),
            -seq_region_name   => $support->seq_region_name,
            -seq_region_start  => $support->qstart,
            -seq_region_end    => $support->qend,
            -seq_region_strand => $support->strand,
            -is_aligned        => 'y'
            );
    }
    elsif ($support->isa('Bio::EnsEMBL::BaseAlignFeature')) {
        $support = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
            -input_seq         => $self->get_input_seq($support->hseqname),
            -seq_region_name   => $support->seq_region_name,
            -seq_region_start  => $support->start,
            -seq_region_end    => $support->end,
            -seq_region_strand => $support->strand,
            -is_aligned        => 'y'
            );
    }
    my $seq_region_id     = $support->seq_region_name;
    my $seq_region_start  = $support->seq_region_start;
    my $seq_region_end    = $support->seq_region_end;
    my $seq_region_strand = $support->seq_region_strand;
    if ($support->isa('Bio::EnsEMBL::Feature')) {
        $support = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
            -input_seq         => $self->get_input_seq($support->hseqname),
            -seq_region_name   => $seq_region_id,
            -seq_region_start  => $seq_region_start,
            -seq_region_end    => $seq_region_end,
            -seq_region_strand => $seq_region_strand,
            -is_aligned        => 'y'
            );
    }
    $name = $support->input_seq->hit_name;
    print STDERR "Starting update $name\n";
    if (defined $seq_region_start and $self->get_track($name.$seq_region_start.$seq_region_end.$seq_region_strand)) {
        print STDERR "has aligned evidence\n";
        $name .= $seq_region_start.$seq_region_end.$seq_region_strand;
    }
    elsif ($self->get_track($name)) {
        print STDERR "has primary evidence\n";
        if (defined $seq_region_start) {
            print STDERR "\thas start\n";
            $self->delete_track($name);
            $name .= $seq_region_start.$seq_region_end.$seq_region_strand;
        }
        else {
            print STDERR 'ZWRT:has evidence aligned: ', $support->is_aligned, "\t", $self->get_track($name)->evidence->is_aligned, "\n";
        }
    }
    else {
        print STDERR "has no evidence\n";
        if (defined $seq_region_start) {
            print STDERR "\thas start\n";
            $name .= $seq_region_start.$seq_region_end.$seq_region_strand;
        }
        else {
            print STDERR 'XKWZ:Is aligned: ', $support->is_aligned, "\n";
        }
    }
    my $evidence_track = Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
          -evidence     => $support,
          -analysis_run => $self->analysis_run,
          -input_id     => $self->input_id
          );
    $self->add_track($name, $evidence_track);
    $self->update_reason($name, $reason_id) if ($reason_id);
}

=head2 write_tracks

 Arg [1]    : $dba, a Bio::EnsEMBL::Analysis::EvidenceTracking::DBAdaptor object;
 Example    : $self->write_tracks($dba);
 Description: It write all the tracks in the database. Call it after or before you write
              the gene/transcripts/...
 Returntype : void
 Exceptions : 


=cut

sub write_tracks {
    my $self = shift;

    return undef unless $self->tracking;
    if (! $self->tracks) {
        warning("No tracking data!");
        return;
    }
    print STDERR 'Writing!!!!!', "\n";
    my $evidencetrack_adaptor = $self->db->get_EvidenceTrackAdaptor;
    foreach my $evidence_track (values %{$self->tracks}) {
        $evidencetrack_adaptor->store($evidence_track);
    }
}

=head2 all_to_noalign

 Example    : $track->all_to_noalign;
 Description: Set to non align all the evidence that are marked unknown
 Returntype : void
 Exceptions : 


=cut

sub all_to_noalign {
    my $self = shift;

    foreach my $evidence_track (values %{$self->tracks}) {
        $self->update_reason($evidence_track, "NoAlignment");
        $evidence_track->evidence->is_aligned('n');
    }
}

=head2 get_track

 Arg [1]    : $name, string
 Example    : my $track = $self->get_track($name);
 Description: Return the track given the specific key, $name or
              $name.$seq_region_name.$seq_region_start.$seq_region_end.$seq_region_strand
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object,
              undef if it does not exists or if the tracking system is off
 Exceptions : 


=cut

sub get_track {
    my $self = shift;
    return undef unless $self->tracking;
    my $name = shift;

    return $self->{'_tracks'}->{$name} if (exists $self->{'_tracks'}->{$name});
    return undef;
}

=head2 delete_track

 Arg [1]    : $name, string
 Example    : my $track = $self->delete_track($name);
 Description: Return the track given the specific key, $name or
              $name.$seq_region_name.$seq_region_start.$seq_region_end.$seq_region_strand
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object,
              undef if it does not exists or if the tracking system is off
 Exceptions : 


=cut

sub delete_track {
    my $self = shift;
    return undef unless $self->tracking;
    my $name = shift;

    delete $self->{'_tracks'}->{$name};
}

=head2 add_track

 Arg [1]    : $track, a Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Example    : $self->add_track($track);
 Description: Add the track to the set of tracks
 Returntype : void
 Exceptions : if not given a Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object 


=cut

sub add_track {
  my $self = shift;
  return undef unless $self->tracking;
  my ($name, $track) = @_;

  if ($track and !$track->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack')) {
      throw('add_track is waiting for a Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack object!');
  }
  $self->{'_tracks'}->{$name} = $track;
}

=head2 update_reason

 Arg [1]    : $evidence, a Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Arg [1]    : $reason_id, int
 Example    : $self->update_reason($evidence, $reason_id);
 Description: Set the reason to the code provided, for now it update
              the first evidence track of the array because in that
              case you are not supposed to have several tracks
 Returntype : void
 Exceptions : 


=cut

sub update_reason {
    my $self = shift;
    return undef unless $self->tracking;
    my ($name, $reason_id) = @_;

    $self->get_track($name)->reason($self->get_reason_by_code($reason_id));
}

=head2 has_tracks

 Example    : my $bool = $self->has_tracks;
 Description: Check if the track got tracks
 Returntype : Boolean
 Exceptions : 


=cut

sub has_tracks {
    my $self = shift;
    return undef unless $self->tracking;
    return exists $self->{'_tracks'};
}

=head2 tracking

 Example    : my $bool = $self->tracking;
 Description: Check if we use the tracking system
 Returntype : Boolean
 Exceptions : 


=cut

sub tracking {
    my $self = shift;

    $self->{'tracking'} = shift if (@_);
    return $self->{'tracking'};
}

=head2  tracks

 Arg [1]    : hashref of Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Example    : $self->tracks($ra_tracks);
 Description: Getter/Setter of the hash of tracks. The key is compose by the name of the
              track if it was not aligned and by the concatenation of the name, seq_region name,
              seq_region start, seq_region end and seq_region strand to have a unique key
 Returntype : hashref of Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Exceptions : if it's not a listref


=cut

sub tracks {
  my $self = shift;
  return undef unless $self->tracking;
  my $ra_tracks = shift;
  if ( $ra_tracks ) {
#      throw('You should provide a reference to an array '.ref($ra_tracks)) unless (ref($ra_tracks) eq 'REF');
      $self->{'_tracks'} = $ra_tracks;
  }
  return $self->{'_tracks'};
}

=head2 _reasons

 Arg [1]    : an listref of Bio::EnsEMBL::Analysis::EvidenceTracking::Reason object
 Example    : $self->_reasons($reason_adaptor->fetch_all);
 Description: Store all the possible reasons, to avoid fetching the object for each evidence tracked
 Returntype : void
 Exceptions : 


=cut

sub _reasons {
  my $self = shift;
  return undef unless $self->tracking;
  my $ra_reasons = shift;
  if ( $ra_reasons ) {
#      throw('You should provide a reference to an array '.ref($ra_reasons)) unless (ref($ra_reasons) eq 'REF');
      foreach my $reason (@{$ra_reasons}) {
          $self->{'_reasons'}->{$reason->code} = $reason;
      }
  }
}

=head2 get_input_seq

 Arg [1]    : $name, a string
 Example    : $input_seq = $track->get_input_seq($name);
 Description: Get the input sequence for an evidence. If it does not already exists, it will
              be created.
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq object
 Exceptions : 


=cut

sub get_input_seq {
  my $self = shift;
  return undef unless $self->tracking;
  my $name = shift;
      print STDERR "name: ", $name, "\n";
  if ($name and !exists $self->{'_input_seq'}->{$name}) {
      my $input_seq = Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq->new(
                        -hit_name => $name
                        );
      $self->{'_input_seq'}->{$name} = $input_seq;
  }
  return $self->{'_input_seq'}->{$name};
}

=head2 get_reason_by_code

 Arg [1]    : $code, a string like NoAlignment, Accepted, LowCoverage,...
 Example    : $reason = $track->get_reason_by_code($code);
 Description: Get a specific reason to set to an evidence by its code. The codes can be 
              found in the reasons table or in the documentation
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::Reason object
 Exceptions : 


=cut

sub get_reason_by_code {
  my $self = shift;
  return undef unless $self->tracking;
  my $code = shift;
  return $self->{'_reasons'}->{$code};
}

=head2 analysis_run

 Arg [1]    : $analysis_run [optional], an Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun object
 Example    : $track->analysis_run($analysis_run);
 Description: Getter/Setter for the analysis_run object that represent the run
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun object
 Exceptions : throw if the object is not of the type Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun


=cut

sub analysis_run {
  my $self = shift;
  return undef unless $self->tracking;
  my $analysis_run = shift;
  if ($analysis_run) {
      throw('Not a Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun object')
        unless $analysis_run->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun');
      $self->{'analysis_run'} = $analysis_run;
  }
  return $self->{'analysis_run'};
}

=head2 input_id

 Arg [1]    : $input_id [optional], a string
 Example    : $track->input_id($input_id);
 Description: Getter/Setter for the input_id of the current job
 Returntype : String
 Exceptions : 


=cut

sub input_id {
  my $self = shift;
  return undef unless $self->tracking;
  $self->{'_input_id'} = shift if ( @_ );
  return $self->{'_input_id'} if (exists $self->{'_input_id'});
}


1;
