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

Bio::EnsEMBL::Analysis::EvidenceTracking::Track - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Track;

use vars qw(@ISA);
use strict;


use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;

#@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

 Arg [1]    : $default_track, a 
 Arg [2]    : $ra_evidences, listref of 
 Example    : $track = ->new(-default_track => $default_track);
 Description: 
 Returntype : 
 Exceptions : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($runnabledb, $tracking, $ra_evidence_names) =
          rearrange([qw(RUNNABLEDB TRACKING EVIDENCE_NAMES)],@args);

  if (!$tracking) {
      $self->tracking(0);
      return $self;
  }
  $self->tracking(1);

  return undef unless $self->tracking;
  $self->db($runnabledb->db);
  my $analysisrun_adaptor = $self->db->get_AnalysisRunAdaptor;
  my $analysis = $runnabledb->analysis;
  my $analysis_id = $analysis->dbID;
  my $current_analysis = $analysisrun_adaptor->get_current_analysis_run_by_analysis_id($analysis_id);
  my $evidence_track = Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
      -analysis_run_id => $current_analysis->dbID,
      -analysis_id     => $analysis_id,
      -reason_id       => 0,
      -is_last         => 1,
      -input_id        => $runnabledb->input_id
      );

  $self->default_track($evidence_track);
  foreach my $name (@{$ra_evidence_names}) {
      my $evidence = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
          -hit_name     => $name,
          -is_aligned   => 0,
          );
      $evidence->add_track($self->default_track);
      $self->update($evidence);
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
    my $db = shift;

    if ( !exists $self->{'db'}) {
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
 Returntype : 
 Exceptions : 


=cut

sub update {
    my $self = shift;
    return undef unless $self->tracking;
    my $support = shift;
    my $reason_id = shift;
    
    my $name;
    my $seq_region_id     = $support->seq_region_name;
    my $seq_region_start  = $support->seq_region_start;
    my $seq_region_end    = $support->seq_region_end;
    my $seq_region_strand = $support->seq_region_strand;
    if ($support->isa('Bio::EnsEMBL::Feature')) {
        $name = $support->hseqname;
        $support = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
            -hit_name          => $name,
            -seq_region_name   => $seq_region_id,
            -seq_region_start  => $seq_region_start,
            -seq_region_end    => $seq_region_end,
            -seq_region_strand => $seq_region_strand,
            -is_aligned        => 1
            );
    }
    $name = $support->hit_name;
    if (! $support->has_tracks) {
        my $evidencetrack = $self->default_track;
        $support->add_track($evidencetrack);
    }
    
    
    
    if ($self->get_evidence($name)) {
        if (defined $seq_region_start) {
            my $evidence = $self->get_evidence($name);
            $evidence->seq_region_name($seq_region_id);
            $evidence->seq_region_start($seq_region_start);
            $evidence->seq_region_end($seq_region_end);
            $evidence->seq_region_strand($seq_region_strand);
            $evidence->is_aligned(1);
            $name .= $seq_region_start.$seq_region_end.$seq_region_strand;
            $self->add_evidence($evidence);
        }
    }
    elsif (defined $seq_region_start and $self->get_evidence($name.$seq_region_start.$seq_region_end.$seq_region_strand)) {
        $name .= $seq_region_start.$seq_region_end.$seq_region_strand;
    }
    else {
        if (defined $seq_region_start) {
            $name .= $seq_region_start.$seq_region_end.$seq_region_strand;
            $self->add_evidence($support);
        }
        else {
            $self->add_evidence($support);
        }
    }
    $self->update_reason($self->get_evidence($name), $reason_id) if ($reason_id);
}

=head2 write_tracks

 Arg [1]    : $dba, a Bio::EnsEMBL::Analysis::EvidenceTracking::DBAdaptor object;
 Example    : $self->write_tracks($dba);
 Description: It write all the tracks in the database. Call it after or before you write
              the gene/transcripts/...
 Returntype : 
 Exceptions : 


=cut

sub write_tracks {
    my $self = shift;

    return undef unless $self->tracking;
    if (! $self->evidences) {
        warning("No tracking data!");
        return;
    }
    my $evidence_adaptor     = $self->db->get_EvidenceAdaptor;
    my $trackevidence_adaptor = $self->db->get_EvidenceTrackAdaptor;
    foreach my $evidence (values %{$self->evidences}) {
        $evidence_adaptor->store($evidence) unless ($evidence->is_stored($self->db));
        $evidence->get_tracks->[0]->evidence_id($evidence->dbID);
        $trackevidence_adaptor->store($evidence->get_tracks->[0]);
    }
}

=head2 get_evidence

 Arg [1]    : $name, string
 Example    : my $evidence = $self->get_evidence($name);
 Description: Return the evidence given the specific key, $name or
              $name.$seq_region_name.$seq_region_start.$seq_region_end.$seq_region_strand
 Returntype : undef if it does not exists
 Exceptions : 


=cut

sub get_evidence {
    my $self = shift;
    return undef unless $self->tracking;
    my $name = shift;

    return $self->{'_evidences'}->{$name} if (exists $self->{'_evidences'}->{$name});
    return undef;
}

=head2 add_evidence

 Arg [1]    : $evidence, a Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Example    : $self->add_evidence($evidence);
 Description: Add the evidence to the 
 Returntype : 
 Exceptions : if not given a Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object 


=cut

sub add_evidence {
  my $self = shift;
  return undef unless $self->tracking;
  my $evidence = shift;

  if ($evidence and !$evidence->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence')) {
      throw('add_evidence is waiting for a Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object!');
  }
  if (defined $evidence->seq_region_start) {
      $self->{'_evidences'}->{$evidence->hit_name.$evidence->seq_region_start
                                .$evidence->seq_region_end.$evidence->seq_region_strand} = $evidence;
  }
  else {
      $self->{'_evidences'}->{$evidence->hit_name} = $evidence;
  }
}

=head2 update_reason

 Arg [1]    : $evidence, a Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Arg [1]    : $reason_id, int
 Example    : $self->update_reason($evidence, $reason_id);
 Description: Set the reason to the code provided, for now it update
              the first evidence track of the array because in that
              case you are not supposed to have several tracks
 Returntype : 
 Exceptions : 


=cut

sub update_reason {
    my $self = shift;
    return undef unless $self->tracking;
    my ($evidence, $reason_id) = @_;

    $evidence->get_tracks->[0]->reason_id($reason_id);
}

=head2 has_evidences

 Example    : my $bool = $self->has_evidences;
 Description: Check if the track got evidences
 Returntype : boolean
 Exceptions : 


=cut

sub has_evidences {
    my $self = shift;
    return undef unless $self->tracking;
    return exists $self->{'_evidences'};
}

=head2 tracking

 Example    : my $bool = $self->tracking;
 Description: Check if we use the tracking system
 Returntype : boolean
 Exceptions : 


=cut

sub tracking {
    my $self = shift;

    $self->{'tracking'} = shift if (@_);
    return exists $self->{'tracking'};
}

=head2 default_track

 Arg [1]    : a Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack
 Example    : my $evidence_track = $self->default_track;
 Description: Getter/Setter for the default track, generic data like the analysis_run_id
              is supposed to be set when you add the EvidenceTrack object. Call it whenever
              you create a new evidence track.
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack
 Exceptions : if not given a Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack object


=cut

sub default_track {
  my $self = shift;
  return undef unless $self->tracking;
  $self->{'default_track'} = shift if ( @_ );

  throw('Not a Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack object')
    unless $self->{'default_track'}->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack');
  return $self->{'default_track'}->clone;
}

=head2  evidences

 Arg [1]    : hashref of Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Example    : $self->evidences($ra_evidences);
 Description: Getter/Setter of the hash of evidences. The key is compose by the name of the
              evidence if it was not aligned and by the concatenation of the name, seq_region name,
              seq_region start, seq_region end and seq_region strand to have a unique key
 Returntype : hashref of Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Exceptions : if it's not a listref


=cut

sub evidences {
  my $self = shift;
  return undef unless $self->tracking;
  my $ra_evidences = shift;
  if ( $ra_evidences ) {
      throw('You should provide a not a reference to an array '.ref($ra_evidences)) unless (ref($ra_evidences) eq 'REF');
      $self->{'_evidences'} = $ra_evidences;
  }
  return $self->{'_evidences'};
}

1;
