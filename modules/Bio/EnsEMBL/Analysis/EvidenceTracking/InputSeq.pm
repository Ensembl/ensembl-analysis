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

Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq - Input sequence to the pipeline

=head1 SYNOPSIS

use Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq;

$inputseq = Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq->new(
    -hit_name =>$name,
    -external_db_id => $external_db_id,
    -molecule_type => 'PROTEIN',
    -submission_date => '2011-03-03',
    );

=head1 DESCRIPTION

  Keep the information about the input sequence. Usually it is added in the database
  automatically when Uniprot is running. Or any time you want to track sequences for
  an analysis you will have to add them in the database

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;
use Data::Dumper;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

 Arg [1]    : $db_id, int
 Arg [2]    : $hit_name, string
 Arg [3]    : $external_db_id, int
 Arg [4]    : $molecule_type, enum, 'PROTEIN', 'MRNA', 'EST'
 Arg [5]    : $submission_date, date YYYY-MM-DD
 Arg [7]    : $adaptor, Bio::EnsEMBL::Analysis::DBSQL::InputSeqAdaptor object
 Example    : $inputseq = Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq->new(
    -hit_name =>$name,
    -external_db_id => $external_db_id,
    -molecule_type => 'PROTEIN',
    -submission_date => '2011-03-03',
    );
 Description: Constructor
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $hit_name, $external_db_id, 
      $molecule_type, $submission_date, $adaptor) =
          rearrange([qw(DBID
                        HIT_NAME       
                        EXTERNAL_DB_ID
                        MOLECULE_TYPE
                        SUBMISSION_DATE
                        ADAPTOR
                        )],@args);

  $self->dbID      ( $id ) if (defined $id);
  $self->hit_name ( $hit_name );
  $self->external_db_id ( $external_db_id );
  $self->molecule_type ( $molecule_type );
  $self->submission_date ( $submission_date );
  $self->adaptor   ( $adaptor );
  return $self; # success - we hope!
}


=head2 hit_name

 Arg [1]    : $hit_name, int [optional]
 Example    : $evidence->hit_name($hit_name);
 Description: Getter/Setter for the hit name
 Returntype : integer, the hit name
 Exceptions : 


=cut

sub hit_name {
  my $self = shift;
  $self->{'hit_name'} = shift if ( @_ );
  return $self->{'hit_name'};
}

=head2 external_db_id

 Arg [1]    : $external_db_id, int [optional]
 Example    : $evidence->external_db_id($external_db_id);
 Description: Getter/Setter for the external db id
 Returntype : integer, the external db id
 Exceptions : 


=cut

sub external_db_id {
  my $self = shift;
  $self->{'external_db_id'} = shift if ( @_ );
  return $self->{'external_db_id'};
}

=head2 molecule_type

 Arg [1]    : $molecule_type, string [optional] 'PROTEIN', 'MRNA', 'EST'
 Example    : $evidence->molecule_type($molecule_type);
 Description: Getter/Setter for the molecule type
 Returntype : string, the molecule type
 Exceptions : 


=cut

sub molecule_type {
  my $self = shift;
  $self->{'molecule_type'} = shift if ( @_ );
  return $self->{'molecule_type'};
}

=head2 submission_date

 Arg [1]    : $submission_date, date [optional] YYYY-MM-DD
 Example    : $evidence->submission_date($submission_date);
 Description: Getter/Setter for the submission date
 Returntype : date, the submission date
 Exceptions : 


=cut

sub submission_date {
  my $self = shift;
  $self->{'submission_date'} = shift if ( @_ );
  return $self->{'submission_date'};
}

=head2 is_stored

 Arg [1]    : $inputseq_adaptor, a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor object
 Example    : $inputseq->is_stored($inputseq_adaptor);
 Description: Test if the object is alreday stored
 Returntype : boolean
 Exceptions : 


=cut

sub is_stored {
  my $self = shift;
  my $db = shift;

  if($db and $db->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor')) {
    if ($db->fetch_by_hit_name($self->hit_name)) {
        print STDERR $self->hit_name, ' is stored!',"\n";
        return 1;
    }
    return 0;
  }
  else {
    throw('db argument must be a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor');
  }
}

1;
