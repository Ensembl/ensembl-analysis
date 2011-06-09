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

Bio::EnsEMBL::Analysis::EvidenceTracking::Reason - Store the reason why we kept/dismissed the sequence

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::EvidenceTracking::Reason;

  my $reason = Bio::EnsEMBL::Analysis::EvidenceTracking::Reason->new(
    -info => 'internal stop codon'
    );

=head1 DESCRIPTION

  Store the reason why we kept/dismissed the sequence.
  If we kept the sequence the code will be lower than 100.
  If we dismissed it the code will be 100 or higher.
  If the code is 0, it means we don't know what happened.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Reason;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::ReasonAdaptor;

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

 Arg [1]    : $dbid, int
 Arg [2]    : $info, string
 Arg [6]    : $adaptor, Bio::EnsEMBL::Analysis::DBSQL::AnalysisRunAdaptor object
 Example    : $reason = ->new(
              -info => 'internal stop codon'
            );
 Description: Constructor
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun
 Exceptions : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $info, $adaptor) = rearrange([qw(DBID INFO ADAPTOR)],@args);

  $self->dbID    ( $id );
  $self->info  ( $info );
  return $self; # success - we hope!
}

=head2 info

 Arg [1]    : $info, string [optional]
 Example    : $reason->info($info);
 Description: Getter/Setter for the reason
 Returntype : string, the reason
 Exceptions : 


=cut

sub info {
  my $self = shift;
  $self->{'info'} = shift if ( @_ );
  return $self->{'info'};
}

1;
