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

Bio::EnsEMBL::Analysis::Tools::Pmatch::ProteinHit - 

=head1 SYNOPSIS


=head1 DESCRIPTION

  holds a pmatch protein hit - simply the name(idenitifier) of the protein
  and a list of ContigHit objects

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Tools::Pmatch::ProteinHit;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

=head2 new

 Title   : new
 Usage   :
 Function:
 Example :
 Returns : 
 Args    : constructor


=cut

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($id) = $self->_rearrange(['ID'], @args);

  $self->throw("No id") unless defined $id;
  $self->id($id);
  
  return $self;

}

=head2 id

 Title   : id
 Usage   :
 Function: get/set for contig id
 Example :
 Returns : 
 Args    :


=cut

sub id {
  my ($self,$id) = @_;
  if ($id) {
    $self->{'id'} = $id;
  }
  return $self->{'id'};
}

=head2 add_ContigHit

 Title   : add_ContigHit
 Usage   :
 Function: adds a ContigHit into $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub add_ContigHit {
  my ($self,$hit) = @_;
  $self->throw('No contig hit') unless defined $hit;
  $self->throw('$contig is not a Bio::EnsEMBL::Analysis::Tools::Pmatch::ContigHit') unless $hit->isa("Bio::EnsEMBL::Analysis::Tools::Pmatch::ContigHit");
  $self->{_contig_hits}{$hit->id()} = $hit;
}

=head2 each_ContigHit

 Title   : each_ContigHit
 Usage   :
 Function: returns all entries in $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub each_ContigHit {
  my ($self) = @_;
  return values %{$self->{_contig_hits}};
}


=head2 get_ContigHit

 Title   : get_ContigHit
 Usage   :
 Function: returns entries for a particular contig in $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub get_ContigHit {
  my ($self,$contig) = @_;
  return ($self->{_contig_hits}{$contig}) if defined $contig;
  return undef;
}

1;
