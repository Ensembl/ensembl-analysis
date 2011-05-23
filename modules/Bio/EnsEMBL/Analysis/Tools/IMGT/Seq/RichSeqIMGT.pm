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

Bio::EnsEMBL::Analysis::Tools::IMGT::Seq::RichSeqIMGT - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Tools::IMGT::Seq::RichSeqIMGT;
use strict;

use base qw(Bio::Seq::RichSeq);


sub new {
  # standard new call..
  my($caller,@args) = @_;
  my $self = $caller->SUPER::new(@args);
  
  my ($data_class) = $self->_rearrange([qw(DATA_CLASS
					    )],
					@args);

  defined $data_class and $self->data_class($data_class);

  return $self;
}


sub data_class {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_data_class'} = $value;
    }
    return $obj->{'_data_class'};

}


1;
