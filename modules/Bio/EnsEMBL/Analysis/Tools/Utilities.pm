
=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Utilities 

- base class which exports utility methods which don't take Bio::XX objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::Utilities qw(shuffle); 

  or 

  use Bio::EnsEMBL::Analysis::Tools:Utilities

  to get all methods

=head1 DESCRIPTION

This is a class which exports Utility methods for genebuilding and
other gene manupulation purposes. 

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::Utilities; 

use strict;
use warnings;
use Exporter;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw( shuffle );


=head2 shuffle 

  Arg [1]   : Reference to Array
  Function  : randomizes the order of an array 
  Returntype: arrayref 
  Exceptions: none
  Example   : 

=cut

sub shuffle {
  my $tref = shift ;
  my $i = @$tref ;
  while ($i--) {
     my $j = int rand ($i+1);
     @$tref[$i,$j] = @$tref[$j,$i];
  }
  return $tref ;
}




1;
