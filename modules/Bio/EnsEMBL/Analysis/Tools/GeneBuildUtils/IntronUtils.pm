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

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils - utilities for transcript objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils qw(get_splice_sites);

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils 

  to get all methods

=head1 DESCRIPTION

All methods in this class should take a Bio::EnsEMBL::Intron
object as their first argument.

The methods provided should carry out some standard 
functionality for said objects such as printing info or
getting splice sites

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(seq_region_coord_string id);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             print_Intron
             Intron_info
             get_splice_sites
             intron_length_less_than_maximum);


=head2 Intron_info

  Arg [1]   : Bio::EnsEMBL::Intron
  Arg [2]   : string, indent of \ts to put infrount of string
  Function  : return a string with coords or intron
  Returntype: string
  Exceptions: throw is not given an argument
  Example   : 

=cut


sub Intron_info{
  my ($intron, $indent) = @_;
  throw("Must be passed an intron") if(!$intron);
  $indent = '' if(!$indent);
  my $coord_string = seq_region_coord_string($intron);
  my $id = id($intron);
  return $indent."INTRON: ".$coord_string;
}


=head2 print_Intron

  Arg [1]   : Bio::EnsEMBL::Intron
  Arg [2]   : string, indent
  Function  : print out info about given Intron
  Returntype: none
  Exceptions: n/a
  Example   : 

=cut



sub print_Intron{
  my ($intron, $indent) = @_;
  print Intron_info($intron, $indent)."\n";
}




=head2 get_splice_sites

  Arg [1]   : Bio::EnsEMBL::Intron
  Function  : get the donor and acceptor splice sites for the
  given intron
  Returntype: array, (two strings in an array the donor and
                      acceptor splice sites)
  Exceptions: none
  Example   : 

=cut



sub get_splice_sites{
  my ($intron) = @_;
  my $intron_slice = $intron->feature_Slice;
  if (defined $intron_slice) {
    print STDERR "Intron slice defined for intron with next exon id ".$intron->next_Exon->dbID." and pre exon id ".$intron->prev_Exon->dbID."\n";
  } else {
    print STDERR "Intron slice NOT defined for intron with next exon id ".$intron->next_Exon->dbID." and pre exon id ".$intron->prev_Exon->dbID."\n";
  }
  my $donor_seq    = uc($intron_slice->subseq(1,2));
  my $acceptor_seq = uc($intron_slice->subseq($intron->length-1,$intron->length));
  return $donor_seq, $acceptor_seq;
}



=head2 intron_length_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Intron
  Arg [2]   : Int, max length of intron
  Function  : check intron length against specificed maximum
  Returntype: boolean 0 for longer than max, 1 for less than
  Exceptions: none
  Example   : 

=cut



sub intron_length_less_than_maximum{
  my ($intron, $max_length) = @_;
  if($intron->length >= $max_length){
    warning("This intron is longer than ".$max_length);
    return 0;
  }
  return 1;
}

1;
