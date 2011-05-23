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

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::HomologyUtils - utilities for gene objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::HomologyUtils qw(clone_Gene);

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::HomologyUtils
  
  to get all methods

=head1 DESCRIPTION

All methods in this class should take a Bio::EnsEMBL::Compara::Homology
object as their first argument.

The methods provided should carry out some standard 
functionality for said objects such as printing info, and 
cloning

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut




package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::HomologyUtils;

use strict;
use warnings;
use Exporter;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             get_gene_obj_out_of_compara_homology_object
            );


use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);

sub get_gene_obj_out_of_compara_homology_object {
  my ( $homology, $species ) = @_ ;

  my $gene ;
  return $gene unless $homology ;

  for my $homology_member_obj ( @{$homology->gene_list}) {
     if ($homology_member_obj->genome_db->name eq $species ) {
       $gene = $homology_member_obj->get_Gene ;
     }
  }
  return $gene ;
}

1;
