# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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

=head1 CONTACT

please send any questions to http://lists.ensembl.org/mailman/listinfo/dev

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
