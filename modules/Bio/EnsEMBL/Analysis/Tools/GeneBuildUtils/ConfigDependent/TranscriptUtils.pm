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

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ConfigDependent::TranscriptUtils - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ConfigDependent::TranscriptUtils - utilities for transcript objects which depend on config files 

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ConfigDependent::TranscriptUtils qw( method );

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ConfigDependent::TranscriptUtils 

  to get all methods

=head1 DESCRIPTION

All methods in this class should take a Bio::EnsEMBL::Transcript
object as their first argument.

The methods provided should carry out some standard 
functionality for said objects such as printing info, and 
cloning and checking phase consistency or splice sites etc

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut

package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ConfigDependent::TranscriptUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump); 

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw( print_peptide );
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw( id );
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg;
use Bio::SeqIO;

use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);

@EXPORT = qw(
             low_complexity_less_than_maximum
            );




=head2 low_complexity_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, maximum low complexity
  Function  : calculate how much low complexity a
  transcripts translation has and check if it is above
  the specificed threshold
  Returntype: boolean, 1 for pass 0 for fail
  Exceptions: none
  Example   : 

=cut



sub low_complexity_less_than_maximum{
  my ($transcript, $complexity_threshold) = @_;
  my $peptide = $transcript->translate;
  my $hit_name = ${$transcript->get_all_supporting_features}[0]->hseqname;
  my $seg = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg->new
    (
     -query => $peptide,
     -analysis => Bio::EnsEMBL::Analysis->new
     (
      -logic_name => 'seg',
      -program_file => 'seg',
     )
    );
  $seg->run;
  my $low_complexity = $seg->get_low_complexity_length;
  logger_info(id($transcript)." ($hit_name) has ".$low_complexity.
              " low complexity sequence");
  #print_peptide($transcript);
  #print id($transcript)." has ".$low_complexity." low complexity sequence compared to".
  #  " ".$complexity_threshold."\n";
  if($low_complexity >= $complexity_threshold){
    warn(id($transcript)."($hit_name)'s low ".
            "complexity (".$low_complexity.") is above ".
            "the threshold ".$complexity_threshold.
            "\n");
    return 0;
  }
  return 1;
}

1;
