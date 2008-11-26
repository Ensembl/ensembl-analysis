
=head1 NAME - Bio::EnsEMBL::Analysis::Tools::Otter::DBSQL::DBAdaptor 

  Inherits from the standard Ensembl Bio::EnsEMBL::DBSQL::DBAdaptor

  However - the get_available_adaptors method is overridden here
  so that we can read in the Otter DnaAlignFeature and
  DnaAlignFeatureHistory adaptors.

=cut

package Bio::EnsEMBL::Analysis::Tools::Otter::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::DBAdaptor);


=head2 get_available_adaptors

  Example    : my %pairs = %{$dba->get_available_adaptors()};
  Description: gets a hash of the available adaptors
  ReturnType : reference to a hash
  Exceptions : none
  Caller     : Bio::EnsEMBL::Utils::ConfigRegistry
  Status     : Stable

=cut 

sub get_available_adaptors {

  my %pairs = (
    # Firstly those that just have an adaptor named after there object
    # in the main DBSQL directory.
    map( { $_ => "Bio::EnsEMBL::DBSQL::${_}Adaptor" } qw(
        AffyFeature              AffyArray            AffyProbe
        Analysis                 ArchiveStableId      Attribute
        AssemblyExceptionFeature AssemblyMapper       CoordSystem
        CompressedSequence       DBEntry            
        DensityFeature           DensityType          Exon
        Gene                     KaryotypeBand        MiscSet
        MiscFeature              OligoArray           OligoFeature
        OligoProbe               PredictionTranscript PredictionExon
        ProteinFeature           ProteinAlignFeature  RepeatConsensus
        RepeatFeature            Sequence             SimpleFeature
        Slice                    SupportingFeature    Transcript
        TranscriptSupportingFeature Translation       UnmappedObject
        UnconventionalTranscriptAssociation
        AssemblySlice
        ) ),
    # Those whose adaptors are in Map::DBSQL
    map( { $_ => "Bio::EnsEMBL::Map::DBSQL::${_}Adaptor" } qw(
        Marker MarkerFeature QtlFeature Qtl Ditag DitagFeature
        ) ),

    # otter ones
    map( { $_ => "Bio::EnsEMBL::Analysis::Tools::Otter::DBSQL::${_}Adaptor" } qw(
        DnaAlignFeature DnaAlignFeatureHistory
        ) ),

    # Finally the exceptions... those that have non-standard mapping
    # between object / adaptor ....
    # 'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
    'MetaCoordContainer' => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
    'MetaContainer'      => 'Bio::EnsEMBL::DBSQL::MetaContainer',
    'SNP'                => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
    # Feature Collections:
    'GeneCollection'       => 'Bio::EnsEMBL::Collection::Gene',
    'TranscriptCollection' => 'Bio::EnsEMBL::Collection::Transcript',
    'ExonCollection'       => 'Bio::EnsEMBL::Collection::Exon',
    'RepeatFeatureCollection' =>
      'Bio::EnsEMBL::Collection::RepeatFeature' );

  return ( \%pairs );
} ## end sub get_available_adaptors


1;
