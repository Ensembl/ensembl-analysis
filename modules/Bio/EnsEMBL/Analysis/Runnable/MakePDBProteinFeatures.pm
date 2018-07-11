=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::MakePDBProteinFeatures -

=head1 SYNOPSIS


=head1 DESCRIPTION

This module makes protein features based on the PDB-UniProt mappings found in the EMBL-EBI PDB SIFTS data and the UniProt-ENSP mappings found in the GIFTS database in order to make the link between PDB and ENSP having a PDB entry as a protein feature for a given ENSP protein.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::MakePDBProteinFeatures;

use warnings;
use strict;

# Bio::DB::HTS::Faidx used in Bio::EnsEMBL::GIFTS::DB needs Perl 5.14.2
use 5.014002;
use parent ('Bio::EnsEMBL::Analysis::Runnable');
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::GIFTS::DB qw(fetch_latest_uniprot_enst_perfect_matches);
use Bio::EnsEMBL::ProteinFeature;

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($core_dba,$gifts_dbc,$pdb_filepath,$species,$cs_version) = rearrange([qw(core_dba gifts_dbc pdb_filepath species cs_version)],@args);

    $self->{'core_dba'} = $core_dba;
    $self->{'gifts_dbc'} = $gifts_dbc;
    $self->{'pdb_filepath'} = $pdb_filepath;
    $self->{'species'} = $species;
    $self->{'cs_version'} = $cs_version;
    
    $self->{'pdb_info'} = undef;
    $self->{'perfect_matches'} = undef;

    return $self;
}

############################################################
#
# Analysis methods
#
############################################################

=head2 run

 Arg [1]    : None
 Description: Parse the PDB file, fetch the perfect matches between Ensembl ENST transcripts and UniProt proteins from the GIFTS database and make the corresponding protein features which link them. 
 Returntype : Integer, 1
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  $self->{'pdb_info'} = $self->parse_pdb_file();
  $self->{'perfect_matches'} = fetch_latest_uniprot_enst_perfect_matches($self->{'gifts_dbc'},$self->{'species'},$self->{'cs_version'});
  $self->make_protein_features();

  return 1;
}

sub parse_pdb_file() {
# Parse the CSV PDB file containing the following 9 columns:
# PDB CHAIN SP_PRIMARY RES_BEG RES_END PDB_BEG PDB_END SP_BEG SP_END
# The SIFTS release date is fetched from the lines starting with '#'.
# Lines starting with 'PDB' are ignored.
# Return reference to an array of hashes with SP_PRIMARY, PDB, CHAIN, RES_BEG, RES_END, SP_BEG, SP_END and SIFTS_RELEASE_DATE as keys.

  my $self = shift;

  open(my $pdb_fh,'<',$self->{'pdb_filepath'}) or die "Cannot open: ".$self->{'pdb_filepath'};
  my @pdb_info = ();
  my $sifts_release_date;
  
  while (my $line = <$pdb_fh>) {
    
    next if $line =~ /^PDB/;

    if ($line =~ /^#/) {
    # parse SIFTS release date
      (undef,$sifts_release_date) = split(/\s+/,$line);
    } else {
    # parse a PDB-UniProt line
      my ($pdb,$chain,$sp_primary,$res_beg,$res_end,undef,undef,$sp_beg,$sp_end) = split(/\s+/,$line);

      push(@pdb_info,{'PDB' => $pdb,
                      'CHAIN' => $chain,
                      'SP_PRIMARY' => $sp_primary,
                      'RES_BEG' => $res_beg,
                      'RES_END' => $res_end,
                      'SP_BEG' => $sp_beg,
                      'SP_END' => $sp_end,
                      'SIFTS_RELEASE_DATE' => $sifts_release_date
                     });
    }
  }
  return \@pdb_info;
}

sub make_protein_features() {
# create the protein features by linking
# the ENSP proteins and the PDB entries
# from 'perfect_matches' and 'pdb_info'
# It returns an array of hash references containing the
# translation db ID as key and the ProteinFeature as value.

  my $self= shift;

  my @pfs = ();

  # get list of transcript stable IDs from the keys of the perfect matches hash
  my @t_sids = keys %{$self->{'perfect_matches'}};

  # loop through all pdb lines and find the corresponding ENSTs (if any)
  # in the perfect matches hash,
  # fetch their proteins and add their protein features
  my $ta = $self->{'core_dba'}->get_TranscriptAdaptor();
  my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'sifts_import',
                                            -display_label => 'SIFTS import',
                                            -displayable => '1',
                                            -description => 'Protein features based on the PDB-UniProt mappings found in the EMBL-EBI PDB SIFTS data and the UniProt-ENSP mappings found in the GIFTS database.',
  );
  
  foreach my $pdb_line (@{$self->{'pdb_info'}}) {
    my $pdb_uniprot = $$pdb_line{'SP_PRIMARY'};    

    if (defined($self->{'perfect_matches'}{$pdb_uniprot})) {
      my @ensts = @{$self->{'perfect_matches'}{$pdb_uniprot}};
      if (scalar(@ensts) > 0) {
        foreach my $enst (@ensts) {
         
          my %pf_translation;
          my $t = $ta->fetch_by_stable_id($enst);

          if ($t) {
            my $translation = $t->translation();
            my $translation_sid = $translation->stable_id();

            if ($$pdb_line{'SP_BEG'} < $$pdb_line{'SP_END'}+1) {
            # ignore complex PDB-UniProt mappings that allow SP_BEG > SP_END
              my $pf = Bio::EnsEMBL::ProteinFeature->new(
                      -start    => $$pdb_line{'SP_BEG'},
                      -end      => $$pdb_line{'SP_END'},
                      -hseqname => $$pdb_line{'PDB'}.".".$$pdb_line{'CHAIN'},
                      -hstart   => $$pdb_line{'RES_BEG'},
                      -hend     => $$pdb_line{'RES_END'},
                      -analysis => $analysis,
                      -hdescription => "Via SIFTS (".$$pdb_line{'SIFTS_RELEASE_DATE'}.
                                       ") UniProt protein ".$$pdb_line{'SP_PRIMARY'}.
                                       " isoform exact match to Ensembl protein $translation_sid"
                   );
              $pf_translation{$translation->dbID()} = $pf;
              push(@pfs,\%pf_translation);
            }
          } # if t
        } # foreach my enst
      } # if scalar
    } # if ensts
  } # foreach my pdb_line
  
  $self->output(\@pfs);
}

1;
