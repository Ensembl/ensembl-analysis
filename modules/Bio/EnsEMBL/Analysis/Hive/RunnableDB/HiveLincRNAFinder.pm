=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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




Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAFinder - 

=head1 SYNOPSIS

Find RNAseq models that don't overlap with protein coding (models) predictions and store them as lincRNA candidates (rnaseq)

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAFinder;

use warnings;
use strict;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils	qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::LincRNA qw(get_genes_of_biotypes_by_db_hash_ref) ;  

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');




=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder
  Function  : instatiates a lincRNAFinder object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder 
  Exceptions: 
  Example   : 

=cut

sub fetch_input {
	my ($self) = @_;

  # This call will set the config file parameters. Note this will set REFGB (which overrides the
  # value in $self->db and OUTDB
	$self->create_analysis;
	my $dba = $self->hrdb_get_dba( $self->param('output_db') );
	$self->hrdb_set_con( $dba, 'output_db' );

	# get cdnas and convert them to single transcript genes. The convert_to_single_transcript_gene located: ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Tools/GeneBuildUtils/GeneUtils.pm
	my $new_cdna = $self->get_genes_of_biotypes_by_db_hash_ref( $self->NEW_SET_1_CDNA );
	my @single_transcript_cdnas =	map { @{ convert_to_single_transcript_gene($_) } } @$new_cdna;  

	# get protein_coding genes and convert them to single transcript genes
	my $new_set_prot =	$self->get_genes_of_biotypes_by_db_hash_ref( $self->NEW_SET_2_PROT );
	my @single_trans_pc =	map { @{ convert_to_single_transcript_gene($_) } } @$new_set_prot;

	# create runnable
	my $runnable = Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder->new(
		-query    => $self->query,
		-analysis => $self->analysis,
	);

  # add hash-keys and hash-values directly to the $runnable hashref. quicker than using constructors...
	$runnable->set_1_cdna_genes( \@single_transcript_cdnas );
	$runnable->set_2_prot_genes( \@single_trans_pc );
	$runnable->ignore_strand( $self->CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND ); # it is not working know, but this need to change! 
	$runnable->find_single_exon_candidates($self->FIND_SINGLE_EXON_LINCRNA_CANDIDATES );
	$runnable->maximum_translation_length_ratio($self->MAXIMUM_TRANSLATION_LENGTH_RATIO );
	$runnable->max_translations_stored_per_gene($self->MAX_TRANSLATIONS_PER_GENE );

	$self->runnable($runnable);
}

sub write_output {
	my ($self) = @_;

	my $adaptor = $self->hrdb_get_con('output_db')->get_GeneAdaptor;
	print "Final output is: \n"; 
	print "have " . @{ $self->output } . " genes to write\n";

  my $analysis = $self->analysis;
  GENE: foreach my $gene ( @{ $self->output } ) {
		if ( !defined $gene->get_all_Transcripts ) {
			$self->throw(" gene does not have any transcripts ....\n");
		}
    empty_Gene($gene);
    attach_Analysis_to_Gene($gene, $analysis);
		$gene->biotype( $self->OUTPUT_BIOTYPE );
    $adaptor->store($gene);
	}
}


#######
#CHECKS
#######


=head2  

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut

sub NEW_SET_1_CDNA {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->param( 'NEW_SET_1_CDNA', $arg );
	}
	return $self->param('NEW_SET_1_CDNA');
}

sub NEW_SET_2_PROT {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->param( 'NEW_SET_2_PROT', $arg );
	}
	return $self->param('NEW_SET_2_PROT');
}

sub OUTPUT_DB {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->param( 'OUTPUT_DB', $arg );
	}
	return $self->param('OUTPUT_DB');
}

sub OUTPUT_BIOTYPE {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'OUTPUT_BIOTYPE', $arg );
	}
	return $self->param('OUTPUT_BIOTYPE');
}

sub DEBUG_OUTPUT_DB {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'DEBUG_OUTPUT_DB', $arg );
	}
	return $self->param('DEBUG_OUTPUT_DB');
}

sub MAX_TRANSLATIONS_PER_GENE {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'MAX_TRANSLATIONS_PER_GENE', $arg );
	}
	return $self->param('MAX_TRANSLATIONS_PER_GENE');
}

sub MAXIMUM_TRANSLATION_LENGTH_RATIO {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'MAXIMUM_TRANSLATION_LENGTH_RATIO', $arg );
	}
	return $self->param('MAXIMUM_TRANSLATION_LENGTH_RATIO');
}

sub WRITE_DEBUG_OUTPUT {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'WRITE_DEBUG_OUTPUT', $arg );
	}
	return $self->param('WRITE_DEBUG_OUTPUT');
}

sub CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND', $arg );
	}
	return $self->param('CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND');
}

sub FIND_SINGLE_EXON_LINCRNA_CANDIDATES {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'FIND_SINGLE_EXON_LINCRNA_CANDIDATES', $arg );
	}
	return $self->param('FIND_SINGLE_EXON_LINCRNA_CANDIDATES');
}


1;
