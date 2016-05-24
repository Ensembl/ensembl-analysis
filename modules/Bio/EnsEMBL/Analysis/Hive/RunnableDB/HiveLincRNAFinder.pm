
=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAFinder;

use warnings;
use vars qw(@ISA);
use strict;
use Data::Dumper;

use Bio::EnsEMBL::Hive::Utils ('destringify');
use Bio::EnsEMBL::Analysis;

# use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
# use Bio::EnsEMBL::Analysis::Config::GeneBuild::lincRNAFinder;
use Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);

# use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils	qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

# @ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild Bio::EnsEMBL::Analysis::RunnableDB);

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

	# $self->query($self->fetch_sequence);
	print "DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 0A !! ...\n\n";

# This call will set the config file parameters. Note this will set REFGB (which overrides the
# value in $self->db and OUTDB
	$self->hive_set_config;

	print	"DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 0B !! ...\n\n";

	# get cdnas and convert them to single transcript genes
	my $new_cdna = $self->get_genes_of_biotypes_by_db_hash_ref( $self->NEW_SET_1_CDNA );

	print "DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 0C !! ...\n\n";
	my @single_transcript_cdnas =
		map { @{ convert_to_single_transcript_gene($_) } } @$new_cdna;

	print "DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 0D !! ...\n\n";

	# get protein_coding genes and convert them to single transcript genes
	my $new_set_prot =
		$self->get_genes_of_biotypes_by_db_hash_ref( $self->NEW_SET_2_PROT );

	print	"DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 0E !! ...\n\n";
	my @single_trans_pc =	map { @{ convert_to_single_transcript_gene($_) } } @$new_set_prot;

	# create runnable
	my $runnable = Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder->new(
		-query    => $self->query,
		-analysis => $self->analysis,
	);

# add hash-keys and hash-values directly to the $runnable hashref. quicker than using constructors...
	print	"DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 1A !! ...\n\n";

	$runnable->set_1_cdna_genes( \@single_transcript_cdnas );
	$runnable->set_2_prot_genes( \@single_trans_pc );

	print	"DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 2o !! ...\n\n";

	$runnable->ignore_strand( $self->CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND );
	$runnable->find_single_exon_candidates($self->FIND_SINGLE_EXON_LINCRNA_CANDIDATES );

	print	"DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 3o !! ...\n\n";

	# $runnable->check_cdna_overlap_with_both_K4_K36($self->CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36 );
	# $runnable->check_cdna_overlap_with_multiple_K36($self->CHECK_CDNA_OVERLAP_WITH_MULTI_K36 );

	# print "DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 4o !! ...\n\n";

	$runnable->maximum_translation_length_ratio($self->MAXIMUM_TRANSLATION_LENGTH_RATIO );
	$runnable->max_translations_stored_per_gene($self->MAX_TRANSLATIONS_PER_GENE );

	print	"DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 5o !! ...\n\n";

	$runnable->efg_clustering_with_cdna_analysis(	$self->create_analysis_object( $self->DEBUG_LG_EFG_CLUSTERING_WITH_CDNA ) );
	$runnable->unclustered_efg_analysis( $self->create_analysis_object( $self->DEBUG_LG_EFG_UNCLUSTERED ) );

	print	"DEBUG_BK::RunnableDB::lincRNAFinder: fetch input :: 6o !! ...\n\n";

	$self->runnable($runnable);
}

sub get_genes_of_biotypes_by_db_hash_ref {
	my ( $self, $href ) = @_;

	my %dbnames_2_biotypes = %$href;
	my @genes_to_fetch;
	foreach my $db_hash_key ( keys %dbnames_2_biotypes ) {

		my @biotypes_to_fetch = @{ $dbnames_2_biotypes{$db_hash_key} };
		my $set_db            = $self->hrdb_get_dba( $self->param($db_hash_key) );
		my $dna_dba = $self->hrdb_get_dba( $self->param('reference_db') );
		if ($dna_dba) {
			$set_db->dnadb($dna_dba);
		}

		my $test_id = $self->param('iid');
		my $slice =	$self->fetch_sequence( $test_id, $set_db, undef, undef, $db_hash_key );

		# implementation of fetch_all_biotypes ....
		my $fetch_all_biotypes_flag;
		foreach my $biotype (@biotypes_to_fetch) {
			if ( $biotype =~ m/fetch_all_biotypes/ ) {
				$fetch_all_biotypes_flag = 1;
			}
		}
		
		if ($fetch_all_biotypes_flag) {
			print "fetching ALL biotypes for slice out of db $db_hash_key :\n";
			my $genes = $slice->get_all_Genes( undef, undef, 1 );
			push @genes_to_fetch, @$genes;
			print scalar(@genes_to_fetch) . " genes fetched in total\n";
		}
		else {
			foreach my $biotype (@biotypes_to_fetch) {
        print "DEBUGG::get_genes_of_biotypes_by_db_hash_ref A2 print flag:: $biotype \n";				
				my $genes = $slice->get_all_Genes_by_type( $biotype, undef, 1 );
				if ( @$genes == 0 ) {
					warning("No genes of biotype $biotype found in $set_db\n");
				}
				print "$db_hash_key [ " . $set_db->dbc->dbname	. " ] Retrieved " . @$genes	. " of type "	. $biotype	. "... in input_id: " . $test_id . "\n";

				if ( $self->param('verbose') ) {
				  print "$db_hash_key [ " . $set_db->dbc->dbname	. " ] Retrieved " . @$genes	. " of type "	. $biotype	. "... in input_id: " . $test_id . "\n";
				}
				push @genes_to_fetch, @$genes;
			}
		}
	}
	return \@genes_to_fetch;
}

sub write_output {
	my ($self) = @_;
	my $dba = $self->hrdb_get_dba( $self->param('lincRNA_output_db') );
	$self->hrdb_set_con( $dba, 'lincRNA_output_db' );

	my $adaptor = $self->hrdb_get_con('lincRNA_output_db')->get_GeneAdaptor;

	print "have " . @{ $self->output } . " genes to write\n";

	my $sucessful_count = 0;

GENE: foreach my $gene ( @{ $self->output } ) {
		print "DEBUG::HiveLincRNAFinder:: write_output gene 1: $gene \n";
		if ( !defined $gene->get_all_Transcripts ) {
			throw(" gene does not have any transcripts ....\n");
		}

		my @tr     = @{ $gene->get_all_Transcripts };
		my $max_ex = 0;
		# print "DEBUG::HiveLincRNAFinder:: write_output gene 2: $gene \n";

		for (@tr) {
			$_->status(undef);
			$_->analysis( $self->analysis );
		}
		# print "DEBUG::HiveLincRNAFinder:: write_output gene 3: $gene \n";

		$gene->biotype( $self->OUTPUT_BIOTYPE );
		$gene->status(undef);
		$gene->analysis( $self->analysis );
		# print "DEBUG::HiveLincRNAFinder:: write_output gene 4: $gene \n";

		eval { $adaptor->store($gene); };
		# print "DEBUG::HiveLincRNAFinder:: write_output gene 5: $gene \n";

		if ($@) {
			print "DEBUG::HiveLincRNAFinder:: write_output gene 6: $gene \n";
			warning( "Failed to write gene " . id($gene) . " "	. coord_string($gene)	. " $@" );
		}
		else {
			# print "DEBUG::HiveLincRNAFinder:: write_output gene 7: $gene \n";
			$sucessful_count++;
			print "STORED LINCRNA GENE " . $gene->dbID;
		}

		# print "DEBUG::HiveLincRNAFinder:: write_output gene 8: $gene \n";

	}

	print "DEBUG::HiveLincRNAFinder:: output:: " . $sucessful_count	. " genes written to " . " @ \n";

	if ( $sucessful_count != @{ $self->output } ) {
		throw("Failed to write some genes");
	}
}



=head2 read_and_check_config

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::lincRNAFinder
  Arg [2]   : hashref from config file
  Function  : call the superclass method to set all the varibles and carry
  out some sanity checking
  Returntype: N/A
  Exceptions: throws if certain key variables aren't set properly
  Example   : 

=cut


#######
#CHECKS
#######
#  foreach my $var(qw(NEW_SET_1_CDNA NEW_SET_2_PROT OUTPUT_DB OUTPUT_BIOTYPE EFG_FEATURE_NAMES EXTEND_EFG_FEATURES )){
#    throw("RunnableDB::lincRNAFinder $var config variable is not defined")
#      unless($self->$var);
#  }

# }

# HIVE check
sub hive_set_config {
	my $self = shift;

	# Throw is these aren't present as they should both be defined
	unless ( $self->param_is_defined('logic_name')
		&& $self->param_is_defined('module') )
	{
		throw(
"You must define 'logic_name' and 'module' in the parameters hash of your analysis in the pipeline config file, "
				. "even if they are already defined in the analysis hash itself. This is because the hive will not allow the runnableDB "
				. "to read values of the analysis hash unless they are in the parameters hash. However we need to have a logic name to "
				. "write the genes to and this should also include the module name even if it isn't strictly necessary"
		);
	}

# Make an analysis object and set it, this will allow the module to write to the output db
	my $analysis = new Bio::EnsEMBL::Analysis(
		-logic_name => $self->param('logic_name'),
		-module     => $self->param('module'),
	);
	$self->analysis($analysis);

# Now loop through all the keys in the parameters hash and set anything that can be set
	my $config_hash = $self->param('config_settings');
	foreach my $config_key ( keys( %{$config_hash} ) ) {
		if ( defined &$config_key ) {
			$self->$config_key( $config_hash->{$config_key} );
		}
		else {
			throw(
"You have a key defined in the config_settings hash (in the analysis hash in the pipeline config) that does "
					. "not have a corresponding getter/setter subroutine. Either remove the key or add the getter/setter. Offending "
					. "key:\n"
					. $config_key );
		}
	}
}

=head2 NEW_SET_1_CDNA 

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

sub EFG_FEATURE_DB {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'EFG_FEATURE_DB', $arg );
	}
	return $self->param('EFG_FEATURE_DB');
}

sub EFG_FEATURE_NAMES {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'EFG_FEATURE_NAMES', $arg );
	}
	return $self->param('EFG_FEATURE_NAMES');
}

sub EXTEND_EFG_FEATURES {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'EXTEND_EFG_FEATURES', $arg );
	}
	return $self->param('EXTEND_EFG_FEATURES');
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

sub DEBUG_LG_EFG_CLUSTERING_WITH_CDNA {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'DEBUG_LG_EFG_CLUSTERING_WITH_CDNA', $arg );
	}
	return $self->param('DEBUG_LG_EFG_CLUSTERING_WITH_CDNA');
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

sub DEBUG_LG_EFG_UNCLUSTERED {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'DEBUG_LG_EFG_UNCLUSTERED', $arg );
	}
	return $self->param('DEBUG_LG_EFG_UNCLUSTERED');
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

sub CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36', $arg );
	}
	return $self->param('CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36');
}

sub CHECK_CDNA_OVERLAP_WITH_MULTI_K36 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'CHECK_CDNA_OVERLAP_WITH_MULTI_K36', $arg );
	}
	return $self->param('CHECK_CDNA_OVERLAP_WITH_MULTI_K36');
}

sub FIND_SINGLE_EXON_LINCRNA_CANDIDATES {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'FIND_SINGLE_EXON_LINCRNA_CANDIDATES', $arg );
	}
	return $self->param('FIND_SINGLE_EXON_LINCRNA_CANDIDATES');
}

sub create_analysis_object {
	my ( $self, $logic_name ) = @_;

	# my $efg_out = $self->get_dbadaptor($self->DB_ANYNAME_REPL);
	# my $aa = $efg_out->get_AnalysisAdaptor() ;
	# my $analysis =  $aa->fetch_by_logic_name($logic_name) ;
	# if ( $analysis ) {
	#   return $analysis ;
	# }
	# need to create analysis object first
	my $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => $logic_name, );
	return $analysis;
}

sub memory_i_am_using {
	print "DEBUG::MEMORY:: ";
	require Carp;
	my $size = `ps -p $$ -h -o size`;
	print "$size \n";
}

# use vars '$AUTOLOAD';
# sub AUTOLOAD {
#  my ($self,$val) = @_;
#  (my $routine_name=$AUTOLOAD)=~s/.*:://; #trim package name
#  $self->{$routine_name}=$val if defined $val ;
#  return $self->{$routine_name} ;
# }
# sub DESTROY {} # required due to AUTOLOAD

1;
