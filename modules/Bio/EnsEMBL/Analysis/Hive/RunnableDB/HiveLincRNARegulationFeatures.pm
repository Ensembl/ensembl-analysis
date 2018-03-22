
=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::HiveLincRNARegulationFeatures - 

=head1 SYNOPSIS

Gets all lincRNA canditates and compare them against histone modification features. (mainly for human and mouse)

=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNARegulationFeatures;

use warnings;
use strict;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::LincRNA qw(get_genes_of_biotypes) ;  

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::HiveLincRNARegulationFeatures
  Function  : instatiates a lincRNAFinder object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::HiveLincRNARegulationFeatures 
  Exceptions: 
  Example   : 

=cut




sub fetch_input {
	my ($self) = @_;

	$self->hive_set_config;

	# checks for parameters if they exist in param! 
	$self->param_required('CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36'); # display_label check can be used, but I am not really checking for both now.
	$self->param_required('CHECK_CDNA_OVERLAP_WITH_MULTI_K36');   # display_label check can be used, but this is the default
	$self->param_required('DEBUG_LG_EFG_CLUSTERING_WITH_CDNA');
	$self->param_required('DEBUG_LG_EFG_UNCLUSTERED');
	$self->param_required('WRITE_DEBUG_OUTPUT'); # need that
	$self->param_required('EXTEND_EFG_FEATURES'); # need that

	my $set_db  = $self->hrdb_get_dba( $self->param('lincRNA_output_db') );
	my $dna_dba = $self->hrdb_get_dba( $self->param('reference_db') );
  $set_db->dnadb($dna_dba);
  $self->hrdb_set_con($set_db, 'lincRNA_output_db');

	my $dba = $self->hrdb_get_dba( $self->param('regulation_reform_db') );
	$self->hrdb_set_con( $dba, 'regulation_reform_db' );
	my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    %{$self->param('regulation_db')},
	);
  $efgdba->dnadb($dna_dba),
  $self->hrdb_set_con($efgdba, 'regulation_db');

	# get data
	my $biotype_tmp = $self->param_required('BIOTYPE_TO_CHECK');
	my $dbName_tmp  = 'lincRNA_output_db';
	my $new_cdna    = $self->get_genes_of_biotypes($biotype_tmp, $dbName_tmp);
	print "fetch input :: number of lincRNA genes: " . scalar(@$new_cdna) . "\n";

	# get simple features
	my $efg_simple_features = $self->get_efg_simple_features();
  $self->param('efg_simple_features', $efg_simple_features);
  $self->param('genes', $new_cdna);
}




sub run {
	my ($self) = @_;
  my $efg_simple_features = $self->param('efg_simple_features');
	$self->output($efg_simple_features) if ($self->param('WRITE_DEBUG_OUTPUT'));

	my $passed_cdnas_dbIDs = $self->score_cdna_overlap_with_multi_K36( $efg_simple_features, $self->param('genes') );

	print "-- passed cdnas:              " . scalar(@$passed_cdnas_dbIDs) . "\n";
	print "-- number of artificial genes "	. scalar(@$efg_simple_features) . "\n";

	my $store_genes_that_overlap_fungen =	1;    # the default is to store genes that overlap regulation elements
	$self->output( $passed_cdnas_dbIDs ) if ($store_genes_that_overlap_fungen);
}


sub write_output {
	my ($self) = @_;

	my $adaptor = $self->hrdb_get_con('regulation_reform_db')->get_GeneAdaptor;
	print "have " . @{ $self->output } . " genes to write\n";
  my $logic_name_to_be = "lincRNA_reg";
  my $analysis = Bio::EnsEMBL::Analysis->new(
                                         -logic_name => $logic_name_to_be,
                                         -displayable => 1
                                         );

	foreach my $gene ( @{ $self->output } ) {
    $gene->status(undef); 
    $gene->analysis($analysis);
		eval { 
			$adaptor->store($gene); 
		};
		if ($@) {
			$self->throw('Can t store this gene');
    }
	}
}


# get histone modification data from fungen database
sub get_efg_simple_features {
	my ($self) = @_;

	#Get fungen db adaptor
	my $efgdba = $self->hrdb_get_con('regulation_db');
	my $fsa = $efgdba->get_FeatureSetAdaptor;
  $self->warning("DBadaptor is not defined! \n") if ( !$fsa ) ;

  print "Fetching EFG domain data\n";
	my ( $fset, @fsets );
	my @feature_sets;    # = ('K562_H3K36me3_ENCODE_Broad_ccat_histone');
	my $all = 'all'; # I want all and then will filter. 

	if ( !defined $all ) {
		foreach my $feature_set (@feature_sets) {
			$fset = $fsa->fetch_by_name($feature_set);
			if ($fset) {
				push @fsets, $fsa->fetch_by_name($feature_set);
			}
			else {
				$self->throw("Could not find Feature Set:\t$feature_set");
			}
		}
	}
	else {
		@fsets = @{ $fsa->fetch_all_by_feature_class('annotated') };
	}
	print "I got : " . scalar(@fsets) . " fsets \n"; 
	return $self->convert_simple_features( \@fsets );
}



# convert and select fungen data to artificial genes. This way, we can find those fungen data that overlap our lincRNA predictions
sub convert_simple_features {
  my ( $self, $sf ) = @_;
	my @converted_sf;
	my @simple_features = @$sf;

	print scalar(@simple_features) . " simple features sets retrieved \n";    # if ($verbose) ;
	my %lg_names;
	my @filtered_sf;
	for my $sf (@simple_features) {
    for my $lg ( @{ $self->EFG_FEATURE_NAMES } ) {
      if ( $sf->analysis->logic_name =~m/^$lg$/) {
		    push @filtered_sf, $sf;
      }
	  }
	}
  print " after filtering by logic names got " . scalar(@filtered_sf) . " features left ... (if you do any filtering! ) \n" ;

	my $offset = $self->param('EXTEND_EFG_FEATURES');
	print "Extending start/end of EFG domains by +-$offset bp.\n";

	#Get fungen db adaptor
	my $efgdba        = $self->hrdb_get_con('regulation_db');
	my $fsa           = $efgdba->get_FeatureSetAdaptor;
	my $slice_adaptor = $efgdba->get_SliceAdaptor;
	my $input_id = $self->param('iid');
	my $slice = $slice_adaptor->fetch_by_name($input_id);

	for my $f (@simple_features) {
		my $f_start = 0;
		my $f_end   = 0;
		my @features = @{ $f->get_Features_by_Slice($slice) };
		for my $f_tmpt (@features) {
			# which I NEED TO PLAY WITH THIS a bit more.  This is hardcoded. I replace it with filter that is happening before. 
			if ( ( $f_tmpt->display_label =~ m/^H3K4me3/ ) or ( $f_tmpt->display_label =~ m/^H3K36me3/ ) ) {  
				$f_start = $f_tmpt->start - $offset;
				$f_end   = $f_tmpt->end + $offset;

				# print_feature($f_tmpt); # print all info of a feature
				my $biotype_to_use = 'efg';
				my $gene           = Bio::EnsEMBL::Gene->new(
					-analysis => $f_tmpt->analysis,
					-biotype  => $biotype_to_use,
					-strand   => $f_tmpt->strand,
				);
				$f_start = 1 if ( $f_start == 0 ); # There are fungen data that start from position 0. I slightly change that, since there is no possibility to have a (artificial) gene starting from 0!
				my $ex = new Bio::EnsEMBL::Exon(
					-START  => $f_start,
					-END    => $f_end,
					-STRAND => "1",    # if you ask for "$f_tmpt->strand", you will get strand 0 (non-applicable, don't know what it is!, so I change that too)
					-SLICE     => $f_tmpt->slice,
					-PHASE     => '0',
					-END_PHASE => '0',
					-ANALYSIS  => $f_tmpt->analysis,
				);
				my $transcript = Bio::EnsEMBL::Transcript->new(
					-analysis => $f_tmpt->analysis,
					-biotype  => $biotype_to_use,
					-strand   => $f_tmpt->strand,
				);

				$transcript->add_Exon($ex);
				$gene->add_Transcript($transcript);
				$gene->status(undef);
				$gene->description( $f_tmpt->display_label );
				$gene->dbID(undef);
				push @converted_sf, $gene;
			}
		}
	}

	print "--  "	. scalar(@simple_features) . " features sets in general ...\n";
	print "--  "	. scalar(@converted_sf) . " features/artificial genes created ...\n";
	return \@converted_sf;
}

=head2 score_cdna_overlap_with_multi_K36

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder
  Arg [2]   : Reference to an array of H3K36me3 gene objects
  Arg [3]   : Reference to an array of cDNA gene objects
  Function  : Selects for cDNA gene objects which overlap with multiple
              H3K36me3 features distributed across a genomic region.
              cDNA genes overlapping with multiple H3K36me3 features
              piling up in the same genomic region will not be 
              selected.
  Returntype: An array dbIDs (integers) of cDNAs which have been
              selected.
  Example   : $self->score_cdna_overlap_with_multi_K36(\@K36_genes, \@cdna_genes);

=cut

sub score_cdna_overlap_with_multi_K36 {
	my ( $self, $K36_genes, $cdna_genes ) = @_;
	print "-- artificial fungen genes:  " . scalar @$K36_genes 
		. "\n-- real genes: " . scalar @$cdna_genes . "\n";
	my ( $clustered, $unclustered ) = @{ simple_cluster_Genes( $K36_genes, "K36", $cdna_genes, "lincRNA" ) };
	my @K36_regions     = (@$K36_genes);
	my $K36_cluster_cnt = scalar(@K36_regions);

  # print "    Got a total of $K36_cluster_cnt clusters after local clustering, ". scalar@$clustered . " of which are K36 clusters and " . scalar@$unclustered . " are K36 singletons.\n"; # this is not really the case. 
	my @cdna_overlapping_multi_K36;
	my %scoring;

	if ( $K36_cluster_cnt > 1 )
	{ # i.e. the H3K36me3 feats don't fall into the same genomic location but distributed across exons
		foreach my $cdna (@$cdna_genes) {
			# my $curr_sr_start = $cdna->start;
			# my $curr_sr_end   = $cdna->end;
			foreach my $K36_region (@K36_regions) {
				# my $K36_clust_start = $curr_sr_start + $K36_region->start;
				# my $K36_clust_end   = $curr_sr_start + $K36_region->end;
				if ( $cdna->seq_region_end >= $K36_region->seq_region_start	&& $cdna->seq_region_start <= $K36_region->seq_region_end ) {    
					 # overlap
					 # if ($cdna->seq_region_end >= $K36_clust_start && $cdna->seq_region_start <= $K36_clust_end) {   # before overlap.... Since I am using the seq_region coordinates, I don't think we need t
					 # print "DEBUG_BK::RunnableDB::lincRNARegulation::score_cdna_overlap_with_multi_K36: cDNAstart: " .  $cdna->start . " cDNAend: " . $cdna->end . " K_start: " . $K36_region->start . " K_end: " . $K36_region->end . " cluster_start: " . $K36_clust_start . " clust_end: " . $K36_clust_end . "\n";
					 # print "DEBUG_BK::RunnableDB::lincRNARegulation::score_cdna_overlap_with_multi_K36: K36_r_start: " .  $K36_region->seq_region_start . " K36_r_end: " .   $K36_region->seq_region_end . "\n" ;
					 # print "    FOUND MATCH in this region for cDNA: " .$cdna->seq_region_start . " -  " . $cdna->seq_region_end . " (dbID " . $cdna->dbID.") " . " description: " . $K36_region->description . "\n";
					if ( $scoring{ $cdna->dbID } ) {
						$scoring{ $cdna->dbID }++;
					}
					else {
					  $scoring{ $cdna->dbID }++;
					  # print " "	. $scoring{ $cdna->dbID } . " "	. $cdna->dbID . "\n";
						$cdna->biotype('overlap_regulation');
						
            # this is because of ZMap. It wants the gene.analysis_id and transcript.analysis_id to be the same !
            # quick: mysql -h  genebuild11   -P 3306 -uXXXXXXX -pXXXXXX  kb15_human_fugn_reform_83_38 -e "UPDATE transcript t  JOIN gene g ON t.gene_id = g.gene_id  SET t.analysis_id = g.analysis_id ;"
					}
				}
			}
		}
	}
	else {
		print "This means only one block of H3K36me3 features are overlapping with cDNA(s).\n";
	}

	print "-- number of genes with fungen elements: " . scalar( keys %scoring ) . "\n";

	foreach my $cdna (@$cdna_genes) {
		if ( defined( $scoring{ $cdna->dbID } ) ) {
			my $biotype_tmp = $self->param_required('OUTPUT_BIOTYPE_OVERLAP');
			$cdna->biotype($biotype_tmp);
		}
		else {
			my $biotype_tmp = $self->param_required('OUTPUT_BIOTYPE_NOT_OVERLAP');
			$cdna->biotype($biotype_tmp);
		}
	}
	return $cdna_genes;
}



=head2 print_feature 

  Function  : print fungen feature for debug perpose
  Returntype: N/A

=cut

sub print_feature {
	my $feature = shift;
	print $feature->display_label . "\t("
		. $feature->seq_region_name . ":"
		. $feature->seq_region_start . "-"
		. $feature->seq_region_end . ")\n";
}




=head2 hive_set_config 

  Function  : loop through parameters and set them. All parameters of basic config (config_settings) will be accessible via param
  Returntype:  

=cut

sub hive_set_config {
	my $self = shift;

	# Throw is these aren't present as they should both be defined
	unless ( $self->param_is_defined('logic_name')
		&& $self->param_is_defined('module') )
	{
		warn("You must define 'logic_name' and 'module' in the parameters hash of your analysis in the pipeline config file, "
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
			warn("You have a key defined in the config_settings hash (in the analysis hash in the pipeline config) that does "
					. "not have a corresponding getter/setter subroutine. Either remove the key or add the getter/setter. Offending "
					. "key:\n"
					. $config_key );
		}
	}
}



=head2 

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut

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

sub OUTPUT_BIOTYPE_OVERLAP {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'OUTPUT_BIOTYPE_OVERLAP', $arg );
	}
	return $self->param('OUTPUT_BIOTYPE_OVERLAP');
}

sub OUTPUT_BIOTYPE_NOT_OVERLAP {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'OUTPUT_BIOTYPE_NOT_OVERLAP', $arg );
	}
	return $self->param('OUTPUT_BIOTYPE_NOT_OVERLAP');
}

sub BIOTYPE_TO_CHECK {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->param( 'BIOTYPE_TO_CHECK', $arg );
	}
	return $self->param('BIOTYPE_TO_CHECK');
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

sub create_analysis_object {
	my ( $self, $logic_name ) = @_;
	# need to create analysis object first
	my $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => $logic_name, );
	return $analysis;
}


=head2 memory_i_am_using 

  Function  : check memory for debug perpose
  Returntype: N/A

=cut

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
