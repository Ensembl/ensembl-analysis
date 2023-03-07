=head1 LICENSE

 Copyright [2021] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessGCA

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessGCA->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module dumps cds coord data from Ensembl dbs to file

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessGCA;

use warnings;
use strict;
use feature 'say';

use File::Spec::Functions;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name is_canonical_splice);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use File::Copy;
use DateTime;
use Bio::EnsEMBL::Utils::Exception qw (warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub param_defaults {
  my ($self) = @_;

  return {
#    %{$self->SUPER::param_defaults},
#    _branch_to_flow_to => 1,
#    use_generic_output_type => 0,
#    generic_output_type => 'cds',
  };
}

=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for cds coord dumping
 Returntype : None

=cut

sub fetch_input {
  my ($self) = @_;

  say "Fetching input";
  # For the combine files option we don't actually need to do anything in fetch input

  my $dirs_to_create    = [];
  my $current_genebuild = $self->param('current_genebuild');
  #my $current_genebuild  = 0;
  my $output_dir_base    = $self->param('base_output_dir');
  my $assembly_accession = $self->param('assembly_accession');
  my $repeatmodeler_library = $self->param('repeatmodeler_library');
  say "assembly_accession " . $assembly_accession;
  my $output_dir         = catdir( $output_dir_base, $assembly_accession );
  my $genome_files_dir   = catdir( $output_dir,      'genome_files' );
  my $short_read_dir     = catdir( $output_dir,      'short_read_fastq' );
  if ( $self->param('use_existing_short_read_dir') and -d $self->param('use_existing_short_read_dir') ) {
    $short_read_dir = $self->param('use_existing_short_read_dir');
  }

  #create the output dir with write permissions for group so that genebuild user can create output dirs here when running BRAKER
  my $od_result = system( 'mkdir -p -m775 ' . $output_dir );
    if ($od_result) {
      $self->throw( "Failed to create dir: " . $output_dir );
    }
  
  my $long_read_dir = catdir( $output_dir, 'long_read_fastq' );
  push( @$dirs_to_create, ($genome_files_dir, $short_read_dir, $long_read_dir ) );

  foreach my $dir (@$dirs_to_create) {
    my $result = system( 'mkdir -p ' . $dir );
    if ($result) {
      $self->throw( "Failed to create dir: " . $dir );
    }
  }

  my $registry_dba = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
    -host   => $self->param('registry_db')->{'-host'},
    -port   => $self->param('registry_db')->{'-port'},
    -user   => $self->param('registry_db')->{'-user'},
    -pass   => $self->param('registry_db')->{'-pass'},
    -dbname => $self->param('registry_db')->{'-dbname'} );

  my $taxonomy_adaptor = new Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor(
    -host   => 'mysql-ens-meta-prod-1',
    -port   => 4483,
    -user   => 'ensro',
    -dbname => 'ncbi_taxonomy' );

  my $pipeline_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $self->param('pipe_db')->{'-host'},
    -port   => $self->param('pipe_db')->{'-port'},
    -user   => $self->param('pipe_db')->{'-user'},
    -pass   => $self->param('pipe_db')->{'-pass'},
    -dbname => $self->param('pipe_db')->{'-dbname'} );
  my $init_file = $self->param('init_config');
  my $output_params = {};
  my $sth;
  my $general_hash = {};
  my ($stable_id_prefix, $clade, $species_taxon_id, $taxon_id, $assembly_name, $common_name, $assembly_refseq_accession, $assembly_date, $species_name, $assembly_group, $stable_id_start);
  if ( -e $init_file ) {
    open( IN, $init_file ) || throw("Could not open $init_file");
    while (<IN>) {
      my $line = $_;

      #$line =~ s/\s//g;
      if ( $line =~ /^(.+)\=(.+)/ ) {
        my $key   = $1;
        my $value = $2;
        say "Found key/value pair: " . $key . " => " . $value;
	if ( $key eq 'stable_id_prefix' ) {	
		#$general_hash->{$key} = $value;
	 $stable_id_prefix = $value;
	 
	 say "dentro if: " . $stable_id_prefix;
	 }
	if ( $key eq 'clade' ) {
	 $clade = $value;
	 }
	if ( $key eq 'species_taxon_id' ) { 
	 $species_taxon_id = $value;
	 }
	if ( $key eq 'taxon_id' ) { 
	 $taxon_id = $value;
 	}
	if ( $key eq 'assembly_name' ) {
	 $assembly_name = $value;
 	}
	if ( $key eq 'common_name' ) {	
       	 $common_name = $value;
        }
	if ( $key eq 'assembly_refseq_accession' ) {
	 $assembly_refseq_accession = $value;
        }
	if ( $key eq 'assembly_date' ) {
	 $assembly_date = $value;
 	}
	if ( $key eq 'species_name' ) {
	 $species_name = $value;
        }
	if ( $key eq 'assembly_group' ) {
	 $assembly_group = $value;
        }
	if ( $key eq 'stable_id_start' ) {
	 $stable_id_start = $value;
        }
	 #        if ( $key eq 'user_w' ) {
	#  $key = 'user';
	#  $general_hash->{$key} = $value;
	#}
	} elsif ( $line eq "\n" ) {
      } else {
        say "Line format not recognised. Skipping line:\n" . $line;
      }
    }
    close IN || throw("Could not close $init_file");

  } else {
    my $sql = "SELECT assembly_id FROM assembly WHERE CONCAT(chain,'.',version) = ?";
    my $sth = $registry_dba->dbc->prepare($sql);
    $sth->bind_param( 1, $assembly_accession );
    $sth->execute();
    my ($assembly_id) = $sth->fetchrow();

    $sql = "SELECT species_prefix,
                 clade,
                 species_id,
                 taxonomy,
                 assembly_name,
                 common_name,
                 refseq_accession,
                 assembly_date,
                 species_name,
                 pri_asm_group,
                 stable_id_space_start
                 FROM assembly JOIN meta as m USING(assembly_id) JOIN stable_id_space USING(stable_id_space_id) WHERE assembly_id=?";
    $sth = $registry_dba->dbc->prepare($sql);

    $sth->bind_param( 1, $assembly_id );
    $sth->execute();

    ( $stable_id_prefix, $clade, $species_taxon_id, $taxon_id, $assembly_name, $common_name, $assembly_refseq_accession, $assembly_date, $species_name, $assembly_group, $stable_id_start ) = $sth->fetchrow();
  }
  say "stable_id_prefix " . $stable_id_prefix;
  my $s = $stable_id_prefix;
  #$s =~ s/(ENS)/BRAKER/m;
  #my $stable_id_prefix = $s;
  $s =~ s/ENS//g;
  my $species_prefix  = uc($s);
  my $scientific_name = $species_name;
  say $species_name;
  $species_name = lc($species_name);
  say $species_name;
  $species_name =~ s/ +/\_/g;
  say $species_name;
  $species_name =~ s/\_$//;                # This cropped up
  say $species_name;
  $species_name =~ /([^\_]+)\_([^\_]+)/;
  say $species_name;
  my $p1                    = $1;
  my $p2                    = $2;
  my $binomial_species_name = $p1 . "_" . $p2;
  my $production_name       = $p1 . "_" . $p2;
  my $max_part_length       = 15;

  unless ( length($production_name) <= ( $max_part_length * 2 ) + 1 ) {
    my $ssp1 = substr( $p1, 0, $max_part_length );
    my $ssp2 = substr( $p2, 0, $max_part_length );
    $production_name = $ssp1 . "_" . $ssp2;
  }

  my $production_gca = $assembly_accession;
  $production_gca =~ s/\./v/;
  $production_gca =~ s/\_//g;
  $production_gca = lc($production_gca);
  $production_name .= "_" . $production_gca;

  # The assembly names for the alt haplotypes from DToL have spaces, probably better to substitute them in the registry
  # as opposed to here, but here is fine too. The paths to the files on the ftp site do not have an spaces so subbing with
  # underscore seems to be the correct thing to do based on the observed paths so far
  $assembly_name =~ s/ /\_/g;

  if ( $self->param('override_clade') ) {
    $clade = $self->param('override_clade');
    say "Clade param set in config, will override registry value";
  }

  my $clade_params = $self->get_clade_params($clade);

  #use taxon id to retrieve genus taxon id
  #genus taxon_id will be used to download genus level rnaseq data
  my $genus_taxon_id;
  my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
  my $taxon_node   = $node_adaptor->fetch_by_taxon_id($taxon_id);
  foreach my $ancestor ( @{ $node_adaptor->fetch_ancestors($taxon_node) } ) {
    #store genus level taxon id
    if ( $ancestor->rank eq 'genus' ) {
      $genus_taxon_id = $ancestor->taxon_id;
    }
    #store taxonomy at species level if species is a subspecies
#    elsif ($ancestor->rank eq 'species'){
#        $species_taxon_id = $ancestor->taxon_id;
#        $assembly_hash->{'species_taxon_id'} = $species_taxon_id;
#      }
    else { }
  }

  my $ensembl_release = $self->param('ensembl_release');

  my $core_db_details = $self->param('core_db');
  my $core_dbname     = $self->param('dbowner') . '_' . $production_gca . '_core_' . $ensembl_release . '_1';
  $core_db_details->{'-dbname'} = $core_dbname;

  my $otherfeatures_db_details = $self->param('otherfeatures_db');
  my $otherfeatures_dbname     = $self->param('dbowner') . '_' . $production_gca . '_otherfeatures_' . $ensembl_release . '_1';
  $otherfeatures_db_details->{'-dbname'} = $otherfeatures_dbname;

  my $rnaseq_summary_file    = catfile( $short_read_dir, $production_name . '.csv' );
  my $long_read_summary_file = catfile( $long_read_dir,  $production_name . '_long_read.csv' );

  my $toplevel_genome_file            = catfile( $output_dir, $species_name . "_toplevel.fa" );
  my $reheadered_toplevel_genome_file = catfile( $output_dir, $species_name . "_reheadered_toplevel.fa" );

  my $protein_file         = $clade_params->{'protein_file'};
  my $busco_protein_file   = $clade_params->{'busco_protein_file'};
  my $rfam_accessions_file = $clade_params->{'rfam_accessions_file'};
  my $busco_group          = $clade_params->{'busco_group'};
  my $max_intron_length    = $clade_params->{'max_intron_length'};

  # Meta details
  my $species_division = $clade_params->{'species_division'};
  my $species_url;
  my $species_display_name;
  my $species_strain       = "reference";
  my $species_strain_group = $production_name;
  my $strain_type          = "strain";
  if ( $assembly_name =~ /alternate_haplotype/ ) {
    $common_name = "alternate haplotype";
    $species_strain = "alternate haplotype";
  }

  if ( !$common_name ) {
    $common_name = "NA";
  }

  $species_display_name = $scientific_name . " (" . $common_name . ") - " . $assembly_accession;
  $species_url          = $scientific_name . "_" . $assembly_accession;
  $species_url =~ s/ /_/g;

  my $anno_commandline = ' --genome_file ' . $reheadered_toplevel_genome_file .
    ' --db_details ' . $core_db_details->{'-dbname'} . ',' .
    $core_db_details->{'-host'} . ',' .
    $core_db_details->{'-port'} . ',' .
    $core_db_details->{'-user'} . ',' .
    $core_db_details->{'-pass'} .
    ' --output_dir ' . $output_dir .
    ' --short_read_fastq_dir ' . $short_read_dir .
    ' --long_read_fastq_dir ' . $long_read_dir .
    ' --max_intron_length ' . $max_intron_length .
    ' --protein_file ' . $protein_file .
    ' --busco_protein_file ' . $busco_protein_file .
    ' --rfam_accessions_file ' . $rfam_accessions_file .
    ' --num_threads ' . $self->param('num_threads');

    $anno_commandline .= ' --repeatmasker_library ' . $repeatmodeler_library if ($repeatmodeler_library);

	$anno_commandline .= ' --run_full_annotation ' .
    ' --load_to_ensembl_db ';

  my $anno_red_commandline = ' --genome_file ' . $reheadered_toplevel_genome_file .
    ' --db_details ' . $core_db_details->{'-dbname'} . ',' .
    $core_db_details->{'-host'} . ',' .
    $core_db_details->{'-port'} . ',' .
    $core_db_details->{'-user'} . ',' .
    $core_db_details->{'-pass'} .
    ' --output_dir ' . $output_dir .
    ' --num_threads ' . $self->param('num_threads') .
    ' --run_masking ' .
    ' --run_repeats ' .
    ' --run_simple_features ' .
    ' --load_to_ensembl_db ';

  if ( $self->param('diamond_validation_db') ) {
    $anno_commandline .= ' --diamond_validation_db ' . $self->param('diamond_validation_db');
  }

  if ( $self->param('validation_type') ) {
    $anno_commandline .= ' --validation_type ' . $self->param('validation_type');
  }

  #Create a local copy of the registry and update the pipeline's resources with the new path
  my $new_registry_file = catfile( $output_dir_base, 'Databases.pm' );
  my $registry_file     = $self->param('registry_file');
  $registry_file =~ s|/+|/|g;

  $sth = $pipeline_db->dbc->prepare("update resource_description set worker_cmd_args=replace(worker_cmd_args, ?, ?);");
  $sth->bind_param( 1, "$registry_file" );
  $sth->bind_param( 2, "$new_registry_file" );
  unless ( $sth->execute() ) {
    throw("Could not update path for registry");
  }

  if ( -e $new_registry_file ) {
    $self->create_registry_entry( $new_registry_file, $core_db_details, $otherfeatures_db_details, $production_name );
  } else {
    system( 'cp ' . $registry_file . ' ' . $new_registry_file );
    $self->create_registry_entry( $new_registry_file, $core_db_details, $otherfeatures_db_details, $production_name );
  }

  #Check genebuild status of assembly
  if ( $current_genebuild == 1 ) {
    $self->update_annotation_status( $registry_dba, $assembly_accession, $current_genebuild );
  } else {
    #it will stop the pipeline if there is already an annotation in progress for this assembly
    $self->check_annotation_status( $registry_dba, $assembly_accession, $current_genebuild );
  }

  #Output
  $output_params->{'core_db'}                         = $core_db_details;
  $output_params->{'core_dbname'}                     = $core_dbname;
  $output_params->{'otherfeatures_db'}                = $otherfeatures_db_details;
  $output_params->{'otherfeatures_dbname'}            = $otherfeatures_dbname;
  $output_params->{'stable_id_start'}                 = $stable_id_start;
  $output_params->{'stable_id_prefix'}                = $stable_id_prefix;
  $output_params->{'clade'}                           = $clade;
  $output_params->{'species_taxon_id'}                = $species_taxon_id;
  $output_params->{'taxon_id'}                        = $taxon_id;
  $output_params->{'genus_taxon_id'}                  = $genus_taxon_id;
  $output_params->{'assembly_name'}                   = $assembly_name;
  $output_params->{'assembly_accession'}              = $assembly_accession;
  $output_params->{'assembly_refseq_accession'}       = $assembly_refseq_accession;
  $output_params->{'assembly_date'}                   = $assembly_date;
  $output_params->{'species_name'}                    = $species_name;
  $output_params->{'assembly_group'}                  = $assembly_group;
  $output_params->{'stable_id_start'}                 = $stable_id_start;
  $output_params->{'ensembl_release'}                 = $ensembl_release;
  $output_params->{'output_path'}                     = $output_dir;
  $output_params->{'genome_files_dir'}                = $genome_files_dir;
  $output_params->{'toplevel_genome_file'}            = $toplevel_genome_file;
  $output_params->{'reheadered_toplevel_genome_file'} = $reheadered_toplevel_genome_file;
  $output_params->{'short_read_dir'}                  = $short_read_dir;
  $output_params->{'long_read_dir'}                   = $long_read_dir;
  $output_params->{'species_url'}                     = $species_url;
  $output_params->{'species_division'}                = $species_division;
  $output_params->{'species_display_name'}            = $species_display_name;
  $output_params->{'species_strain'}                  = $species_strain;
  $output_params->{'species_strain_group'}            = $species_strain_group;
  $output_params->{'strain_type'}                     = $strain_type;
  $output_params->{'production_name'}                 = $production_name;
  $output_params->{'rnaseq_summary_file'}             = $rnaseq_summary_file;               # Problematic
  $output_params->{'long_read_summary_file'}          = $long_read_summary_file;            # Problematic
  $output_params->{'protein_file'}                    = $protein_file;
  $output_params->{'busco_protein_file'}              = $busco_protein_file;
  $output_params->{'busco_group'}                     = $busco_group;
  $output_params->{'rfam_accessions_file'}            = $rfam_accessions_file;
  $output_params->{'anno_commandline'}                = $anno_commandline;
  $output_params->{'anno_red_commandline'}            = $anno_red_commandline;
  $output_params->{'registry_file'}                   = $new_registry_file;
  $output_params->{'species_prefix'}                  = $species_prefix;
  $output_params->{'repeatmodeler_library'}           = $repeatmodeler_library;
  $self->param( 'output_params', $output_params );

}

sub run {
  my ($self) = @_;
}

=head2 write_output
p
 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $output_params = $self->param('output_params');
  $self->dataflow_output_id( $output_params, 1 );
}

sub get_clade_params {
  my ( $self, $clade ) = @_;

  my $clade_params = {};
  $clade_params->{'max_intron_length'} = 100000;

  if ( $clade eq 'lepidoptera' ) {
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/lepidoptera_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/lepidoptera_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_insect_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'lepidoptera_odb10',;
  } elsif ( $clade eq 'trichoptera' ) {
      $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/trichoptera_uniprot_proteins.fa',
	  $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/holometabola_orthodb11v0_proteins.fa',
	  $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_insect_ids.txt',
	  $clade_params->{'species_division'}     = 'EnsemblMetazoa',
	  $clade_params->{'busco_group'}          = 'endopterygota_odb10',;
  }  elsif ( $clade eq 'orthoptera' ) {
      $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/orthoptera_uniprot_proteins.fa',
          $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/hexapoda_orthodb11v0_proteins.fa',
          $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_insect_ids.txt',
          $clade_params->{'species_division'}     = 'EnsemblMetazoa',
          $clade_params->{'busco_group'}          = 'insecta_odb10',;
  } elsif ( $clade eq 'plasmodium' ) {
      $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/plasmodium_uniprot_proteins.fa',
          $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/plasmodium_orthodb_proteins.fa',
          $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_plasmodium_ids.txt',
          $clade_params->{'species_division'}     = 'EnsemblProtists',
          $clade_params->{'busco_group'}          = 'plasmodium_odb10',;
  }   elsif ( $clade eq 'teleostei' ) {
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/teleostei_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/actinopterygii_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_teleostei_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblVertebrates',
      $clade_params->{'busco_group'}          = 'actinopterygii_odb10',;
  } elsif ( $clade eq 'humans' ) {
    # Just temp stuff for testing
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/actinopterygii_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/actinopterygii_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_teleost_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblVertebrates',;
  } elsif ( $clade eq 'hymenoptera' ) {
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/hymenoptera_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/hymenoptera_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_insect_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'hymenoptera_odb10',;
  } elsif ( $clade eq 'diptera' ) {
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/diptera_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/diptera_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_insect_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'diptera_odb10',;
  } elsif ( $clade eq 'coleoptera' ) {
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/endopterygota_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/endopterygota_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_insect_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'endopterygota_odb10',;
  } elsif ( $clade eq 'crustacea' ) {
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/crustacea_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/crustacea_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_crustacea_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'arthropoda_odb10',;
  } elsif ( $clade eq 'mollusca' ) {
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/mollusca_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/mollusca_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_mollusc_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'mollusca_odb10',;
  } elsif ( $clade eq 'hemiptera' ) {
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/hemiptera_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/hemiptera_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_hemiptera_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'hemiptera_odb10',;
  } elsif ( $clade eq 'plants' ) {
    $clade_params->{'max_intron_length'} = 10000;
    $clade_params->{'protein_file'}      = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/viridiplantae_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/viridiplantae_orthodb_proteins_reheader.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_eudicotyledons_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblPlants',
      $clade_params->{'busco_group'}          = 'viridiplantae_odb10',;
  }  elsif ( $clade eq 'eudicotyledons' ) {
      $clade_params->{'max_intron_length'} = 10000;
      $clade_params->{'protein_file'}      = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/plants_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/eudicotyledons_orthodb11v0_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_eudicotyledons_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblPlants',
      $clade_params->{'busco_group'}          = 'eudicots_odb10',;
   } elsif ( $clade eq 'porifera' ) { #sponges
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/metazoa_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/metazoa_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_porifera_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'metazoa_odb10',;
  } elsif ( $clade eq 'viral' ) {
    # Test settings for loading viral dbs
    $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/worm_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/mollusca_orthodb_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_worm_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'metazoa_odb10',;
  }elsif ( $clade eq 'lophotrochozoa' ) {
      $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/lophotrochozoa_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/lophotrochozoa_orthodb11v0_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_lophotrochozoa_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'metazoa_odb10',;
  }elsif ( $clade eq 'cnidaria' ) {
      $clade_params->{'protein_file'} = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/cnidaria_uniprot_proteins.fa',
      $clade_params->{'busco_protein_file'}   = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/protein_sets/cnidaria_orthodb11v0_proteins.fa',
      $clade_params->{'rfam_accessions_file'} = '/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/clade_accessions/rfam_cnidaria_ids.txt',
      $clade_params->{'species_division'}     = 'EnsemblMetazoa',
      $clade_params->{'busco_group'}          = 'metazoa_odb10',;
  } else {
    $self->throw( 'Clade parameters not found for clade: ' . $clade );
  }

  return ($clade_params);
}

sub create_registry_entry {
  my ( $self, $registry_path, $core_db_details, $otherfeatures_db_details, $production_name ) = @_;

  unless ( -e $registry_path ) {
    $self->throw( "A registry file was not found on the path provided. Path:\n" . $registry_path );
  }

  my $core_string = "Bio::EnsEMBL::DBSQL::DBAdaptor->new(\n" .
    "-host => '" . $core_db_details->{'-host'} . "',\n" .
    "-port => '" . $core_db_details->{'-port'} . "',\n" .
    "-dbname => '" . $core_db_details->{'-dbname'} . "',\n" .
    "-user => '" . $core_db_details->{'-user'} . "',\n" .
    "-pass => '" . $core_db_details->{'-pass'} . "',\n" .
    "-species => '" . $production_name . "',\n" .
    "-group => 'core',\n" .
    ");\n";
  my $otherfeatures_string = "Bio::EnsEMBL::DBSQL::DBAdaptor->new(\n" .
    "-host => '" . $otherfeatures_db_details->{'-host'} . "',\n" .
    "-port => '" . $otherfeatures_db_details->{'-port'} . "',\n" .
    "-dbname => '" . $otherfeatures_db_details->{'-dbname'} . "',\n" .
    "-user => '" . $otherfeatures_db_details->{'-user'} . "',\n" .
    "-pass => '" . $otherfeatures_db_details->{'-pass'} . "',\n" .
    "-species => '" . $production_name . "',\n" .
    "-group => 'otherfeatures',\n" .
    ");\n";

  open( my $in, '<', $registry_path )
  	or die "Can't open " . $registry_path . " for reading.\n";
  my @lines = <$in>;
  close $in;

  open( my $out , ">" , $registry_path . ".tmp" )
  	or die "Can't open " . $registry_path . ".tmp for writing.\n";
  foreach my $line (@lines) {
    print $out $line;
    if ( $line =~ /\{/ ) {
      print $out $core_string;
      print $out $otherfeatures_string;
    }
  }
  close $out;

  my $result = system( 'mv ' . $registry_path . ".tmp " . $registry_path );
  if ($result) {
    $self->throw( "Issue overwriting the old registry with the new one. Registry path: " . $registry_path );
  }

}

=pod
=head1 Description of method
This method updates the registry database with the timestamp of when the annotation started. 
It also updates the registry with the status of the annotation as well as the user who started it.
=cut

sub update_annotation_status {
  my ( $self, $registry_dba, $accession, $current_genebuild ) = @_;
  my $dt          = DateTime->now;                                         # Stores current date and time as datetime object
  my $date        = $dt->ymd;
  my $assembly_id = $registry_dba->fetch_assembly_id_by_gca($accession);
  my ( $sql, $sth );
  if ( $current_genebuild == 1 ) {
    say "Updating genebuild status to overwrite";
    $sql = "update genebuild_status set is_current = ? where assembly_accession = ?";
    $sth = $registry_dba->dbc->prepare($sql);
    $sth->bind_param( 1, 0 );
    $sth->bind_param( 2, $accession );
    unless ( $sth->execute() ) {
      throw( "Could not update annotation status for assembly with accession " . $accession );
    }
  }
  $sql = "insert into genebuild_status(assembly_accession,progress_status,date_started,genebuilder,assembly_id,is_current,annotation_source) values(?,?,?,?,?,?,?)";
  $sth = $registry_dba->dbc->prepare($sql);
  $sth->bind_param( 1, $accession );
  $sth->bind_param( 2, 'in progress' );
  $sth->bind_param( 3, $date );
  $sth->bind_param( 4, $ENV{EHIVE_USER} || $ENV{USER} );
  $sth->bind_param( 5, $assembly_id );
  $sth->bind_param( 6, 1 );
  $sth->bind_param( 7, 'pending' );
  say "Accession being worked on is $accession";

  unless ( $sth->execute() ) {
    throw( "Could not update annoation status for assembly with accession " . $accession );
  }
}

=pod
=head1 Description of method
This method checks if there is an existing genebuild entry for the assembly.  
If yes, genebuilder must decide whether to continue annotation or not.
If genebuilder decides to continue, then rerun with option: -current_genebuild 1
This would automatically make this new genebuild the current annotation for tracking purposes
=cut

sub check_annotation_status {
  my ( $self, $registry_dba, $accession, $current_genebuild ) = @_;

  my @status = $registry_dba->fetch_genebuild_status_by_gca($accession);
  if (@status) {
    if ( $status[1] ) {
      throw( "A genebuild entry already exists for this assembly. " . "$accession\nGenebuild status: $status[0]\nDate started: $status[1]\nDate completed: $status[2]\nGenebuilder: $status[3]\nAnnotation source: $status[4]" . "\nTo proceed with this genebuild, re-run script with option: -current_genebuild 1" );
    }
    else {
      $self->update_annotation_status( $registry_dba, $assembly_accession, $current_genebuild );
    }
  }
}

1;
