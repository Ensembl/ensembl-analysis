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
  }
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

  my $dirs_to_create = [];
  my $current_genebuild = $self->param('current_genebuild');
  my $output_dir_base = $self->param('base_output_dir');
  my $assembly_accession = $self->param('assembly_accession');
  my $output_dir = catdir($output_dir_base,$assembly_accession);
  my $genome_files_dir = catdir($output_dir,'genome_files');
  push(@$dirs_to_create,($output_dir,$genome_files_dir));

  foreach my $dir (@$dirs_to_create) {
    my $result = system('mkdir -p '.$dir);
    if($result) {
      $self->throw("Failed to create dir: ".$dir);
    }
  }

  my $registry_dba = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
  -host    => $self->param('registry_db')->{'-host'},
  -port    => $self->param('registry_db')->{'-port'},
  -user    => $self->param('registry_db')->{'-user'},
  -pass   => $self->param('registry_db')->{'-pass'},
  -dbname  => $self->param('registry_db')->{'-dbname'});


  my $taxonomy_adaptor = new Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor(
  -host    => 'mysql-ens-meta-prod-1',
  -port    => 4483,
  -user    => 'ensro',
  -dbname  => 'ncbi_taxonomy');


  my $sql = "SELECT assembly_id FROM assembly WHERE CONCAT(chain,'.',version) = ?";
  my $sth = $registry_dba->dbc->prepare($sql);
  $sth->bind_param(1,$assembly_accession);
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

  $sth->bind_param(1,$assembly_id);
  $sth->execute();
  
  my $output_params = {};
  my ($stable_id_prefix,$clade,$species_taxon_id,$taxon_id,$assembly_name,$common_name,$assembly_refseq_accession,$assembly_date,$species_name,$assembly_group,$stable_id_start) = $sth->fetchrow();
  my $scientific_name = $species_name;
  $species_name = lc($species_name);
  $species_name =~ s/ +/\_/g;
  $species_name =~ s/\_$//; # This cropped up
  $species_name =~ /([^\_]+)\_([^\_]+)/;
  my $p1 = $1;
  my $p2 = $2;
  my $binomial_species_name = $p1."_".$p2;
  my $production_name = $p1."_".$p2;
  my $max_part_length = 15;
  unless(length($production_name) <= ($max_part_length * 2) + 1) {
  
    my $ssp1 = substr($p1,0,$max_part_length);
    my $ssp2 = substr($p2,0,$max_part_length);
    $production_name = $ssp1."_".$ssp2;
  }
  
  my $production_gca = $assembly_accession;
  $production_gca =~ s/\./v/;
  $production_gca =~ s/\_//g;
  $production_gca = lc($production_gca);
  $production_name .= "_".$production_gca;
  
  # The assembly names for the alt haplotypes from DToL have spaces, probably better to substitute them in the registry
  # as opposed to here, but here is fine too. The paths to the files on the ftp site do not have an spaces so subbing with
  # underscore seems to be the correct thing to do based on the observed paths so far
  $assembly_name =~ s/ /\_/g;


  #use taxon id to retrieve genus taxon id
  #genus taxon_id will be used to download genus level rnaseq data
  my $genus_taxon_id;
  my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
  my $taxon_node = $node_adaptor->fetch_by_taxon_id($taxon_id);
  foreach my $ancestor ( @{ $node_adaptor->fetch_ancestors($taxon_node)}){
    #store genus level taxon id
    if ($ancestor->rank eq 'genus'){
      $genus_taxon_id = $ancestor->taxon_id;
    }
    #store taxonomy at species level if species is a subspecies
#    elsif ($ancestor->rank eq 'species'){
#        $species_taxon_id = $ancestor->taxon_id;
#        $assembly_hash->{'species_taxon_id'} = $species_taxon_id;
#      }
     else{}
  }


  my $ensembl_release = $self->param('ensembl_release');

  my $core_db_details = $self->param('core_db');
  my $core_dbname = $production_name.'_core_'.$ensembl_release.'_1';
  if ($self->param_is_defined('dbowner')) {
    $core_dbname = $self->param('dbowner')."_$core_dbname";
  }
  $core_db_details->{'-dbname'} = $core_dbname;


  my $toplevel_genome_file = catfile($output_dir,$species_name."_toplevel.fa");
  my $reheadered_toplevel_genome_file = catfile($output_dir,$species_name."_reheadered_toplevel.fa");


  # Meta details
  my $species_url;
  my $species_display_name;
  if ($assembly_name =~ /alternate_haplotype/) {
    $common_name = "alternate haplotype";
  }

  if(!$common_name) {
    $common_name = "NA";
  }

  $species_display_name = $scientific_name." (".$common_name.") - ".$assembly_accession;
  $species_url = $scientific_name."_".$assembly_accession;
  $species_url =~ s/ /_/g;
  
  # Note this should probably be update so that haps are strain assemblies under a strain group
  # The group should be changed to cut off the GCA from the production name, then the type can
  # be set to alternate haplotype. This needs to be discussed before implementing
  my $species_strain_group = $production_name;
  my $strain_type = "strain";

  my $anno_repeats_commandline = ' --genome_file ' . $reheadered_toplevel_genome_file .
    ' --db_details ' . $core_db_details->{'-dbname'} . ',' .
    $core_db_details->{'-host'} . ',' .
    $core_db_details->{'-port'} . ',' .
    $core_db_details->{'-user'} . ',' .
    $core_db_details->{'-pass'} .
    ' --output_dir ' . $output_dir .
    ' --num_threads ' . $self->param('num_threads') .
    ' --run_repeatmasker ' .
    ' --run_repeats ' .
    ' --load_to_ensembl_db ';

  #Check genebuild status of assembly unless custom loading no need to check registry
  if ( $current_genebuild == 1 ) {
    $self->update_annotation_status( $registry_dba, $assembly_accession, $current_genebuild );
  }
  else {
    #it will stop the pipeline if there is already an annotation in progress for this assembly
    say "UPDATE REGISTRY";
    $self->check_annotation_status( $registry_dba, $assembly_accession, $current_genebuild );
  }  
  $output_params->{'anno_repeats_commandline'} = $anno_repeats_commandline;
  $output_params->{'core_db'} = $core_db_details;
  $output_params->{'core_dbname'} = $core_dbname;
  $output_params->{'stable_id_start'} = $stable_id_start;
  $output_params->{'stable_id_prefix'} = $stable_id_prefix;
  $output_params->{'clade'} = $clade;
  $output_params->{'species_taxon_id'} = $species_taxon_id;
  $output_params->{'taxon_id'} = $taxon_id;
  $output_params->{'genus_taxon_id'} = $genus_taxon_id;
  $output_params->{'assembly_name'} = $assembly_name;
  $output_params->{'assembly_accession'} = $assembly_accession;
  $output_params->{'assembly_refseq_accession'} = $assembly_refseq_accession;
  $output_params->{'assembly_date'} = $assembly_date;
  $output_params->{'species_name'} = $species_name;
  $output_params->{'assembly_group'} = $assembly_group;
  $output_params->{'stable_id_start'} = $stable_id_start;
  $output_params->{'ensembl_release'} = $ensembl_release;
  $output_params->{'output_path'} = $output_dir;
  $output_params->{'genome_files_dir'} = $genome_files_dir;
  $output_params->{'toplevel_genome_file'} = $toplevel_genome_file;
  $output_params->{'reheadered_toplevel_genome_file'} = $reheadered_toplevel_genome_file;
  $output_params->{'species_url'} = $species_url;
  $output_params->{'species_display_name'} = $species_display_name;
  $output_params->{'species_strain_group'} = $species_strain_group;
  $output_params->{'strain_type'} = $strain_type;
  $output_params->{'production_name'} = $production_name;
  $self->param('output_params',$output_params);

  $self->create_registry_entry($self->param('registry_file'),$core_db_details,$production_name);
}




sub run {
  my ($self) = @_;
}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $output_params = $self->param('output_params');
  $self->dataflow_output_id($output_params,1);
}


sub create_registry_entry {
  my ($self,$registry_path,$core_db_details,$production_name) = @_;

  unless(-e $registry_path) {
    $self->throw("A registry file was not found on the path provided. Path:\n".$registry_path);
  }

  my $core_string = "Bio::EnsEMBL::DBSQL::DBAdaptor->new(\n".
                      "-host => '".$core_db_details->{'-host'}."',\n".
                      "-port => '".$core_db_details->{'-port'}."',\n".
                      "-dbname => '".$core_db_details->{'-dbname'}."',\n".
                      "-user => '".$core_db_details->{'-user'}."',\n".
                      "-pass => '".$core_db_details->{'-pass'}."',\n".
                      "-species => '".$production_name."',\n".
                      "-group => 'core',\n".
                      ");\n";
  open(IN,$registry_path);
  my @lines = <IN>;
  close IN;

  open(OUT,">".$registry_path.".tmp");
  foreach my $line (@lines) {
    print OUT $line;
    if($line =~ /\{/) {
      print OUT $core_string;
    }
  }
  close OUT;

  my $result = system('mv '.$registry_path.".tmp ".$registry_path);
  if($result) {
    $self->throw("Issue overwriting the old registry with the new one. Registry path: ".$registry_path);
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
  say "Creating new assembly annotation status in registry...\n";
  $sql = "insert into genebuild_status(assembly_accession,progress_status,date_started,genebuilder,assembly_id,is_current,annotation_source) values(?,?,?,?,?,?,?)";
  $sth = $registry_dba->dbc->prepare($sql);
  $sth->bind_param( 1, $accession );
  $sth->bind_param( 2, 'in progress' );
  $sth->bind_param( 3, $date );
  $sth->bind_param( 4, $ENV{EHIVE_USER} || $ENV{USER} );
  $sth->bind_param( 5, $assembly_id );
  $sth->bind_param( 6, 1 );
  $sth->bind_param( 7, 'pending' );
  say "SQL Successful. Accession being worked on is $accession";

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
      throw( "A genebuild entry already exists for this assembly. " . "$accession\nGenebuild status: $status[0]\nDate started: $status[1]\nDate completed: $status[2]\nGenebuilder: $status[3]\nAnnotation source: $status[4]" . "\nTo proceed with this genebuild, re-run script with option: -current_genebuild 1" );
    }
    else {
      print "Attempting to update annotation status on $accession accession\n";
      $self->update_annotation_status( $registry_dba, $accession, $current_genebuild ); 
    }
}

1;
