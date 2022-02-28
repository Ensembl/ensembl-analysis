=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package HiveLincRNA_conf;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);

use Bio::EnsEMBL::ApiVersion qw/software_version/;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');


sub default_options {
  my ($self) = @_;

  return {
    %{ $self->SUPER::default_options() },

######################################################
#
# Variable settings- You change these!!!
#
######################################################

########################
# Pipe and ref db info
########################

# you will need to have a source with registry creation!
# check https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?pageId=39822808
    pipeline_name => '',

    out_dir => '',  # DO NOT FORGET to set your output dir.

    assembly_name => '',
    species => '', # your species name ie microcebus_murinus

    pipe_db_name => join('_', $self->o('dbowner'), $self->o('species'), $self->o('pipeline_name')),
    pipe_db_server => '', # NOTE! used to generate tokens in the resource_classes sub below
    dna_db_name => '', # what's your dna db name
    dna_db_server => '', # where is your dna db?  NOTE! used to generate tokens in the resource_classes sub below
    user => '',
    password => '',
    user_r => 'ensro',
    pass_r => undef,
    port => ,
    dna_db_port => ,

    cdna_db_host => '', # where is your RNAseq db? Where are your models?
    cdna_db_port => , # where is your RNAseq db? Where are your models?
    cdna_db_name => '', # what's the name of your RNAseq db? Where are your models?

    protein_coding_db_host => '', # where is your core db? Where are your models?
    protein_coding_db_port => , # where is your core db? Where are your models?
    protein_coding_db_name => '',

    # this is the output db. The results of all steps will be stored here!
    # !!! THIS OUTPUT DB NEEDS TO BE IN REGISTRY too !!! 
    lincRNA_db_host => '',
    lincRNA_db_port => ,
    lincRNA_db_name => $self->o('dbowner').'_'.$self->o('species').'_lincrna_3Gen_out_'.$self->o('release_number'),

    # this is for human regulation data
    regulation_db_host => '',
    regulation_db_name => '',

    # this is the output db after regulation data. The final results of regulation step will be stored here!
    regulation_debug_db_host => '',
    regulation_debug_db_name => '',

######################################################
#
# Mostly constant settings
#
######################################################

    'biotype_output' => 'rnaseq',
    'remove_duplicates_script_path' => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'find_and_remove_duplicates.pl'),

########################
# SPLIT PROTEOME File
########################
    'dump_biotypes_to_check_domain' => 'lincRNA_finder_2round', # this should be the same as the output biotype of lincRNA finder
    'max_seqs_per_file' => 20,
    'max_seq_length_per_file' => 20000, # Maximum sequence length in a file
    'max_files_per_directory' => 1000, # Maximum number of files in a directory
    'max_dirs_per_directory'  => $self->o('max_files_per_directory'),
    'file_translations' => catfile($self->o('out_dir'), 'hive_dump_translations.fasta'),

########################
# FINAL Checks parameters - Update biotypes to lincRNA, antisense, sense, problem ...
########################

    'file_for_length' => catfile($self->o('out_dir'), 'check_lincRNA_length.out'),  # list of genes that are smaller than 200bp, if any
    'file_for_biotypes' => catfile($self->o('out_dir'), 'check_lincRNA_need_to_update_biotype_antisense.out'), # mysql queries that will apply or not in your dataset (check update_database) and will update biotypes
     update_database => 'yes', # Do you want to apply the suggested biotypes? yes or no.
    'file_for_introns_support' => catfile($self->o('out_dir'), 'check_lincRNA_Introns_supporting_evidence.out'), # for debug

########################
# Interproscan
########################

    # species       => [],
    interproscan_version => '5.26-65.0',
    interproscan_exe     => 'interproscan.sh',
    run_interproscan     => 1,
    interproscan_applications => ['Pfam'], 

    # Delete rows in tables connected to the existing analysis
    linked_tables => ['protein_feature', 'object_xref'],

    # Retrieve analysis descriptions from the production database;
    # the supplied registry file will need the relevant server details.
    production_lookup => 1,

    # Remove existing analyses; if =0 then existing analyses
    # will remain, with the logic_name suffixed by '_bkp'.
    delete_existing => 1,
    
    run_seg=> 0, 
    
    analyses => [
      {
        'logic_name'    => 'pfam',
        'db'            => 'Pfam',
        'db_version'    => '31.0',
        'ipscan_name'   => 'Pfam',
        'ipscan_xml'    => 'PFAM',
        'ipscan_lookup' => 1,
      },
    ],


########################
# db info
########################
# NOTE! the dbname for each species is generated in the pipeline itself by setup_assembly_loading_pipeline
    'cdna_db' => {
      -host   => $self->o('cdna_db_host'),
      -port   => $self->o('cdna_db_port'),
      -user   => $self->o('user_r'),
      -dbname => $self->o('cdna_db_name'),
    },

    'source_protein_coding_db' => {
      -host   => $self->o('protein_coding_db_host'),
      -port   => $self->o('protein_coding_db_port'),
      -user   => $self->o('user_r'),
      -dbname => $self->o('protein_coding_db_name'),
    },

    'lincRNA_output_db' => {
      -host   => $self->o('lincRNA_db_host'),
      -port   => $self->o('lincRNA_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -dbname => $self->o('lincRNA_db_name'),
      -driver => $self->o('pipeline_db', '-driver'),
    },

    'regulation_db' => {
      -host   => $self->o('regulation_db_host'),
      -port   => $self->o('port'),
      -user   => $self->o('user_r'),
      -dbname => $self->o('regulation_db_name'),
    },

    'regulation_reform_db' => {
      -host   => $self->o('regulation_debug_db_host'),
      -port   => $self->o('port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -dbname => $self->o('regulation_debug_db_name'),
    },
  };
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

###############################################################################
# Set up before lincRNA ANALYSES
##############################################################################
    {
      -logic_name => 'create_lincRNA_output_db',
      -max_retry_count => 0,
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
                       source_db   => $self->o('dna_db'),
                       create_type => 'clone',
                       target_db => $self->o('lincRNA_output_db')
                     },
      -rc_name    => 'default',
      -input_ids => [{}],
      -flow_into => {
        '1->A' => ['AnalysisFactory'],
        'A->1' => ['create_toplevel_slices'],
      },
    },

    { -logic_name        => 'AnalysisFactory',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::ProteinFeatures::AnalysisFactory',
      -max_retry_count   => 1,
      -analysis_capacity => 20,
      -parameters        => {
                              analyses => $self->o('analyses'),
                              run_seg  => $self->o('run_seg'),
                            },
      -flow_into         => {
                              '2' => ['AnalysisSetup'],
                            },
      -meadow_type       => 'LOCAL',
    }, 

    {
      -logic_name        => 'AnalysisSetup',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::Common::AnalysisSetup',
      -max_retry_count   => 0,
      -parameters        => {
                              db_backup_required => 0,
                              db_backup_file     => catfile($self->o('out_dir'), '#species#', 'pre_pipeline_bkp.sql.gz'),
                              delete_existing    => $self->o('delete_existing'),
                              linked_tables      => $self->o('linked_tables'),
                              production_lookup  => $self->o('production_lookup'),
                              program            => $self->o('interproscan_version'),
                              species            => $self->o('species'),
                            },
      -meadow_type       => 'LOCAL',
    },
    

    {
      -logic_name => 'create_toplevel_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
                         target_db => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         feature_constraint => 1,
                         feature_type => 'gene',
                         feature_dbs => [$self->o('source_protein_coding_db'), $self->o('lincRNA_output_db')],
                       },
      -flow_into => {
                       '2' => ['Hive_LincRNARemoveDuplicateGenes'],
                      },
      -rc_name    => 'default',
    },

    {
      -logic_name => 'Hive_LincRNARemoveDuplicateGenes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveDuplicateGenes',
      -max_retry_count => 1,
      -batch_size    => 100,
      -parameters => {
                         logic_name => 'Hive_LincRNARemoveDuplicateGenes_t',
                         module => 'Hive_LincRNARemoveDuplicateGenes_1',
                         config_settings => $self->get_config_settings('HiveRemoveDuplication', 'Hive_Remove_dup_1'),
                         source_protein_coding_db => $self->o('source_protein_coding_db'),
                         lincRNA_output_db => $self->o('lincRNA_output_db'),
                         reference_db => $self->o('dna_db'),
                         cdna_db => $self->o('cdna_db'),
                         biotype_output => $self->o('biotype_output'),
                      },
      -rc_name    => 'normal_4600',
      -flow_into => {
                       '1->A' => ['Hive_LincRNAFinder'],
                       'A->1' => ['Hive_LincRNAEvaluator'],
      }
    },

###############################################################################
# LincRNA ANALYSES
##############################################################################

    {
      -logic_name => 'Hive_LincRNAFinder',
      -batch_size    => 100,
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAFinder',
      -max_retry_count => 1,
      -parameters => {
                       logic_name => 'Hive_LincRNAFinder',
                       module => 'HiveLincRNAFinder_1',
                       config_settings => $self->get_config_settings('HiveLincRNAFinder', 'Hive_lincRNA_1'),
                       source_protein_coding_db => $self->o('source_protein_coding_db'),
                       lincRNA_output_db => $self->o('lincRNA_output_db'),
                       reference_db => $self->o('dna_db'),
                    },
      -rc_name    => 'normal_4600',
      -flow_into => {
                      1 => ['HiveDumpTranslations'],
                     -1 => ['Hive_LincRNAFinder_himem'],
                    },
    },

    {
      -logic_name => 'Hive_LincRNAFinder_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAFinder',
      -batch_size => 5,
      -parameters => {
                       logic_name => 'Hive_LincRNAFinder',
                       module => 'HiveLincRNAFinder_1',
                       config_settings => $self->get_config_settings('HiveLincRNAFinder', 'Hive_lincRNA_1'),
                       source_cdna_db => $self->o('cdna_db'),
                       source_protein_coding_db => $self->o('source_protein_coding_db'),
                       lincRNA_output_db => $self->o('lincRNA_output_db'),
                       reference_db => $self->o('dna_db'),
                    },
      -rc_name    => 'normal_18000',
      -can_be_empty  => 1,
      -flow_into => {
                      1  => ['HiveDumpTranslations'],
                    },
    },


###############################################################################
## INTERPROSCAN
###############################################################################

    
    {
      -logic_name => 'HiveDumpTranslations',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpTranslations',
      -batch_size    => 100,
      -parameters => {
                        dump_translations_script => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'protein', 'dump_translations.pl'),
                        source_db => $self->o('lincRNA_output_db'),
                        dna_db => $self->o('dna_db'),
                        file => $self->o('file_translations'),
                        db_id => '1',
                        biotype => $self->o('dump_biotypes_to_check_domain'),
                     },
      -rc_name    => 'normal_1500_db',
      -flow_into => {
                     2 => ['SplitDumpFiles'],
                    },
     },

     {
       -logic_name => 'SplitDumpFiles',
       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplitFasta',
       -parameters => {
                        fasta_file              => $self->o('file_translations'),
                        out_dir                 => $self->o('out_dir'),
                        max_seqs_per_file       => $self->o('max_seqs_per_file'),
                        max_seq_length_per_file => $self->o('max_seq_length_per_file'),
                        max_files_per_directory => $self->o('max_files_per_directory'),
                        max_dirs_per_directory  => $self->o('max_dirs_per_directory'),
                       },
       -rc_name    => 'normal_1500',
       -flow_into     => {
                        2 => ['InterProScanLocal'],
       }
    },

    {
      -logic_name        => 'InterProScanLocal',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::ProteinFeatures::InterProScan',
      -hive_capacity     => 50,
      -batch_size        => 1,
      -max_retry_count   => 1,
      -can_be_empty      => 1,
      -parameters        =>
      {
        input_file                => '#split_file#',
        run_mode                  => 'local',
        interproscan_exe          => $self->o('interproscan_exe'),
        interproscan_applications => $self->o('interproscan_applications'),
        run_interproscan          => $self->o('run_interproscan'),
        species                   => $self->o('species'),
      },
      -rc_name           => 'normal_7900',
      -flow_into         => ['StoreProteinFeatures'],
    },

    {
      -logic_name        => 'StoreProteinFeatures',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::ProteinFeatures::StoreProteinFeatures',
      -analysis_capacity => 1,
      -batch_size        => 250,
      -max_retry_count   => 1,
      -parameters        => {
                              analyses        => $self->o('analyses'),
                              species         => $self->o('species'),
                            },
      -rc_name           => 'normal_1500_db',
    },


###############################################################################
### LAST PART - CHECK DOMAINs - UPDATE models - DELETE duplications
################################################################################
    {
      -logic_name   => 'Hive_LincRNAEvaluator',
      -module       => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincRNAEvaluator',
      -max_retry_count => 0,
      -batch_size   => 150,
      -parameters   => {
                         logic_name => 'Hive_LincRNAEvaluator_test',
                         module => 'HiveLincRNAEvaluator',
                         config_settings => $self->get_config_settings('HiveLincRNAEvaluator', 'Hive_lincRNAEvaluator_1'),
                         source_protein_coding_db => $self->o('source_protein_coding_db'),
                         lincRNA_output_db => $self->o('lincRNA_output_db'),
                         reference_db => $self->o('dna_db'),
                         source_cdna_db => $self->o('cdna_db'),
                      },
      -rc_name        => 'normal_4600', 
      -flow_into      => {
      	                 2 => ['delete_duplicate_genes'],
                      },
     },

    {
      -logic_name => 'delete_duplicate_genes',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => "perl ".$self->o('remove_duplicates_script_path')
                               ." -dbhost ".$self->o('lincRNA_output_db','-host')
                               ." -dbuser ".$self->o('user')
                               ." -dbpass ".$self->o('password')
                               ." -dbname ".$self->o('lincRNA_output_db','-dbname')
                               ." -dbport ".$self->o('lincRNA_output_db','-port')
                               ." -dnadbhost ".$self->o('dna_db','-host')
                               ." -dnadbuser ".$self->o('user_r')
                               ." -dnadbname ".$self->o('dna_db','-dbname')
                               ." -dnadbport ".$self->o('dna_db','-port'),
                     },
      -max_retry_count => 0,
      -rc_name => 'default',
      -input_ids => [{}],
      -flow_into      => ['Hive_LincRNAAftCheck_pi'],
       -wait_for   => ['Hive_LincRNAEvaluator'],
    },
    
    {
      -logic_name => 'Hive_LincRNAAftCheck_pi',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLincAfterChecks',
      -parameters => {
                         logic_name => 'HiveLincRNA_afterChecks',
                         module => 'HiveLincAfterChecks',
                         config_settings => $self->get_config_settings('HiveLincRNA_checks', 'HiveLincRNA_checks_1'),
                         lincRNA_output_db => $self->o('lincRNA_output_db'),
                         reference_db => $self->o('dna_db'),
                         source_cdna_db => $self->o('cdna_db'),
                         file_l => $self->o('file_for_length'),
                         file_is => $self->o('file_for_introns_support'),
                         file_b => $self->o('file_for_biotypes'),
                         assembly_name => $self->o('assembly_name'),
                         update_database => $self->o('update_database'),
                      },
      -rc_name    => 'normal_7900',
    },
  ];
}


# override the default method, to force an automatic loading of the registry in all workers
sub beekeeper_extra_cmdline_options {
    my $self = shift;
    return '-reg_conf '.$self->o('registry');
}


sub resource_classes {
      my $self = shift;
          return {
            'default' => { 'LSF' => $self->lsf_resource_builder('production-rh7',900,[$self->default_options->{pipe_db_server}]) },
            'normal_1500' => { 'LSF' => $self->lsf_resource_builder('production-rh7',1500,[$self->default_options->{pipe_db_server}]) },
            'normal_1500_db' => { 'LSF' => $self->lsf_resource_builder('production-rh7',1500,[$self->default_options->{pipe_db_server}, $self->default_options->{lincRNA_db_host}, $self->default_options->{dna_db_server}]) },
            'normal_4600' => { 'LSF' => $self->lsf_resource_builder('production-rh7',4600,[$self->default_options->{pipe_db_server}, $self->default_options->{lincRNA_db_host}, $self->default_options->{dna_db_server}, $self->default_options->{cdna_db_host}, $self->default_options->{protein_coding_host}])  },
            'normal_7900' => { 'LSF' => $self->lsf_resource_builder('production-rh7',7900,[$self->default_options->{pipe_db_server}, $self->default_options->{lincRNA_db_host}, $self->default_options->{dna_db_server}, $self->default_options->{cdna_db_host}, $self->default_options->{protein_coding_host}]) },
            'normal_18000' => { 'LSF' => $self->lsf_resource_builder('production-rh7',18000,[$self->default_options->{pipe_db_server}, $self->default_options->{lincRNA_db_host}, $self->default_options->{dna_db_server}, $self->default_options->{cdna_db_host}, $self->default_options->{protein_coding_host}]) },
                }
}


sub get_config_settings {

   # This is a helper sub created to access parameters that historically were held in separate configs in the
   # old pipeline. These are now stored in the master_config_settings sub below this one. In the analyses hashes
   # earlier in the config sets of these param can be requested and stored in the config_settings hash which
   # is them passed in as a parameter to the analysis. The converted analysis modules have code to take the key
   # value pairs from the config_settings hash and assign the values to the getter/setter sub associated with the
   # key.

   # Shift in the group name (a hash that has a collection of logic name hashes and a default hash)
   # Shift in the logic name of the specific analysis
   my $self = shift;
   my $config_group = shift;
   my $config_logic_name = shift;

   # And additional hash keys will be stored in here
   my @additional_configs = @_;

   # Return a ref to the master hash for the group using the group name
   my $config_group_hash = $self->master_config_settings($config_group);
   unless(defined($config_group_hash)) {
     die "You have asked for a group name in master_config_settings that doesn't exist. Group name:\n".$config_group;
   }
   # Final hash is the hash reference that gets returned. It is important to note that the keys added have
   # priority based on the call to this subroutine, with priority from left to right. Keys assigned to
   # $config_logic_name will have most priority, then keys in any additional hashes, then keys from the
   # default hash. A default hash key will never override a $config_logic_name key
   my $final_hash;

   # Add keys from the logic name hash
   my $config_logic_name_hash = $config_group_hash->{$config_logic_name};
   unless(defined($config_logic_name_hash)) {
     die "You have asked for a logic name hash that doesn't exist in the group you specified.\n".
         "Group name:\n".$config_group."\nLogic name:\n".$config_logic_name;
   }

   $final_hash = $self->add_keys($config_logic_name_hash,$final_hash);

   # Add keys from any additional hashes passed in, keys that are already present will not be overriden
   foreach my $additional_hash (@additional_configs) {
     my $config_additional_hash = $config_group_hash->{$additional_hash};
     $final_hash = $self->add_keys($config_additional_hash,$final_hash);
   }

   # Default is always loaded and has the lowest key value priority
   my $config_default_hash = $config_group_hash->{'Default'};
   $final_hash = $self->add_keys($config_default_hash,$final_hash);

   return($final_hash);
}

sub add_keys {
  my ($self,$hash_to_add,$final_hash) = @_;

  foreach my $key (keys(%$hash_to_add)) {
    unless(exists($final_hash->{$key})) {
      $final_hash->{$key} = $hash_to_add->{$key};
    }
  }

  return($final_hash);
}

sub master_config_settings {

  my ($self,$config_group) = @_;
  my $master_config_settings = {

  HiveRemoveDuplication => {
  	Hive_Remove_dup_1 => {
                  # WHERE AND WHAT TO FETCH FOR lincRNA CANDIDATES AND VALIDATION GENES
                  #----------------------------------------------------------------------

                  RNA_DB => {
                                   # Format: database_alias => ['biotype_1', 'biotype_2']
                                   'cdna_db' => ['fetch_all_biotypes'],  # fetch_all_biotypes for all biotypes
                                   # 'cdna_db' => ['rnaseq2_use'],
                                },
  	},
  },

  HiveLincRNAFinder => {
    Hive_lincRNA_1 => {

                  # cDNA / RNASeq model / protein-coding gene databases
                  # ---------------------------------------------------

                  # SOURCE_CDNA_UPDATE_DB is the database holding cDNAs or transcript models (e.g. built from Illumina
                  # RNASeq data) which will be considered as potential lincRNAs

                  # SOURCE_PROTEIN_CODING_DB is the database holding the latest protein_coding genes for the analysed
                  # species.  This database usually is the latest Ensembl core DB available.
                  # to fix biotypes: DB:genebuild11,kb15_human_rna_83_38_linc> update gene set biotype = "human_rnaseq" where biotype = "protein_coding";
                  NEW_SET_1_CDNA => {
                                        'lincRNA_output_db'  => ['fetch_all_biotypes'],
                                        #  'cdna_db'    => ['rnaseq'],
                                    },

                  NEW_SET_2_PROT  => {
                                        'source_protein_coding_db' => ['fetch_all_biotypes'],
                                     },

                  # EFG_FEATURE_DB    => 'regulation_db',
                  # EFG_FEATURE_NAMES => ['H3K4me3','H3K36me3'],
                  # EXTEND_EFG_FEATURES => 350,

                  FIND_SINGLE_EXON_LINCRNA_CANDIDATES => 1000, # I don't want single exon candidates!
                  CDNA_CODING_GENE_CLUSTER_IGNORE_STRAND => 1,
                  MAXIMUM_TRANSLATION_LENGTH_RATIO => 99,
                  MAX_TRANSLATIONS_PER_GENE => 20,

                  OUTPUT_DB => 'lincRNA_output_db',
                  OUTPUT_BIOTYPE => 'lincRNA_finder_2round',

                  WRITE_DEBUG_OUTPUT => 0,     # Set this to "0" to turn off debugging OR to "1000" to set it on.
                  DEBUG_OUTPUT_DB    => 'lincRNA_output_db',    # where debug output (if any) will be written to

                  # logic_name for H3 features which do NOT cluster with cDNA/RNASeq models or protein_coding genes:
                  # DEBUG_LG_EFG_UNCLUSTERED  => 'efg_NO_cdna_update_NO_pc',

                  # logic_name for H3 features which cluster with cDNA/RNASeq models checked to be not overlapping with protein_coding genes:
                  # DEBUG_LG_EFG_CLUSTERING_WITH_CDNA => 'efg_cluster_non_pc_cdna_update',

    },
    Hive_lincRNA_2 => {

    },
  },
  HiveLincRNAEvaluator => {
  	Hive_lincRNAEvaluator_1 => {
                  # WHERE AND WHAT TO FETCH FOR lincRNA CANDIDATES AND VALIDATION GENES
                  #----------------------------------------------------------------------


                  LINCRNA_DB => {
                                   # Format: database_alias => ['biotype_1', 'biotype_2']

                                   # Specify biotypes of genes to be fetched out of the lincRNAFinder output database.
                                   # You have to fetch *all* lincRNA genes, i.e. there is no need to pre-filter the lincRNAFinder output
                                   # to remove lincRNA candidates which contain protein features (pfam/tigfam), because the Evaluator code
                                   # will do the filtering.

                                   # REFERENCE_DB => ['lincRNA_finder_1'],
                                   lincRNA_output_db => ['lincRNA_finder_2round'],
                                   # LINCRNA_OUTPUT_DB => ['fetch_all_biotypes'],  # ,'prot_feat'],
                                },

                  VALIDATION_DBS => {

                                   # Format: database_alias_1 => ['biotype_1', 'biotype_2'],
                                   #         databass_alias_2 => ['biotype_3', 'biotype_4'],

                                   # lincRNAs are "validated" if they do not overlap with Ensembl genes.

                                   # There is usually just one validation DB, and it will be the one with alias "SOURCE_PROTEIN_CODING_DB"
                                   # (containing Ensembl core genes), which was also used in the lincRNAFinder stage.
                                   # Specify the biotype of Ensembl genes to be fetched here.  If all genes are to be fetched from the
                                   # validation DB, put down "fetch_all_biotypes" as the biotype, for it is an alias accepted by BaseGeneBuild
                                   # RunnableDB.

                                   # SOURCE_PROTEIN_CODING_DB => ['rRNA'],
                                   source_protein_coding_db => ['protein_coding'],
                                 },


                  # OPTIONAL EVALUATION CRITERIA
                  #------------------------------

                  # If your lincRNAFinder output contains single-exon lincRNA candidates and you want to remove
                  # them from analysis after all, turn on 'EXCLUDE_SINGLE_EXON_LINCRNAS'

                  EXCLUDE_SINGLE_EXON_LINCRNAS => 1000,

                  # Some single-exon lincRNA candidates are disguised as two-exon ones containing short frameshift
                  # introns.  You can remove these candidates by turning on 'EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS'.
                  # Any introns of length longer than 'MAX_FRAMESHIFT_INTRON_LEN' will be regarded as non-frameshift/sane.

                  MAX_FRAMESHIFT_INTRON_LEN => 9,
                  EXCLUDE_ARTEFACT_TWO_EXON_LINCRNAS => 0,


                  # OUTPUT SETTINGS
                  #-----------------

                  # At the end of the lincRNA pipeline, there is a GeneBuilder step where all lincRNA "genes" (one-gene-
                  # one-transcript) are collapsed into proper Ensembl genes. You can control how many alternative
                  # spliceforms are allowed by changing "MAX_TRANSCRIPTS_PER_CLUSTER". The default value is 3, which is
                  # considerably lower than what's used for a conventional genebuild (set at 15).

                  MAX_TRANSCRIPTS_PER_CLUSTER => 3,

                  FINAL_OUTPUT_BIOTYPE => "lincRNA_pass_Eval_no_pfam",
                  FINAL_OUTPUT_DB      => 'lincRNA_output_db',


                  # DEALING WITH lincRNAs WHICH OVERLAP WITH EXISTING PROCESSED_TRANSCRIPT GENES
                  # ------------------------------------------------------------------------------

                  # Possible actions are:
                  # (1) update the *gene* analysis AND biotype of all processed_transcripts which cluster with newly-identified lincRNAs.
                  #     This is useful if you want to copy your lincRNAs into a 'ready-to-go' core DB where all other types of genes
                  #     are properly annotated (i.e. gene and transcript biotypes and anlaysis logic_names standardised for public release).
                  # AND/OR
                  # (2) write the newly-identified lincRNAs into FINAL_OUTPUT_DB even when they overlap with processed_transcripts
                  #     This is useful for debugging purposes.
                  #
                  # For (1):
                  #
                  # Turn on 'MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB'.  The overlapped proc_trans genes are expected to be in
                  # what's defined by 'UPDATE_SOURCE_DB'.
                  #
                  # Before the update is performed, the code checks that the sanity of processed_transcript genes and requires they
                  # have analysis logic_names containing the string specified in 'PROC_TRANS_HAVANA_LOGIC_NAME_STRING' (usually "havana")
                  # because processed_transcript genes should only be from "havana" or "ensembl_havana_gene" analyses.

                  # The analysis logic_name will be updated to what's defined in 'OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME'.
                  # Make sure this analysis already exists in the DB defined by UPDATE_SOURCE_DB.
                  # The biotype will be updated to 'proc_trans_turned_lincRNA (hard-coded)'.

                  MARK_OVERLAPPED_PROC_TRANS_IN_VALIDATION_DB => 10000, # no validation db check for now. if you say yes here, you have to change the following parameters about validation
                  PROC_TRANS_HAVANA_LOGIC_NAME_STRING => 'havana',
                  OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME => 'ensembl_havana_gene',

                  UPDATE_SOURCE_DB => 'lincRNA_output_db',

                  # For (2):
                  #
                  # Turn on 'WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS'.
                  #
                  # The lincRNAs written will have analysis logic_name just like any other lincRNAs generated by lincRNAEvaluator.
                  # The biotype will be hard-coded as 'lincRNA_clusters_with_proc_trans'.

                  WRITE_LINCRNAS_WHICH_CLUSTER_WITH_PROC_TRANS => 1,


                  # DEALING WITH lincRNAs WHICH OVERLAP WITH EXISTING lincRNA GENES
                  # ----------------------------------------------------------------

                  # Sometimes the VALIDATION DB contains existing lincRNAs from previous work.
                  # The lincRNAs identified in this analysis may overlap with existing lincRNAs. (They do not have to be
                  # identical.)

                  # Possible actions are:
                  # (1) update the *gene* analysis AND biotype of all existing lincRNA genes which cluster with newly-identified lincRNAs;
                  # AND/OR
                  # (2) write the newly-identified lincRNAs into FINAL_OUTPUT_DB even when they overlap with existing lincRNAs

                  # For (1):
                  #
                  # Turn on 'MARK_EXISTING_LINCRNA_IN_VALIDATION_DB'. The overlapped existing lincRNA genes are expected to be in
                  # 'UPDATE_SOURCE_DB' (shared with proc_tran settings).
                  # The analysis logic_name will be updated to what's defined in "OVERLAPPED_GENES_NEW_MERGED_LOGIC_NAME" above.
                  # (Make sure this analysis already exists in the DB defined by "UPDATE_SOURCE_DB".)
                  # The biotype will be changed to "lincRNA_common" (hard-coded).

                  MARK_EXISTING_LINCRNA_IN_VALIDATION_DB => 1,

                  # For (2):
                  #
                  # Turn on 'WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS'
                  # The lincRNAs written will have analysis logic_name just like any other lincRNAs generated by lincRNAEvaluator.
                  # The biotype will be hard-coded as "ncRNA_clusters_with_existing_lincRNA".

                  WRITE_LINCRNAS_WHICH_CLUSTER_WITH_EXISTING_LINCRNAS => 1,


                  # DEBUG OUTPUT WITH REJECTED lincRNAS
                  # ------------------------------------
                  # Turning on this option will write three categories of rejected lincRNAs into the FINAL_OUTPUT_DB:
                  #
                  #   (i) lincRNA candidates which are somehow assoicated with protein domains. They will be given biotype suffix "_reject_prot_dom";
                  #  (ii) lincRNAs which don't associate with protein domains but overlap with Ensembl gene of single biotype. They will be
                  #       given biotype suffix "_reject_single";
                  # (iii) lincRNAs which don't associate with protein domains but overlap with Ensembl genes of multiple biotypes. They will
                  #       be given biotype suffix "_reject_mult".

                  WRITE_REJECTED_NCRNAS => 1,


  	},
  },



  # this could be used only for human and mouse.
  HiveLincRNAReg=> {
    HiveLincRNAReg_1 => {

                  EFG_FEATURE_DB    => 'regulation_db',
                  EFG_FEATURE_NAMES => ['H3K4me3','H3K36me3'],
                  OUTPUT_DB => $self->default_options->{lincRNA_output_db},
                  BIOTYPE_TO_CHECK => 'lincRNA_pass_Eval_no_pfam',

                  EXTEND_EFG_FEATURES => 500, # make efg features longer!

                  CHECK_CDNA_OVERLAP_WITH_BOTH_K4_K36 => 1,

                  CHECK_CDNA_OVERLAP_WITH_MULTI_K36 => 1,

                  WRITE_DEBUG_OUTPUT => 1,     # Set this to "0" to turn off debugging.
                  DEBUG_OUTPUT_DB    => 'regulation_reform_db',    # where debug output (if any) will be written to
                  OUTPUT_BIOTYPE_OVERLAP => 'lincRNA_withReg',
                  OUTPUT_BIOTYPE_NOT_OVERLAP => 'lincRNA_withOutReg',

                  # logic_name for H3 features which do NOT cluster with cDNA/RNASeq models or protein_coding genes:
                  DEBUG_LG_EFG_UNCLUSTERED  => 'efg_NO_cdna_update_NO_pc',

                  # logic_name for H3 features which cluster with cDNA/RNASeq models checked to be not overlapping with protein_coding genes:
                  DEBUG_LG_EFG_CLUSTERING_WITH_CDNA => 'efg_cluster_non_pc_cdna_update',
    },
  },





  HiveLincRNA_checks=> {
    HiveLincRNA_checks_1 => {
    	            Final_BIOTYPE_TO_CHECK  => 'lincRNA_pass_Eval_no_pfam',
    }
  }


  };

  return($master_config_settings->{$config_group});

}

1;
