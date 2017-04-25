=head1 LICENSE

Copyright [2017] EMBL-European Bioinformatics Institute

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

package UniProtAlignment_conf;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions qw(catdir catfile);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');
use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

###########################################
# User parameters: you fill in stuff here #
###########################################

# Path to the fasta file containing the genome sequence
'genome_file'  => $ENV{PWD}.'/data/assembly/genome_dumps/rattus_norvegicus_softmasked_toplevel.fa', # Fill in!!!!

# Path to the fasta file containing the proteins sequences
'protein_file' => $ENV{PWD}.'/data/uniprot_proteins/uniprot_proteins.fa', # Fill in!!!!

# The logic_name of the repeatmasker analysis in the core db
'repeat_masking_logic_names' => ['repeatmasker_repbase_rat'], # Fill in!!!!

# The path to the ensembl-api-folder
    'enscode_root_dir'          => $ENV{PWD}.'/../../../..', #!!!!!!!!!!! git repo checkouts

# Database settings
'pipe_dbname'             => 'test_workshop_ex2_pipe', # Fill in!!!!
'dna_dbname'              => 'test_workshop_rat_core', # Fill in!!!!
'core_dbname'             => $self->o('dna_dbname'), # Fill in!!!!
'uniprot_dbname'   => 'test_workshop_rat_uniprot', # Fill in!!!!
'user_r'                  => $ENV{EHIVE_ROUSER}, # Fill in!!!!

# Exonerate settings
'software_base_path' => '' || $ENV{LINUXBREW_HOME},
'binary_base' => catdir($self->o('software_base_path'), 'bin'),
'exonerate_path'         => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
'exonerate_pid'           => '50', # Fill in!!!!
'exonerate_cov'           => '50', # Fill in!!!!


#######################################################################
# Pre-filled parameters: you don't need to change anything below here #
#######################################################################
'create_type'             => 'clone',
'protein_table_name' => 'protein_sequences',
'exonerate_calculate_coverage_and_pid' => '1',
'exonerate_region_padding'             => 0,
'protein_batch_size' => 1,

#########################################################################
# DB connection info: these databases are used in the pipeline analyses #
#########################################################################

    'pipe_db_server'               => $self->o('host'),
    'pipe_db_port'                  => $self->o('port'), #!!!!!!!!!!!
    'dna_db_server'                => $self->o('host'), #!!!!!!!!!!!
    'dna_db_port'                  => $self->o('port'), #!!!!!!!!!!!
    'reference_dbname'             => $self->o('dna_dbname'),
    'reference_db_server'          => $self->o('dna_db_server'),
    'reference_db_port'            => $self->o('dna_db_port'),
    'uniprot_db_server'                => $self->o('host'), #!!!!!!!!!!!
    'uniprot_db_port'                  => $self->o('port'), #!!!!!!!!!!!
# UniProt output database, this is where gene/transcript models aligning UniProt proteins
# to the genome will be written
'uniprot_db' => {
                    -dbname    => $self->o('uniprot_dbname'),
                    -host      => $self->o('uniprot_db_server'),
                    -port      => $self->o('uniprot_db_port'),
                    -user      => $self->o('user'),
                    -pass      => $self->o('password'),
                    -driver => $self->o('hive_driver'),
},


# Core db from exercise 1
    'reference_db' => {
      -dbname => $self->o('reference_dbname'),
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      $self->db_cmd('CREATE TABLE '.$self->o('protein_table_name').' ('.
                    'accession varchar(50) NOT NULL,'.
                    'biotype varchar(50) NOT NULL,'.
                    'seq text NOT NULL,'.
                    'PRIMARY KEY (accession))'),
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

	{
              -logic_name => 'create_uniprot_output_db',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
              -parameters => {
                               source_db => $self->o('dna_db'),
                               target_db => $self->o('uniprot_db'),
                               create_type => 'clone',
                               enscode_root_dir => $self->o('enscode_root_dir'),
                               user_r => $self->o('user_r'),
                               user_w => $self->o('user'),
                               pass_w => $self->o('password'),
	      },
              -rc_name    => 'default',
              -flow_into => {
                              '1->A' => ['load_protein_file'],
                              'A->1' => ['classify_uniprot_models'],
	      },
              -input_ids  => [ {} ],
	},

	{
              -logic_name => 'load_protein_file',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadProteins',
              -parameters => {
                               sequence_table_name => $self->o('protein_table_name'),
                               protein_file => $self->o('protein_file'),
                               load_biotype => 1,
  	                     },
              -flow_into => {
                              1 => ['generate_protein_jobs'],
	      },
              -rc_name    => 'default',
	},
      {
        -logic_name => 'generate_protein_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'sequence_accession',
                         batch_size => $self->o('protein_batch_size'),
                         sequence_table_name => $self->o('protein_table_name'),
                       },
        -rc_name      => 'default',
        -flow_into => {
                        2 => ['uniprot_alignment'],
                      },
      },

	{
          -logic_name => 'uniprot_alignment',
          -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
          -parameters => {
                           iid_type => 'db_seq',
                           sequence_table_name => $self->o('protein_table_name'),
                           dna_db => $self->o('dna_db'),
                           target_db => $self->o('uniprot_db'),
                           logic_name => 'uniprot_alignment',
                           module     => 'HiveExonerate2Genes',
                           calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
			   %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic',
                                                   'exonerate_protein',
                                                   {FILTER => {
                                                                OBJECT => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                                                                FILTER_PARAMS => {
                                                                                   -coverage => $self->o('exonerate_cov'),
                                                                                   -percent_id => $self->o('exonerate_pid'),
                                                                                   -best_in_genome => 1,
                                                                                   -reject_processed_pseudos => 1,
                                                                                 },
                                                              },
                                                    GENOMICSEQS => $self->o('genome_file'),
                                                    PROGRAM => $self->o('exonerate_path'),
                                                    SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),})}
	                 },
          -rc_name    => 'default',
	},

	{
              -logic_name => 'classify_uniprot_models',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
              -parameters => {
                               classification_type => 'standard',
                               update_gene_biotype => 1,
                               target_db => $self->o('uniprot_db'),
	      },
              -rc_name    => 'default',

	},

    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    };
}

sub resource_classes {
    my $self = shift;
    return {
      'local' => {'LOCAL' => ''},
    }
}

1;
