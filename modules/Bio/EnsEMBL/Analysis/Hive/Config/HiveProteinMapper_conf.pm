=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package HiveProteinMapper_conf;

use warnings;
use strict;
use feature 'say';
use File::Spec::Functions;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;
    return {
      # inherit other stuff from the base class
	    %{ $self->SUPER::default_options() },

'pipeline_name' => 'gifts_uniprot_mapper',

'enscode_root_dir' => '/path/to/enscode/',
'user' => '', # read-only user
'user_r' => '', # read-only user
'user_w' => '', # write user
'password' => '', # write password
'driver' => 'mysql',
'clone_db_script_path' => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/clone_database.ksh', # no need to modify this

'species_name' => 'human',
'repeat_logic_names' => ['repeatmask_repbase_human','dust'],

# database details for the eHive pipe database
'server1' => '',
'port1' => '',
'pipeline_dbname' => '', # this db will be created

# database details for the GIFTS database
'gifts_dbname' => '',
'gifts_dbserver' => '',
'gifts_dbport' => '',
'gifts_dbuser' => '',# read-only user

# database details for the killlist database
'killlist_dbname' => '',
'killlist_dbserver' => '',
'killlist_dbport' => '',

# database details for the other databases
'server2' => '',
'port2' => '',
'core_dbname' => '', # core db containing the genome sequence (dna) for the proteins to be aligned to
'genblast_dbname' => '', # output db which will be created to store the unmapped proteins aligned to the genome

'output_path' => '',
'homology_models_path'       => $self->o('output_path').'/homology_models',

'uniprot_tax_id'             => 9606,
'uniprot_index_name'         => 'uniprot_index',
'uniprot_db_name'            => 'uniprot_db',
'uniprot_query_dir_name'     => 'uniprot_temp',
'uniprot_genblast_batch_size' => 5,
'uniprot_table_name'         => 'uniprot_sequences',

'genblast_path'              => 'genblast',
'genblast_eval'              => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
'genblast_cov'               => '0.5',
'genblast_pid'               => '50',
'genblast_max_rank'          => '5',
'blast_type'                 => 'ncbi',

        'pipeline_db' => {
            # connection parameters
            -dbname => $self->o('pipeline_dbname'),
            -host   => $self->o('server1'),
            -port   => $self->o('port1'),
            -user   => $self->o('user_w'),
            -pass   => $self->o('password'),
            -driver => $self->o('driver'),
        },

        'core_db' =>      {
                            -dbname => $self->o('core_dbname'),
                            -host   => $self->o('server2'),
                            -port   => $self->o('port2'),
                            -user   => $self->o('user'),
                          },

        'genblast_db' =>  {
                            -dbname => $self->o('genblast_dbname'),
                            -host   => $self->o('server2'),
                            -port   => $self->o('port2'),
                            -user   => $self->o('user_w'),
                            -pass   => $self->o('password'),
                          },

        'gifts_db' =>     {
                            -dbname => $self->o('gifts_dbname'),
                            -host   => $self->o('gifts_dbserver'),
                            -port   => $self->o('gifts_dbport'),
                            -user   => $self->o('gifts_dbuser'),
                          },
        'killlist_db' =>  {
                            -dbname => $self->o('killlist_dbname'),
                            -host   => $self->o('killlist_dbserver'),
                            -port   => $self->o('killlist_dbport'),
                            -user   => $self->o('user'),
                          },
    };
  }

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
	    @{$self->SUPER::pipeline_create_commands},
	    $self->hive_data_table('protein', $self->o('uniprot_table_name')),
    ];
  }


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    my %genblast_params = (
      wu    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
      ncbi  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -rl 5000',
      wu_genome    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
      ncbi_genome  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -rl 5000 -softmask',
      wu_projection    => '-P wublast -gff -e #blast_eval# -c #blast_cov# -n 100 -rl 5000 -x 5 ',
      ncbi_projection  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -rl 5000',
      );
    my %commandline_params = (
      'ncbi' => '-num_threads 3 -window_size 40',
      'wu' => '-cpus 3 -hitdist 40',
      'legacy_ncbi' => '-a 3 -A 40',
      );

    return [
      {
        -input_ids => [{}],
        -logic_name => 'dump_softmasked_toplevel',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
        -parameters => {
                         'coord_system_name'    => 'toplevel',
                         'target_db'            => $self->o('core_db'),
                         'output_path'          => $self->o('output_path')."/genome_dumps/",
                         'enscode_root_dir'     => $self->o('enscode_root_dir'),
                         'species_name'         => $self->o('species_name'),
                         'repeat_logic_names'   => $self->o('repeat_logic_names'),
                         'alt_as_scaffolds'     => 1,
                       },
        -flow_into => {
          1 => ['format_softmasked_toplevel'],
        },
        -rc_name    => 'default_himem',
      },

      {
        -logic_name => 'format_softmasked_toplevel',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'    => 'if [ "'.$self->o('blast_type').'" = "ncbi" ]; then convert2blastmask -in '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa -parse_seqids -masking_algorithm repeatmasker -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa.asnb;makeblastdb -in '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa -dbtype nucl -parse_seqids -mask_data '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa.asnb -title "'.$self->o('species_name').'"; else xdformat -n '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa;fi',
                       },
        -rc_name    => 'default_himem_10000',
        -flow_into => { 1 => ['create_genblast_output_db'] },
      },

      {
        -logic_name => 'create_genblast_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('core_db'),
                         target_db => $self->o('genblast_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         '1->A' => ['list_unmapped_uniprot_accessions'],
                         'A->1' => ['classify_genblast'],
                       },
      },

      {
        -logic_name => 'list_unmapped_uniprot_accessions',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
                          db_conn => $self->o('gifts_db'),
                          inputquery => "SELECT uniprot_acc ".
                                        "FROM uniprot_unmapped ".
                                        "WHERE uniprot_tax_id=".$self->o('uniprot_tax_id')." ".
                                        "  AND mapping_history_id=(SELECT max(mapping_history_id) ".
                                                                  "FROM uniprot_unmapped ".
                                                                  "WHERE uniprot_tax_id=".$self->o('uniprot_tax_id').")",
                          step => 100, # 100 is the maximum limit for the Uniprot query which is run in "download_uniprot_files"
                       },
                       -flow_into => { '2->A' => { 'download_uniprot_files' => {'iid' => '#_range_list#'}},
                                       'A->1' => [ 'generate_genblast_jobs' ],
                                     },
                       -rc_name => 'default',
      },

      {
        -logic_name => 'download_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         base_url => "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&format=fasta&accession=",
                         query_url => '#expr(join ",",@{#iid#} )expr#',
                         output_path => $self->o('output_path'),
                         file_name => '#expr(shift @{#iid#})expr#'.".fasta", # cannot use file names consisting of 100 accessions, use the 1st one only
                         dest_dir => $self->o('output_path'),
                       },
        -rc_name => 'default',
        -flow_into => {
                        2 => ['process_uniprot_files'], # this has to be branch 2 because that's the branch the output ids will flow into as defined in the module by the default value of '_branch_to_flow_to'
                      },
      },

      {
        -logic_name => 'process_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
        -parameters => {
                         uniprot_db_name => $self->o('uniprot_db_name'),
                         uniprot_index_name => $self->o('uniprot_index_name'),
                         dest_dir   => $self->o('homology_models_path'),
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                         sequence_table_name => $self->o('uniprot_table_name'),
                      },
        -analysis_capacity => 200,
        -rc_name => 'default',
      },

      {
        -logic_name => 'generate_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'sequence_accession',
                         batch_size => $self->o('uniprot_genblast_batch_size'),
                         sequence_table_name => $self->o('uniprot_table_name'),
                         target_db => $self->o('pipeline_db'),
                       },
        -rc_name      => 'default_himem',
        -flow_into => {
                        2 => ['genblast'],
                      },
      },

      {
        -logic_name => 'genblast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('core_db'),
                         target_db => $self->o('genblast_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa',
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '2h',
                         blast_eval => $self->o('genblast_eval'),
                         blast_cov  => $self->o('genblast_cov'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        -1 => ['split_genblast'],
                        -2 => ['split_genblast'],
                        -3 => ['split_genblast'],
                      },
        -hive_capacity => 1000,
      },

      {
        -logic_name => 'split_genblast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk',
                         batch_size => 1,
                       },
        -rc_name      => 'default',
        -can_be_empty  => 1,
        -flow_into => {
                        2 => ['genblast_retry'],
                      },
      },

      {
        -logic_name => 'genblast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('core_db'),
                         target_db => $self->o('genblast_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa',
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '1h',
                         blast_eval => $self->o('genblast_eval'),
                         blast_cov  => $self->o('genblast_cov'),
                       },
        -rc_name          => 'default_himem_4900',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
        -flow_into => {
                        -1 => ['failed_genblast'],
                        -2 => ['failed_genblast'],
                        -3 => ['failed_genblast'],
                      },
        -hive_capacity => 1000,
      },

      {
        -logic_name => 'failed_genblast',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

      {
        -logic_name => 'classify_genblast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                         target_db => $self->o('genblast_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['dummy'],
                      },
      },

      {
        -logic_name => 'dummy',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {},
        -rc_name          => 'default',
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
      'default' => { LSF => '-M1900 -R"select[mem>1900] rusage[mem=1900]"' },
      'default_himem' => { LSF => '-M2900 -R"select[mem>2900] rusage[mem=2900]"' },
      'default_himem_4900' => { LSF => '-M4900 -R"select[mem>4900] rusage[mem=4900]"' },
      'default_himem_10000' => { LSF => '-M10000 -R"select[mem>10000] rusage[mem=10000]"' },
      'default_20GB' => { LSF => '-M20000 -R"select[mem>20000] rusage[mem=20000]"' },
    }
  }

1;
