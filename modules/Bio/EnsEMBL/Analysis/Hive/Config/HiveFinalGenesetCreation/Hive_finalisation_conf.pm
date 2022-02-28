=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

package Hive_finalisation_conf;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },


##########################################################################
#                                                                        #
# CHANGE STUFF HERE                                                      #
#                                                                        #
##########################################################################
        'output_path'                => '',
        'pipeline_name'              => '',
        'release_number'             => '',
        'species'                    => '',
        'user_r'                     => '',
        'user_w'                     => '',
        'password'                   => '',
        'port'                       => '',

        'pipe_db_name'                => '',
        'reference_db_name'           => '',
        'dna_db_name'                 => '',
        'genblast_db_name'            => '',
        'projection_db_name'          => '', # If you have projection models
        'cdna_db_name'                => '',
        'rnaseq_db_name'              => '',
        'layering_db_name'            => '',
        'utr_source_db_name'          => '',
        'utr_output_db_name'          => '',
        'genebuilder_db_name'         => '',
        'pseudo_db_name'              => '',
        'ncrna_db_name'               => '',
        'final_geneset_db_name'       => '',


        'pipe_db_server'             => '',
        'reference_db_server'        => '',
        'dna_db_server'              => '',
        'genblast_db_server'         => '',
        'projection_db_server'       => '',
        'cdna_db_server'             => '',
        'rnaseq_db_server'           => '',
        'layering_db_server'         => '',
        'utr_source_db_server'       => '',
        'utr_output_db_server'       => '',
        'genebuilder_db_server'      => '',
        'pseudo_db_server'           => '',
        'ncrna_db_server'            => '',
        'final_geneset_db_server'    => '',

        'layer_set_name' => '', # This should correspond to a layering hash in sub layering_set at the bottom of this config

# layering

        'layering_input_gene_dbs' => [
                                       $self->o('genblast_db'),
                                       $self->o('rnaseq_db'),
                                       $self->o('projection_db'),
                                     ],
# UTR
        'utr_gene_dbs' => {
                            'cdna_db'       => $self->o('cdna_db'),
                            'utr_source_db' => $self->o('utr_source_db'),
                            'no_utr_db'     => $self->o('layering_db'),
                          },

        'utr_biotype_priorities'        => { 'best' => 1,
                                             'cdna' => 1,
                                             'cdna_predicted' => 2,
                                           },

# cleaning
'cleaning_blessed_biotypes' => {'pseudogene' => 1, 'processed_pseudogene' => 1}, # These are biotypes not to touch during the cleaning process
'remove_duplicates_script_path' => '/nfs/users/nfs_f/fm2/enscode/ensembl-personal/fm2/scripts/find_and_remove_duplicates.pl',


# Misc stuff

'create_type' => 'clone',
'driver' => 'mysql',
'num_tokens' => 10,
'user' => 'ensro',


# DBs, these are all set with info above
'pipeline_db' => {
                   -dbname => $self->o('pipe_db_name'),
                   -host   => $self->o('pipe_db_server'),
                   -port   => $self->o('port'),
                   -user   => $self->o('user_w'),
                   -pass   => $self->o('password'),
                   -driver => $self->o('driver'),
                 },

'reference_db' => {
                    -dbname => $self->o('reference_db_name'),
                    -host   => $self->o('reference_db_server'),
                    -port   => $self->o('port'),
                    -user   => $self->o('user_r'),
                  },

'dna_db' => {
              -dbname => $self->o('dna_db_name'),
              -host   => $self->o('dna_db_server'),
              -port   => $self->o('port'),
              -user   => $self->o('user_r'),
            },

'genblast_db' => {
  -dbname => $self->o('genblast_db_name'),
  -host   => $self->o('genblast_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},


'projection_db' => {
  -dbname => $self->o('projection_db_name'),
  -host   => $self->o('projection_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'cdna_db' => {
  -dbname => $self->o('cdna_db_name'),
  -host   => $self->o('cdna_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

         'rnaseq_db' => {
           -dbname => $self->o('rnaseq_db_name'),
           -host   => $self->o('rnaseq_db_server'),
           -port   => $self->o('port'),
           -user   => $self->o('user_w'),
           -pass   => $self->o('password'),
        },

        'layering_db' => {
                           -dbname => $self->o('layering_db_name'),
                           -host   => $self->o('layering_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                         },

        'utr_source_db' => {
                           -dbname => $self->o('utr_source_db_name'),
                           -host   => $self->o('utr_source_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_r'),
                        },

        'utr_output_db' => {
                           -dbname => $self->o('utr_output_db_name'),
                           -host   => $self->o('utr_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                        },

        'genebuilder_db' => {
                              -dbname => $self->o('genebuilder_db_name'),
                              -host   => $self->o('genebuilder_db_server'),
                              -port   => $self->o('port'),
                              -user   => $self->o('user_w'),
                              -pass   => $self->o('password'),
                            },


        'pseudo_db' => {
                                -dbname => $self->o('pseudo_db_name'),
                                -host   => $self->o('pseudo_db_server'),
                                -port   => $self->o('port'),
                                -user   => $self->o('user_w'),
                                -pass   => $self->o('password'),
                              },

        'ncrna_db' => {
                                -dbname => $self->o('ncrna_db_name'),
                                -host   => $self->o('ncrna_db_server'),
                                -port   => $self->o('port'),
                                -user   => $self->o('user_w'),
                                -pass   => $self->o('password'),
                              },

        'final_geneset_db' => {
                                -dbname => $self->o('final_geneset_db_name'),
                                -host   => $self->o('final_geneset_db_server'),
                                -port   => $self->o('port'),
                                -user   => $self->o('user_w'),
                                -pass   => $self->o('password'),
                              },

     };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

      {
        -logic_name => 'create_layering_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('layering_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -input_ids => [{}],
        -flow_into => {
                        1 => ['create_utr_db'],
                      },
      },

      {
        -logic_name => 'create_utr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('utr_output_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_genebuilder_db'],
                      },
     },

      {
        -logic_name => 'create_genebuilder_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('genebuilder_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['create_toplevel_slices'],
                        'A->1' => ['create_pseudogene_db'],
                      },
      },

      {
        -logic_name => 'create_toplevel_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'slice',
                         coord_system_name => 'toplevel',
                         include_non_reference => 0,
                         top_level => 1,
                         # These options will create only slices that have a gene on the slice in one of the feature dbs
                         feature_constraint => 1,
                         feature_type => 'gene',
                         feature_dbs => [$self->o('genblast_db'),$self->o('projection_db'),$self->o('rnaseq_db')],
                      },
        -flow_into => {
                       '2' => ['layer_annotation'],
                      },

        -rc_name    => 'default',
      },


      {
        -logic_name => 'layer_annotation',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLayerAnnotation',
        -parameters => {
                         dna_db     => $self->o('dna_db'),
                         logic_name => 'layer_annotation',
                         module     => 'HiveLayerAnnotation',
                         config_settings => $self->get_config_settings('layer_annotation','layers'),
                       },
        -rc_name    => 'layer_annotation',
        -flow_into  => {
                         '1->A' => ['split_slices_on_intergenic'],
                         'A->1' => ['genebuilder'],
                       },
      },

      {
        -logic_name => 'split_slices_on_intergenic',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         input_gene_dbs => $self->o('utr_gene_dbs'),
                         iid_type => 'slice',
                       },
        -batch_size => 100,
        -hive_capacity => 200,
        -rc_name    => 'transcript_finalisation',
        -flow_into => {
                        2 => ['cluster_input_genes'],
                      },
      },


     {
        -logic_name => 'cluster_input_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClusterSourceGenes',
        -parameters => {
                         logic_name => 'cluster_input_genes',
                         dna_db => $self->o('dna_db'),
                         input_gene_dbs => $self->o('utr_gene_dbs'),
                         allowed_input_sets => undef,
                         iid_type => 'slice',
                       },
        -batch_size => 20,
        -hive_capacity => 200,
        -rc_name    => 'transcript_finalisation',
        -flow_into => {
                        2 => ['run_utr_addition'],
                      },

      },

      {
        -logic_name => 'run_utr_addition',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
        -parameters => {
                         logic_name => 'utr_addition',
                         dna_db => $self->o('dna_db'),
                         input_gene_dbs => $self->default_options()->{'utr_gene_dbs'},
                         utr_biotype_priorities => $self->o('utr_biotype_priorities'),
                         utr_output_db => $self->default_options()->{'utr_output_db'},
                         iid_type => 'slice',
                       },
        -batch_size => 20,
        -hive_capacity => 900,
        -rc_name    => 'transcript_finalisation',
     },

      {
        -logic_name => 'genebuilder',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGeneBuilder',
        -parameters => {
                         layering_output_db => $self->o('utr_output_db'),
                         genebuilder_output_db => $self->o('genebuilder_db'),
                         dna_db     => $self->o('dna_db'),
                         logic_name => 'ensembl',
                         module     => 'HiveGeneBuilder',
                         config_settings => $self->get_config_settings('genebuilder','genebuilder_set'),
                       },
        -rc_name    => 'transcript_finalisation',
        -hive_capacity => 900,
      },

      {
        -logic_name => 'create_pseudogene_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('pseudo_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['create_pseudogene_ids'],
                        'A->1' => ['create_final_geneset_db'],
                      },
      },

            {
        -logic_name => 'create_pseudogene_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('genebuilder_db'),
                         iid_type => 'feature_id',
                         feature_type => 'gene',
                      },
        -flow_into => {
                       2 => ['pseudogenes'],
                      },
        -rc_name    => 'default',
      },

      {
        -logic_name => 'pseudogenes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePseudogenes',
        -parameters => {
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->o('reference_db'),
                         output_db => $self->o('pseudo_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'pseudogenes',
                         module     => 'HivePseudogenes',
                         config_settings => $self->get_config_settings('pseudogenes','pseudogenes_set'),
                       },
        -batch_size => 20,
        -hive_capacity => 900,
        -rc_name    => 'transcript_finalisation',
        -flow_into => {
                       2 => ['spliced_elsewhere'],
                      },
      },

      {
        -logic_name => 'concat_blast_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'for i in '.$self->o('output_path').'/pseudogenes/multi_exon_dir/multi_exon_seq*.fasta;'.
                                'do cat $i >> '.$self->o('output_path').'/pseudogenes/multi_exon_dir/all_multi_exon_genes.fasta;'.
                                'done'
                       },
         -rc_name => 'default',
         -wait_for => ['pseudogenes'],
         -input_ids => [{}],
         -flow_into => { 1 => ['format_blast_db'] },
      },

      {
        -logic_name => 'format_blast_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'xdformat -n '.$self->o('output_path').'/pseudogenes/multi_exon_dir/all_multi_exon_genes.fasta'
                       },
         -rc_name => 'default',
      },

      {
        -logic_name => 'spliced_elsewhere',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplicedElsewhere',
        -parameters => {
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->o('reference_db'),
                         output_db => $self->o('pseudo_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'spliced_elsewhere',
                         module     => 'HiveSplicedElsewhere',
                         config_settings => $self->get_config_settings('pseudogenes','pseudogenes_set'),
                       },
        -rc_name          => 'transcript_finalisation',
        -can_be_empty  => 1,
        -wait_for => ['format_blast_db'],
      },

      {
        -logic_name => 'create_final_geneset_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('pseudo_db'),
                         target_db => $self->o('final_geneset_db'),
                         create_type => 'copy',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['run_cleaner'],
                      },
      },


      {
        -logic_name => 'run_cleaner',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanGeneset',
        -parameters => {
                         input_db => $self->o('final_geneset_db'),
                         dna_db => $self->o('dna_db'),
                         output_path => $self->o('output_path').'/clean_genes/',
                         blessed_biotypes => $self->o('cleaning_blessed_biotypes'),
                         flagged_redundancy_coverage_threshold => 95,
                         general_redundancy_coverage_threshold => 95,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['delete_flagged_transcripts'],
                      },
      },


     {
       -logic_name => 'delete_flagged_transcripts',
       -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts',
       -parameters => {
                        dbhost => $self->o('final_geneset_db','-host'),
                        dbname => $self->o('final_geneset_db','-dbname'),
                        dbuser => $self->o('user_w'),
                        dbpass => $self->o('password'),
                        dbport => $self->o('final_geneset_db','-port'),
                        transcript_ids_file => $self->o('output_path').'/clean_genes/transcript_ids_to_remove.txt',
                        delete_transcripts_path => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/genebuild/',
                        delete_genes_path => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/genebuild/',
                        delete_transcripts_script_name => 'delete_transcripts.pl',
                        delete_genes_script_name => 'delete_genes.pl',
                        output_path => $self->o('output_path').'/clean_genes/',
                        output_file_name => 'delete_transcripts.out',
                      },
        -max_retry_count => 0,
        -flow_into => {
                        '1' => ['transfer_ncrnas'],
                      },
      },



     {
       -logic_name => 'transfer_ncrnas',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
                        cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl-analysis/scripts/genebuild/copy_genes.pl'.
                                     ' -sourcehost '.$self->o('ncrna_db','-host').
                                     ' -sourceuser '.$self->o('user_r').
                                     ' -sourceport '.$self->o('ncrna_db','-port').
                                     ' -sourcedbname '.$self->o('ncrna_db','-dbname').
                                     ' -dnauser '.$self->o('user_r').
                                     ' -dnahost '.$self->o('dna_db','-host').
                                     ' -dnaport '.$self->o('dna_db','-port').
                                     ' -dnadbname '.$self->o('dna_db','-dbname').
                                     ' -targetuser '.$self->o('user_w').
                                     ' -targetpass '.$self->o('password').
                                     ' -targethost '.$self->o('final_geneset_db','-host').
                                     ' -targetport '.$self->o('final_geneset_db','-port').
                                     ' -targetdbname '.$self->o('final_geneset_db','-dbname').
                                     ' -all'
                      },
        -rc_name => 'default',
        -flow_into => {
                        '1' => ['delete_duplicate_genes'],
                      },
      },

     {
       -logic_name => 'delete_duplicate_genes',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
                        cmd => "perl ".$self->o('remove_duplicates_script_path')
                               ." -dbhost ".$self->o('final_geneset_db','-host')
                               ." -dbuser ".$self->o('user_w')
                               ." -dbpass ".$self->o('password')
                               ." -dbname ".$self->o('final_geneset_db','-dbname')
                               ." -dbport ".$self->o('final_geneset_db','-port')
                               ." -dnadbhost ".$self->o('dna_db','-host')
                               ." -dnadbuser ".$self->o('user_r')
                               ." -dnadbname ".$self->o('dna_db','-dbname')
                               ." -dnadbport ".$self->o('dna_db','-port'),
                     },
        -max_retry_count => 0,
        -rc_name => 'default',
      },


    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    };
}

# override the default method, to force an automatic loading of the registry in all workers
#sub beekeeper_extra_cmdline_options {
#    my $self = shift;
#    return "-reg_conf ".$self->o("registry");
#}

sub resource_classes {
    my $self = shift;
    return {
      'default' => { LSF => '-q normal -M1900 -R"select[mem>1900] rusage[mem=1900,myens_build11tok=10,myens_build12tok=10,myens_build13tok=10]"' },
      'transcript_finalisation' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,myens_build11tok=10,myens_build12tok=10,myens_build4tok=10,myens_build13tok=10]"' },
      'process_clusters' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,myens_build11tok=10,myens_build12tok=10,myens_build13tok=10]"' },
      'blast' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,myens_build11tok=10,myens_build12tok=10,myens_build13tok=10]"' },
      'layer_annotation' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,myens_build11tok=10,myens_build12tok=10,myens_build13tok=10]"' },
      'filter' => { LSF => '-q normal -M9900 -R"select[mem>9900] rusage[mem=9900,myens_build11tok=10,myens_build12tok=10,myens_build13tok=10]"' },
    }
}

sub get_config_settings {

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

  layer_annotation => {
    Default => {
      TARGETDB_REF => $self->o('layering_db'),
      SOURCEDB_REFS => $self->o('layering_input_gene_dbs'),
      # ordered list of annotation layers. Genes from lower layers
      # are only retained if they do not "interfere" with genes from
      # higher layers. Genes in "Discard" layers are when assessing
      # interference, but are not written to the final database

      # Filtering is using done at the exon-overlap level
      # When no FILTER exists in this file, this is the default behaviour

      # If you would like to filter in a different way, please specify filter
      #FILTER => 'Bio::EnsEMBL::Analysis::Tools::GenomeOverlapFilter',
      #FILTER => 'Bio::EnsEMBL::Analysis::Tools::AllExonOverlapFilter',
      FILTER => 'Bio::EnsEMBL::Analysis::Tools::CodingExonOverlapFilter',
    },

    layers => {
      LAYERS => $self->layering_set(),
    },
  },

  genebuilder => {
    Default => {
               },

    genebuilder_set => {
      INPUT_GENES => { 'input_db' => $self->genebuilder_set()},
      OUTPUT_BIOTYPE => 'ensembl',
      MAX_TRANSCRIPTS_PER_CLUSTER => 10,
      MIN_SHORT_INTRON_LEN => 7, #introns shorter than this seem to be real frame shifts and shoudn't be ignored
      MAX_SHORT_INTRON_LEN => 15,
      BLESSED_BIOTYPES => {
                            'ccds_gene' => 1,
                            'Blessed_UTR_Genes' => 1,
                           },
      #the biotypes of the best genes always to be kept
      MAX_EXON_LENGTH => 20000,
      #if the coding_only flag is set to 1, the transcript clustering into genes is done over coding exons only
      # the current standard way is to cluster only on coding exons
      CODING_ONLY => 1,
    },
  },

                                 pseudogenes => {
                                   Default => {
               # you can set the input- and output database - the names should point to
               # keys in Database.pm
               PS_INPUT_DATABASE  => 'GENEBUILD_DB',
               PS_OUTPUT_DATABASE => 'PSEUDO_DB',

               # configs for the introns in repeats test

               # introns longer than the following are considered "real"
               PS_FRAMESHIFT_INTRON_LENGTH => 9,
               # total length of introns
               PS_MAX_INTRON_LENGTH   => '5000',
               # Types of repeats to run the anaysis with
               PS_REPEAT_TYPES =>  ['LINE','LTR','SINE'],
               # max percent coverage of the introns with the above repeats
               PS_MAX_INTRON_COVERAGE => '80',
               # max allowed exon coverage with the above repeats
               PS_MAX_EXON_COVERAGE   => '99',
               PS_NUM_FRAMESHIFT_INTRONS  => 1,
               PS_NUM_REAL_INTRONS  => 1,
               # biotype of genes to check
               PS_BIOTYPE  => 'protein_coding',

               # Blessed genes dont get called pseudogenes
               # Biotype is a transcript biotype
               BLESSED_BIOTYPES => { 'ccds_gene' => 1 },

               # configs for the spliced elsewhere tests
               # %ID of a tbalstx of the (presumed) retrotransposed query sequence to its
               # homolog that is spliced elsewhere in the genome. hits falling below
               # this cutoff are ignored (80%) is suggested
               PS_PERCENT_ID_CUTOFF   => 40,
               PS_P_VALUE_CUTOFF   => '1.0e-50',
               PS_RETOTRANSPOSED_COVERAGE   => 80,
               PS_ALIGNED_GENOMIC  => 100,
               # logic name to give to pseudogenes
               PS_PSEUDO_TYPE      => 'pseudogene',

               # if a gene is found to be a pseudogene, its gene biotype will be changed to
               # PS_PSEUDO_TYPE. By default, the biotype of its transcript will also be changed
               # to PS_PSEUDO_TYPE.  If you want to keep the original transcript biotype
               # instead (so you can keep track of what type of models actually got turned into a
               # pseudogene), set KEEP_TRANS_BIOTYPE to 1.

               KEEP_TRANS_BIOTYPE  => 0,

               # logic name to give genes with exons covered by repeats
               # if left blank they will just get deleted (recommended)
               PS_REPEAT_TYPE      => '',

               # analysis logic names to run over genes falling into these categories
               SINGLE_EXON      => 'spliced_elsewhere',
               INDETERMINATE    => '',
               RETROTRANSPOSED  => '',
               # if you dont wish to run further tests on retro transposed genes
               # What type would you like to give them?
               # previously set to 'retrotransposed', we change to 'processed_pseudogene' from e70 onwards.
               RETRO_TYPE       => 'processed_pseudogene',

               SPLICED_ELSEWHERE_LOGIC_NAME => 'spliced_elsewhere',
               PSILC_LOGIC_NAME => 'Psilc',
               # SPLICED ELSEWHERE SPECIFIC CONFIG
               # ratio of the spans of the retrotransposed gene vs its spliced homologue
               # spliced / retrotransposed
               # ie: 1 is the same length genes
               # many retrotransposed genes have a ratio > 10
               # used to make retrotransposition decision
               PS_SPAN_RATIO          => 3,
               # mimimum number of exons for the spliced gene to have
               PS_MIN_EXONS           => 4,
               # path of blast db of multi exon genes
               PS_MULTI_EXON_DIR       => "/path/to/my/blast/directory/" ,
               # Chunk size
               PS_CHUNK => '50',
               DEBUG => '1',
                                              },

                                     pseudogenes_set => {
                PS_INPUT_DATABASE  => 'GBUILD_DB',
                 # biotype of genes to check
                 PS_BIOTYPE  => 'protein_coding',
                 # path of blast db of multi exon genes
                 PS_MULTI_EXON_DIR       => $self->o('output_path').'/pseudogenes/multi_exon_dir/',
                                                        },

                spliced_elsewhere_set => {
                                   PS_INPUT_DATABASE  => 'GBUILD_DB',
                 # biotype of genes to check
                 PS_BIOTYPE  => 'ensembl_utr',
                 # path of blast db of multi exon genes
                 PS_MULTI_EXON_DIR       =>$self->o('output_path').'/pseudogenes/multi_exon_dir/',
                },

              }, # end pseudogenes

                                                                  rfam => {
                                   Default => {

                                              },

                                   rfam_set => {
                                                 BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
                                                 PARSER_PARAMS => {
                                                                    -regex => '(\w+)\.\w+',
                                                                    -query_type => 'dna',
                                                                    -database_type => 'dna',
                                                                    -threshold => 0.01,
                                                                  },
                                                BLAST_FILTER => undef,
                                                FILTER_PARAMS => {},
                                                BLAST_PARAMS => {
                                                                  -unknown_error_string => 'FAILED',
                                                                  -type => 'wu',
                                                                }
                                               },
                                         }, # end rfam

                                 mirna => {
                                   Default => {

                                              },

                                   mirna_set => {
                                                  BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::BPliteWrapper',
                                                  PARSER_PARAMS => {
                                                                     -regex => '\w+\s+(\w+)',
                                                                     -query_type => 'dna',
                                                                     -database_type => 'dna',
                                                                   },
                                                  BLAST_FILTER => undef,
                                                  FILTER_PARAMS => {},
                                                  BLAST_PARAMS => {
                                                                    -unknown_error_string => 'FAILED',
                                                                    -type => 'wu',
                                                                  }
                                                },
                                          }, # end mirna
                               };

  return($master_config_settings->{$config_group});
}


sub layering_set {
  my ($self) = @_;

  my $clade = $self->default_options()->{'layer_set_name'};

  if($clade eq 'primates_basic') {
    return ([
             {
              ID         => 'LAYER1',
              BIOTYPES   => [
                             'realign_95',
                             'realign_80',
                             'rnaseq_95',
                             'rnaseq_80',
                             'self_pe12_sp_95',
                             'self_pe12_tr_95',
                             'self_pe12_sp_80',
                             'self_pe12_tr_80',
                             'human_pe12_sp_95',
                             'human_pe12_tr_95',
                             'primates_pe12_sp_95',
                             'primates_pe12_tr_95',
                             'mammals_pe12_sp_95',
                             'mammals_pe12_tr_95',
                            ],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER2',
              BIOTYPES   => [
                             'human_pe12_sp_80',
                             'human_pe12_tr_80',
                             'primates_pe12_sp_80',
                             'primates_pe12_tr_80',
                             'mammals_pe12_sp_80',
                             'mammals_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'primates_pe3_sp_95',
                             'vert_pe12_sp_95',
                             'vert_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'primates_pe3_sp_80',
                             'vert_pe12_sp_80',
                             'vert_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER5',
              BIOTYPES   => [
                              'realign_50',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },


    ]);
  } elsif($clade eq 'rodents_basic') {
    return ([
             {
              ID         => 'LAYER1',
              BIOTYPES   => ['realign_95','realign_80',
                            'rnaseq_95','rnaseq_80',
                            ],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER2',
              BIOTYPES   => [ 'self_pe12_sp_95','self_pe12_sp_80',
                              'mouse_pe12_sp_95','mouse_pe12_sp_80',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,

            },

             {
              ID         => 'LAYER3',
              BIOTYPES   => [ 'self_pe12_tr_95','self_pe12_tr_80',
                              'mouse_pe12_tr_95','mouse_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER4',
              BIOTYPES   => ['rodents_pe12_sp_95','rodents_pe12_sp_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER5',
              BIOTYPES   => ['rodents_pe12_tr_95','rodents_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
              BIOTYPES   => ['human_pe12_sp_95','human_pe12_sp_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => ['human_pe12_tr_95','human_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER8',
              BIOTYPES   => [
                              'mammals_pe12_sp_95','mammals_pe12_sp_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER9',
              BIOTYPES   => [
                              'rodents_pe3_sp_95','rodents_pe3_sp_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER10',
              BIOTYPES   => [
                              'mammals_pe12_tr_95','mammals_pe12_tr_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9'],
              DISCARD    => 0,
            },


             {
              ID         => 'LAYER11',
              BIOTYPES   => [
                              'rodents_pe3_tr_95','rodents_pe3_tr_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9','LAYER10'],
              DISCARD    => 0,
            },


             {
              ID         => 'LAYER12',
              BIOTYPES   => [
                              'vert_pe12_sp_95','vert_pe12_sp_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9','LAYER10','LAYER11'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER13',
              BIOTYPES   => [
                              'vert_pe12_tr_95','vert_pe12_tr_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9','LAYER10','LAYER11','LAYER12'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER14',
              BIOTYPES   => [
                              'realign_50'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9','LAYER10','LAYER11','LAYER12','LAYER13'],
              DISCARD    => 0,
            },
            ]);
  } else {
    die "Unknown clade selected for layering: ".$clade;
  }

}

sub genebuilder_set {
  my ($self) = @_;

  my $clade = $self->default_options()->{'uniprot_set'};

  if($clade eq 'primates_basic') {
    return ([
              'realign_80',
              'realign_95',
              'realign_50',
              'self_pe12_sp_95',
              'self_pe12_sp_80',
              'self_pe12_tr_95',
              'self_pe12_tr_80',
              'rnaseq_95',
              'rnaseq_80',
              'vert_pe12_sp_95',
              'vert_pe12_sp_80',
              'vert_pe12_tr_80',
              'vert_pe12_tr_95',
              'primates_pe12_sp_80',
              'primates_pe12_tr_80',
              'primates_pe12_tr_95',
              'primates_pe12_sp_95',
              'human_pe12_sp_80',
              'human_pe12_sp_95',
              'human_pe12_tr_80',
              'human_pe12_tr_95',
              'primates_pe3_sp_95',
              'primates_pe3_tr_80',
              'primates_pe3_tr_95',
              'primates_pe3_sp_80',
              'primates_pe45_sp_95',
              'primates_pe45_tr_80',
              'primates_pe45_tr_95',
              'primates_pe45_sp_80',
              'mammals_pe12_sp_80',
              'mammals_pe12_sp_95',
              'mammals_pe12_tr_95',
              'mammals_pe12_tr_80',
    ]);
  } elsif($clade eq 'rodents_basic') {
    return ([
              'realign_80',
              'realign_95',
              'realign_50',
              'self_pe12_sp_95',
              'self_pe12_sp_80',
              'self_pe12_tr_95',
              'self_pe12_tr_80',
              'rnaseq_95',
              'rnaseq_80',
              'mouse_pe12_sp_95',
              'mouse_pe12_sp_80',
              'mouse_pe12_tr_95',
              'mouse_pe12_tr_80',
              'vert_pe12_sp_95',
              'vert_pe12_sp_80',
              'vert_pe12_tr_80',
              'vert_pe12_tr_95',
              'rodents_pe12_sp_80',
              'rodents_pe12_tr_80',
              'rodents_pe12_tr_95',
              'rodents_pe12_sp_95',
              'human_pe12_sp_80',
              'human_pe12_sp_95',
              'human_pe12_tr_80',
              'human_pe12_tr_95',
              'rodents_pe3_sp_95',
              'rodents_pe3_tr_80',
              'rodents_pe3_tr_95',
              'rodents_pe3_sp_80',
              'rodents_pe45_sp_95',
              'rodents_pe45_tr_80',
              'rodents_pe45_tr_95',
              'rodents_pe45_sp_80',
              'mammals_pe12_sp_80',
              'mammals_pe12_sp_95',
              'mammals_pe12_tr_95',
              'mammals_pe12_tr_80',
            ]);
  } else {
    die "Unknown clade selected for genebuilding: ".$clade;
  }

}

1;
