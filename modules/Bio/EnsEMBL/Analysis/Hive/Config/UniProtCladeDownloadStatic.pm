=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic

=head1 SYNOPSIS

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub pipeline_analyses {
    my ($self) = @_;

    return [
      {
        -logic_name => 'run_uniprot_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanPep',
        -parameters => {
                         blast_db_path => $self->o('uniprot_blast_db_path'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscanPep',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep')},
                      },
        -flow_into => {
                        -1 => ['run_uniprot_blast_himem'],
                        -2 => ['run_uniprot_blast_long'],
                      },
        -rc_name    => 'blast',
      },
  ];
}

=head1 DESCRIPTION

This is the config file for all analysis downloading clade data from UniProt. You should use it in your Hive configuration file to
specify the parameters of an analysis. You can either choose an existing config or you can create
a new one based on the default hash. 

=head1 METHODS

  _master_config_settings: contains all possible parameters

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key) = @_;

  my $taxon_ids = {
                   'human_taxon_id'    => '9606',
                   'mammals_taxon_id'  => '40674',
                   'mouse_taxon_id'    => '10090',
                   'primates_taxon_id' => '9443',
                   'rodents_taxon_id'  => '9989',
                   'vert_taxon_id'     => '7742',
                   'fish_taxon_id'     => '7898',
                   'aves_taxon_id'     => '8782',
                 };
  my %config = (
    default => {},

    self_patch => {
               self_pe12 => {
                              file_name => 'self_pe12.fasta',
                              taxon_id  => '#taxon_id#',
                              dest_dir  => '#output_path#',
                              compress  => 0,
                              pe_level  => [1,2],
                            },

               self_frag_pe12 => {
                                   file_name => 'self_frag_pe12.fasta',
                                   taxon_id  => '#taxon_id#',
                                   dest_dir  => '#output_path#',
                                   compress  => 0,
                                   fragment  => 1,
                                   pe_level  => [1,2],
                                 },
    },

    primates_basic => {
              self_pe12 =>{
                            file_name => 'self_pe12.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [1,2],
                          },

              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => '#output_path#',
                              compress  => 0,
                              pe_level  => [1,2],
                            },

              primates_pe12 => {
                             file_name => 'primates_pe12.fasta',
                             taxon_id  => $taxon_ids->{'primates_taxon_id'},
                             exclude_id => ['#taxon_id#',$taxon_ids->{'human_taxon_id'}],
                             dest_dir  => '#output_path#',
                             compress  => 0,
                             pe_level  => [1,2],
                           },


               mammals_pe12 => {
                                 file_name  => 'mammals_pe12.fasta',
                                 taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                                 exclude_id => [$taxon_ids->{'primates_taxon_id'}],
                                 dest_dir   => '#output_path#',
                                 compress   => 0,
                                 pe_level   => [1,2],
                               },
             },

    rodents_basic => {
              self_pe12 =>{
                            file_name => 'self_pe12.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [1,2],
                          },

              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => '#output_path#',
                              compress  => 0,
                              pe_level  => [1,2],
                            },

              mouse_pe12 => {
                              file_name  => 'mouse_pe12.fasta',
                              taxon_id   => $taxon_ids->{'mouse_taxon_id'},
                              dest_dir   => '#output_path#',
                              compress   => 0,
                              pe_level   => [1,2],
                            },

              rodents_pe12 => {
                                file_name => 'rodents_pe12.fasta',
                                taxon_id  => $taxon_ids->{'rodents_taxon_id'},
                                exclude_id => [$taxon_ids->{'mouse_taxon_id'},'#taxon_id#'],
                                dest_dir  => '#output_path#',
                                compress  => 0,
                                pe_level  => [1,2],
                              },

              rodents_pe3 => {
                               file_name => 'rodents_pe3.fasta',
                               taxon_id  => $taxon_ids->{'rodents_taxon_id'},
                               dest_dir  => '#output_path#',
                               compress  => 0,
                               pe_level  => [3],
                             },

               mammals_pe12 => {
                                 file_name  => 'mammals_pe12.fasta',
                                 taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                                 exclude_id => [$taxon_ids->{'rodents_taxon_id'},$taxon_ids->{'human_taxon_id'}],
                                 dest_dir   => '#output_path#',
                                 compress   => 0,
                                 pe_level   => [1,2],
                               },

               vert_pe12 => {
                              file_name  => 'vert_pe12.fasta',
                              taxon_id   => $taxon_ids->{'vert_taxon_id'},
                              exclude_id => [$taxon_ids->{'mammals_taxon_id'}],
                              dest_dir   => '#output_path#',
                              compress   => 0,
                              pe_level   => [1,2],
                            },

             },

    mammals_basic => {
              self_pe12 =>{
                            file_name => 'self_pe12.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [1,2],
                          },

              self_pe3 =>{
                            file_name => 'self_pe3.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [3],
                          },

              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => '#output_path#',
                              compress  => 0,
                              pe_level  => [1,2],
                            },

              mouse_pe12 => {
                              file_name  => 'mouse_pe12.fasta',
                              taxon_id   => $taxon_ids->{'mouse_taxon_id'},
                              dest_dir   => '#output_path#',
                              compress   => 0,
                              pe_level   => [1,2],
                            },


               mammals_pe12 => {
                                 file_name  => 'mammals_pe12.fasta',
                                 taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                                 exclude_id => ['#taxon_id#',$taxon_ids->{'mouse_taxon_id'},$taxon_ids->{'human_taxon_id'}],
                                 dest_dir   => '#output_path#',
                                 compress   => 0,
                                 pe_level   => [1,2],
                               },

               vert_pe12 => {
                              file_name  => 'vert_pe12.fasta',
                              taxon_id   => $taxon_ids->{'vert_taxon_id'},
                              exclude_id => [$taxon_ids->{'mammals_taxon_id'}],
                              dest_dir   => '#output_path#',
                              compress   => 0,
                              pe_level   => [1,2],
                            },

             },

    fish_complete => {
              self_pe12 =>{
                            file_name => 'self_pe12.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [1,2],
                          },

              self_pe3 =>{
                            file_name => 'self_pe3.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [3],
                          },

              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => '#output_path#',
                              compress  => 0,
                              pe_level  => [1,2],
                            },

              fish_pe12 => {
                             file_name => 'fish_pe12.fasta',
                             taxon_id  => $taxon_ids->{'fish_taxon_id'},
                             exclude_id => ['#taxon_id#'],
                             dest_dir  => '#output_path#',
                             compress  => 0,
                             pe_level  => [1,2],
                           },

                fish_pe3 => {
                             file_name => 'fish_pe3.fasta',
                             taxon_id  => $taxon_ids->{'fish_taxon_id'},
                             exclude_id => ['#taxon_id#'],
                             dest_dir  => '#output_path#',
                             compress  => 0,
                             pe_level  => [3],
                           },

               mammals_pe12 => {
                                 file_name  => 'mammals_pe12.fasta',
                                 taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                                 exclude_id => [$taxon_ids->{'human_taxon_id'}],
                                 dest_dir   => '#output_path#',
                                 compress   => 0,
                                 pe_level   => [1,2],
                               },

               vert_pe12 => {
                              file_name  => 'vert_pe12.fasta',
                              taxon_id   => $taxon_ids->{'vert_taxon_id'},
                              exclude_id => [$taxon_ids->{'mammals_taxon_id'}, $taxon_ids->{fish_taxon_id}],
                              dest_dir   => '#output_path#',
                              compress   => 0,
                              pe_level   => [1,2],
                            },

             },
    bird_basic => {
              self_pe12 =>{
                            file_name => 'self_pe12.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [1,2],
                          },
              self_pe3 =>{
                            file_name => 'self_pe3.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [3],
                          },
              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => '#output_path#',
                              compress  => 0,
                              pe_level  => [1,2],
                            },
               vert_pe12 => {
                              file_name  => 'vert_pe12.fasta',
                              taxon_id   => $taxon_ids->{'vert_taxon_id'},
                              exclude_id => [$taxon_ids->{'mammals_taxon_id'}, $taxon_ids->{aves_taxon_id}],
                              dest_dir   => '#output_path#',
                              compress   => 0,
                              pe_level   => [1,2],
                            },
               bird_pe12 => {
                              file_name  => 'aves_pe12.fasta',
                              taxon_id   => $taxon_ids->{'aves_taxon_id'},
                              dest_dir   => '#output_path#',
                              exclude_id => ['#taxon_id#'],
                              compress   => 0,
                              pe_level   => [1,2],
                            },
               mammals_pe12 => {
                                 file_name  => 'mammals_pe12.fasta',
                                 taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                                 exclude_id => [$taxon_ids->{'human_taxon_id'}],
                                 dest_dir   => '#output_path#',
                                 compress   => 0,
                                 pe_level   => [1,2],
                               },
             },

    fish_basic => {
              self_pe12 =>{
                            file_name => 'self_pe12.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [1,2],
                          },

              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => '#output_path#',
                              compress  => 0,
                              pe_level  => [1,2],
                            },

              fish_pe12 => {
                             file_name => 'fish_pe12.fasta',
                             taxon_id  => $taxon_ids->{'fish_taxon_id'},
                             exclude_id => ['#taxon_id#'],
                             dest_dir  => '#output_path#',
                             compress  => 0,
                             pe_level  => [1,2],
                           },

               mammals_pe12 => {
                                 file_name  => 'mammals_pe12.fasta',
                                 taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                                 exclude_id => [$taxon_ids->{'human_taxon_id'}],
                                 dest_dir   => '#output_path#',
                                 compress   => 0,
                                 pe_level   => [1,2],
                               },

               vert_pe12 => {
                              file_name  => 'vert_pe12.fasta',
                              taxon_id   => $taxon_ids->{'vert_taxon_id'},
                              exclude_id => [$taxon_ids->{'mammals_taxon_id'}, $taxon_ids->{fish_taxon_id}],
                              dest_dir   => '#output_path#',
                              compress   => 0,
                              pe_level   => [1,2],
                            },

             },
             selenocysteine => {
               query_url => 'taxonomy%3A#taxon_id#+AND+annotation%3A%28type%3Anon_std+Selenocysteine%29+AND+fragment%3Ano&format=fasta&include=yes',
               file_name => '#taxon_id#_seleno.fa',
               dest_dir   => '#output_path#',
             },
      self_isoforms_12 => {
              self_isoforms_12 =>{
                            file_name => 'self_isoforms_12.fasta',
                            taxon_id  => '#taxon_id#',
                            dest_dir  => '#output_path#',
                            compress  => 0,
                            pe_level  => [1,2],
                            isoforms  => 1,
                            format    => 'fasta',
                          },
              },


     havana_human_blast => {

       human_pe12 => {
                       file_name => 'human_pe12.fasta',
                       taxon_id  => $taxon_ids->{'human_taxon_id'},
                       dest_dir  => '#output_path#',
                       compress  => 0,
                       pe_level  => [1,2],
                     },

       primates_pe12 => {
                          file_name => 'primates_pe12.fasta',
                          taxon_id  => $taxon_ids->{'primates_taxon_id'},
                          exclude_id => [$taxon_ids->{'human_taxon_id'}],
                          dest_dir  => '#output_path#',
                          compress  => 0,
                          pe_level  => [1,2],
                        },

       mammals_pe12 => {
                         file_name  => 'mammals_pe12.fasta',
                         taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                         exclude_id => [$taxon_ids->{'primates_taxon_id'}],
                         dest_dir   => '#output_path#',
                        compress   => 0,
                         pe_level   => [1,2],
                       },

     },

  );
  return $config{$key};
}

1;
