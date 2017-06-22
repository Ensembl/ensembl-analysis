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


=head1 DESCRIPTION

This is the config file for all pipeline sanity checks. It can be used for checking
feature counts versus expected values for a clade

=head1 METHODS

  _master_config_settings: contains all possible parameters

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key, $species) = @_;

  my %config = (
  'default' => {},
  'genome_preparation_checks' => {
     'primates_basic' => {
       # repeats
       'dust'                        => [2000000,'repeat'],
       'repeatmask_repbase_primates' => [3000000,'repeat'],
       'trf'                         => [500000,'repeat'],

       # simple features
       'cpg'                         => [15000,'simple'],
       'trnascan'                    => [300,'simple'],
       'eponine'                     => [30000,'simple'],

       # prediction transcripts
       'genscan'                     => [40000,'prediction transcript'],

       # dna align features
       'unigene'                     => [3000000,'dna align'],
       'vertrna'                     => [3000000,'dna align'],

       # protein align features
       'uniprot'                     => [3000000,'protein align'],
     },
  },

  'gene_db_checks' => {
    'primates_basic' => {
      'genblast' => {
        'logic_names' => {
          'genblast'            => 300000,
          'genblast_not_best'   => 300000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 20000,
          'primates_pe12_'      => 50000,
          'primates_pe3_'       => 80000,
          'mammals_pe12_'       => 100000,
          'vert_pe12_'          => 8000,
        }, # biotypes
      }, # genblast
      'projection' => {
        'logic_names' => {
          'project_transcripts' => 50000,
        }, # logic_names
        'biotypes' =>    {
          'projection'          => 50000,
        }, # biotypes
      }, # projection
      'realign' => {
        'logic_names' => {
          # Would actually prefer an upper limit on realign as opposed to a lower limit 
          'project_transcripts' => 20000,
          'realign'             => 1000,
        }, # logic_names
        'biotypes' =>    {
          'realign'             => 40000,
        }, # biotypes
      }, # realign
      'rnaseq_blast' =>  {
        # This one is an issue, logic names, counts are varied and one biotype is a
        # substring of the other. At the moment it's really just a check that theres'
        # some stuff in there
        'logic_names' => {
#          $species.'_merged_rnaseq' => 10000,
        }, # logic_names
        'biotypes' =>    {
          'rnaseq'              => 10000,
         }, # biotypes
      }, # rnaseq_blast
      'layer' => {
        'logic_names' => {
          'genblast'                 => 40000,
          'genblast_not_best'        => 10000,
          'project_transcripts'      => 20000,
#           $species.'_merged_rnaseq' => 10000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 10000,
          'primates_pe12_'      => 15000,
          'primates_pe3_'       => 1,
          'mammals_pe12_'       => 10000,
          'vert_pe12_'          => 1,
          'realign'             => 30000,
          'rnaseq'              => 10000,
        }, # biotypes
      }, # layer
      'genebuilder' => {
        'logic_names' => {
          'ensembl'             => 30000,
        }, # logic_names
        'biotypes' =>    {
          'human_pe12_'         => 2000,
          'primates_pe12_'      => 2000,
          'primates_pe3_'       => 0,
          'mammals_pe12_'       => 1000,
          'vert_pe12_'          => 0,
          'realign'             => 10000,
          'rnaseq'              => 2000,
        }, # biotypes
      }, # genebuilder
      'ncrna' => {
        'logic_names' => {
          'ncrna' => 2000,
        }, # logic_names
        'biotypes' =>    {
          'miRNA'               => 600,
          'misc_RNA'            => 400,
          'ribozyme'            => 0,
          'rRNA'                => 200,
          'scaRNA'              => 0,
          'snoRNA'              => 200,
          'snRNA'               => 400,
        }, # biotypes
      }, # ncrna
    }, # primates_basic
  }, # gene_db_checks
  );
  return $config{$key};
}

1;
