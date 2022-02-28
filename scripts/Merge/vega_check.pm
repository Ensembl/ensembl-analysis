#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

package Merge::vega_check;

use strict;
use warnings;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT_OK = qw(get_combos get_biotype_groups get_actions get_loutre_misc);

$| = 1;

#misc data for loutre QC
my %loutre_misc = (
  'disallowed_gene_trans_biotypes' => {
    'gene'       => ['polymorphic',
                     'pseudogene'],
    'transcript' => ['ambiguous_orf',
                     'disrupted_domain',
                     'non_coding',
                     'ncrna_host'],
  },
  'allowed_transcript_combos' => {},
);

my %allowed_combos = (
     'ensembl' => { # allowed gene-transcript biotypes combination before the merge
       '3prime_overlapping_ncrna'           => ['3prime_overlapping_ncrna',
                                                'retained_intron',
                                                'TEC',
                                                'artifact'],
       'antisense'                          => ['antisense',
                                                'retained_intron',
                                                'TEC',
                                                'artifact'],
       'artifact'                           => ['artifact'],
       'IG_gene'                            => ['IG_gene',
                                                'TEC',
                                                'artifact'],
       'IG_pseudogene'                      => ['IG_pseudogene',
                                                'TEC',
                                                'artifact'],
       'known_ncrna'                        => ['known_ncrna'],
       'lincRNA'                            => ['lincRNA',
                                                'retained_intron',
                                                'TEC',
                                                'artifact'],
       'macro_lncRNA'                       => ['macro_lncRNA',
                                                'processed_transcript'],
       'non_coding'                         => ['non_coding'],
       'polymorphic_pseudogene'             => ['nonsense_mediated_decay',
                                                'non_stop_decay',
                                                'polymorphic_pseudogene',
                                                'protein_coding',
                                                'processed_transcript',
                                                'retained_intron',
                                                'TEC',
                                                'artifact'
                                               ],
       'processed_pseudogene'               => ['processed_pseudogene',
                                                'TEC',
                                                'artifact'
                                               ],
       'processed_transcript'               => ['3prime_overlapping_ncrna',
                                                'antisense',
                                                'macro_lncRNA',
                                                'lincRNA',
                                                'processed_transcript',
                                                'retained_intron',
                                                'sense_intronic',
                                                'sense_overlapping',
                                                'TEC',
                                                'artifact'
                                               ],
       'protein_coding'                     => ['non_stop_decay',
                                                'nonsense_mediated_decay',
                                                'processed_transcript',
                                                'protein_coding',
                                                'retained_intron',
                                                'translated_unprocessed_pseudogene',
                                                'translated_processed_pseudogene',
                                                'TEC',
                                                'artifact'
                                               ],
       'rRNA'                               => ['rRNA',
                                                'TEC',
                                                'artifact'],
       'sense_intronic'                     => ['retained_intron',
                                                'sense_intronic',
                                                'TEC',
                                                'artifact'],
       'sense_overlapping'                  => ['retained_intron',
                                                'sense_overlapping',
                                                'TEC',
                                                'artifact'],
       'snorna'                             => ['snorna'],
       'TR_gene'                            => ['TR_gene',
                                                'TEC',
                                                'artifact'],
       'TR_pseudogene'                      => ['TR_pseudogene'],
       'transcribed_unitary_pseudogene'     => ['transcribed_unitary_pseudogene','processed_transcript','retained_intron','artifact','TEC'],
       'translated_unprocessed_pseudogene'  => ['translated_unprocessed_pseudogene',
                                                'TEC',
                                                'artifact',
                                                'processed_transcript'],
       'translated_processed_pseudogene'    => ['translated_processed_pseudogene',
                                                'TEC',
                                                'artifact',
                                                'processed_transcript'],
       'transcribed_unprocessed_pseudogene' => ['transcribed_unprocessed_pseudogene',
                                                'processed_transcript',
                                                'retained_intron',
                                                'TEC',
                                                'artifact'],
       'transcribed_processed_pseudogene'   => ['transcribed_processed_pseudogene',
                                                'processed_transcript',
                                                'retained_intron',
                                                'TEC',
                                                'artifact'],
       'TEC'                                => ['TEC'],
       'unitary_pseudogene'                 => ['processed_transcript',
                                                'retained_intron',
                                                'unitary_pseudogene',
                                                'TEC',
                                                'artifact'],
       'unprocessed_pseudogene'             => ['unprocessed_pseudogene'],
       'vaultRNA'                           => ['vaultRNA']
     },

     'ensembl_extension' => { # additional allowed gene-transcript biotypes combination after the merge
       'antisense'                          => ['miRNA',
                                                'misc_RNA',
                                                'snoRNA'],
       'lincRNA'                            => ['miRNA',
                                                'misc_RNA',
                                                'rRNA',
                                                'snoRNA',
                                                'snRNA'],
       'miRNA'                              => ['miRNA'],
       'misc_RNA'                           => ['misc_RNA'],
       'non_coding'                         => ['miRNA'],
       'processed_transcript'               => ['miRNA',
                                                'misc_RNA',
                                                'snoRNA',
                                                'antisense',
                                                'lincRNA',
                                                'non_coding',
                                                'sense_intronic',
                                                'sense_overlapping'
                                               ],
       'rRNA'                               => ['rRNA'],
       'sense_intronic'                     => ['miRNA',
                                                'snoRNA',
                                                'snRNA'],
       'sense_overlapping'                  => ['miRNA'],
       'snoRNA'                             => ['snoRNA'],
       'snRNA'                              => ['snRNA'],
       'protein_coding'                     => ['IG_C_gene',
                                                'IG_V_gene',
                                                'IG_LV_gene',
                                                'antisense',
                                                'lincRNA',
                                                'sense_intronic'
                                               ],
       'pseudogene'                         => ['disrupted_domain',
                                                'processed_pseudogene',
                                                'processed_transcript',
                                                'pseudogene',
                                                'retained_intron',
                                                'transcribed_processed_pseudogene',
                                                'transcribed_unprocessed_pseudogene',
                                                'translated_unprocessed_pseudogene',
                                                'translated_processed_pseudogene',
                                                'unitary_pseudogene',
                                                'unprocessed_pseudogene',
                                                'antisense',
                                                'lincRNA',
                                                'TEC',
                                                'artifact'
                                                ],

       'TR_V_gene'                          =>  ['TR_V_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'TR_C_gene'                          =>  ['TR_C_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'TR_D_gene'                          =>  ['TR_D_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'TR_J_gene'                          =>  ['TR_J_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'TR_V_pseudogene'                    =>  ['TR_V_pseudogene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'TR_C_pseudogene'                    =>  ['TR_C_pseudogene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'TR_D_pseudogene'                    =>  ['TR_D_pseudogene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'TR_J_pseudogene'                    =>  ['TR_J_pseudogene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_V_gene'                          =>  ['IG_V_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_V_pseudogene'                    =>  ['IG_V_pseudogene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_C_pseudogene'                    =>  ['IG_C_pseudogene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_D_pseudogene'                    =>  ['IG_D_pseudogene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_J_pseudogene'                    =>  ['IG_J_pseudogene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_LV_gene'                         =>  ['IG_LV_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_J_gene'                          =>  ['IG_J_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_C_gene'                          =>  ['IG_C_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],

       'IG_D_gene'                          =>  ['IG_D_gene',
                                                 'TEC',
                                                 'artifact',
                                                ],
      },

     'loutre' => { # biotype combinations in loutre
       protein_coding                       => ['protein_coding',
                                                'nonsense_mediated_decay',
                                                'non_stop_decay',
                                                'processed_transcript',
                                                'retained_intron',
                                                'artifact',
                                                'tec'],
       polymorphic_pseudogene               => ['protein_coding',
                                                'nonsense_mediated_decay',
                                                'non_stop_decay',
                                                'processed_transcript',
                                                'retained_intron',
                                                'artifact',
                                                'tec',
                                                'polymorphic_pseudogene'],
       processed_transcript                 => ['processed_transcript',
                                                'retained_intron',
                                                'antisense',
                                                'artifact',
                                                'lincrna',
                                                'sense_intronic',
                                                'sense_overlapping',
                                                qq(3'_overlapping_ncrna),
                                                'tr_gene',
                                                'lincrna',
                                                'macro_lncrna',
                                                'tec',
                                                'bidirectional_promoter_lncrna',
                                                'snorna',
                                                'mirna',
                                                'vaultRNA',
                                                'pirna',
                                                'rrna',
                                                'trna',
                                                'snrna',
                                                'sirna',
                                                'scrna',
                                                'tec'],
       processed_pseudogene                 => ['processed_pseudogene'],
       unprocessed_pseudogene               => ['unprocessed_pseudogene'],
       transcribed_processed_pseudogene     => ['transcribed_processed_pseudogene',
                                                'processed_transcript',
                                                'retained_intron',
                                                'artifact',
                                                'tec'],
       transcribed_unprocessed_pseudogene   => ['transcribed_unprocessed_pseudogene',
                                                'processed_transcript',
                                                'retained_intron',
                                                'artifact',
                                                'tec'],
       transcribed_unitary_pseudogene       => [ 'transcribed_unitary_pseudogene',
                                                 'processed_transcript',
                                                 'retained_intron',
                                                 'artifact',
                                                 'tec' ],
       unitary_pseudogene                   => ['unitary_pseudogene',
                                                'artifact',
                                                'tec'],
       translated_processed_pseudogene      => ['translated_processed_pseudogene',
                                                'processed_transcript',
                                                'retained_intron',
                                                'artifact'],
       translated_unprocessed_pseudogene    => ['translated_unprocessed_pseudogene'],
       tec                                  => ['tec'],
       ig_pseudogene                        => ['ig_pseudogene'],
       novel_transcript                     => ['processed_transcript'],
       ig_gene                              => ['ig_gene',],
       tr_gene                              => ['tr_gene'],
       tr_pseudogene                        => ['tr_pseudogene'],
       lincRNA                              => ['retained_intron',
                                                'lincrna',
                                                'artifact'],
       antisense                            => ['retained_intron',
                                                'antisense',
                                                'artifact'],
       sense_intronic                       => ['retained_intron',
                                                'sense_intronic',
                                                'artifact'],
       sense_overlapping                    => ['retained_intron',
                                                'sense_overlapping',
                                                'artifact'],
       qq(3'_overlapping_ncRNA)             => ['retained_intron',
                                                qq(3'_overlapping_ncrna),
                                                'artifact'],
     }
   );

my %biotype_groups = (
     'ensembl' => {
       'gene_coding'                        => ['IG_gene',
                                                'IG_V_gene',
                                                'IG_LV_gene',
                                                'IG_C_gene',
                                                'IG_D_gene',
                                                'IG_J_gene',
                                                'polymorphic_pseudogene',
                                                'protein_coding',
                                                'TR_gene',
                                                'TR_C_gene',
                                                'TR_D_gene',
                                                'TR_J_gene',
                                                'TR_V_gene'],
       'transcript_coding'                  => ['IG_gene',
                                                'IG_C_gene',
                                                'IG_V_gene',
                                                'IG_LV_gene',
                                                'IG_J_gene',
                                                'IG_D_gene',
                                                'non_stop_decay',
                                                'nonsense_mediated_decay',
                                                'polymorphic_pseudogene',
                                                'protein_coding',
                                                'TR_gene',
                                                'TR_C_gene',
                                                'TR_D_gene',
                                                'TR_J_gene',
                                                'TR_V_gene'],
       'gene_non_coding'                    => ['3prime_overlapping_ncrna',
                                                'antisense',
                                                'IG_pseudogene',
                                                'IG_V_pseudogene',
                                                'IG_C_pseudogene',
                                                'IG_D_pseudogene',
                                                'IG_J_pseudogene',
                                                'known_ncrna',
                                                'lincRNA',
                                                'macro_lncRNA',
                                                'miRNA',
                                                'misc_RNA',
                                                'non_coding',
                                                'processed_pseudogene',
                                                'processed_transcript',
                                                'pseudogene',
                                                'rRNA',
                                                'sense_intronic',
                                                'sense_overlapping',
                                                'snoRNA',
                                                'snorna',
                                                'snRNA',
                                                'transcribed_unitary_pseudogene',
                                                'transcribed_unprocessed_pseudogene',
                                                'transcribed_processed_pseudogene',
                                                'translated_processed_pseudogene',
                                                'translated_unprocessed_pseudogene',
                                                'TR_pseudogene',
                                                'TR_V_pseudogene',
                                                'TR_C_pseudogene',
                                                'TR_D_pseudogene',
                                                'TR_J_pseudogene',
                                                'unitary_pseudogene',
                                                'unprocessed_pseudogene',
                                                'TEC',
                                                'artifact',
                                                'vaultRNA'],
       'transcript_non_coding'              => ['3prime_overlapping_ncrna',
                                                'antisense',
                                                'disrupted_domain',
                                                'IG_pseudogene',
                                                'IG_V_pseudogene',
                                                'IG_C_pseudogene',
                                                'IG_D_pseudogene',
                                                'IG_J_pseudogene',
                                                'known_ncrna',
                                                'lincRNA',
                                                'macro_lncRNA',
                                                'miRNA',
                                                'misc_RNA',
                                                'non_coding',
                                                'processed_transcript',
                                                'processed_pseudogene',
                                                'pseudogene',
                                                'rRNA',
                                                'retained_intron',
                                                'sense_intronic',
                                                'sense_overlapping',
                                                'snoRNA',
                                                'snorna',
                                                'snRNA',
                                                'transcribed_unitary_pseudogene',
                                                'transcribed_unprocessed_pseudogene',
                                                'transcribed_processed_pseudogene',
                                                'translated_unprocessed_pseudogene',
                                                'translated_processed_pseudogene',
                                                'TR_pseudogene',
                                                'TR_V_pseudogene',
                                                'TR_C_pseudogene',
                                                'TR_D_pseudogene',
                                                'TR_J_pseudogene',
                                                'unitary_pseudogene',
                                                'unprocessed_pseudogene',
                                                'TEC',
                                                'artifact',
                                                'vaultRNA'] },
     'vega' => {
        'Protein_coding'                     => ['protein_coding'],
        'lncRNAs'                            => ['3prime_overlapping_ncRNA',
                                                 'sense_intronic',
                                                 'sense_overlapping',
                                                 'antisense',
                                                 'non_coding',
                                                 'macro_lncRNA',
                                                 'lincRNA',
                                                 'bidirectional_promoter_lncRNA'],
        'ncRNAs'                             => ['miRNA',
                                                 'rRNA',
                                                 'snoRNA',
                                                 'vaultRNA',
                                                 'piRNA',
                                                 'rRNA',
                                                 'tRNA',
                                                 'snRNA',
                                                 'siRNA',
                                                 'scRNA'],
        'Unclassified_processed_transcripts' => ['processed_transcript'],
        'Pseudogenes'                        => ['polymorphic_pseudogene',
                                                 'processed_pseudogene',
                                                 'unprocessed_pseudogene',
                                                 'transcribed_processed_pseudogene',
                                                 'transcribed_unprocessed_pseudogene',
                                                 'transcribed_unprocessed_pseudogene',
                                                 'transcribed_unitary',
                                                 'unitary_pseudogene',
                                                 'translated_processed_pseudogene',
                                                 'translated_unprocessed_pseudogene',
                                                 'IG_pseudogene',
                                                 'TR_pseudogene'],
        'IG'                                 => ['IG_gene',
                                                 'IG_pseudogene'],
        'TR'                                 => ['TR_gene',
                                                 'TR_pseudogene'],
        'Other'                              => ['TEC'],
        'artifact'                           => ['artifact'], }
   );

# my %actions = (
#   'ensembl' => {
#     'gene' => {
#       'TEC' => {
#         sub {
#           die('die horribly if this biotype is seen');
#           }
#       } } }, {
#     'transcript' => {
#       'polymorphic_pseudogene' => {
#         sub {
#           print('change to protein coding if no stop codons');
#           }
#       } } },
#   'havana' => {
#     'gene' => {
#       'polymorphic' => sub {
#         print('ask again how loutre continues to contain these');
#         }
#     } } );

sub get_combos {
  my $combos_key = shift;
  return $allowed_combos{$combos_key};
}

sub get_biotype_groups {
  my $biotype_groups = shift;
  return $biotype_groups{$biotype_groups};
}

sub get_loutre_misc {
  my $key = shift;
  return $loutre_misc{$key};
}

#sub get_actions {
#  my $actions_key = shift;
#  return $actions{$actions_key};
#}

1;
