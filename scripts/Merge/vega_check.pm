#!/usr/bin/env perl

# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/scripts/Merge/vega_check.pm,v $
# $Revision: 1.2 $
package Merge::vega_check;

use strict;
use warnings;

$| = 1;

my %allowed_combos = (
     'ensembl' => { # allowed gene-transcript biotypes combination before the merge
       '3prime_overlapping_ncrna'           => ['3prime_overlapping_ncrna'],
       'antisense'                          => ['antisense','retained_intron'],
       'IG_gene'                            => ['IG_gene'],
       'IG_pseudogene'                      => ['IG_pseudogene'],
       'lincRNA'                            => ['lincRNA','retained_intron'],
       'non_coding'                         => ['non_coding'],
       'polymorphic_pseudogene'             => ['polymorphic_pseudogene',
                                                'protein_coding',
                                                'processed_transcript',
                                                'nonsense_mediated_decay',
                                                'retained_intron'],
       'processed_pseudogene'               => ['processed_pseudogene',
                                                'processed_transcript',
                                                'transcribed_processed_pseudogene'],
       'processed_transcript'               => ['ambiguous_orf',
                                                'processed_transcript',
                                                'retained_intron',
                                                'TEC'],
       'protein_coding'                     => ['protein_coding',
                                                'nonsense_mediated_decay',
                                                'non_stop_decay',
                                                'ambiguous_orf',
                                                'processed_transcript',
                                                'retained_intron',
                                                'TEC',
                                                '3prime_overlapping_ncrna'],
       'pseudogene'                         => ['disrupted_domain',
                                                'processed_pseudogene',
                                                'processed_transcript',
                                                'pseudogene',
                                                'retained_intron',
                                                'TEC'],
       'sense_intronic'                     => ['sense_intronic','retained_intron'],
       'sense_overlapping'                  => ['sense_overlapping'],
       'TR_gene'                            => ['TR_gene'],
       'TR_pseudogene'                      => ['TR_pseudogene'],
       'translated_processed_pseudogene'    => ['translated_processed_pseudogene']
       'unitary_pseudogene'                 => ['unitary_pseudogene',
                                                'processed_transcript',
                                                'retained_intron'],
       'unprocessed_pseudogene'             => ['processed_transcript',
                                                'retained_intron',
                                                'unprocessed_pseudogene',
                                                'transcribed_unprocessed_pseudogene'] },

     'ensembl_extension' => { # additional allowed gene-transcript biotypes combination after the merge
       'protein_coding'                     => ['IG_C_gene',
                                                'IG_V_gene'],
       'pseudogene'                         => ['unitary_pseudogene',
                                                'unprocessed_pseudogene',
                                                'transcribed_processed_pseudogene',
                                                'transcribed_unprocessed_pseudogene',
                                                'translated_processed_pseudogene',
                                                # temporarily allowed
                                                'antisense',
                                                'lincRNA'],
       # temporarily allowed
       'processed_transcript'               => ['antisense',
                                                'lincRNA',
                                                'non_coding',
                                                'sense_intronic',
                                                'sense_overlapping'],
       'protein_coding'                     => ['antisense',
                                                'lincRNA',
                                                'sense_intronic',
                                                'translated_processed_pseudogene'] },

     'havana' => { # allowed gene-transcript biotypes combination before the merge
       '' => [ '', '' ]
     } );

my %biotype_groups = (
     'ensembl' => {
       'gene_coding'                        => ['IG_gene',
                                                'polymorphic_pseudogene',
                                                'protein_coding',
                                                'TR_gene'],
       'transcript_coding'                  => ['IG_gene',
                                                'IG_C_gene',
                                                'IG_V_gene',
                                                'nonsense_mediated_decay'
                                                'non_stop_decay',
                                                'polymorphic_pseudogene',
                                                'protein_coding',
                                                'TR_gene'],
       'gene_non_coding'                    => ['3prime_overlapping_ncrna',
                                                'antisense',
                                                'IG_pseudogene',
                                                'lincRNA',
                                                'non_coding',
                                                'processed_pseudogene',
                                                'processed_transcript',
                                                'pseudogene',
                                                'sense_intronic',
                                                'sense_overlapping',
                                                'translated_processed_pseudogene'
                                                'TR_pseudogene',
                                                'unitary_pseudogene',
                                                'unprocessed_pseudogene'],
       'transcript_non_coding'              => ['3prime_overlapping_ncrna',
                                                'ambiguous_orf',
                                                'antisense',
                                                'disrupted_domain',
                                                'IG_pseudogene',
                                                'lincRNA',
                                                'non_coding',
                                                'processed_transcript',
                                                'processed_pseudogene',
                                                'pseudogene',
                                                'retained_intron',
                                                'sense_intronic',
                                                'sense_overlapping',
                                                'transcribed_unprocessed_pseudogene',
                                                'transcribed_processed_pseudogene',
                                                'translated_processed_pseudogene',
                                                'TR_pseudogene',
                                                'unitary_pseudogene',
                                                'unprocessed_pseudogene'] },

     'havana' => { }
);

my %actions = (
  'ensembl' => {
    'gene' => {
      'TEC' => {
        sub {
          die('die horribly if this biotype is seen');
          }
      } } }, {
    'transcript' => {
      'polymorphic_pseudogene' => {
        sub {
          print('change to protein coding if no stop codons');
          }
      } } },
  'havana' => {
    'gene' => {
      'polymorphic' => sub {
        print('ask again how loutre continues to contain these');
        }
    } } );

1;
