#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  my $target_db = $self->param('target_db');

  unless($target_db) {
    $self->throw("target_db not passed into param hash. Expected an input id of the format:\n".
                 "{'iid' => 'target_db_connection_hash'}");
  }

  my $target_dba = $self->hrdb_get_dba($target_db);
  $self->hrdb_set_con($target_dba,'target_db');

  my $classification_type = $self->param('classification_type');
  unless($classification_type) {
    $self->warning("The classification_type param has not been passed in so defaulting to 'standard'");
    $classification_type = 'standard';
  }

  $self->classification_type($classification_type);

  return 1;
}

sub run {
  my $self = shift;

  my $type = $self->classification_type();
  $self->run_classification($type);

  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub run_classification {
  my ($self,$type) = @_;

  my $update_gene_biotype = 0;
  if($self->param('update_gene_biotype')) {
    $update_gene_biotype = 1;
  }
  
  my $update_rnaseq_biotype = 0;
  if($self->param('update_rnaseq_biotype')) {
    $update_rnaseq_biotype = 1;
  }

  my $target_db = $self->hrdb_get_con('target_db');

  if($update_rnaseq_biotype) {
      $sth_update = $target_db->dbc->prepare('update gene set biotype="rnaseq" where biotype in ("single","best","other_merged");');
      $sth_update->execute();

      $sth_update = $target_db->dbc->prepare('update gene set biotype="rnaseq_tissue" where biotype != "rnaseq";');
      $sth_update->execute();      

      $sth_update = $target_db->dbc->prepare('update transcript left join gene using(gene_id) set transcript.biotype=gene.biotype;');
      $sth_update->execute();  
  }	

  my $sth_create =  $target_db->dbc->prepare('CREATE table transcript_classify_bak like transcript');
  $sth_create->execute();

  my $sth_copy =  $target_db->dbc->prepare('INSERT INTO transcript_classify_bak SELECT * FROM transcript');
  $sth_copy->execute();

  $sth_create =  $target_db->dbc->prepare('CREATE table gene_classify_bak like gene');
  $sth_create->execute();

  $sth_copy =  $target_db->dbc->prepare('INSERT INTO gene_classify_bak SELECT * FROM gene');
  $sth_copy->execute();

  if($type eq 'standard') {
    my $sth_update = $target_db->dbc->prepare('UPDATE transcript t LEFT JOIN transcript_classify_bak tcb using(transcript_id) '.
                                              'LEFT JOIN transcript_supporting_feature tsf using(transcript_id) '.
                                              'LEFT JOIN protein_align_feature paf on paf.protein_align_feature_id=tsf.feature_id '.
                                              'SET t.biotype=concat(tcb.biotype,"_0") WHERE hcoverage >= 0 AND perc_ident >= 0 '.
                                              'AND tsf.feature_type="protein_align_feature"');
    $sth_update->execute();


    $sth_update = $target_db->dbc->prepare('UPDATE transcript t LEFT JOIN transcript_classify_bak tcb using(transcript_id) '.
                                           'LEFT JOIN transcript_supporting_feature tsf using(transcript_id) '.
                                           'LEFT JOIN protein_align_feature paf on paf.protein_align_feature_id=tsf.feature_id '.
                                           'SET t.biotype=concat(tcb.biotype,"_50") WHERE hcoverage >= 50 AND perc_ident >= 50 '.
                                           'AND tsf.feature_type="protein_align_feature"');
    $sth_update->execute();


    $sth_update = $target_db->dbc->prepare('UPDATE transcript t LEFT JOIN transcript_classify_bak tcb using(transcript_id) '.
                                           'LEFT JOIN transcript_supporting_feature tsf using(transcript_id) '.
                                           'LEFT JOIN protein_align_feature paf on paf.protein_align_feature_id=tsf.feature_id '.
                                           'SET t.biotype=concat(tcb.biotype,"_80") WHERE hcoverage >= 80 AND perc_ident >= 80 '.
                                           'AND tsf.feature_type="protein_align_feature"');
    $sth_update->execute();


    $sth_update = $target_db->dbc->prepare('UPDATE transcript t LEFT JOIN transcript_classify_bak tcb using(transcript_id) '.
                                           'LEFT JOIN transcript_supporting_feature tsf using(transcript_id) '.
                                           'LEFT JOIN protein_align_feature paf on paf.protein_align_feature_id=tsf.feature_id '.
                                           'SET t.biotype=concat(tcb.biotype,"_95") WHERE hcoverage >= 95 AND perc_ident >= 95 '.
                                           'AND tsf.feature_type="protein_align_feature"');
    $sth_update->execute();

    # I am putting this (and the gene backup code) because at the moment the core api doesn't allow us to get transcripts
    # by biotype, only genes. Since LayerAnnotation and GeneBuilder work off gene biotypes, it's very hard to convert them
    # to use transcript biotypes without an equivalent core api method. So for now they will continue to work off gene
    # biotypes
    if($update_gene_biotype) {
      $sth_update = $target_db->dbc->prepare('UPDATE gene g LEFT JOIN transcript t using(gene_id) set g.biotype=t.biotype');
      $sth_update->execute();
    }
  } else {
    $self->throw("Unrecognised classification type, the default is 'standard'. Classification type passed in: ".$type);
  }

}

sub classification_type {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_classification_type',$val);
  }
  return($self->param('_classification_type'));

}


1;
