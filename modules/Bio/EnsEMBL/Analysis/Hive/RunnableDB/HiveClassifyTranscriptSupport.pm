=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: This are the default values:
                update_gene_biotype => 0,
                classification_type => 'standard',
                classification => undef,
                feature_type => 'protein_align_feature',
                skip_analysis => 0,
              feature_type can be either protein_align_feature or dna_align_feature
              classification_type can only be standard at the moment unless you fill
              the classification hash
              classification is the hash that will help build the sql queries. The standard
              classification is:
                '1' => [95, 95],
                '2' => [95, 80],
                '3' => [90, 80],
                '4' => [80, 60],
                '5' => [70, 60],
                '6' => [50, 50],
                '7' => [50, 25],
                '8' => [0, 0],
              The key will be concatenated to the biotype with an underscore, the first element
              is the coverage, the second element is the percentage of identity
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    update_gene_biotype => 0,
    classification_type => 'standard',
    classification => undef,
    feature_type => 'protein_align_feature',
    skip_analysis => 0,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Just checking that you have a 'target_db' and setting the classification hash
              if it's not provided or it is set to 'standard'
 Returntype : None
 Exceptions : Throws if 'target_db' is not set

=cut

sub fetch_input {
  my $self = shift;

  if($self->param('skip_analysis')) {
    $self->complete_early('Skip analysis flag is enabled, so no classification will be carried out');
  }

  my $target_db = $self->param('target_db');
  $self->hrdb_set_con($self->get_database_by_name('target_db'), 'target_db');

  if (!$self->param('classification')) {
    my %classification;
    if ($self->param('classification_type') eq 'standard_old') {
      %classification = (
        '1' => [95, 95],
        '2' => [95, 80],
        '3' => [90, 80],
        '4' => [80, 60],
        '5' => [70, 60],
        '6' => [50, 50],
        '7' => [50, 25],
        '8' => [0, 0],
      );
      $self->param('classification', \%classification);
    } elsif ($self->param('classification_type') eq 'standard' || $self->param('classification_type') eq 'long_read') {
      %classification = (
       '1' => [95, 90],
       '2' => [90, 80],
       '3' => [90, 60],
       '4' => [90, 40],
       '5' => [80, 20],
       '6' => [60, 20],
       '7' => [0, 0],
      );
      $self->param('classification', \%classification);
    } elsif ($self->param('classification_type') eq 'gifts') {
      %classification = (
        '1' => [100,100],
        '2' => [97,97],
        '3' => [95,95],
        '4' => [90,90],
        '5' => [80,90],
        '6' => [80,80],
        '7' => [70,70],
        '8' => [50,50],
        '9' => [0,0]
      );
      $self->param('classification', \%classification);
    } else {
      $self->throw('Unrecognised classification type, the default is "standard". Classification type passed in: '.$self->param('classification_type'));
    }
  }
}


=head2 run

 Arg [1]    : None
 Description: Do nothing, everything is done in write_output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my $self = shift;

  return 1;
}


=head2 write_output

 Arg [1]    : None
 Description: Backup the transcript table, if it already exists it keeps it as it was.
              Update the transcript biotype based on the classification hash.
              If 'update_gene_biotype' is set to 1, it also backs up the gene table and
              update the biotype
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

  my $dbc = $self->hrdb_get_con('target_db')->dbc;

  if($self->param('classification_type') eq 'long_read') {
    $self->set_long_read_biotypes($dbc);
  }

  my $sth = $dbc->prepare('CREATE table transcript_classify_bak like transcript');
  eval {
    $sth->execute();
  };
  if ($@) {
    $self->warning("Table transcript_classify_bak already exists, KEEPING it as it is");
  }
  else {
    $sth = $dbc->prepare('INSERT INTO transcript_classify_bak SELECT * FROM transcript');
    $sth->execute();
  }

  $sth = $dbc->prepare('UPDATE transcript t LEFT JOIN transcript_classify_bak tcb USING(transcript_id) '.
      'LEFT JOIN transcript_supporting_feature tsf USING(transcript_id) '.
      'LEFT JOIN '.$self->param('feature_type').' taf ON taf.'.$self->param('feature_type').'_id = tsf.feature_id '.
      'SET t.biotype = CONCAT(tcb.biotype, ?) WHERE hcoverage >= ? AND perc_ident >= ? '.
      'AND tsf.feature_type = "'.$self->param('feature_type').'"');
  my $classification = $self->param('classification');

  if ($self->param('classification_type') eq 'standard' or
      $self->param('classification_type') eq 'long_read' or
      $self->param('classification_type') eq 'fish' or
      $self->param('classification_type') eq 'gifts') {
    foreach my $key (sort {$b <=> $a} keys %{$classification}) {
      $sth->bind_param(1, '_'.$key);
      $sth->bind_param(2, $classification->{$key}->[0]);
      $sth->bind_param(3, $classification->{$key}->[1]);
      $sth->execute();
    }
  }

# I am putting this (and the gene backup code) because at the moment the core api doesn't allow us to get transcripts
# by biotype, only genes. Since LayerAnnotation and GeneBuilder work off gene biotypes, it's very hard to convert them
# to use transcript biotypes without an equivalent core api method. So for now they will continue to work off gene
# biotypes
  if($self->param('update_gene_biotype')) {
    $sth = $dbc->prepare('CREATE table gene_classify_bak like gene');
    eval {
      $sth->execute();
    };
    if ($@) {
      $self->warning("Table gene_classify_bak already exists, KEEPING it as it is");
    }
    else {
      $sth = $dbc->prepare('INSERT INTO gene_classify_bak SELECT * FROM gene');
      $sth->execute();
    }
    $sth = $dbc->prepare('UPDATE gene g LEFT JOIN transcript t USING(gene_id) SET g.biotype = t.biotype');
    $sth->execute();
  }
}


sub set_long_read_biotypes {
  my ($self,$dbc) = @_;

  # Note that this is somewhat different to other modules where often the transcript biotypes are the most informative
  # Here the gene biotypes are most informative. We want to update things to cdna where the biotype is one of our
  # accepted biotypes. We then update the transcript table via the gene table
  my $sth = $dbc->prepare('UPDATE gene SET biotype="cdna" WHERE biotype NOT LIKE "retained%" '.
                          'AND biotype NOT LIKE "nonsense_mediated_decay%"');
  $sth->execute();

  $sth = $dbc->prepare('UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype=gene.biotype');
  $sth->execute();
}


1;
