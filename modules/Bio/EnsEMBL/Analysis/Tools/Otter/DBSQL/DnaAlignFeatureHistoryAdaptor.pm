# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

# Perl module for Bio::EnsEMBL::Analysis::Tools::Otter::DBSQL::DnaAlignFeatureHistoryAdaptor
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Otter::DBSQL::DnaAlignFeatureHistoryAdaptor

=head1 SYNOPSIS

fetches the DnaAlignFeatureHistory objects from an Otter database

=head1 DESCRIPTION


=head1 CONTACT

  Post questions/comments to the EnsEMBL development list:
  http://lists.ensembl.org/mailman/listinfo/dev

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Tools::Otter::DBSQL::DnaAlignFeatureHistoryAdaptor;

use warnings ;
use Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::Cache;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub fetch_all_by_seq_region {
  my $self = shift;
  my ($seq_region_id, $analysis_id, $daf_id) = @_;

  my ($daf_history, $dbID);
  my $rowHashRef;

  my $sth = $self->prepare(
            "SELECT align_feature_history_id, ".
            "seq_region_id, analysis_id, align_feature_id_start, ".
            "align_feature_id_end, db_version, date ".
            "FROM dna_align_feature_history ".
            "WHERE seq_region_id = ? ".
            "AND analysis_id = ? ".
            "AND align_feature_id_start <= ? ".
            "AND align_feature_id_end >= ? "
            );
  
  $sth->bind_param(1,$seq_region_id ,SQL_INTEGER);
  $sth->bind_param(2,$analysis_id ,SQL_SMALLINT);
  $sth->bind_param(3,$daf_id ,SQL_INTEGER);
  $sth->bind_param(4,$daf_id ,SQL_INTEGER);
  $sth->execute;

  while ($rowHashRef = $sth->fetchrow_hashref) {
    my $daf_history = $self->_objFromHashref($rowHashRef);
    foreach my $num ($daf_history->align_feature_id_start .. $daf_history->align_feature_id_end) {
      my $key = uc(join(':', $daf_history->seq_region_id,$daf_history->analysis->dbID,$daf_id));
      $self->{_dafh_cache}->{$key} = $daf_history;
      #$self->{_dafh_cache}->{$key} = $daf_history;
    }
  }
  $sth->finish;
  
  my @daf_histories = values %{$self->{_dafh_cache}};
  return \@daf_histories;
}

sub clear_cache {
  my ($self) = @_;
  $self->{'_dafh_cache'} = ();
}

sub _tables {
  my $self = shift;
  return (['dna_align_feature_history' , 'dafh']);
}

sub _columns {
  my $self = shift;
  return ('dafh.align_feature_history_id','dafh.seq_region_id',
          'dafh.analysis_id','dafh.align_feature_id_start',
          'dafh.align_feature_id_end','dafh.db_version','dafh.date');
}

sub fetch_by_DnaAlignFeature_info {
  my $self = shift;
  my ($dna_align_feature_id, $seq_region_id, $analysis_id) = @_;

  # each row in dna_align_feature_history is uniquely identified
  my $key = uc(join(':', $seq_region_id,$analysis_id,$dna_align_feature_id));

  #will only use feature_cache if hasn't been no_cache attribute set
  if (!defined $self->db->no_cache || !$self->db->no_cache){
    if (exists($self->{'_dafh_cache'}->{$key}) && defined $self->{'_dafh_cache'}->{$key}) {
      #print STDERR "  Found $key in cache\n";
      return $self->{'_dafh_cache'}->{$key};
    } else {
      #print STDERR "Not found $key , load cache\n";
    }
  }

  #
  # load the cache with some extras
  #
  $self->fetch_all_by_seq_region($seq_region_id, $analysis_id, $dna_align_feature_id);

  my $original_key = uc(join(':', $seq_region_id,$analysis_id,$dna_align_feature_id));
  my $dafh = $self->{'_dafh_cache'}->{$original_key};

  #will only use feature_cache when set attribute no_cache in DBAdaptor
  if (!defined $self->db->no_cache || !$self->db->no_cache){
    # let's try to be clever and make keys for all dafs associated with this dafh
    foreach my $num ($dafh->align_feature_id_start .. $dafh->align_feature_id_end) {
      my $k =  uc(join(':', $seq_region_id,$analysis_id,$num));
      if (!exists $self->{'_dafh_cache'}->{$k} || !defined $self->{'_dafh_cache'}->{$k}) {
        #print STDERR "  created key $k... size hash is ".scalar(keys %{$self->{'_dafh_cache'}})."\n";
        $self->{'_dafh_cache'}->{$k} = $dafh;
      }
    }
  }
  return $dafh;
}

=head2 _objFromHashref

  Arg [1]    : hashref $rowHash
  Description: Private helper function generates an Analysis object from a
               mysql row hash reference.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AnalsisAdaptor::fetch_* methods
  Status     : Stable

=cut

sub _objFromHashref {
  my $self = shift;
  my $rowHash = shift;

  my $aa = $self->db->get_AnalysisAdaptor();
  my $analysis = $aa->fetch_by_dbID($rowHash->{analysis_id});

  my $dna_align_feature_history = Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory->new(
        -id                     => $rowHash->{align_feature_history_id},
        -seq_region_id          => $rowHash->{seq_region_id},
        -analysis               => $analysis,
        -align_feature_id_start => $rowHash->{align_feature_id_start},
        -align_feature_id_end   => $rowHash->{align_feature_id_end},
        -db_version             => $rowHash->{db_version},
        -date                   => $rowHash->{date},
        -adaptor                => $self, 
    );

  return $dna_align_feature_history;
}


=head2 _objFromHashref

  Arg [1]    : $sth 
  Description: Private helper function generates an DnaAlignFeatureHistory object
  Returntype : Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AnalsisAdaptor::fetch_* methods
  Status     : Stable

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my $aa = $self->db->get_AnalysisAdaptor();
  my @out;

  my ($align_feature_history_id, $seq_region_id, $analysis_id, $align_feature_id_start,
      $align_feature_id_end, $db_version, $date);
  $sth->bind_columns(\$align_feature_history_id, \$seq_region_id, \$analysis_id, \$align_feature_id_start,
      \$align_feature_id_end, \$db_version, \$date);

  while($sth->fetch()) {
    my $analysis = $aa->fetch_by_dbID($analysis_id);

    push @out, Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory->new(
        -id => $align_feature_history_id,
        -seq_region_id          => $seq_region_id,
        -analysis               => $analysis,
        -align_feature_id_start => $align_feature_id_start,
        -align_feature_id_end   => $align_feature_id_end,
        -db_version             => $db_version,
        -date                   => $date,
        -adaptor                => $self,
              );
  }
  return \@out;
}


1;
