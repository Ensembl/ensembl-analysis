=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadUpdatecDNAs

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadUpdatecDNAs;

use strict;
use warnings;

use Bio::EnsEMBL::IO::Parser::Fasta;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    sequence_biotype => 'cdna',
    column_names => ['iid'],
    sequence_table_name => 'cdna_sequences',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: If 'strategy' is "update", it first retrieve all aligned and non aligned sequence
              from the cdna database 'old_cdna_db' and load the new sequences into the hive
              database. Then it writes a file of gene_ids 'retire_gene_file' to be able to
              remove from the 'new_cdna_db' any models based on sequences not present in 'cdna_file'.
              If 'strategy' is "complete", it simply loads all the sequences from 'cdna_file'.
 Returntype : None
 Exceptions : Throws if 'strategy' is not set
              Throws if it cannot open or close 'retire_gene_file' when in update mode
              Throws if it cannot query the databases when in update mode

=cut

sub fetch_input {
  my $self = shift;

  my $strategy = $self->param_required('strategy');

  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($self->param_required('cdna_file'));
  my $seq_file = new Bio::SeqIO( -file => $self->param_required('cdna_file'),
                                 -format => "Fasta",
                               );

  my $old_db = $self->get_database_by_name('old_cdna_db');
  my $new_db = $self->get_database_by_name('new_cdna_db');

  # create a hash that has the cdnas that have already been aligned previously as the key
  my %seen_cdna;

  # only need to get the hit names if we're only doing an update and not complete
  if ($strategy eq 'update') {
    my $sth = $old_db->dbc()->prepare('SELECT DISTINCT(hit_name) FROM dna_align_feature') or $self->throw("Sql error\n$!");
    $sth->execute();
    while ( my $cdna  = $sth->fetchrow_array ) {
      $seen_cdna {$cdna} = 1;
    }
    $sth = $old_db->dbc()->prepare('SELECT DISTINCT(identifier) FROM unmapped_object') or $self->throw("Sql error\n$!");
    $sth->execute();
    while ( my $cdna = $sth->fetchrow_array ) {
      $seen_cdna {$cdna} = 1;
    }
    $sth->finish();
  }
  my $biotype = $self->param('sequence_biotype');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param('sequence_table_name'));

  my @iids;
  while($parser->next()) {
    my ($header) = $parser->getHeader() =~ /(\S+).*/;
    my $sequence = $parser->getSequence();

    if (exists $seen_cdna{$header}) {
      delete $seen_cdna{$header};
    } else {
      my $db_row = [{
        'accession'  => $header,
        'seq'        => $sequence,
        'biotype'    => $biotype,
      }];
      $table_adaptor->store($db_row);

      push(@iids, $header);
    }
  }
  # now get a list of the retired genes
  if ($strategy eq 'update') {
    open (GENEIDFILE, '>'.$self->param('retire_gene_file')) || $self->throw('Could not open'.$self->param('retire_gene_file'));

    my $gene_id_query = $new_db->dbc()->prepare('SELECT DISTINCT(g.gene_id) FROM gene g, transcript t, transcript_supporting_feature tsf, dna_align_feature daf WHERE g.gene_id = t.gene_id AND t.transcript_id = tsf.transcript_id AND feature_type = "dna_align_feature" AND tsf.feature_id = daf.dna_align_feature_id AND daf.hit_name = ?') or $self->throw("Sql error\n$!");
    foreach my $key (keys %seen_cdna) {
      $gene_id_query->bind_param(1, $key);
      $gene_id_query->execute();

      while ( my $gene_id  = $gene_id_query->fetchrow_array ) {
        print GENEIDFILE $gene_id, "\n";
      }
    }
    close(GENEIDFILE) || $self->throw('Could not open'.$self->param('retire_gene_file'));
  }
  $self->param('inputlist', \@iids);
}

1;
