
=head1 NAME

ProteinAnnotation.pm - DESCRIPTION of Object

=head1 SYNOPSIS

RunnableDB for copying genes from a source database to a target
database. Very simple at the moment - all genes in source are
copies to target, and source and target database details are 
defined in COPY_SOURCE_DB and COPY_TARGET_DB entries in
Config/GeneBuild/Databases.pm.

Note: make sure analysis tables are synchronised. Not essential,
but ensures that analysis ids are preserved in the copy

=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::CopyGenes;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB
           Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild
           );

################################
sub fetch_input {
  my ($self) = @_;  

  my $source_db = $self->get_dbadaptor('COPY_SOURCE_DB');
  my $slice = $source_db->get_SliceAdaptor->fetch_by_name($self->input_id);

  #
  # total paranoia: fetch everything up front
  #
  my (@genes, %protein_features);

  foreach my $g (@{$slice->get_all_Genes}) {

    foreach my $t (@{$g->get_all_Transcripts}) {
      foreach my $e (@{$t->get_all_Exons}) {
        $e->get_all_supporting_features;
        $e->stable_id;
      }
      my $tr = $t->translation;
      if (defined $tr) {
        $tr->stable_id;
        $tr->get_all_Attributes;
        $tr->get_all_DBEntries;

        # adaptor cascade for gene does not handle protein
        # features. In order to do that here, we have to keep
        # the the translation object so that we can obtain
        # its new dbID after storage
        $protein_features{$tr->dbID} = [$tr];
        foreach my $pf (@{$tr->get_all_ProteinFeatures}) {
          push @{$protein_features{$tr->dbID}}, $pf;
        }
      }
      $t->stable_id;
      $t->get_all_supporting_features;
      $t->get_all_DBEntries;
      $t->display_xref;
      $t->get_all_Attributes;
    }
    $g->stable_id;
    $g->get_all_DBEntries;
    $g->display_xref;
    $g->get_all_Attributes;
  
    push @genes, $g;
  }

  $self->protein_features(\%protein_features);
  $self->output(\@genes);
}


##################################
sub write_output {
  my($self) = @_;
  
  my $target_db = $self->get_dbadaptor('COPY_TARGET_DB');
    
  my $g_adap = $target_db->get_GeneAdaptor;
  my $p_adap = $target_db->get_ProteinFeatureAdaptor;

  foreach my $g (@{$self->output}) {
    $g_adap->store($g);
  }

  # At this point, all translations should have new dbIDs.
  # We can now store the protein features
  foreach my $old_dbid (keys %{$self->protein_features}) {
    my ($trans, @feats) = @{$self->protein_features->{$old_dbid}};

    my $new_dbid = $trans->dbID;
    foreach my $f (@feats) {
      $p_adap->store($f, $new_dbid);
    }
  }
  
  return 1;
}



##########

sub protein_features {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_protein_feature_map} = $val;
  }

  return $self->{_protein_feature_map};
}

1;
