
=head1 NAME

ProteinAnnotation.pm - DESCRIPTION of Object

=head1 SYNOPSIS

RunnableDB for copying genes from a source database to a target
database. By default all the genes in the database COPY_SOURCE_DB are copied into
the database COPY_TARGET_DB

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
use Bio::EnsEMBL::Utils::Argument qw (rearrange);

@ISA = qw (
           Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild
           );

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($source_db, $target_db, $biotype) = rearrange
    (['SOURCE_DB', 'TARGET_DB', 'BIOTYPE'], @args);
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->source_db_name('COPY_SOURCE_DB');
  $self->target_db_name('COPY_TARGET_DB');
  #####################
  my $parameters_hash = $self->parameters_hash;
  $source_db = $parameters_hash->{-source_db} if(!$source_db);
  $target_db = $parameters_hash->{-target_db} if(!$target_db);
  $biotype = $parameters_hash->{-biotype} if(!$biotype);
  $self->source_db_name($source_db);
  $self->target_db_name($target_db);
  $self->biotype($biotype);
  return $self;
}


#getter/setters

sub source_db_name{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'source_db_name'} = $arg;
  }
  return $self->{'source_db_name'};
}

sub target_db_name{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'target_db_name'} = $arg;
  }
  return $self->{'target_db_name'};
}

sub biotype{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'biotype'} = $arg;
  }
  return $self->{'biotype'};
}

################################
sub fetch_input {
  my ($self) = @_;  

  my $source_db = $self->get_dbadaptor($self->source_db_name);
  my $slice = $source_db->get_SliceAdaptor->fetch_by_name($self->input_id);

  #
  # total paranoia: fetch everything up front
  #
  my (@genes, %protein_features);
  my @target_genes;
  if($self->biotype){
    @target_genes = @{$slice->get_all_Genes_by_type($self->biotype)};
  }else{
    @target_genes = @{$slice->get_all_Genes};
  }
  foreach my $g (@target_genes) {

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
  
  my $target_db = $self->get_dbadaptor($self->target_db_name);
    
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
