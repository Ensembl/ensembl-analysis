=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::CopyGenes - 

=head1 SYNOPSIS

RunnableDB for copying genes from a source database to a target
database. By default all the genes in the database COPY_SOURCE_DB are copied into
the database COPY_TARGET_DB

=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::CopyGenes;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($given_source_db, $given_target_db, $given_biotype) = rearrange
    (['SOURCE_DB', 'TARGET_DB', 'BIOTYPE'], @args);

  #### Default values...
  $self->source_db_name("COPY_SOURCE_DB");
  $self->target_db_name("COPY_TARGET_DB");
  $self->biotype("");

  ### ...are over-ridden by parameters given in analysis table...
  my $ph = $self->parameters_hash;
  $self->source_db_name($ph->{-source_db}) if exists $ph->{-source_db};
  $self->target_db_name($ph->{-target_db}) if exists $ph->{-target_db};
  $self->biotype($ph->{-biotype}) if exists $ph->{-biotype};
  
  ### ...which are over-ridden by constructor arguments. 
  $self->source_db_name($given_source_db) if defined $given_source_db;
  $self->target_db_name($given_target_db) if defined $given_target_db;
  $self->biotype($given_biotype) if defined $given_biotype;
 
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
  my (@genes);
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
        $tr->get_all_ProteinFeatures;
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

  $self->output(\@genes);
}


##################################
sub write_output {
  my($self) = @_;
  
  my $target_db = $self->get_dbadaptor($self->target_db_name);
    
  my $g_adap = $target_db->get_GeneAdaptor;

  foreach my $g (@{$self->output}) {
    $g_adap->store($g);
  }
  
  return 1;
}


1;
