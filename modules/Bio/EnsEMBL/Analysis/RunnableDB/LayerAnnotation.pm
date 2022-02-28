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

Bio::EnsEMBL::Analysis::RunnableDB::LayerAnnotation - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::LayerAnnotation;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::LayerAnnotation;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($LAYERANNOTATION_CONFIG_BY_LOGIC);
           
  # other stuff here
  
  return $self;
}




###################################
sub fetch_input {
  my ($self) = @_;  

  foreach my $dbname (@{$self->SOURCEDB_REFS}) {
    my $dbh = $self->get_dbadaptor($dbname);
    my $slice = $dbh->get_SliceAdaptor->fetch_by_name($self->input_id);
    my $tlslice = $dbh->get_SliceAdaptor->fetch_by_region($slice->coord_system->name,
                                                          $slice->seq_region_name);

    foreach my $layer (@{$self->layers}) {
      foreach my $tp (@{$layer->biotypes}) {
        foreach my $g (@{$slice->get_all_Genes_by_type($tp, undef, 1)}) {
          $g = $g->transfer($tlslice);
          push @{$layer->genes}, $g;
        }
      }
    }
    $dbh->dbc->disconnect_when_inactive(0) ;
  }
}


#####################################
sub run {
  my ($self) = @_;

  my (@retained_genes);

  my @layers = @{$self->layers};

  for(my $i=0; $i < @layers; $i++) {
    my $layer = $layers[$i];

    if ($layer->genes) {
      my @layer_genes = sort {$a->start <=> $b->start} @{$layer->genes};

      my @compare_genes;
      
      my %filter_against = map { $_ => 1 } @{$layer->filter_against};

      for(my $j = $i-1; $j>=0; $j--) {
        if (exists $filter_against{$layers[$j]->id} and
            @{$layers[$j]->genes}) {
          push @compare_genes, @{$layers[$j]->genes};
        }
      }
      @compare_genes = sort {$a->start <=> $b->start} @compare_genes;

      if ($layer->filter_object) {
        @layer_genes = @{$layer->filter_object->filter(\@layer_genes, \@compare_genes)};
      } else {
        throw("No filter object");
      }
      
      if (not $layer->discard) {
        push @retained_genes, @layer_genes;
      }

      $layer->genes(\@layer_genes);
    }
  }
  
  $self->output(\@retained_genes);
}


##################################
sub write_output {
  my($self) = @_;
  
  my $target_db = $self->get_dbadaptor($self->TARGETDB_REF);    
  my $g_adap = $target_db->get_GeneAdaptor;

  # fully loading gene is required for the store to work
  # reliably. However, fully loading all genes to be stored
  # up front is expensive in memory. Therefore, load, store and
  # discard one gene at a time

  my $total = 0;
  my $fails = 0;
  my $out = $self->output;
  while (my $g = shift @$out) {
    fully_load_Gene($g);

    # Putting this in to stop the wrong dbIDs being assigned
    empty_Gene($g);
    eval {
      $g_adap->store($g);
    };
    
    if ($@) {
      $self->warning("Unable to store gene".$g->dbID()." ".$g->stable_id()."\n");
      print STDERR "$@\n";
      $fails++;  
    }
    $total++;
  }
  
  if ($fails > 0) {
    throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }

  return 1;
}

####################################
sub layers {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_runnabledb_layers} = $val;
  }

  return $self->{_runnabledb_layers};
}



####################################
sub read_and_check_config {
  my ($self, $hash) = @_;

  $self->SUPER::read_and_check_config($hash);

  # check that all relevant vars have been defined
  foreach my $i (qw(LAYERS
                    SOURCEDB_REFS
                    TARGETDB_REF
                    FILTER)) {
    throw("You must define $i in config")
        if not $self->$i;
  }

  # check types
  throw("Config var LAYERS must be an array reference")
      if ref($self->LAYERS) ne "ARRAY";
  throw("Config var SOURCEDB_REFS must be an array reference")
      if ref($self->SOURCEDB_REFS) ne "ARRAY";
  throw("Config var TARGETDB_REF must be an scalar")
      if ref($self->TARGETDB_REF);
  throw("Config var FILTER must be an scalar")
      if ref($self->FILTER);

  my $filter;
  if ($self->FILTER) {
    $self->require_module($self->FILTER);
    $filter = $self->FILTER->new;
  }     

  my (%biotypes, %layer_ids, @layers);

  # finally, check the integrity of the LAYERS
  foreach my $el (@{$self->LAYERS}) {
    throw("Elements of LAYERS must be hash references")
        if ref($el) ne "HASH";
    
    throw("Each element of LAYERS must contain a layer ID")
        if not exists $el->{ID};

    throw("LAYER " . $el->{ID} . " should a list of BIOTYPES")
        if not exists $el->{BIOTYPES} or ref($el->{BIOTYPES}) ne "ARRAY";

    my $layer_id = $el->{ID};
    my @biotypes = @{$el->{BIOTYPES}};
    my $discard = 0;
    my (@filter_against);

    if (exists $el->{DISCARD} and $el->{DISCARD}) {
      $discard = 1;
    }
    foreach my $tp (@biotypes) {
      if (exists $biotypes{$tp}) {
        throw("biotype $tp occurs more than once");
      }
      $biotypes{$tp} = 1;
    }
    if (exists $el->{FILTER_AGAINST}) {
      throw("In layer $layer_id FILTER_AGAINST must contain a list of layer ids")
          if ref($el->{FILTER_AGAINST}) ne "ARRAY";      
      @filter_against = @{$el->{FILTER_AGAINST}};

      foreach my $id (@filter_against) {
        throw("In FILTER_AGAINST in layer $layer_id, '$id' is not the name ". 
              "of a higher level layer")
            if not exists $layer_ids{$id};
      }
    }      
    

    push @layers, Bio::EnsEMBL::Analysis::RunnableDB::LayerAnnotation::Layer
        ->new(-id => $layer_id,
              -discard => $discard,
              -biotypes =>  \@biotypes,
              -filter_object  => $filter,
              -filter_against => \@filter_against);

    $layer_ids{$layer_id} = 1;
  }

  $self->layers(\@layers);
}


################################################


sub LAYERS {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_layers} = $val;
  }

  return $self->{_layers};
}


sub FILTER {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_filter_object_name} = $val;
  }

  return $self->{_filter_object_name};
}


sub SOURCEDB_REFS {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_inputdb_names} = $val;
  }

  return $self->{_inputdb_names};
}


sub TARGETDB_REF {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_outputdb_name} = $val;
  }

  return $self->{_outputdb_name};
}

##############################################

package Bio::EnsEMBL::Analysis::RunnableDB::LayerAnnotation::Layer;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);

sub new {
  my ($class, @args) = @_;

  my $self = bless {}, $class;

  my ($id,
      $discard,
      $genes,
      $filter_object,
      $filter_against, 
      $biotypes) = rearrange
          (['ID','DISCARD','GENES','FILTER_OBJECT','FILTER_AGAINST','BIOTYPES'],@args);

  throw("Layers must have an id") if not defined $id;
  throw("Layers must have a list of biotypes") if not defined $biotypes;
  $discard = 0 if not defined $discard;
  $genes = [] if not defined $genes;
  $filter_against = [] if not defined $filter_against;

  $self->id($id);
  $self->filter_object($filter_object) if defined $filter_object;
  $self->filter_against($filter_against);
  $self->biotypes($biotypes);
  $self->discard($discard);
  $self->genes($genes);

  return $self;
}


sub id{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_layer_id} = $value;
  }
  return $self->{_layer_id};
}


sub filter_against{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_layer_filter_against} = $value;
  }
  return $self->{_layer_filter_against};
}


sub filter_object{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_layer_filter} = $value;
  }
  return $self->{_layer_filter};
}


sub discard{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_layer_discard} = $value;
    }
  return $self->{_layer_discard};
}


sub biotypes{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_layer_biotypes} = $value;
    }
  return $self->{_layer_biotypes};
}


sub genes{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_layer_genes} = $value;
  }
  return $self->{_layer_genes};
}



1;
