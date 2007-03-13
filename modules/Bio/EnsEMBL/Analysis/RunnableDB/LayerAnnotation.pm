
=head1 NAME

LayerAnnotation - DESCRIPTION of Object

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::LayerAnnotation;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::Analysis::Config::GeneBuild::LayerAnnotation;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);

@ISA = qw (
           Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild
           );

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
        foreach my $g (@{$slice->get_all_Genes_by_type($tp)}) {
          $g = $g->transfer($tlslice);
          push @{$layer->genes}, $g;
        }
      }
    }
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
        @layer_genes = @{$self->generic_filter(\@layer_genes, \@compare_genes)};
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

  my $out = $self->output;
  while(my $g = shift @$out) {
    fully_load_Gene($g);
    $g_adap->store($g);
  }
  
  return 1;
}

#####################################
sub generic_filter {
  my ($self, $these, $others) = @_;

  # interference is judged by overlap at exon level
  # assumption is that @others is sorted by gene start

  my @filtered;

  my $cur_idx = 0;

  foreach my $obj (@$these) {
    my (@genomic_overlap, $left_bound);

  
    for(my $i=$cur_idx; $i < @$others; $i++) {
      my $o_obj = $others->[$i];

      if ($o_obj->end >= $obj->start and not defined $left_bound) {
        $left_bound = $i;
      }
      
      if ($o_obj->end < $obj->start) {
        next;
      } elsif ($o_obj->start > $obj->end) {
        last;
      } else {
        push @genomic_overlap, $o_obj;
      } 
    }
    
    $cur_idx = $left_bound if defined $left_bound;

    my $exon_overlap = 0;
    if (@genomic_overlap) {
      my @exons = @{$obj->get_all_Exons};
      OG: foreach my $o_obj (@genomic_overlap) {
        foreach my $oe (@{$o_obj->get_all_Exons}) {
          foreach my $e (@exons) {
            if ($oe->strand == $e->strand and 
                $oe->end >= $e->start and
                $oe->start <= $e->end) {  
              $exon_overlap = 1;
              last OG;
            }
          }
        }
      }
    }      

    if (not $exon_overlap) {
      push @filtered, $obj;
    }
  }

  return \@filtered;
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
                    TARGETDB_REF)) {
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
    my ($filter, @filter_against);

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
    

    if (exists $el->{FILTER}) {
      $self->require_module($el->{FILTER});
      $filter = $el->{FILTER}->new;
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
