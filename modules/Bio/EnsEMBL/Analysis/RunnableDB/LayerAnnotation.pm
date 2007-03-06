
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

    foreach my $layer (@{$self->LAYERS}) {
      foreach my $tp (@{$layer->{BIOTYPES}}) {
        foreach my $g (@{$slice->get_all_Genes_by_type($tp)}) {
          $g = $g->transfer($tlslice);
          push @{$layer->{GENES}}, $g;
        }
      }
    }
  }

}


#####################################
sub run {
  my ($self) = @_;

  my (@retained_genes);

  my @layers = @{$self->LAYERS};

  for(my $i=0; $i < @layers; $i++) {
    my $layer = $layers[$i];

    if (exists $layer->{GENES}) {
      my @layer_genes = sort {$a->start <=> $b->start} @{$layer->{GENES}};

      if (exists $layer->{FILTER_AGAINST}) {
        my @compare_genes;

        my %filter_against = map { $_ => 1 } @{$layer->{FILTER_AGAINST}};
        for(my $j = $i-1; $j>=0; $j--) {
          if (exists $filter_against{$layers[$j]->{ID}} and
              exists $layers[$j]->{GENES}) {
            push @compare_genes, @{$layers[$j]->{GENES}};
          }
        }
        @compare_genes = sort {$a->start <=> $b->start} @compare_genes;

        if (exists $layer->{FILTER}) {
          @layer_genes = @{$layer->{FILTER}->filter(\@layer_genes, \@compare_genes)};
        } else {
          @layer_genes = @{$self->generic_filter(\@layer_genes, \@compare_genes)};
        }
      }
      
      if (not $layer->{DISCARD}) {
        push @retained_genes, @layer_genes;
      }

      $layer->{GENES} = \@layer_genes;
    }
  }

  $self->output(\@retained_genes);
}


##################################
sub write_output {
  my($self) = @_;
  
  my $target_db = $self->get_dbadaptor($self->TARGETDB_ID);    
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
    
    if (not $exon_overlap) {
      push @filtered, $obj;
    }
  }

  return \@filtered;
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

  my (%biotypes, %layer_ids);

  # finally, check the integrity of the LAYERS
  foreach my $el (@{$self->LAYERS}) {
    throw("Elements of LAYERS must be hash references")
        if ref($el) ne "HASH";
    
    throw("Each element of LAYERS must contain a layer ID")
        if not exists $el->{ID};

    throw("LAYER " . $el->{ID} . " should a list of BIOTYPES")
        if not exists $el->{BIOTYPES} or ref($el->{BIOTYPES}) ne "ARRAY";
  }

  foreach my $el (@{$self->LAYERS}) {
    if (not exists $el->{DISCARD}) {
      $el->{DISCARD} = 0;
    }

    foreach my $tp (@{$el->{BIOTYPES}}) {
      if (exists $biotypes{$tp}) {
        throw("In layer ".$el->{ID} . ", biotype $tp occurs more than once");
      }
      $biotypes{$tp} = 1;
    }
    if (exists $el->{FILTER_AGAINST}) {
      throw("In layer ".$el->{ID} . ", FILTER_AGAINST must contain a list of layer ids")
          if ref($el->{FILTER_AGAINST}) ne "ARRAY";      
      foreach my $id (@{$el->{FILTER_AGAINST}}) {
        throw("In FILTER_AGAINT in layer ". $el->{ID} . ", '$id' is not the name ". 
              "of a higher level layer")
            if not exists $layer_ids{$id};
      }
      
      if (exists $el->{FILTER}) {
        $self->require_module($el->{FILTER});
        $el->{FILTER} = $el->{FILTER}->new();;
      }     
    }

    $layer_ids{$el->{ID}} = 1;
  }
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


1;
