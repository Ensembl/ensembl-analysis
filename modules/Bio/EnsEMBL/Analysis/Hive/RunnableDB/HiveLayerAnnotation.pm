=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLayerAnnotation;

use warnings ;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

###################################
sub fetch_input {
  my ($self) = @_;

  if($self->param('skip_analysis')) {
    $self->complete_early('Skip check flag is enabled, so no check will be carried out');
  }

  # This call will set the config file parameters. Note this will set REFGB (which overrides the
  # value in $self->db and OUTDB
  foreach my $i (qw(LAYERS SOURCEDB_REFS TARGETDB_REF FILTER)) {
     $self->throw("You must define $i in config") unless ($self->param_required($i));
  }
  $self->create_analysis;

  my $target_dba = $self->hrdb_get_dba($self->TARGETDB_REF);
  my $dna_dba;
  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $target_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($target_dba,'target_db');

  my $found_input_genes = 0;
  foreach my $input_db (@{$self->SOURCEDB_REFS}) {
    my $dba = $self->hrdb_get_dba($input_db);
    if($dna_dba) {
      $dba->dnadb($dna_dba);
    }

    my $slice = $dba->get_SliceAdaptor->fetch_by_name($self->param('iid'));

    # Note I have re-written the below to no longer use get_all_Genes_by_Biotype as it causes a massive strain
    # on cpu usage of the servers when many jobs are running at once
    # The below is inefficient in terms of constantly looping through genes, but the overhead is negligable
    # in comparison to the massive savings by not straining the servers. Still should rewrite to be more
    # efficient and possibly consider changing our layering data structures in the static config to reflect this
    my $genes = $slice->get_all_Genes();
    foreach my $layer (@{$self->layers}) {
      foreach my $tp (@{$layer->biotypes}) {
        foreach my $g (@{$genes}) {
          if($g->biotype eq $tp) {
            $found_input_genes = 1;
            push @{$layer->genes}, $g;
	  }
        }
      }
    }
    $dba->dbc->disconnect_if_idle();
  }

  # If there are no input genes then finish and don't flow
  unless($found_input_genes) {
    $self->input_job->autoflow(0);
    $self->complete_early('No genes to process');
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
        $self->throw("No filter object");
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

  my $target_dba = $self->hrdb_get_con('target_db');
  my $g_adap = $target_dba->get_GeneAdaptor;

  # fully loading gene is required for the store to work
  # reliably. However, fully loading all genes to be stored
  # up front is expensive in memory. Therefore, load, store and
  # discard one gene at a time
  my $total = 0;
  my $fails = 0;
  foreach my $g (@{$self->output}) {
    fully_load_Gene($g);
    empty_Gene($g);
    eval {
      $g_adap->store($g);
    };

    if ($@) {
      $self->warning('Unable to store gene'.$g->seq_region_name.' '.$g->seq_region_start.' '.$g->seq_region_end.' '.$g->seq_region_strand.' '.$g->biotype."\n$@");
      ++$fails;
    }
    ++$total;
  }

  if ($fails > 0) {
    $self->throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }

  return 1;
}

####################################
sub layers {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('_runnabledb_layers', $val);
  }
  elsif ($self->LAYERS and !$self->param_is_defined('_runnabledb_layers')) {
    my $filter;
    if ($self->FILTER) {
      $self->require_module($self->FILTER);
      $filter = $self->FILTER->new;
    }

    my (%biotypes, %layer_ids, @layers);

    # finally, check the integrity of the LAYERS
    foreach my $el (@{$self->LAYERS}) {
      $self->throw("Elements of LAYERS must be hash references")
          if ref($el) ne "HASH";

      $self->throw("Each element of LAYERS must contain a layer ID")
          if not exists $el->{ID};

      $self->throw("LAYER " . $el->{ID} . " should a list of BIOTYPES")
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
          $self->throw("biotype $tp occurs more than once");
        }
        $biotypes{$tp} = 1;
      }
      if (exists $el->{FILTER_AGAINST}) {
        $self->throw("In layer $layer_id FILTER_AGAINST must contain a list of layer ids")
            if ref($el->{FILTER_AGAINST}) ne "ARRAY";
        @filter_against = @{$el->{FILTER_AGAINST}};

        foreach my $id (@filter_against) {
          $self->throw("In FILTER_AGAINST in layer $layer_id, '$id' is not the name ".
                "of a higher level layer")
              if not exists $layer_ids{$id};
        }
      }
      push @layers, Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLayerAnnotation::Layer->new(-id => $layer_id,
                                                                                              -discard => $discard,
                                                                                              -biotypes =>  \@biotypes,
                                                                                              -filter_object  => $filter,
                                                                                              -filter_against => \@filter_against);
      $layer_ids{$layer_id} = 1;
    }
    $self->param('_runnabledb_layers', \@layers);
  }

  return $self->param('_runnabledb_layers');
}



################################################


sub LAYERS {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('LAYERS',$val);
  }

  if ($self->param_is_defined('LAYERS')) {
    return $self->param('LAYERS');
  }
  else {
    return;
  }
}


sub FILTER {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('FILTER',$val);
  }

  return $self->param('FILTER');
}


sub SOURCEDB_REFS {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('SOURCEDB_REFS',$val);
  }

  return $self->param('SOURCEDB_REFS');
}


sub TARGETDB_REF {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('TARGETDB_REF',$val);
  }

  return $self->param('TARGETDB_REF');
}

##############################################

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLayerAnnotation::Layer;

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
