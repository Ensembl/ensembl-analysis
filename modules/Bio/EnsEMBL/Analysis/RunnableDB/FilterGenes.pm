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

Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes - 

=head1 SYNOPSIS

my $filtergenes = Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes->new(
                              -db         => $refdb,
			      -analysis   => $analysis_obj,
			      -input_id => $slice_name
			     );

$filtergenes->fetch_input();
$filtergenes->run();
$filtergenes->output();
$filtergenes->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps a genefilter object as defined in its configuration file
The filter object must implement a method called filter_genes which expects
an arrayref of genes as its arguments and returns 2 arrayref of genes, the
first the accepted gene objects, the 2nd the rejected gene objects

These genes are then relabelled as described in the Config file and stored in
the output database

=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils 
  qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Config::GeneBuild::FilterGenes 
  qw(FILTER_CONFIG_BY_LOGIC);

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Function  : This simply checks the configuration is correct
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Exceptions: 
  Example   : 

=cut



sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($FILTER_CONFIG_BY_LOGIC);

  return $self;
}



=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Function  : fetch genes and sequence for the analysis
  Returntype: n/a
  Exceptions: 
  Example   : 

=cut


sub fetch_input{
  my ($self) = @_;

  $self->fetch_genes;
  my $slice = $self->fetch_sequence($self->input_id, $self->db);
  $self->query($slice);
}



=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Function  : run the defined filter object storing the output
  in the object
  Returntype: n/a
  Exceptions: 
  Example   : 

=cut



sub run{
  my ($self) = @_;

  my $filter_object = $self->filter_object;
  my $genes = $self->genes;
  my ($kept, $removed) = $filter_object->filter_genes($genes);
  $self->output($kept);
  $self->removed($removed);
}



=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Function  : store the output and filtered out genes in the
  database, relabelling the biotype as appropriate
  Returntype: n/a
  Exceptions: warns if fails to store a gene, throws if all stores
  fail
  Example   : 

=cut



sub write_output{
  my ($self) = @_;
  my $gene_adaptor = $self->get_dbadaptor($self->OUTPUT_DB)->get_GeneAdaptor;
  my $store_count = 0;
  my $total = @{$self->output};
  logger_info("Storing ".@{$self->output}." genes in ".$self->OUTPUT_DB);
  foreach my $output(@{$self->output}){
    $output->biotype($self->RELABEL_KEPT) if($self->RELABEL_KEPT);
    empty_Gene($output, $self->REMOVE_STABLE_IDS, $self->REMOVE_XREFS);
    eval{
      $gene_adaptor->store($output);
    };
    if($@){
      warning("Failed to store gene ".$@." FilterGenes:write_output");
    }else{
      $store_count++;
    }
  }
  if($self->RELABEL_REMOVED){
    $total += @{$self->removed};
    foreach my $removed(@{$self->removed}){
      $removed->biotype($self->RELABEL_REMOVED);
      empty_Gene($removed, 1, 1);
      eval{
        $gene_adaptor->store($removed);
      };
      if($@){
        warning("Failed to store relabled removed gene ".$@.
                " FilterGenes:write_output");
      }else{
        $store_count++;
      }
    }
  }
  if($store_count == 0){
    throw("Failed to store any genes in FilterGenes:write_output");
  }
}



=head2 fetch_genes

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Function  : gets the genes on the basis of the GENE_SET hash in the
  configuration
  Returntype: arrayref of Bio::EnsEMBL::Gene
  Exceptions: n/a
  Example   : 

=cut


sub fetch_genes{
  my ($self) = @_;

  my $geneset_hash = $self->GENE_SET;
  my %slice_hash;
  my @to_filter;
  foreach my $key(keys(%$geneset_hash)){
    my $hash = $geneset_hash->{$key};
    my $dbname = $hash->{dbname};
    my $biotype_list = $hash->{biotypes};
    my %biotype_hash;
    foreach my $biotype (@$biotype_list){
      $biotype_hash{$biotype} = 1;
    }
    logger_info("Fetching ".$key." geneset from ".$dbname);
    if(!$slice_hash{$dbname}){
      my $db = $self->get_dbadaptor($dbname);
      my $slice = $self->fetch_sequence($self->input_id, $db);
      $slice_hash{$dbname} = $slice;
    }
    my $slice = $slice_hash{$dbname};
    my $genes = $slice->get_all_Genes;
  GENE:foreach my $gene(@$genes){
      next GENE if(keys(%biotype_hash) && !$biotype_hash{$gene->biotype});
      push(@to_filter,$gene);
    }
  }
  $self->genes(\@to_filter);
  return \@to_filter;
}



=head2 genes

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Arg [2]   : arrayref of Bio::EnsEMBL::Gene
  Function  : stores the arrayref
  Returntype: arrayref of Bio::EnsEMBL::Gene
  Exceptions: 
  Example   : 

=cut



sub genes{
  my ($self, $arg) = @_;
  if($arg){
    $self->{genes} = $arg;
  }
  return $self->{genes};
}


=head2 removed

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Arg [2]   : arrayref of Bio::EnsEMBL::Gene
  Function  : stores the arrayref
  Returntype: arrayref of Bio::EnsEMBL::Gene
  Exceptions: 
  Example   : 

=cut



sub removed{
  my ($self, $arg) = @_;
  if($arg){
    $self->{removed} = $arg;
  }
  return $self->{removed};
}



=head2 filter_object

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Arg [2]   : a filter object with the method filter_genes
  Function  : returns said filter object which is created on the
  basis of the config is one doesn't already exist
  Returntype: filter object
  Exceptions: throws if object can't carry out the needed method 
  Example   : 

=cut



sub filter_object{
  my ($self, $arg) = @_;
  if($arg){
    $self->{filter_object} = $arg;
  }
  if(!$self->{filter_object}){
    $self->require_module($self->FILTER_OBJECT);
    my %params = %{$self->FILTER_PARAMS};
    my $filter = $self->FILTER_OBJECT->new(
                                           %params
                                          );
    $self->{filter_object} = $filter;
  }
  throw($self->{filter_object}." must have a method called filter_genes")
    unless(!$self->{filter_object} || 
           $self->{filter_object}->can("filter_genes"));
  return $self->{filter_object};
}



=head2 read_and_check_config

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Function  : calls the superclass method then does some
  sanity checking of the config
  Returntype: n/a
  Exceptions: throws if there are problems with the config
  Example   : 

=cut


sub read_and_check_config {
  my ($self, $config) = @_;

  $self->SUPER::read_and_check_config($config);
  my $logic = $self->analysis->logic_name;
  ##########
  # CHECKS
  ##########
  foreach my $config_var (qw(GENE_SET 
                             OUTPUT_DB
                             FILTER_OBJECT)) {
    throw("You must define $config_var in config for logic '$logic'")
      if not defined $self->$config_var;
  }

  my $gene_hash = $self->GENE_SET;
  my $count = 0;
  foreach my $label(keys(%$gene_hash)){
    $count++;
    my $values_hash = $gene_hash->{$label};
    if(!$values_hash->{dbname} || !$values_hash->{biotypes}){
      throw("The GENE_SET hash must be a hash of hashes, each individual ".
            "hash cotaining 2 key dbname pointing to the database desired ".
            "from Databases.pm and biotypes pointing to an array ref of ".
            "biotypes you want fetched from said databases ".$label." hash ".
            "is missing one or both of these");
    }
  }
  if(!$count){
    throw("There GENE_SET hash must contain keys pointing to other hashes".
          " This hash ".$gene_hash." doesn't");
  }
  
}


=head2 GENE_SET

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FilterGenes
  Arg [2]   : Here are a variety of things, it is mostly strings but there are
  some hash and arrayrefs
  Function  : just to store the defined config variable
  Returntype: what ever is passed in in Arg[2]
  Exceptions: n/a
  Example   : 

=cut



sub GENE_SET{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_SET} = $arg;
  }
  return $self->{GENE_SET};
}

sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_DB} = $arg;
  }
  return $self->{OUTPUT_DB};
}

sub RELABEL_REMOVED{
  my ($self, $arg) = @_;
  if($arg){
    $self->{RELABEL_REMOVED} = $arg;
  }
  return $self->{RELABEL_REMOVED};
}

sub RELABEL_KEPT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{RELABEL_KEPT} = $arg;
  }
  return $self->{RELABEL_KEPT};
}

sub FILTER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_OBJECT} = $arg;
  }
  return $self->{FILTER_OBJECT};
}

sub FILTER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_PARAMS} = $arg;
  }
  return $self->{FILTER_PARAMS};
}

sub REMOVE_STABLE_IDS{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{REMOVE_STABLE_IDS} = $arg;
  }
  return $self->{REMOVE_STABLE_IDS};
}


sub REMOVE_XREFS{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{REMOVE_XREFS} = $arg;
  }
  return $self->{REMOVE_XREFS};
}
