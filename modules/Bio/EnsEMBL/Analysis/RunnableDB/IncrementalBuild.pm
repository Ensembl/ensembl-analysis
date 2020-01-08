=head1 LICENSE
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::IncrementalBuild - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::IncrementalBuild;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::IncrementalBuild;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($INCREMENTAL_HASH);

  return $self;
}



sub fetch_input{
  my ($self) = @_;
  my $secondary_slice = $self->fetch_sequence($self->input_id, 
                                              $self->secondary_database);
  my $secondary_genes = $self->get_Genes($secondary_slice,
                                         $self->SECONDARY_BIOTYPES);
  my $primary_slice = $self->fetch_sequence($self->input_id,
                                            $self->primary_database);
  my $primary_genes = $self->get_Genes($primary_slice,
                                       $self->PRIMARY_BIOTYPES);
  $self->query($secondary_slice);
  $self->secondary_genes($secondary_genes);
  $self->primary_genes($primary_genes);
  $self->primary_slice($primary_slice);
}

sub run{
  my ($self) = @_;

  my $mask_regions = $self->mask_gene_regions($self->primary_genes);
  my @non_overlaps;
  my @overlaps;
  my @genes;
  if($self->PRE_FILTER){
    #print "Running a prefilter\n";
    my $filtered_genes = $self->filter_object->filter_genes($self->secondary_genes);
    push(@genes, @$filtered_genes);
  }else{
    push(@genes, @{$self->secondary_genes});
  }
  #print "Have ".@genes." secondary genes\n";
  foreach my $gene(@genes){
    my $keep_gene = 1;
    my $mask_reg_idx = 0;
    my @test_regions = ({start => $gene->start, end => $gene->end});
    FEAT: foreach my $f (@test_regions) {
        #print "Comparing ".$f->{'start'}." ".$f->{'end'}."\n";
        #print "have ".@$mask_regions." mask regions\n";
        for( ; $mask_reg_idx < @$mask_regions; ) {
        my $mask_reg = $mask_regions->[$mask_reg_idx];
        #print "To ".$mask_reg->{'start'}." ".$mask_reg->{'end'}."\n";
        if ($mask_reg->{'start'} > $f->{'end'}) {
          # no mask regions will overlap this feature
          next FEAT;
		    }
        elsif ( $mask_reg->{'end'} >= $f->{'start'}) {
          # overlap			
          #print "Overlap found\n";
          $keep_gene = 0;
          last FEAT;
        }			
        else {
          $mask_reg_idx++;
        }
      }
    }
    if ($keep_gene) {
      my $gene_to_keep;
      if($self->CALCULATE_TRANSLATION){
        $gene_to_keep = $self->calculate_translation($gene);
      }else{
        $gene_to_keep = $gene;
      }
      push @non_overlaps, $gene_to_keep if($gene_to_keep);
    }else{
      push(@overlaps, $gene);
    }
  }
  #print "Have ".@{$self->primary_genes}." primary genes\n";
  #print "Have ".@non_overlaps." non overlapping genes\n";
  #print "Have ".@overlaps." overlapping genes\n";
  my @non_overlapping_genes;
  if($self->POST_FILTER){
    my $filtered_genes = $self->filter_object->filter_genes(\@non_overlaps);
    push(@non_overlapping_genes, @$filtered_genes);
  }else{
    push(@non_overlapping_genes, @non_overlaps);
  }
  my @output;
  push(@output, @non_overlapping_genes);
  if($self->NEW_BIOTYPE){
    foreach my $output(@output){
      $output->biotype($self->NEW_BIOTYPE);
    }
  }
  my $primary_genes = $self->primary_genes;
  print "Adding primary genes to set\n" if($self->STORE_PRIMARY);
  push(@output, @$primary_genes) if($self->STORE_PRIMARY);
  $self->output(\@output);
  #print "have ".@output." genes to store\n";
  $self->overlapping_genes(\@overlaps);
}

sub write_output{
  my ($self) = @_;

  my $db = $self->output_database;
  my $gene_adaptor = $db->get_GeneAdaptor;
  foreach my $gene (@{$self->output}){
    empty_Gene($gene);
    eval{
      $gene_adaptor->store($gene);
    };
    if($@){
      throw("Failed to store gene $@");
    }
  }
}

sub overlapping_genes{
  my ($self, $value) = @_;
  if($value){
    $self->{'overlapping_genes'} = $value;
  }
  return $self->{'overlapping_genes'};
}

sub secondary_genes{
  my ($self, $value) = @_;
  if($value){
    $self->{'secondary_genes'} = $value;
  }
  return $self->{'secondary_genes'};
}

sub primary_genes{
  my ($self, $value) = @_;
  if($value){
    $self->{'primary_genes'} = $value;
  }
  return $self->{'primary_genes'};
}

sub primary_slice{
  my ($self, $value) = @_;
  if($value){
    $self->{'primary_slice'} = $value;
  }
  return $self->{'primary_slice'};
}


sub secondary_database{
  my ($self, $value) = @_;
  if($value){
    $self->{'secondary_database'} = $value;
  }
  if(!$self->{'secondary_database'}){
    my $db = $self->get_dbadaptor($self->SECONDARY_DB);
    $db->dnadb($self->db);
    $self->{'secondary_database'} = $db;
  }
  return $self->{'secondary_database'};
}

sub primary_database{
  my ($self, $value) = @_;
  if($value){
    $self->{'primary_database'} = $value;
  }
  if(!$self->{'primary_database'}){
    my $db = $self->get_dbadaptor($self->PRIMARY_DB);
    $db->dnadb($self->db);
    $self->{'primary_database'} = $db;
  }
  return $self->{'primary_database'};
}

sub output_database{
  my ($self, $value) = @_;
  if($value){
    $self->{'output_database'} = $value;
  }
  if(!$self->{'output_database'}){
    my $db = $self->get_dbadaptor($self->OUTPUT_DB);
    $db->dnadb($self->db);
    $self->{'output_database'} = $db;
  }
  return $self->{'output_database'};
}


sub PRIMARY_BIOTYPES{
  my ($self,$value) = @_;
  
  if ($value) {
    $self->{'PRIMARY_BIOTYPES'} = $value;
  }

  if (exists($self->{'PRIMARY_BIOTYPES'})) {
    return $self->{'PRIMARY_BIOTYPES'};
  } else {
    return undef;
  }
}


sub PRIMARY_DB{
  my ($self,$value) = @_;
  
  if ($value) {
    $self->{'PRIMARY_DB'} = $value;
  }

  if (exists($self->{'PRIMARY_DB'})) {
    return $self->{'PRIMARY_DB'};
  } else {
    return undef;
  }
}


sub SECONDARY_BIOTYPES{
  my ($self,$value) = @_;
  
  if ($value) {
    $self->{'SECONDARY_BIOTYPES'} = $value;
  }

  if (exists($self->{'SECONDARY_BIOTYPES'})) {
    return $self->{'SECONDARY_BIOTYPES'};
  } else {
    return undef;
  }
}


sub SECONDARY_DB{
  my ($self,$value) = @_;
  
  if ($value) {
    $self->{'SECONDARY_DB'} = $value;
  }

  if (exists($self->{'SECONDARY_DB'})) {
    return $self->{'SECONDARY_DB'};
  } else {
    return undef;
  }
}

sub SECONDARY_FILTER{
  my ($self,$value) = @_;
  
  if ($value) {
    $self->{'SECONDARY_FILTER'} = $value;
  }

  if (exists($self->{'SECONDARY_FILTER'})) {
    return $self->{'SECONDARY_FILTER'};
  } else {
    return undef;
  }
}

sub SECONDARY_PADDING{
  my ($self,$value) = @_;
  
  if ($value) {
    $self->{'SECONDARY_PADDING'} = $value;
  }

  if (exists($self->{'SECONDARY_PADDING'})) {
    return $self->{'SECONDARY_PADDING'};
  } else {
    return undef;
  }
}


sub OUTPUT_DB{
  my ($self,$value) = @_;
  
  if ($value) {
    $self->{'OUTPUT_DB'} = $value;
  }

  if (exists($self->{'OUTPUT_DB'})) {
    return $self->{'OUTPUT_DB'};
  } else {
    return undef;
  }
}

sub PRE_FILTER{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'PRE_FILTER'} = $value;
  }

  if (exists($self->{'PRE_FILTER'})) {
    return $self->{'PRE_FILTER'};
  } else {
    return undef;
  }
}

sub POST_FILTER{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'POST_FILTER'} = $value;
  }

  if (exists($self->{'POST_FILTER'})) {
    return $self->{'POST_FILTER'};
  } else {
    return undef;
  }
}

sub STORE_PRIMARY{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'STORE_PRIMARY'} = $value;
  }

  if (exists($self->{'STORE_PRIMARY'})) {
    return $self->{'STORE_PRIMARY'};
  } else {
    return undef;
  }
}


sub CALCULATE_TRANSLATION{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'CALCULATE_TRANSLATION'} = $value;
  }

  if (exists($self->{'CALCULATE_TRANSLATION'})) {
    return $self->{'CALCULATE_TRANSLATION'};
  } else {
    return undef;
  }
}

sub DISCARD_IF_NO_ORF{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'DISCARD_IF_NO_ORF'} = $value;
  }

  if (exists($self->{'DISCARD_IF_NO_ORF'})) {
    return $self->{'DISCARD_IF_NO_ORF'};
  } else {
    return undef;
  }
}


sub NEW_BIOTYPE{
  my ($self,$value) = @_;
  
  if ($value) {
    $self->{'NEW_BIOTYPE'} = $value;
  }

  if (exists($self->{'NEW_BIOTYPE'})) {
    return $self->{'NEW_BIOTYPE'};
  } else {
    return undef;
  }
}

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($INCREMENTAL_HASH);
  
  #######
  #CHECKS
  #######
  foreach my $config_var (qw(PRIMARY_DB
                             SECONDARY_DB
                             OUTPUT_DB)) {

    throw("You must define $config_var in config for logic '".
          $self->analysis->logicname."'")
        if not defined $self->$config_var;
  }
  

  if($self->SECONDARY_FILTER){
    if($self->SECONDARY_FILTER->{OBJECT}){
      my $module = $self->SECONDARY_FILTER->{OBJECT};
      my $pars = $self->SECONDARY_FILTER->{PARAMETERS};
      (my $class = $module) =~ s/::/\//g;
      eval{
        require "$class.pm";
      };
      throw("Couldn't require ".$class.
            " IncrementalBuild:read_and_check_config $@") if($@);
      
      $self->filter_object($module->new(%{$pars}));
    }
  }
}


sub filter_object{
  my ($self,$value) = @_;
  
  if ( $value) {
    $self->{'filter_object'} = $value;
  }

  if (exists($self->{'filter_object'})) {
    return $self->{'filter_object'};
  } else {
    return undef;
  }
}



sub mask_gene_regions{
  my ($self, $genes) = @_;
  #print "Generating mask regions for ".@$genes." genes\n";
  my @mask_gene_regions;
  foreach my $gene(@$genes){
    push @mask_gene_regions, { start => $gene->start, 
                               end   => $gene->end };
  }
  @mask_gene_regions = sort {$a->{'start'} <=> $b->{'start'} }
    @mask_gene_regions; 
  #print "Returning ".@mask_gene_regions." mask regions\n";
  return \@mask_gene_regions;
}



sub calculate_translation{
  my ($self, $gene) = @_;

  my $gene_to_keep = Bio::EnsEMBL::Gene->new();
  $gene_to_keep->biotype($gene->biotype);
  $gene_to_keep->confidence($gene->confidence);
  $gene_to_keep->analysis($self->analysis);
  $gene_to_keep->source($gene->source);
  $gene_to_keep->stable_id($gene->stable_id);
  
  my @new_transcripts;
 TRANS:foreach my $transcript (@{$gene->get_all_Transcripts}){
    my $seq = $transcript->slice->seq;
    #print "The slice seq of the transcript = ".$seq."\n";
    #print "dbadaptor of transcript ".
    #  $transcript->adaptor->db->dnadb->dbname."\n";
    my $new_translation = Bio::EnsEMBL::Pipeline::Tools::TranslationUtils->return_translation($transcript);
    if(!$new_translation && $self->DISCARD_IF_NO_ORF){
      next TRANS;
    }else{
      $transcript->translation($new_translation);
      push(@new_transcripts, $transcript);
    }
  }
  if(@new_transcripts >= 1){
    foreach my $new_transcript(@new_transcripts){
      $gene_to_keep->add_Transcript($new_transcript);
    }
    return $gene_to_keep;
  }else{
    return undef;
  }
  
}



sub get_Genes{
  my ($self, $slice, $biotypes) = @_;
  my @genes;
  if(@$biotypes){
    foreach my $biotype(@$biotypes){
      push(@genes, @{$slice->get_all_Genes_by_type($biotype)});
    }
  }else{
    push(@genes, @{$slice->get_all_Genes});
  }
  return \@genes;
}


1;
