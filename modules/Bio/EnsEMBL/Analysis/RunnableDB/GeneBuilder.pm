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

Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::GeneBuilder 
  qw(GENEBUILDER_CONFIG_BY_LOGIC);
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(Gene_info attach_Analysis_to_Gene_no_support empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 
  qw(are_strands_consistent are_phases_consistent 
     is_not_folded all_exons_are_valid intron_lengths_all_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(exon_length_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils 
  qw(validate_Translation_coords contains_internal_stops print_Translation print_peptide);
use Bio::EnsEMBL::Analysis::Tools::Logger;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Function  : instatiates a GeneBuilder object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Exceptions: 
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->read_and_check_config($GENEBUILDER_CONFIG_BY_LOGIC);
  return $self;
}



sub fetch_input{
  my ($self) = @_;
  #fetch sequence
  $self->query($self->fetch_sequence); 
  #fetch genes
  $self->get_Genes;
  #print "Have ".@{$self->input_genes}." genes to cluster\n";
  #filter genes
  my @filtered_genes = @{$self->filter_genes($self->input_genes)};
  #print "Have ".@filtered_genes." filtered genes\n";
  #create genebuilder runnable

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder
    ->new(
          -query => $self->query,
          -analysis => $self->analysis,
          -genes => \@filtered_genes,
          -output_biotype => $self->OUTPUT_BIOTYPE,
          -max_transcripts_per_cluster => $self->MAX_TRANSCRIPTS_PER_CLUSTER,
          -min_short_intron_len => $self->MIN_SHORT_INTRON_LEN,
          -max_short_intron_len => $self->MAX_SHORT_INTRON_LEN,
          -blessed_biotypes => $self->BLESSED_BIOTYPES,
          -coding_only => $self->CODING_ONLY,
         );  

  $self->runnable($runnable);
  
};


sub write_output{
  my ($self) = @_;
  my $ga = $self->get_adaptor;
  my $sucessful_count = 0;
  logger_info("WRITE OUTPUT have ".@{$self->output}." genes to write");
  foreach my $gene(@{$self->output}){
    my $attach = 0;
    if(!$gene->analysis){
      my $attach = 1;
      attach_Analysis_to_Gene_no_support($gene, $self->analysis);
    }
    if($attach == 0){
    TRANSCRIPT:foreach my $transcript(@{$gene->get_all_Transcripts}){
        if(!$transcript->analysis){
          attach_Analysis_to_Gene_no_support($gene, $self->analysis);
          last TRANSCRIPT;
        }
      }
    }

   
    foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) 
    {
      $transcript->load ;
      $transcript->dbID(0);
    }

    empty_Gene($gene);
    eval{
      $ga->store($gene);
    };
    if($@){
      warning("Failed to write gene ".id($gene)." ".coord_string($gene)." $@");
    }else{
      $sucessful_count++;
      logger_info("STORED GENE ".$gene->dbID);
    }
  }
  if($sucessful_count != @{$self->output}){
    throw("Failed to write some genes");
  }
}

sub output_db{
  my ($self, $db) = @_;
  if($db){
    $self->{output_db} = $db;
  }
  if(!$self->{output_db}){
    my $db = $self->get_dbadaptor($self->OUTPUT_DB);
    $self->{output_db} = $db;
  }
  return $self->{output_db};
}

sub get_adaptor{
  my ($self) = @_;
  return $self->output_db->get_GeneAdaptor;
}

sub get_Genes{
  my ($self) = @_;
  my @genes;
  foreach my $db_name(keys(%{$self->INPUT_GENES})){
    my $gene_db = $self->get_dbadaptor($db_name);
    my $slice = $self->fetch_sequence($self->input_id, $gene_db);
    my $biotypes = $self->INPUT_GENES->{$db_name};
    foreach my $biotype(@$biotypes){
      my $genes = $slice->get_all_Genes_by_type($biotype);
      print "Retrieved ".@$genes." of type ".$biotype."\n";
      push(@genes, @$genes);
    }
  }
  $self->input_genes(\@genes);
}

sub input_genes{
  my ($self, $arg) = @_;
  if($arg){
    throw("Need to pass input genes an arrayref not ".$arg) 
      if(!ref($arg) || ref($arg) ne 'ARRAY');
    push(@{$self->{input_genes}}, @$arg);
  }
  return $self->{input_genes};
}

sub filter_genes{
  my ($self, $genes) = @_;
  $genes = $self->input_genes if(!$genes);
  print "Have ".@$genes." to filter\n";
  my @filtered;
 GENE:foreach my $gene(@$genes){
    #throw("Genebuilder only works with one gene one transcript structures")
    #  if(@{$gene->get_all_Transcripts} >= 2);
    foreach my $transcript(@{$gene->get_all_Transcripts}){
      if($self->validate_Transcript($transcript)){
        push(@filtered, $gene);
        next GENE;
      }else{
        print Gene_info($gene)." is invalid skipping\n";
        next GENE;
      }
    }
  }
  return \@filtered;
}

sub validate_Transcript{
  my ($self, $transcript) = @_;
  my $slice = $self->query;
  $slice = $transcript->slice if(!$slice);
  my $is_valid = 0;
  #basic transcript validation
  unless(are_strands_consistent($transcript)){
    print "Transcript has inconsistent strands. ";   
    $is_valid++;
  }
  unless(are_phases_consistent($transcript)){
    print "Transcript has inconsistent exon phases. ";
    $is_valid++;
  }
  unless(is_not_folded($transcript)){
    print "Transcript seems to be folded (with secondary structure). ";   
    $is_valid++;
  }
 EXON:foreach my $exon(@{$transcript->get_all_Exons}){
    if(exon_length_less_than_maximum($exon, $self->MAX_EXON_LENGTH)){
      next EXON;
    }else{
      print "Exon in transcript exceeds max length. ";  
      $is_valid++;
      last EXON;
    }
  }
  if(contains_internal_stops($transcript)){
    print "Transcript contains internal stops. ";  
    $is_valid++;
  }
  unless(validate_Translation_coords($transcript)){
    print "Transcript contains invalid translation coords. ";   
    $is_valid++;
  }
  return 0 if($is_valid >= 1);
  return 1;
}



#CONFIG METHODS


=head2 read_and_check_config

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : hashref from config file
  Function  : call the superclass method to set all the varibles and carry
  out some sanity checking
  Returntype: N/A
  Exceptions: throws if certain variables arent set properlu
  Example   : 

=cut

sub read_and_check_config{
  my ($self, $hash) = @_;
  $self->SUPER::read_and_check_config($hash);

  #######
  #CHECKS
  #######
  foreach my $var(qw(INPUT_GENES
                     OUTPUT_DB)){
    throw("RunnableDB::GeneBuilder $var config variable is not defined") 
      unless($self->$var);
  }
  my @keys = keys(%{$self->INPUT_GENES});
  throw("RunnableDB::GeneBuilder INPUT_GENES has needs to contain values") if(!@keys);
  my %unique;
  foreach my $key(@keys){
    my $biotypes = $self->INPUT_GENES->{$key};
    foreach my $biotype(@$biotypes){
      if(!$unique{$biotype}){
        $unique{$biotype} = $key;
      }else{
        if($self->BLESSED_BIOTYPES->{$biotype}){
          throw($biotype." is defined for both ".$key." and ".$unique{$biotype}.
                " and is found in the blessed biotype hash\n".
                "This is likely to cause problems for the filtering done in ".
                "the genebuilder code");
        }else{
          warning($biotype." appears twice in your listing, make sure this ".
                  "isn't for the same database otherwise it will cause issue");
        }
      }
    }
  }
  $self->OUTPUT_BIOTYPE($self->analysis->logic_name) if(!$self->OUTPUT_BIOTYPE);
}


=head2 INPUT_GENES

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GeneBuilder
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions: 
  Example   : 

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::GeneBuilder

sub INPUT_GENES{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'INPUT_GENES'} = $arg;
  }
  return $self->{'INPUT_GENES'};
}

sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'OUTPUT_DB'} = $arg;
  }
  return $self->{'OUTPUT_DB'};
}

sub OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'OUTPUT_BIOTYPE'} = $arg;
  }
  return $self->{'OUTPUT_BIOTYPE'};
}

sub MAX_TRANSCRIPTS_PER_CLUSTER{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MAX_TRANSCRIPTS_PER_CLUSTER'} = $arg;
  }
  return $self->{'MAX_TRANSCRIPTS_PER_CLUSTER'};
}

sub MIN_SHORT_INTRON_LEN{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MIN_SHORT_INTRON_LEN'} = $arg;
  }
  return $self->{'MIN_SHORT_INTRON_LEN'};
}

sub MAX_SHORT_INTRON_LEN{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MAX_SHORT_INTRON_LEN'} = $arg;
  }
  return $self->{'MAX_SHORT_INTRON_LEN'};
} 
sub BLESSED_BIOTYPES{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'BLESSED_BIOTYPES'} = $arg;
  }
  return $self->{'BLESSED_BIOTYPES'};
}



sub MAX_EXON_LENGTH{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'MAX_EXON_LENGTH'} = $arg;
  }
  return $self->{'MAX_EXON_LENGTH'};
}

sub CODING_ONLY{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'CODING_ONLY'} = $arg;
  }
  return $self->{'CODING_ONLY'};
}



1;
