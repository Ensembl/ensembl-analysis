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

Bio::EnsEMBL::Analysis::RunnableDB::Solexa2GenesLiteNewCloud - 

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $refine_genes = Bio::EnsEMBL::Analysis::RunnableDB::Solexa2Genes->new (
                                                    -db      => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $refine_genes->fetch_input();
  $refine_genes->run();
  $refine_genes->write_output(); #writes to DB


=head1 DESCRIPTION


The databases containing the various features to combine is defined in
Bio::EnsEMBL::Analysis::Config::Databases and the configuration for the 
module is defined in Bio::EnsEMBL::Analysis::Config::GeneBuild::Solexa2Genes

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Solexa2GenesLiteNewCloud;

use warnings ;
use strict;

#use Bio::EnsEMBL::Analysis::RunnableDB;
#use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild; 

use Bio::EnsEMBL::Analysis::RunnableDB::Solexa2GenesLiteNew;
use Bio::EnsEMBL::Analysis::Config::ExonerateSolexaCloudConfig;

#use Bio::EnsEMBL::Analysis::Config::GeneBuild::Solexa2GenesLiteNew;
#use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils ;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Time::HiRes qw(gettimeofday);  
use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::Solexa2GenesLiteNew);


# this module is reading from the united alignment db and wrting rough models to the output db 
# the original config is Solexa2GenesLiteNew 
# we need to set : OUTPUT_DB 
#                  ALIGNMENT_DB 
#                  USE_ANALYSIS_LOGIC_NAME_AS_DEFAULT_GENE_OUTPUT_BIOTYPE  

sub new {
  my ( $class, @args ) = @_;

  push @args, ("-ignore_config_file" , 1 ); 
  my $self = $class->SUPER::new(@args);

  # we have to set a few config vars automatically (Solexa2GenesLiteNew) 
  $self->ANALYSIS_BASE_BATCH_CONFIG($ANALYSIS_BASE_BATCH_CONFIG);   

  # this module retrieves an input_id in the format BATCH_NR@chromosome:GRCh37:1:1223:348920:1  
  my ($base_batch, $rest ) = split "@",$self->input_id;
  $self->input_id($rest);
  $self->base_batch($base_batch); 

  $self->OUTPUT_DB($self->rough_model_output_db);  
  $self->ALIGNMENT_DB($self->alignment_db);  
  $self->USE_ANALYSIS_LOGIC_NAME_AS_DEFAULT_GENE_OUTPUT_BIOTYPE(1); 

  print "Aligment data is fetched from " . $self->ALIGNMENT_DB . " (see Bio::EnsEMBL::Analysis::Config::ExonerateSolexaCloudConfig)\n"; 
  print "Rough models will be written to : ". $self->OUTPUT_DB . " (see Bio::EnsEMBL::Analysis::Config::ExonerateSolexaCloudConfig)\n"; 
  print "Biotype for rough models : " . $self->analysis->logic_name ."\n";  

  return $self;
}


sub alignment_db{  
  my ($self) = shift; 
  return ${$self->ANALYSIS_BASE_BATCH_CONFIG}{$self->base_batch}{"STAGE_1_UNITED_OUTPUT_DB"};
}

sub rough_model_output_db {  
  my ($self) = shift; 
  return ${$self->ANALYSIS_BASE_BATCH_CONFIG}{$self->base_batch}{"STAGE_2_OUTPUT_DBNAME"};
}

use vars '$AUTOLOAD';
sub AUTOLOAD {  
 my ($self,$val) = @_;
 (my $routine_name=$AUTOLOAD)=~s/.*:://; #trim package name
 $self->{$routine_name}=$val if $val ; 
 return $self->{$routine_name} ; 
}
sub DESTROY {} # required due to AUTOLOAD



1; 
