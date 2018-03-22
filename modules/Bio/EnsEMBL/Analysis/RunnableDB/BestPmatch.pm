=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::BestPmatch - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::BestPmatch;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Pmatch qw(BESTPMATCH_CONFIG_BY_LOGIC);
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(empty_Object);
use Bio::EnsEMBL::Analysis::Runnable::BestPmatch;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  print "\nReading config : Bio/EnsEMBL/Analysis/Config/GeneBuild/Pmatch.pm\n\n" ; 
  $self->read_and_check_config($BESTPMATCH_CONFIG_BY_LOGIC);
  return $self;
}


sub fetch_input{
  my ($self) = @_;
  my @features ;
  my @a_dbs;
  if (ref($self->INPUT_DB) =~ /ARRAY/) {
      print STDERR 'in ref',"\n";
      @a_dbs = @{$self->INPUT_DB};
  }
  else {
      print STDERR 'else ref',"\n";
      @a_dbs = ($self->INPUT_DB);
  }
  foreach my $dba (@a_dbs) {
      my $pafa = $self->get_dbadaptor($dba)->get_ProteinAlignFeatureAdaptor;  
   
      if ( ref($self->PMATCH_LOGIC_NAME)=~m/ARRAY/ ) {   
         for my  $logic_name  (@{ $self->PMATCH_LOGIC_NAME}) {   
            my @f = @{$pafa->fetch_all_by_logic_name($logic_name)} ; 
            print "Have fetched ".@f." with logic_name : $logic_name from ".$pafa->dbc->dbname."\n"; 
            push @features, @f ;  
         }  
      } else {  
        @features = @{$pafa->fetch_all_by_logic_name($self->PMATCH_LOGIC_NAME)} ; 
        print "Have fetched ".@features." with logic_name : ".$self->PMATCH_LOGIC_NAME." from ".$pafa->dbc->dbname."\n";
      }  
  }  
  $self->pmatch_features(\@features);
}

sub run{
  my ($self) = @_;
  my $pmf2 = Bio::EnsEMBL::Analysis::Runnable::BestPmatch->
    new( 
        '-protein_hits' => $self->pmatch_features,
        '-min_coverage' => $self->MIN_COVERAGE,
        '-analysis' => $self->analysis,
       );
  $pmf2->run;
  my @output = @{$pmf2->output};
  foreach my $output(@output){
    empty_Object($output);
  }
  #my @unique = $self->uniquify(@output);
  $self->output(\@output);
}


sub get_adaptor{
  my ($self) = @_;
  my $output_db = $self->get_dbadaptor($self->OUTPUT_DB);
  return $output_db->get_ProteinAlignFeatureAdaptor;
}

sub pmatch_features{
  my ($self, $features) = @_;
  if($features){
    if(ref($features) eq "ARRAY"){
      push(@{$self->{'pmatch_features'}}, @$features);
    }else{
      push(@{$self->{'pmatch_features'}}, $features);
    }
  }
  return $self->{'pmatch_features'};
}

sub PMATCH_LOGIC_NAME{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PMATCH_LOGIC_NAME'} = $arg;
  }
  return $self->{'PMATCH_LOGIC_NAME'};
}

sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'OUTPUT_DB'} = $arg;
  }
  return $self->{'OUTPUT_DB'};
}



sub INPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'INPUT_DB'} = $arg;
  }
  return $self->{'INPUT_DB'};
}


sub MIN_COVERAGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'MIN_COVERAGE'} = $arg;
  }
  return $self->{'MIN_COVERAGE'} = $arg;
}





sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($BESTPMATCH_CONFIG_BY_LOGIC);
  
  #######
  #CHECKS
  #######
  foreach my $config_var (qw(PMATCH_LOGIC_NAME
                             OUTPUT_DB
                             INPUT_DB)){
    throw("You must define $config_var in config for logic '".
          $self->analysis->logic_name."'")
      if not defined $self->$config_var;
  }
};


