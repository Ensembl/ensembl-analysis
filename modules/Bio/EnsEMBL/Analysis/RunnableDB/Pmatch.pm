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

Bio::EnsEMBL::Analysis::RunnableDB::Pmatch - 

=head1 SYNOPSIS


=head1 DESCRIPTION


Pmatch is a fast alignment program written by Richard Durbin we used
to align species specific proteins to the genome, (We also use it to
align very closely related species proteins sets e.g Mouse to Rat or
Fugu to Tetraodon).

The pmatch source code is available from the
sanger cvs respository in module rd-utils,
(http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/rd-utils/?root=ensembl).
The code to run this process in the ensembl code
base can be found in 2 RunnableDBs and a config
file. Bio::EnsEMBL::Analysis::RunnableDB::Pmatch,
Bio::EnsEMBL::Analysis::RunnableDB::BestPmatch and
Bio::EnsEMBL::Analysis::Config::GeneBuild::Pmatch


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut



package Bio::EnsEMBL::Analysis::RunnableDB::Pmatch;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Pmatch qw(PMATCH_CONFIG_BY_LOGIC);
use Bio::EnsEMBL::Analysis::Runnable::Pmatch;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($PMATCH_CONFIG_BY_LOGIC);

  return $self;
}



sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence($self->input_id, $self->db, 
                                    $self->REPEAT_MASKING);
  $self->query($slice);
  my %parameters = %{$self->parameters_hash};
  my $program = $self->analysis->program_file;
  $program = $self->BINARY_LOCATION if(!$program);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Pmatch
    ->new(
          -query => $self->query,
          -program => $program,
          -analysis => $self->analysis,
          -protein_file => $self->PROTEIN_FILE,
          -max_intron_length => $self->MAX_INTRON_LENGTH,
          -min_coverage => $self->MIN_COVERAGE,
          -options => $self->OPTIONS,
         );
  $self->runnable($runnable);
}


sub PROTEIN_FILE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PROTEIN_FILE'} = $arg;
  }
  return $self->{'PROTEIN_FILE'};
}


sub MIN_COVERAGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'MIN_COVERAGE'} = $arg;
  }
  return $self->{'MIN_COVERAGE'};
}

sub BINARY_LOCATION{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'BINARY_LOCATION'} = $arg;
  }
  return $self->{'BINARY_LOCATION'};
}


sub MAX_INTRON_LENGTH{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'MAX_INTRON_SIZE'} = $arg;
  }
  return $self->{'MAX_INTRON_SIZE'};
}


sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'OUTPUT_DB'} = $arg;
  }
  return $self->{'OUTPUT_DB'};
}

sub REPEAT_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    throw("Runnable::Pmatch ".$arg." must be an array ref of logic names not .".$arg)
      unless(ref($arg) eq 'ARRAY');
    $self->{'REPEAT_MASKING'} = $arg;
  }
  return $self->{'REPEAT_MASKING'};
}

sub OPTIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OPTIONS'} = $value;
  }

  if (exists($self->{'_CONFIG_OPTIONS'})) {
    return $self->{'_CONFIG_OPTIONS'};
  } else {
    return undef;
  }
}



sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($PMATCH_CONFIG_BY_LOGIC);
  
  #######
  #CHECKS
  #######
  foreach my $config_var (qw(PROTEIN_FILE
                             OUTPUT_DB)){
    throw("You must define $config_var in config for logic '".
          $self->analysis->logic_name."'")
      if not defined $self->$config_var;
  }
};


sub get_adaptor{
  my ($self) = @_;
  my $output_db = $self->get_dbadaptor($self->OUTPUT_DB);
  return $output_db->get_ProteinAlignFeatureAdaptor;
}


sub write_output{
  my ($self) = @_;
  my $adaptor = $self->get_adaptor;
  my %unique;
 FEATURE:foreach my $feature(@{$self->output}){
    $feature->analysis($self->analysis);
    $feature->slice($self->query) if(!$feature->slice);
    my $unique_string = $feature->start." ".$feature->end." ".$feature->score." ".$feature->hseqname;
    next FEATURE if($unique{$unique_string});
    $unique{$unique_string} = 1;
    $self->feature_factory->validate($feature);
    eval{
      $adaptor->store($feature);
    };
    if($@){
      throw("RunnableDB:store failed, failed to write ".$feature." to ".
            "the database ".$adaptor->dbc->dbname." $@");
    }
  }
  return 1;
}


1;
