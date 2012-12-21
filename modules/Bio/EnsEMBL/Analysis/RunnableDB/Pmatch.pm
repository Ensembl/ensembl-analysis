=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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



# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/Pmatch.pm,v $
# $Version: $
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
