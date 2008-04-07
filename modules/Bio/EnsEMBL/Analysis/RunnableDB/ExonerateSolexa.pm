
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa->new(
    -db         => $refdb,
    -analysis   => $analysis_obj,
  );

$runnableDB->fetch_input();
$runnableDB->run();
$runnableDB->write_output(); #writes to DB

=head1 DESCRIPTION

Extends Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature to allow
use of compressed dna align features, useful when aligning millions of short 
Solexa reads

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateSolexa;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::ExonerateSolexa;
use Bio::EnsEMBL::Analysis::Tools::Utilities;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::ExonerateAlignFeature);


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($EXONERATE_SOLEXA_CONFIG_BY_LOGIC);

  return $self;
}

=head2 write_output

  Arg [1]   : array ref of Bio::EnsEMBL::DnaDnaAlignFeature
  Function  : Overrides the write_output method from the superclass
              to allow use of compressed DnaAlignFeatures
  Returntype: 1
  Exceptions: Throws if the feature cannot be stored

=cut

sub write_output {
  my ( $self, @output ) = @_;
  # Flag set to 1 = return a pipeline adaptor
  my $outdb = get_db_adaptor_by_string("TEST2_DB",undef,1);
# $outdb = $self->get_dbadaptor($self->OUT_DB,1);

  my $fa;
  if ( $self->COMPRESSION ) {
    $fa = $outdb->get_CompressedDnaAlignFeatureAdaptor;
  } else {
    $fa = $outdb->get_DnaAlignFeatureAdaptor;
  }
  
  foreach my $f (@{$self->output}){
    # calculate the hcoverage and use the evalue feild to store the
    # depth of the features
    
    $f->p_value(1) if $self->COMPRESSION;

    eval{
      $fa->store($f);
    };
    if ($@) {
      $self->throw("Unable to store DnaAlignFeature\n $@");
    }
  }
}

###########################################################
# containers

sub INTRON_OVERLAP {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_INTRON_OVERLAP'} = $value;
  }
  
  if (exists($self->{'_CONFIG_INTRON_OVERLAP'})) {
    return $self->{'_CONFIG_INTRON_OVERLAP'};
  } else {
    return undef;
  }
}

sub BIOTYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BIOTYPE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BIOTYPE'})) {
    return $self->{'_CONFIG_BIOTYPE'};
  } else {
    return undef;
  }
}

sub TRANSDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_TRANSDB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_TRANSDB'})) {
    return $self->{'_CONFIG_TRANSDB'};
  } else {
    return undef;
  }
}

sub OUT_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUT_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OUT_DB'})) {
    return $self->{'_CONFIG_OUT_DB'};
  } else {
    return undef;
  }
}

sub COMPRESSION {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_COMPRESSION'} = $value;
  }
  
  if (exists($self->{'_CONFIG_COMPRESSION'})) {
    return $self->{'_CONFIG_COMPRESSION'};
  } else {
    return undef;
  }
}




1;
