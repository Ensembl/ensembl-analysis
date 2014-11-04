package Bio::EnsEMBL::Analysis::RunnableDB::HiveSubmitAnalysis;

use strict;
use warnings;
use feature 'say';

use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::RunnableDB::HiveBaseRunnable');

sub fetch_input {
  my $self = shift;
  $self->db($self->get_dba($self->param('reference_db')));
  return 1;
}

sub run {
  my $self = shift;

  if (!($self->param('slice')) && !($self->param('single')) && !($self->param('file')) &&
      !($self->param('translation_id')) && !($self->param('hap_pair'))) {
    throw("Must define input as either contig, slice, file, translation_id ".
          "single, seq_level or top_level or hap_pair");
  }

  my $input_id_factory = new Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory
  (
   -db => $self->db(),
   -slice => $self->param('slice'),
   -single => $self->param('single'),
   -file => $self->param('file'),
   -translation_id => $self->param('translation_id'),
   -seq_level => $self->param('seq_level'),
   -top_level => $self->param('top_level'),
   -include_non_reference => $self->param('include_non_reference'),
   -dir => $self->param('dir'),
   -regex => $self->param('regex'),
   -single_name => 'genome', # Don't know why this is set this way
   -logic_name => $self->param('logic_name'),
   -input_id_type => $self->param('input_id_type'),
   -coord_system => $self->param('coord_system_name'),
   -coord_system_version => $self->param('coord_system_version'),
   -slice_size => $self->param('slice_size'),
   -slice_overlaps => $self->param('slice_overlap'),
   -seq_region_name => $self->param('seq_region_name'),
   -hap_pair => $self->param('hap_pair'),
  );

  $input_id_factory->generate_input_ids;
  $self->{'input_id_factory'} = $input_id_factory;

  return 1;
}


sub write_output {
  my $self = shift;

  my $output_ids = $self->{'input_id_factory'}->input_ids();

  unless(scalar(@{$output_ids})) {
    warning("No input ids generated for this analysis!");
  }


  foreach my $id (@{$output_ids}) {

#    unless($id =~ /chromosome\:GRCh38\:6\:14/) {
#      next;
#    }

    my $output_hash = {};
    $output_hash->{'iid'} = $id;
    $self->dataflow_output_id($output_hash,1);
  }

  return 1;
}

sub db {
  my ($self, $value) = @_;
  if($value){
    $self->{'dbadaptor'} = $value;
  }
  return $self->{'dbadaptor'};
}

sub get_dba {
   my ($self,$connection_info) = @_;
   my $dba;

   if (ref($connection_info)=~m/HASH/) {

       $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                  %$connection_info,
                                                );
   }

  $dba->dbc->disconnect_when_inactive(1) ;
  return $dba;

}

sub input_ids {
 my ($self,$value) = @_;

  if (defined $value) {
    $self->{'input_ids'} = $value;
  }

  if (exists($self->{'input_ids'})) {
    return $self->{'input_ids'};
  }

  else {
    return undef;
  }

}

1;
