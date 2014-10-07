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

#  throw("Must define toplevel")
#    unless(defined($self->param('toplevel')));

#  throw("Must define reference db")
#    unless(defined($self->param('reference_db')));


  $self->db($self->get_dba($self->param('reference_db')));

#  $self->toplevel($self->param('toplevel'));

#  $self->coord_system_name($self->param('coord_system_name'));

#  $self->coord_system_version(
#                   $self->param('coord_system_version')
#                 );

#  $self->input_ids(
#                    $self->get_slice_names()
#                  );

}



sub run {
  my $self = shift;

  my $logic_name = $self->param('logic_name');
  my $slice_size = $self->param('slice_size');
  my $slice_overlap = $self->param('slice_overlap');
  my $coord_system_name = $self->param('coord_system_name');
  my $coord_system_version = $self->param('coord_system_version');
  my $slice = $self->param('slice');
  my $input_id_type = $self->param('input_id_type');
  my $file = $self->param('file');
  my $dir = $self->param('dir');
  my $regex = $self->param('regex');
  my $single = $self->param('single');
  my $translation_id = $self->param('translation_id');
  my $seq_region_name = $self->param('seq_region_name');
  my $name = 'genome';
  my $seq_level = $self->param('seq_level');
  my $top_level = $self->param('top_level');
  my $include_non_reference = $self->param('include_non_reference');
  my $hap_pair = $self->param('hap_pair');

  if (!($slice) && !($single) && !($file) && !($translation_id) && !($hap_pair)) {
    throw("Must define input as either contig, slice, file, translation_id ".
          "single, seq_level or top_level or hap_pair");
  }

  my $input_id_factory = new Bio::EnsEMBL::Pipeline::Hive::HiveInputIDFactory
  (
   -db => $self->db(),
   -slice => $slice,
   -single => $single,
   -file => $file,
   -translation_id => $translation_id,
   -seq_level => $seq_level,
   -top_level => $top_level,
   -include_non_reference => $include_non_reference,
   -dir => $dir,
   -regex => $regex,
   -single_name => $name,
   -logic_name => $logic_name,
   -input_id_type => $input_id_type,
   -coord_system => $coord_system_name,
   -coord_system_version => $coord_system_version,
   -slice_size => $slice_size,
   -slice_overlaps => $slice_overlap,
   -seq_region_name => $seq_region_name,
   -hap_pair => $hap_pair,
  );

  $input_id_factory->generate_input_ids;
  $self->{'input_id_factory'} = $input_id_factory;

  return 1;
}


sub write_output {
  my $self = shift;

  say "FM2 Writing output";

  my $output_ids = $self->{'input_id_factory'}->input_ids();
  say "OUTPUT: ".Dumper($output_ids);
  foreach my $id (@{$output_ids}) {
    say "WRITING ID: ".$id;
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
                                                  -dbname => $$connection_info{'-dbname'},
                                                  -user => $$connection_info{'-user'},
                                                  -species => $$connection_info{'-species'},
                                                  -host => $$connection_info{'-host'},
                                                  -port => $$connection_info{'-port'},
                                                );
   }

  $dba->dbc->disconnect_when_inactive(1) ;
  return $dba;

}

sub toplevel {
  my ($self, $value) = @_;
  if($value){
    $self->{'toplevel'} = $value;
  }
  return $self->{'toplevel'};
}

sub coord_system_name {
  my ($self, $value) = @_;
  if($value){
    $self->{'coord_system_name'} = $value;
  }
  return $self->{'coord_system_name'};
}

sub coord_system_version {
  my ($self, $value) = @_;
  if($value){
    $self->{'coord_system_version'} = $value;
  }
  return $self->{'coord_system_version'};
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


sub get_slice_names{
  my ($self) = @_;

  my $csa = $self->db->get_CoordSystemAdaptor();
  my $sa = $self->db->get_SliceAdaptor();

  my $slices;
  $slices = $sa->fetch_all($self->coord_system_name,
                           $self->coord_system_version,
                           0); # left out non-reference

  my @ids;
  foreach my $slice(@$slices){
    push(@ids, $slice->name);
  }

  return \@ids;

}


sub print_params {
  my $self = shift;
}


sub getTopLevelSeqs {
  my $self = shift;
}

1;
