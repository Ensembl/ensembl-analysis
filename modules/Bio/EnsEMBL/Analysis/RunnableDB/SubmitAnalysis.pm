package Bio::EnsEMBL::Analysis::RunnableDB::SubmitAnalysis;

use strict;

use Data::Dumper;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::RunnableDB::HiveBaseRunnable');

sub fetch_input {
  my $self = shift;

  throw("Must define toplevel")
    unless(defined($self->param('toplevel')));

  throw("Must define reference db")
    unless(defined($self->param('reference_db')));


  $self->db($self->get_dba($self->param('reference_db')));

  $self->toplevel($self->param('toplevel'));

  $self->coord_system_name($self->param('coord_system_name'));

  $self->coord_system_version(
                   $self->param('coord_system_version')
                 );

  $self->input_ids(
                    $self->get_slice_names()
                  );

}



sub run {
  my $self = shift;
  return 1;
}


sub write_output {
  my $self = shift;

  foreach my $id (@{$self->input_ids()}) {
    print "WRITING ID: ".$id."\n";
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
