package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB;

use strict;
use Carp;
use Bio::EnsEMBL::Hive::Utils ('stringify');

use parent ('Bio::EnsEMBL::Hive::Process','Bio::EnsEMBL::Analysis::RunnableDB');

=head2 _slurp

Reads the whole content of a file and returns it as a string

=cut

sub parse_hive_input_id {
  my $self = shift;
  my $input_id_string = $self->Bio::EnsEMBL::Hive::Process::input_id;
  unless($input_id_string =~ /.+\=\>.+\"(.+)\"/) {
    throw("Could not parse the value from the input id. Input id string:\n".$input_id_string);
  }

  $input_id_string = $1;
  return($input_id_string);

}

sub _slurp {
  my ($self, $file_name) = @_;
  my $slurped;
  {
    local $/ = undef;
    open(my $fh, '<', $file_name) or $self->throw("Couldnt open file [$file_name]");
    $slurped = <$fh>;
    close($fh);
  }
  return $slurped;
}

1;
