package Bio::EnsEMBL::Analysis::RunnableDB::HiveBaseRunnable;

use strict;
use Carp;
use Bio::EnsEMBL::Hive::Utils ('stringify');

use parent ('Bio::EnsEMBL::Hive::Process','Bio::EnsEMBL::Analysis::RunnableDB');

=head2 _slurp

Reads the whole content of a file and returns it as a string

=cut

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
