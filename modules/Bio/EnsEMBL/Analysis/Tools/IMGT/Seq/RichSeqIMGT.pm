
package Bio::EnsEMBL::Analysis::Tools::IMGT::Seq::RichSeqIMGT;
use strict;

use base qw(Bio::Seq::RichSeq);


sub new {
  # standard new call..
  my($caller,@args) = @_;
  my $self = $caller->SUPER::new(@args);
  
  my ($data_class) = $self->_rearrange([qw(DATA_CLASS
					    )],
					@args);

  defined $data_class and $self->data_class($data_class);

  return $self;
}


sub data_class {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_data_class'} = $value;
    }
    return $obj->{'_data_class'};

}


1;
