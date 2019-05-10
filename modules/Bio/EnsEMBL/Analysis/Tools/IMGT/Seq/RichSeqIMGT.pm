# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Tools::IMGT::Seq::RichSeqIMGT;
use warnings ;
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
