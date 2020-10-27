=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCheckInitialData

check if data for this annotation 

=cut


package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCheckInitialData;

use strict;
use warnings;
use feature 'say';


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    # skip_initial_data_check => 0, # if set to 1, it will skip checks 

  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: fetch the main variable that will be checked
 Returntype : None
 Exceptions : 

=cut

sub fetch_input {
  my $self = shift;

print "DEBUG::" . $self->param_is_defined('skip_initial_data_check') . "\n"; 
  if ($self->param_is_defined('skip_initial_data_check') == 20) {
    $self->complete_early('You asked to skip this analysis');
  }

}


=head2 run

 Arg [1]    : None
 Description: Checks there are enough data and if all files if they are available
 Returntype : None
 Exceptions : 

=cut

sub run {
  my ($self) = @_;
  my $lr_csv = $self->param('lr_csv');
  my $lr_gn_csv = $self->param('lr_gn_csv');
  my $rnaseq_csv = $self->param('rnaseq_csv');
  my $rnaseq_gn_csv = $self->param('rnaseq_gn_csv');
  
  if ( ($self->param_is_defined('skip_rnaseq') == 1 ) 
   and ($self->param_is_defined('skip_long_read') == 1 ) 
   and ($self->param_is_defined('skip_projection') == 1 ) ) {
    $self->throw("You will have very few resources to build models. Set skip_initial_data_check to 1 
      and continue only if you know what you are doing");
  }
   
  if(!-e $lr_csv) {
    say $lr_csv.' does not exist'; 
    if ($self->param_is_defined('skip_long_read') == 0 ) {
    	say 'This is not good'; 
    }
  }
  
  if(!-e $lr_gn_csv) {
    say $lr_gn_csv.' does not exist';
  }  

  if(!-e $rnaseq_csv) {
    say $rnaseq_csv.' does not exist';
    if ($self->param_is_defined('skip_rnaseq') == 0 ) {
      $self->throw("You thought there are rnaseq data, but there are not. Set skip_initial_data_check to 1 
        and continue only if you know what you are doing"); 
    }
  }
  
  if(!-e $rnaseq_gn_csv) {
    say $rnaseq_gn_csv.' does not exist';
  }  
  
}


sub write_output {
  my ($self) = @_;
}
1; 