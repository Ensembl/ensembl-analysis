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

check if data are enough for this annotation 

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
  }
}


=head2 run

 Arg [1]    : None
 Description: It checks if there are enough data and if all files they are available
 Returntype : None
 Exceptions : Throws if you may not have enough data for the genebuild

=cut

sub run {
  my $self = shift;
  
  if ($self->param('skip_initial_data_check') == 1) { 
    $self->complete_early('You asked to skip this analysis'); 
  }
  if ( ($self->param('clade') eq "mammals_basic") or ($self->param('clade') eq "primates_basic") ) { 
  	$self->complete_early('Your clade setting is ok for this genebuild'); 
  }
    
  my $lr_csv = $self->param('lr_csv');
  my $lr_gn_csv = $self->param('lr_gn_csv');
  my $rnaseq_csv = $self->param('rnaseq_csv');
  my $rnaseq_gn_csv = $self->param('rnaseq_gn_csv');
  
  if ( ($self->param('skip_rnaseq') == 1 ) 
  and ($self->param('skip_long_read') == 1 ) 
  and ($self->param('skip_projection') == 1 ) ) {
    $self->throw("You have very few resources to build models. Set skip_initial_data_check to 1 
      and run again if you know what you are doing!");
  }
  
  my $lr_csv_exists = 1; 
  if(!-e $lr_csv) {
    say $lr_csv.' does not exist'; 
    $lr_csv_exists = 0; 
  }

  my $lr_gn_csv_exists = 1; 
  if(!-e $lr_gn_csv) {
    say $lr_gn_csv.' does not exist';
    $lr_gn_csv_exists = 0;
  }  

  if ( ($lr_gn_csv_exists+$lr_csv_exists <1) and ( $self->param_is_defined('skip_long_read') == 0 ) ) {
    say 'Update or check your skip_long_read parameter. It looks like there are no data. \n'; 
  } elsif ( ($lr_gn_csv_exists+$lr_csv_exists >1) and ( $self->param_is_defined('skip_long_read') == 1 ) ) {
  	$self->throw("There are data. Why to skip it?"); 
  }

  my $rnaseq_csv_exists = 1; 
  if (!-e $rnaseq_csv) { 
    say $rnaseq_csv.' does not exist';
    $rnaseq_csv_exists = 0; 
  } 

  my $rnaseq_gn_csv_exists = 1; 
  if (!-e $rnaseq_gn_csv) {
    $rnaseq_gn_csv_exists = 0; 
    say $rnaseq_gn_csv.' does not exist';
  }  

  if ( ($rnaseq_csv_exists+$rnaseq_gn_csv_exists <1) and ( $self->param_is_defined('skip_rnaseq') == 0 ) ) {
    say 'Update or check your skip_rnaseq parameter. It looks like there are no data. \n'; 
  }
  
  if ($lr_gn_csv_exists+$lr_csv_exists+$rnaseq_csv_exists+$rnaseq_gn_csv_exists < 1) {
  	$self->throw("There are not enough data. Missing data. "); 
  }

} 

sub write_output {
  my ($self) = @_;
}

1; 

