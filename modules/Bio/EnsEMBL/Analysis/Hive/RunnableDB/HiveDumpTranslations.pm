=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpTranslations

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpTranslations;


use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpTranslations
  Function  : Dump translations of a database via the script ~/enscode/ensembl-analysis/scripts/protein/dump_translations.pl 
  file
  Returntype: Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpTranslations 
  Exceptions: 
  Example   : 

=cut



use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);

sub param_defaults {
    return {
      dump_translations_path => '~/enscode/ensembl-analysis/scripts/protein/',
      dump_translations_name => 'dump_translations.pl',

      # dump_translations.pl script parameters
      # logic => '', # to set the genes and transcripts analysis (logic names) (optional)
      sourcehost => undef,
      sourceuser => undef,
      sourceport => '3306',
      sourcepass => undef,
      sourcedbname => undef,
      dnahost => undef,
      dnadbname => undef,
      dnauser => undef,
      dnaport => '3306',
      file => undef
    }
}

sub fetch_input {
  my $self = shift;

  $self->param_required('sourcehost');
  $self->param_required('sourceuser');
  $self->param_required('sourceport');
  $self->param_required('sourcepass');
  $self->param_required('sourcedbname');
  $self->param_required('dnahost');
  $self->param_required('dnadbname');
  $self->param_required('dnauser');
  $self->param_required('dnaport');
  $self->param_required('file'); 
  return 1;
}

sub run {

  my $self = shift;
  my $command = "perl ".$self->param('dump_translations_path')
                       .$self->param('dump_translations_name')
                       ." -dbhost ".$self->param('sourcehost')
                       ." -dbport ".$self->param('sourceport')
                       ." -dbname ".$self->param('sourcedbname')
                       ." -dbuser ".$self->param('sourceuser')
                       ." -dnadbhost ".$self->param('dnahost')
                       ." -dnadbname ".$self->param('dnadbname')
                       ." -dnadbuser ".$self->param('dnauser')
                       ." -dnadbport ".$self->param('dnaport')
                       . " -verbose -db_id ";

  $command .= " -dbpass ".$self->param('sourcepass') if ($self->param_is_defined('sourcepass'));
  $command .= " -dnadbpass ".$self->param('dnapass') if ($self->param_is_defined('dnapass'));
  # if ($self->param('logic')) {
  #   $command .= " -logic ".$self->param('logic');
  # }
  my $file = $self->param('file');
  if ($self->param('iid')) {
    $command .= " -slicename ".$self->param('iid');
    $file .= '_'.$self->param('iid').'.fasta';
  }
  $command .= " -file ".$file;
  
  print run_command($command,"Dumping translations...");
  $self->output([$file]);
  return 1;
}

sub write_output {
  my $self = shift;

  foreach my $file (@{$self->output}) {
    $self->dataflow_output_id({fasta_file => $file}, 2);
  }
  return 1;
}

1;






