=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

  {
    -logic_name => 'HiveDumpTranslations',
    -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpTranslations',
    -parameters => {
      dump_translations_script => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'protein', 'dump_translations.pl'),
      source_db => $self->o('lincRNA_output_db'),
      dna_db => $self->o('dna_db'),
      file => $self->o('file_translations') ,
    },
  },

=head1 DESCRIPTION

It dumps the sequence of each translation into a file. It fetches the
genes from a database specified by 'source_db'. If 'source_db' does not
have DNA, you need to specify 'dna_db' which contains DNA. If 'iid' which
should be an Ensembl slice name is given, only the translations from this
region would be dumped.

The output file is specified in 'file'. If 'iid' is used, the filename would
be "'file'_'iid'.'file_extension'"

The filename will be send on branch '_branch_to_flow_to' which is 2 by default
under 'fasta_file'.

=head1 METHODS

param_defaults
fetch_input
run
write_output

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpTranslations;


use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Set default parameters:
               file_extension => 'fasta', # only used when 'iid' is given
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    file_extension => 'fasta',
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Check 'source_db', 'dump_translations_script' and 'file'.
              It creates the command line to be run. It adds a dna_db if
              'dna_db' is set and uses 'iid' as slicename if it is set
 Returntype : None
 Exceptions : Throws is 'source_db' or 'dump_translations_script' or 'file' is not set

=cut

sub fetch_input {
  my $self = shift;

  my $source_db = $self->param_required('source_db');
  my $script = $self->param_required('dump_translations_script');
  my $file = $self->param_required('file');
  my $command = "perl $script ".
                ' -dbhost '.$source_db->{'-host'}.
                ' -dbport '.$source_db->{'-port'}.
                ' -dbuser '.$source_db->{'-user'}.
                ' -dbname '.$source_db->{'-dbname'};
  $command .= ' -dbpass '.$source_db->{'-pass'} if (exists $source_db->{'-pass'} and $source_db->{'-pass'});
  if ($self->param_is_defined('dna_db')) {
    my $dna_db = $self->param('dna_db');
    $command .= ' -dnadbhost '.$dna_db->{'-host'}.
                ' -dnadbport '.$dna_db->{'-port'}.
                ' -dnadbuser '.$dna_db->{'-user'}.
                ' -dnadbname '.$dna_db->{'-dbname'};
    $command .= ' -dnadbpass '.$dna_db->{'-pass'} if (exists $dna_db->{'-pass'} and $dna_db->{'-pass'});
  }
  if ($self->param_is_defined('iid')) {
    $command .= ' -slicename '.$self->param('iid');
    $file .= '_'.$self->param('iid').'.'.$self->param('file_extension');
  }
  $command .= " -file ".$file;
  if ($self->param_is_defined('db_id'))
  {
  	$command .= " -db_id"; 
  }
  if ($self->param_is_defined('stable_id'))
  {
  	$command .= " -stable_id"; 
  }
  if ($self->param_is_defined('biotype')) {
  	$command .= " -biotype " .$self->param('biotype');
  }
  
  $self->output([$file]);
  print $command, "\n";
  $self->param('command', $command);
}


=head2 run

 Arg [1]    : None
 Description: Execute the 'dump_translations_script'
 Returntype : None
 Exceptions : None

=cut

sub run {
  my $self = shift;

  execute_with_wait($self->param('command'));
}


=head2 write_output

 Arg [1]    : None
 Description: Write to '_branch_to_flow_to' the name of the file under 'fasta_file'
              if the file exists and is not empty
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

  foreach my $file (@{$self->output}) {
    $self->dataflow_output_id({fasta_file => $file}, $self->param('_branch_to_flow_to'))
      if (-e $file and -s $file);
  }
}

1;
