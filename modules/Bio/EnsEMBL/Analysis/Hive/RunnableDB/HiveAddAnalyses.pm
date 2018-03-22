=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses

=head1 SYNOPSIS

{
  -logic_name      => 'PrePipelineChecks',
  -module          => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
  -max_retry_count => 0,
  -parameters      => {
                        analyses   => $self->o('analyses'),
                        target_db  => $self->o('lincRNA_output_db'),
                        source_type => 'list',
                     },
  -rc_name         => 'default',
  -flow_into => {
    '1' => ['create_toplevel_slices'],
  },
},

=head1 DESCRIPTION

Add analyses to a database specified in 'target_db' by either giving a list
of logic_names to fetch from a different database. Or by giving a list of
hashes which represent the analysis. The hash will be fed to Bio::EnsEMBL::Analysis.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(empty_Object);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Sets default parameters
               analysis_start => 0, set the value of the first analysis index if not 0
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    analysis_start => 0,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Check that the 'target_db' is set, check 'analyses' is set, check
              that 'source_type' is set. Values allowed are 'db' and 'list'.
              If "db" is set, a 'source_db' should be provided. It will fetch
              the analysis by their logic_name found in 'analyses'.
              If "list" is set, 'analyses' should be an array of hash. These
              hashes should represent an Bio::EnsEMBL::Analysis object like:
              {
                '-logic_name' => 'ensembl',
                '-module' => Genebuilder',
              }
 Returntype : None
 Exceptions : Throws if 'source_type' is not set
              Throws if 'analyses' is not set
              Throws if 'source_db' is not set when 'source_type' is "db"
              Throws if no data was fetched from the database when 'source_type' is "db"
              Throws if elements of the list is not a hash when 'source_type' is "list"

=cut

sub fetch_input {
  my ($self) = @_;

  my $source = $self->param_required('source_type');
  my $analysis_data = $self->param_required('analyses');
  my $target_db = $self->get_database_by_name('target_db');
  $self->hrdb_set_con($target_db, 'target_db');
  if ($source eq 'db') {
    my @analyses;
    my $db = $self->get_database_by_name('source_db');
    my $analysis_adaptor = $db->get_AnalysisAdaptor;
    foreach my $analysis (@$analysis_data) {
      push(@analyses, $analysis_adaptor->fetch_by_logic_name($analysis));
    }
    if (@analyses) {
      $self->output(\@analyses);
    }
    else {
      $self->throw('Could not find any analyses in '.$db->dbc->dbname.' '.$db->dbc->host);
    }
  }
  elsif ($source eq 'list') {
    if (@$analysis_data) {
      foreach my $analysis_hash (@$analysis_data) {
        $self->throw("'$analysis_hash' is not a hashref") unless (ref($analysis_hash) eq 'HASH');
        foreach my $key (keys %$analysis_hash) {
          if ($key !~ /^-/) {
            $analysis_hash->{"-$key"} = $analysis_hash->{$key};
            delete $analysis_hash->{$key};
          }
        }
      }
      $self->param('analysis_data', $analysis_data);
    }
    else {
      $self->throw('Could not find any analyses');
    }
  }
  else {
    $self->throw("'$source' is not supported");
  }
}


=head2 run

 Arg [1]    : None
 Description: Create Bio::EnsEMBL::Analysis object based on the hash given. The hash
              should be similar to {'-logic_name' => 'name', '-module' => 'HiveTest'}
              with any parameters from Bio::EnsEMBL::Analysis.
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  if ($self->param_is_defined('analysis_data')) {
    my @analyses;
    foreach my $hash (@{$self->param('analysis_data')}) {
      push(@analyses, Bio::EnsEMBL::Analysis->new(%$hash));
    }
    $self->output(\@analyses);
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Write new analysis in the target database
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  my $target_adaptor = $self->hrdb_get_con('target_db')->get_AnalysisAdaptor;
  if ($self->param('analysis_start')) {
    my $sth = $target_adaptor->dbc->prepare('ALTER TABLE analysis AUTO_INCREMENT = '.$self->param('analysis_start'));
    $sth->execute();
    $sth->finish();
  }
  foreach my $analysis (@{$self->output}) {
    empty_Object($analysis);
    $target_adaptor->store($analysis);
  }
}

1;
