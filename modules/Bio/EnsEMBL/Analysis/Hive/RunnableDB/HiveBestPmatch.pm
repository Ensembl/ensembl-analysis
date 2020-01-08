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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestPmatch

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestPmatch;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(empty_Object);
use Bio::EnsEMBL::Analysis::Runnable::BestPmatch;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my ($self) = @_;

  $self->create_analysis;
  my @features;
  my $dbs;
  if (ref($self->INPUT_DB) eq 'ARRAY') {
      $dbs = $self->INPUT_DB;
  }
  else {
      $dbs = [$self->INPUT_DB];
  }
  foreach my $dba (@$dbs) {
      my $pafa = $self->hrdb_get_dba($dba)->get_ProteinAlignFeatureAdaptor;

      if (ref($self->PMATCH_LOGIC_NAME) eq 'ARRAY') {
         foreach my $logic_name (@{$self->PMATCH_LOGIC_NAME}) {
            my $f = $pafa->fetch_all_by_logic_name($logic_name);
            print "Have fetched ", scalar(@$f), " with logic_name : $logic_name from ".$pafa->dbc->dbname."\n";
            push(@features, @$f);
         }
      } else {
        @features = @{$pafa->fetch_all_by_logic_name($self->PMATCH_LOGIC_NAME)} ;
        print "Have fetched ".@features." with logic_name : ".$self->PMATCH_LOGIC_NAME." from ".$pafa->dbc->dbname."\n";
      }
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BestPmatch->new(
        -protein_hits => \@features,
        -min_coverage => $self->MIN_COVERAGE,
        -analysis => $self->analysis,
       );
  $self->runnable($runnable);
}

sub run {
  my ($self) = @_;

  my @output;
  foreach my $runnable (@{$self->runnable}) {
    $runnable->run;
    foreach my $obj (@{$runnable->output}) {
      empty_Object($obj);
      push(@output, $obj);
    }
  }
  $self->output(\@output);
}


sub get_adaptor {
  my ($self) = @_;

  my $output_db = $self->hrdb_get_dba($self->OUTPUT_DB);
  return $output_db->get_ProteinAlignFeatureAdaptor;
}


sub PMATCH_LOGIC_NAME {
  my ($self, $arg) = @_;

  if ($arg) {
    $self->param('PMATCH_LOGIC_NAME', $arg);
  }
  return $self->param('PMATCH_LOGIC_NAME');
}


sub OUTPUT_DB {
  my ($self, $arg) = @_;

  if ($arg) {
    $self->param('target_db', $arg);
  }
  return $self->param('target_db');
}


sub INPUT_DB {
  my ($self, $arg) = @_;

  if ($arg) {
    $self->param('source_db', $arg);
  }
  return $self->param('source_db');
}


sub MIN_COVERAGE {
  my ($self, $arg) = @_;

  if ($arg) {
    $self->param('MIN_COVERAGE', $arg);
  }
  return $self->param('MIN_COVERAGE');
}


1;
