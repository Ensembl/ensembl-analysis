=head1 LICENSE

Copyright [2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::FTPSetup

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::FTPSetup;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Registry;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use base qw(Bio::EnsEMBL::Production::Pipeline::FASTA::Base);

=head2 fetch_input

    Description : Implements fetch_input() interface method of Bio::EnsEMBL::Hive::Process that is used to read in parameters and load data.
                  Here we have nothing to do.

=cut

sub fetch_input {
}

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
  my ($self) = @_;

  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($self->param('reg_conf'));
  my $dba = $registry->get_DBAdaptor($self->param('species'), "core");

  my $gca_sql = 'select meta_value from meta where meta_key="assembly.accession"';
  my $gca = $dba->dbc()->sql_helper()->execute_single_result( -SQL => $gca_sql );

  my $repbase_sql = 'select meta_value from meta where meta_key="repeat.analysis" and meta_value like "%repeatmask_repbase_%"';
  my $repbase_library = $dba->dbc()->sql_helper()->execute_single_result( -SQL => $repbase_sql );
  $repbase_library =~ s/repeatmask_repbase_//;

  my $repeatmodeler_library = "/ebi/ftp/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/".$self->param('species')."/".$gca.".repeatmodeler.fa";
  my $softmasked_genome_file = "/hps/nobackup2/production/ensembl/genebuild/production/*/".$self->param('species')."/".$gca."/genome_dumps/".$self->param('species')."_softmasked_toplevel.fa";

  $self->dataflow_output_id(
			    {
        species                => $self->param('species'),
        gca                    => $gca,
        repbase_library        => $repbase_library,
        repeatmodeler_library  => $repeatmodeler_library,
	softmasked_genome_file => $softmasked_genome_file,
    },
  1);

}

=head2 write_output

    Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
                 Here we have nothing to do.

=cut

sub write_output {
}


1;
