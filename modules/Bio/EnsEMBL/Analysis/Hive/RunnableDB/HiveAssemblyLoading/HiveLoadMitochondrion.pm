#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion;

use strict;
use warnings;
use feature 'say';
use File::Fetch;
use File::Spec::Functions;
use File::Path qw(make_path);

use Bio::DB::EUtilities;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters for the module
                ncbi_check => 1, # If mt_accession is not set, search for the RefSeq annotation for the species
                mt_accession => '',
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my $self = shift;

  return {
    %{$self->SUPER::param_defaults},
    ncbi_check => 1,
    mt_accession => '',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Retrieves the mitochondrial genome Genbank file and its RefSeq annotation
              based on the accession provided with 'mt_accession' (usually NC_XXXXXXX).
              If the accession is not provided and 'ncbi_check' is true, it will query
              NCBI to find the accession for the species or a subspecies.
              If multiple entries are found and there is non for the expected taxonomy id,
              it uses the accession with the lowest taxonomy id as it is expected to be the
              lowest in the tree, i.e not a sub species
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my $self = shift;

  if($self->param_is_defined('skip_analysis') && $self->param('skip_analysis')) {
    $self->complete_early('The skip_analysis flag has been set to 1, skipping loading');
  }

  $self->param_required('target_db');
  $self->param_required('enscode_root_dir');
  $self->param_required('species_name');
  my $target_db = $self->get_database_by_name('target_db');

  if (!$self->param('mt_accession') and $self->param('ncbi_check')) {
    my $scientific_name = $target_db->get_MetaContainerAdaptor->get_scientific_name;
    my $taxon_id = $target_db->get_MetaContainerAdaptor->get_taxonomy_id;
    my $email = $ENV{HIVE_EMAIL} || 'genebuild@ebi.ac.uk';
    my $factory = Bio::DB::EUtilities->new(
      -eutil => 'esearch',
      -db     => 'nuccore',
      -term   => '"'.$scientific_name.'"[Organism] AND biomol_genomic[PROP] AND refseq[filter] AND mitochondrion[filter]',
      -email  => $email,
      -retmax => 5,
    );
    if ($factory->get_count) {
      $self->say_with_header(join(' ID: ', $factory->get_ids));
      my @uids = $factory->get_ids;

      $factory->reset_parameters(
        -eutil => 'esummary',
        -db    => 'nuccore',
        -email  => $email,
        -id    => \@uids
      );

      my %mt_accessions;
      while (my $ds = $factory->next_DocSum) {
        $self->say_with_header('ID: '.$ds->get_id.' ACCESSION: '.join(' ', $ds->get_contents_by_name('AccessionVersion')));
        my ($mt_accession) = $ds->get_contents_by_name('AccessionVersion');
        my ($tax_id) = $ds->get_contents_by_name('TaxId');
        $mt_accessions{$tax_id} = $mt_accession;
      }
      if (exists $mt_accessions{$taxon_id}) {
        $self->param('mt_accession', $mt_accessions{$taxon_id});
      }
      else {
        my @tax_ids = sort {$a <=> $b} keys %mt_accessions;
        $self->param('mt_accession', $mt_accessions{$tax_ids[0]});
      }
    }
  }
  if ($self->param('mt_accession')) {
    unless($self->param_is_defined('chromosomes_present')) {
      my $cs_adaptor = $target_db->get_CoordSystemAdaptor;
      my $cs_rank1 = $cs_adaptor->fetch_by_rank(1);
      if ($cs_rank1->name eq 'primary_assembly' or $cs_rank1->name eq 'chromosome') {
        $self->param('chromosomes_present', $cs_rank1->name);
      }
    }

    if (!-e $self->param_required('output_path')) {
      make_path($self->param('output_path'));
    }

    chdir($self->param('output_path'));
    my $fetcher = File::Fetch->new(uri => 'http://www.ncbi.nlm.nih.gov/nuccore/'.$self->param('mt_accession'));
    $self->param('mt_filename', $fetcher->fetch());
  }
  else {
    $self->complete_early('No mitochondria for this species');
  }
}

sub run {
  my $self = shift;

  my $species_name = $self->param('species_name');
  my $output_path = $self->param('output_path');
  my $enscode_dir = $self->param('enscode_root_dir');

  my $mt_filename = $self->param('mt_filename');
  my $mt_accession = $self->param('mt_accession');
  my $toplevel = 'chromosome';
  my $chromosome_flag = "";
  my $scaffold_flag = "";
  if($self->param('chromosomes_present')) {
    $toplevel = 'primary_assembly' if ($self->param('chromosomes_present') eq 'primary_assembly');
    $chromosome_flag = " -chromosome  MT";
    $scaffold_flag = " -scaffold ".$mt_accession;
  } else {
    $toplevel = "scaffold";
    $scaffold_flag = " -scaffold MT";
  }

  my $target_db = $self->param('target_db');
  my $user = $target_db->{'-user'};
  my $pass = $target_db->{'-pass'};
  my $host = $target_db->{'-host'};
  my $port = $target_db->{'-port'};
  my $dbname = $target_db->{'-dbname'};

  my $cmd =  'perl '.catfile($enscode_dir, 'ensembl-analysis', 'scripts', 'refseq', 'load_mitochondria.pl').
             " -dbhost ".$host.
             " -dbuser ".$user.
             " -dbport ".$port.
             " -dbpass ".$pass.
             " -dbname ".$dbname.
             $chromosome_flag.
             " -name MT".
             $scaffold_flag.
             " -toplevel ".$toplevel.
             " -gene_type protein_coding".
             " -trna_type Mt_tRNA".
             " -rrna_type Mt_rRNA".
             " -logic_name mt_genbank_import".
             " -source RefSeq".
             " -accession ".$mt_accession.
             " -genbank_file ".$mt_filename.
             " -non_interactive".
             " -download";

  my $return = system($cmd);
  if($return) {
    $self->throw("The load_mitochondria script returned a non-zero exit code. Commandline used:\n".$cmd);
  }

  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
