# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

=head1 NAME

HiveCopyGenes.pm

=head1 DESCRIPTION

This module copies the genes whose identifiers are listed in a given file from a source database into an output database.

=head1 OPTIONS

-sourcehost     source database host
-sourceuser     source database read-only user name
-sourceport     source database port
-sourcepass     source database read-only user pass
-sourcedbname   source database name
-outhost        destination database host
-outuser        destination database write user name
-outpass        destination database write user pass
-outdbname      destination database name
-outport        destination database port
-dnahost        dna database host (usually same as source)
-dnadbname      dna database name ()
-dnauser        dna database read-only user name
-dnaport        dna database port
-file           file containing the list of gene ids to be copied from the source to the destination database
--merge         The source genes will be copied into existing genes in the output database if they overlap if this parameter is true.
--logic         This option will change the analysis of the genes and transcripts being written to the output database to an analysis with the specified logic_name.
--source        This option will change the source of the genes and transcripts being written to the output database to the source specified.
--biotype       This option will change the biotype of the genes and transcripts being written to the output database to the biotype specified.
--stable_id     Flag for indicating that input file contains stable IDs and not gene IDs.
--verbose       It adds some more prints to the output.

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes -sourcehost HOST -sourceuser READONLY_USER -sourceport PORT -sourcepass READONLY_PASS -sourcedbname SOURCEDB -outhost TARGETHOST -outuser WRITE_USER -outpass WRITE_PASS -outdbname TARGETDB -outport TARGETPORT -dnahost DNAHOST -dnadbname DNADB -dnauser READONLY_USER -dnaport DNAPORT -file gene_ids_to_copy.txt

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene attach_Analysis_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);

sub param_defaults {
    return {
      copy_genes_path => '$ENSCODE/ensembl-analysis/scripts/genebuild/',
      copy_genes_script_name => 'copy_genes.pl',

      # copy_genes.pl script parameters
      source => undef,
      biotype => undef,
      stable_id => undef,
      merge => undef,
      logic => undef, # to set the genes and transcripts analysis (logic names) (optional)
      sourcehost => undef,
      sourceuser => undef,
      sourceport => '3306',
      sourcepass => undef,
      sourcedbname => undef,
      outhost => undef,
      outuser => undef,
      outpass => undef,
      outdbname => undef,
      outport => '3306',
      dnahost => undef,
      dnadbname => undef,
      dnauser => undef,
      dnaport => '3306',
      file => undef,
      copy_genes_directly => 0,
    }
}

sub fetch_input {
  my $self = shift;

  $self->create_analysis;
  return 1;
}

sub run {

  my $self = shift;

  $self->setup_fasta_db;
  if($self->param('copy_genes_directly')) {
    my $input_dba = $self->get_database_by_name('source_db');
    my $output_dba = $self->get_database_by_name('target_db');

    $self->hrdb_set_con($input_dba,'source_db');
    $self->hrdb_set_con($output_dba,'target_db');

    my $input_genes = $self->param('iid');
    $self->copy_genes_directly($input_genes);

  } else {
    $self->param_required('sourcehost');
    $self->param_required('sourceuser');
    $self->param_required('sourceport');
    $self->param_required('sourcepass');
    $self->param_required('sourcedbname');
    $self->param_required('outhost');
    $self->param_required('outuser');
    $self->param_required('outpass');
    $self->param_required('outdbname');
    $self->param_required('outport');
    $self->param_required('dnahost');
    $self->param_required('dnadbname');
    $self->param_required('dnauser');
    $self->param_required('dnaport');
    $self->param_required('file');

    my $command = "perl ".$self->param('copy_genes_path')
                         .$self->param('copy_genes_script_name')
                         ." -sourcehost ".$self->param('sourcehost')
                         ." -sourceport ".$self->param('sourceport')
                         ." -sourcedbname ".$self->param('sourcedbname')
                         ." -outuser ".$self->param('outuser')
                         ." -outpass ".$self->param('outpass')
                         ." -outdbname ".$self->param('outdbname')
                         ." -outport ".$self->param('outport')
                         ." -outhost ".$self->param('outhost')
                         ." -dnahost ".$self->param('dnahost')
                         ." -dnadbname ".$self->param('dnadbname')
                         ." -dnauser ".$self->param('dnauser')
                         ." -dnaport ".$self->param('dnaport')
                         ." -file ".$self->param('file');

  if ($self->param('logic')) {
  	$command .= " -logic ".$self->param('logic');
  }

  if ($self->param('source')) {
    $command .= " -source ".$self->param('source');
  }

  if ($self->param('biotype')) {
    $command .= " -biotype ".$self->param('biotype');
  }

  if ($self->param('stable_id')) {
    $command .= " -stable_id ";
  }

  if ($self->param('merge')) {
    $command .= " -merge ";
  }

  print run_command($command,"Copying genes...");
  }
  return 1;
}

sub write_output {
  my $self = shift;

  if($self->param('copy_genes_directly')) {
    my $output_db = $self->hrdb_get_con('target_db');
    my $output_ga = $output_db->get_GeneAdaptor();
    my $output_genes = $self->output_genes;
    unless(scalar(@{$output_genes})) {
      $self->throw("You have selected to copy genes directly based on the feature id but no genes were present in the output array, so something has went wrong");
    }

    my $copy_biotypes_to_ignore = $self->param('copy_biotypes_to_ignore');
    my $add_analysis = $self->param_is_defined('logic_name') ? $self->param('logic_name') : undef;
    my $add_biotype = $self->param_is_defined('biotype') ? $self->param('biotype') : undef;
    foreach my $gene (@{$output_genes}) {
      if($copy_biotypes_to_ignore && $copy_biotypes_to_ignore->{$gene->biotype}) {
        next;
      }
      if ($add_analysis) {
        attach_Analysis_to_Gene($gene, $self->analysis);
      }
      if ($add_biotype) {
        $gene->biotype($add_biotype);
        foreach my $transcript (@{$gene->get_all_Transcripts}) {
          $transcript->biotype($add_biotype);
        }
      }
      empty_Gene($gene);
      $output_ga->store($gene);
    }
  }

  return 1;
}

sub copy_genes_directly {
  my ($self,$input_genes) = @_;

  my $input_db = $self->hrdb_get_con('source_db');
  my $output_genes = [];

  my $input_ga = $input_db->get_GeneAdaptor();
  foreach my $gene_id (@{$input_genes}) {
    my $gene = $input_ga->fetch_by_dbID($gene_id);
    unless($gene) {
      $self->throw("Problem loading gene with dbID ".$gene_id." from the input database!");
    }

    if($self->param('filter_on_overlap')) {
      $gene = $self->filter_on_overlap($gene);
      unless($gene) {
        next;
      }
    }
    push(@{$output_genes},$gene);
  }

  $self->output_genes($output_genes);
}


sub output_genes {
  my ($self,$output_genes) = @_;

  if($output_genes) {
    $self->param('_output_genes',$output_genes);
  }

  return($self->param('_output_genes'));
}


sub filter_on_overlap {
  my ($self,$gene) = @_;

  my $overlap_filter_type = $self->param_required('overlap_filter_type');
  my $filter_on_strand = $self->param_required('filter_on_strand');
  my $filter_against_biotypes = {};
  if($self->param('filter_against_biotypes')) {
    $filter_against_biotypes = $self->param('filter_against_biotypes');
  }

  my $output_db = $self->hrdb_get_con('target_db');
  my $output_ga = $output_db->get_GeneAdaptor();
  my $output_db_slice = get_feature_slice_from_db($gene,$output_db);
  my @output_db_genes = @{$output_ga->fetch_all_by_Slice($output_db_slice,undef,1)};

  unless(scalar(@output_db_genes)) {
    return($gene);
  }

  foreach my $output_gene (@output_db_genes) {
    unless($filter_against_biotypes->{$output_gene->biotype}) {
      next;
    }

    if($overlap_filter_type eq 'genomic_overlap') {
      if($filter_on_strand && $output_gene->strand == $gene->strand) {
	next;
      }
      if($gene->seq_region_start <= $output_gene->seq_region_end && $output_gene->seq_region_start <= $gene->seq_region_end) {
	return(0);
      }
    } elsif($overlap_filter_type eq 'exon_overlap') {
      if($filter_on_strand && $output_gene->strand == $gene->strand) {
	next;
      }
      if(exon_overlap($gene,$output_gene)) {
	return(0);
      }
    } elsif($overlap_filter_type eq 'coding_exon_overlap') {
      if($filter_on_strand && $output_gene->strand == $gene->strand) {
	next;
      }
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
	foreach my $output_transcript (@{$output_gene->get_all_Transcripts}) {
          if(coding_exon_overlap($transcript,$output_transcript)) {
            return(0);
          }
        }
      }
    } else {
      $self->throw("You have selected an unsupported overlap filter. Filter selected: ".$overlap_filter_type);
    }
  }

  return($gene);
}


sub get_feature_slice_from_db {
  my ( $feature, $db ) = @_;

  # This little helper routine returns a feature slice for a particular
  # region.  The slice will be associated with the given database.

  my $slice = $feature->feature_Slice();

  my @slices = @{
    $db->get_SliceAdaptor()->fetch_by_region_unique(
         $slice->coord_system_name(), $slice->seq_region_name(),
         $slice->start(),             $slice->end(),
         1,            $slice->coord_system()->version(),
         1 ) };

  if ( scalar(@slices) != 1 ) {
    # This will hopefully only happen if the Primary and Secondary
    # databases contain different assemblies.
    die( "!! Problem with projection for feature slice %s\n",
         $slice->name() );
  }

  return $slices[0];
}


1;
