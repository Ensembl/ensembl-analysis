=head1 LICENSE

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::BlastmiRNA - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastmiRNA->new
     (
      -analysis => $analysis,
      -db       => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

Modified blast runnable for specific use with miRNA
Use for running BLASTN of genomic sequence vs miRNAs prior to 
miRNA anaysis
Slice size seems best around 200k

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastPepToGenome;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;
use Bio::EnsEMBL::DnaPepAlignFeature;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlast');



=head2 fetch_input

  Arg [1]   : None
  Function  : fetch sequence out of database, instantiate the filter,
            : parser and finally the blast runnable
  Returntype: None
  Exceptions: none
  Example   : $blast->fetch_input;

=cut

sub fetch_input{
  my ($self) = @_;

  $self->create_analysis();
  $self->analysis->program($self->param('blast_program')) if ($self->param_is_defined('blast_program'));
  $self->analysis->program_file($self->param('blast_exe_path')) if ($self->param_is_defined('blast_exe_path'));
  $self->analysis->parameters($self->param('commandline_params')) if ($self->param_is_defined('commandline_params'));
  $self->analysis->db_file($self->param('blast_db_path')) if ($self->param_is_defined('blast_db_path'));
  $self->analysis->db($self->param('blast_db_name')) if ($self->param_is_defined('blast_db_name'));

  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));

  if($dna_dba) {
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'output_db');

  my $blast_params = $self->param('BLAST_PARAMS');
  my $blast_filter = $self->param('BLAST_FILTER');
  my $blast_parser = $self->param('BLAST_PARSER');
  my $filter_params = $self->param('FILTER_PARAMS');
  my $parser_params = $self->param('PARSER_PARAMS');

  my %blast = %{$blast_params};

  my $parser;
  if($parser_params) {
    $parser = $self->make_parser($parser_params);
  }

  my $filter;
  if($blast_filter){
    $filter = $self->make_filter($filter_params);
  }


  # submit blast module to use via analysis_parameters column of analysis table
  my $options_string ;
  my %options = %{$parser_params};

  my $query_seqs = $self->get_query_seqs($self->input_id);

  foreach my $seq (@{$query_seqs}) {
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Blast->new(
                     -query => $seq,
                     -program => $self->param('blast_exe_path'),
                     -database => $self->param('blast_db_path'),
                     -parser => $parser,
                     -filter => $filter,
                     -analysis => $self->analysis,
                     -options => $self->param('commandline_params'),
                     %blast,
                   );
    $self->runnable($runnable);
  }

  return 1;
}


sub run {
  my ($self) = @_;
  foreach my $runnable(@{$self->runnable}){
    eval {
      $runnable->run;
    }; if ($@) {
      my $error  = $@;
      if($error =~ /VOID/) {
        say "\nError from short sequence, this is okay";
      } else {
       $self->throw("Error running BLAST: ".$error);
      }
    }
    $self->output($runnable->output());
    # This code was added in to deal with File::Temp having issues with having too many files open at once. The files are
    # only deleted when the object is removed. As this module is currently batched on 5MB of 200KB slices, it means occasionally
    # you might get a batch of lots if tiny slices. If the number of tiny slices >= 8185 then the job will die because of File::Temp
    undef($runnable);
  }
  return 1;
}


sub write_output {
  my ($self) = @_;

  # write genes out to a different database from the one we read genes from.
  my $out_dba = $self->hrdb_get_con('output_db');
  my $paf_adaptor = $out_dba->get_ProteinAlignFeatureAdaptor;
  my $slice_adaptor = $out_dba->get_SliceAdaptor;
  foreach my $hit ( @{$self->output} ) {
    my $slice_name = $hit->{'SLICE_NAME'};
    my $slice = $slice_adaptor->fetch_by_name($slice_name);
    $hit->analysis($self->analysis);
    $hit->slice($slice);
    say "FM2 storing hitname: ".$hit->hseqname;
    say "FM2 storing slicename orig: ".$slice_name;
    say "FM2 storing slicename: ".$slice->name;
    $paf_adaptor->store($hit);
  }

  return 1;
} ## end sub write_output


sub get_query_seqs {
  my ($self,$accession_array) = @_;


  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));

  my $biotypes_hash = {};
  my @seqs;

  foreach my $accession (@{$accession_array}) {
    my $db_row = $table_adaptor->fetch_by_dbID($accession);
    unless($db_row) {
      $self->throw("Did not find an entry in the uniprot_sequences table matching the accession. Accession:\n".$accession);
    }

    push(@seqs, Bio::Seq->new(-id => $accession, -seq => $db_row->{seq}));
    $biotypes_hash->{$accession} = $db_row->{'biotype'};
  }

  $self->get_biotype($biotypes_hash);

  return \@seqs;
}


sub get_biotype {
  my ($self,$biotype_hash) = @_;
  if($biotype_hash) {
    $self->param('_biotype_hash',$biotype_hash);
  }
  return($self->param('_biotype_hash'));
}

1;
