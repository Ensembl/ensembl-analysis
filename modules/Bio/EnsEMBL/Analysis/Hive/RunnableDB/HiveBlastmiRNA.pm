=head1 LICENSE

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastmiRNA;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable::BlastmiRNA;
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



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

  #add dna_db

  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                             -program_file => $self->param('blast_exe_path'),
                                             -db_file => $self->param('blast_db_path'),
                                             -parameters => $self->param('commandline_params'),
                                           );
  $self->analysis($analysis);


  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_db;
  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }

    setup_fasta(
                 -FASTA => $self->param('genome_file'),
	       );
  } elsif($self->param('dna_db')) {
    $dna_db = $self->get_database_by_name('dna_db');
    $output_dba->dnadb($dna_db);
  } else {
    $self->throw("You must provide either a flatfile or a dna db to read from");
  }

  $self->hrdb_set_con($output_dba,'output_db');

  my $blast_params = $self->param('BLAST_PARAMS');
  my $blast_filter = $self->param('BLAST_FILTER');
  my $blast_parser = $self->param('BLAST_PARSER');
  my $filter_params = $self->param('FILTER_PARAMS');
  my $parser_params = $self->param('PARSER_PARAMS');

  my %blast = %{$blast_params};

  my $parser = $blast_parser->new(%$parser_params);

  my $filter;
  if($blast_filter){
    $filter = $blast_filter->new(%$filter_params);
  }

  my $input_ids = $self->param('iid');
  foreach my $input_id (@{$input_ids}) {
    # Note that we used to load dust repeats in older versions of the code. But this put a large strain on the servers at scale
    # and is not easy to make load from a flatfile without rewriting a reasonable amount of code. Given that our classifier
    # will later filter on dumped repeats, there is no real need for this and so now the code just fetches an unmasked slice
    my $slice = $self->fetch_sequence($input_id);
    unless ($slice->seq =~ /[CATG]{3}/) {
      say "The following slice is > 3bp after applying repeatmasking, will skip:";
      say $slice->name;
      next;
    }

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastmiRNA->new
    (
     -query    => $slice,
     -program  => $self->analysis->program_file,
     -parser   => $parser,
     -filter   => $filter,
     -database => $self->analysis->db_file,
     -analysis => $self->analysis,
     -options  => $self->param('commandline_params'),
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
       $self->throw("Error running BLAST:\nSlice: ".$runnable->query->name."\nError: ".$error);
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
  my $daf_adaptor = $out_dba->get_DnaAlignFeatureAdaptor;
  foreach my $hit ( @{$self->output} ) {
    $daf_adaptor->store($hit);
  }

  return 1;
} ## end sub write_output

1;
