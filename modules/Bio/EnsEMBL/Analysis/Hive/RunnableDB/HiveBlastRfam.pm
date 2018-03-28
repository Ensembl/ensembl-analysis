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

Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam->
  new(
      -analysis => $analysis,
      -db => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

Modified blast runnable for specific use with RFAMSEQ.
Use for running BLASTN of genomic vs RFAMSEQ prior to 
ncRNA analysis using Infernal.
Slice size seems best around 200k

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRfam;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Runnable::BlastRfam;
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Function  : fetch sequence out of database, instantiate the filter, 
  parser and finally the blast runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;

  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                             -program_file => $self->param('blast_exe_path'),
                                             -db_file => $self->param('blast_db_path'),
                                             -parameters => $self->param('commandline_params'),
                                           );
  $self->analysis($analysis);

  my $repeat_dba = $self->hrdb_get_dba($self->param('repeat_db'));
  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));

  if($dna_dba) {
    $repeat_dba->dnadb($dna_dba);
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($repeat_dba,'repeat_db');
  $self->hrdb_set_con($output_dba,'output_db');

  my $repeat_logic_names = $self->param('repeat_logic_names');

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
    my $slice = $self->fetch_sequence($input_id,$repeat_dba,$repeat_logic_names);

   unless ($slice->seq =~ /[CATG]{3}/) {
    say "The following slice is > 3bp after applying repeatmasking, will skip:";
    say $slice->name;
    next;
  }

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastRfam->new
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
    $runnable->timer($self->param('timer'));
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
