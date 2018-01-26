#!/usr/bin/env perl

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
sub
  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes - 

=head1 SYNOPSIS

my $exonerate2genes = Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes->new(
                              -db         => $refdb,
			      -analysis   => $analysis_obj,
			      -input_id => $chunk_file_name
			     );

$exonerate2genes->fetch_input();
$exonerate2genes->run();
$exonerate2genes->output();
$exonerate2genes->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
It is meant to provide the interface for mapping ESTs to the genome
sequence and writing the results as genes. By the way Exonerate is run
we do not cluster transcripts into genes and only write one transcript per gene.
we then create a dbadaptor for the target database.


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes;

use warnings ;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::KillList::KillList;
use Bio::SeqIO;
use Bio::Seq;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;


  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  if($dna_dba) {
    $dba->dnadb($dna_dba);
  }
  $self->hrdb_set_con($dba,'target_db');


  ##########################################
  # set up the target (genome)
  ##########################################

  $self->create_analysis;
  my @db_files;

  ##########################################
  # set up the query (est/cDNA/protein)
  ##########################################
  my $iid_type = $self->param('iid_type');
  my ($querys, $query_file,$chunk_number,$chunk_total, $query_seq);
  unless($iid_type) {
    $self->throw("You haven't provided an input id type. Need to provide one via the 'iid_type' param");
  }

  if($iid_type eq 'db_seq') {
    $querys = $self->get_query_seqs($self->param('iid'));
    if ($self->param('sequence_table_name') =~ /protein/) {
      $self->peptide_seq($querys->[0]->seq);
    }
    $self->calculate_coverage_and_pid($self->param('calculate_coverage_and_pid'));
    @db_files = ($self->GENOMICSEQS);
  } elsif($iid_type eq 'feature_region') {
    my $feature_region_id = $self->param('iid');
    my ($slice,$accession_array) = $self->parse_feature_region_id($feature_region_id);
    $querys = $self->get_query_seqs($accession_array);
    @db_files = ($slice);
  } elsif($iid_type eq 'feature_id') {
    my $feature_type = $self->param('feature_type');
    if($feature_type eq 'transcript') {
      my $transcript_id = $self->param('iid');
      my $transcript_dba = $self->hrdb_get_dba($self->param('transcript_db'));
      if($dna_dba) {
        $transcript_dba->dnadb($dna_dba);
      }
      $self->hrdb_set_con($transcript_dba,'transcript_db');

      my ($slice,$accession_array) = $self->get_transcript_region($transcript_id);
      #$query_file = $self->output_query_file($accession_array);
      $querys = $self->get_query_seqs($accession_array);
      if ($self->param('sequence_table_name') =~ /protein/) {
        $self->peptide_seq($querys->[0]->seq);
      }
      $self->calculate_coverage_and_pid($self->param('calculate_coverage_and_pid'));
      @db_files = ($slice);
    } else {
      $self->throw("The feature_type you passed in is not supported! Type:\n".$feature_type);
    }
  } elsif($iid_type eq 'projection_transcript_id') {
       my @iid = @{$self->param("iid")};
       my $transcript_dba = $self->hrdb_get_dba($self->param('transcript_db'));
       my $projection_transcript_id = $iid[0];
       my $projection_protein_accession = $iid[1];
       my $padding = $self->param("projection_padding");

       unless(defined($padding)) {
         $padding = 50000;
       }

       if($dna_dba) {
         $transcript_dba->dnadb($dna_dba);
       }
        $self->hrdb_set_con($transcript_dba,'transcript_db');

       my ($slice,$accession_array) = $self->get_transcript_region($projection_transcript_id);
       $query_seq = $self->get_query_seq($accession_array);
       $self->peptide_seq($query_seq->seq);
       @db_files = ($slice);
   } elsif($iid_type eq 'chunk_file') {
    my $query = $self->QUERYSEQS;

    if(-e $query and -d $query) {
      # query seqs is a directory; input id will be a file in that directory
      # As input_id returns a string, I've made it parse out the file name. I don't
      # like this solution but it is the quickest for the moment
      my $input_id = $self->input_id;

      $query_file = "$query/" . $input_id;
      if (not -e $query_file) {
        $self->throw( "Query file '$query_file' does not exist'\n");
      }
      if ($self->USE_KILL_LIST) {
        $query_file = filter_killed_entries($query_file,$self->KILL_TYPE,$self->REFDB,$self->KILLLISTDB,$self->KILL_LIST_FILTER,$self->input_id);
        $self->filtered_query_file($query_file);
      }
    }
    elsif (-e $query and -s $query) {
      # query seqs is a single file; input id will correspond to a chunk number
      $query_file = $query;
      my $iid_regexp = $self->IIDREGEXP;

      $self->throw("When your input ids are not filenames, you must define ".
                   "IIDREGEXP in config to enable inference of chunk number and total")
      if not defined $iid_regexp;

      ($chunk_number, $chunk_total) = $self->input_id =~ /$iid_regexp/;

      ###
      ### DO THE KILL LIST FILTER FOR QUERY FILE. AGAIN THE FILE CAN CONTAIN MULTIPLE ENTIRES
      ###
      if ($self->USE_KILL_LIST) {
        $query_file = filter_killed_entries($query_file,$self->KILL_TYPE,$self->REFDB,$self->KILLLISTDB,$self->KILL_LIST_FILTER);
      }
    } else {
      $self->throw("'$query' refers to something that could not be made sense of\n");
    }
  } elsif($iid_type eq 'filename') {
    @db_files = ($self->GENOMICSEQS);
    $query_file = $self->input_id;
  } else {
   $self->throw("You provided an input id type that was not recoginised via the 'iid_type' param. Type provided:\n".$iid_type);
  }
  ##########################################
  # Annotation file with CDS positions
  ##########################################


  ##########################################
  # setup the runnables
  ##########################################

  my %parameters = %{$self->parameters_hash};
  if (not exists($parameters{-options}) and
      defined $self->OPTIONS) {
    $parameters{-options} = $self->OPTIONS;
  }
  if (not exists($parameters{-coverage_by_aligned}) and
      defined $self->COVERAGE_BY_ALIGNED) {
    $parameters{-coverage_by_aligned} = $self->COVERAGE_BY_ALIGNED;
  }

  # Old code, will leave here as I could activate it again
  if (defined $self->PROGRAM && defined $self->analysis->program_file) {
    if ($self->PROGRAM ne $self->analysis->program_file) {
      # I'm just warning because for debugging it's easier to change just the PROGRAM parameters...
      $self->warning("CONFLICT: You have defined -program in your config file and ".
            "-program_file in your analysis table.");
    }
  }

  my $biotypes_hash = $self->get_biotype();
  foreach my $database ( @db_files ){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
              -program  => $self->PROGRAM ? $self->PROGRAM : $self->analysis->program_file,
              -analysis => $self->analysis,
              -query_type     => $self->QUERYTYPE,
              -annotation_file => $self->QUERYANNOTATION ? $self->QUERYANNOTATION : undef,
              -query_chunk_number => $chunk_number ? $chunk_number : undef,
              -query_chunk_total => $chunk_total ? $chunk_total : undef,
              -biotypes => $biotypes_hash,
              -calculate_coverage_and_pid => $self->param('calculate_coverage_and_pid'),
              %parameters,
              );

      if (ref($database) eq 'Bio::EnsEMBL::Slice') {
        $runnable->target_seqs([$database]);
      }
      else {
        $runnable->target_file($database);
      }
      if (ref($querys) eq 'ARRAY') {
        $runnable->query_seqs($querys);
      }
      else {
        $runnable->query_file($query_file);
      }
      $runnable->_verbose($self->debug) if ($self->debug);
      $self->runnable($runnable);
  }

}


sub run {
  my ($self) = @_;
  my @results;

  $self->throw("Can't run - no runnable objects") unless ($self->runnable);

  foreach my $runnable (@{$self->runnable}){
    # This is to catch the closing exonerate errors, which we currently have no actual solution for
    # It seems to be more of a problem with the exonerate code itself
    $runnable->run;

    push ( @results, @{$runnable->output} );
  }
  if ($self->USE_KILL_LIST) {
    unlink $self->filtered_query_file;
    # print "Removed temporary query file ".$self->filtered_query_file."\n";
  }
  if ($self->filter) {
    my $filtered_transcripts = $self->filter->filter_results(\@results);
    @results = @$filtered_transcripts;
  }

  my @genes = $self->make_genes(@results);
  $self->param('output_genes',\@genes);
}


sub write_output {
  my ($self) = @_;

  my $outdb = $self->hrdb_get_con('target_db');
  my $gene_adaptor = $outdb->get_GeneAdaptor;

  my @output = @{$self->param('output_genes')};
  $self->param('output_genes',undef);
  my $fails = 0;
  my $total = 0;

  foreach my $gene (@output){
    empty_Gene($gene);
    eval {
      $gene_adaptor->store($gene);
    };
    if ($@){
      $self->warning("Unable to store gene!!\n$@");
      $fails++;
    }
    $total++;
  }
  if ($fails > 0) {
    $self->throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }


}


sub get_transcript_region {
  my ($self,$transcript_id) = @_;

  my $transcript_dba = $self->hrdb_get_con('transcript_db');

  my $transcript = $transcript_dba->get_TranscriptAdaptor()->fetch_by_dbID($transcript_id);
  my $tsf = $transcript->get_all_supporting_features();


  my $feature_pair = ${$tsf}[0];
  my $accession = $feature_pair->hseqname();

  if($self->param('use_genblast_best_in_genome')) {
    my $logic_name = $transcript->analysis->logic_name();
    if($logic_name =~ /_not_best$/) {
      $self->best_in_genome_transcript(0);
    } else {
      $self->best_in_genome_transcript(1);
    }
  }

  my $padding = $self->param('region_padding');
  unless(defined($self->param('region_padding'))) {
    $self->warning("You didn't pass in any value for padding. Defaulting to 10000");
    $padding = 10000;
  }

  my $start = $transcript->seq_region_start;
  my $end = $transcript->seq_region_end;
  my $strand = $transcript->strand;
  my $slice = $transcript->slice();
  my $slice_length = $slice->length();
  if($padding) {
    $start = $start - $padding;
    if($start < 1) {
      $start = 1;
    }
    $end = $end + $padding;
    my $slice = $transcript->slice();
    my $slice_length = $slice->length();
    if($end > $slice_length) {
      $end = $slice_length;
    }
  }

  my @slice_array = split(':',$slice->name());
  $slice_array[3] = $start;
  $slice_array[4] = $end;
  $slice_array[5] = $strand;
  my $new_slice_name = join(':',@slice_array);

  my $sa = $transcript_dba->get_SliceAdaptor();
  my $transcript_slice = $sa->fetch_by_name($new_slice_name);

  return($transcript_slice,[$accession]);
}


sub best_in_genome_transcript {
   my ($self,$val) = @_;

   if(defined($val) && $val==1) {
     $self->param('best_in_genome_transcript', 1);
   } elsif(defined($val) && $val==0) {
     $self->param('best_in_genome_transcript', 0);
   }

   if ($self->param_is_defined('best_in_genome_transcript')) {
     return $self->param('best_in_genome_transcript');
   }
   else {
     return;
   }
}


sub make_genes{
  my ($self,@transcripts) = @_;

  my (@genes);

  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  my %genome_slices;

  foreach my $tran ( @transcripts ){
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->biotype($self->analysis->logic_name);
    $tran->analysis($self->analysis);

    if(defined($self->best_in_genome_transcript()) && $self->best_in_genome_transcript() == 0) {
      my $analysis = $self->analysis;
      my $logic_name = $analysis->logic_name."_not_best";
      $analysis->logic_name($logic_name);
      $gene->analysis($analysis);
      $gene->biotype($logic_name);
      $tran->analysis($analysis);
    }

    ############################################################
    # put a slice on the transcript
    my $slice_id = $tran->start_Exon->seqname;
    if (not exists $genome_slices{$slice_id}) {
      # assumes genome seqs were named in the Ensembl API Slice naming
      # convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
      $genome_slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
    }
    my $slice = $genome_slices{$slice_id};

    foreach my $exon (@{$tran->get_all_Exons}){
      $exon->slice($slice);
      foreach my $evi (@{$exon->get_all_supporting_features}){
        $evi->slice($slice);
        $evi->analysis($self->analysis);
      }
    }
    foreach my $evi (@{$tran->get_all_supporting_features}) {
      $evi->slice($slice);
      $evi->analysis($self->analysis);
    }

    if (!$slice){
        my ($sf);

        if (@{$tran->get_all_supporting_features}) {
          ($sf) = @{$tran->get_all_supporting_features};
        } else {
          my @exons = @{$tran->get_all_Exons};
          ($sf) = @{$exons[0]->get_all_supporting_features};
        }
        print $sf->hseqname."\t$slice_id\n";
    }

    $self->throw("Have no slice") if(!$slice);
    $tran->slice($slice);
    my $accession = $tran->{'accession'};
    my $transcript_biotype = $self->get_biotype->{$accession};
    $tran->biotype($transcript_biotype);

   if($self->calculate_coverage_and_pid) {
      $self->realign_translation($tran);
    }

    $gene->add_Transcript($tran);
    push( @genes, $gene);
  }
  return @genes;
}

############################################################

sub get_chr_names{
  my ($self) = @_;
  my @chr_names;
  my @chromosomes;

  my $chr_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  #also fetching non-reference regions like DR52 for human by default.
  #specify in Exonerate2Genes config-file.
  if(defined($self->NONREF_REGIONS)){
    @chromosomes = @{$chr_adaptor->fetch_all('toplevel', undef, 1)};
  }
  else{
    @chromosomes = @{$chr_adaptor->fetch_all('toplevel')};
  }

  foreach my $chromosome ( @chromosomes ){
    push( @chr_names, $chromosome->seq_region_name );
  }

  return @chr_names;
}

sub get_output_db {
  my ($self) = @_;

  my $outdb;

  if ($self->OUTDB) {
    if ( ref($self->OUTDB)=~m/HASH/) {

      $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{$self->OUTDB},
                                                -dnadb => $self->hrdb_get_con('target_db'));
    }else{
      $outdb = $self->get_dbadaptor($self->OUTDB);
    }
  } else {
    $outdb = $self->hrdb_get_con('target_db');
  }
  return $outdb;
}


sub parse_feature_region_id {
  my ($self,$feature_region_id) = @_;

  my $dba = $self->hrdb_get_con('target_db');
  my $sa = $dba->get_SliceAdaptor();

  unless($feature_region_id =~ s/\:([^\:]+)$//) {
    $self->throw("Could not parse the accession from the feature region id. Expecting a normal slice id, with an extra colon ".
                 "followed by the accession. Offending feature_region_id:\n".$feature_region_id);
  }

  my $slice_name = $feature_region_id;
  my $accession = $1;

  my $slice = $sa->fetch_by_name($slice_name);

  return($slice,[$accession]);

}
############################################################
#
# get/set methods
#
############################################################

sub QUERYSEQS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('QUERYSEQS',$value);
  }

  if ($self->param_is_defined('QUERYSEQS')) {
    return $self->param('QUERYSEQS');
  } else {
    return undef;
  }
}

sub QUERYTYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('QUERYTYPE',$value);
  }

  if ($self->param_is_defined('QUERYTYPE')) {
    return $self->param('QUERYTYPE');
  } else {
    return undef;
  }
}

sub QUERYANNOTATION {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('QUERYANNOTATION',$value);
  }

  if ($self->param_is_defined('QUERYANNOTATION')) {
    return $self->param('QUERYANNOTATION');
  } else {
    return undef;
  }
}

sub GENOMICSEQS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('GENOMICSEQS',$value);
  }

  if ($self->param_is_defined('GENOMICSEQS')) {
    return $self->param('GENOMICSEQS');
  } else {
    return undef;
  }
}

sub IIDREGEXP {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('IIDREGEXP',$value);
  }

  if ($self->param_is_defined('IIDREGEXP')) {
    return $self->param('IIDREGEXP');
  } else {
    return undef;
  }
}

sub REFDB {
  my ($self,$value) = @_;

  if (defined $value) {
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                  %$value
                                                );
    $self->param('REFDB',$dba);
    # Set this to override the default dbc which is inherited from Process and is to the Hive db
    #$self->db($dba);
  }

  if ($self->param_is_defined('REFDB')) {
    return $self->param('REFDB');

  } else {
    return undef;
  }
}

sub OUTDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('OUTDB',$value);
  }

  if ($self->param_is_defined('OUTDB')) {
    return $self->param('OUTDB');

  } else {
    return undef;
  }
}

sub KILLLISTDB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('KILLLISTDB',$value);
  }

  if ($self->param_is_defined('KILLLISTDB')) {
    return $self->param('KILLLISTDB');
  } else {
    return undef;
  }
}

sub COVERAGE_BY_ALIGNED {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('COVERAGE_BY_ALIGNED',$value);
  }

  if ($self->param_is_defined('COVERAGE_BY_ALIGNED')) {
    return $self->param('COVERAGE_BY_ALIGNED');
  } else {
    return undef;
  }
}


sub OPTIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('OPTIONS',$value);
  }

  if ($self->param_is_defined('OPTIONS')) {
    return $self->param('OPTIONS');
  } else {
    return undef;
  }
}

sub NONREF_REGIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('NONREF_REGIONS',$value);
  }

  if ($self->param_is_defined('NONREF_REGIONS')) {
    return $self->param('NONREF_REGIONS');
  } else {
    return undef;
  }
}

sub PROGRAM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('PROGRAM',$value);
  }

  if ($self->param_is_defined('PROGRAM')) {
    return $self->param('PROGRAM');
  } else {
    return undef;
  }
}

sub USE_KILL_LIST {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('USE_KILL_LIST',$value);
  }

  if ($self->param_is_defined('USE_KILL_LIST')) {
    return $self->param('USE_KILL_LIST');
  } else {
    return undef;
  }
}

sub KILL_TYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('KILL_TYPE',$value);
  }

  if ($self->param_is_defined('KILL_TYPE')) {
    return $self->param('KILL_TYPE');
  } else {
    return undef;
  }
}

sub KILL_LIST_FILTER {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('KILL_LIST_FILTER',$value);
  }

  if ($self->param_is_defined('KILL_LIST_FILTER')) {
    return $self->param('KILL_LIST_FILTER');
  } else {
    return undef;
  }
}

sub SOFT_MASKED_REPEATS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('SOFT_MASKED_REPEATS',$value);
  }

  if ($self->param_is_defined('SOFT_MASKED_REPEATS')) {
    return $self->param('SOFT_MASKED_REPEATS');
  } else {
    return undef;
  }
}

sub SEQFETCHER_PARAMS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('SEQFETCHER_PARAMS',$value);
  }

  if ($self->param_is_defined('SEQFETCHER_PARAMS')) {
    return $self->param('SEQFETCHER_PARAMS');
  } else {
    return undef;
  }
}

sub SEQFETCHER_OBJECT {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('SEQFETCHER_OBJECT',$value);
  }

  if ($self->param_is_defined('SEQFETCHER_OBJECT')) {
    return $self->param('SEQFETCHER_OBJECT');
  } else {
    return undef;
  }
}

sub query_seqs {
  my ($self, @seqs) = @_;
  if( @seqs ) {
    unless ($seqs[0]->isa("Bio::PrimarySeqI") || $seqs[0]->isa("Bio::SeqI")){
      $self->throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    push( @{$self->param('_query_seqs')}, @seqs);
  }
  return @{$self->param('_query_seqs')};
}

############################################################

sub genomic {
  my ($self, $seq) = @_;
  if ($seq){
    unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")){
      $self->throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->param('_genomic',$seq);
  }
  return $self->param('_genomic');
}


############################################################

sub database {
  my ($self, $database) = @_;
  if ($database) {
    $self->param('_database',$database) = $database;
  }
  return $self->param('_database');
}

############################################################

sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->param('_transcript_filter',$val);
  }
  elsif (!$self->param_is_defined('_transcript_filter')
    and $self->param_is_defined('FILTER')
    and exists $self->param('FILTER')->{OBJECT}) {
    my $module = $self->require_module($self->param('FILTER')->{OBJECT});
    $self->param('_transcript_filter', $module->new(%{$self->param('FILTER')->{PARAMETERS}}));
  }
  if ($self->param_is_defined('_transcript_filter')) {
    return $self->param('_transcript_filter');
  }
  else {
    return;
  }
}

############################################################

sub filtered_query_file {
  my ($self, $val) = @_;
  if ($val) { 
    $self->param('_filtered_query_file',$val);
  }
  return $self->param('_filtered_query_file');
}


###############################################
###     end of config
###############################################

sub filter_killed_entries {
  my ($orig_query_filename, $mol_type,$input_db,$killlist_db,$filter_params,$inputID) = @_;
  my $kill_list_object = Bio::EnsEMBL::KillList::KillList
      ->new(-TYPE => $mol_type, -GB_REF_DB => $input_db, -KILL_LIST_DB => $killlist_db, -FILTER_PARAMS => $filter_params);
  my %kill_list = %{ $kill_list_object->get_kill_list() };

  my $seqin  = new Bio::SeqIO(-file   => "<$orig_query_filename",
                            -format => "Fasta",
                          );


  my @sequences;
  while( my $query_entry = $seqin->next_seq ){
    my $display_id  = $query_entry->display_id;
    my $no_ver_id;
    # Depending on the display ID's format, strip off the
    # version number because the kill_list hash keys are
    # without version numbers

    if ($display_id =~/\w+\.\d/) {
      ($no_ver_id) = $display_id =~/(\w+)\.\d/;
    } elsif ($display_id =~/\w+\-\d/) {
      ($no_ver_id) = $display_id =~/(\w+)\-\d/;
    } elsif ($display_id =~/\w+/ ) {
      ($no_ver_id) = $display_id;
    }

    if ( !$kill_list{$no_ver_id} ) {
      push(@sequences, $query_entry);
    } elsif ( $kill_list{$no_ver_id} ) {
      print "$mol_type $display_id is in the kill_list. Discarded from analysis.\n";
    }
  }
  return \@sequences;
}

sub get_query_seqs {
  my ($self, $accession_array) = @_;

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));

  my $biotypes_hash = {};
  my @query_sequences;
  foreach my $accession (@{$accession_array}) {
    my $db_row = $table_adaptor->fetch_by_dbID($accession);
    unless($db_row) {
      $self->throw('Did not find an entry in the '.$self->param('sequence_table_name')." table matching the accession. Accession:\n".$accession);
    }

    my $seq = $db_row->{'seq'};
    $biotypes_hash->{$accession} = $db_row->{'biotype'};

    push(@query_sequences, Bio::Seq->new(-display_id => $accession, -seq => $seq));
  }

  $self->get_biotype($biotypes_hash);

  return \@query_sequences;
}


sub get_biotype {
  my ($self,$biotype_hash) = @_;
  if($biotype_hash) {
    $self->param('_biotype_hash',$biotype_hash);
  }
  return($self->param('_biotype_hash'));
}


sub peptide_seq {
  my ($self, $value) = @_;
  if($value){
    $self->param('_peptide_seq', $value);
  }
  return $self->param('_peptide_seq');
}


sub calculate_coverage_and_pid {
  my ($self, $value) = @_;
  if($value){
    $self->param('calculate_coverage_and_pid', $value);
  }
  return $self->param('calculate_coverage_and_pid');
}


sub realign_translation {
  my ($self,$transcript) = @_;

  my $query_seq = $self->peptide_seq;
  my $translation = $transcript->translate->seq();

  my $align_input_file = "/tmp/exonerate_align_".$$.".fa";
  my $align_output_file = "/tmp/exonerate_align_".$$.".aln";

  open(INPUT,">".$align_input_file);
  say INPUT ">query";
  say INPUT $query_seq;
  say INPUT ">target";
  say INPUT $translation;
  close INPUT;

  my $align_program_path = 'muscle';

  my $cmd = $align_program_path." -in ".$align_input_file." -out ".$align_output_file;
  my $result = system($cmd);

  if($result) {
    $self->throw("Got a non-zero exit code from alignment. Commandline used:\n".$cmd);
  }

  my $file = "";
  open(ALIGN,$align_output_file);
  while(<ALIGN>) {
    $file .= $_;
  }
  close ALIGN;

  unless($file =~ /\>.+\n(([^>]+\n)+)\>.+\n(([^>]+\n)+)/) {
    $self->throw("Could not parse the alignment file for the alignment sequences. Alignment file: ".$align_output_file);
  }

  my $aligned_query_seq = $1;
  my $aligned_target_seq = $3;

  $aligned_query_seq =~ s/\n//g;
  $aligned_target_seq =~ s/\n//g;

  say "Aligned query:\n".$aligned_query_seq;
  say "Aligned target:\n".$aligned_target_seq;

  `rm $align_input_file`;
  `rm $align_output_file`;

  # Work out coverage
  my $coverage;
  my $temp = $aligned_target_seq;
  my $target_gap_count = $temp =~ s/\-//g;
  my $ungapped_query_seq = $aligned_query_seq;
  $ungapped_query_seq  =~ s/\-//g;

  if(length($ungapped_query_seq) == 0) {
    $coverage = 0;
  } else {
    $coverage = 100 - (($target_gap_count/length($ungapped_query_seq)) * 100);
  }

  # Work out percent identity
  my $match_count = 0;
  my $aligned_positions = 0;
  for(my $j=0; $j<length($aligned_query_seq); $j++) {
    my $char_query = substr($aligned_query_seq,$j,1);
    my $char_target = substr($aligned_target_seq,$j,1);
    if($char_query eq '-' || $char_target  eq '-') {
      next;
    }
    if($char_query eq $char_target) {
      $match_count++;
    }
    $aligned_positions++;
  }

  unless($aligned_positions) {
    $self->throw("Pairwise alignment between the query sequence and the translation shows zero aligned positions. Something has gone wrong");
  }

  my $percent_id = ($match_count / $aligned_positions) * 100;

  # Get all exons and transcript supporting features
  my $transcript_supporting_features = $transcript->get_all_supporting_features();
  my $exons = $transcript->get_all_Exons();

  # Now clean these out
  $transcript->flush_Exons();
  $transcript->flush_supporting_features();

  # Loop through the TSFs and add the coverage and pid, then add back into transcript
  foreach my $transcript_supporting_feature (@{$transcript_supporting_features}) {
    $transcript_supporting_feature->hcoverage($coverage);
    $transcript_supporting_feature->percent_id($percent_id);
    $transcript->add_supporting_features($transcript_supporting_feature);
  }

  # Loop through exons, get supporting features for each, flush existing SF, add coverage and pid, add back to exon, add exon to transcript
  foreach my $exon (@{$exons}) {
    my $exon_supporting_features = $exon->get_all_supporting_features();
    $exon->flush_supporting_features();
    foreach my $exon_supporting_feature (@{$exon_supporting_features}) {
      $exon_supporting_feature->hcoverage($coverage);
      $exon_supporting_feature->percent_id($percent_id);
      $exon->add_supporting_features($exon_supporting_feature);
    }
    $transcript->add_Exon($exon);
  }
}

sub output_db_file {
  my ($self) = @_;

  $self->throw('Method '.ref($self).'::output_db_file has not been implemeted. Only the Runnable should write temporary files. You should pass the databases as Bio::Seq');
}

1;
