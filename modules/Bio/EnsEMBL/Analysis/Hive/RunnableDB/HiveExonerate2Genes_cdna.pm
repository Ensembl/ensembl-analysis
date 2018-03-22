#!/usr/bin/env perl

=head1 LICENSE

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna;

use warnings ;
use strict;
use feature 'say';

use Bio::EnsEMBL::Hive::DBSQL::DataflowRuleAdaptor;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::UnmappedObject;
use Bio::SeqIO;
use Bio::Seq;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    unaligned_to_unmapped => 0,
  }
}


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
  my @target_list = $self->GENOMICSEQS;


  foreach my $target (@target_list){ 

    if(ref $target eq 'ARRAY'){
      #check to see if we have multiple files or directories:

      my $dir = 0;
      foreach my $alt_target (@$target){
        if (-d $alt_target){
            $dir = 1;
            last;
        }
      }

      # genome is in multiple directories;  the order of directories determines
      # which file is used in case of duplicates. New versions should therefore
      # be in the directory listed first.

      if ($dir) {

        foreach my $chr_name ($self->get_chr_names) {
          my $found = 0;
          DIRCHECK:
          foreach my $alt_target (@$target){
            if (-s "$alt_target/$chr_name.fa") {
              push @db_files, "$alt_target/$chr_name.fa";
              $found = 1;
              last DIRCHECK;
            }
          }
          if(!$found){
	        $self->warning( "Could not find fasta file for '$chr_name' in directories:\n".
            join("\n\t", @$target)."\n");
          }
        }
      }else{
        foreach my $alt_target (@$target){
          if (-s $alt_target){
             push @db_files, $alt_target;
          }
        }
      }
    } # // end target is a directory
    else {
      $target =~s/^\s+//;
      if (-e $target and -d $target) {
        # genome is in a directory; the directory must contain the complete
        # genome else we cannot do best-in-genome filtering.
        # We would like to use exonerate's ability to accept a directory as
        # target (because bestn then works), but we must check that the directory
        # contains only toplevel sequence files

        my %dir_contents;
        opendir DIR, $target;
        while(my $entry = readdir DIR) {
          if ($entry ne '.' and $entry ne '..') {
            $dir_contents{$entry} = 0;
          }
        }
        closedir(DIR);

        foreach my $chr ($self->get_chr_names) {
          my $seq_fname = "$chr.fa"; 
          if (-s "$target/$seq_fname") {
            $dir_contents{$seq_fname}++;
            push @db_files, "$target/$seq_fname";
          } else {
            $self->warning( "Could not find fasta file for '$chr' in '$target'\n");
          }
        }

        # if all files in dir were expected, we can revert to having
        # the whole directory as target
        if (not grep { $dir_contents{$_} == 0 } keys %dir_contents) {
          @db_files = ($target);
        }
      }
      elsif (-e $target and -s $target) {
        # genome sequence is in a single file
        @db_files = ($target);
      } else {
        $self->throw("'$target' refers to something that could not be made sense of");
      }
    }
  }

  ##########################################
  # set up the query (est/cDNA/protein)
  ##########################################
  my $iid_type = $self->param('iid_type');
  my $timer = $self->param('timer');
  #$timer = parse_timer($timer);
  my ($query_file,$chunk_number,$chunk_total);
  unless($iid_type) {
    $self->throw("You haven't provided an input id type. Need to provide one via the 'iid_type' param");
  }

  if($iid_type eq 'db_seq') {
    $query_file = $self->output_query_file();
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

#  my $transcript_biotype = $self->transcript_biotype();
  my $biotypes_hash = $self->get_biotype();
  foreach my $database ( @db_files ){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
              -program  => $self->PROGRAM ? $self->PROGRAM : $self->analysis->program_file,
              -analysis => $self->analysis,
              -target_file    => $database,
              -query_type     => $self->QUERYTYPE,
              -annotation_file => $self->QUERYANNOTATION ? $self->QUERYANNOTATION : undef,
              -query_chunk_number => $chunk_number ? $chunk_number : undef,
              -query_chunk_total => $chunk_total ? $chunk_total : undef,
              -biotypes => $biotypes_hash,
              -timer => $timer,
               %parameters,
              );

      if (ref($query_file) eq 'ARRAY') {
        $runnable->query_seqs($query_file);
      }
      else {
        $runnable->query_file($query_file);
      }
      $runnable->_verbose($self->debug);
      $self->runnable($runnable);
  }

  $self->throw("Can't run - no runnable objects") unless ($self->runnable);

}


sub run {
  my ($self) = @_;
  my @results;

  foreach my $runnable (@{$self->runnable}){
    eval {
      $runnable->run;
    };

    if($@) {
      my $except = $@;
      if($except =~ /still running after your timer/) {
        $self->dataflow_output_id({iid => $self->param('iid')}, Bio::EnsEMBL::Hive::DBSQL::DataflowRuleAdaptor::branch_name_2_code('RUNLIMIT'));
        $self->complete_early("Exonerate took longer than the timer limit (".$self->param('timer')."), will dataflow input id on branch 'RUN_LIMIT'. Exception:\n".$except);
      } else {
        $self->throw("Issue with running Exonerate:\n".$except);
      }
    } else {
      push ( @results, @{$runnable->output} );
    }
  }
  if ($self->USE_KILL_LIST) {
    unlink $self->filtered_query_file;
  }
  if ($self->filter and @results) {
    my ($filtered_transcripts, $unmapped) = $self->filter->filter_results(\@results);
    @results = @$filtered_transcripts;
    if (@$unmapped) {
      my @unmapped_objects;
      foreach my $unmapped_object (@$unmapped) {
        my $id = $unmapped_object->identifier;
        push(@unmapped_objects, $unmapped_object) unless (grep {$_->{accession} eq $id} @results);
      }
      $self->param('unmapped_objects', \@unmapped_objects) if (@unmapped_objects);
    }
  }

  $self->output($self->make_genes(\@results));
}


sub write_output {
  my ($self) = @_;

  my $outdb = $self->hrdb_get_con('target_db');
  my $gene_adaptor = $outdb->get_GeneAdaptor;
  my $analysis = $self->analysis;
  my $output = $self->output;
  my $fails = 0;
  my $total = 0;

  if (scalar(@$output) == 0){
    my $iids;
    if ($self->param_is_defined('unmapped_objects')) {
      my $unmapped = $self->param('unmapped_objects');
      foreach my $iid (@{$self->param('iid')}) {
        push(@$iids, $iid) unless (grep {$_->identifier eq $iid} @$unmapped);
      }
    }
    else {
      $iids = $self->param('iid');
    }
    $self->dataflow_output_id({iid => $iids}, $self->param('_branch_to_flow_to')) if ($iids);
  } else {
    foreach my $gene (@$output){
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
    if ($fails){
      $self->throw("Not all genes could be written successfully " .
        "($fails fails out of $total)");
    }
  }
  if ($self->param_is_defined('unmapped_objects')) {
    my $unmapped_adaptor = $outdb->get_UnmappedObjectAdaptor;
    foreach my $unmapped_object (@{$self->param('unmapped_objects')}) {
      $unmapped_adaptor->store($unmapped_object);
    }
  }
}


sub make_genes{
  my ($self, $transcripts) = @_;

  my (@genes);

  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  my %genome_slices;

  foreach my $tran ( @$transcripts ){
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->biotype($self->analysis->logic_name);

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
    $tran->analysis($self->analysis);
    my $accession = $tran->{'accession'};
    my $transcript_biotype = $self->get_biotype->{$accession};
    $tran->biotype($transcript_biotype);
    $gene->add_Transcript($tran);
    push( @genes, $gene);
  }
  return \@genes;
}

############################################################

sub runnable_failed {
  my ($self,$runnable_failed) = @_;
  if (defined $runnable_failed) {
    $self->param('_runnable_failed',$runnable_failed);
  }
  return ($self->param('_runnable_failed'));
}

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
  $self->hrdb_get_con('target_db')->dbc->disconnect_when_inactive(1) ;
  $outdb->dbc->disconnect_when_inactive(1) ;
  return $outdb;
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
    $self->param('COVERAGE',$value);
  }

  if ($self->param_is_defined('COVERAGE')) {
    return $self->param('COVERAGE');
  } else {
    return undef;
  }
}

sub FILTER {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('FILTER',$value);
  }

  if ($self->param_is_defined('FILTER')) {
    return $self->param('FILTER');
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
    $self->param('_SOFT_MASKED_REPEATS',$value);
  }

  if ($self->param_is_defined('_SOFT_MASKED_REPEATS')) {
    return $self->param('_SOFT_MASKED_REPEATS');
  } else {
    return undef;
  }
}

sub SEQFETCHER_PARAMS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('_SEQFETCHER_PARAMS',$value);
  }

  if ($self->param_is_defined('_SEQFETCHER_PARAMS')) {
    return $self->param('_SEQFETCHER_PARAMS');
  } else {
    return undef;
  }
}

sub SEQFETCHER_OBJECT {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('_SEQFETCHER_OBJECT',$value);
  }

  if ($self->param_is_defined('_SEQFETCHER_OBJECT')) {
    return $self->param('_SEQFETCHER_OBJECT');
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

#############################################################
# Declare and set up config variables
#############################################################

#sub read_and_check_config {
#  my $self = shift;

#  $self->SUPER::read_and_check_config($EXONERATE_CONFIG_BY_LOGIC);

  ##########
  # CHECKS
  ##########
#  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
#  foreach my $config_var (qw(QUERYSEQS 
#                             QUERYTYPE
#                             GENOMICSEQS)) {

#   throw("You must define $config_var in config for logic '$logic'")
#        if not defined $self->$config_var;
#  }
  
#  throw("QUERYANNOTATION '" . $self->QUERYANNOTATION . "' in config must be readable")
#      if $self->QUERYANNOTATION and not -e $self->QUERYANNOTATION;

  # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
#  if ($self->FILTER) {
#    if (not ref($self->FILTER) eq "HASH" or
#        not exists($self->FILTER->{OBJECT}) or
#        not exists($self->FILTER->{PARAMETERS})) {

#      throw("FILTER in config fo '$logic' must be a hash ref with elements:\n" . 
#            "  OBJECT : qualified name of the filter module;\n" .
#            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
#    } else {
#      my $module = $self->FILTER->{OBJECT};
#      my $pars   = $self->FILTER->{PARAMETERS};
      
#      (my $class = $module) =~ s/::/\//g;
#      eval{
#        require "$class.pm";
#      };
#      throw("Couldn't require ".$class." Exonerate2Genes:require_module $@") if($@);
#    
#      $self->filter($module->new(%{$pars}));
#    }
#  }
#}



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

sub output_query_file {
  my ($self) = @_;

  my $accession_array = $self->param('iid');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();

  # table name here should probably be changed to something more general
  $table_adaptor->table_name('cdna_sequences');


#  my $output_dir = $self->param('query_seq_dir');
#
#  # Note as each accession will occur in only one file, there should be no problem using the first one
#  my $outfile_name = "exonerate_".${$accession_array}[0].".fasta";
#  my $outfile_path = $output_dir."/".$outfile_name;

  my $biotypes_hash = {};

#  unless(-e $output_dir) {
#    `mkdir $output_dir`;
#  }
#
#  if(-e $outfile_path) {
#    $self->warning("Found the query file in the query dir already. Overwriting. File path:\n".$outfile_path);
#  }

  my @query_sequences;
  foreach my $accession (@{$accession_array}) {
    my $db_row = $table_adaptor->fetch_by_dbID($accession);
    unless($db_row) {
      $self->throw('Did not find an entry in the '.$self->param('query_table_name')." table matching the accession. Accession:\n".$accession);
    }

    my $seq = $db_row->{'seq'};
    $biotypes_hash->{$accession} = $db_row->{'biotype'};

    push(@query_sequences, Bio::Seq->new(-id => $accession, -seq => $seq));
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




1;
