=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::KillList::KillList;
use Bio::SeqIO;
use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


############################################################
sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($EXONERATE_CONFIG_BY_LOGIC);

  return $self;
}


sub fetch_input {
  my( $self) = @_;

  my $logic = $self->analysis->logic_name;

  ##########################################
  # set up the target (genome)
  ##########################################

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
      
      if ($dir){  
  
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
	        warning( "Could not find fasta file for '$chr_name' in directories:\n".
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
    else{
      $target =~s/^\s+//;  
      if (-e $target and -d $target) {
        # genome is in a directory; the directory must contain the complete
        # genome else we cannot do best-in-genome filtering. 
        # 
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
            warning( "Could not find fasta file for '$chr' in '$target'\n");
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
        throw("'$target' refers to something that could not be made sense of");
      }
    }
  }

  ##########################################
  # set up the query (est/cDNA/protein)
  ##########################################

  my ($query_file, $chunk_number, $chunk_total);

  my $query = $self->QUERYSEQS;

  if (-e $query and -d $query) {
    # query seqs is a directory; input id will be a file in that directory
    $query_file = "$query/" . $self->input_id;
    if (not -e $query_file) {
      throw( "Query file '$query_file' does not exist'\n");
    }
    if ($self->USE_KILL_LIST) {
      $query_file = filter_killed_entries($query_file, $self->KILL_TYPE, $self->input_id);
      $self->filtered_query_file($query_file);
    }  
  }
  elsif (-e $query and -s $query) {
    # query seqs is a single file; input id will correspond to a chunk number
    $query_file = $query;
    my $iid_regexp = $self->IIDREGEXP;
    
    throw("When your input ids are not filenames, you must define ".
          "IIDREGEXP in config to enable inference of chunk number and total")
        if not defined $iid_regexp;

    ($chunk_number, $chunk_total) = $self->input_id =~ /$iid_regexp/;

    ###
    ### DO THE KILL LIST FILTER FOR QUERY FILE. AGAIN THE FILE CAN CONTAIN MULTIPLE ENTIRES
    ###
    if ($self->USE_KILL_LIST) {
      $query_file = filter_killed_entries($query_file, $self->KILL_TYPE);
    }
  } else {
    throw("'$query' refers to something that could not be made sense of\n");
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
    $parameters{-options} = $self->OPTIONS
  }
  if (not exists($parameters{-coverage_by_aligned}) and
      defined $self->COVERAGE_BY_ALIGNED) {
    $parameters{-coverage_by_aligned} = $self->COVERAGE_BY_ALIGNED;
  }

  if (defined $self->PROGRAM && defined $self->analysis->program_file) {
    if ($self->PROGRAM ne $self->analysis->program_file) {
# I'm just warning because for debugging it's easier to change just the PROGRAM parameters...
      warning("CONFLICT: You have defined -program in your config file and ".
            "-program_file in your analysis table.");
    }
  }

  foreach my $database ( @db_files ){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
        ->new(
              -program  => $self->PROGRAM ? $self->PROGRAM : $self->analysis->program_file,
              -analysis => $self->analysis,
              -target_file    => $database,
              -query_type     => $self->QUERYTYPE,
              -query_file     => $query_file,
              -annotation_file => $self->QUERYANNOTATION ? $self->QUERYANNOTATION : undef,
              -query_chunk_number => $chunk_number ? $chunk_number : undef,
              -query_chunk_total => $chunk_total ? $chunk_total : undef,
              %parameters,
              );
    $self->runnable($runnable);
  }

}

############################################################

sub run{
  my ($self) = @_;
  my @results;
  
  throw("Can't run - no runnable objects") unless ($self->runnable);
  
  foreach my $runnable (@{$self->runnable}){
    
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
  
  $self->output(\@genes);
}


############################################################

sub write_output{
  my ($self,@output) = @_;

  my $outdb = $self->get_output_db;
  my $gene_adaptor = $outdb->get_GeneAdaptor;  

  unless (@output){
    @output = @{$self->output};
  }
  
  my $fails = 0;
  my $total = 0;
  foreach my $gene (@output){

    eval {
      $gene_adaptor->store($gene);
    };    
    if ($@){
      warning("Unable to store gene!!\n$@");
      $fails++;
    }
    $total++;
  }
  if ($fails > 0) {
    throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }
}

############################################################

sub make_genes{
  my ($self,@transcripts) = @_;
  
  my (@genes);

  my $slice_adaptor = $self->db->get_SliceAdaptor;
  
  my %genome_slices;

  foreach my $tran ( @transcripts ){
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
 
    throw("Have no slice") if(!$slice);
    $tran->slice($slice);
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

  my $chr_adaptor = $self->db->get_SliceAdaptor;
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


############################################################

sub get_output_db {
  my ($self) = @_;

  my $outdb;

  if ($self->OUTDB) {
    if ( ref($self->OUTDB)=~m/HASH/) {

      $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{$self->OUTDB},
                                                -dnadb => $self->db);
    }else{
      $outdb = $self->get_dbadaptor($self->OUTDB);
    }
  } else {
    $outdb = $self->db;
  }
  $self->db->dbc->disconnect_when_inactive(1) ;
  $outdb->dbc->disconnect_when_inactive(1) ;
  return $outdb;
}


############################################################
#
# get/set methods
#
############################################################


sub query_seqs {
  my ($self, @seqs) = @_;
  if( @seqs ) {
    unless ($seqs[0]->isa("Bio::PrimarySeqI") || $seqs[0]->isa("Bio::SeqI")){
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    push( @{$self->{_query_seqs}}, @seqs);
  }
  return @{$self->{_query_seqs}};
}

############################################################

sub genomic {
  my ($self, $seq) = @_;
  if ($seq){
    unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")){
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_genomic} = $seq ;
  }
  return $self->{_genomic};
}


############################################################

sub database {
  my ($self, $database) = @_;
  if ($database) {
    $self->{_database} = $database;
  }
  return $self->{_database};
}

############################################################

sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->{_transcript_filter} = $val;
  }
  return $self->{_transcript_filter};
}

############################################################

sub filtered_query_file {
  my ($self, $val) = @_;
  if ($val) { 
    $self->{_filtered_query_file} = $val;
  }
  return $self->{_filtered_query_file};
}

#############################################################
# Declare and set up config variables
#############################################################

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($EXONERATE_CONFIG_BY_LOGIC);

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
  foreach my $config_var (qw(QUERYSEQS 
                             QUERYTYPE
                             GENOMICSEQS)) {

   throw("You must define $config_var in config for logic '$logic'")
        if not defined $self->$config_var;
  }
  
  throw("QUERYANNOTATION '" . $self->QUERYANNOTATION . "' in config must be readable")
      if $self->QUERYANNOTATION and not -e $self->QUERYANNOTATION;

  # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
  if ($self->FILTER) {
    if (not ref($self->FILTER) eq "HASH" or
        not exists($self->FILTER->{OBJECT}) or
        not exists($self->FILTER->{PARAMETERS})) {
          
      throw("FILTER in config fo '$logic' must be a hash ref with elements:\n" . 
            "  OBJECT : qualified name of the filter module;\n" .
            "  PARAMETERS : anonymous hash of parameters to pass to the filter");
    } else {
      my $module = $self->FILTER->{OBJECT};
      my $pars   = $self->FILTER->{PARAMETERS};
      
      (my $class = $module) =~ s/::/\//g;
      eval{
        require "$class.pm";
      };
      throw("Couldn't require ".$class." Exonerate2Genes:require_module $@") if($@);
    
      $self->filter($module->new(%{$pars}));
    }
  }
}

sub QUERYSEQS {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_QUERYSEQS'} = $value;
  }

  if (exists($self->{'_CONFIG_QUERYSEQS'})) {
    return $self->{'_CONFIG_QUERYSEQS'};
  } else {
    return undef;
  }
}

sub QUERYTYPE {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_QUERYTYPE'} = $value;
  }

  if (exists($self->{'_CONFIG_QUERYTYPE'})) {
    return $self->{'_CONFIG_QUERYTYPE'};
  } else {
    return undef;
  }
}


sub QUERYANNOTATION {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_QUERYANNOTATION'} = $value;
  }

  if (exists($self->{'_CONFIG_QUERYANNOTATION'})) {
    return $self->{'_CONFIG_QUERYANNOTATION'};
  } else {
    return undef;
  }
}



sub GENOMICSEQS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_GENOMICSEQS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_GENOMICSEQS'})) {
    return $self->{'_CONFIG_GENOMICSEQS'};
  } else {
    return undef;
  }
}


sub IIDREGEXP {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_IIDREGEXP'} = $value;
  }

  if (exists($self->{'_CONFIG_IIDREGEXP'})) {
    return $self->{'_CONFIG_IIDREGEXP'};
  } else {
    return undef;
  }
}

sub OUTDB {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_OUTDB'} = $value;
  }

  if (exists($self->{'_CONFIG_OUTDB'})) {
    return $self->{'_CONFIG_OUTDB'};

  } else {
    return undef;
  }
}

sub COVERAGE_BY_ALIGNED {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_COVERAGE'} = $value;
  }

  if (exists($self->{'_CONFIG_COVERAGE'})) {
    return $self->{'_CONFIG_COVERAGE'};
  } else {
    return undef;
  }
}


sub FILTER {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_FILTER'} = $value;
  }
  
  if (exists($self->{'_CONFIG_FILTER'})) {
    return $self->{'_CONFIG_FILTER'};
  } else {
    return undef;
  }
}

sub OPTIONS {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{'_CONFIG_OPTIONS'} = $value;
  }

  if (exists($self->{'_CONFIG_OPTIONS'})) {
    return $self->{'_CONFIG_OPTIONS'};
  } else {
    return undef;
  }
}

sub NONREF_REGIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_NONREF_REGIONS'} = $value;
  }

  if (exists($self->{'_CONFIG_NONREF_REGIONS'})) {
    return $self->{'_CONFIG_NONREF_REGIONS'};
  } else {
    return undef;
  }
}

sub PROGRAM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_PROGRAM'} = $value;
  }

  if (exists($self->{'_CONFIG_PROGRAM'})) {
    return $self->{'_CONFIG_PROGRAM'};
  } else {
    return undef;
  }
}

sub USE_KILL_LIST {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_USE_KILL_LIST'} = $value;
  }

  if (exists($self->{'_CONFIG_USE_KILL_LIST'})) {
    return $self->{'_CONFIG_USE_KILL_LIST'};
  } else {
    return undef;
  }
}

sub KILL_TYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_KILL_TYPE'} = $value;
  }

  if (exists($self->{'_CONFIG_KILL_TYPE'})) {
    return $self->{'_CONFIG_KILL_TYPE'};
  } else {
    return undef;
  }
}

sub SOFT_MASKED_REPEATS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_SOFT_MASKED_REPEATS'} = $value;
  }

  if (exists($self->{'_SOFT_MASKED_REPEATS'})) {
    return $self->{'_SOFT_MASKED_REPEATS'};
  } else {
    return undef;
  }
}

sub SEQFETCHER_PARAMS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_SEQFETCHER_PARAMS'} = $value;
  }

  if (exists($self->{'_SEQFETCHER_PARAMS'})) {
    return $self->{'_SEQFETCHER_PARAMS'};
  } else {
    return undef;
  }
}

sub SEQFETCHER_OBJECT {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_SEQFETCHER_OBJECT'} = $value;
  }

  if (exists($self->{'_SEQFETCHER_OBJECT'})) { 
    return $self->{'_SEQFETCHER_OBJECT'};
  } else {
    return undef;
  }
}


###############################################
###     end of config
###############################################

sub filter_killed_entries {
  my ($orig_query_filename, $mol_type, $inputID) = @_;
  my $kill_list_object = Bio::EnsEMBL::KillList::KillList
      ->new(-TYPE => $mol_type);
  my %kill_list = %{ $kill_list_object->get_kill_list() };

  my $seqin  = new Bio::SeqIO(-file   => "<$orig_query_filename",
                            -format => "Fasta",
                          );

  my $filtered_seqout_filename = "/tmp/$inputID"."_filtered";
  print "Filename for my filtered sequence: $filtered_seqout_filename.\n";

  my $seqout = new Bio::SeqIO(-file   => ">$filtered_seqout_filename",
                              -format => "Fasta"
                           );

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
      $seqout->write_seq($query_entry);
    } elsif ( $kill_list{$no_ver_id} ) {
      print "$mol_type $display_id is in the kill_list. Discarded from analysis.\n";
    }
  }
  return $filtered_seqout_filename; 
}    

1;
