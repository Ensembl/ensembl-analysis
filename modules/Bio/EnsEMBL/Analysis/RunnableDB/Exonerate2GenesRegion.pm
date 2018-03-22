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

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2GenesRegion - 

=head1 SYNOPSIS

my $exonerate2genes = Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2GenesRegion->new(
                              -db         => $refdb,
			      -analysis   => $analysis_obj,
			      -input_id => $slice_and_query_seq_id
			     );

$exonerate2genes->fetch_input();
$exonerate2genes->run();
$exonerate2genes->output();
$exonerate2genes->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
It is meant to provide the interface for mapping nucleotide or protein sequences to 
regions of the genome
sequence and writing the results as genes. By the way Exonerate is run
we do not cluster transcripts into genes and only write one transcript per gene.
we then create a dbadaptor for the target database.
It expects an input id consisting of an ensembl slice name followed by :: and then either
an accession/id that can be located in the specified index or the name of a file located 
in the QUERYSEQS dir.
Example input ID: chromosome:GRCh37:10:100005443:100030007:1::NM_032211.6
Reads config from Config/Exonerate2Genes.pm


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2GenesRegion;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (parse_config create_file_name write_seqfile);
use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;
use Bio::EnsEMBL::Analysis::Config::General qw ( ANALYSIS_WORK_DIR ) ;
#Needed for subroutine &Bio::DB::Flat::OBDAIndex::FileHandle
#called at bioperl-live/Bio/DB/Flat/OBDAIndex.pm
use FileHandle;
use Bio::SeqIO;
use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


############################################################
sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $genomicseqs,       $querytype,           $queryseqs,
       $queryannotation,   $iidregexp,           $outdb,
       $filter,            $coverage_by_aligned, $options,
       $nonref_regions,    $program,             $seqfetcher_object,
       $seqfetcher_params, $soft_masked_repeats, $use_kill_list,
       $kill_type ) =

      rearrange( [ 'GENOMICSEQS',       'QUERYTYPE',
                   'QUERYSEQS',         'QUERYANNOTATION',
                   'IIDREGEXP',         'OUTDB',
                   'FILTER',            'COVERAGE_BY_ALIGNED',
                   'OPTIONS',           'NONREF_REGIONS',
                   'PROGRAM',           'SEQFETCHER_OBJECT',
                   'SEQFETCHER_PARAMS', 'SOFT_MASKED_REPEATS',
                   'USE_KILL_LIST',     'KILL_TYPE'  
                 ],
                 @args );

  #config precedence is default, param hash, constructor and finally logic_name config
  #read default config entries and do checks
  parse_config($self, $EXONERATE_CONFIG_BY_LOGIC, 'DEFAULT');

  #parse input ID
  my($genomic_slice_name, $query_acc) = split /:_:/,$self->input_id;
  throw "Must have genomic slice as first half of input_id split by ::" if not defined $genomic_slice_name;
  throw "Must have query acc or file name as second half of input_id split by ::" if not defined $query_acc;
  $self->genomic_slice_name($genomic_slice_name);
  $self->query_acc($query_acc);

  #Defaults are over-ridden by parameters given in analysis table...
  my $ph = $self->parameters_hash;
  $self->GENOMICSEQS($ph->{-genomicseqs})                   if $ph->{-genomicseqs};
  $self->QUERYTYPE($ph->{-querytype})                       if $ph->{-querytype};
  $self->QUERYSEQS($ph->{-queryseqs})                       if $ph->{-queryseqs};
  $self->QUERYANNOTATION($ph->{-queryannotation})           if $ph->{-queryannotation};
  $self->IIDREGEXP($ph->{-iidregexp})                       if $ph->{-iidregexp};
  $self->OUTDB($ph->{-outdb})                               if $ph->{-outdb};
  $self->FILTER($ph->{-filter})                             if $ph->{-filter};
  $self->COVERAGE_BY_ALIGNED($ph->{-coverage_by_aligned})   if $ph->{-coverage_by_aligned};
  $self->OPTIONS($ph->{-options})                           if $ph->{-options};
  $self->NONREF_REGIONS($ph->{-nonref_regions})             if $ph->{-nonref_regions};
  $self->PROGRAM($ph->{-program})                           if $ph->{-program};
  $self->SEQFETCHER_OBJECT($ph->{-seqfetcher_object})       if $ph->{-seqfetcher_object};
  $self->SEQFETCHER_PARAMS($ph->{-seqfetcher_params})       if $ph->{-seqfetcher_params};
  $self->SOFT_MASKED_REPEATS($ph->{-soft_masked_repeats})   if $ph->{-soft_masked_repeats};
  $self->USE_KILL_LIST($ph->{-use_kill_list})               if $ph->{-use_kill_list};
  $self->KILL_TYPE($ph->{-kill_type})                       if $ph->{-kill_type};

  #...which are over-ridden by constructor arguments.
  $self->GENOMICSEQS($genomicseqs);
  $self->QUERYTYPE($querytype);
  $self->QUERYSEQS($queryseqs);
  $self->QUERYANNOTATION($queryannotation);
  $self->IIDREGEXP($iidregexp);
  $self->OUTDB($outdb);
  $self->FILTER($filter);
  $self->COVERAGE_BY_ALIGNED($coverage_by_aligned);
  $self->OPTIONS($options);
  $self->NONREF_REGIONS($nonref_regions);
  $self->PROGRAM($program);
  $self->SEQFETCHER_OBJECT($seqfetcher_object);
  $self->SEQFETCHER_PARAMS($seqfetcher_params);
  $self->SOFT_MASKED_REPEATS($soft_masked_repeats);
  $self->USE_KILL_LIST($use_kill_list);
  $self->KILL_TYPE($kill_type);

  #Finally, analysis specific config
  #use uc as parse_config call above switches logic name to upper case
  $self->GENOMICSEQS(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{GENOMICSEQS});
  $self->QUERYTYPE(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{QUERYTYPE});
  $self->QUERYSEQS(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{QUERYSEQS});
  $self->QUERYANNOTATION(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{QUERYANNOTATION});
  $self->IIDREGEXP(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{IIDREGEXP});
  $self->OUTDB(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{OUTDB});
  $self->FILTER(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{FILTER});
  $self->COVERAGE_BY_ALIGNED(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{COVERAGE_BY_ALIGNED});
  $self->OPTIONS(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{OPTIONS});
  $self->NONREF_REGIONS(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{NONREF_REGIONS});
  $self->PROGRAM(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{PROGRAM});
  $self->SEQFETCHER_OBJECT(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{SEQFETCHER_OBJECT});
  $self->SEQFETCHER_PARAMS(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{SEQFETCHER_PARAMS});
  $self->SOFT_MASKED_REPEATS(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{SOFT_MASKED_REPEATS});
  $self->USE_KILL_LIST(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{USE_KILL_LIST});
  $self->KILL_TYPE(${$EXONERATE_CONFIG_BY_LOGIC}{uc($self->analysis->logic_name)}{KILL_TYPE});

  $self->read_and_check_config($EXONERATE_CONFIG_BY_LOGIC);

  #NOTE: use of kill list is not currently implemented in this module
  #the parameters are read from the config file but are not used
  if($self->USE_KILL_LIST or $self->KILL_TYPE){
    throw("Use of kill list not implemented in this module.\n");
  } 

  return $self;
}


sub fetch_input {
  my( $self) = @_;

  my $logic = $self->analysis->logic_name;

  ##########################################
  # set up the target (genome)
  ##########################################

  #write genomic seq to temp file
  my $genomic_file = create_file_name("genomic_" ,"fa" , $self->ANALYSIS_WORK_DIR);

  #repeat masking logic names
  throw "Repeat logic names are not in an array\n" if(!(ref($self->SOFT_MASKED_REPEATS) eq "ARRAY"));

  foreach my $repeat_logic_name ( @{ $self->SOFT_MASKED_REPEATS } ) {
    my $repeat_analysis =
      $self->db->get_AnalysisAdaptor->fetch_by_logic_name($repeat_logic_name);
    if ( !$repeat_analysis ) {
      throw(   "Failed to find an analysis with the logic_name "
             . $repeat_logic_name . ". Cannot continue." );
    }
  }

  my $genomic_seq = $self->db->get_SliceAdaptor->fetch_by_name($self->genomic_slice_name)->
  get_repeatmasked_seq($self->SOFT_MASKED_REPEATS,1);

  write_seqfile($genomic_seq, $genomic_file, 'fasta'); 

  ##########################################
  # set up the query (est/cDNA/protein)
  ##########################################
  my $query_file;
  my $delete_query;
  # check if QUERYSEQ dir exists and file exists
  if(defined $self->QUERYSEQS){
    if(-d $self->QUERYSEQS && -e $self->QUERYSEQS.'/'.$self->query_acc){
      # read use this file as seq file 
      print "Using existing file\n";
      $query_file = $self->QUERYSEQS.'/'.$self->query_acc;
      $delete_query = 0;
    }
  }
  #retrieve and write seq 
  else{
    my $seqfetcher = $self->seqfetcher;
    my $query_seq  = $seqfetcher->get_Seq_by_acc( $self->query_acc );
    if(!$query_seq){
      throw("No entry in sequence index for ".$self->query_acc."\n");
    }
    $query_file    = create_file_name( "query_", "fa", $self->ANALYSIS_WORK_DIR );

    write_seqfile($query_seq, $query_file, 'fasta');
    $delete_query = 1;
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
      throw("CONFLICT: You have defined -program in your config file and ".
            "-program_file in your analysis table.");
    }
  }

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
        ->new(
              -program  => $self->PROGRAM ? $self->PROGRAM : $self->analysis->program_file,
              -analysis => $self->analysis,
              -target_file    => $genomic_file,
              -query_type     => $self->QUERYTYPE,
              -query_file     => $query_file,
              -annotation_file => $self->QUERYANNOTATION ? $self->QUERYANNOTATION : undef,
              %parameters,
              );
    $runnable->files_to_delete($genomic_file);
    if($delete_query){
      $runnable->files_to_delete($query_file);
    }
    $self->runnable($runnable);

}

############################################################

sub run {

  my ($self) = @_;
  my @results;

  throw("Can't run - no runnable objects") unless ( $self->runnable );

  foreach my $runnable ( @{ $self->runnable } ) {

    $runnable->run;
    push( @results, @{ $runnable->output } );
  }

  if ( $self->filter ) {
    my $filtered_transcripts = $self->filter->filter_results( \@results );
    @results = @$filtered_transcripts;
  }

  my @genes = $self->make_genes(@results);

  $self->output( \@genes );
}


############################################################

sub write_output {
  my ( $self, @output ) = @_;

  my $outdb        = $self->get_output_db;
  my $gene_adaptor = $outdb->get_GeneAdaptor;

  unless (@output) {
    @output = @{ $self->output };
  }

  my $fails = 0;
  my $total = 0;
  foreach my $gene (@output) {

    eval { $gene_adaptor->store($gene); };
    if ($@) {
      warning("Unable to store gene!!\n$@");
      $fails++;
    }
    $total++;
  }
  if ( $fails > 0 ) {
    throw(   "Not all genes could be written successfully "
           . "($fails fails out of $total)" );
  }
} ## end sub write_output

############################################################

sub make_genes {
  my ( $self, @transcripts ) = @_;

  my (@genes);

  my $slice_adaptor = $self->db->get_SliceAdaptor;

  my %genome_slices;

  foreach my $tran (@transcripts) {
    $tran->analysis( $self->analysis );
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis( $self->analysis );
    $gene->biotype( $self->analysis->logic_name );

    ############################################################
    # put a slice on the transcript

    my $slice_id = $tran->start_Exon->seqname;
    if ( not exists $genome_slices{$slice_id} ) {
      # assumes genome seqs were named in the Ensembl API Slice naming
      # convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
      $genome_slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
    }
    my $slice = $genome_slices{$slice_id};

    foreach my $exon ( @{ $tran->get_all_Exons } ) {
      $exon->slice($slice);
      foreach my $evi ( @{ $exon->get_all_supporting_features } ) {
        $evi->slice($slice);
        $evi->analysis( $self->analysis );
      }
    }
    foreach my $evi ( @{ $tran->get_all_supporting_features } ) {
      $evi->slice($slice);
      $evi->analysis( $self->analysis );
    }

    if ( !$slice ) {
      my ($sf);

      if ( @{ $tran->get_all_supporting_features } ) {
        ($sf) = @{ $tran->get_all_supporting_features };
      } else {
        my @exons = @{ $tran->get_all_Exons };
        ($sf) = @{ $exons[0]->get_all_supporting_features };
      }
      print $sf->hseqname . "\t$slice_id\n";
    }

    throw("Have no slice") if ( !$slice );
    $tran->slice($slice);
    $gene->add_Transcript($tran);
    push( @genes, $gene );
  } ## end foreach my $tran (@transcripts)
  return @genes;
} ## end sub make_genes

############################################################

sub get_chr_names {
  my ($self) = @_;
  my @chr_names;
  my @chromosomes;

  my $chr_adaptor = $self->db->get_SliceAdaptor;
  #also fetching non-reference regions like DR52 for human by default.
  #specify in E2G config-file.
  if ( defined( $self->NONREF_REGIONS ) ) {
    @chromosomes = @{ $chr_adaptor->fetch_all( 'toplevel', undef, 1 ) };
  } else {
    @chromosomes = @{ $chr_adaptor->fetch_all('toplevel') };
  }

  foreach my $chromosome (@chromosomes) {
    push( @chr_names, $chromosome->seq_region_name );
  }

  return @chr_names;
}


############################################################

sub seqfetcher {
  my ( $self, $arg ) = @_;

  if ($arg) {
    throw(   "RunnableDB::Exonerate2GenesRegion " . $arg
           . " must have a method get_Seq_by_acc" )
      unless ( $arg->can("get_Seq_by_acc") );
    $self->{seqfetcher} = $arg;
  }
  if ( !$self->{seqfetcher} ) {
    $self->require_module( $self->SEQFETCHER_OBJECT );
    my %params = %{ $self->SEQFETCHER_PARAMS };
    print $params{-db}->[0], "\n";
    $self->{seqfetcher} = $self->SEQFETCHER_OBJECT->new( %params, );
  }
  return $self->{seqfetcher};
}


sub get_output_db {
  my ($self) = @_;

  my $outdb;

  if ( $self->OUTDB ) {
    if ( ref( $self->OUTDB ) =~ m/HASH/ ) {

      $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $self->OUTDB },
                                                   -dnadb => $self->db );
    } else {
      $outdb = $self->get_dbadaptor( $self->OUTDB );
    }
  } else {
    $outdb = $self->db;
  }

  return $outdb;
}



############################################################
#
# get/set methods
#
############################################################


sub query_seqs {
  my ( $self, @seqs ) = @_;
  if (@seqs) {
    unless ( $seqs[0]->isa("Bio::PrimarySeqI") || $seqs[0]->isa("Bio::SeqI") )
    {
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    push( @{ $self->{_query_seqs} }, @seqs );
  }
  return @{ $self->{_query_seqs} };
}

############################################################

sub genomic {
  my ( $self, $seq ) = @_;
  if ($seq) {
    unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI") ) {
      throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_genomic} = $seq;
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


#############################################################
# Declare and set up config variables
#############################################################

sub read_and_check_config {
  my $self = shift;

  ##########
  # CHECKS
  ##########
  my $logic = $self->analysis->logic_name;

  # check that compulsory options have values
  foreach my $config_var (qw(QUERYTYPE
                             GENOMICSEQS)) {

   throw("You must define $config_var in config for logic '$logic'")
        if not defined $self->$config_var;
  }

  my %index_params = %{$self->SEQFETCHER_PARAMS};

  if ( ( not -d $index_params{-db}->[0] ) and ( not -d $self->QUERYSEQS ) ) {
    throw(   "Need to define SEQFETCHER_PARAMS -db dir or QUERYSEQS dir "
           . "(with query seq file name or acc in input id after ::).\n" );
  }

  #check that, if specified, query annotation is readable and contains an entry for the
  #current query seq. If it tries to use annotation and there isn't an entry, Exonerate 
  #may produce models on both strands.
  if($self->QUERYANNOTATION and not -e $self->QUERYANNOTATION){
    throw("QUERYANNOTATION '" . $self->QUERYANNOTATION . "' in config must be readable");
  }
  elsif($self->QUERYANNOTATION){
    my $acc = $self->query_acc;
    open(F, $self->QUERYANNOTATION) or throw("Could not open supplied annotation file for reading");
    my $unmatched = 1;
    LINE: while(<F>) {
      my @fields = split;
      $unmatched = 0 if $fields[0] eq $acc;
      last LINE if !$unmatched;
    }
    close(F);
    throw ("No entry for ".$self->query_acc." in supplied query annotation 
    file:".$self->QUERYANNOTATION."\n") if $unmatched;
  }

  # filter does not have to be defined, but if it is, it should
  # give details of an object and its parameters
  if ($self->FILTER) {
    if (not ref($self->FILTER) eq "HASH" or
        not exists($self->FILTER->{OBJECT}) or
        not exists($self->FILTER->{PARAMETERS})) {

      throw( "FILTER in config fo '$logic' must be a hash ref with elements:\n"
           . "  OBJECT : qualified name of the filter module;\n"
           . "  PARAMETERS : anonymous hash of parameters to pass to the filter"
      );

    } else {
      my $module = $self->FILTER->{OBJECT};
      my $pars   = $self->FILTER->{PARAMETERS};

      (my $class = $module) =~ s/::/\//g;
      eval{
        require "$class.pm";
      };
      throw("Couldn't require ".$class." E2G:require_module $@") if($@);

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

sub genomic_slice_name {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_CONFIG_genomic_slice_name'} = $value;
  }

  if (exists($self->{'_CONFIG_genomic_slice_name'})) {
    return $self->{'_CONFIG_genomic_slice_name'};
  } else {
    return undef;
  }
}

sub query_acc {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_CONFIG_query_acc'} = $value;
  }

  if (exists($self->{'_CONFIG_query_acc'})) {
    return $self->{'_CONFIG_query_acc'};
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

sub SEQFETCHER_OBJECT {
  my ($self, $value) = @_;
  if($value){
    $self->{'SEQFETCHER_OBJECT'} = $value;
  }
  return $self->{'SEQFETCHER_OBJECT'};
}

sub SEQFETCHER_PARAMS {
  my ($self, $value) = @_;
  if($value){
    $self->{'SEQFETCHER_PARAMS'} = $value;
  }
  return $self->{'SEQFETCHER_PARAMS'};
}

sub SOFT_MASKED_REPEATS {
  my ($self, $value) = @_;
  if($value){
    $self->{'SOFT_MASKED_REPEATS'} = $value;
  }
  return $self->{'SOFT_MASKED_REPEATS'};
}

sub ANALYSIS_WORK_DIR {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_ana_dir} = $val;
  }
  return $self->{_ana_dir};
}


###############################################
###     end of config
###############################################




1;
