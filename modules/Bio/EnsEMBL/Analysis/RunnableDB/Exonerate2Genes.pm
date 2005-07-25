
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes;




=head1 SYNOPSIS

my $exonerate2genes = Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes->new(
                              -db         => $refdb,
			      -analysis   => $analysis_obj,
			      -database   => $EST_GENOMIC,
			      -query_seqs => \@sequences,
			      -query_type => 'dna',
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


=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Gene;

use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);


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
  my $target = $self->GENOMICSEQS;

  if(ref $target eq 'ARRAY'){
    # genome is in multiple directories;  the order of directories determines
    # which file is used in case of duplicates. New versions should therefore
    # be in the directory listed first.
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
  }
  else{
    if (-e $target and -d $target) {
      # genome is in a directory; the directory must contain the complete
      # genome else we cannot do best-in-genome filtering. 
      foreach my $chr_name ($self->get_chr_names) {
	if (-s "$target/$chr_name.fa") {
	  push @db_files, "$target/$chr_name.fa";
	} else {
	  warning( "Could not find fasta file for '$chr_name' in '$target'\n");
	}
      }
    }
    elsif (-e $target and -s $target) {
      # genome sequence is in a single file
      @db_files = ($target);
    } else {
      throw("'$target' refers to something that could not be made sense of");
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
  }
  elsif (-e $query and -s $query) {
    # query seqs is a single file; input id will correspond to a chunk number
    $query_file = $query;
    my $iid_regexp = $self->IIDREGEXP;
    
    throw("When your input ids are not filenames, you must define ".
          "IIDREGEXP in config to enable inference of chunk number and total")
        if not defined $iid_regexp;

    ($chunk_number, $chunk_total) = $self->input_id =~ /$iid_regexp/;
  } else {
    throw("'$query' refers to something that could not be made sense of\n");
  }

  ##########################################
  # load and create filterm is there is one
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


  foreach my $database ( @db_files ){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript
        ->new(
              -program  => $self->analysis->program_file,
              -analysis => $self->analysis,
              -target_file    => $database,
              -query_type     => $self->QUERYTYPE,
              -query_file     => $query_file,
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
    $gene->type($self->analysis->logic_name);
    
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
    $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{$self->OUTDB},
                                                -dnadb => $self->db);
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

  # output db does not have to be defined, but if it is, it should be a hash
  throw("OUTDB in config for '$logic' must be a hash ref of db connection pars.")
      if $self->OUTDB and ref($self->OUTDB) ne "HASH";

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

###############################################
###     end of config
###############################################




1;
