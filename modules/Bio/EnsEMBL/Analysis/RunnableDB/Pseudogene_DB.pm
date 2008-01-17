# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB.pm

=head1 SYNOPSIS

my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB->new(
						-db => $db_adaptor,
						-input_id => $slice_id,		
					        -analysis => $analysis,
								       );

$runnabledb->fetch_input();
$runnabledb->run();
my @array = @{$runnabledb->output};
$runnabledb->write_output();

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Analysis::Runnable::Pseudogene.pm 

Opens connections to 2 dbs:
1 for repeat sequences (GB_DB)
1 for fetching genes from (GB_FINAL)

fetches all the genes on the slice, all the repeats associtaed with each gene and 
collects alignment feature evidence for single exon genes and passes them to the 
runnable.

This module forms the base class of the pseudogene analysis for the gene build,
it will identify obvious pseudogenes but will also flag genes that look
interesting for analysis by either:
Bio::EnsEMBL::Analysis::RunnableDB::Spliced_elsewhere or
Bio::EnsEMBL::Analysis::RunnableDB::PSILC
PSILC will work on any gene with a PFAM domain and will attempt to predict if
it is a pseudogene - it underpredicts pseudogenes but is useful for sequences
that look  dodgy but have no overt evidence of being pseudogenes.
Spliced_elsewhere tests for retrotransposition and tends to be run over single
exon genes.
Configuration for all three of these modules is here:
zBio::EnsEMBL::Analysis::Config::Pseudogene


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 APPENDIX



=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB;

use strict;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::Pseudogene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor;
use Bio::EnsEMBL::Pipeline::Flag;
#use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases qw(DATABASES DNA_DBNAME);
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Blessed;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Combined;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning 
                                      stack_trace);
use Data::Dumper;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);



=head2 fetch_input

  Title   :   fetch_input
  Usage   :   $self->fetch_input
  Function:   Fetches input data for Pseudogene.pm from the database
  Returns :   none
  Args    :   none

=cut

sub fetch_input {
  my( $self) = @_;
  throw("No input id") unless defined($self->input_id);

  my $results = [];		# array ref to store the output
  my %repeat_blocks;
  my %homolog_hash;
  my @transferred_genes;

#  print "Loading database ".$GB_DBNAME.":".$GB_DBHOST."\n";
#  my $rep_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
#    (
#     '-host'   => $GB_DBHOST,
#     '-user'   => $GB_DBUSER,
#     '-dbname' => $GB_DBNAME,
#     '-pass'   => $GB_DBPASS,
#     '-port'   => $GB_DBPORT,
#    );
  #now using Analysis:Databases
  print "Loading reference database : REFERENCE_DB.\n";
  my $rep_db = $self->get_dbadaptor("REFERENCE_DB") ; 
  #my $rep_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $$DATABASES{REFERENCE_DB} });

  #store repeat db internally
  $self->rep_db($rep_db);
  my $rsa = $rep_db->get_SliceAdaptor;

  #genes come from final genebuild database
#  my $genes_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
#    (
#     '-host'   => $GB_FINALDBHOST,
#     '-user'   => $GB_FINALDBUSER,
#     '-dbname' => $GB_FINALDBNAME,
#     '-pass'   => $GB_FINALDBPASS,
#     '-port'   => $GB_FINALDBPORT,
#    );
  #now using Analysis:Databases
  print "Loading genes database : PS_INPUT_DATABASE => $PS_INPUT_DATABASE ( defined in Databases.pm ) \n";
  #my $genes_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $$DATABASES{GENEBUILD_DB} });
  my $genes_db = $self->get_dbadaptor("$PS_INPUT_DATABASE") ; 

  $self->gene_db($genes_db);
  #genes are written to the pseudogene database
  # genes_slice holds all genes on input id slice
 
  my $genedb_sa = $self->gene_db->get_SliceAdaptor; 
  print  "DB NAME: ".$self->db->dbc->dbname."\n";
  my $genes_slice = $genedb_sa->fetch_by_name($self->input_id);
  $self->query($genes_slice);

  my $genes = $genes_slice->get_all_Genes;
  print  $genes_slice->name."\t".
    $genes_slice->start."\n";
  GENE: foreach my $gene (@{$genes}) {
    # Ignore genes that are not protein_coding
    unless ( $gene->biotype eq $PS_BIOTYPE ) {
      $self->ignored_genes($gene);
      next GENE;
    }


    ############################################################################
    # transfer gene coordinates to entire chromosome to prevent problems arising
    # due to offset with repeat features 
    my $chromosome_slice = $rsa->fetch_by_region(
						 'toplevel',
						 $genes_slice->chr_name,
						);

    my $transferred_gene = $gene->transfer($chromosome_slice);
    $self->lazy_load($transferred_gene);
    push @transferred_genes,$transferred_gene;

    # repeats come from core database
    # repeat slice only covers gene to avoid sorting repeats unnecessarily
    my $rep_gene_slice = $rsa->fetch_by_region(
					       'toplevel',
					       $genes_slice->chr_name,
					       $transferred_gene->start,
					       $transferred_gene->end,
					      );
    # get repeat blocks
    my @feats = @{$rep_gene_slice->get_all_RepeatFeatures};
    @feats = map { $_->transfer($chromosome_slice) } @feats;
    my $blocks = $self->get_all_repeat_blocks(\@feats);
    # make hash of repeat blocks using the gene as the key
    $repeat_blocks{$transferred_gene} = $blocks;
  }

  $self->genes(\@transferred_genes);
  $self->repeat_blocks(\%repeat_blocks);
  $self->make_runnable;

  return 1;
}


sub make_runnable {
  my ($self) = @_;
      
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Pseudogene->new
      ( 
        '-analysis' => $self->analysis,
        '-genes' => $self->genes,
        '-repeat_features' => $self->repeat_blocks,
        );
  $self->runnable($runnable);
}


=head2 get_all_repeat_blocks

  Args       : none
  Description: merges repeats into blocks for each gene
  Returntype : array of Seq_Feature blocks;

=cut 

sub get_all_repeat_blocks {
  my ($self,$repeat_ref) = @_;
  my @repeat_blocks;
  my @repeats = @{$repeat_ref};
  @repeats = sort {$a->start <=> $b->start} @repeats;
  my $curblock = undef;

 REPLOOP: foreach my $repeat (@repeats) {
    my $rc = $repeat->repeat_consensus;
    my $use = 0;
    foreach my $type (@$PS_REPEAT_TYPES){
      if ($rc->repeat_class =~ /$type/) {
	$use = 1;
	last;
      }
    }
    next REPLOOP unless $use;
    if ($repeat->start <= 0) { 
      $repeat->start(1); 
    }
    if (defined($curblock) && $curblock->end >= $repeat->start) {
      if ($repeat->end > $curblock->end) { 
	$curblock->end($repeat->end); 
      }
    } else {
      $curblock = Bio::EnsEMBL::Feature->new(
						-START => $repeat->start,
						-END => $repeat->end, 
						-STRAND => $repeat->strand
					    );
      push (@repeat_blocks,$curblock);
    }
  }
    @repeat_blocks = sort {$a->start <=> $b->start} @repeat_blocks;
  return\@repeat_blocks;
}

=head2 write_output

  Args       : none
  description: writes gene array into db specified in Bio::EnsEMBL::Config::GeneBuild::Databases.pm
  exception  : warns if it cannot write gene
  Returntype : none

=cut 


sub write_output {
my($self) = @_;
  my $genes = $self->output;
  push @{$genes},@{$self->ignored_genes} if $self->ignored_genes;
  my %feature_hash;
  #  empty_Analysis_cache();
  # write genes out to a different database from the one we read genes from.

#  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
#					      '-host'   => $PSEUDO_DBHOST,
#					      '-user'   => $PSEUDO_DBUSER,
#					      '-dbname' => $PSEUDO_DBNAME,
#					      '-pass'   => $PSEUDO_DBPASS,
#					      '-port'   => $PSEUDO_DBPORT,
#					     );

# my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $$DATABASES{PSEUDO_DB} });
  my $db = $self->get_dbadaptor("$PS_OUTPUT_DATABASE") ; 

  #now using Analysis:Databases
  print "Writing to database PS_OUTPUT_DATABASE : $PS_OUTPUT_DATABASE\n"; 


  # sort out analysis
  my $analysis = $self->analysis;
  unless ($analysis){
    throw("an analysis logic name must be defined in the command line");
  }
  my $gene_adaptor = $db->get_GeneAdaptor;
  foreach my $gene (@{$genes}) {
    # store gene
    eval {
      $gene_adaptor->store($self->lazy_load($gene));
      print STDERR  "wrote gene " . $gene->dbID . " to database ".
        $gene_adaptor->db->dbname."\n";
    };
    if ( $@ ) {
      warning("UNABLE TO WRITE GENE:\n$@");
    }
  }

  return 1;
} 


=head2 run

  Args       : none
  Description: overrides runnableDb run method to allow gene objects to be validated 
before runnning the runnable
  Returntype : scalar

=cut 

sub run  {
  my ($self) = @_;
  foreach my $runnable (@{$self->runnable}) {
    throw("Runnable module not set") unless ($runnable->isa("Bio::EnsEMBL::Analysis::Runnable"));
    $runnable->run();
    $self->output($runnable->output);
    if ($SINGLE_EXON){
      $self->store_ids($runnable->single_exon_genes,$SPLICED_ELSEWHERE_LOGIC_NAME);
    }    
    if ($INDETERMINATE){
      $self->store_ids($runnable->indeterminate_genes,$PSILC_LOGIC_NAME);
    }
  }
  return 0;
}

=head2 store_ids

Arg [none] :
  Description: stores gene dbIDS of genes being held back for further analyses
  Returntype : scalar
  Exceptions : throws if it cannot open the file
  Caller     : general

=cut

sub store_ids {
  my ($self, $id_list,$analysis_name) = @_;

  my $flag_adaptor = Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor->new($self->db);
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
  my $analysis = $analysis_adaptor->fetch_by_logic_name($analysis_name);
  unless ($analysis) {  
     throw( "analysis " . $analysis_name  . " can't be found in " . $self->db->dbname ." @ " . $self->db->host.  "\n")  ; 
  } 
# What do you do if the analysis thing isnt found?
  return 0 unless (scalar(@{$id_list} >0));
  foreach my $id(@{$id_list}){
    my $flag = Bio::EnsEMBL::Pipeline::Flag->new(
						 '-type'         => 'gene',
						 '-ensembl_id'   => $id,
						 '-goalAnalysis' => $analysis,
						);
    $flag_adaptor->store($flag);
  }
  return;
}

=head2 lazy_load

  Arg [1]    : Bio::EnsEMBL::Gene
  Description: forces lazy loading of transcripts etc.s
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub lazy_load {
  my ($self, $gene) = @_;
  if ($gene){
    unless ($gene->isa("Bio::EnsEMBL::Gene")){
      throw("gene is not a Bio::EnsEMBL::Gene, it is a $gene");
    }
    foreach my $trans(@{$gene->get_all_Transcripts}){
      my $transl = $trans->translation; 
       $trans->get_all_supporting_features() ; 
      if ($transl){
	$transl->get_all_ProteinFeatures;
      }
    }
  }
  return $gene;
}

=head2 remove_transcript_from_gene

  Args       : Bio::EnsEMBL::Gene object , Bio::EnsEMBL::Transcript object
  Description: steves method for removing unwanted transcripts from genes
  Returntype : scalar

=cut 

sub _remove_transcript_from_gene {
  my ($self, $gene, $trans_to_del)  = @_;
  # check to see if it is a blessed transcript first
  foreach my $blessed ( @{$GB_BLESSED_GENETYPES} ) {
    if ( $trans_to_del->biotype eq $blessed->{'type'} or 
	 $trans_to_del->biotype eq $GB_BLESSED_COMBINED_GENETYPE ) {
      # transcript is blessed dont delete it
      return 'BLESSED';
    }
  }

  my @newtrans;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans != $trans_to_del) {
      push @newtrans,$trans;
    }
  }

  # The naughty bit!
  $gene->{_transcript_array} = [];

  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }

  return;
}

=head2 transcript_to_keep

  Args       : Bio::EnsEMBL::Transcript object
  Description: removes the translation provided it is not a blessed transcript
  Returntype : scalar

=cut 


sub transcript_to_keep {
  my ($self, $trans_to_keep)  = @_;
  foreach my $blessed ( @{$GB_BLESSED_GENETYPES} ) {
    if ( $trans_to_keep->biotype eq $blessed->{'type'} or 
	 $trans_to_keep->biotype eq $GB_BLESSED_COMBINED_GENETYPE ) {
      # transcript is blessed dont delete the translation
      return;
    }
  }
  $trans_to_keep->translation(undef);
  return;
}


############################################################################
# container methods

=head2 gene_db

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set gene db adaptor
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub gene_db {
  my ($self, $gene_db) = @_;
  if ($gene_db){
    unless ($gene_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
      throw("gene db is not a Bio::EnsEMBL::DBSQL::DBAdaptor, it is a $gene_db");
    }
    $self->{'_gene_db'} = $gene_db;
  }
  return $self->{'_gene_db'};
}

=head2 rep_db

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set gene db adaptor
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub rep_db {
  my ($self, $rep_db) = @_;
  if ($rep_db){
    unless ($rep_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
      throw("gene db is not a Bio::EnsEMBL::DBSQL::DBAdaptor, it is a $rep_db");
    }
    $self->{'_rep_db'} = $rep_db;
  }
  return $self->{'_rep_db'};
}

=head2 genes

  Arg [1]    : array ref
  Description: get/set genescript set to run over
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub genes {
  my ($self, $genes) = @_;
  if ($genes) {
    foreach my $gene (@{$genes}) {
      unless  ($gene->isa("Bio::EnsEMBL::Gene")){
	throw("Input isn't a Bio::EnsEMBL::Gene, it is a $gene\n$@");
      }
    }
    $self->{'_genes'} = $genes;
  }
  return $self->{'_genes'};
}


=head2 repeat_blocks

  Arg [1]    : array ref
  Description: get/set genescript set to run over

=cut

sub repeat_blocks {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_repeat_blocks} = $val;
  }
  return $self->{_repeat_blocks};
}



=head2 pseudo_genes

  Arg [1]    : Bio::EnsEMBL::Gene
  Description: get/set for pseudogenes 
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub pseudo_genes {
  my ($self, $pseudo_gene) = @_;
  if ($pseudo_gene) {
    unless ($pseudo_gene->isa("Bio::EnsEMBL::Gene")){
      throw("pseudo gene is not a Bio::EnsEMBL::Gene, it is a $pseudo_gene");
    }
    push @{$self->{'_pseudo_gene'}},$self->lazy_load($pseudo_gene);
  }
  return $self->{'_pseudo_gene'};
}

=head2 real_genes

Arg [1]    : Bio::EnsEMBL::Gene
 Description: get/set for 'functional' genes
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub real_genes {
  my ($self, $real_gene) = @_;
  if ($real_gene) {
    unless ($real_gene->isa("Bio::EnsEMBL::Gene")){
      throw("real gene is not a Bio::EnsEMBL::Gene, it is a $real_gene");
    }
    push @{$self->{'_real_gene'}},$self->lazy_load($real_gene);
  }
  return $self->{'_real_gene'};
}


=head2 ignored_genes

Arg [1]    : Bio::EnsEMBL::Gene
 Description: get/set for genes that pseudogene does not check
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub ignored_genes {
  my ($self, $ignored_gene) = @_;
  if ($ignored_gene) {
    unless ($ignored_gene->isa("Bio::EnsEMBL::Gene")){
      throw("ignored gene is not a Bio::EnsEMBL::Gene, it is a $ignored_gene");
    }
    push @{$self->{'_ignored_gene'}},$self->lazy_load($ignored_gene);
  }
  return $self->{'_ignored_gene'};
}



1;
