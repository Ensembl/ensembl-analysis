=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::PSILC - 

=head1 SYNOPSIS 
my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::PSILC->new
(
-db => $db_adaptor,
-input_id => flag start:flag end,
-analysis => $analysis,
);

$runnabledb->fetch_input();
$runnabledb->run();
my @array = @{$runnabledb->output};
$runnabledb->write_output();

=head1 DESCRIPTION

Prepares and runs PSILC.
PSILC is a pseudogene predicting program that compares DNA alignments of 
transcripts from the organism of interest and homologous sequences from 3 
other organisms. It uses regions of conservation within the seqence alignment 
defined by PFAM domains to make its predictions.
It uses flagged seqences as input ids that were identified for further 
analysis by Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB.pm or 
Bio::EnsEMBL::Analysis::RunnableDB::Spliced_elsewhere.
Takes start and end flag ids as input, fetches all flags that lie within the 
specified range that are associated with the psilc analysis, retrieves the 
associated gene objects. Identifies all the transcripts in the set that contain 
PFAM domains (requires protein annotationto be run before the pseudogene tests)
Runs a BLASTP search to identify orthologs from the 3 comparison databases, 
then clusters the results using clustalw and passed the alignment to PSILC to 
run over.
Parses the resuls and modifies gene objects to make them into pseudogenes if 
required and then writes the genes to the final database.

There is a script ensembl-personal/sw4/Scripts/Pseudogenes/prepare_PSILC.pl
that prepares the module for pipelining, makes input ids etc.
A second script:
ensembl-personal/sw4/Scripts/Pseudogenes/make_PSILC_bsubs.pl
which will prepare the blast databases to use from the 3 species.

The module is configured through:
Bio::EnsEMBL::Analysis::Config::Pseudogene.pm

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::PSILC;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB;
use Bio::EnsEMBL::Analysis::Runnable::PSILC_BlastP;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Analysis::Runnable::PSILC; 
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


my $runnable;

=head2 fetch_input

  Arg [1]    : none
  Description: Opens and stores connections to the 3 ortholog databases and the genebuild databases, parses the 
flag input id and fetches the genes to run PSILC on. Determines which transcripts contain
PFAM domains and stores the domains internally in a hash that ties the transcript id to the
PFAM id
  Returntype : none
  Exceptions : throws if the input ids are not recognised.
  Caller     : general

=cut


sub fetch_input{
  my ($self)=@_;
  # open the dbs, get the adaptors and store them
  my ($start, $end);
  my $count=0;
  my @genes; 
  my $subjectdb = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $self->PSILC_SUBJECT_DBHOST,
     '-user'   => 'ensro',
     '-dbname' => $self->PSILC_SUBJECT_DBNAME,
     '-port'   => $self->PSILC_SUBJECT_DBPORT,
    );
  my $orth1db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $self->PSILC_ORTH1_DBHOST,
     '-user'   => 'ensro',
     '-dbname' => $self->PSILC_ORTH1_DBNAME,
     '-port'   => $self->PSILC_ORTH1_DBPORT,
    );  
  my $orth2db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $self->PSILC_ORTH2_DBHOST,
     '-user'   => 'ensro',
     '-dbname' => $self->PSILC_ORTH2_DBNAME,
     '-port'   => $self->PSILC_ORTH2_DBPORT,
    ); 

  $self->species_db($self->SUBJECT,$subjectdb);
  $self->species_db($self->ORTH1,$orth1db);
  $self->species_db($self->ORTH2,$orth2db);

  #store repeat db internally
  my $dna_db = $self->get_dbadaptor($self->DNA_DBNAME) ;
  $self->rep_db($dna_db);

  #genes come from final genebuild database
  my $genes_db = $self->get_dbadaptor("GENEBUILD_DB");
  $self->gene_db($genes_db);
 
  my $ga = $genes_db->get_GeneAdaptor;
  my $fa = Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor->new($self->db);
  my $ids = $fa->fetch_by_analysis($self->analysis);
  $self->throw("No flags found for analysis $self->PSILC_LOGIC_NAME\n")  unless (scalar(@{$ids}>0));
  if ($self->input_id =~ /(\d+):(\d+)/) {
    $start = $1;
    $end = $2;
  } else {
    $self->throw("Input id not recognised\n");
  }
  # get ids
  foreach my $flag (@{$ids}) {
    my $gene;
    if ($flag->dbID >= $start && $flag->dbID <= $end) {
      $count++;
      if ($flag->ensembl_id =~ /\w+/){
	$gene = $ga->fetch_by_stable_id($flag->ensembl_id);
      }
      else {
	$gene = $ga->fetch_by_dbID($flag->ensembl_id);
      }
      push @genes, $self->lazy_load($gene);
    }
  }
  $self->genes(\@genes);

  # get all the genes out of your core database with pfam domain hits
  # If $REP_TRANSCRIPT flag is set only select one transcript per gene 
  # the longest translation that has a pfam domain

  my $pfam_transcripts =0;
  my %domains;

 GENE: foreach my $gene (@genes) {  
    my %already_seen;
    my @rep_transcripts;
    my $domain;
  TRANS: foreach my $trans (@{$gene->get_all_Transcripts}) {
      my %already_seen;
      if ($trans->translation) {
	my $pep = $trans->translation;
	#If the transcript has a pfam domain store translation object
	foreach my $pf (@{$pep->get_all_ProteinFeatures}) {
	  if ($pf->analysis->logic_name eq 'Pfam') {
	    # remove the transcript from the array unless it has a pfam domain
	    push @rep_transcripts,$trans;
	    next TRANS;
	  }
	}
      }
    }
    next GENE unless (scalar (@rep_transcripts > 0));
    if (scalar(@rep_transcripts > 1) && $self->REP_TRANSCRIPT){
      # arrange transcripts by length so you can pick the longest as rep if needed
      @rep_transcripts = sort {$a->length <=> $b->length} @rep_transcripts;
      my @temp = pop @rep_transcripts;
      print $temp[0]->stable_id."\n";
      @rep_transcripts = @temp;
    }
    # othewise store domains found in the rep transcript
    foreach my $rep_transcript(@rep_transcripts){
      foreach my $pf (@{$rep_transcript->translation->get_all_ProteinFeatures}) {
	unless ($already_seen{$pf->hseqname}) {
	    if ($pf->analysis->logic_name eq 'Pfam') {
	      my $domain_id = $pf->hseqname;
	      $domain_id =~ s/\.\d+//;
	      $domain .= "$domain_id\t$domain_id\n";
	      $already_seen{$pf->hseqname} =1;
	    }
	  }
      }
      $pfam_transcripts++;
      $domains{$rep_transcript->dbID}=$domain;
    }
  }
  print STDERR "Found $pfam_transcripts Pfam transcripts out of ".scalar(@genes)." genes \n";
  $self->domains(\%domains);
  return 1;
}

=head2 run 

  Arg [1]    : none
  Description: Runs a BLASTP search on each transcript with a PFAM domain, then clusters
the sequences using clustalw follwed by running PSILC
  Returntype : none
  Exceptions : 
  Caller     : general

=cut

sub run{
  my ($self) = @_;
  my $ga = $self->gene_db->get_GeneAdaptor;
  $self->make_dir;
  my @genes = @{$self->genes};
  my $pfam_transcript;
  my %domains = %{$self->domains};

 PFAM:  foreach my $gene(@genes){
    my $pseudo_trans = 0;
    my $pfam_trans = 0;
  TRANSCRIPT: foreach my $transcript (@{$gene->get_all_Transcripts}){
      my $PSILC_results;
      if ($domains{$transcript->dbID}){
        $pfam_transcript = $transcript;
	$pfam_trans++;
      }
      else {
	next TRANSCRIPT;
      }
      if ($pfam_transcript->stable_id) {
	print STDERR "Running analysis for ".$pfam_transcript->stable_id."\n";
      } else {
	print STDERR "Running analysis for ".$pfam_transcript->dbID."\n";
      }
      # run blast comparison against the three databases
      my $blast_results = $self->run_blast($pfam_transcript);
      # fetch the sequences of the homologs identified in blast
      my @homologs = @{$self->fetch_trans($blast_results)};
      if (scalar(@homologs) >= 2) {
	$PSILC_results = $self->run_PSILC($pfam_transcript,\@homologs);
      }
      unless ($PSILC_results) {
	print STDERR "Not enough sequences to continue for transcript ".$pfam_transcript->dbID;
	print STDERR " could only find ".scalar(@homologs)."\n";
	$self->real_genes($gene);
	next PFAM;
      }
      if ($PSILC_results->{'postNMax'} > 0.95){
	$pseudo_trans++;
      }
    }
    if ($pfam_trans > 0 && $pseudo_trans == $pfam_trans){
      $self->pseudo_genes($gene)
    }
    else {
      $self->real_genes($gene);
    }
  }
  $self->output($self->pseudo_genes);
  $self->output($self->real_genes);
  return 1;
}

=head2 fetch_trans

  Arg [1]    : HASH ref
  Description: fetches transcripts from the 3 ortholog databases, ensures sequences
are non-redundant and limits the transcripts to a pre set number of sequences from each species
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : throws unless passed a Bio::EnsEMBL::DBSQL::DBAdaptor, throws if the species name is not recognised
  Caller     : general

=cut

sub fetch_trans{
  my ($self,$blast_results) =@_;
  my %transcript_hash;
  my $transcript_count;
  my @homologs;
  my %already_seen;
 SPECIES: foreach my $species (keys %{$blast_results}){
    my $ta = $self->species_db($species)->get_TranscriptAdaptor;
    print "Fetching $species transcripts\n";
    $transcript_count=0;
  TRANS:   foreach my $trans_id (@{$blast_results->{$species}}){;
      next TRANS if ($already_seen{$trans_id});
      $transcript_count++;
      next SPECIES unless ($transcript_count <= $self->PS_SPECIES_LIMIT);
      my $transcript = $ta->fetch_by_translation_stable_id($trans_id);
      $already_seen{$trans_id} = 1;
      unless ($transcript && $transcript->isa("Bio::EnsEMBL::Transcript")){
	$self->warn("Cannot find transcript $trans_id for species ".$transcript_hash{$trans_id}." I get a $transcript\n$@\n");
      }
      next TRANS unless ($transcript->translateable_seq);
      push @homologs,$transcript;
    }
  }
  return \@homologs;
}

=head2 run_blast

  Arg [1]    : Bio::EnsEMBL::Transcript
  Description: runs the BLASTP runnable
  Returntype : Hash reference
  Exceptions : none
  Caller     : general

=cut

sub run_blast{
  my ($self,$pfam_transcript) = @_;
  my $blast = Bio::EnsEMBL::Analysis::Runnable::PSILC_BlastP->new 
    (
     '-trans'    => $pfam_transcript,
     '-analysis' => $self->analysis,
    );
  $blast->run;
  return $blast->output;
}

=head2 species_db

  Arg [1]    : Bio::EnsEMBL::Transcript, Array contining homolog transcripts
  Description: runs the PSILC runnable
  Returntype : Hash reference
  Exceptions : none
  Caller     : general

=cut

sub run_PSILC{
  my ($self,$pfam_transcript,$homologs) = @_;
  my %domains = %{$self->domains};
  # Run PSILC analysis
  my $PSILC = Bio::EnsEMBL::Analysis::Runnable::PSILC->new 
    (
     '-trans'     => $pfam_transcript,
     '-homologs'  => $homologs,
     '-analysis'  => $self->analysis,
     '-domain'    => \$domains{$pfam_transcript->dbID},
     '-input_id'  => \$self->input_id,
     '-psilc_work_dir' => $self->PSILC_WORK_DIR,
    );
  eval{
    $PSILC->run;
  };
  if ($@){
    $self->warn("PSILC FAILED TO RUN $@\n");
    return 0;
  }
  return $PSILC->output;
}

=head2 species_db

  Arg [1]    : scalar species name, Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set db adaptors for the three species chosen
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : throws unless passed a Bio::EnsEMBL::DBSQL::DBAdaptor, throws if the species name is not recognised
  Caller     : general

=cut

sub species_db{
  my ($self,$species,$db)= @_;
  if ($db){
    unless ($db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
      $self->throw("$species object is not a Bio::EnsEMBL::DBSQL::DBAdaptor, it is a $db$@\n");
    }
    if ($species eq $self->SUBJECT){
      $self->{'_db_subject'} = $db;
    }
    if ($species eq $self->ORTH1){
      $self->{'_db_orth1'} = $db;
    }
    if ($species eq $self->ORTH2){
      $self->{'_db_orth2'} = $db;
    }
  }
  if ($species eq $self->SUBJECT){
    return $self->{'_db_subject'};
  }
  if ($species eq $self->ORTH1){
    return $self->{'_db_orth1'};
  }
  if ($species eq $self->ORTH2){
    return $self->{'_db_orth2'};
  }
  $self->throw("Species not recognised : $species\n");
  return;
}

=head2 make_dir

  Arg [1]    : none
  Description: creates output directory for storing PSILC results
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub make_dir{
  my ($self) = @_;
  my $input_id = $self->input_id;  
  if ($self->PSILC_WORK_DIR){
    system("mkdir $self->PSILC_WORK_DIR/$input_id");
  }
  else{
    $self->throw("Cannot make output directory\n");
  }
  return 1;
}



#####################################################################
# Containers

=head2 domains

  Arg [1]    : Hash
  Description: get/set for hash tying pfam domain identifiers to transcript dbID
  Returntype : Hash
  Exceptions : none
  Caller     : general

=cut


sub domains {
  my ($self, $domains) = @_;
  if ($domains) {
    $self->{'_domains'} = $domains;
  }
  return $self->{'_domains'};
}

=head2 transcripts

  Arg [1]    : Bio::EnsEMBL::Transcript
  Description: get/set for transcripts
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut

sub transcripts {
  my ($self, $transcript) = @_;
  if ($transcript) {
    push @{$self->{'_transcripts'}}, $transcript;
  }
  return $self->{'_transcripts'};
}



1;
