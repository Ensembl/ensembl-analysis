
package Bio::EnsEMBL::Analysis::RunnableDB::Spliced_elsewhere;


use strict;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


sub fetch_input{
  my ($self)=@_;
  # get all the genes out of your core database with pfam domain hits
  my @genes;
  my %parameters;
  if ($self->parameters_hash) {
    %parameters = %{$self->parameters_hash};
  }
  my $runname = "Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere";

  my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_DBHOST,
     '-user'   => $GB_DBUSER,
     '-dbname' => $GB_DBNAME,
     '-pass'   => $GB_DBPASS,
     '-port'   => $GB_DBPORT,
    );
  #store repeat db internally
  $self->dna_db($dna_db);

  #genes come from final genebuild database
  my $genes_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_FINALDBHOST,
     '-user'   => $GB_FINALDBUSER,
     '-dbname' => $GB_FINALDBNAME,
     '-pass'   => $GB_FINALDBPASS,
     '-port'   => $GB_FINALDBPORT,
     '-dnadb'  => $dna_db,
    );

  $self->gene_db($genes_db);

  my $ga = $genes_db->get_GeneAdaptor;

  open(GENES,$PS_GENE_DIR.$self->input_id) or $self->throw("Cannot open file ".$PS_GENE_DIR.$self->input_id."\n");
  while (<GENES>){
    chomp;
    my $gene = $ga->fetch_by_transcript_stable_id($_);
    unless ($gene->isa("Bio::EnsEMBL::Gene")){
      $self->warn("Unable to find gene object by stable id $_\n$@\n");
    }
    else{
      push @genes,$ga->fetch_by_transcript_stable_id($_);
    }
  }
  my $runnable = $runname->new (
				'-genes' => \@genes,
				'-analysis' => $self->analysis,
			       );
  $self->runnable($runnable);
return 1;
}

sub run {
  my($self) = @_;
  foreach my $runnable (@{$self->runnable}) {
    $self->throw("Runnable module not set") unless ($runnable);
    $runnable->run;
    $self->parse_results($runnable->output);
  }
  return 1;
}

sub  parse_results{
  my ($self,$results)=@_;
 RESULT: foreach my $result_array (@{$results}) {
    my $ta = $self->gene_db->get_TranscriptAdaptor;
    my $ga = $self->gene_db->get_GeneAdaptor;
    my $retro_trans = @$result_array[0];
    my $retro_span = $retro_trans->cdna_coding_end- $retro_trans->cdna_coding_start;
    my @dafs =  @$result_array[1];
    @dafs = sort {$b->p_value <=> $a->p_value} @dafs;
  DAF: foreach my $daf (@{$dafs[0]}) {
      # is the percent id above threshold?
      next DAF unless ($daf->percent_id > $PS_PERCENT_ID_CUTOFF);
      my $coverage = int($daf->length/$retro_trans->translate->length*100);
      # is the coverage above the threshold?
      next DAF unless ($coverage > $PS_PERCENT_ID_CUTOFF);

      my $real_trans;
      # Warn if transcript cannot be found
      eval{
	$real_trans =   $ta->fetch_by_translation_stable_id($daf->hseqname);
      };
      if ($@) {
	$self->warn("Unable to find translation $daf->hseqname \n$@\n");
	next;
      }
      # real span is the genomic span of the region of the transcript
      # that aligns to the retro_transcript
      # + and - 2 are designed to chuck things that lie right on the 
      # edge of an exon and so need to overlap each exon by at least 3
      # before you consider it tp be ok.
      # thats a crap number really should be altered but not just yet...

      my $real_span;
      my @genomic_coords = $real_trans->pep2genomic($daf->hstart+3,$daf->hend-3);
      @genomic_coords = sort {$a->start <=> $b->start} @genomic_coords;
      $real_span = $genomic_coords[$#genomic_coords]->end - $genomic_coords[0]->start;


      # Is the span higher than the allowed ratio?
      if($real_span / $retro_span > $PS_SPAN_RATIO ) {
	print STDERR $retro_trans->stable_id." matches " . $daf->hseqname . " at "
	  . $daf->percent_id . " %ID and $coverage % coverage with pseudogene span of $retro_span vs real gene span of $real_span\n";
	# Gene is pseudo, store it internally
	$self->pseudo_genes($ga->fetch_by_transcript_stable_id($retro_trans->stable_id));
	next RESULT;
      }
    }
    $self->real_genes($ga->fetch_by_transcript_stable_id($retro_trans->stable_id));
  }
  return 1; 
}

sub write_output {
  my($self) = @_;
  unless ($self->real_genes){
    $self->throw("Where have all the genes gone?\n$@\n");
  }
  my @genes = @{$self->real_genes};
  push @genes,@{$self->pseudo_genes};
  unless (@genes){
    print STDERR "No genes to store\n";
    return 1;
  }

  # write genes out to a different database from the one we read genes from.
  my $dbname = $PSEUDO_DBNAME;
  my $dbhost = $PSEUDO_DBHOST;
  my $dbuser = $PSEUDO_DBUSER;
  my $dbpass = $PSEUDO_DBPASS;
  my $dbport = $PSEUDO_DBPORT;
  
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					      '-host'   => $dbhost,
					      '-user'   => $dbuser,
					      '-dbname' => $dbname,
					      '-pass'   => $dbpass,
					      '-port'   => $dbport,
					     );
  # sort out analysis
  my $analysis = $self->analysis;
  print STDERR $analysis->logic_name."\n";
  unless ($analysis){
    $self->throw("an analysis logic name must be defined in the command line");
  }

  my $gene_adaptor = $db->get_GeneAdaptor;
  foreach my $gene (@genes) { 
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->translation;
    }

    # store
    eval {
      $gene_adaptor->store($gene);
      print STDERR  "wrote gene " . $gene->dbID . " to database ".
	$gene_adaptor->db->dbname."\n";
    };
    if ( $@ ) {
      $self->warn("UNABLE TO WRITE GENE:\n$@");
    }
  }
}



=head2 remove_transcript_from_gene

  Args       : Bio::EnsEMBL::Gene object , Bio::EnsEMBL::Transcript object
  Description: steves method for removing unwanted transcripts from genes
  Returntype : scalar

=cut 


sub _remove_transcript_from_gene {
  my ($self, $gene, $trans_to_del)  = @_;

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

  return scalar(@newtrans);
}



########################################################
# Containers

=head2 dna_db

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set gene db adaptor
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub dna_db {
  my ($self, $dna_db) = @_;
  
  unless ($dna_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
  $self->throw("gene db is not a Bio::EnsEMBL::DBSQL::DBAdaptor, it is a $dna_db");
}
  $self->{'_dna_db'} = $dna_db;
  return $self->{'_dna_db'};
}

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
      $self->throw("gene db is not a Bio::EnsEMBL::DBSQL::DBAdaptor, it is a $gene_db");
    }
    $self->{'_gene_db'} = $gene_db;
  }
  return $self->{'_gene_db'};
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
  if ($pseudo_gene){
    unless ($pseudo_gene->isa("Bio::EnsEMBL::Gene")){
      $self->throw("pseudo gene is not a Bio::EnsEMBL::Gene, it is a $pseudo_gene");
    }
    $pseudo_gene->type('pseudogene');
    my @pseudo_trans = @{$pseudo_gene->get_all_Transcripts};
    # should only be 1 transcript but you never know;
    unless (scalar(@pseudo_trans) == 1){
      $self->warn("WARNING\nsingle exon gene has the wrong number of transcripts, this is likely to be a mistake\n");
    }
    @pseudo_trans = sort {$a->length <=> $b->length} @pseudo_trans;
    my $only_transcript_to_keep = pop  @pseudo_trans;
    $only_transcript_to_keep->translation(undef);
    foreach my $pseudo_transcript (@pseudo_trans) {
      $pseudo_transcript->translation(undef);
      $self->_remove_transcript_from_gene($pseudo_gene,$pseudo_transcript);
    }
    push @{$self->{'_pseudo_gene'}},$pseudo_gene;
  }
  return $self->{'_pseudo_gene'};
}

sub real_genes {
  my ($self, $real_gene) = @_;
  if ($real_gene){
    unless ($real_gene->isa("Bio::EnsEMBL::Gene")){
      $self->throw("real gene is not a Bio::EnsEMBL::Gene, it is a $real_gene");
    }
    push @{$self->{'_real_gene'}},$real_gene;
  }
  return $self->{'_real_gene'};
}
1;
