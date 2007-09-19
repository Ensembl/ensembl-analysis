package Bio::EnsEMBL::Analysis::RunnableDB::ExonerateForGenewise;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(Transcript_info set_stop_codon set_start_codon attach_Slice_to_Transcript evidence_coverage low_complexity_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(contains_internal_stops);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(Gene_info print_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise qw(GENEWISE_CONFIG_BY_LOGIC);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name write_seqfile);

@ISA = qw (
           Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
           );



sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($GENEWISE_CONFIG_BY_LOGIC);
  $self->extra_sanity_check();
 
  return $self;
}

sub create_bmg_runnables{
  my ($self, $hits_hash, $seqfetcher, $analysis) = @_;
  throw("RunnableDB::BlastMiniGenewise Can't create BlastMiniGenewise ".
        "runnables without a hash of hits") unless($hits_hash);
  $seqfetcher = $self->seqfetcher if(!$seqfetcher);
  $analysis = $self->analysis if(!$analysis);
  my %exonerate_parameters = %{$self->EXONERATE_PARAMETERS};
  my %params_hash = %{$self->parameters_hash};
  if($self->LIMIT_TO_FEATURE_RANGE == 0){
    my $genomic_file = write_seqfile($self->query, 
                                     create_file_name("genomic", "fa",
                                                      "/tmp/"));
    $self->genomic_file($genomic_file);
    my @ids = sort keys(%$hits_hash);
    my $seqs = $self->get_protein_sequences(\@ids);
    my $peptide_file = write_seqfile($seqs, 
                                     create_file_name("peptide", "fa",
                                                      "/tmp/"));
    $self->peptide_file($peptide_file);
    my $r = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new
      (
       -analysis => $self->analysis,
       -target_file => $genomic_file,
       -query_type => 'protein',
       -query_file => $peptide_file,
       -program => $self->analysis->program_file,
       %exonerate_parameters,
       %params_hash,
      );
    $self->runnable($r);
  }else{
    foreach my $id(keys(%$hits_hash)){
      my $features = $hits_hash->{$id};
      foreach my $feature(@$features){
        print "MAKING RUNNABLE with ".$id." protein\n";
        #This should be made a non overlapping range
        my $start = $feature->seq_region_start - $self->FEATURE_RANGE_PADDING;
        my $end = $feature->seq_region_end + $self->FEATURE_RANGE_PADDING;
        $start = 1 if($start < 1);
        $end = $self->query->seq_region_length 
          if($end > $self->query->seq_region_length);
        my $name = ($self->query->coord_system->name.":".
                    $self->query->coord_system->version.":".
                    $self->query->seq_region_name.":".
                    $start.":".$end.":".
                    $self->query->strand);
        my $db = $self->get_dbadaptor($self->GENE_SOURCE_DB);
        my $query = $self->fetch_sequence($name, $db, $self->REPEATMASKING);
        logger_info("Creating ExonerateTranscript Runnable over a limited range with "
                    .$id." and ".$seqfetcher." to run on ".$query->name);
        my $genomic_file = write_seqfile($query, 
                                            create_file_name("genomic", "fa",
                                                             "/tmp/"));
        $self->genomic_file($genomic_file);
        my $seqs = $self->get_protein_sequences([$id]);
        my $peptide_file = write_seqfile($seqs, 
                                            create_file_name("peptide", "fa",
                                                             "/tmp/"));
        $self->peptide_file($peptide_file);
        my $r = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new
          (
           -analysis => $self->analysis,
           -target_file => $genomic_file,
           -query_type => 'protein',
           -query_file => $peptide_file,
           %exonerate_parameters,
           %params_hash,           
          );
        $self->runnable($r);
      }
    }
  }
}


sub process_transcripts{
  my ($self, $transcripts) = @_;
  my @complete_transcripts;
  foreach my $transcript(@$transcripts){
    my $query_name = $transcript->start_Exon->seqname;
    my $unmasked_slice = $self->fetch_sequence($query_name, 
                                               $self->output_db);
    attach_Slice_to_Transcript($transcript, $unmasked_slice);
    my $start_t = set_start_codon($transcript);
    my $end_t = set_stop_codon($start_t);
    if(contains_internal_stops($transcript)){
      warning(Transcript_info($transcript)." contains internal stops will be thrown ".
              "This will be thrown away");
    }
    push(@complete_transcripts, $end_t);
  }
  foreach my $file($self->genomic_file, @{$self->peptide_file}){
    #unlink $file;
  }
  return \@complete_transcripts
}

=head2 filter_genes

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : arrayref of Genes
  Function  : filters the gene sets using the defined filter
  Returntype: n/a, stores the genes in the output array
  Exceptions: 
  Example   : 

=cut


sub filter_genes{
  my ($self, $genes) = @_;
  my $output;
  my $rejected_set;
  my $filtered_set;
  if($self->filter_object){
    ($filtered_set, $rejected_set) = 
      $self->filter_object->filter_genes($genes);
    #$self->rejected_set($rejected_set);
    foreach my $gene(@{$filtered_set}){
      my $transcripts = $gene->get_all_Transcripts;
      my $transcript = $transcripts->[0];
      my $coverage = evidence_coverage($transcript);
      if(low_complexity_less_than_maximum($transcript, $coverage)){
        push(@$output, $gene);
      }else{
        $gene->description("low_complexity is greater than evidence coverage");
        push(@$rejected_set, $gene);
      }
    }
  }else{
    $self->output($genes);
  }
  $self->output($output);
  $self->rejected_set($rejected_set);
}

sub genomic_file{
  my ($self) = @_;
  if(@_){
    $self->{genomic_file} = shift;
  }
  return $self->{genomic_file};
}

sub peptide_file{
  my ($self) = @_;
  if(@_){
    my $file = shift;
    push(@{$self->{peptide_file}}, $file);
  }
  return $self->{peptide_file};
}


sub get_protein_sequences{
  my ($self, $ids) = @_;
  my $seqfetcher = $self->seqfetcher;
  my $seqs;
  my $seq;
  foreach my $id(@$ids){
    eval {
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };
    if ($@) {
      throw("Problem fetching sequence for [$id]: [$@]\n");
    }
    if(!$seq){
      throw("Can't fetch the sequence $id with ".$seqfetcher." - it's not in the index\n"); 
    }
    push(@$seqs, $seq);
  }
  return $seqs;
}


sub extra_sanity_check{
  my ($self) = @_;
  throw("ExonerateAlignFeatures needs an exonerate parameters hash defined")
    if(!$self->EXONERATE_PARAMETERS);
}


1;
