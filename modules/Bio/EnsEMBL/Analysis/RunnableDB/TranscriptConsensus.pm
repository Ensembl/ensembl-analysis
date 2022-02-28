=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus -

=head1 SYNOPSIS

  my $runnabledb =
    Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus->new(
                                                 -query    => $query,
                                                 -analysis => $analysis,
    )

=head1 DESCRIPTION

TranscriptConsensus is an extension of TranscriptCoalescer that combines
protein coding and EST transcripts to identify the best supported
transcript models.  The initial gene sets are clustered, then collapsed
into a non-redundant set of exons and introns which are assigned
scores according to the amount of supporting evidence they have.
The similarity genes have UTR added using the est transcripts and
the resulting models are assigned a score by summing the individual
exon and intron scores.  The transcripts are sorted by score and the
highest scoring models are made into gene objets and written to the
TranscriptConsensus database.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus;

use strict;
use warnings;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus;

use Bio::EnsEMBL::Utils::Exception
  qw(verbose throw warning info stack_trace_dump );
use Bio::EnsEMBL::Analysis::Tools::Utilities
  qw ( merge_config_details get_evidence_set convert_prediction_transcripts_to_genes );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);




=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus
  Function  : instatiates a TranscriptConsensus object and reads and checks the config
  file
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus
  Exceptions:
  Example   :

=cut


sub new {
  my ($class,@args) = @_ ;
  my $self = $class->SUPER::new(@args) ;
  $self->read_and_check_config($TRANSCRIPT_CONSENSUS_CONFIG_BY_LOGIC) ;

  return $self ;
  print STDERR " TC-db : new runnable created\n" ;
}


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus
  Function  : fetch input (gene-objects) out of different databases specified
              in the TranscriptConsensus config.
  Returntype: 1
  Exceptions: none
  Example   :

=cut



sub fetch_input {
  my ($self) = @_ ;

  $self->throw("No input id") unless defined($self->input_id) ;

  my $slice = $self->fetch_sequence($self->input_id, $self->db ) ;
  my $introns ;

  $self->query($slice);
  # fetch any RNASeq introns first
  if ( @{$self->RNASEQ_INTRON_DB} ) {
    foreach my $intron_db ( @{$self->RNASEQ_INTRON_DB}) {
      my $rnaseq_db = $self->get_dbadaptor($intron_db);
      my $rna_slice = $self->fetch_sequence($self->input_id, $rnaseq_db );

      my %intron_name_hash = %{$self->RNASEQ_INTRON_NAME};
      foreach my $ln ( @{$intron_name_hash{$intron_db}} ) {
        print STDERR "Found logic name ".$ln." for database ".$intron_db."\n";
        push @$introns , @{$rna_slice->get_all_DnaAlignFeatures($ln)};
      }
    }

    # make sure they are sorted
    @$introns = sort {$a->start <=> $b->start} @$introns;
    print STDERR "Found " . scalar(@$introns) ." unique intron features from the RNASeq db\n";
  }

  $self->{evidence_sets} = {
                'est'      => $self->EST_SETS,
                'simgw'    => $self->SIMGW_SETS,
                'abinitio' => $self->ABINITIO_SETS,
               };

  my %databases = %{$self->INPUT_GENES} ;

  # check if every biotype/logic_name belongs to an EVIDENCE set
  $self->_check_config() ;

  # now looping through all keys of %databases defined in the config
  my %biotypes_to_genes ;
  for my $db (keys %databases) {
    next unless (exists $databases{$db}) ;

    my $dba = $self->get_dbadaptor($db) ;
    $dba->dnadb($self->db) ;
    my $slice = $self->fetch_sequence($self->input_id, $dba );

    # organise genes into "EVIDENCE_SETS" depending on their underlying source
    # these sets are specified in the config file TranscriptConsensus.pm

    for my $biotype ( @{$databases{$db} }) {
      my $genes = $slice->get_all_Genes_by_type($biotype) ;
      # lazy load
      foreach my $gene (@$genes){
    foreach my $trans (@{$gene->get_all_Transcripts}){
      if ($trans->translateable_seq){
        $trans->translation->seq;
      }
      $trans->get_all_supporting_features;
    }
    foreach my $exon ( @{$gene->get_all_Exons} ){
      $exon->get_all_supporting_features;
    }
      }
      my ( $set ) = $self->get_evidence_set ( $biotype ) ;

      if (scalar( @{$genes} ) > 0 ) {
        # store genes of specific biotype in hash which holds them in an array
        # add specific information to exons / re-bless them

        my @multi_exon_genes ;
        for my $g (@$genes){
          my $single_exon_gene ;

          for my $t (@{$g->get_all_Transcripts}){
            throw ("gene has more than one transcript - only processing 1-gene-1-transcript-genes")
          if (@{$g->get_all_Transcripts}>1) ;

            bless $t, "Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended" ;
        $t->ev_set($set);
            my @all_exons = @{ $t->get_all_Exons } ;

            $single_exon_gene = 1 if (@all_exons == 1 ) ;

            map { bless $_,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended" } @all_exons ;

            for (my $int=0; $int < @all_exons ; $int++) {

              my $e = $all_exons[$int] ;
              $e->biotype($g->biotype) ;
              $e->ev_set($set) ;
              $e->transcript($t) ;
              $e->analysis($g->analysis) ;

              # setting pointer to previous exon
              $e->next_exon($all_exons[$int+1]);

              # setting pointer to previous exon
              if ($int == 0 ) {
                $e->prev_exon(0);
              } else {
                $e->prev_exon($all_exons[$int-1]);
              }
            }
          }
      # allow single exon genes
      push @multi_exon_genes, $g ;
    }
        push @{$biotypes_to_genes{$biotype} } , @multi_exon_genes ;
      }
    }

    # PREDICTION TRANSCRIPT PROCESSING
    #
    # get all PredictionTranscripts by their logic_name
    # these transcripts are converted to Genes with biotype eq logicname of analysis

    for my $logic_name_becomes_biotype ( @{ $self->AB_INITIO_LOGICNAMES }) {

      # get all PredictionTranscripts and convert them to Genes, set biotype to logic_name
      my $pt = $slice->get_all_PredictionTranscripts( $logic_name_becomes_biotype ) ;

      # get ev-set
      my $result_set_name  = $self->get_evidence_set( $logic_name_becomes_biotype ) ;
      my $ab_initio_genes = $self->convert_prediction_transcripts_to_genes( $pt,$logic_name_becomes_biotype,$result_set_name ) ;
      info ( scalar (@{$ab_initio_genes}) . "\t$logic_name_becomes_biotype-genes\n") ;
      if ( scalar(@$ab_initio_genes)>0){
    push @{$biotypes_to_genes{$logic_name_becomes_biotype} } , @{$ab_initio_genes} ;
      }
    }
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus->new
    (
     -query                => $self->query,
     -analysis             => $self->analysis,
     -all_genes            => \%biotypes_to_genes , # ref to $hash{biotype_of_gene} = @all_genes_of_this_biotype
     -evidence_sets        => $self->{evidence_sets} ,
     -verbose              => $self->VERBOSE,
     -filter_singletons    => $self->FILTER_SINGLETONS,
     -filter_non_consensus => $self->FILTER_NON_CONSENSUS,
     -filter_ests          => $self->FILTER_ESTS,
     -add_utr              => $self->ADD_UTR,
     -min_consensus        => $self->MIN_CONSENSUS,
     -utr_penalty          => $self->UTR_PENALTY,
     -end_exon_penalty     => $self->END_EXON_PENALTY,
     -est_overlap_penalty  => $self->EST_OVERLAP_PENALTY,
     -short_intron_penalty => $self->SHORT_INTRON_PENALTY,
     -short_exon_penalty   => $self->SHORT_EXON_PENALTY,
     -good_percent         => $self->GOOD_PERCENT,
     -good_biotype         => $self->GOOD_BIOTYPE,
     -small_biotype        => $self->SMALL_BIOTYPE,
     -bad_biotype          => $self->BAD_BIOTYPE,
     -rnaseq_introns       => $introns,
    );
  $self->runnable($runnable);

  return 1;
}


sub write_output {
  my ($self) = @_;
  my $out_dba = $self->get_dbadaptor($self->OUTPUT_DATABASE) ;
  my $gene_a = $out_dba->get_GeneAdaptor() ;
  info ("trying to write output") ;

  foreach my $gene (@{$self->output}) {
    info("STORED GENE $gene" ) ;
    eval {
      $gene_a->store($gene) ;
    } ;
    if($@) {
      $self->throw("Failed to write gene ". $@) ;
    }
  }
  return ;
}




=head2  _check_config

   Arg[0] : Hash-reference to $TRANSCRIPT_CONSENSUS_BY_LOGIC - Hash out of TranscriptConsensus.pm
   Arg[1] : href with evidence sets
   Func   : Checks if every logic_name / biotype of gene belongs to an
         Evidence set
   Return : 1 if all is ok

=cut

sub _check_config {
  my ( $self ) = @_ ;
  my $ev_sets = $self->{evidence_sets} ;

  # check that every biotype / logic_name in TranscriptConsensus.pm
  # belongs to an ev-set
  my %database_definition ;
  for my $db_class (keys %{$self->INPUT_GENES}) {
    map $database_definition{$_}=(),@{${$self->INPUT_GENES}{$db_class}} ;
  }

  my %ev_set_definitions;

  for my $set_name (keys %{ $ev_sets } ) {
    map $ev_set_definitions{$_}=(), @{$$ev_sets{$set_name}} ;
  }


  # check that every biotype / logicname mentioned in the evidence_sets
  # has a database where it's fetched from

  for my $ln ( keys %ev_set_definitions ) {
    throw ("\n\tError in config-file Analysis/Config/GeneBuild/TranscriptConsensus.pm\n ".
           "\tType $ln has an evidence set but there's NO database for this type defined in Config".
             " TranscriptConsensus.pm : $ln\n")
            unless ( exists $database_definition{$ln} ) ;

  }

  # check if every ev-set has an entry in the db-section as well

  for my $ln ( keys %database_definition) {
    throw ("\n\tError in config-file Analysis/Config/GeneBuild/TranscriptConsensus.pm\n".
           "\tType $ln has an entry in a database-section but there's NO evidence-set defined in Config TranscriptConsensus.pm for this type : $ln \n" )
            unless ( exists $ev_set_definitions{$ln} ) ;
  }
  return 1 ;
}



=head2 INPUT_GENES

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype:
  Exceptions:
  Example   :

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus


sub INPUT_GENES {
  my ($self, $arg) = @_ ;
  if (defined $arg) {
    $self->{'INPUT_GENES'} = $arg ;
  }
  return $self->{'INPUT_GENES'} ;
}

sub ABINITIO_SETS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'ABINITIO_SETS'} = $arg ;
  }
  return $self->{'ABINITIO_SETS'} ;
}

sub SIMGW_SETS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'SIMGW_SETS'} = $arg ;
  }
  return $self->{'SIMGW_SETS'} ;
}

sub EST_SETS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'EST_SETS'} = $arg ;
  }
  return $self->{'EST_SETS'} ;
}

sub AB_INITIO_LOGICNAMES {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'AB_INITIO_LOGICNAMES'} = $arg ;
  }
  return $self->{'AB_INITIO_LOGICNAMES'} ;
}

sub VERBOSE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'VERBOSE'} = $arg ;
  }
  return $self->{'VERBOSE'} ;
}

sub OUTPUT_DATABASE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'OUTPUT_DATABASE'} = $arg ;
  }
  return $self->{'OUTPUT_DATABASE'} ;
}

sub FILTER_SINGLETONS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'FILTER_SINGLETONS'} = $arg ;
  }
  return $self->{'FILTER_SINGLETONS'} ;
}

sub FILTER_NON_CONSENSUS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'FILTER_NON_CONSENSUS'} = $arg ;
  }
  return $self->{'FILTER_NON_CONSENSUS'} ;
}

sub FILTER_ESTS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'FILTER_ESTS'} = $arg ;
  }
  return $self->{'FILTER_ESTS'} ;
}

sub ADD_UTR {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'ADD_UTR'} = $arg ;
  }
  return $self->{'ADD_UTR'} ;
}

sub MIN_CONSENSUS {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'MIN_CONSENSUS'} = $arg ;
  }
  return $self->{'MIN_CONSENSUS'} ;
}

sub UTR_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'UTR_PENALTY'} = $arg ;
  }
  return $self->{'UTR_PENALTY'} ;
}

sub END_EXON_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'END_EXON_PENALTY'} = $arg ;
  }
  return $self->{'END_EXON_PENALTY'} ;
}

sub EST_OVERLAP_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'EST_OVERLAP_PENALTY'} = $arg ;
  }
  return $self->{'EST_OVERLAP_PENALTY'} ;
}

sub SHORT_INTRON_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'SHORT_INTRON_PENALTY'} = $arg ;
  }
  return $self->{'SHORT_INTRON_PENALTY'} ;
}

sub SHORT_EXON_PENALTY {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'SHORT_EXON_PENALTY'} = $arg ;
  }
  return $self->{'SHORT_EXON_PENALTY'} ;
}

sub GOOD_PERCENT {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'GOOD_PERCENT'} = $arg ;
  }
  return $self->{'GOOD_PERCENT'} ;
}

sub GOOD_BIOTYPE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'GOOD_BIOTYPE'} = $arg ;
  }
  return $self->{'GOOD_BIOTYPE'} ;
}

sub BAD_BIOTYPE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'BAD_BIOTYPE'} = $arg ;
  }
  return $self->{'BAD_BIOTYPE'} ;
}

sub SMALL_BIOTYPE {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'SMALL_BIOTYPE'} = $arg ;
  }
  return $self->{'SMALL_BIOTYPE'} ;
}

sub RNASEQ_INTRON_DB {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'RNASEQ_INTRON_DB'} = $arg ;
  }
  return $self->{'RNASEQ_INTRON_DB'} ;
}

sub RNASEQ_INTRON_NAME {
  my ($self, $arg) = @_ ;
  if(defined $arg) {
    $self->{'RNASEQ_INTRON_NAME'} = $arg ;
  }
  return $self->{'RNASEQ_INTRON_NAME'} ;
}



1;
