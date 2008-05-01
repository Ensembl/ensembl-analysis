=pod

=head1 NAME Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus->new
     (
     -query            => $query,
     -analysis        => $analysis,
     )

=head1 DESCRIPTION

TranscriptConsensus is an extension of TranscriptCoalescer that combines protein coding and est
transcripts to identify the best supported transcript models.
The initial gene sets are clustered, then collapsed into a non-redundant set of exons
and introns which are assigned scores according to the amount of supporting evidence they have.
The similarity genes have UTR added using the est transcripts and the resulting models are assigned
a score by summing the individual exon and intron scores.
The transcripts are sorted by score and the highest scoring models are made into
gene objets and written to the TranscriptCoalescer database.

=head1 CONTACT

Post questions to the EnsEMBL developer list: <ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus;

use strict;
use warnings;
use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes; 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptCoalescer;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info stack_trace_dump );
use Bio::EnsEMBL::Analysis::RunnableDB::TranscriptCoalescer;
use Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::TranscriptCoalescer);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::TranscriptConsensus
  Function  : fetch input (gene-objects) out of different databases specified
              in the TranscriptCoalescer config.
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;

  $self->throw("No input id") unless defined($self->input_id);

  my $slice = $self->fetch_sequence($self->input_id, $self->db );
  my $solexa_slice;
  my (%est_genes, %simgw_genes, %pred_transcripts )  ;  

  $self->query($slice);


  # the config for TranscriptCoalescer.pm is distributed over 3 files 
  #    
  # - Databases.pm 
  # - Exonerate2Genes.pm
  # - TranscriptCoalescer.pm 
  #
  # we merge the configs such that we have all in one hash $database 
  #
  
  my %databases = %{$self->merge_config_details( $DATABASES, $EXONERATE_CONFIG_BY_LOGIC, $TRANSCRIPT_COALESCER_DB_CONFIG)} ;   

  # check if every biotype/logic_name belongs to an EVIDENCE set
  _check_config(\%databases, $TRANSCRIPT_COALESCER_DB_CONFIG, $self->{evidence_sets} ) ; 
 
  # now looping through all keys of %databases definend in the (merged) configs 

  my %biotypes_to_genes ; 
  
  for my $database_class (keys %databases ) { 
  
    # skip all databases out of Exonerate2Genes.pm & Databases.pm which
    # have no configuration in TranscriptCoalescer.pm
    # Unless its the solexa db in the trascript consensus config

    if ($SOLEXA && exists $$TRANSCRIPT_CONSENSUS_DB_CONFIG{$database_class}) {
      my $dba = $self->get_connection($database_class,\%databases);

      # attach refdb as dnadb  
      $dba->dnadb($self->db) ; 
      $solexa_slice = $self->fetch_sequence($self->input_id, $dba );
     }

    next unless (exists $$TRANSCRIPT_COALESCER_DB_CONFIG{$database_class}) ;  


    my $dba = $self->get_connection($database_class,\%databases);
    # attach refdb as dnadb  
    $dba->dnadb($self->db) ; 
    
    my $slice = $self->fetch_sequence($self->input_id, $dba );
    # 
    # organise genes into "EVIDENCE_SETS" depending on their underlying source 
    # these sets are specified in the config file TranscriptCoalescer.pm 
    #
    #   
    for my $biotype ( @{$databases{$database_class}{BIOTYPES} }) {
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
      my ( $set ) = $self->_get_evidence_set ( $biotype ) ;

      my $db_used = "${$databases{$database_class}}{db}{-dbname}\@".
        "${$databases{$database_class}}{db}{-host}::".
          "${$databases{$database_class}}{db}{-port}"; 
      info( "$db_used :" ); 
      info(scalar (@{$genes}) . "\t$biotype-genes fetched \n") ;



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
  # 
    # PREDICTION TRANSCRIPT PROCESSING 
    #
    # get all PredictionTranscripts by their logic_name 
    # these transcripts are converted to Genes with biotype eq logicname of analysis 
    #  

    for my $logic_name_becomes_biotype ( @{$databases{$database_class}{AB_INITIO_LOGICNAMES} }) { 

     # get all PredictionTranscripts and convert them to Genes, set biotype to logic_name  
     my $pt = $slice->get_all_PredictionTranscripts( $logic_name_becomes_biotype ) ;  
     
     # get ev-set 
     my $result_set_name  = $self->_get_evidence_set( $logic_name_becomes_biotype ) ; 

     my $ab_initio_genes = convert_prediction_transcripts_to_genes(
                            $pt,$logic_name_becomes_biotype,$result_set_name ) ; 

     info ( scalar (@{$ab_initio_genes}) . "\t$logic_name_becomes_biotype-genes\n") ; 

     if ( scalar(@$ab_initio_genes)>0){ 
       push @{$biotypes_to_genes{$logic_name_becomes_biotype} } , @{$ab_initio_genes} ;
     }
   }
  } 
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus->new
    (
     -query    => $self->query,
     -solexa   => $solexa_slice,
     -analysis => $self->analysis,
     -all_genes   => \%biotypes_to_genes , # ref to $hash{biotype_of_gene} = @all_genes_of_this_biotype
     -evidence_sets   => $self->{evidence_sets} , 
     -dnadb => $self->db , 
     -utils_verbosity => $self->{utils_verbosity}, 
    );
  $runnable->solexa_slice($solexa_slice) if ($SOLEXA && $solexa_slice);
  $self->runnable($runnable);



  return 1;
}



=head2  _check_config

   Arg[0] : Hash-reference to $TRANSCRIPT_COALESCER_DB_CONFIG-Hash out of TranscriptCoalescer.pm 
   Arg[1] : href with evidence sets 
   Func   : Checks if every logic_name / biotype of gene belongs to an
         Evidence set 
   Return : 1 if all is ok
 
=cut  

sub _check_config{ 
  my ( $databases, $TRANSCRIPT_COALESCER_DB_CONFIG, $ev_sets) = @_ ; 

  # check TranscriptCoalescer.pm:  
  for my $db_class (keys %{$TRANSCRIPT_COALESCER_DB_CONFIG}) { 
    # check that each db_class (like REFERENCE_DB ...) in TranscriptCoalescer.pm has 
    # also a database in Databases.pm / Exonerate2Genes.pm 
    throw ( "\n\tConfig-Error: There is configuration for \"$db_class\" in Analysis/Config/GeneBuild/TranscriptCoalescer.pm " .
     "but no Database defined in Exonerate2Genes.pm or Databases.pm") 
     unless exists $$databases{$db_class} ; 
  }


  # check that every biotype / logic_name in TranscriptCoalescer.pm 
  # belongs to an ev-set 
  
  my %database_definition ; 
  for my $db_class (keys %{$TRANSCRIPT_COALESCER_DB_CONFIG}) {
    for my $key ( keys %{ $$TRANSCRIPT_COALESCER_DB_CONFIG{$db_class} } ) { 
       map $database_definition{$_}=(),@{ ${$$TRANSCRIPT_COALESCER_DB_CONFIG{$db_class}}{$key}};  
    }
  }
   
  my %ev_set_definitions; 

  for my $set_name (keys %{ $ev_sets} ) { 
    map $ev_set_definitions{$_}=(), @{$$ev_sets{$set_name}} ; 
  }

 
  # check that every biotype / logicname mentioned in the evidence_sets 
  # has a database where it's fetched from 
  
  for my $ln ( keys %ev_set_definitions ) { 
    throw ("\n\tError in config-file Analysis/Config/GeneBuild/TranscriptCoalescer.pm\n ". 
           "\tType $ln has an evidence set but there's NO database for this type defined in Config". 
             " TranscriptCoalescer.pm : $ln\n")
            unless ( exists $database_definition{$ln} ) ; 
     
  } 

  # check if every ev-set has an entry in the db-section as well 
  #
  
  for my $ln ( keys %database_definition) { 
    throw ("\n\tError in config-file Analysis/Config/GeneBuild/TranscriptCoalescer.pm\n". 
           "\tType $ln has an entry in a database-section but there's NO evidence-set defined in Config TranscriptCoalescer.pm for this type : $ln \n" )
            unless ( exists $ev_set_definitions{$ln} ) ; 
     
  } 
  return 1 ;
}


sub get_connection {
  my ($self,$database_class,$d) = @_;
  my %databases = %$d;
#  print "Connection  " . $database_class ."\n";
  unless ( $self->{_db_conn}{$database_class} ) {
    $self->{_db_conn}{$database_class} = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ ${$databases{$database_class}}{db}} ) ; 
  }
  return $self->{_db_conn}{$database_class};
}
1;
