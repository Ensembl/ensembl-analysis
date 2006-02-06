# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::TranscriptCoalescer
#
# Copyright (c) 2005 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::TranscriptCoalescer

=head1 SYNOPSIS

  my $cond = Bio::EnsEMBL::Analysis::RunnableDB::TranscriptCoalescer->new
             (
               -analysis => $analysis,
               -db => $db,
               -input_id => 'chromosome:JGI2:chr_04q:70000:80000:1'
              );
  $cond->fetch_input;
  $cond->run;
  $cond->write_output; 

=head1 DESCRIPTION

  This module acts as an intermediate between the runnable and the
  core database. It reads configuration and uses information from the analysis
  object to setup the runnable and then write the results back to the 
  database specified in the config file.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::TranscriptCoalescer;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;

use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes; 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptCoalescer;

use Bio::EnsEMBL::Analysis::Runnable::TranscriptCoalescer;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended ;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info stack_trace_dump );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->verbose(0) ;  
  my @est_biotypes = @{ $EST_SETS } ; 
  my @simgw_biotypes = @{ $SIMGW_SETS } ; 
  my @abinitio_logicnames = @{ $ABINITIO_SETS } ; 
  #
  #  dont' change the key-names of the $tmp hash - they are 
  #  used in the runnable to distinguish between the differnt sets !
  #
  $self->{evidence_sets}=  {
                             'est' => \@est_biotypes , 
                             'simgw' => \@simgw_biotypes , 
                             'abinitio' => \@abinitio_logicnames , 
                            } ; 
  # verbosity-level can be controlled by test_RunnableDB -utils_verbosity INFO 
  # see Bio::EnsEMBL::Utils::Exception
 
  &verbose($self->{utils_verbosity})  ; 
  return $self;
}



=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::TranscriptCoalescer
  Function  : fetch input (gene-objects) out of different databases specified 
              in the TranscriptCoalescer config 
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;

  $self->throw("No input id") unless defined($self->input_id);

  my $slice = $self->fetch_sequence($self->input_id, $self->db );

  my (%est_genes, %simgw_genes, %pred_transcripts )  ;  

  $self->query($slice);


  # the config for TranscriptCoalescer.pm is distributed over 3 files 
  # - Databases.pm 
  # - Exonerate2Genes.pm
  # - TranscriptCoalescer.pm 
  #
  # we merge the configs such that we have all in one hash $database 
  my %databases = %{merge_database_hashes ( $DATABASES, $DATA, undef ) }  ;
  %databases = %{merge_database_hashes ( $EXONERATE_CONFIG_BY_LOGIC, $DATA, \%databases ) }  ; 
 
  # check if every biotype/logic_name belongs to an EVIDENCE set
  check_config(\%databases, $DATA, $self->{evidence_sets} ) ; 
 
  # now looping through all keys of %databases definend in the (merged) configs 

  my %biotypes_to_genes ; 
  
  for my $database_class (keys %databases ) { 

    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ ${$databases{$database_class}}{db}} ) ; 
    my $slice = $self->fetch_sequence($self->input_id, $dba );
    # 
    # organise genes into "EVIDENCE_SETS" depending on their underlying source 
    # these sets are specified in the config file TranscriptCoalescer.pm 
    #
    #   
    for my $biotype ( @{$databases{$database_class}{BIOTYPES} }) { 
      my $genes = $slice->get_all_Genes_by_type($biotype) ; 

      my ( $set ) = $self->get_evidence_set ( $biotype ) ; 

      my $db_used = "${$databases{$database_class}}{db}{-dbname}\@".
        "${$databases{$database_class}}{db}{-host}::".
          "${$databases{$database_class}}{db}{-port}"; 
      info( "$db_used :" ); 
      info(scalar (@{$genes}) . "\t$biotype-genes fetched \n") ; 


      if (scalar( @{$genes} ) > 0 ) { 
        # store genes of specific biotype in hash which holds them in an array 
        # add specific information to exons / re-bless them 
        
        for my $g (@$genes){ 

          for my $t (@{$g->get_all_Transcripts}){
            bless $t, "Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended" ; 

            my @all_exons = @{ $t->get_all_Exons } ; 
            map { bless $_,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended" } @all_exons ; 

            for (my $int=0; $int < @all_exons ; $int++) {

              my $e = $all_exons[$int] ; 
              $e->biotype($g->biotype) ;
              $e->ev_set($set) ;
              $e->transcript($t) ; 
              $e->number_exons($#all_exons) ; 

              $e->analysis($g->analysis) ; 

              $e->source_transcript ($t) ; # this refers to original transcript the exon is from  
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
        }  
        push @{$biotypes_to_genes{$biotype} } , @{$genes} ;
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
     my $result_set_name  = $self->get_evidence_set( $logic_name_becomes_biotype ) ; 

     my $ab_initio_genes = convert_prediction_transcripts_to_genes(
                            $pt,$logic_name_becomes_biotype,$result_set_name ) ; 

     info ( scalar (@{$ab_initio_genes}) . "\t$logic_name_becomes_biotype-genes\n") ; 

     if ( scalar(@$ab_initio_genes)>0){ 
       push @{$biotypes_to_genes{$logic_name_becomes_biotype} } , @{$ab_initio_genes} ;
     }
   }
  } 

  #
  # remove redundant transcripts which have all same exons to speed computation up 
  # 
  my $non_redundant_genes = remove_redundant_genes(\%biotypes_to_genes) ;  


  my $runnable = Bio::EnsEMBL::Analysis::Runnable::TranscriptCoalescer->new
    (
     -query    => $self->query,
     -analysis => $self->analysis,
     -all_genes   => $non_redundant_genes , # ref to $hash{biotype_of_gene} = @all_genes_of_this_biotype
     -new_biotype     =>  $NEW_BIOTYPE,
     -evidence_sets   => $self->{evidence_sets} , 
    );
  $self->runnable($runnable);

  return 1;
}


sub merge_database_hashes {
  my ($databases, $data, $new_data  ) = @_ ; 

  my %db = %{$databases  } ; 
  my %tc = %{ $data } ; 

  my %result ;   
  
  %result = %$new_data if $new_data ; ; 

  # COALESCER_DB is NOT merged !
  
  for my $db_class ( keys %db ) {  
    #print "DB_CLASS $db_class\n" ; 
    if ( exists $tc{$db_class} ) {   
      # there is the same logic name in TranscriptCoalescer.pm and a config-file 
      # (Exonerate2Genes.pm or Databases.pm)
      
      # get database-config out of Exonerate2Genes.pm 
      if ( exists  ${$db{$db_class}}{OUTDB} ) { 
        # TranscriptCoaleser.pm and Exonerate2Genes share same db_class key 

         unless ($db_class=~/DEFAULT/){
         # we are not processing any entry in the DEFAULT section of Exonerate2Genes.pm
        
          # get OUTDB anonymous hash and hand it over to new config   
          $result{$db_class}{db} = ${$db{$db_class}}{OUTDB} ;
                   
        } 
      }
     $result{$db_class}{db}=$db{$db_class} unless $result{$db_class}{db}; 
     # there is some configuration (biotype/logic_name) for this db_class
     # in Databases.pm as well as in TrancriptCoalescer.pm 
     
     my %tc_conf = %{ $tc{$db_class }} ; 

     for my $tc_name ( keys %tc_conf ) { 
       # tc_name is either AB_INITIO_LOGICNAME or BIOTYPE
       $result{$db_class}{$tc_name}=$tc_conf{$tc_name} ;  
     }
    }
  } 
  return \%result ; 
} 


sub remove_redundant_genes {
  my ($biotype2genes ) = @_ ; 
  info ("removing redundant Transcripts........\n\n" ); 
  
  my %biotype2genes = %{$biotype2genes} ; 
  my %result ;
  my %tmp ;  
  for my $biotype ( keys %biotype2genes) {
    my ( %non_redundant_genes) ;  
    my @genes = @{$biotype2genes{$biotype}} ; 
    push @{$tmp{$biotype}}, scalar(@genes) ; 

    info ( "RunnableDB:  $biotype: " . scalar( @genes ) . "\n" ); 

    my $red = 0 ;      
    for my $g (@genes ) {
      $red++;
      my @all_tr = @{ $g->get_all_Transcripts } ; 

      throw ("Only processing genes with 1gene:1transcript relation")  if (@all_tr > 1 ) ; 

      my $tr_hashkey = "" ;  
      for my $t ( @all_tr ) { 
        for my $e ( @{ $t->get_all_Exons} ) { 
          $e->end_phase(0) unless ($e->end_phase); 
          $tr_hashkey .=$e->hashkey; 
        }
        $non_redundant_genes{$tr_hashkey}=$g; 
      }
    }
    
    for my $key (keys %non_redundant_genes) {
      push @{ $result{$biotype}},  $non_redundant_genes{$key} ;
    }
    push @{$tmp{$biotype}}, scalar(@{$result{$biotype}}) ;
  }
  for (keys %tmp) { 
    info("having ${$tmp{$_}}[1] non-redundant genes ( out of ${ $tmp{$_} }[0] of biotype $_)\n") ; 
  }
  return \%result ; 
}



=head2  check_config

   Arg[0] : Hash-reference to $DATA-Hash out of TranscriptCoalescer.pm 
   Arg[1] : href with evidence sets 
   Func   : Checks if every logic_name / biotype of gene belongs to an
         Evidence set 
   Return : 1 if all is ok
 
=cut  

sub check_config{ 
  my ( $databases, $DATA , $ev_sets) = @_ ; 

  # check TranscriptCoalescer.pm:  
  for my $db_class (keys %{$DATA}) { 
    # check that each db_class (like REFERENCE_DB ...) in TranscriptCoalescer.pm has 
    # also a database in Databases.pm / Exonerate2Genes.pm 
    throw ( "\n\tConfig-Error: There is configuration for \"$db_class\" in GeneBuild/TranscriptCoalescer.pm " .
     "but no Database defined in Exonerate2Genes.pm or Databases.pm") 
     unless exists $$databases{$db_class} ; 
  }


  # check that every biotype / logic_name in TranscriptCoalescer.pm 
  # belongs to an ev-set 
  
  my %database_definition ; 
  for my $db_class (keys %{$DATA}) {
    for my $key ( keys %{ $$DATA{$db_class} } ) { 
       map $database_definition{$_}=(),@{ ${$$DATA{$db_class}}{$key}};  
       #push @database_definition, @{ ${$$DATA{$db_class}}{$key}}; 
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



=head2 convert_prediction_transcripts_to_genes 

  Arg [0]   : reference to an array of Bio::EnsEMBL::PredictionTranscript-objects
  Arg [1]   : String describing the logic_name  
  Arg [2]   : String describing name of the evidence-set 
  Function  : Creates Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended-Objects from a 
              set of Bio::EnsEMBL::PredictionTranscript-Objects and adds further information to these Transcripts.
              The PredictionExons are re-blessed to Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended-objects
  Returntype: Ref. to arrray of  Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended-objects
  Example   : 

=cut


sub  convert_prediction_transcripts_to_genes {
  my ($pt,$logic_name_becomes_biotype,$ev_set_name ) = @_ ; 
  my @new_genes ;  
  for my $pt (@$pt) { 
    # conversion 
    my $gene_from_pt = Bio::EnsEMBL::Gene->new( 
                       -start => $pt->start , 
                       -end => $pt->end , 
                       -strand => $pt->strand ,  
                       -slice =>$pt->slice ,  
                       -biotype => $logic_name_becomes_biotype,
                       -analysis=>$pt->analysis, 
                       ) ;  

    my $new_tr = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended->new(
                    -BIOTYPE => $logic_name_becomes_biotype ,  
                    -ANALYSIS => $pt->analysis , 
                 ) ; 
 
    my @pt_exons  = @{$pt->get_all_Exons} ; 

    for (my $i=0 ; $i<scalar(@pt_exons) ; $i++) { 
      
      # converting Bio::EnsEMBL::PredictionExon into ExonExtened (ISA Bio::EnsEMBL::Exon)  
      my $pte =$pt_exons[$i] ;  
      bless $pte,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended" ; 
      $pte->biotype($logic_name_becomes_biotype) ;  
      $pte->ev_set($ev_set_name) ;  
      $pte->end_phase(0); 
      $pte->phase(0); 
      $pte->next_exon($pt_exons[$i+1]) ; 
      $pte->prev_exon($pt_exons[$i-1]) ; 
      $pte->transcript($new_tr) ; 
      $pte->analysis($pt->analysis) ; 
    } ;
     
    #
    # Extending the Bio::EnsEMBL::Transcript object by ev_set methods 
    #
    for (@pt_exons) { 
      $new_tr->add_Exon($_); 
    }   

    $gene_from_pt->add_Transcript($new_tr) ; 

    push @new_genes , $gene_from_pt ; 
  }
  return \@new_genes ;  
}


=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::TranscriptCoalescer
  Function  : go through the runnables, run each one and check for errors
  and push the output on to the output array
  Returntype: 1;
  Exceptions: throws if blast run fails
  Example   : 

=cut



sub run {
  my ($self) = @_;

  my @runnables = @{$self->runnable};
  foreach my $runnable(@runnables){
    eval{
      $runnable->run;
    };
    #
    # checking errors from Runnable 
    #
    if(my $err = $@){
      chomp $err;
      print $err ;  
      # only match '"ABC_DEFGH"' and not all possible throws
      if ($err =~ /^\"([A-Z_]{1,40})\"$/i) {
        my $code = $1;
        if ($code ne 'VOID') {
          $self->failing_job_status($1);          
          throw("TranscriptCoalescer::run failed $@");
        }
      }
    } else { 
      print $@; 
      print STDERR "RunnableDB $self run was OK !\n\n " ;
    }  
   $self->output($runnable->output);
   #print $runnable->output."\n" ; 
  }
  1;
}






=head2 write_output

  Arg [1]   : $sself
  Function  : get appropriate adaptors and write output to database
       after validating and attaching sequence and analysis objects
  Returntype: undef
  Exceptions: throws if the store fails
  Example   : 

=cut


sub write_output{
  my ($self) = @_;

  my $out_dba = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{$$DATABASES{COALESCER_DB}} ) ; 

  my $gene_a = $out_dba->get_GeneAdaptor() ; 
  print "trying to write output\n" ;  
  foreach my $gene (@{$self->output}){
     print "storing $gene\n" ; 
    $gene_a->store($gene) ; 
  }  
  return ;
}


=head2 get_evidence_set ($logic_name_or_biotype)

  Name     : get_evidence_set( $logic_name_or_biotype )
  Arg      : String 
  Func     : returns the name of the evidence_set of a genee / PredictionTranscript 
  Returnval: String describing evidence_set_name

=cut 

sub get_evidence_set {
  my ($self, $logic_name_or_biotype) = @_ ; 

  my %ev_sets = %{ $self->{evidence_sets} } ;
  my $result_set_name ; 
  for my $set_name (keys %ev_sets){
    my @logic_names = @{$ev_sets{$set_name}} ; 
    for my $ln (@logic_names ) { 
       if ($ln eq $logic_name_or_biotype) { 
         $result_set_name = $set_name ;
      }
    }
  }
  return $result_set_name ; 
}     


1;
