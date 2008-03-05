# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::BestTargetted

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $cond = Bio::EnsEMBL::Analysis::RunnableDB::BestTargetted
    ->new (-db         => $pipelinedb,
           -input_id   => $input_id
           -analysis   => $analysis );
  $cond->fetch_input();
  $cond->run();
  $cond->write_output(); 




=head1 DESCRIPTION

  This module acts as an intermediate between the runnable and the
  core database. It reads configuration and uses information from the analysis
  object to setup the runnable and then write the results back to the 
  database specified in the config file.

=head1 CONTACT

  Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BestTargetted;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BestTargetted;

use Bio::EnsEMBL::Analysis::Runnable::BestTargetted;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info stack_trace_dump );





use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->read_and_check_config($BEST_TARGETTED_CONFIG);
  
  
  my $seqfetcher = $self->make_seqfetcher($self->SEQFETCHER_DIR, 
                                          $self->SEQFETCHER_OBJECT);
  $self->seqfetcher($seqfetcher);
  
  
  return $self; 
}



sub fetch_input{
  
  my ($self) = @_;

  $self->throw("No input id") unless defined($self->input_id);
  
  my $dnadb = $self->get_dbadaptor('REFERENCE_DB');
  my $db = $self->get_dbadaptor($self->DB_NAME);
  $db->dnadb($dnadb);
  
  my $slice = $self->fetch_sequence($self->input_id, $db);
  $self->query($slice);
  
  my @genes;
  foreach my $biotype (@{$self->BIOTYPES}) {
    push @genes, @{ $slice->get_all_Genes_by_type($biotype) }  ;
  }
  print "\nGot ".scalar(@genes)." genes\n";
  $self->genes(\@genes);
#  my @genes = @{ $slice->get_all_Genes($self->PRIMARY_LOGICNAME) }  ;
#  print "\nGot ".scalar(@genes)." genes from primary logic_name\n";
#  push @genes, @{ $slice->get_all_Genes($self->SECONDARY_LOGICNAME) }  ;
#  print "Now have ".scalar(@genes)." genes in total";
#  $self->genes(\@genes);
  

  return 1;
}



sub run {
  my ($self) = @_;

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BestTargetted->new
    (
     -query             => $self->query,
     -analysis          => $self->analysis,
     -seqfetcher        => $self->seqfetcher,
     -biotypes          => $self->BIOTYPES,
     -program           => $self->EXONERATE_PROGRAM,
     -verbose           => $self->VERBOSE,
     -keep_single_analysis => $self->KEEP_SINGLE_ANALYSIS,
     -genes             => $self->genes,
     
    );

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
        throw("BestTargetted::run failed $@");
      }
    }
  } else { 
    print $@; 
  }  
  
  $self->output($runnable->output);

  1;
}



sub write_output{
  my ($self) = @_;

  my $out_dba = $self->get_dbadaptor($self->OUT_DB_NAME);

  my $gene_a = $out_dba->get_GeneAdaptor() ; 
  info ("trying to write output") ;  
  foreach my $gene (@{$self->output}){
     info("storing $gene" ); 
    $gene_a->store($gene) ; 
  }  
  return ;
}


=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: if $index exists, 
           returns a Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs, otherwise throws
 Example :
 Returns : Bio::DB::RandomAccessI
 Args    : $indexname - string


=cut



sub make_seqfetcher {
  my ( $self, $index, $seqfetcher_class ) = @_;

  my $seqfetcher;
  
  (my $class = $seqfetcher_class) =~ s/::/\//g;

  throw ("Configuration-error !! There's no entry for SEQFETCHER_OBJECT in Targetted.pm\n") if (length($class)==0)  ;

  require "$class.pm";

  if(defined $index && $index ne ''){
    my @db = ( $index );
    
    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => \@db, );
  }
  else{
    $self->throw("can't make seqfetcher\n");
  }

  return $seqfetcher;
}


####################################
# config variable holders
####################################



sub read_and_check_config {
  my ($self, $hash) = @_;

  $self->SUPER::read_and_check_config($hash);
 
  my $logic = $self->analysis->logic_name;

  foreach my $var (qw(BIOTYPES 
                      DB_NAME
                      SEQFETCHER_DIR
                      SEQFETCHER_OBJECT)) {

    throw("You must define $var in config for logic '$logic'" . 
          " or in the DEFAULT entry")
        if not $self->$var;
        
  }
  
  if (not $self->OUT_DB_NAME) {
    warn ("\n\tOUT_DB_NAME has not been provided in config-file Analysis/Config/GeneBuild/BestTargetted.pm\n ".
          "\toutput will be written to " . $self->DB_NAME . "\n");
    $self->OUT_DB_NAME($self->DB_NAME);   
  }
  
}


  
   






##############
# containers
##############

sub seqfetcher {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_seq_fetcher} = $value;
  }
  return $self->{_seq_fetcher};
}


sub genes {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_genes} = $value;
  }
  return $self->{_genes};
}



################
# CONFIG VARS
################

sub BIOTYPES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_biotypes} = $value;
  }
  return $self->{_biotypes};
}


sub SEQFETCHER_DIR {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_seqfetcher_dir} = $value;
  }
  return $self->{_seqfetcher_dir};
}


sub SEQFETCHER_OBJECT {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_seqfetcher_object} = $value;
  }
  return $self->{_seqfetcher_object};
}


sub DB_NAME {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_dbname} = $value;
  }
  return $self->{_dbname};
}


sub OUT_DB_NAME {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_outdb_name} = $value;
  }
  return $self->{_outdb_name};
}


sub EXONERATE_PROGRAM {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_exon_prog} = $value;
  }
  return $self->{_exon_prog};
}

sub VERBOSE {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_verbose} = $value;
  }
  return $self->{_verbose};
}

sub KEEP_SINGLE_ANALYSIS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_keep_single_analysis} = $value;
  }
  return $self->{_keep_single_analysis};
}

1;
