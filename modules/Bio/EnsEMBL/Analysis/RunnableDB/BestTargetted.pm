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

Bio::EnsEMBL::Analysis::RunnableDB::BestTargetted - 

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

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BestTargetted;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BestTargetted;

use Bio::EnsEMBL::Analysis::Runnable::BestTargetted;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info stack_trace_dump );
use FileHandle;




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
  
  my $slice = $self->fetch_sequence($self->input_id, $dnadb);
  $self->query($slice);
 
  my @genes; 
  my @all_bt_for_clustering;  

  # fetch genes from different databases  
  foreach my $db_alias ( keys  %{ $self->INPUT_DATA_FROM_DBS } ) { 
     print "fetching out of DB : " . $db_alias ." : " ;  
     my @biotypes_to_fetch = @{${$self->INPUT_DATA_FROM_DBS }{$db_alias}};
     for ( @biotypes_to_fetch ) { 
       print $_ . " " ;  
     }
     print "\n";  
     # get db adaptor for db 
     my $input_db = $self->get_dbadaptor($db_alias) ;
     # get slice 
     my $input_slice = $self->fetch_sequence($self->input_id, $input_db);
     # get biotypes  
     #
     for my $bt ( @biotypes_to_fetch ) {  
        my @genes_fetched = @{ $input_slice->get_all_Genes_by_type($bt,undef,1) }  ; 
        print "-> $bt    ".scalar(@genes_fetched) . " genes fetched \n" ;  
        push @genes, @genes_fetched; 
        push @all_bt_for_clustering , $bt ; 
     }
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
     -cluster_on_coding_exons => $self->CLUSTER_ON_CODING_EXONS, 
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

  foreach my $var (qw(SEQFETCHER_DIR
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


sub INPUT_DATA_FROM_DBS  {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_input_data_from_dbs} = $value;
  }
  return $self->{_input_data_from_dbs};
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

sub CLUSTER_ON_CODING_EXONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{_cluster_on_coding_exons} = $value;
  }
  return $self->{_cluster_on_coding_exons};
}

1;
