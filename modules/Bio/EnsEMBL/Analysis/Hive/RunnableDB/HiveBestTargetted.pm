=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestTargetted

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestTargetted;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BestTargetted;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
    my $self = shift;

    return {
        %{$self->SUPER::param_defaults},
        keep_single_analysis => 1, # If set to 0, do not write model created by only one analysis
        cluster_on_coding_exons => 1,
        exonerate_program => 'exonerate', # You should use version 0.9.0
        verbose => 1,
    }
}


sub fetch_input{
  my ($self) = @_;

  $self->create_analysis;
  $self->throw("No input id") unless defined($self->input_id);
  my $seqfetcher = $self->make_seqfetcher($self->SEQFETCHER_DIR,
                                          $self->SEQFETCHER_OBJECT);
  $self->seqfetcher($seqfetcher);
  my $dnadb = $self->get_database_by_name('dna_db');

  $self->hrdb_set_con($self->get_database_by_name($self->OUT_DB_NAME, $dnadb), 'target_db');

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
     my $input_db = $self->hrdb_get_dba($self->param('source_db')->{$db_alias}, $dnadb) ;
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
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BestTargetted->new
    (
     -query             => $self->query,
     -analysis          => $self->analysis,
     -seqfetcher        => $self->seqfetcher,
     -biotypes          => $self->BIOTYPES,
     -program           => $self->EXONERATE_PROGRAM,
     -verbose           => $self->VERBOSE,
     -keep_single_analysis => $self->KEEP_SINGLE_ANALYSIS,
     -genes             => \@genes,
     -cluster_on_coding_exons => $self->CLUSTER_ON_CODING_EXONS,
    );
    if ($self->param_is_defined('protein_min_coverage')) {
      $runnable->min_coverage($self->param('protein_min_coverage'));
    }
    if ($self->param_is_defined('protein_min_identity')) {
      $runnable->min_identity($self->param('protein_min_identity'));
    }
  $self->runnable($runnable);
#  my @genes = @{ $slice->get_all_Genes($self->PRIMARY_LOGICNAME) }  ;
#  print "\nGot ".scalar(@genes)." genes from primary logic_name\n";
#  push @genes, @{ $slice->get_all_Genes($self->SECONDARY_LOGICNAME) }  ;
#  print "Now have ".scalar(@genes)." genes in total";
#  $self->genes(\@genes);


  return 1;
}



sub write_output{
  my ($self) = @_;

  my $out_dba = $self->hrdb_get_con($self->OUT_DB_NAME);

  my $gene_a = $out_dba->get_GeneAdaptor() ;
  my $analysis = $self->analysis;
  foreach my $gene (@{$self->output}){
    $gene->analysis($analysis);
    $gene_a->store($gene) ;
  }
  return ;
}


=head2 make_seqfetcher

 Arg [1]    : Arrayref of String, indexes to be used
 Arg [2]    : String, class of the indexer, usually Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher
 Description: Returns the initialised indexer
 Returntype : Bio::DB::RandomAccessI
 Exceptions : Throws if Arg[1] does not exist
              Throws if Arg[2] does not exist

=cut

sub make_seqfetcher {
  my ( $self, $index, $seqfetcher_class ) = @_;

  my $seqfetcher;

  (my $class = $seqfetcher_class) =~ s/::/\//g;

  $self->throw ("Configuration-error !! There's no entry for SEQFETCHER_OBJECT in Targetted.pm\n") if (length($class)==0)  ;

  require "$class.pm";

  if(defined $index && $index ne ''){
    my $dbs;
    if (ref($index) eq 'ARRAY') {
      $dbs = $index;
    }
    else {
      $dbs = [ $index ];
    }

    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => $dbs, );
  }
  else{
    $self->throw("can't make seqfetcher\n");
  }

  return $seqfetcher;
}


##############
# containers
##############

sub seqfetcher {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('_seq_fetcher', $value);
  }
  return $self->param('_seq_fetcher');
}



################
# CONFIG VARS
################

sub BIOTYPES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('BIOTYPES', $value);
  }
  return $self->param('BIOTYPES');
}


sub INPUT_DATA_FROM_DBS  {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('INPUT_DATA_FROM_DBS', $value);
  }
  return $self->param('INPUT_DATA_FROM_DBS');
}



sub SEQFETCHER_DIR {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('SEQFETCHER_DIR', $value);
  }
  return $self->param('SEQFETCHER_DIR');
}


sub SEQFETCHER_OBJECT {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('SEQFETCHER_OBJECT', $value);
  }
  return $self->param('SEQFETCHER_OBJECT');
}


sub DB_NAME {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('dbname', $value);
  }
  return $self->param('dbname');
}


sub OUT_DB_NAME {
  my ($self,$value) = @_;

#  if (defined $value) {
#    $self->param('target_db', $value);
#  }
#  return $self->param('target_db');
  return 'target_db';
}


sub EXONERATE_PROGRAM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('exonerate_program', $value);
  }
  return $self->param('exonerate_program');
}

sub VERBOSE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('verbose', $value);
  }
  return $self->param('verbose');
}

sub KEEP_SINGLE_ANALYSIS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('keep_single_analysis', $value);
  }
  return $self->param('keep_single_analysis');
}

sub CLUSTER_ON_CODING_EXONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('cluster_on_coding_exons', $value);
  }
  return $self->param('cluster_on_coding_exons');
}

1;
