=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::miRNA - 

=head1 SYNOPSIS

     my $runnableDB = Bio::EnsEMBL::Analysis::RunnableDB::miRNA->new(
            	-db       => $db_adaptor,
		-input_id => 'analysis logic name',
		-analysis => $analysis,
	        );
    $runnabledb->fetch_input();
    $runnabledb->run();
    $runnabledb->write_output();

=head1 DESCRIPTION

RunnableDB to provide database access for miRNA detection.
Runs as an accumulator job on miRNA blast hits found by Bio::EnsEMBL::RunnableDB::BlastmiRNA
Takes an analysis logic name as an input id and uses it to fetch all dna align features associated with 
that analysis.
It then groups the dna align features by miRNA families and ignores families with > 50
members as there is a high probability that these are hitting repetitive sequences.
Creates and runs the miRNA runnable.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::miRNA;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild ; 
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::miRNA;
use Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


=head2 fetch_input

  Title      : fetch_input
  Usage      : $miRNA->fetch_input();
  Function   : opens and stores connections to databases in Bio::EnsEMBL::Analysis::Config::Databases
             : fetches all dna align features by analysis logic name specified in the 
             : input id 
  Returns    : Hash reference
  Exceptions : throws if the analysis object is not found or no dna align features are retrieved
  Args       : None

=cut

sub fetch_input{
  my ($self) = @_;

  # dna database
  my $dna_db = $self->get_dbadaptor($DNA_DBNAME) ;

  # if you want to write the final genes into the pipeline database need 
  # to catch it first and store the $self->db as the genes->db otherwise the
  # registry will cause problems 
  if ( $$DATABASES{'GENEBUILD_DB'}{'-dbname'} eq $self->db->dbc->dbname &&
       $$DATABASES{'GENEBUILD_DB'}{'-port'} == $self->db->dbc->port &&
       $$DATABASES{'GENEBUILD_DB'}{'-host'} eq $self->db->dbc->host){
         $self->gene_db($self->db);
  } else { 
    my $genes_db = $self->get_dbadaptor("GENEBUILD_DB");
    $self->gene_db($genes_db);
  }
  $self->db->dnadb($dna_db);
  my $aa = $self->db->get_AnalysisAdaptor;
  my $analysis = $aa->fetch_by_logic_name($self->input_id);
  $self->throw("Analysis BlastmiRNA not found $@\n") unless $analysis;
  my $dafa = $self->db->get_DnaAlignFeatureAdaptor;
  my @flags;
  my @dafs;
  my $fa = Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor->new($self->db);
  print "Fetching features\n";
  eval{
    @flags = @{$fa->fetch_by_analysis($self->analysis)};
  };
  foreach my $flag (@flags){
    if ($flag->goalAnalysis->logic_name eq $self->analysis->logic_name){
      my $daf = $dafa->fetch_by_dbID($flag->ensembl_id);
      push @dafs, $daf;
    }
  }
  $self->throw("No dna align features found ") unless (scalar(@dafs) >=1);
  print scalar(@dafs)." dafs found\n";
  my %families = %{$self->family(\@dafs)};
  # empty the array
  @dafs = ();
  
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::miRNA->new
    (
     -queries => \%families,
     -analysis => $self->analysis,
    );
  $self->runnable($runnable);
}


=head2 family

  Title      : family
  Usage      : my %families = %{$self->family(\@dafs)};
  Function   : order dafs by family, removes families with blast hits to repeats
  Returns    : Hash reference
  Exceptions : None
  Args       : Array ref of Bio::EnsEMBL::DnaDnaAlignFeatures

=cut

sub family{
  my ($self,$dafs_ref) = @_;
  my %families;
  foreach my $daf (@$dafs_ref){
    push @{$families{$daf->hseqname}},$daf;
  }
  my %filtered_fam;
  foreach my $key (keys %families){
  if (scalar @{$families{$key}} <= 50){
    $filtered_fam{$key} = $families{$key};
    } else {
      # take top scoring 50 hits
      my @array = sort {$a->p_value <=> $b->p_value} @{$families{$key}};
      my @filtered_array =  splice(@array,0,50);
      $filtered_fam{$key} = \@filtered_array;
    }
 }
  return \%filtered_fam;
}

=head2 write_output

  Args       : none
  Description: Writes the single exon miRNA genes into the final genebuild database,
             : also stores attributes associated with the transcript
  Exceptions : Throws if gene or transcript attribute fail to write to the database
  Returntype : scalar

=cut

sub write_output{
  my ($self) = @_;
  my $adaptor = $self->gene_db->get_GeneAdaptor;
  my $aa = $self->gene_db->get_AttributeAdaptor;
  my @attributes; 
  my $xref;
  foreach my $gene_hash (@{$self->output}){
    my $gene = $gene_hash->{'gene'};
    @attributes = @{$gene_hash->{'attrib'}};
    $gene->analysis($self->analysis);
    $gene->status('PREDICTED');
    foreach my $trans (@{$gene->get_all_Transcripts}){
      $trans->analysis($self->analysis);
      $trans->status('PREDICTED');
    }
    $gene->slice($self->query) if(!$gene->slice);
    $self->feature_factory->validate($gene);
    eval{
      $adaptor->store($gene);
 	print STDERR "Attemting to store in ".$adaptor->db->dbname."\n";
    };
    if($@){
      $self->throw("miRNA:store failed, failed to write ".$gene." to ".
		   "the database $@");
    }
    foreach my $trans (@{$gene->get_all_Transcripts}){
      eval{
	$aa->store_on_Transcript($trans->dbID,\@attributes);
	$self->gene_db->get_TranscriptAdaptor->update($trans);
	$self->gene_db->get_GeneAdaptor->update($gene);	
      };
      if($@){
	$self->throw("miRNA:store failed, failed to write ".@attributes." on transcript ".
		     $trans." in the database $@");
      }
    }
  }
  return 1;
}


#########################################################
# Containers

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

1;
