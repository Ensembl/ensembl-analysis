# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::miRNA
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::miRNA

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

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::miRNA;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::miRNA;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


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

  my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_DBHOST,
     '-user'   => $GB_DBUSER,
     '-dbname' => $GB_DBNAME,
     '-pass'   => $GB_DBPASS,
     '-port'   => $GB_DBPORT,
    );

  my $genes_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_FINALDBHOST,
     '-user'   => $GB_FINALDBUSER,
     '-dbname' => $GB_FINALDBNAME,
     '-pass'   => $GB_FINALDBPASS,
     '-port'   => $GB_FINALDBPORT,
    );
  $self->gene_db($genes_db);
  # add dna_db 
  $self->db->dnadb($dna_db);

  my $aa = $self->db->get_AnalysisAdaptor;
  my $analysis = $aa->fetch_by_logic_name($self->input_id);
  $self->throw("Analysis ".$self->input_id." not found $@\n") unless $analysis;
  my $dafa = $self->db->get_DnaAlignFeatureAdaptor;
  my @dafs = @{$dafa->generic_fetch(" analysis_id = ".$analysis->dbID)};
  $self->throw("No dna align features found ") unless (scalar(@dafs) >=1);
  my %families = %{$self->family(\@dafs)};
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
    $filtered_fam{$key} = $families{$key} if (scalar @{$families{$key}} <= 50);

# This bit would take top 50 instaed of ignoring the familly ....
#    my @array = sort {$a->p_value <=> $b->p_value} @{$families{$key}};
#    if ($filtered_fam{$key}){
#      print "Familly $key has ".scalar(@{$filtered_fam{$key}})." members\n";
#    } else { 
#      print "Familly $key rejected had  ".scalar(@{$families{$key}})." members\n";
#    }
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
  my $dbea = $self->gene_db->get_DBEntryAdaptor;
  my @attributes; 
  my $xref;
  foreach my $gene_hash (@{$self->output}){
    my $gene = $gene_hash->{'gene'};
    @attributes = @{$gene_hash->{'attrib'}};
    $xref = $gene_hash->{'xref'};
    $gene->analysis($self->analysis);
    $gene->slice($self->query) if(!$gene->slice);
    $self->feature_factory->validate($gene);
    eval{
      $adaptor->store($gene);
    };
    if($@){
      $self->throw("miRNA:store failed, failed to write ".$gene." to ".
		   "the database $@");
    }
    foreach my $trans (@{$gene->get_all_Transcripts}){
      eval{
	$aa->store_on_Transcript($trans,\@attributes);
	$dbea->store($xref, $trans->dbID, 'Transcript') if $xref;
	$self->gene_db->get_TranscriptAdaptor->update($trans);
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
