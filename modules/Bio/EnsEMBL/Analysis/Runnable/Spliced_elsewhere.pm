# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere ->new
  (
   '-trans' => \@transcript_array_ref,
   '-analysis' => $analysis_object,
  );
    $runnable->run;
    $output = $runnable->output

=head1 DESCRIPTION

Runnable for Bio::EnsEMBL::Analysis::RunnableDB::Spliced_elsewhere
Does the blast analysis and returns the results to the runnable_db for
parsing.
Uses Bio::EnsEMBL::Analysis::Config::Pseudogene for config

=head1 CONTACT

Simon White

sw4@sanger.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;


use strict;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


=head2 new

  Args       : various
  Description: Runnable constructor
  Returntype : Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;
  Caller     : general

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_trans'} = [];	#array of transcripts to test;  

  my($trans) = $self->_rearrange([qw(
				     TRANS
				    )], @args);
  if ($trans) {
    $self->trans($trans);
  }
  # Path to blast database
  $self->db("$PS_MULTI_EXON_DB"."multi_exon_blastdb.fasta");
  return $self;
}

=head2 run

Arg [none] :
  Description: calls run_blast for each transcript in turn
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run {
  my ($self)=@_;
  my @trans = @{$self->trans};
  foreach my $trans(@trans){
    if ($trans->translateable_seq){
      $self->run_blast($trans);
    }
  }
  return 1;
}

=head2 run_blast

Arg [none] :
  Description: runs the blast analysis, uses BPliteWapper as parser
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub  run_blast{
  my ($self,$trans)=@_;
  my $bplitewrapper = Bio::EnsEMBL::Analysis::Tools::BPliteWrapper-> new
    (
     -query_type => 'pep',
     -database_type => 'pep',
    );

  my $query = $trans->translate;
  my $blast =  Bio::EnsEMBL::Analysis::Runnable::Blast->new 
    ('-query'     => $query,
     '-program'   => 'blastp',
     '-database'  => $self->db,
     '-threshold' => 1e-6,
     '-parser'    => $bplitewrapper,
     '-options'   => 'V=10',
     '-analysis'  => $self->analysis,
    );

  $blast->run();
  my @array =($trans, $blast->output);
  push @{$self->output},\@array;
  return 1;
}

########################################################
# Containers

=head2 trans

  Arg [1]    : array ref
  Description: get/set transcript set to run over
  Returntype : array ref to Bio::EnsEMBL::Transcript objects
  Exceptions : throws if not a Bio::EnsEMBL::Transcript
  Caller     : general

=cut

sub trans {
  my ($self, $trans) = @_;
  if ($trans) {
    foreach my $tran (@{$trans}) {
      unless  ($tran->isa("Bio::EnsEMBL::Transcript")){
	$self->throw("Input isn't a Bio::EnsEMBL::Transcript, it is a $tran\n$@");
      }
    }
    $self->{'_trans'} = $trans;
  }
  return $self->{'_trans'};
}

=head2 db

  Arg [1]    : scalar
  Description: get/set path to blast database
  Returntype : scalar
  Exceptions : none
  Caller     : general

=cut

sub db {
  my ($self, $db) = @_;
  if ($db) {
    $self->find_file($db);
    $self->{'_db'} = $db;
  }
  return $self->{'_db'};
}

1;
