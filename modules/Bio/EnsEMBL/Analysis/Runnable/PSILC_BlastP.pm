# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::PSILC;

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::PSILC_BlastP->new 
    (
     '-trans'    => $transcript,
     '-analysis' => $analysis,
    );

  $blast->run;
  my $output = $blast->output;


=head1 DESCRIPTION

BLASTP runnable for Bio::EnsEMBL::Analysis::RunnableDB::PSILC
Runs BLASTP for the supplied transcript against a multi species blast database.
Filters the results and returns a hash ref of homologous transcript identifiers

=head1 CONTACT

Simon White

sw4@sanger.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::PSILC_BlastP;
use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 new

  Arg [1]   : array of Bio::EnsEMBL::Transcript
  Function  : cretes PSILC BLASTP runnable
  Returntype: Bio::EnsEMBL::Analysis::Runnable::PSILC_BlastP
  Exceptions: none

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_trans'} = {};	#gene to test;  

  my( $trans) = $self->_rearrange([qw(
				     TRANS
				    )], @args);
  if ($trans) {
    $self->trans($trans);
  }
  return $self;
}

=head2 run

  Arg [1]   : none
  Function  : runs Bio::EnsEMBL::Analysis::Runnable::PSILC_BlastP
  Returntype: none
  Exceptions: none

=cut

sub run {
  my ($self)=@_;
  my $trans = $self->trans;
  $self->run_blast($trans);
}

=head2 run

  Arg [1]   : none
  Function  : runs BLAST on transcript against multi species blast database
  Returntype: Hash containing transcript ids of homologs that pass criteria:
target is non identical to query has at least 80% coverage and 50% sequence identity
  Exceptions: throws if it cannot recognise the transcript identifier

=cut

sub  run_blast{
  my ($self,$trans)=@_;
  my %transcript_hash;

  print STDERR "Blast search of ".$trans->dbID.".\n";
  # might want to use minimal blast module
  my $bplitewrapper = Bio::EnsEMBL::Analysis::Tools::BPliteWrapper-> new
    (
     #  -regex => '^\w+\s+(\w+)'
     -query_type => 'dna',
     -database_type => 'dna',
    );
  # make a query object containing the cds

  my $query = Bio::Seq->new(
			    '-display_id' => $trans->dbID,
			    '-seq'        => $trans->translate->seq,
			    );

  my $blast =  Bio::EnsEMBL::Analysis::Runnable::Blast->new 
    ('-query'     => $query,
     '-program'   => 'blastp',
     '-database'  => "$PSILC_BLAST_DB/multi_species.fasta",
     '-threshold' => 1e-6,
     '-parser'    => $bplitewrapper,
     '-options'   => 'V=10 -cpus=1',
     '-analysis'  => $self->analysis,
    );

  $blast->run();
  my $results = $blast->output; 

  # get dnadnaalignfeatures for each HSP
  # include in set if > 50% ID
  # dont include if > 99%id and 80% coverage
  # Coverage based on cds not cdna length

  foreach my $daf (@{$results}) {
    my $coverage = $daf->length/length($query->seq);
    if ($daf->percent_id <= 99 && 
	$coverage > 0.8  && 
	$daf->percent_id > 50 ){    
      # organise the transcripts in to a hash to avoid redundancy and
      # group them by species 
      if ($daf->hseqname =~ /(\w+)_(\w+)/) {
	push @{$transcript_hash{$2}},$1;
      } else {
	$self->throw("Cannot recognise trans_id ".$daf->hseqname."  $@\n");
      }
    }
  }
  $self->output(\%transcript_hash);
}



#######################################
# Containers


=head2 trans

Arg [1]    : array ref
  Description: get/set trans set to run over
  Returntype : array ref to Bio::EnsEMBL::Transcript objects
  Exceptions : throws if not a Bio::EnsEMBL::Transcript
  Caller     : general

=cut

sub trans {
  my ($self, $trans) = @_;
  if ($trans) {
    unless  ($trans->isa("Bio::EnsEMBL::Transcript")){
      $self->throw("Input isn't a Bio::EnsEMBL::Transcript, it is a $trans\n$@");
    }
    $self->{'_trans'} = $trans;
  }
  return $self->{'_trans'};
}

=head2 output

  Arg [1]    : hash ref
  Description: over-ridden output in runnable.pm to allow hash refs rather than array refs
  Returntype : array ref to Bio::EnsEMBL::Transcript objects
  Exceptions : none
  Caller     : general

=cut

sub output{
  my ($self, $trans) = @_;
  if ($trans) {
    $self->{'_output'} = $trans;
  }
  return $self->{'_output'};
}



1;
