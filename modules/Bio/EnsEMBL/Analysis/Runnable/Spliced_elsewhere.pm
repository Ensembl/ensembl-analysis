# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere ->new
  (
   '-genes' => \@genes_array_ref,
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

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;


use strict;
use Bio::EnsEMBL::Analysis::Runnable::Pseudogene;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Pseudogene);


=head2 new

  Args       : various
  Description: Runnable constructor
  Returntype : Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;
  Caller     : general

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_genes'} = [];	#array of genescripts to test;  

  my($genes) = rearrange([qw( GENES)], @args);

  if ($genes) {
    $self->genes($genes);
  }
  # Path to blast database
  $self->db("$PS_MULTI_EXON_DIR"."/all_multi_exon_genes.fasta");
  return $self;
}

=head2 run

Arg [none] :
  Description: calls run_blast for each genescript in turn
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run {
  my ($self)=@_;
  my @genes = @{$self->genes};
 GENE: foreach my $gene(@genes){
    $self->run_blast($gene);
    next GENE;
  }
  return 1;
}

=head2 run_blast

Arg [none] :
  Description: runs the blast analysis, uses BPliteWapper as parser
  Returntype : returns the blast result in a hash keyed by the transcript dbid
also returns the gene object in the hash for refence
  Exceptions : none
  Caller     : general

=cut

sub  run_blast{
  my ($self,$gene)=@_;
  my $bplitewrapper = Bio::EnsEMBL::Analysis::Tools::BPliteWrapper-> new
    (
     -query_type => 'dna',
     -database_type => 'dna',
    );
  my %output_hash;
  foreach my $trans (@{$gene->get_all_Transcripts}){ 
    next unless ($trans->translateable_seq);
    my $query = $trans->feature_Slice->expand(1000,1000)->get_repeatmasked_seq;
    my $test = $trans->feature_Slice->get_repeatmasked_seq->seq;
  #  print "BEFORE  Ns have gone " . length($test) . "\n"; 
    $test =~ s/N//g;
#    print "AFTER Ns have gone " . length($test) . "\n";
    unless ( length($test) > 5 ) {
      print STDOUT "ignoring " . $trans->display_id. " all Ns \n";
      $output_hash{$trans->dbID}= 'REPEATS';
      $self->output(\%output_hash);
    return 1;
    }
    print  $trans->stable_id if ( $trans->stable_id ) ; 
    # need to mask out other coding sequences here
    my $blast =  Bio::EnsEMBL::Analysis::Runnable::Blast->new 
      ('-query'     => $query,
       '-program'   => 'blastn',
       '-database'  => $self->db,
       '-threshold' => 1e-50,
       '-parser'    => $bplitewrapper,
       '-options'   => 'V=10',
       '-analysis'  => $self->analysis,
      );
    eval {
      $blast->run();
    };
    if ($@) {
      $self->throw("Problem with Blast $@ \n ");
    } else {
      if ( $blast->output ) {
	$output_hash{$trans->dbID}= $blast->output;
      } else {
	$output_hash{$trans->dbID}= 'NONE';
      }
    }
  }
  $self->output(\%output_hash);
  return 1;
}



########################################################
# Containers


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

=head2 output

  Arg [1]    : hasref
  Description: overrides output array
  Returntype : array
  Exceptions : none
  Caller     : general

=cut

sub output {
  my ($self, $hash_ref) = @_;
  if ($hash_ref) {
    push @{$self->{'_output'}},$hash_ref;
  }
  return $self->{'_output'};
}
