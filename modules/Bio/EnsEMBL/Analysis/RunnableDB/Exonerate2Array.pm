#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array->new(
					     -dbobj     => $db,
					     -input_id  => $id,
					     -analysis   => $analysis
                                             );
    $obj->fetch_input();
    $obj->run();

    my @newfeatures = $obj->output();
    
    $obj->write_output();

=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Array;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateArray;
use Bio::EnsEMBL::Analysis::Config::General;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Array
  Function  : fetch data out of fasta files and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 


=cut

sub fetch_input {
  my( $self) = @_;
  print STDERR "Fetching input \n";
  
  my $input_id = $self->input_id;
  my $analysis = $self->analysis;
  my $program = $analysis->program_file;
  my $query_type = 'dna';
  my $target_type = 'dna';
  my $query_file = $ANALYSIS_INPUT_DIR.$input_id;
  my $target_dir =$ANALYSIS_TARGET_DIR;
  my $verbose = "all";
  
  my $options = "--showalignment no --bestn 100 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 25 --dnawordthreshold 11 --querytype $query_type --targettype $target_type  --target $target_dir --query " ;

  verbose($verbose);
  #$target_dir .= "22.fa"; ##only for testing

  ###Runnable::ExonerateArray take a array of query_seq_obj, so it's need to be generated here###

  my @query_seqs;

  my $in = Bio::SeqIO->newFh(
			     -FILE => $query_file,
			     -FORMAT => 'Fasta',
			    );

  while (my $seq = <$in>) {
    push (@query_seqs, $seq);
  }


  # prepare runnable
  
  throw("Can't run Exonerate without both query and target sequences") 
    unless (defined($query_file) && defined($target_dir));
  
  info("exonerate is '$program', target_dir is $target_dir, query_file is $ query_file\n");

  #my $target_file = $target_dir . "*";###exonerate-0.8.2 can use both file and dir
  my $runnable = new Bio::EnsEMBL::Analysis::Runnable::ExonerateArray(
								      '-db'           => $self->db,
								      '-query_seqs'   => \@query_seqs,
								      '-program'      => $program,
								      '-options'      => $options,
								      '-verbose'      => $verbose,
								     );
  $self->runnable($runnable);
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->{_output} into $self->dbobj
    Returns :   1
    Args    :   None

=cut

sub write_output {

  my($self) = @_;
  
  my @misc_features = @{$self->output()}; 
  
  my $mfa = $self->db->get_MiscFeatureAdaptor();
  $mfa->store( @misc_features );
  
  return 1;
}

1;
