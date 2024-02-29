=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Array - 

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


=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Array;

use warnings ;
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
  
  my $options = "--showalignment no --bestn 100 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 25 --dnawordthreshold 11 --querytype $query_type --targettype $target_type  --target $target_dir --query " ;

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
                      '-analysis'     => $self->analysis,
								     );
  $self->runnable($runnable);
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->output into $self->dbobj
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
