# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Blast
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::BlastRfam

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::BlastRfam->
  new(
      -query => $slice,
      -program => 'wublastn',
      -database => 'embl_vertrna',
      -options => 'hitdist=40 -cpus=1',
      -parser => $bplitewrapper,
      -filter => $featurefilter,
     );
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Runnable::BlastRfam;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Blast);





=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Blast
  Function  : override results parsing to allow for 
coverage calculations
  Returntype: none
  Exceptions: none
  Example   : 

=cut


sub parse_results{
  my ($self) = @_;
  my $results = $self->results_files;
  my @daf_coverage_results;
  my $filtered_output;
  my $bplite = $self->parser->get_parsers($results);
  foreach my $blast (@{$bplite}){
      while( my $subject = $blast->nextSbjct){
	 while (my $hsp = $subject->nextHSP) {
	   my @daf_results;
	   my $hsp_length = $hsp->length."\t";
	   my $subject_length = $subject->{'LENGTH'};
	   $subject_length = 1 unless ($subject_length);
	   my $coverage = $hsp_length/$subject_length*100;
	   $coverage =~ s/\.\d+//;
#	   print "subject $subject_length hsp $hsp_length coverage $coverage\n";
	   unless ($coverage > 70){
	     next;
	   }
	   $subject->name =~ /^(\S+)\/\S+\s+(\w+);\w+/;
	   my $name = $2."-".$1;
	   push @daf_results, $self->parser->split_hsp($hsp,$name);
	   # add coverage into daf score?
	   foreach my $daf(@daf_results){
	     $daf->score($coverage);
	     push  @daf_coverage_results, $daf;
	   }
	 }
      }
    }
#  print "Before clustering ".scalar( @daf_coverage_results)."\n";
  return undef unless  ( @daf_coverage_results);
  my $output = $self->cluster(\@daf_coverage_results);
  $self->output($output);
}


sub cluster{
  my ($self,$dafs_ref)=@_;
  my @dafs = @$dafs_ref;
  my $start =0;
  my @representative_sequence;
 DAFS: foreach my $daf (@dafs){
    $start ++;
    next DAFS unless($daf);
    my @temp_array;
    push @temp_array,$daf;
  MATCHES:  for (my $index = $start ; $index <= $#dafs ; $index ++){
      next MATCHES unless ($dafs[$index]);
      if ($daf->overlaps($dafs[$index])){
	push @temp_array,$dafs[$index];
	$dafs[$index] = undef;
      }
    }
    # get highest scoring alignment to be representative
    @temp_array = sort{$a->p_value <=> $b->p_value} @temp_array;
    push @representative_sequence,shift @temp_array;;
  }
  return \@representative_sequence;
}



1;
