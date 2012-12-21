=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/Runnable/ProteinAnnotation.pm,v $
# $Revision: 1.7 $
package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use File::Copy qw(mv cp);
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::ProteinFeature;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);    

  my ($database) = rearrange(['DATABASE'], @args);

  if(!$database){
    $database = $self->analysis->db_file;
  }
  $self->database($database);

  my ($threshold_file) = rearrange(['THRESHOLDS'], @args);
  $self->thresholds($threshold_file);

  
  if(!$self->query){
    throw("need to have a query defined in order to run protein ".
          "annotation");
  }
  if(!$self->analysis){
    throw("need to have a analysis defined in order to run protein ".
          "annotation");
  }
  if(not $self->program){
    if ($self->analysis->program_file) {
      $self->program($self->analysis->program_file);
    } elsif ($self->analysis->program) {
      $self->program($self->analysis->program);
    }
  }
    
  return $self; # success - we hope!
}


sub run{
  my ($self, $dir) = @_;
  
  $self->workdir ('/tmp') unless ($self->workdir($dir));
  $self->checkdir;

  my @input_files;
    
  if(-s $self->query and not $self->multiprotein){
    # The input is a sequence file. but we have to break up 
    # the file into single-sequence entries for this analysis
    my %files = %{$self->get_individual_protein_files};
    foreach my $file(keys(%files)){
      my $id = $files{$file};
      $self->queryfile($file);
      $self->run_analysis();
      $self->parse_results($id);
    }
  } else {
    if (-s $self->query) {
      $self->queryfile($self->query);
      $self->run_analysis;
      $self->parse_results;
    }
    elsif (ref($self->query) and
           $self->query->isa("Bio::PrimarySeqI")) {
      #The input is a sequence object
      # write sequence to file
      my $filename = $self->write_seq_file($self->query);
      $self->files_to_delete($filename);
      $self->queryfile($filename);
      $self->run_analysis;
      $self->parse_results($self->query->id);
    } else {
      throw("Can't run if ".$self->query." isn't either a Bio::PrimarySeq  " .
            " or a file which has a size greater than 0");
    }

  }

  $self->delete_files;
}

sub get_individual_protein_files {
  my ($self) = @_;
  
  if(!$self->{_protein_files}){
    $self->{_protein_files} = {};
  }

  my $in  = Bio::SeqIO->new(-file => $self->query, '-format' =>'Fasta');
  while(my $tmpseq = $in->next_seq()){
    
    my $stub = $self->analysis->logic_name.".".$tmpseq->display_id.".$$";
    my $filename = $self->create_filename($self->analysis->logic_name, 
                                          $stub . "seq");
    $filename = $self->write_seq_file($tmpseq, $filename);
    $self->{_protein_files}{$filename} = $tmpseq->display_id;
    $self->files_to_delete($filename);
  }

  return $self->{_protein_files};
}


#######################################
sub create_protein_feature{
  my ($self, $start, $end, $score, $seqname, $hstart, 
      $hend, $hseqname,$analysis, $p_value, $percent_id) = @_;

  my $fp = Bio::EnsEMBL::ProteinFeature->new(
                                             -start    => $start,
                                             -end      => $end,
                                             -hstart   => $hstart,
                                             -hend     => $hend,
                                             -percent_id => $percent_id,
                                             -score    => $score,
                                             -p_value  => $p_value,
                                             -hseqname => $hseqname,
                                             -seqname => $seqname,
                                             -analysis => $analysis,
                                            );

  return $fp;
}

###################################
sub query {
  my ($self, $seq) = @_;
  
  if ($seq) {
    if (ref($seq) and $seq->isa("Bio::PrimarySeqI")) {
      $self->{_sequence} = $seq;
    } elsif (-e $seq) {
      $self->{_sequence} = $seq;
    } else {
      throw("You must provide either a Bio::Seq or a file which ".
            " exists not $seq");
    }
  }
  
  return $self->{_sequence};
}


##################################
sub queryfile{
  my ($self, $filename) = @_;

  if($filename){
    $self->{_queryfile} = $filename;
  }
  if(not exists $self->{_queryfile}){
    $self->{_queryfile} = $self->create_filename($self->analysis->logic_name . '.seq', 
                                                 'fa');
    $self->files_to_delete($self->{_query_file});
  }

  return $self->{_queryfile};
}

##############################
sub resultsfile{
  my ($self, $filename) = @_;
  
  if($filename){
    $self->{_resultsfile} = $filename;
  }
  if(not exists $self->{_resultsfile}){
    $self->{_resultsfile} = $self->create_filename($self->analysis->logic_name . '.results',
                                                   'out');
    $self->files_to_delete($self->{_resultsfile});
  }

  return $self->{_resultsfile};
}

################################
sub database {
  my ($self, $database) = @_;

  if (defined $database) {
    $self->{_database} = $database;
  }

  return $self->{_database};
} 

##################################
sub multiprotein{
  my ($self) = @_;
  return 1;
}

sub thresholds {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_thresholds} = $val;
  }

  return $self->{_thresholds};
}





1;
