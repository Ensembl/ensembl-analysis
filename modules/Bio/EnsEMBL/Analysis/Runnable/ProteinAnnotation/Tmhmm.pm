=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Tmhmm - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut


=pod 

=head1 NAME

 Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Tmhmm

=head1 SYNOPSIS

 my $seqstream = Bio::SeqIO->new ( -file => $clonefile,
                                   -fmt => 'Fasta',
                                 );
 $seq = $seqstream->next_seq;

 my $tmhmm = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Tmhmm->new(-analysis => $ana,
                                                                             -query    => $seq,
                                                                             -modelfile => $mf,
                                                                             -optionsfile => $o);
 $tmhmm->run;
 my @results = $tmhmm->output;

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Tmhmm;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;


@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);


sub multiprotein{
  my ($self) = @_;
  return 1;
}


sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($options_file, 
      $model_file) =  rearrange(['OPTIONSFILE',
                                 'MODELFILE'
                                 ], @args);
  
  throw("You must supply a model file for Tmhmm")
      if not defined $model_file;
  
  $self->options_file($options_file)
      if defined $options_file;
  
  $self->model_file($model_file);

  return $self;
}


sub run_analysis {
  my ($self) = @_;
  # run program
  my $cmd = "cat " . $self->queryfile . " | " .  $self->program;
  if ($self->options_file) {
    $cmd .= " -f " . $self->options_file;
  }
  if ($self->options) {
    $cmd .= " " .  $self->options;
  }
  $cmd .=  " -modelfile " . $self->model_file . " > " . $self->resultsfile;
  
  print "Running $cmd\n";

  # decodeanhmm returns a non-zero exit status, even when successful
  # So we just have to run it and wait until we parse the output
  # file to see whether it ran successfully
  system($cmd);
}


sub parse_results {
  my ($self) = @_;

  my ($fh, $id, @pfs);

  my $resfile = $self->resultsfile;
  
  if (-e $resfile) {
    # it's a filename
    if (-z $resfile) {  
      return;
    } else {
      open ($fh, "<$resfile") or throw("Error opening $resfile");
    }
  } else {
    # it'a a filehandle
    $fh = $resfile;
  }
  
  my $id_count = 0;
  while (<$fh>) {
    chomp;
    next if /^$/;
    if (/^\>(\S+)/) {
      $id = $1;
      $id_count++;
    }
    elsif (/^\%pred/) {
      my ($junk, $values) = split /:/;
      my @tm = split (/,/, $values);
      foreach (@tm) {
        /(\w+)\s+(\d+)\s+(\d+)/;
        my $orien = $1;
        my $start = $2;
        my $end = $3;
        $orien = uc ($orien);
        if ($orien eq "M") {
          my $fp = $self->create_protein_feature($start, $end, 0, $id, 0, 
                                                 0, 'Tmhmm', 
                                                 $self->analysis, 0, 0);
          push @pfs, $fp;
        }
      }
    }
  }
  close($fh);

  throw("Something went wrong when running decodeanhmm; no ids found in output")
      if not $id_count;

  $self->output(\@pfs);
}


sub model_file {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_model_file} = $val;
  }

  return $self->{_model_file};
}


sub options_file {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_options_file} = $val;
  }

  return $self->{_options_file};
}



1;
