# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie

=head1 SYNOPSIS

  my $nessie = Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'nessie',
  );
  $nessie->run;
  my @annotated_features = @{$nessie->output};

=head1 DESCRIPTION

Nessie expects to run the program Nessie and predicts features which 
can be stored in the predicted_feature table in the eFG database

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Funcgen);


=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Nessie
  Arg [2]     : string, program name
  Usage       : 
  Description : 
  Returns     : 
  Exceptions  : 

=cut

sub run_analysis {
        
    my ($self, $program) = @_;

    if(!$program){
        $program = $self->program;
    }

    throw($program." is not executable Nessie::run_analysis ") 
        unless($program && -x $program);
    
    (my $outfile = $self->infile()) =~ s/\.dat$/.out/;
    my @fields = (0..2,4);
    $self->output_fields(\@fields);

    (my $resultsfile = $self->infile()) =~ s/\.dat$/_peaks.out/;
    
    my $replicates = scalar(keys %{$self->result_features});
    if ($replicates == 1) { $replicates=2 }

    my $command = $self->program.' --data="'.$self->infile().'" '.
        '--replicates='.$replicates.
        $self->analysis->parameters. ' > '.$outfile;

    warn("Running analysis " . $command . "\n");
    eval { system($command) };
    throw("FAILED to run $command: ", $@) if ($@);

    my $parser = 'oligo_peaks.pl --gap 500 --sang --split 2500 '.
        '--bridge --mindist 400 --minpost 0.99 --autonorm -l '.
        $outfile.' '.$self->infile.' > '.$resultsfile; 
        
    warn("Running parser " . $parser . "\n");
    eval { system($parser) };
    throw("FAILED to run $parser: ", $@) if ($@);

    throw("No peak resultfile.") if (! -e $resultsfile);
    $self->resultsfile($resultsfile);
    
}

=head2 write_infile

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Nessie
    Arg [2]     : filename
    Description : 
    Returntype  : 
    Exceptions  : 
    Example     : 

=cut

sub write_infile {

    my ($self, $filename) = @_;

    if (! $filename) {
        $filename = $self->infile();
    }


    foreach my $rset_name (sort keys %{$self->result_features})
    {

        (my $datafile = $filename) =~ s,\.dat,.$rset_name.dat,;

        open(F, ">".$datafile)
            or throw("Can't open file $datafile.");

        foreach (keys %{$self->result_features}) {
            next unless ($datafile =~ m/$_/);
            foreach my $ft (@{$self->result_features->{$_}}) {
                #print join("\t", @$ft), "\n";
                print F join("\t", @$ft), "\n";
            }
            last;
        }
        
        close F;

        $self->datafiles($datafile);

    }
    #print Dumper $self->datafiles;

    $self->write_filelist();
    map {$self->files_to_delete($_)} @{$self->datafiles};

    return $filename;

}

=head2 infile

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Nessie
  Arg [2]     : filename (string)
  Description : will hold a given filename or if one is requested but none
  defined it will use the create_filename method to create a filename
  Returntype  : string, filename
  Exceptions  : none
  Example     : 

=cut


sub infile{

  my ($self, $filename) = @_;

  if($filename){
    $self->{'infile'} = $filename;
  }
  if(!$self->{'infile'}){
    $self->{'infile'} = $self->create_filename('nessie', 'dat');
  }

  return $self->{'infile'};

}

=head2 datafiles

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Nessie
  Arg [2]     : filename (string)
  Description : will hold/add a list of filename
  Returntype  : listref of filenames
  Exceptions  : none
  Example     : 

=cut


sub datafiles{

  my ($self, $filename) = @_;

  if(! $self->{'datafiles'}){
      $self->{'datafiles'} = [];
  }

  if($filename){
      push @{$self->{'datafiles'}}, $filename;
  }

  return \@{$self->{'datafiles'}};

}

=head2 write_filelist

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Nessie
    Description : writes the nessie specific wrapper file containing a 
                  list of replicate data file names
    Returntype  : none
    Exceptions  : 
    Example     : 

=cut

sub write_filelist {

    my ($self) = shift;

    throw("No infile found") if (! $self->infile());
    
    open(F, ">".$self->infile())
        or throw("Can't open file $self->infile.");

    foreach my $f (@{$self->datafiles}) {
        print F $f, "\n";
    }

    # nessie hack: if there is only one replicate data set we duplicate this one ...
    if (scalar(@{$self->datafiles}) == 1) {
        warn("Duplicating input data file.");
        print F ${$self->datafiles}[0], "\n";
    }

    close F;

}

1;
