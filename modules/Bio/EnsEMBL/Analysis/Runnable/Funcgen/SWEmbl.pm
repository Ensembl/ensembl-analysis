# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'program.pl',
  );
  $runnable->run;
  my @features = @{$runnable->output};

=head1 DESCRIPTION

SWEmbl expects bed or maq mapview files as input and predicts features which 
can be stored in the annotated_feature table in the eFG database

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Funcgen);

=head2 run

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::SWEmbl
    Arg [2]     : filename
    Description : 
    Returntype  : 
    Exceptions  : 
    Example     : 

=cut

sub run {

    print "Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl::run\n";
    my ($self, $dir) = @_;

    $self->run_analysis;


    print "Parsing results ... ";
    $self->parse_results();
    print "done!\n";

}


=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::SWEmbl
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
    throw($program." is not executable SWEmble::run_analysis ") 
        unless($program && -x $program);

    my @fields = (0..2,5);
    $self->output_fields(\@fields);

    (my $resultsfile = $self->infile());
    $self->resultsfile($resultsfile);
    warn("RESULTS FILE: ".$resultsfile);
    
    my $command = $self->program . ' -B -z -i ' . $self->query . ' ' . 
        $self->options . ' -o ' . $resultsfile;
    
    warn("Running analysis " . $command);

    eval { system($command) };
    throw("FAILED to run $command: ", $@) if ($@);

}

#=head2 infile
#
#  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::SWEmble
#  Arg [2]     : filename (string)
#  Description : will hold a given filename or if one is requested but none
#  defined it will use the create_filename method to create a filename
#  Returntype  : string, filename
#  Exceptions  : none
#  Example     : 
#
#=cut
#
#
#sub infile{
#
#  my ($self, $filename) = @_;
#
#  if($filename){
#    $self->{'infile'} = $filename;
#  }
#  if(!$self->{'infile'}){
#    $self->{'infile'} = $self->create_filename($self->analysis->logic_name, 'dat');
#  }
#
#  return $self->{'infile'};
#
#}


sub query{
  my $self = shift;
  $self->{'query'} = shift if(@_);

  throw("file ".$self->{'query'}. " doesn't exist") if (! -e $self->{'query'});

  return $self->{'query'};
}


sub config_file {
    my $self = shift;
    $self->{'config_file'} = shift if(@_);
    return $self->{'config_file'};
}

1;
