# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen::SWEmbl

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

=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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

    throw($program." is not executable") if(! ($program && -x $program));

    my @fields = (0..2,5,9);
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
