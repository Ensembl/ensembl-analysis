
=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Protein::Hmmpfam

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Pipeline::Runnable::ProteinAnnotation::Superfamily->new 
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => '/superfamily/db');
  $hmm->run;
  my @results = @{$hmm->output};

=head1 DESCRIPTION

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Superfamily;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam;
use Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam);


#################################
sub parse_results {
  my ($self) = @_;

  $self->SUPER::parse_results;
  
  my $filter = Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter->new();  
  my $out = $filter->filter_results($self->output);
  #my $out = $self->output;

  $self->output($out);
}


sub output {
  my ($self, $output) = @_;

  if (defined $output) {
    throw("Must pass Runnable:output an arrayref not a ".$output)
        unless(ref($output) eq 'ARRAY');
    $self->{output} = $output;
  }

  if (not exists $self->{output}) {
    $self->{output} = [];
  }

  return $self->{output};
}

1;
