
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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam;
use Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Hmmpfam);


#################################

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($scop_map_file,
      $evalue) = rearrange(['SCOPMAP',
                            'EVALUE',
                            ], @args);
    
  $self->scop_map($scop_map_file)
      if defined $scop_map_file;

  if (defined $evalue) {
    $self->evalue_cutoff($evalue);
  } else {
    # set it to something sensible
    warning("No evalue cutoff given; defaulting to 0.02");
    $self->evalue_cutoff("0.02");
  }

  my $opts = defined($self->options) ? $self->options : "";
  $opts =~ s/\-E\s+\S+//; 
  $self->options("$opts -E " . $self->evalue_cutoff);
  
  return $self;
}


sub parse_results {
  my ($self) = @_;

  $self->SUPER::parse_results;

  # although the search will have been done with the
  # required evalur cutoff, we have to filter again here
  # on the basis of per-domain evalue (as opposed to 
  # whole sequence evalue)

  my @pass_eval;
  foreach my $hit (@{$self->output}) {
    if (not defined $self->evalue_cutoff or 
        $hit->p_value <= $self->evalue_cutoff) {
      push @pass_eval, $hit;
    }
  }
  
  my $filter = Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter->new();  
  my $out = $filter->filter_results(\@pass_eval);
  #my $out = $self->output;

  # finally, change the names to scop ids if we have the mappings
  if ($self->scop_map) {
    foreach my $hit (@$out) {
      if (exists($self->scop_map->{$hit->hseqname})) {
        my $scopname = $self->scop_map->{$hit->hseqname};
        $hit->hseqname( $scopname );
      }
    }
  }

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

sub scop_map {
  my ($self, $file) = @_;

  my %scop_map;

  if (defined $file) {
    open SCOPMAP, $file or throw("Could not open SCOP mapfile '$file'");
    while(<SCOPMAP>) {
      /^(\S+)\s+(\S+)/ and $scop_map{$1} = $2;
    }
    close(SCOPMAP);
    
    $self->{_scop_map} = \%scop_map;
  }

  return $self->{_scop_map};
}

sub evalue_cutoff {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_evalue_cutoff} = $val;
  }

  if (not exists $self->{_evalue_cutoff}) {
    return undef;
  } else {
    return $self->{_evalue_cutoff};
  }
}


1;
