=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Minimap

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Runnable::Minimap;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Analysis::Runnable::Samtools;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($target_file, $use_threading, $samtools, $verbose) = rearrange (['TARGETFILE', 'USE_THREADING', 'SAMTOOLS', 'VERBOSE'], @args);
  throw('You need a file with your sequences') unless ($self->queryfile);
  throw('Target file does not exists '.$target_file) unless (-e $target_file);
  $self->target_file($target_file);
  $self->use_threading($use_threading || 0);
  $self->samtools($samtools || 'samtools', $verbose);
  return $self;
}


sub run {
  my ($self) = @_;

  throw('No query file '.$self->queryfile) unless (-s $self->queryfile);
  $self->resultsfile($self->queryfile.'_unsorted.bam');
  my $cmd = $self->program;
  if ($self->options) {
    $cmd .= ' '.$self->options;
  }
  $cmd .= ' '.$self->target_file.' '.$self->queryfile;
  my $to_bam_cmd = $self->samtools->make_commandline('view', '-b', $self->resultsfile, '-');
  $cmd .= ' | '.$to_bam_cmd;
  execute_with_wait($cmd, 'Could not run minimap');
  throw("Something went wrong for '$cmd'") unless (-s $self->resultsfile > 5000);
  $self->samtools->sort($self->queryfile, $self->resultsfile);
  $self->samtools->index($self->queryfile.'.bam');
  unlink $self->resultsfile || throw('Could not delete '.$self->resultsfile);
  $self->output([$self->queryfile.'.bam']);
}

=head2 use_threading

 Arg [1]    : Boolean, 0 or 1
 Example    : $sef->use_threading(1);
 Description: Getter/Setter, set this to 1 if you want picard to use threads, is ~20% faster but use ~20% more memory
 Returntype : Boolean
 Exceptions : None

=cut

sub use_threading {
    my ($self, $use_threading) = @_;

    if ($use_threading) {
        $self->{_use_threading} = $use_threading;
    }
    return $self->{_use_threading};
}


=head2 samtools

 Arg [1]    : String, '/usr/bin/samtools'
 Example    : $self->samtools('/usr/bin/samtools');
 Description: Creates a object to run samtools command
 Returntype : Bio::EnsEMBL::Analysis::Runnable::Samtools
 Exceptions : None

=cut

sub samtools {
    my ($self, $samtools, $verbose) = @_;

    if ($samtools) {
        $self->{_samtools} = Bio::EnsEMBL::Analysis::Runnable::Samtools->new(
                              -program => $samtools,
                              -verbose => $verbose,
                              -use_threading => $self->use_threading
                              );
    }
    return $self->{_samtools};
}

1;
