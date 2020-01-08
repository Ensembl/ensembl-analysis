=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Runnable::Aspera

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Runnable::Aspera;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(locate_executable execute_with_wait);

use File::Spec::Functions qw(catfile splitpath);


=head2 new

 Arg [PROGRAM]      : String
 Arg [USE_THREADING]: Integer
 Arg [VERBOSE]      : Integer
 Description: Creates a Bio::EnsEMBL::Analysis::Runnable::Samtools object to run samtools command
 Returntype : Bio::EnsEMBL::Analysis::Runnable::Samtools
 Exceptions : None

=cut

sub new {
  my ($class,@args) = @_;
  my $self = bless {},$class;
  my ($program, $sshkey, $verbose) = rearrange( [ qw( PROGRAM SSHKEY VERBOSE) ], @args);


  $program = 'ascp' unless ($program);
  $self->program($program);
  if (!$sshkey and exists $ENV{LINUXBREW_HOME}) {
    $sshkey = catfile($ENV{LINUXBREW_HOME}, 'etc', 'asperaweb_id_dsa.openssh');
  }
  $self->sshkey($sshkey);
  logger_verbosity($verbose) if ($verbose);

  return $self;
}

=head2 program

 Arg [1]    : String (optional) $program
 Description: Getter/setter for the binary. It tries to locate properly the executable
 Returntype : String
 Exceptions : Throw if it can't find the executable Arg[1]

=cut

sub program {
    my ($self, $file) = @_;

    if ($file) {
        $self->{_program} = locate_executable($file);
    }
    return $self->{_program};
}


=head2 sshkey

 Arg [1]    : String path to sshkey file
 Description: Getter/setter for the ssh key file
 Returntype : String path
 Exceptions : Throws if the file does not exists

=cut

sub sshkey {
  my ($self, $file) = @_;

  if ($file) {
    if (-e $file) {
      $self->{_sshkey_file} = $file;
    }
    else {
      throw('Could not find ssh key file '.$file);
    }
    if ($self->options) {
      $self->options =~ s/-i \S+/-i $file/;
    }
    else {
      $self->options('-i '.$file);
    }
  }
  return $self->{_sshkey_file};
}


sub options {
  my ($self, $options) = @_;

  if ($options) {
    $self->{_options} = $options;
    if ($self->sshkey and $options !~ /-i \S+/) {
      $self->{_options} .= ' -i '.$self->sshkey;
    }
  }
  return $self->{_options};
}


sub source {
  my ($self, $source) = @_;

  if ($source) {
    $self->{_source} = $source;
  }
  return $self->{_source};
}


sub target {
  my ($self, $target) = @_;

  if ($target) {
    $self->{_target} = $target;
  }
  return $self->{_target};
}

sub fetch {
  my ($self, $source, $target) = @_;

  $source = $self->source unless ($source);
  $target = $self->target unless ($target);
  throw("Could not find $target") unless (-e $target);
  my $options = $self->options;
  my $cmd = $self->program.' '.$self->options.' '.$source.' '.$target;
  execute_with_wait($cmd);
  my (undef, undef, $file) = splitpath($source);
  return catfile($target, $file);
}


sub put {
  my ($self, $source, $target) = @_;

  $source = $self->source unless ($source);
  $target = $self->target unless ($target);
  throw("Could not find $source") unless (-e $source);
  my $cmd = $self->program.' '.$self->options.' '.$source.' '.$target;
  execute_with_wait($cmd);
}

1;
