=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Runnable::Indicate

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Runnable::Indicate;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(locate_executable execute_with_wait);
use File::Spec::Functions qw(catfile splitpath);

sub new {
  my ($class,@args) = @_;
  my $self = bless {},$class;
  my ($program, $dbname, $verbose) = rearrange( [ qw( PROGRAM DBNAME VERBOSE) ], @args);


  $program = 'indicate' unless ($program);
  $self->program($program);
  $self->dbname($dbname);
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


sub dbname {
  my ($self, $dbname) = @_;

  if ($dbname) {
    $self->{_dbname} = $dbname;
  }
  return $self->{_dbname};
}

sub build_index {
  my ($self, @files) = @_; #I try to mimic Bio::DB::Flat::BinarySearch so I don't have to do to many tests

  my $cmd = $self->program.' --index '.catfile($self->directory, $self->dbname).' --parser '.$self->primary_parser;
  foreach my $file (@files) {
    my (undef, $dir, $name) = splitpath($file);
    $cmd .= ' -d '.$dir.' -f '.$file;
    execute_with_wait($cmd);
  }
}

1;
