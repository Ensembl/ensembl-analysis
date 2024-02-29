# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Test::Checker


=head1 SYNOPSIS
Module to define the interface for and implement the basic functionality of
the various Checker classes

=head1 DESCRIPTION

=head1 CONTACT

  Steve Searle <searle@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Checker;
use vars qw(@ISA $AUTOLOAD);
use strict;


sub ignorewarnings {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_ignorewarnings} = $arg;
  }
  return $self->{_ignorewarnings};
}


=head2 add_Error

 Title   : add_Error
 Usage   : $obj->add_Error($newval)
 Function:
 Returns : value of errors


=cut

sub add_Error {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      push @{$self->{_errors}},$value;
    }
    return @{$self->{_errors}};
}

sub get_all_Errors {
   my $self = shift;
   if (!defined($self->{_errors})) {
     @{$self->{_errors}} = ();
   }
   return @{$self->{_errors}};
}

sub add_Warning {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      push @{$self->{_warnings}},$value;
    }
    return @{$self->{_warnings}};
}

sub get_all_Warnings {
   my $self = shift;
   if (!defined($self->{_warnings})) {
     @{$self->{_warnings}} = ();
   }
   return @{$self->{_warnings}};
}

sub has_Errors {
  my $self = shift;

  if (scalar($self->get_all_Errors) ||
      (scalar($self->get_all_Warnings) && !$self->ignorewarnings)) {
    return 1;
  }
  return 0;
}

sub output {
  my $self = shift;

  if (scalar($self->get_all_Errors)) {
    print "\nErrors:\n";
    foreach my $error ($self->get_all_Errors()) {
      print "  " . $error;
    }
  }
  if (scalar($self->get_all_Warnings)) {
    print "\nWarnings:\n";
    foreach my $warning ($self->get_all_Warnings) {
      print "  " . $warning;
    }
  }
}


sub check {
  my $self = shift;

  $self->throw("Object did not provide the check method of the Checker interface");
}

1;
