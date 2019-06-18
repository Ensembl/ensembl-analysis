=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::MergeBamFiles -

=head1 SYNOPSIS

  [mergebamfiles]
  module=MergeBamFiles

=head1 DESCRIPTION

  Merge the bam files to create the "super" bam file needed by
  Bam2Genes, the rough model step.
  It uses samtools to check but it uses picard for merging,
  no config file as everything is fetched from the analysis table
  The module is looking for files named: *_sorted.bam

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::MergeBamFiles;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Analysis::Config::MergeBamFiles qw (MERGE_BAM_FILES_BY_LOGIC);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


sub new{
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->read_and_check_config($MERGE_BAM_FILES_BY_LOGIC);
  return $self;
}

sub fetch_input {
    my ($self) = @_;

    if (scalar(@{$self->INPUT_FILES}) == 0) {
        throw('You did not specify input files for '.$self->analysis->logic_name);
    }
    elsif (scalar(@{$self->INPUT_FILES}) == 1 and $self->OPTIONS !~ /-b /) {
        $self->input_is_void(1);
    }
    if (defined $self->PICARD_LIB) {
        $self->require_module('Bio::EnsEMBL::Analysis::Runnable::PicardMerge');
        $self->runnable(Bio::EnsEMBL::Analysis::Runnable::PicardMerge->new(
            -program => $self->JAVA || 'java',
            -java_options => $self->JAVA_OPTIONS,
            -lib => $self->PICARD_LIB,
            -options => $self->OPTIONS,
            -analysis => $self->analysis,
            -output_file => $self->OUTPUT_FILE,
            -input_files => $self->INPUT_FILES,
            -use_threading => $self->USE_THREADING,
            -samtools => $self->SAMTOOLS || 'samtools',
            ));
    }
    else {
        $self->require_module('Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge');

        $self->runnable(Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge->new(
            -program => $self->SAMTOOLS || 'samtools',
            -options => $self->OPTIONS,
            -analysis => $self->analysis,
            -output_file => $self->OUTPUT_FILE,
            -input_files => $self->INPUT_FILES,
            -use_threading => $self->USE_THREADING,
            ));
    }
}

sub run {
    my ($self) = @_;

    foreach my $runnable (@{$self->runnable}) {
        $runnable->run;
        $runnable->check_output_file;
    }
    return 1;
}

sub write_output {
    my ($self) = @_;

    return 1;
}



sub OUTPUT_FILE {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_OUTPUT_FILE'} = $value;
    }

    if (exists($self->{'_CONFIG_OUTPUT_FILE'})) {
        return $self->{'_CONFIG_OUTPUT_FILE'};
    } else {
        return undef;
    }
}


sub INPUT_FILES {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_INPUT_FILES'} = $value;
    }

    if (exists($self->{'_CONFIG_INPUT_FILES'})) {
        return $self->{'_CONFIG_INPUT_FILES'};
    } else {
        return undef;
    }
}


sub PICARD_LIB {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_PICARD_LIB'} = $value;
    }

    if (exists($self->{'_CONFIG_PICARD_LIB'})) {
        return $self->{'_CONFIG_PICARD_LIB'};
    } else {
        return undef;
    }
}


sub OPTIONS {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_OPTIONS'} = $value;
    }

    if (exists($self->{'_CONFIG_OPTIONS'})) {
        return $self->{'_CONFIG_OPTIONS'};
    } else {
        return undef;
    }
}


sub SAMTOOLS {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_SAMTOOLS'} = $value;
    }

    if (exists($self->{'_CONFIG_SAMTOOLS'})) {
        return $self->{'_CONFIG_SAMTOOLS'};
    } else {
        return undef;
    }
}


sub GENES_DB {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_GENES_DB'} = $value;
    }

    if (exists($self->{'_CONFIG_GENES_DB'})) {
        return $self->{'_CONFIG_GENES_DB'};
    } else {
        return undef;
    }
}


sub JAVA_OPTIONS {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_JAVA_OPTIONS'} = $value;
    }

    if (exists($self->{'_CONFIG_JAVA_OPTIONS'})) {
        return $self->{'_CONFIG_JAVA_OPTIONS'};
    } else {
        return undef;
    }
}


sub USE_THREADING {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_USE_THREADING'} = $value;
    }

    if (exists($self->{'_CONFIG_USE_THREADING'})) {
        return $self->{'_CONFIG_USE_THREADING'};
    } else {
        return undef;
    }
}


sub JAVA {
    my ($self,$value) = @_;

    if (defined $value) {
        $self->{'_CONFIG_JAVA'} = $value;
    }

    if (exists($self->{'_CONFIG_JAVA'})) {
        return $self->{'_CONFIG_JAVA'};
    } else {
        return undef;
    }
}

1;
