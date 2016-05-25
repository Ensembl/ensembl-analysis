# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::Config::KillListFilter

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::KillListFilter;

=head1 DESCRIPTION

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first declared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Exonerate> hash is asked to be set.

Since The RunnableDB that this config controls can be used for 
inferring transcript structures from (different sets of) EST, cDNA 
and proteins, and several uses may be required in the same pipeline, 
this Config contains one primary config variable, EXONERATE_TRANSCRIPT_CONFIG.
This is hash keyed off logic name, each entry of which is a hash
containing the variable that affect the behaviour of the RunnableDB.
When the RunnableDB instance is created, the correct entry is identified
by logic name and value for a corresponding set of local variables are
set.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Analysis::Config::GeneBuild::KillListFilter;

use strict;
use vars qw( %Config );

# Hash containing config info
%Config = (
  KILL_LIST_CONFIG_BY_LOGIC => {
    DEFAULT => {
      FILTER_PARAMS => {
        -before_date => undef,
        -for_analyses => [],
        -for_external_db_ids => [],
        -for_species => [],
        -having_status => undef,
        -only_mol_type => 'protein',
        -reasons => [],
        -source_species => undef,
        -user_name => undef
      },
      GB_REF_DB => 'REFERENCE_DB',
      KILL_LIST_DB => 'KILL_LIST_DB'
    },
    EST => {
      FILTER_PARAMS => {
        -before_date => undef,
        -for_analyses => [],
        -for_external_db_ids => [],
        -for_species => [],
        -from_source_species => undef,
        -having_status => undef,
        -only_mol_type => 'EST',
        -reasons => [],
        -user_id => undef
      }
    },
    PROTEIN => {
      FILTER_PARAMS => {
        -before_date => undef,
        -for_analyses => [],
        -for_external_db_ids => [],
        -for_species => [],
        -from_source_species => undef,
        -having_status => undef,
        -only_mol_type => 'protein',
        -reasons => [],
        -user_id => undef
      }
    },
    cDNA => {
      FILTER_PARAMS => {
        -before_date => undef,
        -for_analyses => [],
        -for_external_db_ids => [],
        -for_species => [],
        -from_source_species => undef,
        -having_status => undef,
        -only_mol_type => 'cDNA',
        -reasons => [],
        -user_id => undef
      }
    },
    cdna_update => {
      FILTER_PARAMS => {
        -before_date => undef,
        -for_analyses => [
          'cdna_update'
        ],
        -for_external_db_ids => [],
        -for_species => [
          9606
        ],
        -from_source_species => undef,
        -having_status => undef,
        -only_mol_type => 'cDNA',
        -reasons => [],
        -user_id => undef
      }
    }
  }
);


sub import {
  my ($callpack) = caller(0); # Name of the calling package
  my $pack = shift; # Need to move package off @_

  # Get list of variables supplied, or else everything
  my @vars = @_ ? @_ : keys( %Config );
  return unless @vars;
  
  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $Config{$_} ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$Config{ $_ };
	} else {
	    die "Error: Config: $_ not known\n";
	}
    }
}

1;
