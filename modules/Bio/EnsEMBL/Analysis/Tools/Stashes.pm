# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Tools::Stashes;

use warnings;
use strict;
no strict "refs";

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Exporter;
use vars qw (@ISA @EXPORT %stash $alias @alias %alias);

@ISA    = qw(Exporter);
@EXPORT = qw( package_stash );

sub package_stash {
  my ($packageName) = @_;

  my %result;

  local (*alias);
  *stash = *{"${packageName}::"};

  while ( my ( $varName, $globValue ) = each %stash ) {
    # only return the config hash
    next if $varName =~ m/BEGIN/;
    next if $varName =~ m/import/;

    *alias = $globValue;
    $result{$varName} = $alias  if ( defined($alias) );
    $result{$varName} = \@alias if ( *alias{ARRAY} );
    $result{$varName} = \%alias if ( *alias{HASH} );
  }

  if ( scalar( keys %result > 1 ) ) {
    throw( "Have more than one item exported from " .
           "$packageName - you'll run into trouble\n" );
  }
  my $hash_name = ( keys %result )[0];
  return [ $result{$hash_name}, $hash_name ];
} ## end sub package_stash

1;
