# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Config::GeneBuild::Solexa2GenesLiteNew;

use strict;
use vars qw( %Config );

# Hash containing config info
%Config = (

SOLEXA2GENES_CONFIG_BY_LOGIC => 
     {
    
      DEFAULT => {
                  # databases are defined as hash keys from Bio::EnsEMBL::Analysis::Config::Databases
                  OUTPUT_DB    => '',
	          ALIGNMENT_DB => '',
		  
		  
                  # logic_names for dna_align_features to be fetched out of ALIGNMENT_DB
                  LOGIC_NAMES  => [],
		  
		  # logic_name for repeats will be used to merge exons separated by repeats 
                  REPEAT_LN  => 'repeatmask',	

                  # options for filtering out small gene models
                  MIN_LENGTH => 100,
                  MIN_EXONS  =>   1,
		  
		  #�genes with a span < MIN_SPAN will also be considered single exon
		  MIN_SINGLE_EXON_LENGTH => 1000,
		  
                  # 'span' = genomic extent / cdna length
                  MIN_SPAN   =>   1.5,
		  
                  # number of read pairs needed to confirm an end exon
		  END_EXON_COVERAGE => 2,
		  
                  # exclude exons where the coverage is lower than 
		  # *percentage* of the *average* exon coverage for the transcript
		  # ie a value of 10 means remove all exon pairings where the
		  # coverage is < 10% of the average 
		  EXCLUDE_EXONS => 0.1,
		
		  # Maximum number of reads to process at a time
		  # slows down the processing but uses less memory
		  # fetches all reads if left blank
		  BATCH_FETCH  => ,
                  },
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
