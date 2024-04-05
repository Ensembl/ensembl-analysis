#!/usr/bin/env perl
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

=head1 NAME

ensembl-analysis/scripts/recalculate_exon_phases.pl

=head1 DESCRIPTION

The script fetches all transcripts from the given database,
runs calculate_exon_phases and prints the result of the calculate exon_phases sub
if there is any difference between them and the existing phases.

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)
    -dbport    what port to connect to (port= in locator)
    -dbname    what name to connect to (dbname= in locator)
    -dbuser    what username to connect as (user= in locator)
    -dbpass    what password to use (pass= in locator)

=cut

use warnings ;
use strict;

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';

GetOptions( 'dbhost|host|h:s'               => \$host,
            'dbport|port|P:n'               => \$port,
            'dbname|db|D:s'                 => \$dbname,
            'dbuser|user|u:s'               => \$dbuser,
            'dbpass|pass|p:s'               => \$dbpass,
);

if (!$host || !$dbname || !$dbuser) {
  my $message =
    "Need -dbhost '$host' -dbuser '$dbuser' and -dbname '$dbname' to run. ".
    'Use -help for more detailed documentation.';
  throw($message);
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
(
   -dbname          => $dbname,
   -host            => $host,
   -user            => $dbuser,
   -port            => $port,
   -pass            => $dbpass,
);

my $sa = $db->get_SliceAdaptor;
foreach my $s (@{$sa->fetch_all('toplevel',undef,1)}) {
  foreach my $t (@{$s->get_all_Transcripts}) {
    if ($t->translate()) {
      my @ep;

      foreach my $e (@{$t->get_all_Exons()}) {
        push(@ep, [$e->phase(),$e->end_phase()]);
      }
  
      my $start_phase = $t->get_all_translateable_Exons()->[0]->phase();
      calculate_exon_phases($t,$start_phase);
  
      my @recalculated_exon_phases = ();
      my $i = 0;
      foreach my $e (@{$t->get_all_Exons()}) {
        if ($ep[$i]->[0] != $e->phase or $ep[$i]->[1] != $e->end_phase) {
 
          print($t->seq_region_name()." "); 

          if ($t->stable_id()) {
            print($t->stable_id());
          }
          else {
            print($t->dbID());
          }
          print(" ");
          if ($e->stable_id()) {
            print($e->stable_id());
          }
          else {
            print($e->dbID());
          }
          
          print(
             " original phases: ".
             $ep[$i]->[0]." ".
             $ep[$i]->[1].
             " recalculated phases: ".
             $e->phase()." ".
             $e->end_phase().
             " start phase used: ".
             $start_phase.
             "\n"
          );
        } # end if different phases
        ++$i;
      } #end foreach e
    } #end if t translate
  } #end foreach t
} #end foreach s

1;
