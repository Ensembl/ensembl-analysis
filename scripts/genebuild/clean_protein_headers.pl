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

#!/usr/bin/env perl

use strict; 
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::PolyA;

$| = 1; # disable buffering

my $protfile;
my $seqoutfile;
my $min_length = 10;
my $species_name;
&GetOptions( 
	    'protfile:s'     => \$protfile,
	    'outfile:s'     => \$seqoutfile,
	    'min_length:n'    => \$min_length,
      'short_name:s'  => \$species_name, 
	   );

# usage
if(!defined $protfile    ||
   !defined $seqoutfile
  ){
  print  "USAGE: clean_protein_headers.pl -protfile protfile -outfile outfile".
    " -min_length <min_prot_length> -short_name <species_short_name>\n";
  exit(1);
}

if(!defined $species_name
    ){
    print  "No short species name specified, must specify using -short_name flag: e.g. SORAR for common shrew, PIG for pig etc.\n";
    print 
  exit(1);
}


my $seqin  = new Bio::SeqIO(-file   => "<$protfile",
			    -format => "Fasta",
			  );

my $seqout = new Bio::SeqIO(-file   => ">$seqoutfile", 
			    -format => "Fasta"
			   );

my %hadhash;
my $in = 0;
my $out = 0;
my $too_short = 0;
my $unrecognised_header = 0;
my $unrecognised_description = 0;
my $duplicate = 0;
my $np = 0;
my $xp = 0;
my $yp = 0;
my $ap = 0;
my $sp = 0;
my $tr = 0;

SEQFETCH:
while( my $prot = $seqin->next_seq ){
  $in++;
  
  if ($prot->length <= $min_length) {
    print STDERR "Length < $min_length: rejecting ".$prot->display_id." with length ".$prot->length."\n";
    $too_short++;
    next;
  }
  
  my $display_id  = $prot->display_id;
  #print "$display_id\n";
  my $description = $prot->desc;
  #print "$description\n";

  if ( $display_id =~/gi\|\S+\|ref\|(NP\_\d+\.\d+)\|/) {
    $display_id = $1;
    $np++;
  } elsif ( $display_id =~/gi\|\S+\|ref\|(XP\_\d+\.\d+)\|/) {
    #print STDERR "Skipping $1\n";
    $xp++;
    next SEQFETCH;
  } elsif ( $display_id =~/gi\|\S+\|ref\|(AP\_\d+\.\d+)\|/) {
    $display_id = $1;
    $ap++;
  } elsif ( $display_id =~/gi\|\S+\|ref\|(YP\_\d+\.\d+)\|/) {
    #print STDERR "Skipping $1\n";
    $yp++;
    next SEQFETCH;
  } elsif ( $display_id =~/sp\|(\S+)\|\S+\_$species_name/ ){
    $display_id = $1;
    if ( $prot->desc =~/\s+SV=(\d+)$/ ){
      $display_id.= ".".$1;
    }
    $sp++;
  } elsif ( $display_id =~/tr\|(\S+)\|\S+\_$species_name/ ){
    $display_id = $1;
    # >tr|B0LXP6|B0LXP6_PIG IkB kinase-a OS=Sus scrofa GN=IKK-alpha PE=2 SV=1
    if ( $prot->desc =~/\s+SV=(\d+)$/ ){
      $display_id.= ".".$1;
    }
    $tr++;
  } 
  elsif ( $display_id =~/sp\|(\S+)\|\S+/ ){
    $display_id = $1;
    if ( $prot->desc =~/\s+SV=(\d+)$/ ){
      $display_id.= ".".$1;
    }
    $sp++;
  } elsif ( $display_id =~/tr\|(\S+)\|\S+/ ){
    $display_id = $1;
    if ( $prot->desc =~/\s+SV=(\d+)$/ ){
      $display_id.= ".".$1;
    }
    $tr++;
  } else {
    die "Unrecognised format ".$display_id."\n"; 
    $unrecognised_header++;
  }
  
  if (exists($hadhash{$display_id})) {
    print STDERR "Skipping $display_id as already had one with that name\n";
    $duplicate++;
    next SEQFETCH;
  }

  $hadhash{$display_id} = 1;

  $prot->display_id($display_id);
  $prot->desc("");
  
  
  if($@){
    warn("can't parse sequence for [$description]:\n$@\n");
    $unrecognised_description++;
    next SEQFETCH;
  }

  my $new_prot;

  $prot->seq(uc($prot->seq));

  # write sequence
  $out++;
  $seqout->write_seq($prot);
}

print "\nRead in $in sequences and wrote out $out sequences\n".
      "too_short $too_short\n".
      "Unrecognised header $unrecognised_header\n".
      "unrecognised_description $unrecognised_description\n".
      "duplicate $duplicate\n".
      "np $np\n".
      "xp $xp\n".
      "ap $ap\n".
      "yp $yp\n".
      "sp $sp\n".
      "tr $tr\n" ;
