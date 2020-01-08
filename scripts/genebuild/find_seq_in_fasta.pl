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

#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my ($fasta_file,$id_file);
my $id;
my $prefix;

&GetOptions(
            'id:s'           => \$id,
            'fasta_file:s'         => \$fasta_file, 
            'id_file:s'         => \$id_file, 
             'prefix:s'       => \$prefix,
           );

if (!defined($fasta_file) || (!defined($id) && !defined($prefix) && !defined($id_file))) {
  die "ERROR: Must at least set file (-fasta_file) and full id (-id)\n" .
      "       or prefix (-prefix)\n";
}

my %ids ;
if ( $id_file ) {  
  open(I,"$id_file") || die ( "Cant read file : $id_file\n") ; 
  while(my $line=<I>){ 
    chomp($line); 
    $ids{$line} = 1;  
  } 
} elsif ( $id ) {  
  $ids{$id} = 1;  
} 


my $inputer = Bio::SeqIO->new(-file => "<" . $fasta_file , '-format' => 'Fasta') ;
my $outputer = Bio::SeqIO->new(-file => ">-" , '-format' => 'Fasta') ;

while (my $seq = $inputer->next_seq) {
   #print $seq->id ."\n" ; 
  if (exists $ids{$seq->id}) {
    $outputer->write_seq($seq);
  }
  if (defined($prefix)) {
    if ($seq->id =~ /^$prefix/) {
      $outputer->write_seq($seq);
    }
  }
}
$inputer->close;
$outputer->close;
