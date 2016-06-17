#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePreparecDNAs;

use strict;
use warnings;
use feature 'say';


use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my ($self) = shift;

  my $query_hash = $self->param('prepare_seqs');
  my $output_path = $query_hash->{'dest_dir'};
  my $gss_file = $query_hash->{'gss_file'};
  my $embl_file = $query_hash->{'embl_file'};
  my $refseq_file = $query_hash->{'refseq_file'};
  my $polyA_script = $query_hash->{'polyA_script'};
  my $cdna_file = $query_hash->{'cdna_file'};
  my $species = $query_hash->{'species'};

  $self->fix_headers($embl_file,$refseq_file,$output_path,$cdna_file,$species);

  $self->remove_kill_list_object($cdna_file,$output_path,$gss_file);
  say "Finished removing kill-list objects";

  print $output_path, "\n";
  $self->polyA_clipping($polyA_script,$cdna_file,$output_path,$gss_file);
  say "Finished clipping polyA tails";

  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

sub fix_headers {
  my ($self,$embl_file,$refseq_file,$output_path,$cdna_file,$species) = @_;
  local $/ = "\n>";
  open( CF, ">", $output_path . "/" . $cdna_file ) or die "can't create $cdna_file\n";
  #open( EF, "<", $output_path . "/" . $embl_file ) or die "can't read $embl_file\n";
  
  my $header;

  # read EMBL file
  #while ( my $entry = <EF> ) {
  #  # Need this to include the first record when using $/='\n>'
  #  $entry =~ s/^>//;
  #  # Extract & save id
  #  $entry =~ s/^([\w\.\d]+)\s.*\n{1}?/$1\n/;
  #  if ( !$1 ) {
  #     say "\nunmatched id pattern:\n$entry";
  #  }
  #  # Re-write fasta entry
  #  $entry =~ s/\>//g;
  #  print CF '>' . $entry;
  #}
  #close(EF);
  # Just to avoid a missing \n
  #print CF "\n";
  #print("\nRead EMBL file.\n");

  # Read RefSeq file
  open( RF, "<", $output_path . "/" . $refseq_file ) or die "can't read $refseq_file\n";
  while ( my $entry = <RF> ) {
    # Need this to include the first record when using $/='\n>'
    # we're not using 'predicted' XM entries for now
    $entry =~ s/^>//;
    if ( $entry =~ m/^gi.+ref\|(NM_.+)\| Homo sapiens.*/ ) {
      $header = $1;
    } elsif ( $entry =~ m/^gi.+ref\|(NR_.+)\| Homo sapiens.*/ ) {
      $header = $1;
    } else {
      next;
    }
    $entry =~ s/\>//g;
    if ($header) {
      # Reduce header to accession number
      $entry =~ s/^gi.+\n{1}?/$header\n/g;
      print CF '>' . $entry;
    }
  }
  print "read RefSeq file.\n";
  close(RF);
  close(CF);
}

sub remove_kill_list_object {
  my ($self,$newfile,$output_path,$GSS) = @_;
  require Bio::EnsEMBL::KillList::KillList;
  my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new( -TYPE => 'cdna_update' );
  my %kill_list = %{ $kill_list_object->get_kill_list() };

  open( LIST, "<", $GSS ) or die "can't open gss list $GSS";
  my %gss;
  while (<LIST>) {
    my @tmp = split /\s+/, $_;
    $gss{ $tmp[1] } = 1;
  }
  close LIST;

  # Go through file removing any seqs which appear on the kill list
  local $/ = "\n>";

  my $newfile2 = $newfile . ".seqs";
  open( SEQS, "<", $output_path . "/" . $newfile )
    or die "Can't open seq file $newfile";
  open( OUT, ">", $output_path . "/" . $newfile2 )
    or die"Can't open seq file $newfile2";
  while (<SEQS>) {
    s/>//g;

    my @tmp = split /\n/, $_;

    # Store the accession number
    my $acc;
    if ( $tmp[0] =~ /(\w+)\./ ) {
      $acc = $1;
    }
    if ( ( !exists $kill_list{$acc} ) && ( !exists $gss{$acc} ) ) {
      print OUT ">$_";
    }
  }
  local $/ = "\n";
  close OUT;
  close SEQS;
  return $newfile2;
} 

# Clipping the polyA tail
sub polyA_clipping {
  my ($self,$POLYA_CLIPPING,$trim_file,$output_path) = @_;

  # Clip ployA tails
  print("\nPerforming polyA clipping...\n");
  my $newfile3 = $output_path . "/" . $trim_file. ".clipped";
  my $cmd = "perl \$" . $POLYA_CLIPPING . " " ;
  $cmd.="-errfile $output_path/polyA.err ";
  #if ( $MIN_LENGTH ) {
  #   $cmd.="-min_length $MIN_LENGTH ";
  #}
  $cmd .=  $output_path . "/" . $trim_file . " " . $newfile3;

  $cmd = 'bsub -I -q yesterday -M1000 -R"select[mem>1000] rusage[mem=1000]" "'.$cmd.'"';
  print $cmd, "\n";

  system($cmd);
  #  die"Couldn't clip file.$@\n";
  #}

  # Split fasta files, store into CHUNKDIR
#    print("Splitting fasta file.\n");
#    $cmd = "$FASTA_SPLIT $newfile3 $CHUNK $chunkDIR";
#    if ( system($cmd) ) {
#        die "Couldn't split file.$@\n";
#    }

    # Isolate biggest sequences
#    check_chunksizes();

#    print "\nChopped up the file.\n";
}

1;

