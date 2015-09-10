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
  $self->remove_kill_list_object('cdna_update');
  say "Finished removing kill-list objects";

  $self->polyA_clipping('cdna_update.seqs');
  say "Finished clipping polyA tails";
  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

sub remove_kill_list_object {
    my ($self,$newfile) = @_;
    require Bio::EnsEMBL::KillList::KillList;
    my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new( -TYPE => 'cdna_update' );
    my %kill_list = %{ $kill_list_object->get_kill_list() };

    my $DATA_DIR = '/lustre/scratch109/ensembl/dm15/hive_cdna/';

    my $GSS = '/nfs/users/nfs_d/dm15/cvs_checkout_head/ensembl-personal/genebuilders/cDNA_update/gss_acc.txt';
    #my $newfile = '/lustre/scratch109/ensembl/dm15/hive_cdna/cdna_update';

    open( LIST, "<", $GSS ) or die "can't open gss list $GSS";
    my %gss;
    while (<LIST>) {
        my @tmp = split /\s+/, $_;
        $gss{ $tmp[1] } = 1;
    }
    close LIST;

#    # Go through file removing any seqs which appear on the kill list
    local $/ = "\n>";
#    if ( !( -e $DATA_DIR . "/" . $newfile ) ) {
#        print("\tNewfile not here so need to create it first.\n");
#        write_to_file();
#    }

    my $newfile2 = $newfile . ".seqs";
    open( SEQS, "<", $DATA_DIR . "/" . $newfile )
      or die "Can't open seq file $newfile";
    open( OUT, ">", $DATA_DIR . "/" . $newfile2 )
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
    my ($self,$trim_file) = @_;

    my $DATA_DIR = '/lustre/scratch109/ensembl/dm15/hive_cdna/';
    my $POLYA_CLIPPING = '~/enscode/ensembl-pipeline/scripts/EST/new_polyA_clipping.pl';
    # Clip ployA tails
    print("\nPerforming polyA clipping...\n");
    my $newfile3 = $DATA_DIR . "/" . $trim_file. ".clipped";
    my $cmd = "perl " . $POLYA_CLIPPING . " " ;
    $cmd.="-errfile $DATA_DIR/polyA.err ";
    #if ( $MIN_LENGTH ) {
    #   $cmd.="-min_length $MIN_LENGTH ";
    #}
      $cmd .=  $DATA_DIR . "/" . $trim_file . " " . $newfile3;

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

