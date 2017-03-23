#!/usr/env perl

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;


# Connection to the target DB
my $host   = '';
my $port   = '3306';
my $user   = 'ensadmin';
my $pass   = '';
my $dbname = '';
my $write;

&GetOptions (
            'host=s'       => \$host,
            'port=s'       => \$port,
            'user=s'       => \$user,
            'pass=s'       => \$pass,
            'dbname=s'     => \$dbname,
            'write!'       => \$write,
        );

sub check_slice {
    my ($slice) = @_;

    my @genes_to_deletes;
    my @duplicates_genes;
    my @genes = sort ($a->seq_region_start <=> $b->seq_region_start) @{$slice->get_all_Genes};
    for (my $i = 0; $i+1 < scalar(@genes); $i++) {
        if (scalar(@{$genes[$i]->get_all_Transcripts}) == 0) {
            push(@genes_to_deletes, $genes[$i]);
            next;
        }
        elsif (scalar(@{$genes[$i]->get_all_Transcripts}) == 0) {
            push(@genes_to_deletes, $genes[$i]);
            next;
        }
        if ($genes[$i]->seq_region_start == $genes[$i+1]->seq_region_start and $genes[$i]->seq_region_end == $genes[$i+1]->seq_region_end) {
            push(@duplicates_genes, $genes[$i], $genes[$i+1]);
        }
    }
}
