#!/usr/env perl

use strict;
use warnings;

use Getopt::Long;
use Term::ANSIColor;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::KaryotypeBand;


# Connection to the target DB
my $host;
my $port   = 3306;
my $user;
my $pass;
my $dbname;
my $file;
my $write = 0;
my $verbose = 0;
my $help = 0;

my @colors = (
    'ON_RED',
    'ON_GREEN',
    'ON_YELLOW',
    'ON_BLUE',
    'ON_MAGENTA',
    'ON_CYAN',
    'BOLD RED',
    'BOLD GREEN',
    'BOLD YELLOW',
    'BOLD BLUE',
    'BOLD MAGENTA',
    'BOLD CYAN',
    'RED',
    'GREEN',
    'YELLOW',
    'BLUE',
    'MAGENTA',
    'CYAN',
    );

&GetOptions (
            'host=s'   => \$host,
            'port=i'   => \$port,
            'user=s'   => \$user,
            'pass=s'   => \$pass,
            'dbname=s' => \$dbname,
            'file=s'   => \$file,
            'write!'   => \$write,
            'verbose!' => \$verbose,
            'help!'    => \$help,
        );

if ($help) {
    Usage();
    exit(0);
}

die("You need to provide database information for host: $host user: $user password: $pass port: $port dbname: $dbname\n Or the file to load the karyotype bands from: $file\n") unless ($host and $dbname and $user and $file);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
            -host => $host,
            -port => $port,
            -user => $user,
            -pass => $pass,
            -dbname => $dbname,
);

my %seq_regions;
my %density_colors;
open(RF, $file) || die("Could not open $file\n");
while (my $line = <RF>) {
    next if ($line =~ /^\s*#/);
    my ($region, $arm, $band, $start, $end, $density) = $line =~ /^(\S+)\s+(\w)\s+([0-9\.]+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s(\w+)/;
    push(@{$seq_regions{$region}}, [$arm.$band, $density, $start, $end]);
    if ($verbose and !exists $density_colors{$density}) {
        $density_colors{$density} = pop @colors;
    }
}
close(RF) || die("Could not close $file\n");

my $slice_adaptor = $db->get_SliceAdaptor;
my @karyotype_bands;
foreach my $slice (@{$slice_adaptor->fetch_all('toplevel')}) {
    if (exists $seq_regions{$slice->seq_region_name}) {
        foreach my $info_ref (@{$seq_regions{$slice->seq_region_name}}) {
            my $karyotype_band = Bio::EnsEMBL::KaryotypeBand->new(
                -start => $info_ref->[2],
                -end => $info_ref->[3],
                -slice => $slice,
                -name => $info_ref->[0],
                -stain => $info_ref->[1],
                );
            push(@karyotype_bands, $karyotype_band);
        }
    }
}

if ($verbose) {
    foreach my $karyotype_band (@karyotype_bands) {
        print STDOUT $karyotype_band->slice->seq_region_name, "\t", $karyotype_band->start, "\t", $karyotype_band->end, "\t", $karyotype_band->name, "\t", colored( $karyotype_band->stain, $density_colors{$karyotype_band->stain}), "\n";
    }
}

if ($write) {
    my $karyotype_adaptor = $db->get_KaryotypeBandAdaptor();
    $karyotype_adaptor->store(@karyotype_bands);
}

sub Usage {
    print STDOUT <<EOF

$0 -host <host> -user <user> -dname <dname> -file <karyotype file> [-pass <pass>] [-port <port>] [-write] [-verbose] [-help]
  -host    The server where your database is
  -user    The user, preferably with write permissions
  -dname   The name of the database
  -file    The file containing the bands, downloaded from the NCBI: ftp://ftp.ncbi.nlm.nih.gov/pub/gdp/
  -pass    Password of the user
  -port    Port to connect to
  -write   Write the karyotype to your database
  -verbose Print each of your band with a color for each line, similar bands have the same color
  -help    This help

EOF
}
