#!/usr/local/ensembl/bin/perl -w 

# Dumps markers from the marker table into a file, in the format expected by
# the EPCR Runnable


use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $host   = 'ecs4';
my $user   = 'ensro';
my $pass   = undef;
my $dbname = 'steve_mouse_nt_fix';
my $port   = 3350;
my $outfile  = 'markers.dat';

my $path = 'NCBIM33';

$| = 1;

&GetOptions(
  'host:s'   => \$host,
  'user:s'   => \$user,
  'dbname:s' => \$dbname,
  'port:n'   => \$port,
  'pass:s'   => \$pass,
  'path:s'   => \$path,
  'outfile:s' => \$outfile,
);


# Open database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -pass   => $pass,
  -port   => $port,
  -dbname => $dbname,
  -path   => $path,
);


my $sts = $db->get_MarkerAdaptor->fetch_all;
die "No markers in database" unless @{$sts};


dump_sts_file($outfile, $sts);

sub dump_sts_file {
    my ($dest, $sts) = @_;

    open DEST, "> $dest";
    foreach my $m (@{$sts}) {
        unless (ref $m && $m->isa("Bio::EnsEMBL::Map::Marker")) {
            die "Object not a Bio::EnsEMBL::Map::Marker: [$m]";
        }
        next if $m->max_primer_dist == 0;
        next unless length($m->left_primer) > 0;
        next unless length($m->right_primer) > 0;
        print DEST join("\t",
            $m->dbID,
            $m->left_primer,
            $m->right_primer,
            join("-", $m->min_primer_dist, $m->max_primer_dist),
        ), "\n";
    }
    close DEST;
}

