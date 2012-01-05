#!/usr/local/ensembl/bin/perl

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;

# Connection to the target DB
my $host   = '';
my $port   = '3306';
my $user   = 'ensadmin';
my $pass   = '';
my $dbname = '';
my $output;
my $verbose;
my %h_mapped;
my %h_count;
my %h_unmapped;
my %h_seen;

&GetOptions (
            'host|dbhost=s'  => \$host,
            'port|dbport=s'  => \$port,
            'user|dbuser=s'  => \$user,
            'pass|dbpass=s'  => \$pass,
            'dbname=s'       => \$dbname,
            'output=s'       => \$output,
            'verbose!'      => \$verbose,
        );
if ($output) {
    open(STDOUT, '>'.$output) || die('Could not open file for writing: '.$output);
}

my $dbe = Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor->new(
    -host => $host,
    -port => $port,
    -user => $user,
    -pass => $pass,
    -dbname => $dbname);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -dbconn => $dbe->dbc);
my $evidencetrack_adaptor = $dbe->get_EvidenceTrackAdaptor;

my $t = localtime;
print STDERR 'Fetching tracks... ', $t, "\n" if ($verbose);
my @a_evidencetracks = @{$evidencetrack_adaptor->fetch_all()};
$t = localtime;
print STDERR 'Get aligned sequences: ', $t, "\n" if ($verbose);
foreach my $evidencetrack (@a_evidencetracks) {
    my $reason_id = $evidencetrack->reason->dbID;
    my $code = $evidencetrack->code;
    my $name = $evidencetrack->name;
    print STDERR $name, "\t", $code, "\n";
    if ($reason_id == 0) {
        print STDERR 'Problem with: ', $evidencetrack->name, ' ', $evidencetrack->code, ' ', $evidencetrack->analysis_run->dbID, ' ', $evidencetrack->input_id, "\n";
    }
    elsif ($reason_id < 100) {
        $h_mapped{$name} = $code;
        $h_count{$code}++ unless (exists $h_seen{$name});
        $h_seen{$name}++;
    }
    else {
#        next if ($code eq "Rejected");
        $h_unmapped{$name}{$code}++;
    }
}
$t = localtime;
print STDERR 'Look at unaligned sequences: ', $t, "\n" if ($verbose);
foreach my $name (keys %h_unmapped) {
    next if (exists $h_mapped{$name});
    foreach my $code (keys %{$h_unmapped{$name}}) {
        $h_count{$code}++;
    }
}

$t = localtime;
print STDERR 'Printing stats! ', $t, "\n" if ($verbose);
foreach my $key (keys %h_count) {
    print STDOUT $key, "\t", $h_count{$key}, "\n";
}
