#!/usr/local/ensembl/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Tools::TrackUtils qw(get_date get_id_trembl get_id_swiss get_id_dna get_id_rna get_id_embl get_id_aa );
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor;

my $file;
my $mhost = 'cbi5d';
my $mport = '3306';
my $muser = 'genero';
my $host;
my $port = 3306;
my $user;
my $pass = '';
my $dbname;
my $molecule_type = 'MRNA';

&GetOptions (
            'dbhost|host=s'     => \$host,
            'dbport|port=s'     => \$port,
            'dbuser|user=s'     => \$user,
            'dbpass|pass=s'     => \$pass,
            'dbname=s'          => \$dbname,
            'file=s'       => \$file,
            'molecule_type=s'          => \$molecule_type,
        );
my $adaptor = Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor->new(
    -host => $host,
    -dbname => $dbname,
    -pass => $pass,
    -user => $user,
    -port => $port
);
my $input_seq_adaptor = $adaptor->get_InputSeqAdaptor;
my $eentry_adaptor;
my $nentry_adaptor;
my $pentry_adaptor;
if ($molecule_type eq 'PROTEIN') {
    my $pdbname = connect_and_retrieve_from_db($mhost, 'mm_ini', $mport, 'SELECT database_name FROM ini WHERE current = "yes" AND available = "yes" AND database_category = "refseq"');
    my $uniprot_db = new Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor(
            -host    => $mhost,
            -user    => $muser,
            -port    => $mport,
            -dbname  => $pdbname
            );
    $pentry_adaptor = $uniprot_db->get_EntryAdaptor;
}
else {
    my $edbname = connect_and_retrieve_from_db($mhost, 'mm_ini', $mport, 'SELECT database_name FROM ini WHERE current = "yes" AND available = "yes" AND database_category = "emblrelease"');
    my $ndbname = connect_and_retrieve_from_db($mhost, 'mm_ini', $mport, 'SELECT database_name FROM ini WHERE current = "yes" AND available = "yes" AND database_category = "emblnew"');
    my $embl_db = new Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor(
            -host    => $mhost,
            -user    => $muser,
            -port    => $mport,
            -dbname  => $edbname
            );
    my $emblnew_db = new Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor(
            -host    => $mhost,
            -user    => $muser,
            -port    => $mport,
            -dbname  => $ndbname
            );

    $eentry_adaptor = $embl_db->get_EntryAdaptor;
    $nentry_adaptor = $emblnew_db->get_EntryAdaptor;
}
my $mdbname = connect_and_retrieve_from_db($mhost, 'mm_ini', $mport, 'SELECT database_name FROM ini WHERE current = "yes" AND available = "yes" AND database_category = "refseq"');
my $refseq_db = new Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor(
        -host    => $mhost,
        -user    => $muser,
        -port    => $mport,
        -dbname  => $mdbname
        );
my $rentry_adaptor = $refseq_db->get_EntryAdaptor;
open(IF, $file) || die('Could not open '.$file);
while (my $line = <IF>) {
    next unless $line =~ /^>/;
    chomp $line;
    my ($accession, $id) = $line =~ />((\w+)\.\d+)/;
    my $entry;
    my $date;
    my $ext_id;
    if ($molecule_type eq 'PROTEIN') {
        if ($id =~ /[A-Z]{2}_\d/) {
            $entry = $rentry_adaptor->fetch_by_accession($id);
            if (!$entry) {
                print STDERR 'No Entry for ', $id, "\n";
                $date = "0000-00-00";
            }
            else {
                $date = $entry->last_updated;
            }
            $ext_id = &get_id_aa;
        }
        else {
            my $cmd = 'mfetch -f lau -i acc:'.$id;
            open INF, "$cmd |" || die ('Could not execute mfetch');
            while (<INF>) {
                if (/(UniProt\S+)/) {
                    $date = &get_date($_);
                    $ext_id = ($_ =~ /Swiss-Prot/o) ? &get_id_swiss : &get_id_trembl;
                    last;
                }
            }
            close(INF);
        }
    }
    elsif ($molecule_type eq 'MRNA') {
        if ($id =~ /[A-Z]{2}_\d/) {
            $entry = $rentry_adaptor->fetch_by_accession($id);
            if (!$entry) {
                print STDERR 'No Entry for ', $id, "\n";
                $date = "0000-00-00";
            }
            else {
                $date = $entry->last_updated;
            }
            $ext_id = &get_id_rna;
        }
        else {
            $entry = $eentry_adaptor->fetch_by_accession($id);
            $entry = $nentry_adaptor->fetch_by_accession($id) unless ($entry);
            if (!$entry) {
                print STDERR 'No Entry for ', $id, "\n";
                $date = "0000-00-00";
            }
            else {
                $date = $entry->first_submitted;
            }
            $ext_id = &get_id_embl;
        }
    }
    my $input_seq = Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq->new(
                -hit_name        => $accession,
                -molecule_type   => $molecule_type,
                -submission_date => $date,
                -external_db_id  => $ext_id
            );
    $input_seq_adaptor->store($input_seq);
}

sub connect_and_retrieve_from_db {
    my ($host, $db, $port, $query) = @_;
    my $dbh = DBI->connect('DBI:mysql:database='.$db.';host='.$host.';port='.$port, 'genero', undef, {RaiseError => 1, AutoCommit => 0}) || die("Could not connect to the $db on ".$host." with port $port!!\n");
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my ($result) = $sth->fetchrow_array();
    $sth->finish();
    $dbh->disconnect();
    return $result;
}
