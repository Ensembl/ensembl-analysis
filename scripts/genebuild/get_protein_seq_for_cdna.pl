use strict;

=head1 NAME

  get_protein_seq_for_cdna.pl

=head1 SYNOPSIS

  get_protein_seq_for_cdna.pl

=head1 DESCRIPTION

  get_protein_seq_for_cdna.pl makes a file of protein 
  fasta sequences, with the corresponding cDNA accession
  followed by the protein accession in the header, e.g.
  
  >NM_XXXXXX.1  NP_XXXXXX.1
  MSALPGSKLSERVRTVGWQISRPYFCHFFPIRITAPPATCSANKGFPELEHARPCPKRCP
  GSISQAIHVGKMAAVQVAASLPCEQPREAPRELSLEQNNGFRRLSARLRALQPDDSTVSR
  ......

  The output sequence file is often used at the BestTargetted
  step in genebuild when gene models from Exonerate cdna2genome 
  analyses need to be considered.

  The script will read in a list of cDNA accessions and
  then will go to Mole to look for the protein accession.
  Sequences will be pfetched and then written to an outfile.

=head1 OPTIONS

  -dbnames
   database names for Mole databases eg. embl_91,refseq_23. list
   database in descending chronological order (i.e. newest comes
   first).
  
  -infile
   path to the list of cDNA accessions 

  -outfile
   path to the fasta file

=cut

use Bio::SeqIO;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::ExternalData::Mole::DBXref;
use Bio::EnsEMBL::ExternalData::Mole::Entry;
use Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor;

$| = 1; # disable buffering

my (
        @dbnames,
        $infile,
        $outfile,
);
my $dbuser;
my $dbhost;
my $dbport;


&GetOptions(
        'dbnames=s'              => \@dbnames,
        'dbuser=s'               => \$dbuser,
        'dbhost=s'               => \$dbhost,
        'dbport=s'               => \$dbport,
        'infile=s'               => \$infile,
        'outfile=s'              => \$outfile,
);

# check commandline
if (!defined($dbhost) || !defined($dbuser) || !scalar(@dbnames)) {
  throw ("Please set dbhost (-hdbost), dbport (-dbport), dbnames (-dbnames) and dbuser (-dbuser)");
}

if (!defined($infile) || !defined($outfile)) {
  throw ("Please enter a file name to read (-infile) and an outfile name (-outfile)");
}

if (scalar(@dbnames)) {
  @dbnames = split(/,/,join(',',@dbnames));
}

# connect to databases
my @dbs;
foreach my $dbname (@dbnames) {
  my $db = Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor->new(
        '-dbname' => $dbname,
        '-host'   => $dbhost,
        '-user'   => $dbuser,
        '-port'   => $dbport,
  );
  push @dbs, $db;
}

#outfile, to print the protein sequences to
my $seqout = new Bio::SeqIO( -file   => ">$outfile",
                             -format => "Fasta"
                           );
 
# read in the cDNA accession list
open(INFILE, "<", "$infile") or die("can t read $infile $!\n");
while (<INFILE>) {
  chomp;
  my $acc = $_;
  my ($protein_id,$protein_seq) = get_protein_from_mole($acc);
  if ($protein_id && $protein_seq) {
    my $to_print = new Bio::Seq(-id => "$acc\t$protein_id", -seq => "$protein_seq");
    $seqout->write_seq($to_print);
  } else {
    print STDERR "Unable to find corresponding protein seq for cDNA $acc\n";
  }
}
close INFILE;

=head2 get_protein_from_mole 

  Arg [1]   : String - accession of cDNA, with version 
  Function  : Reads the protein accession out of Mole 
  Returntype: Array of ($protein_id,$protein_seq). First element a string
              that's the protein accession. Second element a string that's
              the protein sequence.
  Exceptions:
  Example   : my ($protein_id,$protein_seq) = get_protein_from_mole($cdna_accession);

=cut
sub get_protein_from_mole {
  my ($id)  = @_;
  my ($in_db, $entry, $protein_id, $protein_seq);

  if ($id =~ m/\.\d/) {
    #accession has a version
    foreach my $db (@dbs) { 
      $entry = $db->get_EntryAdaptor->fetch_by_accession_version($id);
      $in_db = $db;
      last if defined $entry;
    }
  }
  if (defined $entry) {
   # my ($protein_id, $protein_seq) = fetch_protein_info($entry,$in_db);
    ($protein_id) = fetch_protein_info($entry,$in_db);
  } else {
    warn("Could_not_find_cDNA_in_mole $id\n");
  }
  if (!$protein_id || $protein_id eq '-') {
    warn("Could_not_find_proteinID_for $id\n");
  }
  $protein_seq = pfetch_seq($protein_id);
  if (!$protein_seq) {
    warn("Could_not_find_protein_seq_for $protein_id\n");
  }
  return ($protein_id, $protein_seq); 
}  
  

=head2 fetch_protein_info 

  Arg [1]   : Entry object from Mole database 
  Arg [2]   : Arrayref of Mole databases 
  Function  : Reads the protein_id out of primary_id 
  Returntype: String, protein's accession 
  Exceptions: 
  Example   : my ($protein_acc) = fetch_protein_info($entry, $db); 

=cut
sub fetch_protein_info {
  my ($entry, $db) = @_;
  my $primary_id;

  # Get necessary information...
  my @dbxref_objs = @{$db->get_DBXrefAdaptor->fetch_all_by_entry_id($entry->dbID)}; 

  foreach my $dbxo (@dbxref_objs) {
    if (defined $dbxo->database_id && defined $dbxo->primary_id) {
      if ($dbxo->database_id eq 'EMBL-CDS' || $dbxo->database_id eq 'RefSeq-CDS') {
        $primary_id = $dbxo->primary_id;
      }
    }
  }
  if ($primary_id eq '-') {
    # we can't get the protein id from mole
    # maybe we can get it by doing a full pfetch
    my $command = "pfetch -F ".$entry->accession_version;
    open (OUTFILE, $command." | ") or die "couldn't open ".$command;
    while(<OUTFILE>){
      chomp;
      my $line = $_;
      if ($line =~ /\/protein\_id=\"(.*)\"/) {
        # eg. /protein_id="AAK84156.1"
        $primary_id = $1;
      }
    }
    close(OUTFILE);
  }
#  if ($primary_id eq '-' && $entry->accession_version =~ m/^NM/) { 
#    $primary_id = $entry->accession_version;
#    $primary_id =~ s/M/P/;
#  }
  return ($primary_id);
}

sub pfetch_seq {
  my ($id) = @_;
  my $command = "pfetch -q ".$id;
  my $seq;
  open (OUTFILE, $command." | ") or die "couldn't open ".$command;
  while(<OUTFILE>){
    chomp;
    my $line = $_;
    if ($line ne "no match") {
      $seq .= $line;
    }
  }
  close (OUTFILE);
  return ($seq);
}

