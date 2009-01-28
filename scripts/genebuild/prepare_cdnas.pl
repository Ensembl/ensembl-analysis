use strict;

=head1 NAME

  prepare_cdnas.pl

=head1 SYNOPSIS

  prepare_cdnas.pl

=head1 DESCRIPTION

  prepare_cdnas.pl makes a file of fasta sequences 
  and also prints out an annotation file for use
  by exonerate's cdna2genome 

  The script reads in each cDNA and clips it using 
  Bio::EnsEMBL::Analysis::Tools::PolyAClipping. It then 
  looks for the cDNA in Mole and finds its coordinates. 
  If the coordinates are not found or can't be parsed, 
  a message is sent to STDERR and the cDNA is not 
  written to the outfile. If the cDNA is found, we check 
  that it has not been clipped into the CDS and then 
  write the cDNA's fasta to outfile (/path/to/cdna.clip)
  and its annotation to annotation_file 
  (/path/to/cdna.annot). This annotation file contains 
  additional information for the use to know what has 
  happened during clipping. In order for this annotation 
  file to be used by exonerate's cdna2genome model, the 
  annotation file with have to be edited to contain only 
  the first 4 columns of each line
  eg. less /path/to/cdna.annot | awk '{print $1"\t"$2"\t"$3"\t"$4}' 

=head1 OPTIONS

  -dbnames
   database names for Mole databases eg. embl_91,refseq_23
  
  -infile
   path to the original downloaded cDNA fasta file

  -outfile
   path to the clipped cDNA fasta file to be written

  -annotation
   path to the annotation file

=cut

use Bio::EnsEMBL::Analysis::Tools::PolyAClipping;
use Bio::SeqIO;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Mole::DBXref;
use Bio::EnsEMBL::Mole::Entry;
use Bio::EnsEMBL::Mole::DBSQL::DBAdaptor;

$| = 1; # disable buffering

my (
        @dbnames,
        $infile,
        $outfile,
        $annotation_file,
);
my $buffer = 10;
my $window_size = 3;
my $dbuser = 'genero';
my $dbhost = 'cbi3';
my $dbport = 3306;

@dbnames = ('embl_91','refseq_23', 'refseq_24', 'embl_90', 'embl_89');

&GetOptions(
        'dbnames=s'              => \@dbnames,
        'dbuser=s'               => \$dbuser,
        'dbhost=s'               => \$dbhost,
        'dbport=s'               => \$dbport,
        'infile=s'               => \$infile,
        'outfile=s'              => \$outfile,
        'annotation=s'           => \$annotation_file,
        'buffer=s'               => \$buffer,
        'window=s'               => \$window_size,
);

# check commandline
if (!defined($dbhost) || !defined($dbuser) || !scalar(@dbnames)) {
  throw ("Please set dbhost (-hdbost), dbport (-dbport), dbnames (-dbnames) and dbuser (-dbuser)");
}

if (!defined($infile) || !defined($annotation_file) || !defined($outfile)) {
  throw ("Please enter a file name to read (-infile), a file name to write to (-outfile) and a file name for annotation (-annotation). The infile must contain a list of accessions whose CDS coordinates you want to find");
}

if (scalar(@dbnames)) {
  @dbnames = split(/,/,join(',',@dbnames));
}

# connect to databases
my @dbs;
foreach my $dbname (@dbnames) {
  my $db = Bio::EnsEMBL::Mole::DBSQL::DBAdaptor->new(
        '-dbname' => $dbname,
        '-host'   => $dbhost,
        '-user'   => $dbuser,
        '-port'   => $dbport,
  );
  push @dbs, $db;
}

# open the downloaded cdna file
my $seqin  = new Bio::SeqIO( -file   => "<$infile",
                             -format => "Fasta",
                           );

my $seqout = new Bio::SeqIO( -file   => ">$outfile",
                             -format => "Fasta"
                           );
 
SEQFETCH:
while ( my $cdna = $seqin->next_seq ) {
  my ($clipped, $clip_end, $num_bases_removed) = clip_if_necessary($cdna, $buffer, $window_size);
  my $id = $clipped->id;
  #  next unless in_mole;
  # check whether in Mole

  my ($found, $substr) = in_mole($id, $clip_end, $num_bases_removed, $clipped->length);
  if ($substr) {
    # the clipping cut off too much. We must get it back.
    if ($clip_end eq "head") {
      $clipped->seq(substr($cdna->seq, $substr));
    } elsif ($clip_end eq "tail") {
      $clipped->seq(substr($cdna->seq, 0, $substr));
    }
  }
  if ($found) {  
    $seqout->write_seq($clipped);
  } else {
   # warning("Not in Mole $id.");
  }
}


=head2 in_mole 

  Arg [1]   : String - accession of cDNA, with version 
  Arg [2]   : String - "head" or "tail" depending on which end was clipped
  Arg [3]   : Integer - number of bases removed by PolyAClipping
  Arg [4]   : Integer - length of clipped cDNA sequence
  Function  : Reads the coordinates out of tertiary_id
  Returntype: Array of ($write_cdna, $do_substr). First element boolean whether
              cDNA should be written to file or not. Second element an Integer
              that signals whether we need to adjust how much was clipped 
  Exceptions:
  Example   : my ($found, $substr) = in_mole($id, $clip_end, $num_bases_removed, $clipped->length);

=cut
sub in_mole {
  my ($id, $clip_end, $num_bases_removed, $length)  = @_;
  my $entry;
  my $write_cdna = 0;
  my $do_substr;

  my $in_db;
  if ($id =~ m/\.\d/) {
    #accession has a version
    foreach my $db (@dbs) { 
      $entry = $db->get_EntryAdaptor->fetch_by_accession_version($id);
      $in_db = $db;
      last if defined $entry;
    }
  }
  if (defined $entry) {
    my ($strand, $start, $end, $coords) = fetch_coords($entry, $in_db); 
    if (defined $strand && defined $start && defined $end && defined $coords) {
      my $cdslength = $end - $start + 1;
      my $string = "$id\t$strand\t$start\t$cdslength";
      if ($clip_end eq "tail") {
        if ($length < $end) {
          warning("Clipped off too many bases from the tail: $id");
          $do_substr = $end;
          $string .= "\t| $clip_end | $end | substr_tail";
        } else {
          $string .= "\t| $clip_end | $num_bases_removed | do_nothing";
        }
      } elsif ($clip_end eq "head") {
        if ($start-$num_bases_removed < 0) {
          warning("Clipped off too many bases from the head: $id");
          $do_substr = $start-1;
          $string = "$id\t$strand\t1\t$cdslength"; 
          $string .= "\t| $clip_end | oldstart $start oldend $end | substr_head";
        } else {
          $string = "$id\t$strand\t".($start-$num_bases_removed)."\t$cdslength";
          $string .= "\t| $clip_end | $num_bases_removed | oldstart $start oldend $end";
          #eg. AF067164.1      +       1075    2946    | head | 12 | oldstart 1087 oldend 2958
        }
      }
      open (ANNOTATION, ">>$annotation_file") or die "Cannot open $annotation_file\n";
      print ANNOTATION "$string\n";
      close ANNOTATION;
      $write_cdna = 1;
    } else {
      print STDERR "Parse_problem $id ($strand\t$start\t$end\t$coords) db ".$in_db->dbc->dbname." \n";
    } 
  } else {
    print STDERR "Not_in_mole $id\n";
  }
  return ($write_cdna, $do_substr); 
}  
  

=head2 fetch_coords 

  Arg [1]   : Entry object from Mole database 
  Arg [2]   : Arrayref of Mole databases 
  Function  : Reads the coordinates out of tertiary_id  
  Returntype: Array of strand, start, end, unprocessed tertiary_id 
  Exceptions: 
  Example   : my ($strand, $start, $end, $coords) = fetch_coords($entry, \@dbs); 

=cut
sub fetch_coords {
  my ($entry, $db) = @_;
  my $strand;
  my $start;
  my $end;
  my $tertiary_id;

  # Get necessary information...
  my @dbxref_objs = @{$db->get_DBXrefAdaptor->fetch_all_by_entry_id($entry->dbID)}; 

  foreach my $dbxo (@dbxref_objs) {
    if (defined $dbxo->database_id && defined $dbxo->tertiary_id) {
      if ($dbxo->database_id eq 'EMBL-CDS' || $dbxo->database_id eq 'RefSeq-CDS') {
        # The entry will look like this: 245..1273 
        $tertiary_id = $dbxo->tertiary_id;
        if ($dbxo->tertiary_id =~ m/^complement\((\d+)\.\.(\d+)\)/i) {
          $strand = "-";
          $start = $1;
          $end = $2;
        } elsif ($dbxo->tertiary_id =~ m/^(\d+)\.\.(\d+)/) { 
          $strand = "+";
          $start = $1;
          $end = $2; 
        } 
      }
    }
  }
  return ($strand, $start, $end, $tertiary_id);
}
