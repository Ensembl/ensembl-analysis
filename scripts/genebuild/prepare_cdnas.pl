#!/usr/local/ensembl/bin/perl

use strict;
use warnings;

=head1 NAME

  prepare_cdnas.pl

=head1 SYNOPSIS

  prepare_cdnas.pl

=head1 DESCRIPTION

  prepare_cdnas.pl makes a file of fasta sequences 
  and also prints out an annotation file for use
  by exonerate's cdna2genome 

  The script reads in each cDNA, checks the cDNA against
  the Ensembl kill_list DB if required, and clips it using 
  Bio::EnsEMBL::Analysis::Tools::PolyAClipping. It then 
  looks for the cDNA in Mole and finds its coordinates. 

  If the coordinates are not found or can't be parsed, 
  a message is sent to STDERR and the cDNA is not 
  written to the outfile. (Parsing problem often caused by
  the cDNA not having exact start/end position data in the 
  EMBL-CDS DB.) If the cDNA is found, we check that it has
  not been clipped into the CDS and then write the cDNA's 
  fasta to outfile (/path/to/cdna.clip) and its annotation 
  to annotation_file (/path/to/cdna.annot). 

  The annotation file contains additional information for 
  the user to know what has happened during clipping. In order
  for this annotation file to be used by exonerate's 
  cdna2genome model, the annotation file will have to be 
  edited to contain only the first 4 columns of each line
  eg. less /path/to/cdna.annot | awk '{print $1"\t"$2"\t"$3"\t"$4}' 

=head1 OPTIONS

  -dbnames
   database names for Mole databases separated by commas without
   any whitespaces, eg. embl_91,refseq_23
  
  -infile
   path to the original downloaded cDNA fasta file

  -outfile
   path to the clipped cDNA fasta file to be written

  -min_length
   minimum length of cDNA after its polyA tail has been clipped.
   set to be 0 by default.

  -annotation
   path to the annotation file

=head1 EXAMPLE

  perl prepare_cdnas.pl -infile all_DL_cDNAs.fa -outfile clipped_annotated_cDNAs.fa\
  -annotation cDNAs.annotate -min_length 60 -dbnames embl_91,refseq_23


=cut

use Bio::EnsEMBL::Analysis::Tools::PolyAClipping;
use Bio::SeqIO;
use Bio::EnsEMBL::KillList::KillList;
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
        $annotation_file,
);
my $buffer = 10;
my $window_size = 3;
my $dbuser = 'genero';
my $dbhost = 'cbi3';
my $dbport = 3306;
my $min_length = 0;

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
        'min_length:s'           => \$min_length,
);

# check commandline
if (!defined($dbhost) || !defined($dbuser) || !scalar(@dbnames)) {
  throw ("Please set dbhost (-hdbost), dbport (-dbport), dbnames (-dbnames) and dbuser (-dbuser)");
}

if ( !defined($infile) || !defined($annotation_file) || !defined($outfile) ) {
  throw( "Please enter a file name to read (-infile), "
       . "a file name to write to (-outfile) "
       . "and a file name for annotation (-annotation). "
       . "The infile must contain a list of accessions "
       . "whose CDS coordinates you want to find" );
}

@dbnames = split(/,/,join(',',@dbnames));

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

# get kill_list object and a killed_objects hash

my $kill_list_object = Bio::EnsEMBL::KillList::KillList
     ->new(-TYPE => "cDNA");
my %kill_list = %{ $kill_list_object->get_kill_list() };

# open the downloaded cdna file
my $seqin  = new Bio::SeqIO( -file   => "<$infile",
                             -format => "Fasta",
                           );

my $seqout = new Bio::SeqIO( -file   => ">$outfile",
                             -format => "Fasta"
                           );

# Various counter to keep track of cDNAs lost/discarded in the process:

my $total_cDNAs = 0;
my $XM_count = 0;
my $killed_count = 0;
my $parse_problem = 0;
my $not_in_mole = 0;
my $too_short_after_clip = 0;
my $AT_only_count = 0;
my $cdna_written = 0;

# Now work on cDNAs one by one

SEQFETCH:
while ( my $cdna = $seqin->next_seq ) {
  $total_cDNAs ++;
  if ($cdna->display_id =~ /^gi\|\d+\|ref\|(XM_.+)\|/ ) { 
    # if cDNA is XM_ RefSeq (predicted), we don't even bother clipping.
    print STDERR "XM RefSeq mRNA, skipped: $1.\n";
    $XM_count ++;
    next SEQFETCH;
  }
  
  my ($clipped, $clip_end, $num_bases_removed) = clip_if_necessary($cdna, $buffer, $window_size);
  if (defined $clipped) {  # $clipped would have been returned undef if the entire cDNA seq seems to be polyA/T tail/head
    my $id = $clipped->id;
    # next SEQFETCH unless in_mole;
    # check whether in Mole
    if ($kill_list{$id}) {
      print STDERR "$id is in kill list, discarded.\n";
      $killed_count ++;
      next SEQFETCH;
    }
   
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
      if (length($clipped->seq) >= $min_length) {
        $cdna_written ++;
        $seqout->write_seq($clipped);
      } else {
      print STDERR "Clipped sequence for $id is shorter than $min_length bp and is excluded from output.\n";
      $too_short_after_clip ++;
      }
    }
  } else {
    $AT_only_count ++;
  }
}

my $expected_output_nr = $total_cDNAs-$XM_count-$killed_count-$too_short_after_clip-$AT_only_count-$parse_problem-$not_in_mole;
if ($cdna_written != $expected_output_nr) {
  print "BEWARE, NUMBERS DON'T ADD UP!  EXPECTED TO HAVE $expected_output_nr ENTRIES IN OUTPUT FILE, ".
        "BUT YOU ONLY HAVE $cdna_written.\n";
}
 
print "\n";
print "Number of cDNAs processed: $total_cDNAs.\n";
print "XM RefSeq mRNA discarded: $XM_count.\n";
print "cDNAs present in kill_list DB: $killed_count.\n";
print "cDNAs too short after polyA clipped: $too_short_after_clip.\n";
print "Entire cDNA sequence seems to be of polyA or polyT: $AT_only_count.\n";
print "Can't parse the head/tail info in EMBL file: $parse_problem.\n";
print "Can't find the accession in mole DB: $not_in_mole.\n";
print "cDNAs written to output: $cdna_written.\n";


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

      if ( defined $clip_end ) {
        if ( $clip_end eq "tail" ) {
          if ( $length < $end ) {
            print STDERR "Clipped off too many bases from the tail: $id\n";
            $do_substr = $end;
            $string .= "\t| $clip_end | $end | substr_tail";
          } else {
            $string .= "\t| $clip_end | $num_bases_removed | do_nothing";
          }
        } elsif ( $clip_end eq "head" ) {
          if ( $start - $num_bases_removed < 0 ) {
            print STDERR "Clipped off too many bases from the head: $id\n";
            $do_substr = $start - 1;
            $string    = "$id\t$strand\t1\t$cdslength";
            $string   .= "\t| $clip_end | oldstart $start oldend $end | substr_head";
          } else {
            $string  = "$id\t$strand\t" . ( $start - $num_bases_removed ) . "\t$cdslength";
            $string .= "\t| $clip_end | $num_bases_removed | oldstart $start oldend $end";
            #eg. AF067164.1      +       1075    2946    | head | 12 | oldstart 1087 oldend 2958
          }
        }
      }
      open (ANNOTATION, ">>$annotation_file") or die "Cannot open $annotation_file\n";
      print ANNOTATION "$string\n";
      close ANNOTATION;
      $write_cdna = 1;
    } else {
      print STDERR "Parse_problem $id in db ".$in_db->dbc->dbname." \n";
      $parse_problem ++;
    } 
  } else {
    print STDERR "Not_in_Mole $id.\n";
    $not_in_mole ++;
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
