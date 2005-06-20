#!/usr/local/ensembl/bin/perl -w

# Script to check the integrity of some or all of the genes in an Ensembl 
# database 

# any questions please send to ensembl-dev@ebi.ac.uk

=head1 NAME

check_GBGenes.pl

=head1 SYNOPSIS

perl check_GBGenes.pl -chromosome 1 -coordsystem chromosome -transcripts

This will run the transcript checks on chromosome 1 on the final database and
reference database specified in Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases

=head1 DESCRIPTION

This script runs a series of tests on genes in the database it is pointed
too. These tests are described in more detail in the various perl
modules in this directory

=head1 OPTIONS

  -dbhost database host
  -dbuser database user
  -dbpass database password
  -dbport database port
  -dbname database name
  -dnahost host for dna database
  -dnaport port for dna database
  -dnadbname name for dna database

  These settings by default are taken from:
  Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases

  -ignorewarnings flag to specify whether to ignore the warnings in the 
                  code

  This is take by default from:
  Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder


  -chromosome name of seq region to get genes from
  -coordsystem name of the coordinate system the seq region belongs to

  these are both obligatory options

  -chrstart the start coordinate of the piece of seq region to fetch
  -chrend the end coordinate of the piece of seq region to fetch

  without these the whole seq region if fetched

  -duplicates check for duplicate exons (off by default)
  -transcripts check the transcripts (on by default)

  -help print the perl docs

=head1 EXAMPLES

perl check_GBGenes.pl -chromosome 1 -coordsystem chromosome -transcripts

runs the standard checks on chromosomes using the databases defined in
configuration


perl check_GBGenes.pl -chromosome 3 -coordsystem chromosome -duplicates
-host yourhost -user youruser -port 3306 -dbname yourdb -dnahost yourhost
-dnaport 3306 -dnadbname yourdnadb

runs the standard checks and checks for duplicate exons on the database
specified on the commandline

=head1 NOTES

Note this script uses several settings from the 
Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder

Also this scripts uses several modules which can be found in this directory
so you need to make sure you run the script in this directory or 
alternatively put this directory in your PERL5LIB

=cut

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use TranscriptChecker;
use ContigGenesChecker;
use GeneChecker;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder qw (
					 GB_MAXSHORTINTRONLEN
					 GB_MINSHORTINTRONLEN
					 GB_MINLONGINTRONLEN
					 GB_MAX_EXONSTRANSCRIPT
					 GB_MAXSHORTEXONLEN
					 GB_MINSHORTEXONLEN
					 GB_MINLONGEXONLEN 
					 GB_MAXTRANSCRIPTS
					 GB_MINTRANSLATIONLEN
					 GB_IGNOREWARNINGS
					 GB_MAXGENELEN
					 );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
					 GB_FINALDBHOST
					 GB_FINALDBNAME
					 GB_FINALDBUSER
					 GB_FINALDBPASS
					 GB_FINALDBPORT
					 );
$| = 1;

my $host = $GB_FINALDBHOST || undef;
my $dbname = $GB_FINALDBNAME || undef;
my $user = $GB_FINALDBUSER || 'ensro';
my $pass = $GB_FINALDBPASS || ''; 
my $port = $GB_FINALDBPORT || ''; 

# default path comes out of database
my $path = undef; 

my $maxshortintronlen  =  50;
if (defined($GB_MAXSHORTINTRONLEN)) { $maxshortintronlen  = $GB_MAXSHORTINTRONLEN; }

my $minshortintronlen  =  3;
if (defined($GB_MINSHORTINTRONLEN)) { $minshortintronlen  = $GB_MINSHORTINTRONLEN; }

my $minlongintronlen   =  100000;
if (defined($GB_MINLONGINTRONLEN )) { $minlongintronlen   = $GB_MINLONGINTRONLEN; }

my $maxexonstranscript =  150;
if (defined($GB_MAX_EXONSTRANSCRIPT)) { $maxexonstranscript = $GB_MAX_EXONSTRANSCRIPT; }

my $maxshortexonlen    =  10;
if (defined($GB_MAXSHORTEXONLEN)) { $maxshortexonlen    = $GB_MAXSHORTEXONLEN; }

my $minshortexonlen    =  3;
if (defined($GB_MINSHORTEXONLEN)) { $minshortexonlen    = $GB_MINSHORTEXONLEN; }

my $minlongexonlen     =  50000;
if (defined($GB_MINLONGEXONLEN )) { $minlongexonlen     = $GB_MINLONGEXONLEN; }

my $maxtranscripts     =  10; 
if (defined($GB_MAXTRANSCRIPTS)) { $maxtranscripts     = $GB_MAXTRANSCRIPTS; }

my $mintranslationlen  =  10; 
if (defined($GB_MINTRANSLATIONLEN)) { $mintranslationlen  = $GB_MINTRANSLATIONLEN; }

my $ignorewarnings     =  0; 
if (defined($GB_IGNOREWARNINGS)) { $ignorewarnings     = $GB_IGNOREWARNINGS; }

my $maxgenelen     =  2_000_000; 
if (defined($GB_MAXGENELEN)) { $maxgenelen     = $GB_MAXGENELEN; }

my @chromosomes;

my $specstart = 1;
my $specend   = undef;

my $dnadbname = "";
my $dnahost   = "";
my $dnaport   = "";

my $exon_dup_check = 0;

my $check_transcripts = 1;
my $schema = 20;
my $coordsystem = 'chromosome';
my $help;
my @genetypes;

&GetOptions(
            'dbhost|host:s'    => \$host,
            'dbuser|user:s'    => \$user,
            'dbpass|pass:s'    => \$pass,
            'dbport|port:n'    => \$port,
            'dbname:s'         => \$dbname,
            'dnahost:s'        => \$dnahost,
            'dnaport:n'        => \$dnaport,
            'dnadbname:s'      => \$dnadbname,
            'path:s'           => \$path,
            'ignorewarnings!'  => \$ignorewarnings,
            'chromosomes:s'    => \@chromosomes,
            'coordsystem:s'    => \$coordsystem,
            'chrstart:n'       => \$specstart,
            'chrend:n'         => \$specend,
            'schema:n'         => \$schema,
            'genetypes:s'      => \@genetypes,
            'duplicates!'      => \$exon_dup_check,
            'transcripts!'     => \$check_transcripts,
            'help!'            => \$help,
           ) or perldocs("Failed to get options");

if (!defined($host) || !defined($dbname)) {
  die "ERROR: Must at least set host (-host), dbname (-dbname)\n" .
      "       (options can also be set ".
        "Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases)\n";
}

if($help){
  perldocs();
}

if (scalar(@chromosomes)) {
  @chromosomes = split(/,/,join(',',@chromosomes));
}

if (scalar(@genetypes)) {
  @genetypes = split(/,/,join(',',@genetypes));
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
                                            -user => $user,
                                            -port => $port,
                                            -dbname => $dbname,
                                            -pass => $pass);
                                            

if ($path) {
  $db->assembly_type($path);
}

if ($dnadbname ne "") {

  if ($dnahost eq "") {
    $dnahost = $host;
  }

  my $dnadbase = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dnahost,
                                                    -user   => $user,
                                                    -port   => $dnaport,
                                                    -dbname => $dnadbname,
                                                    -pass   => $pass,
                                                   );
  $db->dnadb($dnadbase);
}

my $sa = $db->get_SliceAdaptor();

#Not practical to do any other way
if ($exon_dup_check) {
  print "Performing exon duplicate check for ALL exons\n";
  find_duplicate_exons($db);
  print "Done duplicate check\n";
}

my $chrhash;

if ($schema == 20) {
  $chrhash = get_chrlengths($db, $path,$coordsystem);
} else {
  $chrhash = get_chrlengths_19($db, $path);
}

#filter to specified chromosome names only 
if (scalar(@chromosomes)) {
  foreach my $chr (@chromosomes) {
    my $found = 0;
    foreach my $chr_from_hash (keys %$chrhash) {
      if ($chr_from_hash =~ /^${chr}$/) {
        $found = 1;
        last;
      }
    }
    if (!$found) {
      die "Didn't find chromosome named $chr in database $dbname\n";
    }
  }
  HASH: foreach my $chr_from_hash (keys %$chrhash) {
    foreach my $chr (@chromosomes) {
      if ($chr_from_hash =~ /^${chr}$/) {next HASH;}
    }
    delete($chrhash->{$chr_from_hash});
  }
}

#print "Start $specstart End $specend\n";


my @failed_transcripts;

my $total_transcripts_with_errors = 0;
my $total_genes_with_errors = 0;
my $total_genes = 0;
my $total_transcripts = 0;

if (!defined($path)) {
  $path = $db->assembly_type;
}

# Begin testing genes
foreach my $chr (sort bychrnum keys %$chrhash) {

  my $chrstart = $specstart;
  my $chrend = (defined ($specend) && $specend < $chrhash->{$chr}) ? $specend :
               $chrhash->{$chr};

  my $slice;
  my $slicename;
  if ($schema == 20) {
    $slicename = "$coordsystem:$path:$chr:$chrstart:$chrend:1";
    print "Slice = $slicename\n";
    $slice = $sa->fetch_by_name($slicename);
  } else {
    $slicename = "$chr:$chrstart-$chrend";
    print "Slice = $slicename\n";
    $slice = $sa->fetch_by_chr_start_end($chr,$chrstart,$chrend);
  }
  
  my $genes;

  if (!scalar(@genetypes)) {
    $genes = $slice->get_all_Genes;
  } else {

    foreach my $type (@genetypes) {
      push @$genes, @{$slice->get_all_Genes_by_type($type)};
    }
  }

  my $cgc = new ContigGenesChecker(
                                     -slice          => $slice,
                                     -genes          => $genes,
                                     -ignorewarnings => $ignorewarnings,
                                     -adaptor        => $db, 
                                     );
  $cgc->check;

  if ($cgc->has_Errors()) {
    print "--------------------------------------\n";
    print "Slice with errors: $slicename\n";
    $cgc->output;
  }

  GENE: foreach my $gene (@$genes) {
    $total_genes++;
  
    my $gc = new GeneChecker(-gene          => $gene,
                             -maxtransgene  => $maxtranscripts,
                             -maxgenelen    => $maxgenelen,
                             -ignorewarnings=> $ignorewarnings,
                             -adaptor       => $db, 
                             -slice         => $slice,
                            );
    $gc->check;
    my $nwitherror = 0;

    if ($gc->has_Errors()) {
      $gc->output;

      $total_genes_with_errors++;
      $nwitherror=1; 
    }
  

    if ($check_transcripts) {
      my @trans = @{$gene->get_all_Transcripts()};

      TRANSCRIPT: foreach my $transcript (@trans) {
        $total_transcripts++;
        my $tc = new TranscriptChecker(
                                       -transcript         => $transcript,
                                       -minshortintronlen  => $minshortintronlen,
                                       -maxshortintronlen  => $maxshortintronlen,
                                       -minlongintronlen   => $minlongintronlen,
                                       -minshortexonlen    => $minshortexonlen,
                                       -maxshortexonlen    => $maxshortexonlen,
                                       -minlongexonlen     => $minlongexonlen,
                                       -mintranslationlen  => $mintranslationlen,
                                       -maxexonstranscript => $maxexonstranscript,
                                       -ignorewarnings     => $ignorewarnings,
                                       -adaptor            => $db, 
                                       -slice              => $slice,
                                       );
        $tc->check;
  
        if ($tc->has_Errors()) {
          $total_transcripts_with_errors++;
          if (!$nwitherror) {
            $gc->output;
            $total_genes_with_errors++;
          }
          $tc->output;
          # Don't store for now!!!  push @failed_transcripts,$tc;
          $nwitherror++;
        }
      }
    }
  }
}

print "Summary:\n";
print "Number of genes checked           = $total_genes\n";
print "Number of transcripts checked     = $total_transcripts\n\n";
print "Number of transcripts with errors = $total_transcripts_with_errors\n";
print "Number of genes with errors       = $total_genes_with_errors\n\n";

#End of main

sub get_chrlengths {
  my $db   = shift;
  my $type = shift;
  my $coordsystem = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my $query = "select seq_region.name, seq_region.length as mce from seq_region,coord_system where" .
              " seq_region.coord_system_id=coord_system.coord_system_id and" .
              " coord_system.version = '" . $type . "' and coord_system.name='" . $coordsystem ."'";

  my $sth = $db->prepare($query);

  $sth->execute;

  my %chrhash;

  my $hashref;
  while (($hashref = $sth->fetchrow_hashref) && defined($hashref)) {
    $chrhash{$hashref->{'name'}} = $hashref->{mce};
#    print $hashref->{'name'} . " " . $hashref->{'mce'} . "\n";
  }
  return \%chrhash;
}

sub get_chrlengths_19 {
  my $db   = shift;
  my $type = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my $query = "select chromosome.name,max(chr_end) as mce from assembly,chromosome where" .
              " assembly.chromosome_id=chromosome.chromosome_id and" .
              " assembly.type = '" . $type . "' group by assembly.chromosome_id";

  my $sth = $db->prepare($query);

  $sth->execute;

  my %chrhash;

  my $hashref;
  while (($hashref = $sth->fetchrow_hashref) && defined($hashref)) {
    $chrhash{$hashref->{'name'}} = $hashref->{mce};
    # print $hashref->{'name'} . " " . $hashref->{'mce'} . "\n";
  }
  return \%chrhash;
}



sub bychrnum {

  my @awords = split /_/,$a;
  my @bwords = split /_/,$b;

  my $anum = $awords[0];
  my $bnum = $bwords[0];

#  if ($anum !~ /^chr/ || $bnum !~ /^chr/) {
#    die "Chr name doesn't begin with chr for $a or $b";
#  }
   
  $anum =~ s/chr//;
  $bnum =~ s/chr//;

  if ($anum !~ /^[0-9]*$/) {
    if ($bnum !~ /^[0-9]*$/) {
      return $anum cmp $bnum;
    } else {
      return 1;
    }
  }
  if ($bnum !~ /^[0-9]*$/) {
    return -1;
  }

  if ($anum <=> $bnum) {
    return $anum <=> $bnum;
  } else {
    if ($#awords == 0) {
      return -1;
    } elsif ($#bwords == 0) {
      return 1;
    } else {
      return $awords[1] cmp $bwords[1];
    }
  }
}

sub find_duplicate_exons {
  my $db = shift;
  
  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "find_duplicate_exons should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my $q = qq( SELECT e1.exon_id, e2.exon_id 
              FROM exon e1, exon e2 
              WHERE e1.exon_id<e2.exon_id AND e1.seq_region_start=e2.seq_region_start AND 
                    e1.seq_region_end=e2.seq_region_end AND e1.seq_region_id=e2.seq_region_id AND 
                    e1.seq_region_strand=e2.seq_region_strand AND e1.phase=e2.phase AND
                    e1.end_phase=e2.end_phase
              ORDER BY e1.exon_id 
            ); 
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
 
  my $ndup = 0;
  while( my ($exon1_id, $exon2_id) = $sth->fetchrow_array) {
    print "ERROR: Exon duplicate pair: $exon1_id and $exon2_id\n"; 
    $ndup++;
  }
  print "Total number of duplicate pairs = $ndup\n";
}

sub perldocs{
  my ($msg) = @_;
  print $msg."\n" if($msg);
  exec('perldoc', $0);
  exit(0);
}
