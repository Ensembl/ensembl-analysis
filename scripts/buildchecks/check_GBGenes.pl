#!/usr/local/ensembl/bin/perl -w

# Script to check the integrity of some or all of the genes in an Ensembl 
# database 

# Maintained by:  Steve Searle (searle@sanger.ac.uk)

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use TranscriptChecker;
use ContigGenesChecker;

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

my @chromosomes;

my $specstart = 1;
my $specend   = undef;

my $dnadbname = "";
my $dnahost   = "";
my $dnaport   = "";

my $exon_dup_check = 0;

my $check_transcripts = 1;
my $schema = 20;

&GetOptions(
            'host:s'           => \$host,
            'user:s'           => \$user,
            'pass:s'           => \$pass,
            'port:n'           => \$port,
            'dbname:s'         => \$dbname,
            'dnahost:s'        => \$dnahost,
            'dnaport:n'        => \$dnaport,
            'dnadbname:s'      => \$dnadbname,
            'path:s'           => \$path,
            'ignorewarnings!'  => \$ignorewarnings,
            'chromosomes:s'    => \@chromosomes,
            'chrstart:n'       => \$specstart,
            'chrend:n'         => \$specend,
            'schema:n'         => \$schema,
            'duplicates!'      => \$exon_dup_check,
            'transcripts!'     => \$check_transcripts,
           );

if (!defined($host) || !defined($dbname)) {
  die "ERROR: Must at least set host (-host), dbname (-dbname)\n" .
      "       (options can also be set in GeneConf.pm)\n";
}

if (scalar(@chromosomes)) {
  @chromosomes = split(/,/,join(',',@chromosomes));
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
  $chrhash = get_chrlengths($db, $path);
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
    $slicename = "chromosome:$path:$chr:$chrstart:$chrend:1";
    print "Slice = $slicename\n";
    $slice = $sa->fetch_by_name($slicename);
  } else {
    $slicename = "$chr:$chrstart-$chrend";
    print "Slice = $slicename\n";
    $slice = $sa->fetch_by_chr_start_end($chr,$chrstart,$chrend);
  }
  
  my $genes = $slice->get_all_Genes;

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
  
    my @trans = @{$gene->get_all_Transcripts()};
    $total_genes++;
  
    my $nwitherror = 0;
    if (scalar(@trans) == 0) {
      print_geneheader($gene);
      print "ERROR: Gene " . $gene->dbID . " has no transcripts\n";
      $total_genes_with_errors++;
      $nwitherror=1; 
    } elsif (scalar(@trans) > $maxtranscripts) {
      print_geneheader($gene);
      print "ERROR: Gene " . $gene->dbID .  
            " has an unexpected large number of transcripts (" .  scalar(@trans) . ")\n";
      $total_genes_with_errors++;
      $nwitherror=1; 
    }
  
    if ($check_transcripts) {
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
            print_geneheader($gene);
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

sub print_geneheader {
  my $gene = shift;

  print "\n++++++++++++++++++++++++++++\n";
  print "Gene " . $gene->dbID . " type " . $gene->type . "\n";
}

sub get_chrlengths {
  my $db   = shift;
  my $type = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my $query = "select seq_region.name, seq_region.length as mce from seq_region,coord_system where" .
              " seq_region.coord_system_id=coord_system.coord_system_id and" .
              " coord_system.version = '" . $type . "' and coord_system.name='chromosome'";

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
    print $hashref->{'name'} . " " . $hashref->{'mce'} . "\n";
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
