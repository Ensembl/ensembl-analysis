#!/usr/local/bin/perl

#
# Source: schema 20+ ensembl db with a 'new' assembly (eg NCBI35)
# Target: schema 20+ ensembl db with an 'old' assembly (eg NCBI34)
# Action: script copies the new assembly into the db with the old assembly.
# Purpose: to have two assemblies co-exist in the same ensembl db,
#   so you can move annotations between them (for instance...) 
#
# WARNING: This script is not 'polished' - it is provided here so you
# can have a starting point if you want to do this again (I find it easier.
# to copy other people's scripts) but don't assume it will just work for you...
# -- in particular, it will assume your coordinate systems look 'human' - 
# chromosome/clone/contig etc. I made on attempt to generalise this. 

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $source_host    = 'ecs4';
my $source_user    = 'ensro';
my $source_pass    = '';
my $source_dbname  = 'val_homo_sapiens_23_35_ref';
my $source_port    = 3353;

my $target_host    = 'ecs2';
my $target_user    = 'ensro';
my $target_pass    = undef;
#my $target_user    = 'ensadmin';
#my $target_pass    = 'ensembl';
my $target_port    = 3361;
my $target_dbname  = 'vivek_homo_sapiens_core_24_34';

my $ref_host = 'ecs4';
my $ref_user = 'ensro';
my $ref_port = '3353';
my $ref_dbname = 'val_homo_sapiens_23_35_ref';

#my $chr      = '2';
#my $chrstart = 1;
#my $chrend   = 243615958;
#my $version     = 'NCBI34';

my $chr      = undef;
my $chrstart = undef;
my $chrend   = undef;
my $version     = 'NCBI35';

&GetOptions( 
  'source_host:s'    => \$source_host,
  'source_user:s'    => \$source_user,
  'source_pass:s'    => \$source_pass,
  'source_port:s'    => \$source_port,
  'source_dbname:s'  => \$source_dbname,
  'target_host:s'=> \$target_host,
  'target_user:s'=> \$target_user,
  'target_pass:s'=> \$target_pass,
  'target_port:s'=> \$target_port,
  'target_dbname:s'  => \$target_dbname,
  'chr:s'     => \$chr,
  'chrstart:n'=> \$chrstart,
  'chrend:n'  => \$chrend,
  'version:s'  => \$version,
);

my $source_db = 
  new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $source_host,
    -user=>$source_user,
    -pass => $source_pass,
    -port => $source_port,
    -dbname => $source_dbname
  );

my $target_db = 
  new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $target_host,
    -user => $target_user,
    -pass => $target_pass,
    -port => $target_port,
    -dbname => $target_dbname
  );

my $ref_db = 
  new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $ref_host,
    -user => $ref_user,
    -port => $ref_port,
    -dbname => $ref_dbname
  );

# Find the 'contig' coordinate system id in the target db
# so that we can check for duplicates
my $coordinate_system_statement =
  $target_db->prepare(
    qq{
      select coord_system_id from coord_system where name = 'contig'
    }
  );
$coordinate_system_statement->execute;
my $coordinate_system_hashref = $coordinate_system_statement->fetchrow_hashref;
if(not defined $coordinate_system_hashref){
  die 'couldnt find contig coordinate system in target database';
}
my $target_contig_coord_system_id = $coordinate_system_hashref->{'coord_system_id'};

#Find 'chromosome' coordinate system in target db.
$coordinate_system_statement =
  $target_db->prepare(
    qq{
      select coord_system_id from coord_system where name = 'chromosome' and version = '$version'
    }
  );
# print STDERR $coordinate_system_statement->{'Statement'}."\n";
$coordinate_system_statement->execute;
my $hashref = $coordinate_system_statement->fetchrow_hashref;
if(not defined $hashref){
  die 'couldnt find chromosome coordinate system in target database';
}
my $target_chromosome_coord_system_id = $hashref->{'coord_system_id'};

# Find the chromosome of interest by name in the target db.
my $chromosome_statement = 
  $source_db->prepare( 
    qq {
     select seq_region_id from seq_region 
     where seq_region.name='$chr' and coord_system_id = 1
    }
  );

# print STDERR $chromosome_statement->{'Statement'}."\n";
$chromosome_statement->execute;
my $chromosome_hashref = $chromosome_statement->fetchrow_hashref;
if (not defined $chromosome_hashref) {
  die "Couldn't find chromosome $chr in source database\n";
}
my $source_chromosome_id = $chromosome_hashref->{'seq_region_id'};

print STDERR "Source chromosome id = $source_chromosome_id\n";

# Find the chromosome of interest by name in the target db.
$chromosome_statement = 
  $target_db->prepare( 
    qq {
     select seq_region_id from seq_region 
     where seq_region.name='$chr' and
     seq_region.coord_system_id = $target_chromosome_coord_system_id
    }
  );

# print STDERR $chromosome_statement->{'Statement'}."\n";
$chromosome_statement->execute;
$chromosome_hashref = $chromosome_statement->fetchrow_hashref;
if (not defined $chromosome_hashref) {
  die "Couldn't find chromosome $chr in target database\n";
}
my $target_chromosome_id = $chromosome_hashref->{'seq_region_id'};

print STDERR "Target chromosome id = $target_chromosome_id\n";

# First check for duplicate clones used in contigs
my $dup_sth = $source_db->prepare( qq {
  select component_region.name from 
         assembly, seq_region component_region, dna, seq_region assembled_region
  where  assembly.cmp_seq_region_id = assembled_region.seq_region_id and
         assembled_region.name='$chr' and 
         component_region.seq_region_id = assembly.cmp_seq_region_id and 
         assembly.asm_start >= $chrstart and
         assembly.asm_end <=$chrend and
         dna.seq_region_id = component_region.seq_region_id
  order by assembly.cmp_start 
});
        
$dup_sth->execute;

my %clones;
my %duplicated_clones;
my $number_duplicate = 0;
my $hashref;
while ($hashref = $dup_sth->fetchrow_hashref()) {
  my $contigname = $hashref->{'name'};
  my @bits = split /\./,$contigname;

  my $clone = $bits[0] . "." . $bits[1];

  if (exists($clones{$clone})) {
    print "Duplicate clone$clone \n";
    $duplicated_clones{$clone} = 1;
    $number_duplicate++;
  }

  $clones{$clone} = $_;
}

print STDERR "Have found ".scalar(keys %duplicated_clones)." duplicated clones\n";

my $string = 
  qq {
    select 
      assembly.asm_seq_region_id,
      assembly.cmp_seq_region_id,
      assembly.asm_start,
      assembly.asm_end,
      assembly.cmp_start,
      assembly.cmp_end,
      assembly.ori,
      component_region.name
    from 
      assembly, seq_region component_region
    where
      assembly.asm_seq_region_id = $source_chromosome_id and
      assembly.asm_start >=  $chrstart and
      assembly.asm_end <= $chrend and
      assembly.cmp_seq_region_id = component_region.seq_region_id
    order by assembly.asm_start 
  };

print STDERR $string."\n";

my $statement = $source_db->prepare($string);
        
$statement->execute;

my $source_slice_adaptor = $source_db->get_SliceAdaptor();
my $target_slice_adaptor = $target_db->get_SliceAdaptor();
my $ref_slice_adaptor = $ref_db->get_SliceAdaptor();


my $hashref;
my $ndiffseq = 0;
my $totfailed = 0;
my $nnotfound = 0;
my $nwritten = 0;
my $number_duplicate_contig = 0;

while ($hashref = $statement->fetchrow_hashref()) {

  # print STDERR "Processing row: \n";
  my $old_asm_seq_region_id = $hashref->{'asm_seq_region_id'};
  my $old_cmp_seq_region_id = $hashref->{'cmp_seq_region_id'};
  my $asm_start = $hashref->{'asm_start'};
  my $asm_end = $hashref->{'asm_end'};
  my $cmp_start = $hashref->{'cmp_start'};
  my $cmp_end = $hashref->{'cmp_end'};
  my $ori = $hashref->{'ori'};
  my $contigname = $hashref->{'name'};

  my @bits = split /\./,$contigname;
  my $clonename_plus_accession = $bits[0] . "." . $bits[1];

  if (!exists($duplicated_clones{$clonename_plus_accession})) {

    my $contig = $target_slice_adaptor->fetch_by_region(undef, $contigname);
    if($contig){
      print STDERR "Fetched contig: ".$contig->dbID." from target slice\n";
    }
   
    # check for contigs with name without version (AC.V.ST.ED -> AC.ST.ED)
    if (!defined($contig)) {
      $contig = find_contig_by_accession($clonename_plus_accession, $target_slice_adaptor);
    }

    # check for contigs with without start-end (AC.V.ST.ED -> AC.V)
    if (!defined($contig)) {
      $contig = find_contig_with_clone($contigname, $target_db);
    }

    if (defined($contig)) {
      print STDERR "Fetching source contig by slicename: $contigname\n";
      my $source_contig= $source_slice_adaptor->fetch_by_region(undef, $contigname);
      my $source_contig_sequence = $source_contig->subseq($cmp_start,$cmp_end);
      #$contig is the 'target' contig
      my $target_contig_sequence = $contig->subseq($cmp_start,$cmp_end);

      if ($source_contig_sequence ne $target_contig_sequence) {
        print STDERR "start = " . $cmp_start . " end " . $cmp_end . "\n";
        print STDERR "NOTE Contig subseqs sequences different for $clonename_plus_accession\n";
        compare_seqs($source_contig_sequence, $target_contig_sequence);
        $ndiffseq++; 
        $totfailed++;


      } else {

        my $contig_id = $target_slice_adaptor->get_seq_region_id($contig);
        print STDERR "Found contig " . $contig->seq_region_name. "($contig_id)\n";

        my $insert_statement = 
          $target_db->prepare( 
            qq:insert into assembly(asm_seq_region_id,asm_start,asm_end,cmp_seq_region_id,cmp_start,cmp_end,ori) values($target_chromosome_id, $asm_start, $asm_end, $contig_id, $cmp_start, $cmp_end, $ori):
          );

        print STDERR $insert_statement->{Statement} . "\n";
        #$insert_statment->execute;
        $nwritten++;
      }
    } else {
      print STDERR "Didn't find $contigname in target db\n";
      $nnotfound++; 
      $totfailed++;
    } 
  } else {
    print STDERR "Clone with multiple contigs $clonename_plus_accession\n";
    $totfailed++;
    $number_duplicate_contig++;
  }
}

print STDERR "Number of assembly elements written = $nwritten\n";
print  STDERR "Total failures = $totfailed\n";
print  STDERR "No. with different seq = $ndiffseq\n";
print  STDERR "No. not found in target db = $nnotfound\n";
print  STDERR "No. of clones with duplicates = " . scalar(keys(%duplicated_clones)) . " (n contig  = $number_duplicate_contig)\n";


sub compare_seqs {
  my ($seq1, $seq2) = @_;

  $seq1 =~ s/(.{80})/$1\n/g;
  $seq2 =~ s/(.{80})/$1\n/g;

  # print "Chr = $chrstr\n";
  # print "Contig = " . $contigsubstr . "\n";

  if ($seq1 ne $seq2) {
    my $ndiffline = 0;
    my @seq2lines = split /\n/,$seq2;
    my @seq1lines = split /\n/,$seq1;
    for (my $linenum = 0; $linenum<scalar(@seq2lines); $linenum++) {
      if ($seq1lines[$linenum] ne $seq2lines[$linenum]) {
        $ndiffline++;
      }
    }
    print "N diff line = $ndiffline N line = " . scalar(@seq2lines)."\n";
    if ($ndiffline > 0.95*scalar(@seq2lines)) {
#        print "Chr = $chrstr\n";
#        print "Contig = " . $seq1 . "\n";
      for (my $linenum = 0; $linenum<scalar(@seq2lines); $linenum++) {
        if ($seq1lines[$linenum] eq $seq2lines[$linenum]) {
          print "Matched line in very different: $seq1lines[$linenum]\n";
        }
      }
    }
  }
}

sub find_contig_by_accession{
  my ($clonename_plus_accession, $target_slice_adaptor) = @_;
  print STDERR "Trying for contig with $clonename_plus_accession\n";
  my $contig = $target_slice_adaptor->fetch_by_region(undef, $clonename_plus_accession);
  if(defined($contig)){
    print STDERR "Looked for $clonename_plus_accession successful\n";
    return $contig;
  }else{
    print STDERR "Couldnt find contig\n";
    return;
  }
}

sub find_contig_with_clone{
  my ($contigname, $target_db) = @_;
  print STDERR "Missing contig " . $contigname . " - trying with just clone name\n";
  my $new_contigname = $contigname;
  $new_contigname =~ s/\.[0-9]*\.[0-9]*\.[0-9]*$//;
  if($new_contigname){
    my $contig_search_statement = $target_db->prepare("select name from seq_region where name like '$new_contigname%'");
    print STDERR $contig_search_statement->{Statement};
    $contig_search_statement->execute;
    my $hashref2 = $contig_search_statement->fetchrow_hashref;
  
    if (defined($hashref2)) {
      print STDERR "Looking for $new_contigname did find " . $hashref2->{'name'} . "\n";
      my $contig = $target_slice_adaptor->fetch_by_region(undef, $hashref2->{'name'});
      return $contig;
    }else{
      print STDERR "Couldnt find alternate contigname\n";
      return;
    }
  }else{
    print STDERR "Couldnt create an alternate contigname\n";
    return;
  }
}
