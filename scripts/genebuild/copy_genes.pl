#!/usr/local/ensembl/bin/perl

=head1 NAME

  copy_genes.pl

=head1 DESCRIPTION

 This script takes database options and a file of gene ids or
 stable_ids and copies them between two databases. It can, if asked,
 split multi-transcript genes into single genes.

 When using the in and out_config_name options
 it reads the equivalent database from the
 Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases file.

=head1 OPTIONS

  -all               This option will copy all genes from $sourcedbname to the
                     output database

  -split             This option will split multi-transcript genes into
                     single-transcript genes

  -logic             This option will change the analysis of the genes being
                     written to the output database to an analysis
                     with the specified logic_name

  -remove_xrefs      Using this flag will remove stable IDs from
                     genes/transcripts/translations

  -remove_stable_ids Using this flag will remove stable_ids from
                     genes/transcripts/translations/exons

  -transform_to      Transforms the genes from one coordinate system to the given one
                     cooord_system_name:version 

  -stable_id         Flag for indicating that input file contains stable IDs and
                     not gene IDs

=head1 EXAMPLE

  perl copy_genes.pl -in_config_name COALESCER_DB -out_config_name UTR_DB \
    -file gene_ids_to_copy

  or

  perl copy_genes.pl -sourcehost host -sourceuser ensro -sourcedbname est_db -outhost host1 \
    -outuser user -outpass **** -outdbname utr_db -file gene_ids_to_copy -stable_id

=cut

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw ( throw warning ) ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
use Getopt::Long;
use Bio::EnsEMBL::Analysis::Tools::Utilities;

my $sourcehost      = '';
my $sourceuser      = '';
my $sourcedbname    = '';
my $sourceport      = 3306;

my $outhost   = '';
my $outuser   = '';
my $outpass   = '';
my $outdbname = undef;
my $outport   = 3306;

my $dnahost   = '';
my $dnauser   = 'ensro';
my $dnapass   = '';
my $dnadbname = undef;
my $dnaport   = 3306;

my $in_config_name;
my $out_config_name;

my $logic;
my $split = 0;
my $infile;
my $all;
my $remove_xrefs;
my $remove_stable_ids;
my $transform_to;
my $stable_id;
my $verbose; 
GetOptions( 'sourcehost:s'        => \$sourcehost,
            'sourceuser:s'        => \$sourceuser,
            'sourcedbname:s'      => \$sourcedbname,
            'sourceport:n'        => \$sourceport,
            'in_config_name:s'    => \$in_config_name,
            'out_config_name:s'   => \$out_config_name,
            'outhost:s'           => \$outhost,
            'outuser:s'           => \$outuser,
            'outpass:s'           => \$outpass,
            'outdbname:s'         => \$outdbname,
            'outport:n'           => \$outport,
            'dnahost:s'           => \$dnahost,
            'dnauser:s'           => \$dnauser,
            'dnadbname:s'         => \$dnadbname,
            'dnaport:n'           => \$dnaport,
            'logic:s'             => \$logic,
            'split!'              => \$split,
            'all'                 => \$all,
            'remove_xrefs'        => \$remove_xrefs,
            'remove_stable_ids'   => \$remove_stable_ids,
            'transform_to:s'      => \$transform_to,
            'verbose'            => \$verbose,
            'stable_id'           => \$stable_id,
            'file:s'              => \$infile );


my $transform_to_version;
if($transform_to){ 
  if ( $transform_to=~m/:/) { 
    ( $transform_to,  $transform_to_version ) = split /:/,$transform_to ; 
  }
}

if ($all && $infile) {
  throw("Specify either -all or -infile");
} elsif (!$all && !$infile) {
  throw("Specify -all (to copy all genes from $sourcedbname) or -infile (to copy a list of gene_ids)");
}

#print "Connecting to ".$sourcedbname." at ".$sourcehost." ".$sourceuser."\n";
my $sourcedb;
if ($in_config_name) {
  $sourcedb = get_db_adaptor_by_string($in_config_name);
} else {

  $sourcedb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $sourcehost,
                                        -user   => $sourceuser,
                                        -port   => $sourceport,
                                        -dbname => $sourcedbname );
}

my $dnadb;
if ( !$dnadbname ) {
  my $dna_query = q(SELECT count(1) FROM dna);
  my $dna_sth = $sourcedb->dbc->prepare($dna_query);
  $dna_sth->execute();
  my $dna_count = $dna_sth->fetchrow;

  if ( !$dna_count ) {
    croak( "\nYour database does not contain DNA. "
         . "Please provide a database with genomic sequences.\n");
  }
} else {
  $dnadb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                        -user   => $dnauser,
                                        -port   => $dnaport,
                                        -dbname => $dnadbname );

  $sourcedb->dnadb($dnadb);
}

my $ga = $sourcedb->get_GeneAdaptor;

my @genes;
my @copy_genes;

if ($remove_xrefs && $remove_stable_ids) {
  print STDERR "Fetching genes, removing xrefs and stable ids. (Adaptors and dbIDs also removed.)\n";
} elsif ($remove_xrefs) {
  print STDERR "Fetching genes, removing xrefs, keeping stable ids. (Adaptors and dbIDs also removed.)\n";
} elsif ($remove_stable_ids) {
  print STDERR "Fetching genes, removing stable ids, keeping xrefs. (Adaptors and dbIDs also removed.)\n";
} else {
  print STDERR "Fetching genes, keeping xrefs and stable IDs. (Adaptors and dbIDs removed.)\n";
}

if ($infile) {
  open(INFILE, "<$infile") or die ("Can't read $infile $! \n");

  my $i = 0 ; 
  while(<INFILE>){
    chomp;
    $i++;
    my $gene_id = $_;
    my $gene;

    if ($stable_id) {
      $gene = $ga->fetch_by_stable_id($gene_id);
      $gene->load();
    } else {
      $gene = $ga->fetch_by_dbID($gene_id);
      $gene->load();
    }
    print "fetched $i genes\n" if $verbose ;
    empty_Gene($gene, $remove_stable_ids, $remove_xrefs);
    push( @copy_genes, $gene );
  }
  close(INFILE);
} elsif ($all) {
  my @genes = @{$ga->fetch_all()};   
  my $i = 0 ; 
  foreach my $gene (@genes ) {  
    $gene->load();
    $i++;
    print "fetched $i / " . scalar(@genes) . " \n" if $verbose ; 
    empty_Gene($gene, $remove_stable_ids, $remove_xrefs);
    push( @copy_genes, $gene );
  }
}

my $outdb;
if ($out_config_name) {
  $outdb = get_db_adaptor_by_string($out_config_name);
} else {
  $outdb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $outhost,
                                        -user   => $outuser,
                                        -pass   => $outpass,
                                        -port   => $outport,
                                        -dbname => $outdbname );
}

if ($split) {
  foreach my $gene (@copy_genes) {
    push( @genes,
          @{convert_to_genes( $gene->get_all_Transcripts, $gene->analysis,
                              $gene->biotype) }
        );
  }
} else {
  @genes = @copy_genes;
}
print STDERR "Fetched ".scalar(@genes)." genes\n" if $verbose ; 

my $outga = $outdb->get_GeneAdaptor;

if (defined $logic) {
  my $analysis; 
  $analysis = $outdb->get_AnalysisAdaptor->fetch_by_logic_name($logic);
  if ( !defined $analysis ) {
    $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => $logic, );
  }


  foreach my $g (@genes) {
    $g->analysis($analysis);
    foreach my $t (@{$g->get_all_Transcripts}) {
      $t->analysis($analysis);
    }
  }
}

my $si = 0 ; 
foreach my $gene (@genes) { 
  $si++; 
  my $old_stable_id = $gene->stable_id ; 
  print "transforming $old_stable_id\n" if $verbose ; 
  $gene->load();#fully_load_Gene($gene);
  if ($transform_to) {
    my $transformed_gene = $gene->transform( $transform_to , $transform_to_version );
    $gene = $transformed_gene ;
  } 
  if ( $gene ) { 
    empty_Gene($gene, $remove_stable_ids, $remove_xrefs);
    $outga->store($gene);   
    print "stored $si / " . scalar(@genes) . " \n" if $verbose ; 
  } else { 
     print STDERR "gene $old_stable_id did not tranform\n";
  } 
}
