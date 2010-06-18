#!/usr/local/ensembl/bin/perl

=head1 NAME

  copy_genes.pl

=head1 DESCRIPTION

 This script takes database options and a file of gene ids and copies
 them between two databases. It can if asked split multi transcript
 genes into single genes.

 When using the in and out_config_name options
 it reads the equivalent database from the
 Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases file.

=head1 EXAMPLE

 perl copy_genes.pl -in_config_name COALESCER_DB -out_config_name UTR_DB \
   -file gene_ids_to_copy

 or

 perl copy_genes.pl -dbhost host -dbuser ensro -dbname est_db -outhost host1 \
   -outuser user -outpass **** -outdbname utr_db -file gene_ids_to_copy

=cut

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
use Getopt::Long;
use Bio::EnsEMBL::Analysis::Tools::Utilities;

my $host      = '';
my $user      = '';
my $pass      = undef;
my $dbname    = '';
my $port      = 3306;

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

GetOptions( 'dbhost:s'          => \$host,
            'dbuser:s'          => \$user,
            'dbpass:s'          => \$pass,
            'dbname:s'          => \$dbname,
            'dbport:n'          => \$port,
            'in_config_name:s'  => \$in_config_name,
            'out_config_name:s' => \$out_config_name,
            'outhost:s'         => \$outhost,
            'outuser:s'         => \$outuser,
            'outpass:s'         => \$outpass,
            'outdbname:s'       => \$outdbname,
            'outport:n'         => \$outport,
            'dnahost:s'         => \$dnahost,
            'dnauser:s'         => \$dnauser,
            'dnadbname:s'       => \$dnadbname,
            'dnaport:n'         => \$dnaport,
            'logic:s'           => \$logic,
            'split!'            => \$split,
            'file:s'            => \$infile );

#print "Connecting to ".$dbname." at ".$host." ".$user."\n";
my $db;
if ($in_config_name) {
  $db = get_db_adaptor_by_string($in_config_name);
} else {

  $db =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                        -user   => $user,
                                        -pass   => $pass,
                                        -port   => $port,
                                        -dbname => $dbname );
}

my $dnadb;
if ( !$dnadbname ) {
  my $dna_query = q(SELECT count(1) FROM dna);
  my $dna_sth = $db->dbc->prepare($dna_query);
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

  $db->dnadb($dnadb);
}


my $ga = $db->get_GeneAdaptor;

my @genes;
my @copy_genes;

open(INFILE, "<$infile") or die ("Can't read $infile $! \n");

while(<INFILE>){
  #print "$_";
  chomp;
  my $gene_id = $_;
  my $gene    = $ga->fetch_by_dbID($gene_id);
  empty_Gene($gene);
  push( @copy_genes, $gene );
}
close(INFILE);

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
print STDERR "Fetched ".scalar(@genes)." genes\n";

my $outga = $outdb->get_GeneAdaptor;

if (defined $logic) {
  my $analysis; 
  $analysis = $outdb->get_AnalysisAdaptor->fetch_by_logic_name($logic);
  if (!defined $analysis) {
    $analysis = Bio::EnsEMBL::Analysis->new(
                  -logic_name  => $logic,
                                           ); 
  }


  foreach my $g (@genes) {
    $g->analysis($analysis);
    foreach my $t (@{$g->get_all_Transcripts}) {
      $t->analysis($analysis);
    }
  }
}

foreach my $gene (@genes) {
  fully_load_Gene($gene);
  empty_Gene($gene);
  $outga->store($gene);
}
