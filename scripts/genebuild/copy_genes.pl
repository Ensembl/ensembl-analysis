#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/scripts/genebuild/copy_genes.pl,v $
# $Revision: 1.30 $

=head1 NAME

  copy_genes.pl

=head1 DESCRIPTION

This script takes database options and a file of gene ids or stable_ids
and copies them between two databases.  It can, if asked, split
multi-transcript genes into single genes.

The source database (from which you copy genes FROM) must contain DNA,
or else the script will throw.  If your source database does not contain
DNA, then provide details of a DNA database using the options --dnauser,
--dnahost, --dnadbname etc.

When using the in and out_config_name options it reads the equivalent
database from the Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases
file.

=head1 OPTIONS

  --all                 This option will copy all genes from
                        $sourcedbname to the output database.

  --split               This option will split multi-transcript genes
                        into single-transcript genes.

  --logic               This option will change the analysis of the
                        genes being written to the output database to an
                        analysis with the specified logic_name.

  --remove_xrefs        Using this flag will remove xrefs from genes,
                        transcripts, and from translations.

  --remove_stable_ids   Using this flag will remove stable IDs from
                        genes, transcripts, translations, and from
                        exons.

  --transform_to        Transforms the genes from one coordinate system
                        to the given one.  The value must be on the
                        format "coord_system_name:version".

  --stable_id           Flag for indicating that input file contains
                        stable IDs and not gene IDs.

  --file                Read gene dbIDs out of a supplied file and only
                        copy these genes.

=head1 EXAMPLE

  perl copy_genes.pl \
    --source_config_name=CONSENSUS_DB \
    --target_config_name=UTR_DB \
    --file=gene_ids_to_copy

  or

  perl copy_genes.pl \
    --sourcehost=srchost --sourceuser=ensro \
    --sourcedbname=est_db \
    --targethost=trghost --targetuser=user --targetpass=XXX \
    --targetdbname=utr_db \
    --file=gene_ids_to_copy --stable_id


In the options, 'source' is synonymous with 'in' (--in_config_name is
the same option as --source_config_name etc.) and 'target' is synonymous
with 'out' (--targetdbname is the same as --outdbname).

=cut

use strict;
use warnings;

use Carp;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning verbose);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
use Bio::EnsEMBL::Analysis::Tools::Utilities;

my $sourcehost   = '';
my $sourceuser   = 'ensro';
my $sourcedbname = '';
my $sourceport   = 3306;

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

verbose('EXCEPTION');

GetOptions( 'inhost|sourcehost:s'                  => \$sourcehost,
            'inuser|sourceuser:s'                  => \$sourceuser,
            'indbname|sourcedbname:s'              => \$sourcedbname,
            'inport|sourceport:n'                  => \$sourceport,
            'in_config_name|source_config_name:s'  => \$in_config_name,
            'out_config_name|target_config_name:s' => \$out_config_name,
            'outhost|targethost:s'                 => \$outhost,
            'outuser|targetuser:s'                 => \$outuser,
            'outpass|targetpass:s'                 => \$outpass,
            'outdbname|targetdbname:s'             => \$outdbname,
            'outport|targetport:n'                 => \$outport,
            'dnahost:s'                            => \$dnahost,
            'dnauser:s'                            => \$dnauser,
            'dnadbname:s'                          => \$dnadbname,
            'dnaport:n'                            => \$dnaport,
            'logic:s'                              => \$logic,
            'split!'                               => \$split,
            'all!'                                 => \$all,
            'remove_xrefs!'                        => \$remove_xrefs,
            'remove_stable_ids!' => \$remove_stable_ids,
            'transform_to:s'     => \$transform_to,
            'verbose!'           => \$verbose,
            'stable_id!'         => \$stable_id,
            'file:s'             => \$infile ) ||
  throw("Error while parsing command line options");

if ($verbose) {
  verbose('WARNING');
}

my $transform_to_version;

my $attach_dna_db = 0;
if ($transform_to) {
  $attach_dna_db = 1;
  if ( $transform_to =~ m/:/ ) {
    ( $transform_to, $transform_to_version ) = split /:/, $transform_to;
  }
}

if ( $all && defined($infile) ) {
  throw("Specify either -all or -file");
}
elsif ( !$all && !defined($infile) ) {
  throw( "Specify -all (to copy all genes from $sourcedbname) " .
         "or -file (to copy a list of gene_ids)" );
}

#print "Connecting to ".$sourcedbname." at ".$sourcehost." ".$sourceuser."\n";

my $sourcedb;
my $dnadb;

if (defined($in_config_name)) {
  $sourcedb =
    get_db_adaptor_by_string( $in_config_name,
                              0, 0,
                              { -do_not_attach_dna_db => !$attach_dna_db
                              } );

  if ( $attach_dna_db && $dnadbname ) {
    $dnadb =
      new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                          -user   => $dnauser,
                                          -port   => $dnaport,
                                          -dbname => $dnadbname );

    $sourcedb->dnadb($dnadb);
  }
}
else {
  $sourcedb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $sourcehost,
                                        -user   => $sourceuser,
                                        -port   => $sourceport,
                                        -dbname => $sourcedbname );
  if ($attach_dna_db) {
    if ( !defined($dnadbname) ) {
      my $dna_query = q(SELECT count(1) FROM dna);
      my $dna_sth   = $sourcedb->dbc()->prepare($dna_query);
      $dna_sth->execute();

      my $dna_count = $dna_sth->fetchrow();

      if ( !$dna_count ) {
        croak( "\nYour source database does not contain DNA. " .
               "Please provide a database with genomic sequences.\n" );
      }
    }
    else {
      $dnadb =
        new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                            -user   => $dnauser,
                                            -port   => $dnaport,
                                            -dbname => $dnadbname );

      $sourcedb->dnadb($dnadb);
    }
  }
} ## end else [ if ($in_config_name) ]

my $outdb;
if ($out_config_name) {
  $outdb = get_db_adaptor_by_string( $out_config_name, 0, 0,
                         { -do_not_attach_dna_db => !$attach_dna_db } );
}
else {
  $outdb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $outhost,
                                        -user   => $outuser,
                                        -pass   => $outpass,
                                        -port   => $outport,
                                        -dbname => $outdbname );
}

my $ga = $sourcedb->get_GeneAdaptor();

my @genes;
my @copy_genes;

if ( $remove_xrefs && $remove_stable_ids ) {
  print STDERR "Fetching genes, removing xrefs and stable ids. " .
    "(Adaptors and dbIDs also removed.)\n";
}
elsif ($remove_xrefs) {
  print STDERR "Fetching genes, removing xrefs, keeping stable ids. " .
    "(Adaptors and dbIDs also removed.)\n";
}
elsif ($remove_stable_ids) {
  print STDERR "Fetching genes, removing stable ids, keeping xrefs. " .
    "(Adaptors and dbIDs also removed.)\n";
}
else {
  print STDERR "Fetching genes, keeping xrefs and stable IDs. " .
    "(Adaptors and dbIDs removed.)\n";
}

if ( defined($infile) ) {
  open( INFILE, "<$infile" ) or die("Can't read $infile $! \n");

  my $i = 0;
  while (<INFILE>) {
    chomp;
    $i++;
    my $gene_id = $_;
    my $gene;

    if ($stable_id) {
      $gene = $ga->fetch_by_stable_id($gene_id);
    }
    else {
      $gene = $ga->fetch_by_dbID($gene_id);
    }

    if ($verbose) {
      print "fetched $i genes\n";
    }

    push( @copy_genes, $gene );
  }
  close(INFILE);
} ## end if ( defined($infile) )
elsif ($all) {
  @copy_genes = @{ $ga->fetch_all() };
}

if ($split) {
  foreach my $gene (@copy_genes) {
    push( @genes,
          @{convert_to_genes( $gene->get_all_Transcripts(),
                              $gene->analysis(),
                              $gene->biotype() ) } );
  }
}
else {
  @genes = @copy_genes;
}

print STDERR "Fetched " . scalar(@genes) . " genes\n" if $verbose;

my $outga = $outdb->get_GeneAdaptor();

my $analysis;
if ( defined($logic) ) {
  $analysis = $outdb->get_AnalysisAdaptor->fetch_by_logic_name($logic);
  if ( !defined($analysis) ) {
    $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => $logic, );
  }
}

my $si = 0;
foreach my $gene (@genes) {
  $si++;
  my $old_stable_id = $gene->stable_id();

  print "Loading gene $old_stable_id from source DB\n" if $verbose;

  $gene->load();    # fully_load_Gene($gene);

  if ( defined($logic) ) {
    $gene->analysis($analysis);
    foreach my $t ( @{ $gene->get_all_Transcripts() } ) {
      $t->analysis($analysis);
    }
  }

  if ($transform_to) {
    print "transforming $old_stable_id\n" if $verbose;

    my $transformed_gene =
      $gene->transform( $transform_to, $transform_to_version );

    # only check transform if transform is successful
    if ( defined($transformed_gene) ) {
      check_transform( $gene, $transformed_gene,
                       $transform_to_version );
    }
    $gene = $transformed_gene;
  }

  if ( defined($gene) ) {
    empty_Gene( $gene, $remove_stable_ids, $remove_xrefs );
    $outga->store($gene);
    print "stored $si / " . scalar(@genes) . " \n" if $verbose;
  }
  else {
    print STDERR "gene $old_stable_id did not transform\n";
  }
} ## end foreach my $gene (@genes)

sub check_transform {
  my ( $old_gene, $new_gene, $new_assembly_version ) = @_;

  my $old_transcripts = $old_gene->get_all_Transcripts();
  my $new_transcripts = $new_gene->get_all_Transcripts();

  if ( scalar(@$old_transcripts) != scalar(@$new_transcripts) ) {
    print STDERR "TRANSFORM_CHECK: old gene " .
      $old_gene->stable_id() . " has " .
      ( scalar(@$old_transcripts) ) . " transcripts but new gene has " .
      ( scalar(@$new_transcripts) ) . " transcripts\n";
  }

  # loop through and compare
  foreach my $old_transc ( @{$old_transcripts} ) {
    foreach my $new_transc ( @{$new_transcripts} ) {
      if ( $old_transc->stable_id() ne $new_transc->stable_id() ) {
        next;
      }

      # check number of exons, not sure if this helps
      if ( scalar( @{ $old_transc->get_all_Exons() } ) !=
           scalar( @{ $new_transc->get_all_Exons() } ) )
      {
        print STDERR "TRANSFORM_CHECK: old transcript " .
          $old_transc->stable_id() .
          " has " . ( scalar( @{ $old_transc->get_all_Exons() } ) ) .
          " exons but new transcript has " .
          ( scalar( @{ $new_transc->get_all_Exons() } ) ) . " exons\n";
      }

      # check translation
      if ( defined( $old_transc->translation() ) &&
           !defined($new_assembly_version) )
      {
        # We don't want to do have to deal with this if transforming
        # between assembly _VERSIONS_.
        my $new_translation = $new_transc->translate->seq();
        my $old_translation = $old_transc->translate->seq();

        if ( !defined( $new_transc->translation() ) ||
             $old_translation ne $new_translation )
        {
          print "TRANSFORM_CHECK: " .
            "old translation does not match new translation\n" .
            ">old_" . $old_transc->stable_id() .
            "\n" . $old_translation . "\n" . ">new_" .
            $new_transc->stable_id() . "\n" . $new_translation . "\n";
        }
      }
    } ## end foreach my $new_transc ( @{...})
  } ## end foreach my $old_transc ( @{...})

  return 1;
} ## end sub check_transform
